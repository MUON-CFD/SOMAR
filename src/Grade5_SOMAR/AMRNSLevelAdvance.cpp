/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2024 Thomas Jefferson University and Arizona State University
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
 *  USA
 *
 *  For up-to-date contact information, please visit the repository homepage,
 *  https://github.com/MUON-CFD/SOMAR.
 ******************************************************************************/
#include "AMRNSLevel.H"
#include "AMRNSLevelF_F.H"
#include "FiniteDiffF_F.H"
#include "memusage.H"
#include "ProblemContext.H"
#include "Convert.H"
#include "Integral.H"
#include "Subspace.H"
#include "IO.H"
#include "Comm.H"
#include "Masks.H"
#include "Debug.H"
#include "MiscUtils.H"

#include <fstream>
#include <iomanip>
#include <chrono>


// -----------------------------------------------------------------------------
// Computes the maximum stable timestep using the current state.
// -----------------------------------------------------------------------------
Real
AMRNSLevel::computeDt()
{
    BEGIN_FLOWCHART();

    // Declare variables.
    const int verbThresh = 0;

    const DisjointBoxLayout& grids     = m_levGeoPtr->getBoxes();
    DataIterator             dit       = grids.dataIterator();
    const RealVect&          dXi       = m_levGeoPtr->getDXi();
    const ProblemContext*    ctx       = ProblemContext::getInstance();
    const Real               maxDt     = ctx->time.maxDt;
    const Real               dtMult    = ctx->time.dtMult;
    const Real               maxDtGrow = ctx->time.maxDtGrow;
    Real newDt;

    // Let the user know what is happening.
    if (s_verbosity >= verbThresh) {
        // Format::pushFlags(pout());
        pout() << Format::indent() << "computeDt on level " << m_level << ":"
               << Format::indent() << Format::scientific << endl;
    }

    // Initialize to max allowed dt.
    {
        if (s_verbosity >= verbThresh) {
            pout() << "dt limit maxDt = " << maxDt << endl;
        }
        newDt = maxDt;
        CH_assert(newDt > 0.0);
    }

    // Advective stability limit...
    {
        Real advDt = 1.0e100;
        for (dit.reset(); dit.ok(); ++dit) {
            const GeoSourceInterface& geoSrc    = m_levGeoPtr->getGeoSource();
            const FluxBox&            cartVelFB = m_statePtr->vel[dit];
            const Box&                ccValid   = grids[dit];

            FArrayBox ccDXidxFAB(ccValid, SpaceDim);
            for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
                geoSrc.fill_dXidx(ccDXidxFAB, fcDir, fcDir, dXi);
            }

            FORT_COMPUTEADVECTIVEDT (
                CHF_REAL(advDt),
                CHF_CONST_FRA1(cartVelFB[0],0),
                CHF_CONST_FRA1(cartVelFB[1],0),
                CHF_CONST_FRA1(cartVelFB[SpaceDim-1],0),
                CHF_CONST_FRA(ccDXidxFAB),
                CHF_BOX(ccValid),
                CHF_CONST_REALVECT(dXi));
        }

        Comm::reduce(advDt, MPI_MIN);
        advDt *= m_parkPtr->ERKStabilityIm();

        if (s_verbosity >= verbThresh) {
            pout() << "dt advective limit = " << advDt << endl;
        }
        newDt = min(newDt, advDt);
    } // end advective limiting.


    // Viscous and diffusive stability limits...
    if (!ctx->rhs.doImplicitDiffusion &&
        (ctx->rhs.doViscousForcing || ctx->rhs.doTemperatureDiffusion ||
         ctx->rhs.doSalinityDiffusion || ctx->rhs.doScalarDiffusion)) {
        Real viscDt = 1.0e100;
        Real diffDt = 1.0e100;

        for (dit.reset(); dit.ok(); ++dit) {
            const FArrayBox& nuTFAB = m_statePtr->eddyNu[dit];
            const Box&       valid  = grids[dit];

            // Compute cell widths.
            Tuple<FArrayBox, CH_SPACEDIM> dxFAB;
            for (int d = 0; d < SpaceDim; ++d) {
                const IntVect e = BASISV(d);

                Box region = Subspace::flattenBox(valid, e);
                region.surroundingNodes(d);
                FArrayBox xFAB(region, 1);
                m_levGeoPtr->getGeoSource().fill_physCoor(xFAB, 0, d, dXi);

                region.enclosedCells(d);
                dxFAB[d].define(region, 1);
                const Real dummyDx = 1.0;

                FiniteDiff::partialD(dxFAB[d], 0, region, xFAB, 0, d, dummyDx);
            }

            // Compute delta2 = 1 / (1/dx^2 + 1/dy^2 + 1/dz^2).
            FArrayBox delta2FAB(valid, 1);
            FORT_COMPUTEDELTA2FORDIFFUSIVEDTLIMIT(
                CHF_FRA1(delta2FAB, 0),
                CHF_CONST_FRA1(dxFAB[0], 0),
                CHF_CONST_FRA1(dxFAB[1], 0),
                CHF_CONST_FRA1(dxFAB[SpaceDim - 1], 0),
                CHF_BOX(valid));

            // Do viscous limiting.
            if (ctx->rhs.doViscousForcing) {
                const Real nu           = ctx->rhs.nu;
                const Real dummyPrandtl = 1.0;

                FORT_COMPUTEDIFFUSIVEDT(
                    CHF_REAL(viscDt),
                    CHF_FRA1(delta2FAB, 0),
                    CHF_CONST_REAL(nu),
                    CHF_CONST_FRA1(nuTFAB, 0),
                    CHF_CONST_REAL(dummyPrandtl),
                    CHF_BOX(valid));
            }

            // Do diffusive limiting -- Temperature.
            if (ctx->rhs.doTemperatureDiffusion) {
                const Real kappa        = ctx->rhs.TKappa;
                const Real eddyPrandtl  = ctx->rhs.eddyPrandtlT;

                FORT_COMPUTEDIFFUSIVEDT(
                    CHF_REAL(diffDt),
                    CHF_FRA1(delta2FAB, 0),
                    CHF_CONST_REAL(kappa),
                    CHF_CONST_FRA1(nuTFAB, 0),
                    CHF_CONST_REAL(eddyPrandtl),
                    CHF_BOX(valid));
            }

            // Do diffusive limiting -- Salinity.
            if (ctx->rhs.doSalinityDiffusion) {
                const Real kappa        = ctx->rhs.SKappa;
                const Real eddyPrandtl  = ctx->rhs.eddyPrandtlS;

                FORT_COMPUTEDIFFUSIVEDT(
                    CHF_REAL(diffDt),
                    CHF_FRA1(delta2FAB, 0),
                    CHF_CONST_REAL(kappa),
                    CHF_CONST_FRA1(nuTFAB, 0),
                    CHF_CONST_REAL(eddyPrandtl),
                    CHF_BOX(valid));
            }

            // Do diffusive limiting -- User-defined scalars.
            if (ctx->rhs.doScalarDiffusion && this->numScalars() > 0) {
                const unsigned int numComps = (unsigned int)this->numScalars();

                for (unsigned int comp = 0; comp < numComps; ++comp) {
                    const Real sKappa = ctx->rhs.getScalarsKappa(comp);
                    const Real eddyPr = ctx->rhs.getEddyPrandtlScalars(comp);

                    FORT_COMPUTEDIFFUSIVEDT(
                        CHF_REAL(diffDt),
                        CHF_FRA1(delta2FAB, 0),
                        CHF_CONST_REAL(sKappa),
                        CHF_CONST_FRA1(nuTFAB, 0),
                        CHF_CONST_REAL(eddyPr),
                        CHF_BOX(valid));
                }
            }
        } // dit

        Comm::reduce(viscDt, MPI_MIN);
        Comm::reduce(diffDt, MPI_MIN);

        viscDt *= m_parkPtr->ERKStabilityRe();
        diffDt *= m_parkPtr->ERKStabilityRe();

        if (s_verbosity >= verbThresh) {
            if (viscDt < 1.0e100) {
                pout() << "dt viscous limit = " << viscDt << endl;
            }

            if (diffDt < 1.0e100) {
                pout() << "dt diffusive limit = " << diffDt << endl;
            }
        }
        newDt = min(newDt, viscDt);
        newDt = min(newDt, diffDt);
    } // end viscous / diffusive limiting.


    // // Internal wave speed limit
    // const bool doGravityForcing = ctx->rhs.doGravityForcing;
    // if (doGravityForcing && m_hasStrat) {
    //     const LevelData<FluxBox>& advVel = m_statePtr->vel;

    //     // Loop over grids and compute minDt.
    //     Real iwsDt = 1.0e100;
    //     for (dit.reset(); dit.ok(); ++dit) {
    //         const Box& valid = grids[dit];

    //         D_TERM(
    //         FArrayBox c0iFAB(valid, SpaceDim);,
    //         m_levGeoPtr->getGeoSource().fill_dXidx(c0iFAB, 0, 0, dXi, m_c1);,
    //         m_levGeoPtr->getGeoSource().fill_dXidx(c0iFAB, 1, 1, dXi, m_c1);)
    //         c0iFAB.setVal(0.0, SpaceDim - 1);

    //         FArrayBox ccVelFAB(valid, SpaceDim);
    //         for (int dir = 0; dir < SpaceDim; ++dir) {
    //             Convert::Simple(ccVelFAB, dir, valid, advVel[dit][dir], 0);
    //         }

    //         FORT_COMPUTEMINBVDT(
    //             CHF_REAL(iwsDt),
    //             CHF_CONST_FRA(c0iFAB),
    //             CHF_CONST_FRA(ccVelFAB),
    //             CHF_CONST_REALVECT(dXi),
    //             CHF_BOX(valid));

    //         CH_assert(iwsDt >= 0.0);
    //     }

    //     newDt = min(newDt, iwsDt);

    //     if (s_verbosity >= verbThresh) {
    //         pout() << "dt internal wave speed limit = " << iwsDt << endl;
    //     }
    // }

    // Embedded RK controllers
    if (ctx->time.useElementaryController ||
        ctx->time.usePIController ||
        ctx->time.usePIDController)
    {
        const auto velNorm = Analysis::pNorm(*m_velPtr, 0);
        const Real velScale = std::max({D_DECL(velNorm[0], velNorm[1], velNorm[2])});
        // const Real tol = std::max(ctx->time.absTol, ctx->time.relTol * velScale);
        const Real tol = ctx->time.absTol + ctx->time.relTol * velScale;

        const Real controllerDt = m_parkPtr->controllerDt(
            tol,
            ctx->rhs.doImplicitDiffusion,
            ctx->time.useElementaryController,
            ctx->time.usePIController,
            ctx->time.usePIDController);

        if (s_verbosity >= verbThresh) {
            pout() << Format::scientific
                    << "dt controller limit = " << controllerDt
                    << endl;
        }

        newDt = std::min(newDt, controllerDt);

// CHECKPOINT();
// if (m_time > 0.0 && controllerDt < 1000.0) {
//     pout() << "Forcing controller dt.\n";
//     newDt = controllerDt;
// } else {
//     IO::tout(0) << Format::hired << "Not forcing controller dt.\n" << Format::none;
// }
    }

    // Communicate: compute the minimum over all procs.
    Comm::reduce(newDt, MPI_MIN);

    // Report
    if (s_verbosity >= verbThresh) {
        pout() << "New dt (before dtMult applied) = " << newDt
               << Format::unindent << Format::unindent << endl;
        // Format::popFlags(pout());
    }

    // Scale by dtMult.
    newDt *= dtMult;

    // Ensure dt is not growing too fast.
    if (m_dt > smallReal) {
        if (newDt > maxDtGrow * m_dt) {
            newDt = maxDtGrow * m_dt;
            pout() << "Dt growth restricted by maxDtGrow = " << maxDtGrow << endl;
        }
    }


    CH_assert(newDt > -smallReal);
    if (newDt < 0.0) newDt = smallReal; // Let AnisotropicAMR deal with this.

    CH_assert(newDt <= maxDt);


    // Now that dt has been calculated to keep things stable and accurate,
    // check if we will overstep the simulation's end time. If so, reduce dt.
    // NOTE: I add smallReal to ensure the AMR driver knows we are at the last
    // timestep without running into rounding errors.
    if (m_time + newDt > ctx->time.stopTime) {
        newDt = ctx->time.stopTime - m_time + smallReal;

        if (s_verbosity >= verbThresh) {
            pout() << "\tNew dt reduced to " << newDt
                   << " to match final simulation time." << endl;
        }
    }

    // static std::vector<Real> vdt(1, newDt);
    // static Real lastTime = 0.0;
    // if (RealCmp::neq(m_time, lastTime)) {
    //     vdt.push_back(newDt);
    // } else {
    //     vdt.back() = newDt;
    // }
    // POUT(vdt);

    return newDt;
}


// -----------------------------------------------------------------------------
// advance by one timestep (returns maximum allowable timestep)
// -----------------------------------------------------------------------------
Real
AMRNSLevel::advance()
{
    BEGIN_FLOWCHART();

    // Has the halt flag been thrown? If so, abort!
    {
        ifstream ifile("halt");
        if (ifile) {
            remove("halt");
            MayDay::Abort("Halt flag was set. Terminating run...");
        }
    }

    const ProblemContext*    ctx   = ProblemContext::getInstance();
    const DisjointBoxLayout& grids = this->getBoxes();
    DataIterator             dit   = grids.dataIterator();

    // Notify user and set pout indentation.
    pout() << Format::fixed;
    if (s_verbosity >= 2) {
        pout() << "Timestep on level " << m_level << ":\n";
    }
    pout() << Format::indent() << flush;

    const auto start = std::chrono::high_resolution_clock::now();

    // Set time variables.
    const Real dt = m_dt;
    const Real oldTime = m_time;
    const Real newTime = m_time + m_dt;
    CH_verify(!RealCmp::isZero(dt));

    // Initialize the flux register between this and the finer level.
    if (m_velFluxRegPtr) m_velFluxRegPtr->setToZero();
    if (m_qFluxRegPtr) m_qFluxRegPtr->setToZero();


    // Compute and apply the sponge forcing via Forward Euler.
    // We cannot do this within setExplicitRHS because RK methods (including
    // RK3) do not guarantee dt = 0 at each stage.
    //
    // Warning: Do this before calling PARK::advance or else the final
    // divergence will not be zero.
    if (ctx->rhs.doSpongeForcing) {
        debugCheckValidFaceOverlap(*m_velPtr);

        LevelData<FluxBox> spongekvel;
        LevelData<FArrayBox> spongekq;

        spongekvel.define(grids, m_velPtr->nComp(), m_velPtr->ghostVect());
        spongekq.define(grids, m_qPtr->nComp(), m_qPtr->ghostVect());

        setValLevel(spongekvel, 0.0);
        setValLevel(spongekq, 0.0);

        if (!RealCmp::isZero(dt)) {
            Real invTimeScale = 1.0 / (dt * ctx->rhs.spongeTimeCoeff);
            this->addSpongeForcing(spongekvel,
                                   spongekq,
                                   m_statePtr->vel,
                                   m_statePtr->p,
                                   m_statePtr->q,
                                   oldTime,
                                   invTimeScale);

            // static int fidx = 0;
            // char fn[80];
            // sprintf(fn, "force_%04d.hdf5", fidx++);
            // IO::writeHDF5(fn, spongekq, *m_levGeoPtr);
            // barrier();
            // CH_verify(false);

            for (dit.reset(); dit.ok(); ++dit) {
                for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                    (*m_velPtr)[dit][velComp].plus(spongekvel[dit][velComp], dt);
                }
                (*m_qPtr)[dit].plus(spongekq[dit], dt);
            }
        }

        this->setBC(*m_statePtr, oldTime);
    }

    debugCheckValidFaceOverlap(*m_velPtr);

    // Do the RK timestep.
    m_parkPtr->advance(*m_velPtr, *m_pPtr, *m_qPtr, oldTime, dt, this);
    this->time(newTime);

    debugCheckValidFaceOverlap(*m_velPtr);

    // {
    //     const RealVect& dx = m_levGeoPtr->getDXi();
    //     const Real nu = ctx->rhs.nu;

    //     LevelData<FluxBox> vel(grids, 1);
    //     m_velPtr->copyTo(vel);

    //     Real epsilon = Python.PythonReturnFunction<Real>(
    //         "Diagnostic", "epsilon", vel, nu, dx);
    //     pout() << "\nepsilon = " << epsilon << endl;
    // }


    // Report memory usage & timing info
    if (s_verbosity >= 2) {
        const Real memory = get_memory_usage_from_OS();
        pout() << "Current memory usage = " << memory << "MB" << endl;

        const auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        pout() << "Timestep computation time = " << elapsed.count() << " s\n";
    }

    // Restore pout indentation.
    pout() << Format::unindent << flush;

    return m_dt; // This is actually not used by the AMR driver.
}


// -----------------------------------------------------------------------------
// Things to do after a timestep.
// a_step is the step number that is being completed.
// Unlike in AnisotropicAMR.cpp, a_step starts counting from 1.
// Called from finest -> coarsest synced levels.
// -----------------------------------------------------------------------------
void
AMRNSLevel::postTimeStep(int a_step)
{
    BEGIN_FLOWCHART();

    if (s_verbosity >= 3) {
        pout() << "postTimeStep on level " << m_level
            << " with dt = " << this->getLevel(0)->dt() << endl;
    }

    TODONOTE(
        "We should average and project more carefully since this is called "
        "from top to bottom.");

    const ProblemContext* ctx = ProblemContext::getInstance();

    // Check if this is the coarsest synchronized level, lbase.
    Real crseTime = -1.0;
    {
        const AMRNSLevel* crsePtr = this->crseNSPtr();
        if (crsePtr) crseTime = crsePtr->m_time;
    }
    const bool atBaseLevel = RealCmp::neq(crseTime, m_time);

    if (!this->isFinestLevel() && atBaseLevel) {
        // 1. Refluxing
        {
            // Gather AMR state.
            Vector<LevelData<FluxBox>*>   amrVel;
            Vector<LevelData<FArrayBox>*> amrQ;
            this->allocateAndAliasVel(
                amrVel, m_time, m_statePtr->vel.interval(), m_level);
            this->allocateAndAliasScalars(
                amrQ, m_time, m_statePtr->q.interval(), m_level);

            // Reflux.
            this->simpleRefluxing(amrVel, amrQ);
            // this->explicitRefluxing(amrVel, amrQ);

            // Free memory
            this->deallocate(amrQ);
            this->deallocate(amrVel);
        }

        // 2. Sync projection
        if (ctx->proj.doSyncProj) {
            // Okay, I bet you're wondering why this is false!
            // Well, it shouldn't be. This is just temporary to that the
            // pressure is not updated. In the future, I think we need to store
            // the sync pressure increment so that we can set CFI BCs properly.
            const bool sync = false;

            pout() << "\nSynchronization projection:" << endl;
            this->projectDownToThis(sync);
            pout() << endl;
        }

        // 3. Average down
        this->averageDownToThis();

        // 4. Reset BCs
        this->setBCsDownToThis();
    } // end if not finest level.

    // Finally, write diagnostic info to terminal.
    this->printDiagnostics(a_step, false);
}


// -----------------------------------------------------------------------------
// Refluxes via Q = Q + dQ where
//  dQ are the contents of the flux registers.
// You must call this from lbase = the coarsest refluxed level.
// -----------------------------------------------------------------------------
void
AMRNSLevel::simpleRefluxing(Vector<LevelData<FluxBox>*>&   a_amrVel,
                            Vector<LevelData<FArrayBox>*>& a_amrQ)
{
    const ProblemContext* ctx   = ProblemContext::getInstance();
    const int             lbase = m_level;
    const int             lmax  = a_amrVel.size() - 1;

    const AMRNSLevel* levPtr = this;
    for (int l = lbase; l < lmax; ++l) {
        if (ctx->rhs.doMomAdvRefluxing && levPtr->m_velFluxRegPtr) {
            levPtr->m_velFluxRegPtr->reflux(*a_amrVel[l], 1.0);
        }

        if (levPtr->m_qFluxRegPtr) {
            levPtr->m_qFluxRegPtr->reflux(*a_amrQ[l], 1.0);
        }

        levPtr = levPtr->fineNSPtr();
    }
}


// -----------------------------------------------------------------------------
// Refluxes via Q = Q + [1 + dt*D](dQ) where
//  D is the diffusive operator and
//  dQ are the contents of the flux register.
// You must call this from lbase = the coarsest refluxed level.
// -----------------------------------------------------------------------------
void
AMRNSLevel::explicitRefluxing(Vector<LevelData<FluxBox>*>&   a_amrVel,
                              Vector<LevelData<FArrayBox>*>& a_amrQ)
{
    const int             lbase      = m_level;
    const int             lmax       = a_amrVel.size() - 1;
    const Real            refluxDt   = m_dt;
    const Real            refluxTime = m_time;
    const ProblemContext* ctx        = ProblemContext::getInstance();

    for (int l = lbase; l < lmax; ++l) {
        const AMRNSLevel*        levPtr = this->getLevel(l);
        const DisjointBoxLayout& grids  = levPtr->getBoxes();
        DataIterator             dit    = grids.dataIterator();
        const LevelGeometry&     levGeo = *(levPtr->m_levGeoPtr);

        if (levPtr->m_velFluxRegPtr) {
            LevelData<FluxBox>& cartVel     = *a_amrVel[l];
            const int           numVelComps = cartVel.nComp();
            const IntVect&      ghostVect   = cartVel.ghostVect();

            // Put reflux increment into dedicated holder.
            LevelData<FluxBox> dvel(grids, numVelComps, ghostVect);
            setValLevel(dvel, 0.0);
            levPtr->m_velFluxRegPtr->reflux(dvel, 1.0);
            levPtr->setVelBC(dvel, refluxTime, true);

            // Diffuse the momentum refluxing increment, if necessary.
            if (ctx->rhs.doViscousForcing) {
                const RealVect              nu = RealVect::Unit * ctx->rhs.nu;
                const LevelData<FArrayBox>& eddyNu = levPtr->m_statePtr->eddyNu;

                // Diffuse.
                StaggeredFluxLD    velFlux(grids);
                LevelData<FluxBox> velDiv(grids, 1, IntVect::Unit);
                levPtr->computeMomentumDiffusion(
                    velDiv, velFlux, dvel, nu, eddyNu);

                // Add to dvel.
                for (dit.reset(); dit.ok(); ++dit) {
                    for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                        dvel[dit][velComp].plus(velDiv[dit][velComp], refluxDt);
                    }
                } // dit

                levPtr->setVelBC(dvel, refluxTime, true);
            } // if diffusing

            // Update velocity.
            for (dit.reset(); dit.ok(); ++dit) {
                for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                    levGeo.divByJ(dvel[dit][velComp], dit());
                    cartVel[dit][velComp].plus(dvel[dit][velComp], 1.0);
                }
            } // dit
        } // if refluxing the momentum

        if (levPtr->m_qFluxRegPtr) {
            const int                   numQComps = a_amrQ[l]->nComp();
            const IntVect&              ghostVect = a_amrQ[l]->ghostVect();
            const LevelData<FArrayBox>& eddyNu = levPtr->m_statePtr->eddyNu;

            // Define workspace.
            LevelData<FluxBox>   qFlux(grids, 1);
            LevelData<FArrayBox> qDiv(grids, 1);

            // Put reflux increment into dedicated holder.
            LevelData<FArrayBox> dq(grids, numQComps, ghostVect);
            setValLevel(dq, 0.0);
            levPtr->m_qFluxRegPtr->reflux(dq, 1.0);
            levPtr->setScalarBC(dq, refluxTime, true);

            // Diffuse the temperature refluxing increment.
            if (ctx->rhs.doTemperatureDiffusion) {
                const Real      kappa        = ctx->rhs.TKappa;
                const Real      eddyPrandtlT = ctx->rhs.eddyPrandtlT;
                const Interval& ivl          = levPtr->m_statePtr->TInterval;
                const int       comp         = levPtr->m_statePtr->TComp;

                LevelData<FArrayBox> dqComp;
                aliasLevelData(dqComp, &dq, ivl);

                levPtr->setTemperatureBC(dqComp, refluxTime, true);
                levPtr->computeScalarDiffusion(
                    qDiv, qFlux, dqComp, kappa, eddyNu, eddyPrandtlT);

                for (dit.reset(); dit.ok(); ++dit) {
                    levGeo.divByJ(qDiv[dit], dit());
                    dq[dit].plus(qDiv[dit], refluxDt, 0, comp, 1);
                }
            }

            // Diffuse the salinity refluxing increment.
            if (ctx->rhs.doSalinityDiffusion) {
                const Real      kappa        = ctx->rhs.SKappa;
                const Real      eddyPrandtlS = ctx->rhs.eddyPrandtlS;
                const Interval& ivl          = levPtr->m_statePtr->SInterval;
                const int       comp         = levPtr->m_statePtr->SComp;

                LevelData<FArrayBox> dqComp;
                aliasLevelData(dqComp, &dq, ivl);

                levPtr->setSalinityBC(dqComp, refluxTime, true);
                levPtr->computeScalarDiffusion(
                    qDiv, qFlux, dqComp, kappa, eddyNu, eddyPrandtlS);

                for (dit.reset(); dit.ok(); ++dit) {
                    levGeo.divByJ(qDiv[dit], dit());
                    dq[dit].plus(qDiv[dit], refluxDt, 0, comp, 1);
                }
            }

            // Diffuse the scalar refluxing increments.
            if (ctx->rhs.doScalarDiffusion && levPtr->numScalars() > 0) {
                const Interval& ivl       = levPtr->m_statePtr->scalarsInterval;
                const size_t    startComp = ivl.begin();
                const size_t    endComp   = ivl.end();
                CH_assert(ivl.size() == levPtr->numScalars());

                {
                    LevelData<FArrayBox> dqComps;
                    aliasLevelData(dqComps, &dq, ivl);
                    levPtr->setScalarBC(dqComps, refluxTime, true);
                }

                for (size_t comp = startComp; comp <= endComp; ++comp) {
                    const Real sKappa = ctx->rhs.getScalarsKappa(comp);
                    const Real eddyPr = ctx->rhs.getEddyPrandtlScalars(comp);

                    LevelData<FArrayBox> dqComp;
                    aliasLevelData(dqComp, &dq, Interval(comp, comp));

                    levPtr->computeScalarDiffusion(
                        qDiv, qFlux, dqComp, sKappa, eddyNu, eddyPr);

                    for (dit.reset(); dit.ok(); ++dit) {
                        levGeo.divByJ(qDiv[dit], dit());
                        dqComp[dit].plus(qDiv[dit], refluxDt);
                    }
                } // end loops over scalar comps (comp)
            }

            // Reflux q.
            for (dit.reset(); dit.ok(); ++dit) {
                (*a_amrQ[l])[dit].plus(dq[dit], 1.0);
            }
        } // end if flux reg defined
    } // end loop over refluxed levels
}


// -----------------------------------------------------------------------------
// Writes basic simulation info to the terminal.
// -----------------------------------------------------------------------------
void
AMRNSLevel::printDiagnostics(const int  a_step,
                             const bool a_writeBanner) const
{
    BEGIN_FLOWCHART();

    // Banner
    const int stepWidth  = 6;
    const int levelWidth = 6;
    const int timeWidth  = 10;
    const int dtWidth    = 10;
    const int momWidth   = 15;
    const int eWidth     = 15;
    const int elWidth    = 15;

    if (a_writeBanner) {
        const std::string energyStr = m_hasStrat ? "KE + APE" : "KE + GPE";

        IO::tout(0) << Format::pushFlags;
        IO::tout(0) << setfill(' ') << std::right << Format::hiwhite
                    << std::setw(stepWidth) << "Step"
                    << Format::white << "│" << Format::hiwhite
                    << std::setw(levelWidth) << "Level"
                    << Format::white << "│" << Format::hiwhite
                    << std::setw(timeWidth) << "time" << std::setw(dtWidth)
                    << "dt"
                    << Format::white << "│" << Format::hiwhite
                    D_TERM(<< std::setw(momWidth) << "x_mom",
                           << std::setw(momWidth) << "y_mom",
                           << std::setw(momWidth) << "z_mom")
                    << Format::white << "│" << Format::hiwhite
                    << std::setw(eWidth) << energyStr
                    << std::setw(elWidth) << "\u0394E/E"
                    << endl;
        IO::tout(0) << Format::none << "──────"
                    << "┼"
                    << "──────"
                    << "┼"
                    << "──────────"
                    << "──────────"
                    << "┼" D_TERM(<< "───────────────",
                                  << "───────────────",
                                  << "───────────────")
                    << "┼"
                    << "───────────────"
                    << "───────────────" << endl;
        IO::tout(0) << Format::none << Format::popFlags;
    }

    // Diagnostics
    if (m_level > 0) {
        IO::tout(0) << Format::pushFlags;
        IO::tout(0) << setfill(' ') << std::right
                    << Format::fixed << std::setw(stepWidth) << ""
                    << Format::white << "│" << Format::hiblack;

        IO::tout(0) << setfill(' ') << std::right
                    << Format::fixed
                    << Format::hiblack
                    << std::setw(levelWidth) << m_level
                    << Format::white
                    << "│";

        if (a_step == 0 || (1.0e-5 < abs(m_time) && abs(m_time) < 1.0e5)) {
            IO::tout(0) << Format::fixed
                        << std::setprecision(5)
                        << std::setw(timeWidth) << m_time;
        } else {
            IO::tout(0) << Format::scientific
                        << std::setprecision(3)
                        << std::setw(timeWidth) << m_time;
        }

        if (a_step == 0 || (1.0e-5 < abs(m_dt) && abs(m_dt) < 1.0e5)) {
            IO::tout(0) << Format::fixed
                        << std::setprecision(5)
                        << std::setw(dtWidth) << m_dt;
        } else {
            IO::tout(0) << Format::scientific
                        << std::setprecision(3)
                        << std::setw(dtWidth) << m_dt;
        }

        IO::tout(0) << Format::white
                    << "│"
                    D_TERM(
                    << Format::scientific << std::setprecision(5) << std::setw(momWidth) << "",
                    << Format::scientific << std::setprecision(5) << std::setw(momWidth) << "",
                    << Format::scientific << std::setprecision(5) << std::setw(momWidth) << "")
                    << "│"
                    << Format::scientific << std::setprecision(5) << std::setw(eWidth) << ""
                    << std::setw(elWidth) << ""
                    << endl;
        IO::tout(0) << Format::none << Format::popFlags;

    } else if (m_level == 0) {
        // Gather:
        //  amrVel     = composite CC momentum
        //  amrE       = composite total energy
        //  amrT       = composite temperature
        //  amrS       = composite salinity
        //  amrScalars = composite scalars (if any)
        Vector<LevelData<FArrayBox>*> amrVel(0), amrE(0), amrT(0), amrS(0),
            amrScalars(0);
        const int nScal = this->numScalars();
        {
            const AMRNSLevel* levPtr = this;
            while (levPtr) {
                // Allocate CC holders
                amrVel.push_back(new LevelData<FArrayBox>(
                    levPtr->getBoxes(), SpaceDim, IntVect::Zero));

                amrE.push_back(new LevelData<FArrayBox>(
                    levPtr->getBoxes(), 1, IntVect::Zero));

                amrT.push_back(new LevelData<FArrayBox>(
                    levPtr->getBoxes(), 1, IntVect::Zero));

                amrS.push_back(new LevelData<FArrayBox>(
                    levPtr->getBoxes(), 1, IntVect::Zero));

                if (nScal > 0) {
                    amrScalars.push_back(new LevelData<FArrayBox>);
                }

                // bdrySum[vel]
                LevelData<FluxBox> vel(levPtr->getBoxes(), 1);
                levPtr->sendToAdvectingVelocity(vel, levPtr->m_statePtr->vel);
                Real bdrySum = Integral::bdrySum(vel, *levPtr->m_levGeoPtr);
                levPtr->sendToCartesianVelocity(vel, vel);

                if (s_verbosity >= 2) {
                    pout() << "Level " << Format::fixed << levPtr->m_level
                           << " bdry momentum = " << Format::scientific
                           << bdrySum << endl;
                }

                // Send to CC vel.
                levPtr->m_levGeoPtr->multByJ(vel);
                Convert::FacesToCells(*amrVel.back(), vel);

                // Compute energy
                levPtr->computeTotalEnergy(*amrE.back(), m_time);
                // setValLevel(*amrE.back(), 0.0);
                // levPtr->addKE(*amrE.back(), vel);
                vel.clear();

                // Alias T, S, scalars
                aliasLevelData(
                    *amrT.back(), &levPtr->m_statePtr->T, Interval(0, 0));

                aliasLevelData(
                    *amrS.back(), &levPtr->m_statePtr->S, Interval(0, 0));

                if (nScal > 0) {
                    aliasLevelData(*amrScalars.back(),
                                   &levPtr->m_statePtr->scalars,
                                   Interval(0, nScal - 1));
                }

                // Move to next level
                levPtr = levPtr->fineNSPtr();
            }
        }

        // Integrate the momentum and energy
        RealVect totalMom;
        for (int dir = 0; dir < SpaceDim; ++dir) {
            // Remember, vel is alredy scaled by J.
            const bool scaleByJ = false;
            totalMom[dir] = Integral::sum(amrVel, *m_levGeoPtr, scaleByJ, dir);
        }

        // RealVect bdryMom(D_DECL(0., 0., 0.));
        // PhysBdryIter it(m_levGeoPtr->getBoxes());
        // for (it.reset(); it.ok(); ++it) {
        //     const DataIndex&     di         = it->di;
        //     const int            isign      = it->isign;
        //     const Real           rsign      = it->rsign;
        //     const int            bdryDir    = it->dir;
        //     const Box&           fcBdryBox  = it->fcBdryBox;

        //     const IntVect    eb   = ((isign > 0) ? IntVect::Zero : BASISV(bdryDir));
        //     const FArrayBox& vFAB = m_statePtr->vel[di][bdryDir];
        //     const RealVect&  dXi  = m_levGeoPtr->getDXi();
        //     const Real       dA   = dXi.product() / dXi[bdryDir];

        //     bdryMom[bdryDir] += rsign * vFAB.sum(fcBdryBox, 0, 1) * dA;
        // }
        // for (int dir = 0; dir < SpaceDim; ++dir) {
        //     Comm::reduce(bdryMom[dir], MPI_SUM);
        //     pout() << "Boundary momentum sum[" << dir << "] = " << bdryMom[dir] << endl;
        // }
        // pout() << "Net boundary momentum = " << bdryMom.sum() << endl;

        const Real totalE = Integral::sum(amrE, *m_levGeoPtr);
        static Real lastTotalE = 0.0;
        Real dEonE = (totalE - lastTotalE) / lastTotalE;
        lastTotalE = totalE;

        // Integrate T and S.
        {
            Real totalT = Integral::sum(amrT, *m_levGeoPtr, true, 0);
            Real totalS = Integral::sum(amrS, *m_levGeoPtr, true, 0);

            if (s_verbosity >= 2) {
                pout() << "Sum[T] = " << totalT << endl;
                pout() << "Sum[S] = " << totalS << endl;
            }
        }

        // Integrate the scalars.
        Vector<Real> totalScalars(nScal, -1.0);
        for (int comp = 0; comp < nScal; ++comp) {
            totalScalars[comp] =
                Integral::sum(amrScalars, *m_levGeoPtr, true, comp);

            if (s_verbosity >= 2) {
                pout() << "Sum[" << this->getScalarName(comp)
                    << "] = " << totalScalars[comp] << endl;
            }
        }

        // Free memory
        for (unsigned int idx = 0; idx < amrVel.size(); ++idx) {
            if (nScal > 0) {
                delete amrScalars[idx];
                amrScalars[idx] = nullptr;
            }

            delete amrS[idx];
            amrS[idx] = nullptr;

            delete amrT[idx];
            amrT[idx] = nullptr;

            delete amrE[idx];
            amrE[idx] = nullptr;

            delete amrVel[idx];
            amrVel[idx] = nullptr;
        }
        amrScalars.resize(0);
        amrE.resize(0);
        amrVel.resize(0);


        // Write!
        IO::tout(0) << Format::pushFlags;
        IO::tout(0) << setfill(' ') << std::right
                    << Format::fixed
                    << std::setw(stepWidth) << a_step
                    << "│";

        IO::tout(0) << setfill(' ') << std::right
                    << Format::fixed
                    << Format::hiblack
                    << std::setw(levelWidth) << m_level
                    << Format::white
                    << "│";

        if (a_step == 0 || (1.0e-5 < abs(m_time) && abs(m_time) < 1.0e5)) {
            IO::tout(0) << Format::fixed
                        << std::setprecision(5)
                        << std::setw(timeWidth) << m_time;
        } else {
            IO::tout(0) << Format::scientific
                        << std::setprecision(3)
                        << std::setw(timeWidth) << m_time;
        }

        if (a_step == 0 || (1.0e-5 < abs(m_dt) && abs(m_dt) < 1.0e5)) {
            IO::tout(0) << Format::fixed
                        << std::setprecision(5)
                        << std::setw(dtWidth) << m_dt;
        } else {
            IO::tout(0) << Format::scientific
                        << std::setprecision(3)
                        << std::setw(dtWidth) << m_dt;
        }

        IO::tout(0) << "│"
                    << Format::scientific
                    D_TERM(
                    << std::setprecision(5) << std::setw(momWidth) << totalMom[0],
                    << std::setprecision(5) << std::setw(momWidth) << totalMom[1],
                    << std::setprecision(5) << std::setw(momWidth) << totalMom[2])
                    << "│"
                    << std::setprecision(5) << std::setw(eWidth) << totalE;
        if (a_step > 0) {
            if (dEonE <= 0.0) {
                IO::tout(0) << std::setprecision(3)
                            << std::setw(elWidth) << dEonE;
            } else {
                IO::tout(0) << Format::hipurple
                            << std::setprecision(3)
                            << std::setw(elWidth) << dEonE
                            << Format::none;
            }
        }
        IO::tout(0) << endl << Format::popFlags;
    }
}


