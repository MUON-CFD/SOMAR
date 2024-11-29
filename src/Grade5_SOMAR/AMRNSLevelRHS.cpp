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
#include "Debug.H"
#include "ProblemContext.H"
#include "SetValLevel.H"
#include "Convert.H"
#include "FABAlgebra.H"
#include "MiscUtils.H"
#include "LayoutTools.H"
#include "NodeFArrayBox.H"
#include "EdgeDataBox.H"
#include "BoxIterator.H"
#include "BCTools.H"
#include "Subspace.H"
#include <chrono>
#include "GNUC_Extensions.H"

// Viscous solver stuff
#include "LevelSolver.H"
#include "ViscousOp.H"
#include "DiffusiveOp.H"
#include "VelBC.H"
#include "ScalarBC.H"


#ifndef NDEBUG
    // Debug mode
#   define nanCheck(x) checkForValidNAN(x)
#else
    // Release mode
#   define nanCheck(x)
#endif



// -----------------------------------------------------------------------------
// ARKRHS override.
// The explicit Navier-Stokes forcing. Be sure to set the comps with no
// forcing to zero. Q is non-const so that you can set BCs.
//
// Both a_vel and a_kvel are in the Cartesian basis. If you need
// a_vel in the mapped basis, use LevelGeometry::sendToMappedBasis on
// a local copy.
// -----------------------------------------------------------------------------
void
AMRNSLevel::setExplicitRHS(LevelData<FluxBox>&   a_kvel,
                           LevelData<FArrayBox>& a_kq,
                           LevelData<FluxBox>&   a_vel,
                           LevelData<FArrayBox>& a_p,
                           LevelData<FArrayBox>& a_q,
                           const Real            a_time,
                           const Real            a_refluxDt)
{
    // Collect references
    const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();
    DataIterator             dit   = grids.dataIterator();
    const ProblemContext*    ctx   = ProblemContext::getInstance();

    // Prepare state variables
    LevelData<FluxBox> cartVel;
    LevelData<FluxBox> advVel;
    LevelData<FArrayBox> T, S;
    LevelData<FArrayBox> Tpert, Spert, bpert;
    LevelData<FArrayBox> scalars;
    {
        nanCheck(a_vel);
        nanCheck(a_p);
        nanCheck(a_q);

        // cartVel
        aliasLevelData(cartVel, &a_vel, a_vel.interval());
        this->setVelBC(cartVel, a_time, false);

        // advVel
        advVel.define(grids, 1, IntVect::Unit);
        this->sendToAdvectingVelocity(advVel, cartVel);

        // T, Tpert
        aliasLevelData(T, &a_q, m_statePtr->TInterval);
        this->setTemperatureBC(T, a_time, false);
        Tpert.define(grids, 1, IntVect::Unit);
        for (dit.reset(); dit.ok(); ++dit) {
            Tpert[dit].copy(T[dit]);
        }
        Subspace::addHorizontalExtrusion(Tpert, 0, *m_TbarPtr, 0, 1, -1.0);

        // S, Spert
        aliasLevelData(S, &a_q, m_statePtr->SInterval);
        this->setSalinityBC(S, a_time, false);
        Spert.define(grids, 1, IntVect::Unit);
        for (dit.reset(); dit.ok(); ++dit) {
            Spert[dit].copy(S[dit]);
        }
        Subspace::addHorizontalExtrusion(Spert, 0, *m_SbarPtr, 0, 1, -1.0);

        // bpert
        bpert.define(grids, 1, IntVect::Unit);
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox zFAB(bpert[dit].box(), 1);
            m_levGeoPtr->fill_physCoor(zFAB, 0, SpaceDim - 1);
            this->equationOfState(bpert[dit], T[dit], S[dit], zFAB);
        }
        Subspace::addHorizontalExtrusion(bpert, 0, *m_bbarPtr, 0, 1, -1.0);

        // scalars
        if (m_statePtr->numScalars > 0) {
            aliasLevelData(scalars, &a_q, m_statePtr->scalarsInterval);
            this->setScalarBC(scalars, a_time, false);
        }

        nanCheck(cartVel);
        nanCheck(advVel);
        debugCheckValidFaceOverlap(cartVel);
        debugCheckValidFaceOverlap(advVel);
        nanCheck(a_p);
        nanCheck(T);
        nanCheck(Tpert);
        nanCheck(S);
        nanCheck(Spert);
        nanCheck(bpert);
    }

    // Create workspace
    StaggeredFluxLD      momentumFlux(grids);
    LevelData<FluxBox>   qFlux(grids, 1);
    LevelData<FArrayBox> qDiv(grids, 1);

    // Start clean
    setValLevel(a_kvel, 0.0);
    setValLevel(a_kq, 0.0);

    // Momentum advection
    if (ctx->rhs.doMomentumAdvection) {
        this->computeMomentumAdvection(a_kvel, momentumFlux, cartVel, advVel);

        // Refluxing
        if (ctx->rhs.doMomAdvRefluxing) {
            this->incrementFluxRegisters(momentumFlux, a_refluxDt);
        }

        nanCheck(a_kvel);
    }


    // Temperature advection
    if (ctx->rhs.doTemperatureAdvection) {
        LevelData<FArrayBox> kT;
        aliasLevelData(kT, &a_kq, m_statePtr->TInterval);

        // Linear eqs.
        // LevelData<FArrayBox> Tbar(grids, 1, IntVect::Unit);
        // setValLevel(Tbar, 0.0);
        // Subspace::addHorizontalExtrusion(Tbar, 0, *m_TbarPtr, 0, 1, 1.0);
        // this->computeScalarAdvection(kT, qFlux, Tbar, advVel);

        this->computeScalarAdvection(kT, qFlux, T, advVel);

        if (ctx->rhs.doTemperatureAdvRefluxing) {
            this->incrementFluxRegisters(
                qFlux, m_statePtr->TInterval, a_refluxDt);
        }
    } // end temperature advection


    // Salinity advection
    if (ctx->rhs.doSalinityAdvection) {
        LevelData<FArrayBox> kS;
        aliasLevelData(kS, &a_kq, m_statePtr->SInterval);

        this->computeScalarAdvection(kS, qFlux, S, advVel);

        if (ctx->rhs.doSalinityAdvRefluxing) {
            this->incrementFluxRegisters(
                qFlux, m_statePtr->SInterval, a_refluxDt);
        }
    } // end salinity advection


    // Scalar advection
    if (ctx->rhs.doScalarAdvection && this->numScalars() > 0) {
        const int startComp = m_statePtr->scalarsInterval.begin();
        const int endComp   = m_statePtr->scalarsInterval.end();

        for (int comp = startComp; comp <= endComp; ++comp) {
            const Interval ivl(comp, comp);

            LevelData<FArrayBox> kscalar, scalar;
            aliasLevelData(kscalar, &a_kq, ivl);
            aliasLevelData(scalar, &a_q, ivl);

            this->computeScalarAdvection(kscalar, qFlux, scalar, advVel);

            if (ctx->rhs.doScalarAdvRefluxing) {
                this->incrementFluxRegisters(qFlux, ivl, a_refluxDt);
            }
        }
    } // end scalars advection


    // Eddy viscosity and diffusivity
    LevelData<FArrayBox> eddyNu;
    aliasLevelData(eddyNu, &a_q, m_statePtr->eddyNuInterval);
    this->computeEddyNu(eddyNu, cartVel, a_time);

    // // Viscous forcing
    // if (ctx->rhs.doViscousForcing) {
    //     const RealVect nu = RealVect::Unit * ctx->rhs.nu;
    //     const Real primaryScale   = ctx->rhs.doImplicitDiffusion ? 0.0 : 1.0;
    //     const Real transposeScale = 1.0;

    //     this->computeMomentumDiffusion(a_kvel,
    //                                    momentumFlux,
    //                                    cartVel,
    //                                    nu,
    //                                    eddyNu,
    //                                    primaryScale,
    //                                    transposeScale);

    //     // Refluxing
    //     if (ctx->rhs.doViscousRefluxing) {
    //         this->incrementFluxRegisters(momentumFlux, a_refluxDt);
    //     }
    // } // end viscous forcing

    // Viscous forcing (w/ debugging)
    if (ctx->rhs.doViscousForcing) {
        const RealVect nu = RealVect::Unit * ctx->rhs.nu;
        const Real primaryScale   = ctx->rhs.doImplicitDiffusion ? 0.0 : 1.0;
        const Real transposeScale = 1.0;

        // nu
        LevelData<FluxBox> viscForce(grids, a_kvel.nComp(), a_kvel.ghostVect());
        setValLevel(viscForce, 0.0);

        LevelData<FArrayBox> zeroEddyNu(grids, 1, eddyNu.ghostVect());
        setValLevel(zeroEddyNu, 0.0);

        this->computeMomentumDiffusion(viscForce,
                                       momentumFlux,
                                       cartVel,
                                       nu,
                                       zeroEddyNu,
                                       primaryScale,
                                       transposeScale);

        auto maxViscForce = Analysis::pNorm(viscForce, 0);
        pout() << "max |visc force| = " << Format::scientific << maxViscForce << '\n';

        for (dit.reset(); dit.ok(); ++dit) {
            a_kvel[dit] += viscForce[dit];
        }


        // eddyVisc
        setValLevel(viscForce, 0.0);
        const RealVect zeroNu = RealVect::Zero;

        this->computeMomentumDiffusion(viscForce,
                                       momentumFlux,
                                       cartVel,
                                       zeroNu,
                                       eddyNu,
                                       primaryScale,
                                       transposeScale);

        maxViscForce = Analysis::pNorm(viscForce, 0);
        pout() << "max |eddy visc force| = " << Format::scientific << maxViscForce << '\n';

        for (dit.reset(); dit.ok(); ++dit) {
            a_kvel[dit] += viscForce[dit];
        }

        // Refluxing
        if (ctx->rhs.doViscousRefluxing) {
            this->incrementFluxRegisters(momentumFlux, a_refluxDt);
        }
    } // end viscous forcing


    // Temperature diffusion
    if (ctx->rhs.doTemperatureDiffusion && !ctx->rhs.doImplicitDiffusion) {
        const Real kappa        = ctx->rhs.TKappa;
        const Real eddyPrandtlT = ctx->rhs.eddyPrandtlT;

        LevelData<FArrayBox> kT;
        aliasLevelData(kT, &a_kq, m_statePtr->TInterval);

        this->computeScalarDiffusion(
            qDiv, qFlux, Tpert, kappa, eddyNu, eddyPrandtlT);

        for (dit.reset(); dit.ok(); ++dit) {
            kT[dit].plus(qDiv[dit], 1.0);
        }

        if (ctx->rhs.doTemperatureDiffusiveRefluxing) {
            this->incrementFluxRegisters(
                qFlux, m_statePtr->TInterval, a_refluxDt);
        }
    } // end temperature diffusion


    // Salinity diffusion
    if (ctx->rhs.doSalinityDiffusion && !ctx->rhs.doImplicitDiffusion) {
        const Real kappa = ctx->rhs.SKappa;
        const Real eddyPrandtlS = ctx->rhs.eddyPrandtlS;

        LevelData<FArrayBox> kS;
        aliasLevelData(kS, &a_kq, m_statePtr->SInterval);

        this->computeScalarDiffusion(
            qDiv, qFlux, Spert, kappa, eddyNu, eddyPrandtlS);

        for (dit.reset(); dit.ok(); ++dit) {
            kS[dit].plus(qDiv[dit], 1.0);
        }

        if (ctx->rhs.doSalinityDiffusiveRefluxing) {
            this->incrementFluxRegisters(
                qFlux, m_statePtr->SInterval, a_refluxDt);
        }
    } // end salinity diffusion


    // Scalar diffusion (overwrites ks)
    if (ctx->rhs.doScalarDiffusion && this->numScalars() > 0 && !ctx->rhs.doImplicitDiffusion) {
        const int startComp = m_statePtr->scalarsInterval.begin();
        const int endComp  = m_statePtr->scalarsInterval.end();

        for (int comp = startComp; comp <= endComp; ++comp) {
            const Interval ivl(comp, comp);

            const Real sKappa = ctx->rhs.getScalarsKappa(comp - startComp);
            const Real eddyPr = ctx->rhs.getEddyPrandtlScalars(comp - startComp);

            LevelData<FArrayBox> kscomp, scomp;
            aliasLevelData(kscomp, &a_kq, ivl);
            aliasLevelData(scomp, &a_q, ivl);

            this->computeScalarDiffusion(
                qDiv, qFlux, scomp, sKappa, eddyNu, eddyPr);

            for (dit.reset(); dit.ok(); ++dit) {
                kscomp[dit].plus(qDiv[dit], 1.0);
            }

            if (ctx->rhs.doScalarDiffusiveRefluxing) {
                this->incrementFluxRegisters(qFlux, ivl, a_refluxDt);
            }
        } // end loops over scalar comps (comp)
    } // end scalars diffusion


    // ------------------------------------------------------------------
    // All of the forces above this line need to be scaled by 1/J.
    for(dit.reset(); dit.ok(); ++dit) {
        // The CC fields are easy.
        const FArrayBox& ccJinvFAB = m_levGeoPtr->getCCJinv()[dit];
        const Box        ccValid   = grids[dit];
        for (int comp = 0; comp < a_kq.nComp(); ++comp) {
            a_kq[dit].mult(ccJinvFAB, ccValid, 0, comp, 1);
        }

        // The FC fields require a harmonic average of Jinv.
        // That is, first average J to FC, then reciprocate.
        const FArrayBox& ccJFAB = m_levGeoPtr->getCCJ()[dit];
        for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
            const Box fcValid = surroundingNodes(ccValid, fcDir);

            FArrayBox fcJFAB(fcValid, 1);
            Convert::Simple(fcJFAB, ccJFAB);
            a_kvel[dit][fcDir].divide(fcJFAB, fcValid, 0, 0, 1);
        }
    }

    // ------------------------------------------------------------------
    // All of the forces below this line DO NOT need to be scaled by 1/J.

    // Gravity forcing
    if (ctx->rhs.doGravityForcing) {
        // Original version...
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox&       kwFAB    = a_kvel[dit][SpaceDim - 1];
            const FArrayBox& bpertFAB = bpert[dit];
            const FArrayBox& ccJFAB   = m_levGeoPtr->getCCJ()[dit];
            const Box        fcValid  = grids[dit].surroundingNodes(SpaceDim - 1);

            FORT_ADDEXPLICITGRAVITYFORCING (
                CHF_FRA1(kwFAB, 0),
                CHF_CONST_FRA1(bpertFAB, 0),
                CHF_CONST_FRA1(ccJFAB, 0),
                CHF_BOX(fcValid));
        }

        // // 4th order version...
        // LevelData<FArrayBox> b2(grids, 1, 2*IntVect::Unit);
        // for (dit.reset(); dit.ok(); ++dit) {
        //     b2[dit].copy(bpert[dit]);
        // }
        // BCTools::extrapAllGhosts(b2, 2, IntVect::Unit);
        // b2.exchange();

        // for (dit.reset(); dit.ok(); ++dit) {
        //     FArrayBox&       kwFAB    = a_kvel[dit][SpaceDim - 1];
        //     const FArrayBox& bpertFAB = b2[dit];
        //     const FArrayBox& ccJFAB   = m_levGeoPtr->getCCJ()[dit];
        //     const Box        fcValid  = grids[dit].surroundingNodes(SpaceDim - 1);

        //     FArrayBox fcbFAB(fcValid, 1);
        //     Convert::fourthOrderFV_Unmapped(fcbFAB, Interval(0,0), fcValid, bpertFAB, Interval(0,0));
        //     kwFAB.plus(fcbFAB, fcValid, fcValid, -1.0, 0, 0, 1);
        // }

        nanCheck(a_kvel);
    }

    // Tidal forcing -- Cartesian only
    // Applying this removes the barotropic component from the "pressure."
    if (ctx->rhs.doTidalForcing) {
        static const RealVect& tidalU0        = ctx->rhs.tidalU0;
        static const RealVect& tidalOmega     = ctx->rhs.tidalOmega;
        static const RealVect& tidalInitPhase = ctx->rhs.tidalInitPhase;
        Real arg, tidalF;

        for (int dir = 0; dir < SpaceDim; ++dir) {
            arg = tidalOmega[dir] * a_time + tidalInitPhase[dir];
            tidalF = tidalU0[dir] * tidalOmega[dir] * cos(arg);

            if (RealCmp::isZero(tidalF)) continue;

            for (dit.reset(); dit.ok(); ++dit) {
                a_kvel[dit][dir].plus(tidalF);
            }
        }
    }

    // Coriolis forcing -- Cartesian only
    if (ctx->rhs.doCoriolisForcing &&
        ctx->rhs.coriolisF.vectorLength() > smallReal)
    {
#if CH_SPACEDIM == 2
        BUG("Coriolis forcing has not been tested in 2D.");
#endif
        const Real coriolisFy = ctx->rhs.coriolisF[SpaceDim - 2];
        const Real coriolisFz = ctx->rhs.coriolisF[SpaceDim - 1];

        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox&       kuFAB = a_kvel[dit][0];
            FArrayBox&       kvFAB = a_kvel[dit][1];
            FArrayBox&       kwFAB = a_kvel[dit][SpaceDim - 1];
            const FArrayBox& uFAB  = a_vel[dit][0];
            const FArrayBox& vFAB  = a_vel[dit][1];
            const FArrayBox& wFAB  = a_vel[dit][SpaceDim - 1];
            const Box&       kuBox = surroundingNodes(grids[dit], 0);
            const Box&       kvBox = surroundingNodes(grids[dit], 1);
            const Box&       kwBox = surroundingNodes(grids[dit], SpaceDim - 1);

            FORT_ADDEXPLICITCORIOLISFORCING3D(
                CHF_FRA1(kuFAB, 0),
                CHF_FRA1(kvFAB, 0),
                CHF_FRA1(kwFAB, 0),
                CHF_CONST_FRA1(uFAB, 0),
                CHF_CONST_FRA1(vFAB, 0),
                CHF_CONST_FRA1(wFAB, 0),
                CHF_CONST_REAL(coriolisFy),
                CHF_CONST_REAL(coriolisFz),
                CHF_BOX(kuBox),
                CHF_BOX(kvBox),
                CHF_BOX(kwBox));
        }
    }

    // debugCheckValidFaceOverlap(a_kvel);
    nanCheck(a_kvel);
    nanCheck(a_kq);

    // Let the user modify the NS forces.
    this->addExplicitRHS(
        a_kvel, a_kq, a_vel, a_p, a_q, advVel, a_time, a_refluxDt);

    // Small discrepancies at overlying faces on different grids can lead to
    // large discrepancies after many timesteps.If left unchecked, they magnify
    // and take down some simulations. Averaging the two overlying values seems
    // to help, but we need to find the source of the discrepancies.
    // BUG("Figure out why we have discrepancies at overlapping faces of a_kvel.");
    LayoutTools::averageOverlappingValidFaces(a_kvel);

    // debugCheckValidFaceOverlap(a_kvel);
    nanCheck(a_kvel);
    nanCheck(a_kq);
}


// -----------------------------------------------------------------------------
// ARKRHS override.
// The implicit Navier-Stokes forcing. Be sure to set the comps with no
// forcing to zero. Q is non-const so that you can set BCs.
//
// Both a_vel and a_kvel are in the Cartesian basis. If you need
// a_vel in the mapped basis, use LevelGeometry::sendToMappedBasis on
// a local copy.
// -----------------------------------------------------------------------------
void
AMRNSLevel::setImplicitRHS(LevelData<FluxBox>&   a_kvel,
                           LevelData<FArrayBox>& a_kq,
                           LevelData<FluxBox>&   a_vel,
                           LevelData<FArrayBox>& a_p,
                           LevelData<FArrayBox>& a_q,
                           const Real            a_time,
                           const Real            a_refluxDt)
{
    setValLevel(a_kvel, 0.0);
    setValLevel(a_kq, 0.0);

    const ProblemContext* ctx = ProblemContext::getInstance();
    if (!ctx->rhs.doImplicitDiffusion) return;

    const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();

    LevelData<FluxBox> fcFluxes(grids, 1); // workspace

    LevelData<FArrayBox> eddyNu;
    aliasLevelData(eddyNu, &a_q, m_statePtr->eddyNuInterval);
    this->computeEddyNu(eddyNu, a_vel, a_time);

    // Viscous forcing
    if (ctx->rhs.doViscousForcing) {
        this->setVelBC(a_vel, a_time, false);

        StaggeredFluxLD momentumFlux(grids);
        const RealVect nu = RealVect::Unit * ctx->rhs.nu;
        const Real primaryScale   = 1.0;
        const Real transposeScale = 0.0; // was done explicitly

        this->computeMomentumDiffusion(a_kvel,
                                       momentumFlux,
                                       a_vel,
                                       nu,
                                       eddyNu,
                                       primaryScale,
                                       transposeScale);
    }

    // Temperature diffusion
    if (ctx->rhs.doTemperatureDiffusion) {
        LevelData<FArrayBox> T, kT;
        aliasLevelData(T, &a_q, m_statePtr->TInterval);
        aliasLevelData(kT, &a_kq, m_statePtr->TInterval);

        // Set BCs on T
        this->setTemperatureBC(T, a_time, false);

        // Tpert
        LevelData<FArrayBox> Tpert(grids, 1, IntVect::Unit);
        for (DataIterator dit(grids); dit.ok(); ++dit) {
            Tpert[dit].copy(T[dit]);
        }
        Subspace::addHorizontalExtrusion(Tpert, 0, *m_TbarPtr, 0, 1, -1.0);

        // Compute kT
        this->computeScalarDiffusion(kT,
                                     fcFluxes,
                                     Tpert,
                                     ctx->rhs.TKappa,
                                     eddyNu,
                                     ctx->rhs.eddyPrandtlT);
    }

    // Salinity diffusion
    if (ctx->rhs.doSalinityDiffusion) {
        LevelData<FArrayBox> S, kS;
        aliasLevelData(S, &a_q, m_statePtr->SInterval);
        aliasLevelData(kS, &a_kq, m_statePtr->SInterval);

        // Set BCs on S
        this->setSalinityBC(S, a_time, false);

        // Spert
        LevelData<FArrayBox> Spert(grids, 1, IntVect::Unit);
        for (DataIterator dit(grids); dit.ok(); ++dit) {
            Spert[dit].copy(S[dit]);
        }
        Subspace::addHorizontalExtrusion(Spert, 0, *m_SbarPtr, 0, 1, -1.0);

        // Compute kS
        this->computeScalarDiffusion(kS,
                                     fcFluxes,
                                     Spert,
                                     ctx->rhs.SKappa,
                                     eddyNu,
                                     ctx->rhs.eddyPrandtlS);
    }

    // Scalars diffusion
    if (ctx->rhs.doScalarDiffusion && (this->numScalars() > 0)) {
        LevelData<FArrayBox> scalars, kScalars;
        aliasLevelData(scalars, &a_q, m_statePtr->scalarsInterval);
        aliasLevelData(kScalars, &a_kq, m_statePtr->scalarsInterval);

        // Set BCs on S
        this->setScalarBC(scalars, a_time, false);

        for (int comp = 0; comp < this->numScalars(); ++comp) {
            LevelData<FArrayBox> S, kS;
            aliasLevelData(S, &scalars, Interval(comp, comp));
            aliasLevelData(kS, &kScalars, Interval(comp, comp));

            // Compute kS
            this->computeScalarDiffusion(kS,
                                         fcFluxes,
                                         S,
                                         ctx->rhs.getScalarsKappa(comp),
                                         eddyNu,
                                         ctx->rhs.getEddyPrandtlScalars(comp));
        }
    }

    // All forces need to be scaled by 1/J.
    m_levGeoPtr->divByJ(a_kvel);
    m_levGeoPtr->divByJ(a_kq);
}


// -----------------------------------------------------------------------------
// ARKRHS override.
// Solves [1 - gammaDt * D] Q^{n+1} = Q^{n} in place.
//
// a_vel is in the Cartesian basis.
// -----------------------------------------------------------------------------
void
AMRNSLevel::solveImplicit(LevelData<FluxBox>&   a_vel,
                          LevelData<FArrayBox>& a_q,
                          const Real            a_gammaDt,
                          const Real            a_time,
                          const Real            /*a_refluxDt*/)
{
    const ProblemContext* ctx = ProblemContext::getInstance();
    if (!ctx->rhs.doImplicitDiffusion) return;

    const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();

    LevelData<FArrayBox> eddyNu;
    aliasLevelData(eddyNu, &a_q, m_statePtr->eddyNuInterval);
    this->computeEddyNu(eddyNu, a_vel, a_time);

    // Compute dz values at z-faces.
    const std::shared_ptr<FArrayBox> dzFABPtr = [&]() {
        const Box dzBox = Subspace::verticalDataBox(m_levGeoPtr->getDomainBox())
                              .surroundingNodes(SpaceDim - 1);

        std::shared_ptr<FArrayBox> ptr(new FArrayBox(dzBox, 1));

        for (IntVect iv = dzBox.smallEnd();
             iv[SpaceDim - 1] <= dzBox.bigEnd(SpaceDim - 1);
             ++iv[SpaceDim - 1]) {
            (*ptr)(iv) = m_levGeoPtr->getFaceDx(SpaceDim - 1, iv, SpaceDim - 1);
        }

        return ptr;
    }();

    // Viscous solve
    if (ctx->rhs.doViscousForcing) {
        for (DataIterator dit(grids); dit.ok(); ++dit) {
            eddyNu[dit] += ctx->rhs.nu;
        }

        this->solveMomentumDiffusion(a_vel, eddyNu, a_gammaDt, a_time);

        for (DataIterator dit(grids); dit.ok(); ++dit) {
            eddyNu[dit] -= ctx->rhs.nu;
        }
    }

    // Temperature diffusion
    if (ctx->rhs.doTemperatureDiffusion) {
        LevelData<FArrayBox> T;
        aliasLevelData(T, &a_q, m_statePtr->TInterval);

        // crseTPtr
        LevelData<FArrayBox> crseT;
        if (m_level > 0) {
            const auto* crsePtr          = this->crseNSPtr();
            const auto& crseGrids        = crsePtr->getBoxes();
            const bool  setBCs           = false;
            const bool  removeBackground = true;

            crseT.define(crseGrids, 1);
            crsePtr->fillTemperature(crseT, a_time, setBCs, removeBackground);
        }
        const LevelData<FArrayBox>* crseTPtr = m_level ? &crseT : nullptr;

        CH_assert(m_TbarPtr);
        const std::shared_ptr<BCTools::BCFunction> bcFuncPtr(
            new ScalarBC::BackgroundScalarWrapper(this->temperaturePhysBC(),
                                                  m_TbarPtr,
                                                  dzFABPtr));

        this->solveScalarDiffusion(T,
                                   crseTPtr,
                                   {ctx->rhs.TKappa},
                                   eddyNu,
                                   {ctx->rhs.eddyPrandtlT},
                                   a_gammaDt,
                                   a_time,
                                   bcFuncPtr,
                                   "Temperature");
    }

    // Salinity diffusion
    if (ctx->rhs.doSalinityDiffusion) {
        LevelData<FArrayBox> S;
        aliasLevelData(S, &a_q, m_statePtr->SInterval);

        // crseSPtr
        LevelData<FArrayBox> crseS;
        if (m_level > 0) {
            const auto* crsePtr          = this->crseNSPtr();
            const auto& crseGrids        = crsePtr->getBoxes();
            const bool  setBCs           = false;
            const bool  removeBackground = true;

            crseS.define(crseGrids, 1);
            crsePtr->fillSalinity(crseS, a_time, setBCs, removeBackground);
        }
        const LevelData<FArrayBox>* crseSPtr = m_level ? &crseS : nullptr;

        CH_assert(m_SbarPtr);
        const std::shared_ptr<BCTools::BCFunction> bcFuncPtr(
            new ScalarBC::BackgroundScalarWrapper(this->salinityPhysBC(),
                                                  m_SbarPtr,
                                                  dzFABPtr));

        this->solveScalarDiffusion(S,
                                   crseSPtr,
                                   {ctx->rhs.SKappa},
                                   eddyNu,
                                   {ctx->rhs.eddyPrandtlS},
                                   a_gammaDt,
                                   a_time,
                                   bcFuncPtr,
                                   "Salinity");
    }

    // Scalars diffusion
    if (ctx->rhs.doScalarDiffusion && (this->numScalars() > 0)) {
        LevelData<FArrayBox> S;
        aliasLevelData(S, &a_q, m_statePtr->scalarsInterval);

        // crseSPtr
        LevelData<FArrayBox> crseS;
        if (m_level > 0) {
            const auto* crsePtr          = this->crseNSPtr();
            const auto& crseGrids        = crsePtr->getBoxes();
            const bool  setBCs           = false;

            crseS.define(crseGrids, this->numScalars());
            crsePtr->fillScalar(crseS, a_time, setBCs);
        }
        const LevelData<FArrayBox>* crseSPtr = m_level ? &crseS : nullptr;

        const std::vector<Real> vkappa = [&]() {
            std::vector<Real> v(S.nComp(), 0.0);
            for (int comp = 0; comp < S.nComp(); ++comp)
                v[comp] = ctx->rhs.getScalarsKappa(comp);
            return v;
        }();

        const std::vector<Real> vEddyPrandtl = [&]() {
            std::vector<Real> v(S.nComp(), 0.0);
            for (int comp = 0; comp < S.nComp(); ++comp)
                v[comp] = ctx->rhs.getEddyPrandtlScalars(comp);
            return v;
        }();

        this->solveScalarDiffusion(S,
                                   crseSPtr,
                                   vkappa,
                                   eddyNu,
                                   vEddyPrandtl,
                                   a_gammaDt,
                                   a_time,
                                   this->scalarsPhysBC(),
                                   "Custom scalars");
    }
}


// -----------------------------------------------------------------------------
// Computes viscous force and fluxes for the FC momentum (aka velocity).
// a_momentumFlux will be overwritten. a_kvel will be added to.
// Flux registers will not be updated. That is left to the caller.
// Result will not be scaled by 1/J. That is left to the caller.
// -----------------------------------------------------------------------------
void
AMRNSLevel::computeMomentumDiffusion(LevelData<FluxBox>&         a_kvel,
                                     StaggeredFluxLD&            a_momentumFlux,
                                     const LevelData<FluxBox>&   a_cartVel,
                                     const RealVect&             a_nu,
                                     const LevelData<FArrayBox>& a_eddyNu,
                                     const Real                  a_primaryScale,
                                     const Real                  a_transposeScale) const
{
    CH_assert(a_kvel.nComp() == 1);
    CH_assert(a_cartVel.nComp() == 1);
    CH_assert(a_eddyNu.nComp() == 1);

    CH_assert(a_kvel        .getBoxes() == m_levGeoPtr->getBoxes());
    CH_assert(a_momentumFlux.getBoxes() == m_levGeoPtr->getBoxes());
    CH_assert(a_cartVel     .getBoxes() == m_levGeoPtr->getBoxes());

    const RealVect&          dXi   = m_levGeoPtr->getDXi();
    const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();
    DataIterator             dit   = grids.dataIterator();

    // Compute J*grad[u] components. Last index will be Cartesian-based.
    m_finiteDiffPtr->levelVectorGradient(a_momentumFlux, a_cartVel);

    for (dit.reset(); dit.ok(); ++dit) {
        // In what follows, we are computing du^a/dt.
        // a is the Cartesian velocity component that we are updating.

        for (int a = 0; a < SpaceDim; ++a) {
            // 1. Convert J*grad[u] into 2*J*S^{ia}.

            // Easy, diagonal elements first.
            {
                FArrayBox& SaaFAB = a_momentumFlux[a][a][dit];
                SaaFAB *= (a_primaryScale + a_transposeScale);
                checkForNAN(SaaFAB, SaaFAB.box());
            }

            // Off-diagonal comps need to be symmetrized.
            for (int i = a + 1; i < SpaceDim; ++i) {
                FArrayBox& SiaFAB = a_momentumFlux[i][a][dit];
                FArrayBox& SaiFAB = a_momentumFlux[a][i][dit];
                const Box& ecBox  = SiaFAB.box();

                // Send last index to mapped basis.
                FArrayBox dxdXiFAB(ecBox, 2);
                m_levGeoPtr->getGeoSource().fill_dxdXi(dxdXiFAB, 0, a, dXi);
                m_levGeoPtr->getGeoSource().fill_dxdXi(dxdXiFAB, 1, i, dXi);
                SiaFAB.divide(dxdXiFAB, 0, 0, 1);
                SaiFAB.divide(dxdXiFAB, 1, 0, 1);

                // Symmetrize
                FArrayBox tmpFAB(ecBox, 1);
                tmpFAB.copy(SiaFAB);

                SiaFAB *= a_primaryScale;
                SiaFAB.plus(SaiFAB, a_transposeScale);

                SaiFAB *= a_primaryScale;
                SaiFAB.plus(tmpFAB, a_transposeScale);

                // Send last index back to Cartesian basis.
                SiaFAB.mult(dxdXiFAB, 0, 0, 1);
                SaiFAB.mult(dxdXiFAB, 1, 0, 1);

                checkForNAN(SiaFAB, ecBox);
                checkForNAN(SaiFAB, ecBox);
            } // i


            // 2. Convert 2*S^{ia} into the viscous stresses for each i,
            //    T^{ia} = (nu + eddyNu) * 2*J*S^{ia}.
            for (int i = 0; i < SpaceDim; ++i) {
                FArrayBox&       TiaFAB    = a_momentumFlux[i][a][dit];
                const Real       nuDir     = a_nu[i];
                const FArrayBox& eddyNuFAB = a_eddyNu[dit];
                const Box&       region    = TiaFAB.box();

                if (a == i) {
                    FABAlgebra::CCmultCC(
                        TiaFAB, 0, region, eddyNuFAB, 0, nuDir);
                } else {
                    FABAlgebra::ECmultCC(
                        TiaFAB, 0, region, eddyNuFAB, 0, nuDir);
                }
                checkForNAN(TiaFAB, region);
            }


            // 3. Compute the viscous force and add it to a_kvel.
            FArrayBox& FaFAB = a_kvel[dit][a];
            const Box  fcBox = surroundingNodes(grids[dit], a);
            const Real scale = 1.0;

            D_TERM(
            const int i0 = a;,
            const int i1 = (a + 1) % SpaceDim;,
            const int i2 = (a + 2) % SpaceDim;)

            D_TERM(
            const FArrayBox& T0aFAB = a_momentumFlux[i0][a][dit];,
            const FArrayBox& T1aFAB = a_momentumFlux[i1][a][dit];,
            const FArrayBox& T2aFAB = a_momentumFlux[i2][a][dit];)

            checkForNAN(FaFAB, fcBox);
            D_TERM(
            checkForNAN(T0aFAB, T0aFAB.box());,
            checkForNAN(T1aFAB, T1aFAB.box());,
            checkForNAN(T2aFAB, T2aFAB.box());)

#if CH_SPACEDIM == 2
            FORT_ADDVISCOUSFORCECOMP2D(
                CHF_FRA1(FaFAB, 0),
                CHF_CONST_FRA1(T0aFAB, 0),
                CHF_CONST_FRA1(T1aFAB, 0),
                CHF_BOX(fcBox),
                CHF_CONST_REALVECT(dXi),
                CHF_CONST_INT(i0),
                CHF_CONST_INT(i1),
                CHF_CONST_REAL(scale));
#else
            FORT_ADDVISCOUSFORCECOMP3D(
                CHF_FRA1(FaFAB, 0),
                CHF_CONST_FRA1(T0aFAB, 0),
                CHF_CONST_FRA1(T1aFAB, 0),
                CHF_CONST_FRA1(T2aFAB, 0),
                CHF_BOX(fcBox),
                CHF_CONST_REALVECT(dXi),
                CHF_CONST_INT(i0),
                CHF_CONST_INT(i1),
                CHF_CONST_INT(i2),
                CHF_CONST_REAL(scale));
#endif

            checkForNAN(FaFAB, fcBox);
        } // a
    } // dit

    nanCheck(a_kvel);
    LayoutTools::averageOverlappingValidFaces(a_kvel);
}


// -----------------------------------------------------------------------------
// Computes diffusive force and fluxes for any CC scalar.
// a_kq and a_qFlux will be overwritten.
// Flux registers will not be updated. That is left to the caller.
// Result will not be scaled by 1/J. That is left to the caller.
// -----------------------------------------------------------------------------
void
AMRNSLevel::computeScalarDiffusion(
    LevelData<FArrayBox>&       a_kq,
    LevelData<FluxBox>&         a_qFlux,
    const LevelData<FArrayBox>& a_q,
    const Real                  a_kappa,
    const LevelData<FArrayBox>& a_eddyNu,
    const Real                  a_eddyPrandtl) const
{
    CH_assert(a_kq.nComp() == 1);
    CH_assert(a_qFlux.nComp() == 1);
    CH_assert(a_q.nComp() == 1);

    const DisjointBoxLayout& grids = a_q.getBoxes();
    DataIterator dit = a_q.dataIterator();

    debugInitLevel(a_kq);
    debugInitLevel(a_qFlux);

    m_finiteDiffPtr->levelGradientMAC(a_qFlux, a_q);

    for (dit.reset(); dit.ok(); ++dit) {
        for (int dir = 0; dir < SpaceDim; ++dir) {
            const Box fcValid = surroundingNodes(grids[dit], dir);
            FABAlgebra::FCmultCC(a_qFlux[dit][dir],
                                 0,
                                 fcValid,
                                 a_eddyNu[dit],
                                 0,
                                 a_kappa,
                                 1.0 / a_eddyPrandtl);
        }
    }
    nanCheck(a_qFlux);

    const bool scaleByJinv = false;
    m_finiteDiffPtr->levelDivergenceMAC(a_kq, a_qFlux, scaleByJinv);
    nanCheck(a_kq);
}


// -----------------------------------------------------------------------------
void
AMRNSLevel::solveMomentumDiffusion(LevelData<FluxBox>&         a_cartVel,
                                   const LevelData<FArrayBox>& a_totalNu,
                                   const Real                  a_gammaDt,
                                   const Real                  a_time) const
{
    const ProblemContext* ctx = ProblemContext::getInstance();
    CH_assert(ctx->rhs.doImplicitDiffusion);
    CH_assert(ctx->rhs.doViscousForcing);

    const DisjointBoxLayout& grids       = m_levGeoPtr->getBoxes();
    const bool               writeToPout = (s_verbosity >= 1);
    const int                normType    = ctx->proj.normType;

    // Begin indentation...
    if (writeToPout) {
        pout() << "Viscous solve:"
               << Format::pushFlags
               << Format::indent()
               << '\n';
    }

    // crseVelPtr & set cartVel BCs
    LevelData<FluxBox> crseVel;
    if (m_level > 0) {
        const auto* crsePtr   = this->crseNSPtr();
        const auto& crseGrids = crsePtr->getBoxes();
        const bool  setBCs    = true;

        crseVel.define(crseGrids, 1, IntVect::Unit);
        crsePtr->fillVelocity(crseVel, a_time, setBCs);
    }
    LevelData<FluxBox>* crseVelPtr = m_level ? &crseVel : nullptr;
    this->setVelBC(a_cartVel, a_time, false, crseVelPtr);

    // Create op
    std::shared_ptr<Elliptic::ViscousOp> opPtr(new Elliptic::ViscousOp(
        m_levGeoPtr, 1.0, -a_gammaDt, a_totalNu, this->velPhysBC()));

    // Create solver
    // Elliptic::MGSolver<LevelData<FluxBox>> solver;
    // {
    //     auto opts = solver.getOptions();
    //     // opts.verbosity = 10;
    //     opts.bottomOptions.verbosity = 0;
    //     opts.bottomOptions.normType = normType;
    //     // opts.absTol = 1.0e-6;
    //     opts.relTol = 1.0e-6;
    //     opts.numSmoothDown = 0;
    //     opts.numSmoothUp = 0;
    //     opts.numSmoothPrecond = 0;
    //     opts.numSmoothBottom = 0;
    //     // opts.numSmoothUpFMG = 16;
    //     opts.maxDepth = 0;
    //     opts.numCycles = 1;
    //     opts.normType = normType;
    //     solver.define(*opPtr, opts);
    // }

    Elliptic::BiCGStabSolver<LevelData<FluxBox>> solver;
    solver.define(opPtr);
    solver.options().verbosity = 2;
    solver.options().normType = normType;

    // Elliptic::RelaxSolver<LevelData<FluxBox>> solver;
    // solver.define(opPtr);
    // solver.options().normType = normType;

    // Set up RHS
    LevelData<FluxBox> rhs0(grids, 1, IntVect::Unit);
    opPtr->assignLocal(rhs0, a_cartVel);
    m_levGeoPtr->multByJ(rhs0);

    // Set up initial guess.
    {
        LevelData<FluxBox> kvel(grids, 1, IntVect::Unit);
        StaggeredFluxLD momentumFlux(grids);
        opPtr->setToZero(kvel);
        this->computeMomentumDiffusion(
            kvel, momentumFlux, a_cartVel, RealVect::Zero, a_totalNu, 1.0, 0.0);
        opPtr->incr(a_cartVel, kvel, a_gammaDt);
    }

    // Solve!
    {
        const auto start = std::chrono::high_resolution_clock::now();
        constexpr bool homogBCs  = false;
        constexpr bool setToZero = false;
        solver.solve(a_cartVel, crseVelPtr, rhs0, a_time, homogBCs, setToZero);
        const auto finish = std::chrono::high_resolution_clock::now();

        // Write solve time.
        if (writeToPout) {
            std::chrono::duration<double> elapsed = finish - start;
            pout() << "Solve time: " << elapsed.count() << " s\n";
        }
    }

    // Restore indentation state.
    if (writeToPout) {
        pout() << Format::unindent
               << Format::popFlags
               << std::flush;
    }
}


// -----------------------------------------------------------------------------
void
AMRNSLevel::solveScalarDiffusion(
    LevelData<FArrayBox>&                       a_phi,
    const LevelData<FArrayBox>*                 a_crsePhiPtr,
    const std::vector<Real>&                    a_vKappa,
    const LevelData<FArrayBox>&                 a_eddyNu,
    const std::vector<Real>&                    a_vEddyPrandtl,
    const Real                                  a_gammaDt,
    const Real                                  a_time,
    const std::shared_ptr<BCTools::BCFunction>& a_bcFuncPtr,
    const std::string                           a_scalarName) const
{
    const ProblemContext* ctx = ProblemContext::getInstance();
    CH_assert(ctx->rhs.doImplicitDiffusion);
    CH_assert(ctx->rhs.doTemperatureDiffusion);

    const DisjointBoxLayout& grids       = m_levGeoPtr->getBoxes();
    const bool               writeToPout = (s_verbosity >= 1);
    const int                normType    = ctx->proj.normType;

    // Begin indentation...
    if (writeToPout) {
        pout() << Format::pushFlags;
        pout() << a_scalarName << " diffusive solve:"
                << Format::indent() << std::endl;
    }

    // Op = D[JgupKappa * G[T]], not scaled by 1/J!
    std::shared_ptr<Elliptic::DiffusiveOp> opPtr(
        new Elliptic::DiffusiveOp(*m_levGeoPtr,
                                  this->getCrseGridsPtr(),
                                  1.0,        // alpha
                                  -a_gammaDt, // beta
                                  a_vKappa,
                                  a_eddyNu,
                                  a_vEddyPrandtl,
                                  a_bcFuncPtr,
                                  a_phi.nComp()));

    // Solver
    // Elliptic::MGSolver<LevelData<FArrayBox>> solver;
    // {
    //     auto opts = solver.getOptions();
    //     // opts.verbosity = 10;
    //     opts.bottomOptions.verbosity = 0;
    //     opts.bottomOptions.normType = normType;
    //     opts.absTol = 1.0e-6;
    //     opts.relTol = 1.0e-6;
    //     // opts.numSmoothDown = 1;
    //     // opts.numSmoothUp = 6;
    //     // opts.numSmoothUpFMG = 16;
    //     // opts.maxDepth = 0;
    //     // opts.numCycles = -1;
    //     opts.normType = normType;
    //     solver.define(*opPtr, opts);
    // }

    Elliptic::BiCGStabSolver<LevelData<FArrayBox>> solver;
    solver.define(opPtr);
    solver.options().verbosity = 2;
    solver.options().normType = normType;

    // Elliptic::RelaxSolver<LevelData<FArrayBox>> solver;
    // solver.define(opPtr);
    // solver.options().normType = normType;

    // Now, we set up the eq...
    // rhs = J * phi
    LevelData<FArrayBox> rhs(grids, a_phi.nComp(), IntVect::Unit);
    opPtr->assignLocal(rhs, a_phi);
    m_levGeoPtr->multByJ(rhs);

    // Set phi BCs. RK doesn't seem to need this, but SDC does.
    opPtr->applyBCs(a_phi, a_crsePhiPtr, a_time, false, false);

    // Set up initial guess.
    {
        CH_verify(a_eddyNu.nComp() == 1);
        CH_verify((size_t)a_phi.nComp() == a_vKappa.size());
        CH_verify((size_t)a_phi.nComp() == a_vEddyPrandtl.size());

        LevelData<FArrayBox> kq(grids, 1, IntVect::Unit);
        LevelData<FluxBox> qFlux(grids, 1);

        for (int comp = 0; comp < a_phi.nComp(); ++comp) {
            const Interval ivl(comp, comp);

            LevelData<FArrayBox> q;
            aliasLevelData(q, &a_phi, ivl);

            opPtr->setToZero(kq);
            this->computeScalarDiffusion(
                kq, qFlux, q, a_vKappa[comp], a_eddyNu, a_vEddyPrandtl[comp]);
            opPtr->incr(q, kq, a_gammaDt);
        }
    }

    // Solve Op[phi] = rhs w/ inhomog BCs.
    {
        constexpr bool homogBCs  = false;
        constexpr bool setToZero = false;

        const auto start = std::chrono::high_resolution_clock::now();
        solver.solve(a_phi, a_crsePhiPtr, rhs, a_time, homogBCs, setToZero);
        const auto finish = std::chrono::high_resolution_clock::now();

        // Write solve time.
        if (writeToPout) {
            std::chrono::duration<double> elapsed = finish - start;
            pout() << "Solve time: " << elapsed.count() << " s\n";
        }

        // opPtr->applyBCs(a_phi, a_crsePhiPtr, a_time, homogBCs, homogBCs);
    }

    // Restore indentation state.
    if (writeToPout) {
        pout() << Format::unindent << std::flush;
        pout() << Format::popFlags;
    }
}


// -----------------------------------------------------------------------------
// Increments both the registers on this and the coarser level as needed.
// This function can be used to reflux the FC momentum, but not the
// CC scalars.
// a_flux must have the form a_flux[derivDir][velComp] where
//   d(u^{velComp})/dt = Div_{derivDir}(a_flux[derivDir][velComp]).
// -----------------------------------------------------------------------------
void
AMRNSLevel::incrementFluxRegisters(const StaggeredFluxLD& a_flux,
                                   const Real             a_refluxDt)
{
    if (RealCmp::isZero(a_refluxDt)) return;

    CH_assert(a_flux.getBoxes() == m_levGeoPtr->getBoxes());
    DataIterator dit(a_flux.getBoxes());

    // Increment flux register between this and the finer level.
    if (m_velFluxRegPtr) {
        const RealVect& crseDXi = m_levGeoPtr->getDXi();
        for (dit.reset(); dit.ok(); ++dit) {
            for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                for (int derivDir = 0; derivDir < SpaceDim; ++derivDir) {
                    const Real scale = a_refluxDt / crseDXi[derivDir];
                    for (SideIterator sit; sit.ok(); ++sit) {
                        m_velFluxRegPtr->incrementCoarse(
                            a_flux[derivDir][velComp][dit],
                            scale,
                            dit(),
                            Interval(0,0),
                            Interval(0,0),
                            velComp,
                            derivDir,
                            sit());
                    } // sit
                } // derivDir
            } // velComp
        } // dit
    }

    // Increment flux register between this and the coarser level.
    if (m_level > 0) {
        FluxRegisterFace* crseFRPtr = this->crseNSPtr()->m_velFluxRegPtr;
        const RealVect& crseDXi = this->crseNSPtr()->m_levGeoPtr->getDXi();

        for (dit.reset(); dit.ok(); ++dit) {
            for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                for (int derivDir = 0; derivDir < SpaceDim; ++derivDir) {
                    const Real scale = a_refluxDt / crseDXi[derivDir];
                    for (SideIterator sit; sit.ok(); ++sit) {
                        crseFRPtr->incrementFine(
                            a_flux[derivDir][velComp][dit],
                            scale,
                            dit(),
                            Interval(0,0),
                            Interval(0,0),
                            velComp,
                            derivDir,
                            sit());
                    } // sit
                } // derivDir
            } // velComp
        } // dit
    }
}


// -----------------------------------------------------------------------------
// Increments both the registers on this and the coarser level as needed.
// This function can be used to reflux any of the CC scalars, but not the
// FC momentum.
// -----------------------------------------------------------------------------
void
AMRNSLevel::incrementFluxRegisters(const LevelData<FluxBox>& a_flux,
                                   const Interval&           a_regInterval,
                                   const Real                a_refluxDt)

{
    if (RealCmp::isZero(a_refluxDt)) return;

    CH_assert(a_flux.getBoxes() == m_levGeoPtr->getBoxes());
    DataIterator dit = a_flux.dataIterator();

    // Increment flux register between this and the finer level.
    if (m_qFluxRegPtr) {
        const RealVect& crseDXi = m_levGeoPtr->getDXi();
        for (dit.reset(); dit.ok(); ++dit) {
            for (int dir = 0; dir < SpaceDim; ++dir) {
                const Real scale = a_refluxDt / crseDXi[dir];
                m_qFluxRegPtr->incrementCoarse(
                    a_flux[dit][dir],
                    scale,
                    dit(),
                    Interval(0,0),
                    a_regInterval,
                    dir);
            }
        }
    }

    // Increment flux register between this and the coarser level.
    if (m_level > 0) {
        AnisotropicFluxRegister* crseFRPtr = this->crseNSPtr()->m_qFluxRegPtr;
        const RealVect& crseDXi = this->crseNSPtr()->m_levGeoPtr->getDXi();
        for (dit.reset(); dit.ok(); ++dit) {
            for (int dir = 0; dir < SpaceDim; ++dir) {
                const Real scale = a_refluxDt / crseDXi[dir];
                crseFRPtr->incrementFine(
                    a_flux[dit][dir],
                    scale,
                    dit(),
                    Interval(0,0),
                    a_regInterval,
                    dir);
            }
        }
    }
}
