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
#include "LoadBalance.H"
#include "SetValLevel.H"
#include "Subspace.H"
#include "Debug.H"
#include <chrono>


// -----------------------------------------------------------------------------
// Tag cells during runtime.
// -----------------------------------------------------------------------------
void
AMRNSLevel::tagCells(IntVectSet& a_tags)
{
    BEGIN_FLOWCHART();
    TODONOTE("Clean up tagCells and create vorticity tagging.");

    // Gather and prepare data structures.
    const ProblemContext* ctx       = ProblemContext::getInstance();
    const int             growTags  = ctx->amr.growTags;

    const ProblemDomain&     domain     = m_levGeoPtr->getDomain();
    const Box                domBox     = domain.domainBox();
    const Box                growDomBox = grow(domBox, IntVect::Unit);
    const DisjointBoxLayout& grids      = m_levGeoPtr->getBoxes();
    const IntVect&           refRatio   = this->getFineRefRatio();

    // Find the domain's interior cells.
    Box domInteriorBox = domBox;
    for (int dir = 0; dir < SpaceDim; ++dir) {
        if (!domain.isPeriodic(dir)) {
            domInteriorBox.grow(-BASISV(dir));
        }
    }

    // Set all BCs.
    this->setBC(*m_statePtr, m_time);


    // // Specialized tagging...

    // // Tag bottom boundary.
    // if (0) {
    //     Box tagBox(domBox);
    //     tagBox.setRange(SpaceDim - 1, domBox.smallEnd()[SpaceDim - 1], 4);
    //     a_tags |= tagBox;
    // }

    // // Tag focusing region
    // if (0) {
    //     const IntVect& sm = domBox.smallEnd();
    //     // const IntVect& bg = domBox.bigEnd();
    //     const IntVect& sz = domBox.size();
    //     const IntVect  md = sm + sz / 2;
    //     const int width = sz[0] / 8;

    //     Box tagBox(domBox);
    //     tagBox.setRange(0, md[0] - width / 2, width);
    //     if (SpaceDim == 3) {
    //         tagBox.setRange(1, md[1] - width / 2, width);
    //     }
    //     tagBox.setRange(SpaceDim - 1, sm[SpaceDim - 1], sz[SpaceDim - 1] / 4);

    //     a_tags |= tagBox;
    // }

    // // Random tagging
    // if (0) {
    //     int levelTaggingLuck = rand() % 10 + 1;

    //     if (levelTaggingLuck > 3) {
    //         const IntVect& smallEnd = m_problem_domain.domainBox().smallEnd();
    //         const IntVect& size = m_problem_domain.size();

    //         const int maxLevelTags = 4*numProc();
    //         for (int idx = 0; idx < maxLevelTags; ++idx) {
    //             IntVect levelTag;
    //             for (int dir = 0; dir < SpaceDim; ++dir) {
    //                 levelTag[dir] = rand() % size[dir] + smallEnd[dir];
    //             }
    //             a_tags |= levelTag;
    //         }
    //     }
    // }
    // // a_tags &= domBox;
    // // return;

    // // L-shaped tag in middle of domain
    // if (0) {
    //     Box region;

    //     // for (int dir = 0; dir < SpaceDim; ++dir) {
    //     for (int dir = 1; dir <= 1; ++dir) {
    //         region = m_problem_domain.domainBox();

    //         for (int growDir = 0; growDir < SpaceDim; ++growDir) {
    //             region.growLo(growDir, -m_problem_domain.size(growDir) / 2.6);
    //             if (growDir != dir) {
    //                 region.growHi(growDir, -m_problem_domain.size(growDir) / 2.6);
    //             }
    //         }
    //         region.shift(-8*IntVect::Unit);

    //         Box restrictBox = m_problem_domain.domainBox();
    //         // restrictBox.grow(-16);
    //         region &= restrictBox;

    //         denseTags |= region;
    //     }
    // }

    // // Plus-shaped tag in middle of domain
    // if (0) {
    //     const IntVect& smallEnd = m_problem_domain.domainBox().smallEnd();
    //     const IntVect& bigEnd   = m_problem_domain.domainBox().bigEnd();
    //     const IntVect  middle   = (smallEnd + bigEnd) / 2;
    //     // const IntVect& size     = m_problem_domain.domainBox().size();

    //     for (int dir = 0; dir < 2; ++dir) {
    //         const int width = 4;
    //         const int sm = middle[dir] - width;
    //         const int bg = middle[dir] + width;

    //         Box region = m_problem_domain.domainBox();
    //         region.setSmall(dir, sm);
    //         region.setBig(dir, bg);

    //         Box restrictBox = m_problem_domain.domainBox();
    //         // restrictBox.grow(-16);
    //         region &= restrictBox;

    //         denseTags |= region;
    //     }
    // }


    // Standard tagging...

    TODONOTE("Should we check differences in dirs that we aren't refining?");

    // Tag on velocity differences
    const Real velTagTol = ctx->amr.velTagTol;
    if (velTagTol > 0.0) {
        for (DataIterator dit(grids); dit.ok(); ++dit) {
            const Box ccTagRegion = grids[dit] & domInteriorBox;

            for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                const IntVect    ef     = BASISV(velComp);
                const FArrayBox& velFAB = m_statePtr->vel[dit][velComp];
                const Box fcTagRegion = surroundingNodes(ccTagRegion, velComp);
                Real vel0, dif;

                if (refRatio[velComp] > 1) {
                    for(BoxIterator bit(ccTagRegion); bit.ok(); ++bit) {
                        const IntVect& cc = bit();
                        dif = abs(velFAB(cc+ef) - velFAB(cc));
                        if (dif >= velTagTol) {
                            a_tags |= cc;
                        }
                    } // bit
                }

                for (int diffOffset = 1; diffOffset < SpaceDim; ++diffOffset) {
                    const int diffDir = (velComp + diffOffset) % SpaceDim;
                    const IntVect ed = BASISV(diffDir);

                    if (refRatio[diffDir] == 1) continue; // Skip non-refined dirs.

                    for(BoxIterator bit(fcTagRegion); bit.ok(); ++bit) {
                        const IntVect& fc = bit();
                        vel0 = velFAB(fc,0);

                        dif = abs(velFAB(fc+ed) - vel0);
                        if (dif >= velTagTol) {
                            a_tags |= fc      + ed;
                            a_tags |= fc          ;
                            a_tags |= fc - ef + ed;
                            a_tags |= fc - ef     ;
                        }

                        dif = abs(vel0 - velFAB(fc-ed,0));
                        if (dif >= velTagTol) {
                            a_tags |= fc          ;
                            a_tags |= fc      - ed;
                            a_tags |= fc - ef     ;
                            a_tags |= fc - ef - ed;
                        }
                    } // bit
                } // diffOffset, diffDir
            } // velComp
        } // dit
    }

    // Function to do tagging on generic CC data holder.
    auto doQTagging = [&domInteriorBox, &a_tags, &refRatio](
                          const LevelData<FArrayBox>& a_q,
                          const Real                  a_tagTol) {
        if (a_tagTol > 0.0) {
            for (DataIterator dit = a_q.dataIterator(); dit.ok(); ++dit) {
                const FArrayBox& qFAB = a_q[dit];
                const Box tagRegion = a_q.getBoxes()[dit] & domInteriorBox;

                IntVect left, right;
                Real dif;

                for (int dir = 0; dir < SpaceDim; ++dir) {
                    if (refRatio[dir] == 1) continue; // Skip non-refined dirs.

                    for (BoxIterator bit(tagRegion); bit.ok(); ++bit) {
                        right = bit();
                        left  = bit() - BASISV(dir);

                        dif = abs(qFAB(right) - qFAB(left));

                        if (dif >= a_tagTol) {
                            a_tags |= left;
                            a_tags |= right;
                        }
                    } // end loop over box
                } // end loop over difference dirs
            } // dit
        } // end if tagTol > 0
    };

    // Tag on T.
    LevelData<FArrayBox> T(grids, 1, IntVect::Unit);
    this->fillTemperature(T, m_time, true);
    doQTagging(T, ctx->amr.TTagTol);

    // Tag on S.
    LevelData<FArrayBox> S(grids, 1, IntVect::Unit);
    this->fillSalinity(S, m_time, true);
    doQTagging(S, ctx->amr.STagTol);

    // Tag on b.
    LevelData<FArrayBox> b(grids, 1, IntVect::Unit);
    for (DataIterator dit(grids); dit.ok(); ++dit) {
        FArrayBox zFAB(b[dit].box(), 1);
        m_levGeoPtr->fill_physCoor(zFAB, 0, SpaceDim - 1);
        this->equationOfState(b[dit], T[dit], S[dit], zFAB);
    }
    doQTagging(b, ctx->amr.bTagTol);

    // Tag on Tpert.
    Subspace::addHorizontalExtrusion(T, 0, *m_TbarPtr, 0, 1, -1.0);
    doQTagging(T, ctx->amr.TpertTagTol);

    // Tag on Spert.
    Subspace::addHorizontalExtrusion(S, 0, *m_SbarPtr, 0, 1, -1.0);
    doQTagging(S, ctx->amr.SpertTagTol);

    // Tag on bpert.
    Subspace::addHorizontalExtrusion(b, 0, *m_bbarPtr, 0, 1, -1.0);
    doQTagging(b, ctx->amr.bpertTagTol);

    // Tag on scalars
    const int N = this->numScalars();
    if (N > 0 && ctx->amr.scalarsTagTol.size() > 0) {
        if ((int)ctx->amr.scalarsTagTol.size() != N) {
            MAYDAYERROR(
                "ctx->amr.scalarsTagTol.size() must equal this->numScalars()."
                << "\nctx->amr.scalarsTagTol.size() = "
                << ctx->amr.scalarsTagTol.size()
                << "\nthis->numScalars() = " << N);
        }

        LevelData<FArrayBox> scalars(grids, N, IntVect::Unit);
        this->fillScalar(scalars, m_time, true);

        for (int comp = 0; comp < this->numScalars(); ++comp) {
            const Real tol = ctx->amr.scalarsTagTol[comp];
            if (tol > 0.0) {
                LevelData<FArrayBox> scalComp;
                aliasLevelData(scalComp, &scalars, Interval(comp, comp));
                doQTagging(scalComp, tol);
            }
        }
    }

    // Finalize.
    a_tags.grow(growTags);
    a_tags &= domBox;
}



// -----------------------------------------------------------------------------
// Perform any pre-regridding ops -- lBase is the finest unchanged level.
// This function will be called from the finest level downward to a_lBase.
// NOTE: If a new level is about to be created, this function will be run on
// that level as well...even if it doesn't exist yet!
// -----------------------------------------------------------------------------
void
AMRNSLevel::preRegrid(int                        a_lBase,
                      const Vector<Vector<Box>>& a_newGrids)
{
    BEGIN_FLOWCHART();

    IO::tout(0) << Format::pushFlags << Format::yellow
                << "AMRNSLevel::preRegrid" << Format::none << Format::popFlags
                << endl;

    // Save all data that will be regridded.
    if (m_level > a_lBase && !a_newGrids.empty()) {
        const DisjointBoxLayout& grids = this->getBoxes();

        m_oldVelPtr.reset(new LevelData<FluxBox>(grids, m_velPtr->nComp(), m_velPtr->ghostVect()));
        m_oldPPtr.reset(new LevelData<FArrayBox>(grids, m_pPtr->nComp(), m_pPtr->ghostVect()));
        m_oldQPtr.reset(new LevelData<FArrayBox>(grids, m_qPtr->nComp(), m_qPtr->ghostVect()));

        for (DataIterator dit(grids); dit.ok(); ++dit) {
            (*m_oldVelPtr)[dit].copy((*m_velPtr)[dit]);
            (*m_oldPPtr)[dit].copy((*m_pPtr)[dit]);
            (*m_oldQPtr)[dit].copy((*m_qPtr)[dit]);
        }

        // Convert vel to advecting velocity to help when refining.
        this->sendToAdvectingVelocity(*m_oldVelPtr, *m_oldVelPtr);
        debugCheckValidFaceOverlap(*m_oldVelPtr);
    }
}


// -----------------------------------------------------------------------------
// Set up data on this level after regridding.
// Called from coarsest regridded level --> finest regridded level
// Load balancing occurs here.
// -----------------------------------------------------------------------------
void
AMRNSLevel::regrid(const Vector<Box>& a_new_grids)
{
    BEGIN_FLOWCHART();

    IO::tout(0) << Format::pushFlags << Format::yellow << "AMRNSLevel::regrid"
                << Format::none << Format::popFlags << endl;

    if (s_verbosity >= 6) {
        pout() << "AMRNSLevel::regrid " << m_level
               << " a_new_grids.size() =" << a_new_grids.size()
               << endl;
    }

    // Check if level is inactive.
    if (a_new_grids.size() == 0) {
        this->deactivateLevel();
    } else {
        // Reallocate data holders, interpolators, etc...
        // WARNING: This deactivates all levels above *this!
        this->activateLevel(a_new_grids);

        // Interpolate coarse-level data.
        {
            const AMRNSLevel*        crsePtr   = this->crseNSPtr();
            const DisjointBoxLayout& crseGrids = crsePtr->getBoxes();

            debugInitLevel(*m_velPtr);
            debugInitLevel(*m_pPtr);
            debugInitLevel(*m_qPtr);

            // Get coarse vel and interpolate up.
            // Leave m_velPtr as advecting vel for now.
            {
                LevelData<FluxBox> crseVel(crseGrids,
                                           m_velPtr->nComp(),
                                           IntVect::Unit);
                crsePtr->fillVelocity(crseVel, m_time);
                crsePtr->sendToAdvectingVelocity(crseVel, crseVel);

                m_cfInterpPtr->refine(*m_velPtr, crseVel);
            }

            // Get coarse P and interpolate up.
            // This will be made more accurate in postRegrid().
            {
                LevelData<FArrayBox> crseP(crseGrids,
                                           m_pPtr->nComp(),
                                           IntVect::Unit);
                crsePtr->fillPressure(crseP, m_time, false);
                m_cfInterpPtr->refine(*m_pPtr, crseP);
            }

            // Get coarse Q and interpolate up.
            {
                LevelData<FArrayBox> crseQ(crseGrids,
                                           m_qPtr->nComp(),
                                           IntVect::Unit);
                setValLevel(crseQ, 0.0); // eddyNu = 0 for now on new level.

                if (this->numScalars() > 0) {
                    LevelData<FArrayBox> crseS;
                    const Interval& sivl = crsePtr->m_statePtr->scalarsInterval;
                    aliasLevelData(crseS, &crseQ, sivl);
                    crsePtr->fillScalar(crseS, m_time);
                }

                LevelData<FArrayBox> crseT;
                aliasLevelData(crseT, &crseQ, crsePtr->m_statePtr->TInterval);
                crsePtr->fillTemperature(crseT, m_time);

                LevelData<FArrayBox> crseS;
                aliasLevelData(crseS, &crseQ, crsePtr->m_statePtr->SInterval);
                crsePtr->fillSalinity(crseS, m_time);

                m_cfInterpPtr->refine(*m_qPtr, crseQ);
            }
        } // end interpolate up

        // Copy fine-level data where it already existed and free old storage.
        if (m_oldVelPtr) {
            Copier oldToNewCopier;
            {
                const DisjointBoxLayout& oldGrids = m_oldVelPtr->getBoxes();
                const DisjointBoxLayout& newGrids = m_velPtr->getBoxes();
                const ProblemDomain& domain = newGrids.physDomain();
                oldToNewCopier.define(oldGrids, newGrids, domain);
            }

            debugCheckValidFaceOverlap(*m_oldVelPtr);
            m_oldVelPtr->copyTo(*m_velPtr);
            this->sendToCartesianVelocity(*m_velPtr, *m_velPtr);
            debugCheckValidFaceOverlap(*m_velPtr);

            m_oldPPtr->copyTo(*m_pPtr, oldToNewCopier);
            m_oldQPtr->copyTo(*m_qPtr, oldToNewCopier);
        }

        // Pressure estimates and composite projection done in postRegrid().

        // Free memory
        m_oldVelPtr.reset();
        m_oldPPtr.reset();
        m_oldQPtr.reset();
    }
}


// -----------------------------------------------------------------------------
// Perform any post-regridding ops -- lbase is the finest unchanged level
// Called from finest regridded level --> a_lbase (coarsest regridded level - 1)
// -----------------------------------------------------------------------------
void
AMRNSLevel::postRegrid(int a_lbase)
{
    BEGIN_FLOWCHART();

    if (s_verbosity >= 3) {
        pout() << "AMRNSLevel::postRegrid on level " << m_level
               << " with a_lbase = " << a_lbase << endl;
    }
    const ProblemContext* ctx = ProblemContext::getInstance();

    // Validate the ops and solvers.
    // This is a bit trickier than usual because a level may have been created,
    // so be safe -- re-validate everything from lbase --> finest.
    if (this->isFinestLevel()) {
        AMRNSLevel* levPtr = this;
        while (levPtr) {
            levPtr->m_amrOpsValid       = false;
            levPtr->m_levelSolversValid = false;
            levPtr->m_amrSolversValid   = false;
            levPtr = levPtr->crseNSPtr();
        }
        this->coarsestNSPtr()->validateOpsAndSolvers();
    }

    // Project from lBase + 1 to finest level.
    if (m_level == a_lbase + 1 && ctx->proj.doInitProj) {
        if (s_verbosity >= 2) {
            pout() << "Initialization projection:" << endl;
        }
        this->projectDownToThis();
        pout() << endl;
    }

    // Pressure estimation from lBase to finest level.
    if (m_level == a_lbase && ctx->rhs.computeInitPressure) {
        this->initAMRPressure(0.5);
    }
}


// -----------------------------------------------------------------------------
// Create a load-balanced DisjointBoxLayout from a collection of boxes
// -----------------------------------------------------------------------------
DisjointBoxLayout
AMRNSLevel::loadBalance(const Vector<Box>& a_grids)
{
    if (s_verbosity >= 6) {
        pout() << "AMRNSLevel::loadBalance " << m_level
               << "  boxes=" << a_grids.size() << endl;
    }

    Vector<int> proc_map;
    LoadBalance(proc_map, a_grids);

    if (s_verbosity >= 5) {
        pout() << "AMRNSLevel::loadBalance: processor map: " << endl;
        for (unsigned int igrid = 0; igrid < a_grids.size(); ++igrid) {
            pout() << igrid << ": " << proc_map[igrid] << "\t" << a_grids[igrid]
                   << endl;
        }
        pout() << endl;
    }

    return DisjointBoxLayout(a_grids, proc_map, m_problem_domain);
}


// -----------------------------------------------------------------------------
void
AMRNSLevel::initLevelPressure(const Real a_dt)
{
    if (s_verbosity >= 2) {
        pout() << "Pressure initialization on level " << m_level << ":"
               << Format::indent() << endl;
    }

    // Sanity checks
    CH_assert(!RealCmp::isZero(a_dt));
    CH_assert(m_levelProjSolverPtr);

    // Gather references
    const DisjointBoxLayout& grids = m_statePtr->grids;
    LevelData<FluxBox>&      vel   = m_statePtr->vel;
    LevelData<FArrayBox>&    q     = m_statePtr->q;
    LevelData<FArrayBox>&    p     = m_statePtr->p;

    // Save current data
    LevelData<FluxBox> saveVel(grids, vel.nComp(), vel.ghostVect());
    LevelData<FArrayBox> saveQ(grids, q.nComp(), q.ghostVect());
    for (DataIterator dit(grids); dit.ok(); ++dit) {
        saveVel[dit].copy(vel[dit]);
        saveQ[dit].copy(q[dit]);
    }
    Real saveTime = m_time;
    Real saveDt = m_dt;

    // Set temporary timestep.
    this->dt(a_dt);

    // I've tried to make this as cheap as possible by forcing the
    // level projector to run just one iteration.
    // I've also swapped advance() with FEadvance(), which is less
    // accurate, but much faster.
    // Keep in mind that this only needs to produce a 1st order accurate
    // (in time) pressure to produce a 2nd order accurate pressure in the
    // main timestepper.
    if constexpr (1) {
        // Cheap version...
        const auto saveOpts = m_levelProjSolverPtr->getOptions();
        m_levelProjSolverPtr->modifyOptionsExceptMaxDepth(LevelProjSolver::getQuickAndDirtyOptions());
        m_parkPtr->FEadvance(vel, p, q, m_time, m_dt, this);
        m_levelProjSolverPtr->modifyOptionsExceptMaxDepth(saveOpts);

    } else {
        // Accurate version...
        m_parkPtr->advance(vel, p, q, m_time, m_dt, this);
    }

    // Restore data
    this->time(saveTime);
    this->dt(saveDt);
    for (DataIterator dit(grids); dit.ok(); ++dit) {
        vel[dit].copy(saveVel[dit]);
        q[dit].copy(saveQ[dit]);
    }

    // Set BCs on new pressure estimate.
    this->setPressureBC(p, m_time, false);

    if (s_verbosity >= 2) {
        pout() << Format::unindent;
    }
}

// -----------------------------------------------------------------------------
void
AMRNSLevel::initAMRPressure(const Real a_dtMult)
{
    if (s_verbosity >= 2) {
        pout() << "Pressure initialization on levels " << m_level
                << " to "
                << this->finestNSPtr()->m_level
                << ":" << Format::indent() << endl;
    }
    auto start = std::chrono::high_resolution_clock::now();

    // Compute small timestep
    Real smallDt = this->finestNSPtr()->dt();
    if (RealCmp::isZero(smallDt)) {
        smallDt = this->finestNSPtr()->computeDt();
    }
    CH_assert(!RealCmp::isZero(smallDt));

    CH_assert(a_dtMult > 0.0);
    smallDt *= a_dtMult;

    // Climb up the levels, calling the level init routine.
    AMRNSLevel* levPtr = this;
    while (levPtr) {
        levPtr->initLevelPressure(smallDt);
        levPtr = levPtr->fineNSPtr();
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    if (s_verbosity >= 2) {
        pout() << "AMR pressure init took " << elapsed.count() << " s\n"
                << Format::unindent << endl;
    }
}
