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
#include "Subspace.H"


// -----------------------------------------------------------------------------
// Static utility
// Splits the domain into its sponge layers and the interior. The locations of
// the splitting are given as face indices. If the domain is not split, these
// indices will lie outside of the domain.
// You almost certainly do NOT want to override this.
// -----------------------------------------------------------------------------
void AMRNSLevel::computeSpongeRegions(Tuple<Box, 2>&       a_spongeBox,
                                      Tuple<int, 2>&       a_splitFaceIndex,
                                      Box&                 a_interior,
                                      const int            a_dir,
                                      const LevelGeometry& a_levGeo)
{
    // Sanity checks
    CH_assert(0 <= a_dir);
    CH_assert(a_dir < SpaceDim);

    // Gather some needed info
    const ProblemContext* ctx    = ProblemContext::getInstance();
    const Box             domBox = a_levGeo.getDomainBox();
    const RealVect&       dXi    = a_levGeo.getDXi();
    Box                   interiorBox[2];

    // Lo side
    int s = 0;
    Real spongeWidth = ctx->rhs.spongeWidth[a_dir][s];
    if (spongeWidth > 0.0) {
        const int numCells = ceil(spongeWidth / dXi[a_dir]);
        a_splitFaceIndex[s] = domBox.smallEnd(a_dir) + numCells;
        a_spongeBox[s] = domBox;
        interiorBox[s] = a_spongeBox[s].chop(a_dir, a_splitFaceIndex[s]);
    }

    // Hi side
    s = 1;
    spongeWidth = ctx->rhs.spongeWidth[a_dir][s];
    if (spongeWidth > 0.0) {
        const int numCells = ceil(spongeWidth / dXi[a_dir]);
        a_splitFaceIndex[s] = domBox.bigEnd(a_dir) - numCells + 1;
        interiorBox[s] = domBox;
        a_spongeBox[s] = interiorBox[s].chop(a_dir, a_splitFaceIndex[s]);
    }

    // Compute the interior region
    a_interior = interiorBox[0] & interiorBox[1];
}


// -----------------------------------------------------------------------------
// By default, this does nothing. You can override this in your physics
// class if you need to perform global computations before filling targets.
// -----------------------------------------------------------------------------
void
AMRNSLevel::prepareSpongeTargets(const LevelData<FluxBox>&   /*a_vel*/,
                                 const LevelData<FArrayBox>& /*a_p*/,
                                 const LevelData<FArrayBox>& /*a_q*/,
                                 const Real                  /*a_time*/)
{
    // Do nothing.
}


// -----------------------------------------------------------------------------
// By default, this throws an error for all BCs except TIDAL.
// You can override this to fit your own needs.
// -----------------------------------------------------------------------------
void
AMRNSLevel::fillVelSpongeTarget(FArrayBox&       a_targetFAB,
                                const int        a_targetFABComp,
                                const int        a_velComp,
                                const FArrayBox& a_physCoorFAB[[maybe_unused]],
                                const Box&       a_spongeBox,
                                const DataIndex& a_di[[maybe_unused]],
                                const Real       a_time,
                                const int        a_bdryDir,
                                const Side::LoHiSide& a_side) const
{
    // Sanity checks
    CH_assert(0 <= a_targetFABComp);
    CH_assert(a_targetFABComp < a_targetFAB.nComp());

    CH_assert(0 <= a_velComp);
    CH_assert(a_velComp < SpaceDim);

    CH_assert(a_physCoorFAB.nComp() == SpaceDim);

    CH_assert(0 <= a_bdryDir);
    CH_assert(a_bdryDir < SpaceDim);

    CH_assert(a_spongeBox.type() == a_targetFAB.box().type());

    if (m_problem_domain.isPeriodic(a_bdryDir)) {
        MAYDAYERROR("Cannot use sponge at a periodic boundary.");
    }

    // Gather BC type.
    const int             iside = int(a_side);
    const ProblemContext* ctx   = ProblemContext::getInstance();
    const auto            btype = ctx->rhs.velBCType[a_bdryDir][iside];

    // Where are we setting the target?
    const Box fillBox = a_spongeBox & a_targetFAB.box();
    if (fillBox.isEmpty()) return;

    // Gather BC values.
    const int        numComps  = m_statePtr->vel.nComp();
    const FArrayBox& stateFAB  = m_statePtr->vel[a_di][a_velComp];
    FArrayBox alphaFAB(fillBox, numComps);
    FArrayBox betaFAB(fillBox, numComps);
    FArrayBox bcFAB(fillBox, numComps);

    m_velBCPtr->operator()(alphaFAB,
                           betaFAB,
                           bcFAB,
                           stateFAB,
                           a_physCoorFAB,
                           a_di,
                           a_bdryDir,
                           a_side,
                           a_time,
                           false); // not homog

    // If BCs are Neumann or Robin, throw error.
    if (!RealCmp::isZero(betaFAB.norm(2))) {
        MAYDAYERROR(
            "You must override AMRNSLevel::fillVelSpongeTarget for "
            "velBCType["
            << a_bdryDir << "][" << iside << "] = " << btype
            << ". This seems to be a Neumann or Robin BC.");
    }

    // Set Dirichlet BCs.
    a_targetFAB.copy(bcFAB, a_targetFABComp, a_targetFABComp, 1);
    a_targetFAB.divide(alphaFAB, a_targetFABComp, a_targetFABComp, 1);
}


// -----------------------------------------------------------------------------
// By default, this sets T and S to the stratification and does nothing
// to the other q comps. You can override this to fit your own needs.
// -----------------------------------------------------------------------------
void
AMRNSLevel::fillQSpongeTarget(
    FArrayBox&            a_targetFAB,
    const int             a_targetFABComp,
    const int             a_qComp,
    const FArrayBox&      a_physCoorFAB[[maybe_unused]],
    const Box&            a_spongeBox,
    const DataIndex&      a_di[[maybe_unused]],
    const Real            a_time[[maybe_unused]],
    const int             a_bdryDir,
    const Side::LoHiSide& a_side[[maybe_unused]]) const
{
    // Sanity checks
    CH_assert(0 <= a_targetFABComp);
    CH_assert(a_targetFABComp < a_targetFAB.nComp());

    CH_assert(0 <= a_qComp);
    CH_assert(a_qComp < m_statePtr->q.nComp());

    CH_assert(a_physCoorFAB.nComp() == SpaceDim);

    CH_assert(0 <= a_bdryDir);
    CH_assert(a_bdryDir < SpaceDim);

    CH_assert(a_spongeBox.type() == a_targetFAB.box().type());

    // Where are we setting the target?
    const Box fillBox = a_spongeBox & a_targetFAB.box();
    if (fillBox.isEmpty()) return;

    if (a_bdryDir < SpaceDim - 1) {
        if (a_qComp == m_statePtr->TComp) {
            // Temperature. Set target to background.
            Subspace::horizontalExtrusion(
                a_targetFAB, a_targetFABComp, *m_TbarPtr, 0, 1);
        } else if (a_qComp == m_statePtr->SComp) {
            // Salinity. Set target to background.
            Subspace::horizontalExtrusion(
                a_targetFAB, a_targetFABComp, *m_SbarPtr, 0, 1);
        } else {
            // Do nothing for all other comps.
        }
    } else {
        MAYDAYERROR(
            "AMRNSLevel::fillQSpongeTarget is not yet written for the vertical "
            "boundaries.");
    }
}


// -----------------------------------------------------------------------------
// Adds ramp * (target - state) * invTimeScale to the forcing, a_k*.
// invTimeScale is typically 1.0 / (c * dt), where c ~ 10 or so.
// By default, this only operates on vel and b. Override if you like.
// -----------------------------------------------------------------------------
void
AMRNSLevel::addSpongeForcing(LevelData<FluxBox>&         a_kvel,
                             LevelData<FArrayBox>&       a_kq,
                             const LevelData<FluxBox>&   a_vel,
                             const LevelData<FArrayBox>& a_p,
                             const LevelData<FArrayBox>& a_q,
                             const Real                  a_time,
                             const Real                  a_invTimeScale)
{
    // Gather references, etc.
    const int                numQComps = a_q.nComp();
    const RealVect&          dXi       = m_levGeoPtr->getDXi();
    const Box&               domBox    = this->getDomainBox();
    const DisjointBoxLayout& grids     = this->getBoxes();
    DataIterator             dit       = grids.dataIterator();
    const ProblemContext*    ctx       = ProblemContext::getInstance();

    // Sanity checks
    CH_assert(a_kvel.getBoxes() == grids);
    CH_assert(a_kq  .getBoxes() == grids);
    CH_assert(a_vel .getBoxes() == grids);
    CH_assert(a_q   .getBoxes() == grids);
    CH_assert(a_kvel.nComp() == 1);
    CH_assert(a_kvel.nComp() == a_vel.nComp());
    CH_assert(a_kq  .nComp() == a_q  .nComp());

    // Calculate the sponge regions and the interior
    Tuple<Tuple<Box, 2>, CH_SPACEDIM> spongeBoxes;
    Tuple<Tuple<int, 2>, CH_SPACEDIM> splitFaceIndices;
    Tuple<Box, CH_SPACEDIM> interiors;
    for (int bdryDir = 0; bdryDir < SpaceDim; ++bdryDir) {
        this->computeSpongeRegions(spongeBoxes[bdryDir],
                                   splitFaceIndices[bdryDir],
                                   interiors[bdryDir],
                                   bdryDir,
                                   *m_levGeoPtr);
    }

    // Get physical coordinates.
    LevelData<FluxBox> fcPhysCoor(grids, SpaceDim, IntVect::Unit);
    LevelData<FArrayBox> ccPhysCoor(grids, SpaceDim, IntVect::Unit);
    for (dit.reset(); dit.ok(); ++dit) {
        for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
            m_levGeoPtr->fill_physCoor(fcPhysCoor[dit][fcDir]);
        }
        m_levGeoPtr->fill_physCoor(ccPhysCoor[dit]);
    } // dit

    // Allocate workspace.
    LevelData<FluxBox> velWS(grids, 1);
    LevelData<FArrayBox> qWS(grids, numQComps);

    // Prepare target values globally (just in case the user needs to perform
    // global computations before filling the targets).
    this->prepareSpongeTargets(a_vel, a_p, a_q, a_time);

    for (int bdryDir = 0; bdryDir < SpaceDim; ++bdryDir) {
        for (SideIterator sit; sit.ok(); ++sit) {
            const Side::LoHiSide& side      = sit();
            const int             iside     = int(side);
            const Box&            spongeBox = spongeBoxes[bdryDir][iside];
            const Real spongeWidth      = ctx->rhs.spongeWidth[bdryDir][iside];
            const Real spongeWidthCells = spongeWidth / dXi[bdryDir];

            // Fill targets.
            for (dit.reset(); dit.ok(); ++dit) {
                // Start clean.
                velWS[dit].setVal(0.0);
                qWS[dit].setVal(0.0);

                // vel...
                for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                    this->fillVelSpongeTarget(
                        velWS[dit][velComp],
                        0,
                        velComp,
                        fcPhysCoor[dit][velComp],
                        surroundingNodes(spongeBox, velComp),
                        dit(),
                        a_time,
                        bdryDir,
                        side);
                }

                // q...
                for (int qComp = 0; qComp < numQComps; ++qComp) {
                    this->fillQSpongeTarget(qWS[dit],
                                            qComp,
                                            qComp,
                                            ccPhysCoor[dit],
                                            spongeBox,
                                            dit(),
                                            a_time,
                                            bdryDir,
                                            side);
                }
            }

            // Add sponge force.
            Real endPos;
            if (iside == 0) {
                endPos = Real(domBox.smallEnd(bdryDir));
            } else {
                endPos = Real(domBox.bigEnd(bdryDir) + 1);
            }

            for (dit.reset(); dit.ok(); ++dit) {
                // vel...
                for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                    FArrayBox&       kvelFAB   = a_kvel[dit][velComp];
                    const FArrayBox& targetFAB = velWS[dit][velComp];
                    const FArrayBox& stateFAB  = a_vel[dit][velComp];
                    const Box fcSpongeBox = surroundingNodes(spongeBox, velComp);
                    const Box region      = fcSpongeBox & velWS[dit][velComp].box();
                    Real  pos, r;

                    Real offset, invLength;
                    if (bdryDir == velComp) {
                        offset = 0.0;
                        invLength = 1.0 / (spongeWidthCells + 1.0);
                    } else {
                        offset = 0.5;
                        invLength = 1.0 / spongeWidthCells;
                    }

                    for (BoxIterator bit(region); bit.ok(); ++bit) {
                        const IntVect& fc = bit();
                        pos = Real(fc[bdryDir]) + offset;
                        r = 1.0 - abs(pos - endPos) * invLength;
                        r = this->spongeRampFunction(r);

                        kvelFAB(fc)
                            += r * a_invTimeScale
                             * (targetFAB(fc) - stateFAB(fc));
                    }
                }

                // Do these scalar comps too...
                Vector<int> vComps(2);
                vComps[0] = m_statePtr->TComp; // Temperature
                vComps[1] = m_statePtr->SComp; // Salinity

                {
                    FArrayBox&       kqFAB     = a_kq[dit];
                    const FArrayBox& targetFAB = qWS[dit];
                    const FArrayBox& stateFAB  = a_q[dit];
                    const Box        region    = spongeBox & qWS[dit].box();
                    Real  pos, r;

                    Real invLength = 1.0 / (spongeWidthCells + 1.0);

                    TODONOTE("Move function call out of loop and use fortran.");
                    for (const int qComp : vComps) {
                        for (BoxIterator bit(region); bit.ok(); ++bit) {
                            const IntVect& cc = bit();
                            pos = Real(cc[bdryDir]) + 0.5;
                            r = 1.0 - abs(pos - endPos) * invLength;
                            r = this->spongeRampFunction(r);

                            kqFAB(cc, qComp)
                                += r * a_invTimeScale
                                 * (targetFAB(cc, qComp) - stateFAB(cc, qComp));
                        }
                    } // qComp
                }
            } // dit
        } // sit
    } // bdryDir
}

