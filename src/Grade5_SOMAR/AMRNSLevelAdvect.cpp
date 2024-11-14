/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2022
 *    Thomas Jefferson University and
 *    University of North Carolina at Chapel Hill
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
 *  https://github.com/UNC-CFD/LES-SOMAR.
 ******************************************************************************/
#include "AMRNSLevel.H"
#include "AMRNSLevelF_F.H"
#include "Convert.H"
#include "FourthOrder.H"
#include "Debug.H"


// =============================================================================
// Momentum advection

// -----------------------------------------------------------------------------
void
AMRNSLevel::computeMomentumAdvection(LevelData<FluxBox>&       a_kvel,
                                     StaggeredFluxLD&          a_momentumFlux,
                                     const LevelData<FluxBox>& a_cartVel,
                                     const LevelData<FluxBox>& a_advVel) const
{
    CH_assert(a_kvel.nComp() == 1);
    CH_assert(a_cartVel.nComp() == 1);
    CH_assert(a_advVel.nComp() == 1);

    debugCheckValidFaceOverlap(a_kvel);
    debugCheckValidFaceOverlap(a_cartVel);
    debugCheckValidFaceOverlap(a_advVel);

    static const RHSParameters& rhs = ProblemContext::getInstance()->rhs;
    static const auto           velReconstruction = rhs.velReconstruction;
    static const Real           momAdvSkewness    = rhs.momAdvSkewness;

    // Compute fluxes
    switch (velReconstruction) {
        case RHSParameters::Reconstruction::SecondOrder:
            this->momAdvFlux_Cons_Form_SecondOrder(
                a_momentumFlux, a_cartVel, a_advVel);
            break;
        case RHSParameters::Reconstruction::FourthOrder:
            this->momAdvFlux_Cons_Form_FourthOrder(
                a_momentumFlux, a_cartVel, a_advVel);
            break;
    };

    // Compute divergence
    if (momAdvSkewness < 1.0) {
        // Has a conservative form component.
        const bool scaleByJinv = false; // AMRNSLevelRHS will do this.
        const bool addTokvel   = true;
        m_finiteDiffPtr->levelVectorDivergence(a_kvel,
                                               a_momentumFlux,
                                               scaleByJinv,
                                               addTokvel,
                                               1.0 - momAdvSkewness);
        // m_finiteDiffPtr->levelVectorDivergence_4thOrder(a_kvel,
        //                                                 a_momentumFlux,
        //                                                 scaleByJinv,
        //                                                 addTokvel,
        //                                                 1.0 - momAdvSkewness);
    }

    if (momAdvSkewness > 0.0) {
        // Has an advective form component.
        this->momAdvForce_AdvForm_SecondOrder(a_kvel,
                                              a_cartVel,
                                              a_advVel,
                                              momAdvSkewness);
    }

    // debugCheckValidFaceOverlap(a_kvel);
    LayoutTools::averageOverlappingValidFaces(a_kvel);
}


// -----------------------------------------------------------------------------
void
AMRNSLevel::momAdvFlux_Cons_Form_SecondOrder(
    StaggeredFluxLD&          a_momentumFlux,
    const LevelData<FluxBox>& a_cartVel,
    const LevelData<FluxBox>& a_advVel) const
{
    const GeoSourceInterface& geoSrc = m_levGeoPtr->getGeoSource();
    const RealVect&           dXi    = m_levGeoPtr->getDXi();
    const DisjointBoxLayout&  grids  = m_levGeoPtr->getBoxes();

    for (DataIterator dit(grids); dit.ok(); ++dit) {
        // velComp = advectED vel component.
        for (int velComp = 0; velComp < SpaceDim; ++velComp) {
            const FArrayBox& uaFAB = a_cartVel[dit][velComp];

            checkForNAN(uaFAB, uaFAB.box());

            FArrayBox JFAB(uaFAB.box(), 1);
            geoSrc.fill_J(JFAB, 0, dXi);

            // derivDir = advectING vel component and derivative dir.
            for (int derivDir = 0; derivDir < SpaceDim; ++derivDir) {
                const FArrayBox& JubFAB = a_advVel[dit][derivDir];
                FArrayBox& fluxFAB    = a_momentumFlux[derivDir][velComp][dit];
                const Box& fluxRegion = fluxFAB.box();
                const Real dXiDir     = dXi[derivDir];

                checkForNAN(JubFAB, JubFAB.box());
                debugInit(fluxFAB);

                FORT_MOMADVFLUX_SECONDORDER(
                    CHF_FRA1(fluxFAB, 0),
                    CHF_BOX(fluxRegion),
                    CHF_CONST_FRA1(JubFAB, 0),
                    CHF_CONST_FRA1(uaFAB, 0),
                    CHF_CONST_FRA1(JFAB, 0),
                    CHF_CONST_INT(derivDir),
                    CHF_CONST_INT(velComp),
                    CHF_CONST_REAL(dXiDir));
            } // derivDir
        } // velComp
    } // dit
}


#if 1
// Old version
// -----------------------------------------------------------------------------
void
AMRNSLevel::momAdvFlux_Cons_Form_FourthOrder(
    StaggeredFluxLD&          a_momentumFlux,
    const LevelData<FluxBox>& a_cartVel,
    const LevelData<FluxBox>& a_advVel) const
{
    const GeoSourceInterface& geoSrc      = m_levGeoPtr->getGeoSource();
    const RealVect&           dXi         = m_levGeoPtr->getDXi();
    const DisjointBoxLayout&  grids       = m_levGeoPtr->getBoxes();

    // Prepare source data.
    LevelData<FluxBox> cartVelGrow(grids, 1, 2 * IntVect::Unit);
    LevelData<FluxBox> advVelGrow(grids, 1, 2 * IntVect::Unit);
    {
        debugInitLevel(cartVelGrow);
        debugInitLevel(advVelGrow);

        for (DataIterator dit(grids); dit.ok(); ++dit) {
            cartVelGrow[dit].copy(a_cartVel[dit]);
            advVelGrow[dit].copy(a_advVel[dit]);
        }

        const bool    extrapOrder = 4;
        const IntVect skipGhosts  = IntVect::Unit;
        BCTools::extrapAllGhosts(cartVelGrow, extrapOrder, skipGhosts);
        BCTools::extrapAllGhosts(advVelGrow, extrapOrder, skipGhosts);

        cartVelGrow.exchange();
        advVelGrow.exchange();
    }

    for (DataIterator dit(grids); dit.ok(); ++dit) {
        // velComp = advectED vel component.
        for (int velComp = 0; velComp < SpaceDim; ++velComp) {
            const FArrayBox& uaFAB  = cartVelGrow[dit][velComp];

            checkForNAN(uaFAB, uaFAB.box());

            FArrayBox JFAB(uaFAB.box(), 1);
            geoSrc.fill_J(JFAB, 0, dXi);

            // derivDir = advectING vel component and derivative dir.
            for (int derivDir = 0; derivDir < SpaceDim; ++derivDir) {
                const FArrayBox& JubFAB = advVelGrow[dit][derivDir];
                FArrayBox& fluxFAB    = a_momentumFlux[derivDir][velComp][dit];
                const Box& fluxRegion = fluxFAB.box();
                const Real dXiDir     = dXi[derivDir];

                checkForNAN(JubFAB, JubFAB.box());
                debugInit(fluxFAB);

                FORT_MOMADVFLUX_FOURTHORDER(
                    CHF_FRA1(fluxFAB, 0),
                    CHF_BOX(fluxRegion),
                    CHF_CONST_FRA1(JubFAB, 0),
                    CHF_CONST_FRA1(uaFAB, 0),
                    CHF_CONST_FRA1(JFAB, 0),
                    CHF_CONST_INT(derivDir),
                    CHF_CONST_INT(velComp),
                    CHF_CONST_REAL(dXiDir));

            } // derivDir
        } // velComp
    } // dit
}

#else

// Fully fourth order verion
// -----------------------------------------------------------------------------
void
AMRNSLevel::momAdvFlux_Cons_Form_FourthOrder(
    StaggeredFluxLD&          a_momentumFlux,
    const LevelData<FluxBox>& a_cartVel,
    const LevelData<FluxBox>& a_advVel) const
{
    const GeoSourceInterface& geoSrc      = m_levGeoPtr->getGeoSource();
    const RealVect&           dXi         = m_levGeoPtr->getDXi();
    const DisjointBoxLayout&  grids       = m_levGeoPtr->getBoxes();

    // Prepare advectED velocity.
    LevelData<FluxBox> cartVelGrow(grids, 1, 2 * IntVect::Unit);
    {
        debugInitLevel(cartVelGrow);
        for (DataIterator dit(grids); dit.ok(); ++dit) {
            cartVelGrow[dit].copy(a_cartVel[dit]);
        }

        BCTools::extrapAllGhosts(cartVelGrow, 4, IntVect::Unit);
        // cartVelGrow.exchange();

        for (DataIterator dit(grids); dit.ok(); ++dit) {
            // velComp = advectED vel component.
            for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                FArrayBox& uaFAB  = cartVelGrow[dit][velComp];
                for (int convDir = 0; convDir < SpaceDim; ++convDir) {
                    FourthOrder::nodalToAvg(uaFAB, uaFAB.box(), convDir);
                }
            }
        }
        cartVelGrow.exchange();
    }

    // Prepare advectING velocity.
    LevelData<FluxBox> advVelGrow(grids, 1, 2 * IntVect::Unit);
    {
        debugInitLevel(advVelGrow);
        for (DataIterator dit(grids); dit.ok(); ++dit) {
            advVelGrow[dit].copy(a_advVel[dit]);
        }
        BCTools::extrapAllGhosts(advVelGrow, 4, IntVect::Unit);
        FourthOrder::nodalToAvg(advVelGrow);
        advVelGrow.exchange();
    }

    for (DataIterator dit(grids); dit.ok(); ++dit) {
        // velComp = advectED vel component.
        for (int velComp = 0; velComp < SpaceDim; ++velComp) {
            const FArrayBox& uaFAB  = cartVelGrow[dit][velComp];

            checkForNAN(uaFAB, uaFAB.box());

            FArrayBox JFAB(uaFAB.box(), 1);
            geoSrc.fill_J(JFAB, 0, dXi);

            // derivDir = advectING vel component and derivative dir.
            for (int derivDir = 0; derivDir < SpaceDim; ++derivDir) {
                const FArrayBox& JubFAB = advVelGrow[dit][derivDir];
                FArrayBox& fluxFAB    = a_momentumFlux[derivDir][velComp][dit];
                const Box  fluxRegion = fluxFAB.box(); //grids[dit].convert(fluxFAB.box().type());
                const Real dXiDir     = dXi[derivDir];

                checkForNAN(JubFAB, JubFAB.box());
                debugInit(fluxFAB);

                if (velComp == derivDir) {
                    // {
                    //     FArrayBox flux2FAB(fluxRegion, 1);
                    //     FORT_MOMADVFLUX_SECONDORDER(
                    //         CHF_FRA1(flux2FAB, 0),
                    //         CHF_BOX(fluxRegion),
                    //         CHF_CONST_FRA1(JubFAB, 0),
                    //         CHF_CONST_FRA1(uaFAB, 0),
                    //         CHF_CONST_FRA1(JFAB, 0),
                    //         CHF_CONST_INT(derivDir),
                    //         CHF_CONST_INT(velComp),
                    //         CHF_CONST_REAL(dXiDir));

                    //     fluxFAB.copy(flux2FAB);
                    //     FiniteDiff::addSecondDifference(fluxFAB,
                    //                                     0,
                    //                                     flux2FAB,
                    //                                     0,
                    //                                     fluxRegion,
                    //                                     velComp,
                    //                                     -1.0 / 32.0);
                    // }

                    // Old version. 4th-order FV recentering, 2nd-order flux calc.
                    FORT_MOMADVFLUX_FULLYFOURTHORDER(
                        CHF_FRA1(fluxFAB, 0),
                        CHF_BOX(fluxRegion),
                        CHF_CONST_FRA1(JubFAB, 0),
                        CHF_CONST_FRA1(uaFAB, 0),
                        CHF_CONST_FRA1(JFAB, 0),
                        CHF_CONST_INT(derivDir),
                        CHF_CONST_INT(velComp),
                        CHF_CONST_REAL(dXiDir));

                } else {
                    // {
                    //     FArrayBox flux2FAB(fluxRegion, 1);
                    //     FORT_MOMADVFLUX_SECONDORDER(
                    //         CHF_FRA1(flux2FAB, 0),
                    //         CHF_BOX(fluxRegion),
                    //         CHF_CONST_FRA1(JubFAB, 0),
                    //         CHF_CONST_FRA1(uaFAB, 0),
                    //         CHF_CONST_FRA1(JFAB, 0),
                    //         CHF_CONST_INT(derivDir),
                    //         CHF_CONST_INT(velComp),
                    //         CHF_CONST_REAL(dXiDir));

                    //     fluxFAB.copy(flux2FAB);
                    //     FiniteDiff::addSecondDifference(fluxFAB,
                    //                                     0,
                    //                                     flux2FAB,
                    //                                     0,
                    //                                     fluxRegion,
                    //                                     velComp,
                    //                                     1.0 / 96.0);
                    // }

                    FORT_MOMADVFLUX_FULLYFOURTHORDER(
                        CHF_FRA1(fluxFAB, 0),
                        CHF_BOX(fluxRegion),
                        CHF_CONST_FRA1(JubFAB, 0),
                        CHF_CONST_FRA1(uaFAB, 0),
                        CHF_CONST_FRA1(JFAB, 0),
                        CHF_CONST_INT(derivDir),
                        CHF_CONST_INT(velComp),
                        CHF_CONST_REAL(dXiDir));
                }
            } // derivDir
        } // velComp
    } // dit
}

#endif


// -----------------------------------------------------------------------------
void
AMRNSLevel::momAdvForce_AdvForm_SecondOrder(
    LevelData<FluxBox>&       a_kvel,
    const LevelData<FluxBox>& a_cartVel,
    const LevelData<FluxBox>& a_advVel,
    const Real                a_scale) const
{
    const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();
    const RealVect&          dXi   = m_levGeoPtr->getDXi();

    for (DataIterator dit(grids); dit.ok(); ++dit) {
        for (int velComp = 0; velComp < SpaceDim; ++velComp) { // advectED
            const FArrayBox& uFAB    = a_cartVel[dit][velComp];
            const Box        fcValid = grids[dit].surroundingNodes(velComp);

            FArrayBox kuFAB(fcValid, 1);

            for (int derivDir = 0; derivDir < SpaceDim; ++derivDir) { // advectING
                const FArrayBox& JuFAB  = a_advVel[dit][derivDir];
                const Real       dXiDir = dXi[derivDir];

                debugInit(kuFAB);

                FORT_MOMADVFORCE_ADVFORM_SECONDORDER(
                    CHF_FRA1(kuFAB, 0),
                    CHF_BOX(fcValid),
                    CHF_CONST_FRA1(JuFAB, 0),
                    CHF_CONST_FRA1(uFAB, 0),
                    CHF_CONST_INT(derivDir),
                    CHF_CONST_INT(velComp),
                    CHF_CONST_REAL(dXiDir)
                );

                a_kvel[dit][velComp].plus(kuFAB, a_scale);
            } // derivDir
        } // velComp
    } // dit
}


// =============================================================================
// Scalar advection

// -----------------------------------------------------------------------------
void
AMRNSLevel::computeScalarAdvection(LevelData<FArrayBox>&       a_kq,
                                   LevelData<FluxBox>&         a_qFlux,
                                   const LevelData<FArrayBox>& a_q,
                                   const LevelData<FluxBox>&   a_advVel) const
{
    CH_assert(a_kq.nComp() == 1);
    CH_assert(a_qFlux.nComp() == 1);
    CH_assert(a_q.nComp() == 1);

    debugInitLevel(a_kq);
    debugInitLevel(a_qFlux);

    static const RHSParameters& rhs = ProblemContext::getInstance()->rhs;
    static const auto scalReconstruction = rhs.scalReconstruction;

    // Compute fluxes
    switch (scalReconstruction) {
        case RHSParameters::Reconstruction::SecondOrder:
            this->scalAdvFlux_Cons_Form_SecondOrder(a_qFlux, a_q, a_advVel);
            break;
        case RHSParameters::Reconstruction::FourthOrder:
            this->scalAdvFlux_Cons_Form_FourthOrder(a_qFlux, a_q, a_advVel);
            break;
    };

    // Compute divergence
    const bool scaleByJinv = false; // AMRNSLevelRHS will do this.
    m_finiteDiffPtr->levelDivergenceMAC(a_kq, a_qFlux, scaleByJinv);
    // m_finiteDiffPtr->levelDivergenceMAC_4thOrder(a_kq, a_qFlux, scaleByJinv);
    checkForValidNAN(a_kq);
}


// -----------------------------------------------------------------------------
void
AMRNSLevel::scalAdvFlux_Cons_Form_SecondOrder(
    LevelData<FluxBox>&         a_qFlux,
    const LevelData<FArrayBox>& a_q,
    const LevelData<FluxBox>&   a_advVel) const
{
    const GeoSourceInterface& geoSrc = m_levGeoPtr->getGeoSource();
    const RealVect&           dXi    = m_levGeoPtr->getDXi();
    const DisjointBoxLayout&  grids  = a_q.getBoxes();

    for (DataIterator dit(grids); dit.ok(); ++dit) {
        FArrayBox JqFAB(a_q[dit].box(), 1);
        geoSrc.fill_J(JqFAB, 0, dXi);
        JqFAB.mult(a_q[dit]);

        Convert::CellsToAllFaces(a_qFlux[dit], JqFAB);
        m_levGeoPtr->divByJ(a_qFlux[dit], dit());

        a_qFlux[dit].mult(a_advVel[dit], grids[dit], 0, 0, 1);
        a_qFlux[dit].negate();
    }
    checkForValidNAN(a_qFlux);
}


// -----------------------------------------------------------------------------
#if 1
void
AMRNSLevel::scalAdvFlux_Cons_Form_FourthOrder(
    LevelData<FluxBox>&         a_qFlux,
    const LevelData<FArrayBox>& a_q,
    const LevelData<FluxBox>&   a_advVel) const
{
    const DisjointBoxLayout& grids = a_q.getBoxes();

    // Create q that we can modify.
    LevelData<FArrayBox> Jq(grids, 1, IntVect::Unit);
    for (DataIterator dit(grids); dit.ok(); ++dit) {
        Jq[dit].copy(a_q[dit]);
        m_levGeoPtr->multByJ(Jq[dit], dit());
    }

    // Construct 2nd order FC Jq.
    Convert::CellsToAllFaces(a_qFlux, Jq);

    // Upgrade to 4th order...
    // Compute slopes.
    LevelData<FluxBox> deltaJq(grids, 1, IntVect::Unit);
    for (DataIterator dit(grids); dit.ok(); ++dit) {
        for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
            FArrayBox&       deltaJqFAB = deltaJq[dit][fcDir];
            const FArrayBox& JqFAB      = Jq[dit];
            const Box        fcValid    = grids[dit].surroundingNodes(fcDir);
            constexpr Real   dummyDXi   = 1.0;

            deltaJqFAB.setVal(quietNAN);
            FiniteDiff::partialD(
                deltaJqFAB, 0, fcValid, JqFAB, 0, fcDir, dummyDXi);
        }
    }
    constexpr bool extrapOrder = 2;
    BCTools::extrapAllGhosts(deltaJq, extrapOrder);
    deltaJq.exchange();

    // Promote FC Jq to 4th order.
    for (DataIterator dit(grids); dit.ok(); ++dit) {
        for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
            FArrayBox&       JqFAB      = a_qFlux[dit][fcDir];
            const FArrayBox& deltaJqFAB = deltaJq[dit][fcDir];
            const IntVect    e          = BASISV(fcDir);
            constexpr Real   coeff      = 1.0 / 12.0;

            const Box upgradeBox =
                grids[dit].surroundingNodes(fcDir);

            for (BoxIterator bit(upgradeBox); bit.ok(); ++bit) {
                const IntVect& fc = bit();
                JqFAB(fc) -= coeff * (deltaJqFAB(fc + e) - deltaJqFAB(fc - e));
            }
        }
    }

    // Compute FC flux.
    for (DataIterator dit(grids); dit.ok(); ++dit) {
        m_levGeoPtr->divByJ(a_qFlux[dit], dit());
        a_qFlux[dit].mult(a_advVel[dit], grids[dit], 0, 0, 1);
        a_qFlux[dit].negate();
    }
    checkForValidNAN(a_qFlux);
}

#else

void
AMRNSLevel::scalAdvFlux_Cons_Form_FourthOrder(
    LevelData<FluxBox>&         a_qFlux,
    const LevelData<FArrayBox>& a_q,
    const LevelData<FluxBox>&   a_advVel) const
{
    const DisjointBoxLayout& grids = a_q.getBoxes();

    {
        // Prepare advecting velocity
        LevelData<FluxBox> advVelGrow(grids, 1, 2 * IntVect::Unit);
        {
            debugInitLevel(advVelGrow);
            for (DataIterator dit(grids); dit.ok(); ++dit) {
                advVelGrow[dit].copy(a_advVel[dit]);
            }
            BCTools::extrapAllGhosts(advVelGrow, 4, IntVect::Unit);
            FourthOrder::nodalToAvg(advVelGrow);
            advVelGrow.exchange();
        }


        // Prepare advected scalar.
        LevelData<FArrayBox> ccJq(grids, 1, 2 * IntVect::Unit);
        for (DataIterator dit(grids); dit.ok(); ++dit) {
            ccJq[dit].copy(a_q[dit], a_q[dit].box());
            m_levGeoPtr->multByJ(ccJq[dit], dit());
        }
        BCTools::extrapAllGhosts(ccJq, 4, IntVect::Unit);
        ccJq.exchange();
        for (DataIterator dit(grids); dit.ok(); ++dit) {
            FourthOrder::nodalToAvg(ccJq[dit], ccJq[dit].box());
        }

        // Construct 2nd order FC Jq.
        LevelData<FluxBox> fcJq(grids, 1);
        Convert::CellsToAllFaces(fcJq, ccJq); // Redundant!

        for (DataIterator dit(grids); dit.ok(); ++dit) {
            for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
                FArrayBox&       fcJqFAB = fcJq[dit][fcDir];
                const FArrayBox& ccJqFAB = ccJq[dit];
                const Box        fcValid = grids[dit].surroundingNodes(fcDir);
                const IntVect    e       = BASISV(fcDir);

                for (BoxIterator bit(fcValid); bit.ok(); ++bit) {
                    const IntVect& iv = bit();
                    fcJqFAB(iv) =
                        (7.0 / 12.0) * (ccJqFAB(iv) + ccJqFAB(iv - e)) -
                        (1.0 / 12.0) * (ccJqFAB(iv + e) + ccJqFAB(iv - 2 * e));
                }
            }
        }
        fcJq.exchange();

        // Compute FC flux.
        for (DataIterator dit(grids); dit.ok(); ++dit) {
            m_levGeoPtr->divByJ(fcJq[dit], dit());

            // a_qFlux[dit].mult(a_advVel[dit], grids[dit], 0, 0, 1);

            for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
                FArrayBox&       fluxFAB = a_qFlux[dit][fcDir];
                FArrayBox&       qFAB    = fcJq[dit][fcDir];
                const FArrayBox& JuFAB   = a_advVel[dit][fcDir];
                const Box        fcValid = grids[dit].surroundingNodes(fcDir);

                D_TERM(
                FArrayBox JuSlopeFAB(fcValid, SpaceDim - 1);,
                m_finiteDiffPtr->computeSlopes(JuSlopeFAB, 0, JuFAB, 0, fcValid, (fcDir + 1) % SpaceDim, 1.0);,
                m_finiteDiffPtr->computeSlopes(JuSlopeFAB, 1, JuFAB, 0, fcValid, (fcDir + 2) % SpaceDim, 1.0);)

                D_TERM(
                FArrayBox qSlopeFAB(fcValid, SpaceDim - 1);,
                m_finiteDiffPtr->computeSlopes(qSlopeFAB, 0, qFAB, 0, fcValid, (fcDir + 1) % SpaceDim, 1.0);,
                m_finiteDiffPtr->computeSlopes(qSlopeFAB, 1, qFAB, 0, fcValid, (fcDir + 2) % SpaceDim, 1.0);)

                for (BoxIterator bit(fcValid); bit.ok(); ++bit) {
                    const IntVect& fc = bit();

                    fluxFAB(fc) = JuFAB(fc) * qFAB(fc)
                        + (D_TERM(,
                           + JuSlopeFAB(fc,0) * qSlopeFAB(fc,0),
                           + JuSlopeFAB(fc,1) * qSlopeFAB(fc,1))) / 12.0;
                }
            }
            a_qFlux[dit].negate();
        }
        checkForValidNAN(a_qFlux);
    }

    // // Upwind PLM
    // {
    //     // Create q that we can modify.
    //     LevelData<FArrayBox> Jq(grids, 1, IntVect::Unit);
    //     for (DataIterator dit(grids); dit.ok(); ++dit) {
    //         Jq[dit].copy(a_q[dit]);
    //         m_levGeoPtr->multByJ(Jq[dit], dit());
    //     }

    //     // Compute 2nd order slopes
    //     LevelData<FArrayBox> deltaJq2(grids, SpaceDim, IntVect::Unit);
    //     for (DataIterator dit(grids); dit.ok(); ++dit) {
    //         for (int d = 0; d < SpaceDim; ++d) {
    //             m_finiteDiffPtr->computeSlopes(deltaJq2[dit], d, Jq[dit], 0, Jq[dit].box(), d, 1.0);
    //         }
    //     }
    //     BCTools::extrapAllGhosts(deltaJq2, 2);
    //     deltaJq2.exchange();

    //     // Compute 4th order slopes
    //     LevelData<FArrayBox> deltaJq4(grids, SpaceDim, IntVect::Unit);
    //     for (DataIterator dit(grids); dit.ok(); ++dit) {
    //         const Box region = grids[dit];
    //         for (int d = 0; d < SpaceDim; ++d) {
    //             FArrayBox        d4FAB(Interval(d, d), deltaJq4[dit]);
    //             const FArrayBox  d2FAB(Interval(d, d), deltaJq2[dit]);
    //             const FArrayBox& qFAB = Jq[dit];
    //             const IntVect    e     = BASISV(d);

    //             for (BoxIterator bit(region); bit.ok(); ++bit) {
    //                 const IntVect& cc = bit();
    //                 const IntVect  ccl = cc - e;
    //                 const IntVect  ccr = cc + e;

    //                 d4FAB(cc) = (2.0/3.0) * (
    //                     (qFAB(ccr) - 0.25*d2FAB(ccr)) -
    //                     (qFAB(ccl) + 0.25*d2FAB(ccl))
    //                 );
    //             }
    //         }
    //     }
    //     BCTools::extrapAllGhosts(deltaJq4, 2);
    //     deltaJq4.exchange();

    //     // Left and right states
    //     LevelData<FArrayBox> qL(grids, SpaceDim, IntVect::Unit);
    //     LevelData<FArrayBox> qR(grids, SpaceDim, IntVect::Unit);
    //     for (DataIterator dit(grids); dit.ok(); ++dit) {
    //         const FArrayBox& dqFAB = deltaJq4[dit];

    //         for (int d = 0; d < SpaceDim; ++d) {
    //             qL[dit].copy(Jq[dit], 0, d, 1);
    //             qL[dit].plus(dqFAB, -0.5, d, d, 1);

    //             qR[dit].copy(Jq[dit], 0, d, 1);
    //             qR[dit].plus(dqFAB, 0.5, d, d, 1);
    //         }
    //         m_levGeoPtr->divByJ(qL[dit], dit());
    //         m_levGeoPtr->divByJ(qR[dit], dit());
    //     }
    //     BCTools::extrapAllGhosts(qL, 2);
    //     qL.exchange();
    //     BCTools::extrapAllGhosts(qR, 2);
    //     qR.exchange();

    //     // Compute fluxes
    //     for (DataIterator dit(grids); dit.ok(); ++dit) {
    //         for (int d = 0; d < SpaceDim; ++d) {
    //             FArrayBox&       fluxFAB = a_qFlux[dit][d];
    //             const FArrayBox  qLFAB(Interval(d, d), qL[dit]);
    //             const FArrayBox  qRFAB(Interval(d, d), qR[dit]);
    //             const FArrayBox& JuFAB      = a_advVel[dit][d];
    //             const Box        fluxRegion = surroundingNodes(grids[dit], d);
    //             const IntVect    e          = BASISV(d);

    //             for (BoxIterator bit(fluxRegion); bit.ok(); ++bit) {
    //                 const IntVect  fc    = bit();
    //                 const IntVect  ccl   = fc - e;
    //                 const IntVect& ccr   = fc;

    //                 const Real vel = JuFAB(fc);
    //                 const Real qL  = qRFAB(ccl);
    //                 const Real qR  = qLFAB(ccr);

    //                 constexpr Real small = 1.0e-8;

    //                 if (vel > small) {
    //                     fluxFAB(fc) = -vel * qL;
    //                 } else if (vel < -small) {
    //                     fluxFAB(fc) = -vel * qR;
    //                 } else {
    //                     fluxFAB(fc) = -vel * 0.5 * (qL + qR);
    //                 }
    //             }
    //         }
    //     }
    // }

    // // Upwind PPM
    // {
    //     // Create q that we can modify.
    //     LevelData<FArrayBox> Jq(grids, 1, IntVect::Unit);
    //     for (DataIterator dit(grids); dit.ok(); ++dit) {
    //         Jq[dit].copy(a_q[dit]);
    //         m_levGeoPtr->multByJ(Jq[dit], dit());
    //     }

    //     // Compute 2nd order slopes
    //     LevelData<FArrayBox> deltaJq2(grids, SpaceDim, IntVect::Unit);
    //     for (DataIterator dit(grids); dit.ok(); ++dit) {
    //         for (int d = 0; d < SpaceDim; ++d) {
    //             m_finiteDiffPtr->computeSlopes(deltaJq2[dit], d, Jq[dit], 0, Jq[dit].box(), d, 1.0);
    //         }
    //     }
    //     BCTools::extrapAllGhosts(deltaJq2, 2);
    //     deltaJq2.exchange();

    //     // Compute 4th order slopes
    //     LevelData<FArrayBox> deltaJq4(grids, SpaceDim, IntVect::Unit);
    //     for (DataIterator dit(grids); dit.ok(); ++dit) {
    //         const Box region = grids[dit];
    //         for (int d = 0; d < SpaceDim; ++d) {
    //             FArrayBox        d4FAB(Interval(d, d), deltaJq4[dit]);
    //             const FArrayBox  d2FAB(Interval(d, d), deltaJq2[dit]);
    //             const FArrayBox& qFAB = Jq[dit];
    //             const IntVect    e     = BASISV(d);

    //             for (BoxIterator bit(region); bit.ok(); ++bit) {
    //                 const IntVect& cc = bit();
    //                 const IntVect  ccl = cc - e;
    //                 const IntVect  ccr = cc + e;

    //                 d4FAB(cc) = (2.0/3.0) * (
    //                     (qFAB(ccr) - 0.25*d2FAB(ccr)) -
    //                     (qFAB(ccl) + 0.25*d2FAB(ccl))
    //                 );
    //             }
    //         }
    //     }
    //     BCTools::extrapAllGhosts(deltaJq4, 2);
    //     deltaJq4.exchange();

    //     // Left and right states
    //     LevelData<FArrayBox> qL(grids, SpaceDim, IntVect::Unit);
    //     LevelData<FArrayBox> qR(grids, SpaceDim, IntVect::Unit);
    //     for (DataIterator dit(grids); dit.ok(); ++dit) {
    //         const Box region = grids[dit];

    //         for (int d = 0; d < SpaceDim; ++d) {
    //             FArrayBox        qLFAB(Interval(d, d), qL[dit]);
    //             FArrayBox        qRFAB(Interval(d, d), qR[dit]);
    //             const FArrayBox  dqFAB(Interval(d, d), deltaJq4[dit]);
    //             const FArrayBox& qFAB = Jq[dit];
    //             const IntVect    e    = BASISV(d);

    //             for (BoxIterator bit(region); bit.ok(); ++bit) {
    //                 const IntVect& cc = bit();
    //                 const IntVect  ccl = cc - e;
    //                 const IntVect  ccr = cc + e;

    //                 qRFAB(cc) = 0.5 * (qFAB(cc) + qFAB(ccr)) + (1.0/6.0) * (dqFAB(cc) - dqFAB(ccr));
    //                 qLFAB(cc) = 0.5 * (qFAB(cc) + qFAB(ccl)) - (1.0/6.0) * (dqFAB(cc) - dqFAB(ccl));
    //             }
    //         }
    //         m_levGeoPtr->divByJ(qL[dit], dit());
    //         m_levGeoPtr->divByJ(qR[dit], dit());
    //     }
    //     BCTools::extrapAllGhosts(qL, 2);
    //     qL.exchange();
    //     BCTools::extrapAllGhosts(qR, 2);
    //     qR.exchange();

    //     // PPM limiter
    //     for (DataIterator dit(grids); dit.ok(); ++dit) {
    //         const Box region = grids[dit];

    //         for (int d = 0; d < SpaceDim; ++d) {
    //             FArrayBox        qLFAB(Interval(d, d), qL[dit]);
    //             FArrayBox        qRFAB(Interval(d, d), qR[dit]);
    //             const FArrayBox  dqFAB(Interval(d, d), deltaJq4[dit]);
    //             const FArrayBox& qFAB = Jq[dit];

    //             for (BoxIterator bit(region); bit.ok(); ++bit) {
    //                 const IntVect& cc = bit();

    //                 Real alphaR = qRFAB(cc) - qFAB(cc);
    //                 Real alphaL = qLFAB(cc) - qFAB(cc);

    //                 if (alphaL * alphaR < 0.0) {
    //                     const Real s = (alphaR > alphaL) ? 1.0 : -1.0;
    //                     if (alphaR*alphaR > alphaL * alphaL) {
    //                         alphaR = s * std::min(s * alphaR, -2.0 * s * alphaL);
    //                     } else {
    //                         alphaL = s * std::min(s * alphaL, -2.0 * s * alphaR);
    //                     }
    //                 } else {
    //                     alphaL = 0.0;
    //                     alphaR = 0.0;
    //                 }

    //                 qRFAB(cc) = qFAB(cc) + alphaR;
    //                 qLFAB(cc) = qFAB(cc) + alphaL;
    //             }
    //         }
    //     }
    //     BCTools::extrapAllGhosts(qL, 2);
    //     qL.exchange();
    //     BCTools::extrapAllGhosts(qR, 2);
    //     qR.exchange();

    //     // Compute fluxes
    //     for (DataIterator dit(grids); dit.ok(); ++dit) {
    //         for (int d = 0; d < SpaceDim; ++d) {
    //             FArrayBox&       fluxFAB = a_qFlux[dit][d];
    //             const FArrayBox  qLFAB(Interval(d, d), qL[dit]);
    //             const FArrayBox  qRFAB(Interval(d, d), qR[dit]);
    //             const FArrayBox& JuFAB      = a_advVel[dit][d];
    //             const Box        fluxRegion = surroundingNodes(grids[dit], d);
    //             const IntVect    e          = BASISV(d);

    //             for (BoxIterator bit(fluxRegion); bit.ok(); ++bit) {
    //                 const IntVect  fc    = bit();
    //                 const IntVect  ccl   = fc - e;
    //                 const IntVect& ccr   = fc;

    //                 const Real vel = JuFAB(fc);
    //                 const Real qL  = qRFAB(ccl);
    //                 const Real qR  = qLFAB(ccr);

    //                 constexpr Real small = 1.0e-8;

    //                 if (vel > small) {
    //                     fluxFAB(fc) = -vel * qL;
    //                 } else if (vel < -small) {
    //                     fluxFAB(fc) = -vel * qR;
    //                 } else {
    //                     fluxFAB(fc) = -vel * 0.5 * (qL + qR);
    //                 }
    //             }
    //         }
    //     }
    // }

    checkForValidNAN(a_qFlux);
}
#endif
