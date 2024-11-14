/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2018
 *    Jefferson (Philadelphia University + Thomas Jefferson University) and
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
 *  https://github.com/somarhub.
 ******************************************************************************/
#include "AMRNSLevel.H"
#include "AMRNSLevelF_F.H"
#include "Convert.H"


// AMRNSLevel::rateOfStrain(StaggeredFluxLD&          a_Sia,
//                          const LevelData<FluxBox>& a_cartVel,
//                          const bool                a_multByJ,
//                          const bool                a_makeCovariant,
//                          const Real                a_primaryScale,
//                          const Real                a_transposeScale) const
// {
//     CH_assert(a_Sia    .getBoxes() == m_levGeoPtr->getBoxes());
//     CH_assert(a_cartVel.getBoxes() == m_levGeoPtr->getBoxes());

//     const RealVect&          dXi   = m_levGeoPtr->getDXi();
//     const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();
//     DataIterator             dit   = grids.dataIterator();

//     // Compute J*grad[u] components. Last index will be Cartesian-based.
//     m_finiteDiffPtr->levelVectorGradient(a_Sia, a_cartVel);

//     // Convert J*grad[u] into S^{ia}.
//     for (dit.reset(); dit.ok(); ++dit) {
//         for (int a = 0; a < SpaceDim; ++a) {
//             // Easy, diagonal elements first.
//             {
//                 FArrayBox& SaaFAB = a_Sia[a][a][dit];
//                 SaaFAB *= (a_primaryScale + a_transposeScale);
//                 checkForNAN(SaaFAB, SaaFAB.box());

//                 if (!a_multByJ) {
//                     m_levGeoPtr->divByJ(SaaFAB, dit());
//                 }

//                 if (a_makeCovariant) {
//                     const Box& ccBox = SaaFAB.box();
//                     FArrayBox gdnFAB(ccBox, 1);
//                     m_levGeoPtr->getGeoSource().fill_gdn(gdnFAB, 0, a, dXi);
//                     SaaFAB.mult(gdnFAB, ccBox, 0, 0, 1);
//                 }
//             }

//             // Off-diagonal comps need to be symmetrized.
//             for (int i = a + 1; i < SpaceDim; ++i) {
//                 FArrayBox& SiaFAB = a_Sia[i][a][dit];
//                 FArrayBox& SaiFAB = a_Sia[a][i][dit];
//                 const Box& ecBox  = SiaFAB.box();

//                 // Send last index to mapped basis.
//                 FArrayBox dxdXiFAB(ecBox, 2);
//                 m_levGeoPtr->getGeoSource().fill_dxdXi(dxdXiFAB, 0, a, dXi);
//                 m_levGeoPtr->getGeoSource().fill_dxdXi(dxdXiFAB, 1, i, dXi);
//                 SiaFAB.divide(dxdXiFAB, 0, 0, 1);
//                 SaiFAB.divide(dxdXiFAB, 1, 0, 1);

//                 // Symmetrize
//                 FArrayBox tmpFAB(ecBox, 1);
//                 tmpFAB.copy(SiaFAB);

//                 SiaFAB *= a_primaryScale;
//                 SiaFAB.plus(SaiFAB, a_transposeScale);

//                 SaiFAB *= a_primaryScale;
//                 SaiFAB.plus(tmpFAB, a_transposeScale);

//                 // Send last index back to Cartesian basis.
//                 SiaFAB.mult(dxdXiFAB, 0, 0, 1);
//                 SaiFAB.mult(dxdXiFAB, 1, 0, 1);

//                 checkForNAN(SiaFAB, ecBox);
//                 checkForNAN(SaiFAB, ecBox);

//                 if (!a_multByJ) {
//                     m_levGeoPtr->divByJ(SiaFAB, dit());
//                     m_levGeoPtr->divByJ(SaiFAB, dit());
//                 }

//                 if (a_makeCovariant) {
//                     FArrayBox gdnFAB(ecBox, 1);
//                     m_levGeoPtr->getGeoSource().fill_gdn(gdnFAB, 0, a, dXi);
//                     SiaFAB.mult(gdnFAB, ecBox, 0, 0, 1);
//                     SaiFAB.mult(gdnFAB, ecBox, 0, 0, 1);
//                 }
//             } // i
//         } // a
//     } // dit
// }


// // -----------------------------------------------------------------------------
// void
// AMRNSLevel::viscousStress(StaggeredFluxLD&            a_Tia,
//                           const StaggeredFluxLD&      a_JSia,
//                           const RealVect&             a_nu,
//                           const LevelData<FArrayBox>& a_eddyNu) const
// {
//     CH_assert(a_Tia   .getBoxes() == m_levGeoPtr->getBoxes());
//     CH_assert(a_JSia  .getBoxes() == m_levGeoPtr->getBoxes());
//     CH_assert(a_eddyNu.getBoxes() == m_levGeoPtr->getBoxes());

//     const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();
//     DataIterator             dit   = grids.dataIterator();

//     // Copy T^{ia} = J*S^{ia}, if necessary.
//     if (&a_Tia != &a_JSia) {
//         for (int a = 0; a < SpaceDim; ++a) {
//             for (int i = 0; i < SpaceDim; ++i) {
//                 for (dit.reset(); dit.ok(); ++dit) {
//                     a_Tia[i][a][dit].copy(a_JSia[i][a][dit]);
//                 }
//             }
//         }
//     }

//     // Convert T^{ia} = J*S^{ia} into T^{ia} = (nu + eddyNu) * 2*J*S^{ia}.
//     for (dit.reset(); dit.ok(); ++dit) {
//         for (int a = 0; a < SpaceDim; ++a) {
//             for (int i = 0; i < SpaceDim; ++i) {
//                 FArrayBox&       TiaFAB    = a_Tia[i][a][dit];
//                 const Real       nuDir     = a_nu[i];
//                 const FArrayBox& eddyNuFAB = a_eddyNu[dit];
//                 const Box&       region    = TiaFAB.box();

//                 if (a == i) {
//                     FABAlgebra::CCmultCC(
//                         TiaFAB, 0, region, eddyNuFAB, 0, nuDir);
//                 } else {
//                     FABAlgebra::ECmultCC(
//                         TiaFAB, 0, region, eddyNuFAB, 0, nuDir);
//                 }

//                 TiaFAB *= 2.0;

//                 checkForNAN(TiaFAB, region);
//             }
//         } // a
//     } // dit
// }

// -----------------------------------------------------------------------------
void
AMRNSLevel::computeEddyNu(LevelData<FArrayBox>&     a_eddyNu,
                          const LevelData<FluxBox>& a_cartVel,
                          const Real                a_time) const
{
    const auto* ctx            = ProblemContext::getInstance();
    const int   eddyViscMethod = ctx->rhs.eddyViscMethod[m_level];

    if (eddyViscMethod <= 0) {
        setValLevel(a_eddyNu, 0.0);
        return;

    } else if (eddyViscMethod == RHSParameters::EddyViscMethods::DUCROS) {
        const int numFilterSweeps = 3;
        const RealVect& eddyScale = ctx->rhs.eddyScale;
        this->SGSModel_Ducros(
            a_eddyNu, a_cartVel, a_time, numFilterSweeps, eddyScale);

    } else if (eddyViscMethod == RHSParameters::EddyViscMethods::STATIC_SMAG) {
        UNDEFINED_FUNCTION();

    } else if (eddyViscMethod == RHSParameters::EddyViscMethods::DYNAMIC_SMAG) {
        UNDEFINED_FUNCTION();

    } else {
        MAYDAYERROR("eddyViscMethod = " << eddyViscMethod << " not recognized.");
    }

    // Set all BCs on eddy viscosity.
    // For now, just extrapolate and exchange.
    // Someday, maybe interp at CFI if needed.
    BCTools::extrapAllGhosts(a_eddyNu, 2);
    a_eddyNu.exchange(m_statePtr->qExCopier);
    a_eddyNu.exchange(m_statePtr->qExCornerCopier);
}


// -----------------------------------------------------------------------------
void
AMRNSLevel::LaplacianFilter(LevelData<FArrayBox>& a_ccCartVel,
                            const int             a_numFilterSweeps,
                            const RealVect&       a_dirScale) const
{
    CH_assert(a_ccCartVel.getBoxes() == this->getBoxes());
    CH_assert(a_ccCartVel.ghostVect() >= IntVect::Unit);
    CH_assert(a_numFilterSweeps > 0);

    const DisjointBoxLayout& grids = this->getBoxes();

    LevelData<FArrayBox>* filtVelPtr = &a_ccCartVel;
    LevelData<FArrayBox> _ccVel2(grids, SpaceDim, IntVect::Unit);
    LevelData<FArrayBox>* origVelPtr = &_ccVel2;

    // At this point,
    //   *origVelPtr = bogus data.
    //   *filtVelPtr = velocity we want to filter.
    // But this will soon change.
    //
    // At the beginning of each filter sweep, there is a call to std::swap.
    // After std::swap,
    //   *origVelPtr = velocity we want to filter.
    //   *filtVelPtr = bogus data / to be filled with the filtered velocity,
    // and this is how you should think about these pointers.
    //
    // As we loop and begin each new sweep, *origVelPtr will be swapped with
    // *filtVelPtr, so that *origVelPtr = velocity we want to filter.
    //
    // When this block is complete, *filtVelPtr will contain the filtered
    // velocity.
    for (int iter = 0; iter < a_numFilterSweeps; ++iter) {
        // Put filtered velocity into *origVelPtr.
        // *filtVelPtr will be overwritten in this iter.
        std::swap(filtVelPtr, origVelPtr);
        debugInitLevel(*filtVelPtr);

        // Filter once, filling *filtVelPtr.
        for (DataIterator dit(grids); dit.ok(); ++dit) {
            FArrayBox&       filtFAB = (*filtVelPtr)[dit];
            const FArrayBox& origFAB = (*origVelPtr)[dit];
            const Box&       valid   = grids[dit];

            FORT_FILTER_LAPLACIAN (
                CHF_FRA(filtFAB),
                CHF_CONST_FRA(origFAB),
                CHF_BOX(valid),
                CHF_CONST_REALVECT(a_dirScale));
        } // end loop over grids (dit)

        BCTools::extrapAllGhosts(*filtVelPtr, 2);
        filtVelPtr->exchange();
    } // end loop over filter iters (iter)

    // Copy back to user's holder, if needed.
    if (filtVelPtr != &a_ccCartVel) {
        for (DataIterator dit(grids); dit.ok(); ++dit) {
            FArrayBox&       destFAB = a_ccCartVel[dit];
            const FArrayBox& srcFAB  = (*filtVelPtr)[dit];

            destFAB.copy(srcFAB);
        }
    }
}


// -----------------------------------------------------------------------------
void
AMRNSLevel::SGSModel_Ducros(LevelData<FArrayBox>&     a_nuT,
                            const LevelData<FluxBox>& a_cartVel,
                            const Real                /*a_time*/,
                            const int                 a_numFilterSweeps,
                            const RealVect&           a_dirScale) const
{
    // Sanity checks
    CH_assert(a_nuT.nComp() == 1);
    CH_assert(a_nuT.getBoxes() == this->getBoxes());
    CH_assert(a_cartVel.nComp() == 1);
    CH_assert(a_cartVel.getBoxes() == this->getBoxes());
    CH_assert(a_numFilterSweeps > 0);

    // Gather data structures.
    const RealVect&          dXi   = m_levGeoPtr->getDXi();
    const DisjointBoxLayout& grids = this->getBoxes();

    // Create CC vel.
    LevelData<FArrayBox> ccVel(grids, SpaceDim, IntVect::Unit);
    for (DataIterator dit(grids); dit.ok(); ++dit) {
        for (int dir = 0; dir < SpaceDim; ++dir) {
            FArrayBox fcVelCompFAB(a_cartVel[dit][dir].box(), 1);
            fcVelCompFAB.copy(a_cartVel[dit][dir]);

            m_levGeoPtr->multByJ(fcVelCompFAB, dit());
            Convert::Simple(ccVel[dit], dir, grids[dit], fcVelCompFAB, 0);
            m_levGeoPtr->divByJ(ccVel[dit], dit(), dir);
        }
    }
    BCTools::extrapAllGhosts(ccVel, 2);  // Properly setting BCs increases nuT at the CFI.
    ccVel.exchange();

    // Filter the CC velocity.
    this->LaplacianFilter(ccVel, a_numFilterSweeps, a_dirScale);

    // Compute the eddy viscosity.
    for (DataIterator dit(grids); dit.ok(); ++dit) {
        FArrayBox&       nuTFAB  = a_nuT[dit];
        const FArrayBox& velFAB  = ccVel[dit];
        const FArrayBox& ccJFAB  = m_levGeoPtr->getCCJ()[dit];
        const Box&       ccValid = grids[dit];

        // In 2D, only vel comps 0 and 1 will be used.
        FORT_SGSMODEL_DUCROS (
            CHF_FRA1(nuTFAB, 0),
            CHF_CONST_FRA1(velFAB, 0),
            CHF_CONST_FRA1(velFAB, 1),
            CHF_CONST_FRA1(velFAB, SpaceDim - 1),
            CHF_CONST_FRA1(ccJFAB, 0),
            CHF_BOX(ccValid),
            CHF_CONST_REALVECT(dXi),
            CHF_CONST_REALVECT(a_dirScale));
    }

    BCTools::extrapAllGhosts(a_nuT, 2);
    a_nuT.exchange();
}

