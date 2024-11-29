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
#include "Integral.H"
#include "IntegralF_F.H"
#include "Masks.H"
#include "Comm.H"


// // -----------------------------------------------------------------------------
// // Returns the integral of phi over the valid region.
// // This version does not require a full LevelGeometry.
// // This function can only handle cell-centered data.
// // -----------------------------------------------------------------------------
// Real
// Integral::sum(Real&                       a_vol,
//               const LevelData<FArrayBox>& a_phi,
//               const DisjointBoxLayout*    a_finerGridsPtr,
//               const IntVect&              a_fineRefRatio,
//               const RealVect&             a_dXi,
//               const LevelData<FArrayBox>& a_CCJinv,
//               const int                   a_comp)
// {
//     // Sanity checks
//     CH_assert(0 <= a_comp);
//     CH_assert(a_comp < a_phi.nComp());
//     CH_assert(a_CCJinv.getBoxes().compatible(a_phi.getBoxes()));

//     // Gather geometric structures
//     const Real               dXiScale = a_dXi.product();
//     const DisjointBoxLayout& grids    = a_phi.getBoxes();
//     DataIterator             dit      = grids.dataIterator();

//     // Loop over this level's grids, adding to sum.
//     Real localSum = 0.0;
//     Real localVol = 0.0;
//     for (dit.reset(); dit.ok(); ++dit) {
//         // Create references for convenience / allocate calucation space.
//         const FArrayBox& phiFAB  = a_phi[dit];
//         const FArrayBox& JinvFAB = a_CCJinv[dit];
//         const Box&       valid   = grids[dit];

//         // If needed, we can adopt other centerings later.
//         CH_assert(phiFAB.box().type() == IntVect::Zero);

//         // Phi should be equal to or larger than the valid region.
//         CH_assert(phiFAB.box().contains(valid));

//         // Create a grid mask. 0 = covered by finer grid, 1 = valid data.
//         BaseFab<int> maskFAB(valid, 1);
//         Masks::zeroInvalid(maskFAB, a_finerGridsPtr, a_fineRefRatio);

//         // Add to localSum
//         FORT_INTEGRAL_SUMWITHJINV(CHF_REAL(localSum),
//                                   CHF_REAL(localVol),
//                                   CHF_CONST_FRA1(phiFAB, a_comp),
//                                   CHF_CONST_FRA1(JinvFAB, 0),
//                                   CHF_CONST_FIA1(maskFAB, 0),
//                                   CHF_BOX(valid),
//                                   CHF_CONST_REAL(dXiScale));
//     }  // end loop over this level's grids

//     // Compute global sum (this is where the MPI communication happens)
// #ifdef CH_MPI
//     Real globalSum = 0.0;
//     int result = MPI_Allreduce(
//         &localSum, &globalSum, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

//     if (result != MPI_SUCCESS) {
//         MayDay::Error(
//             "Sorry, but I had a communication error in computeMappedSum");
//     }

//     Real globalVol = 0.0;
//     result = MPI_Allreduce(
//         &localVol, &globalVol, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

//     if (result != MPI_SUCCESS) {
//         MayDay::Error(
//             "Sorry, but I had a communication error in computeMappedSum");
//     }

// #else
//     Real globalSum = localSum;
//     Real globalVol = localVol;
// #endif

//     a_vol += globalVol;
//     return globalSum;
// }


// -----------------------------------------------------------------------------
// Returns the integral of phi over the valid region.
// This function can only handle cell-centered data.
// -----------------------------------------------------------------------------
Real
Integral::sum(Real&                       a_vol,
              const LevelData<FArrayBox>& a_phi,
              const LevelGeometry&        a_levGeo,
              const bool                  /*a_sumJPhi*/,
              const int                   a_comp)
{
    // Assertions are done in this function call.
    Real localVol = 0.0;
    Real localSum = localMappedSum(localVol,
                                   a_phi,
                                   a_levGeo.getDXi(),
                                   a_levGeo.getFineRefRatio(),
                                   a_levGeo.getFineGridsPtr(),
                                   a_levGeo.getCCJ(),
                                   a_comp);

    Comm::reduce(localSum, MPI_SUM);
    Comm::reduce(localVol, MPI_SUM);

    a_vol = localVol;
    return localSum;
}


// -----------------------------------------------------------------------------
// Returns the integral of phi over the valid region.
// a_levGeo can be any levGeo in the hierarchy.
// If a_sumJPhi == false, the we will sum a_phi directly, but still use J to
// compute cell volumes.
// This function can only handle cell-centered data.
// -----------------------------------------------------------------------------
Real
Integral::sum(Real&                                a_vol,
              const Vector<LevelData<FArrayBox>*>& a_phi,
              const LevelGeometry&                 a_levGeo,
              const bool                           a_sumJPhi,
              const int                            a_comp,
              const int                            a_lBase)
{
    // Sanity check on a_lBase
    const int vectorSize = a_phi.size();
    CH_assert(0 <= a_lBase);
    CH_assert(a_lBase < vectorSize);

    // We need at least one defined level.
    CH_assert(a_phi[a_lBase] != NULL);

    // Find highest defined level.
    int lev;
    for (lev = a_lBase; lev < vectorSize; ++lev) {
        if (a_phi[lev] == NULL) break;
        if (!a_phi[lev]->isDefined()) break;
    }
    const int topLevel = lev - 1;

    // Get every level's levGeo.
    Vector<const LevelGeometry*> vLevGeos = a_levGeo.getAMRLevGeos();

    // Loop over levels and add to total sum.
    Real localSum = 0.0;
    Real localVol = 0.0;
    for (lev = a_lBase; lev <= topLevel; ++lev) {
        const LevelData<FArrayBox>& levelPhi = *a_phi[lev];

        const DisjointBoxLayout* finerGridsPtr = NULL;
        if (lev < topLevel) {
            finerGridsPtr = &(a_phi[lev + 1]->getBoxes());
        }

        const LevelGeometry& levGeoRef = *vLevGeos[lev];
        CH_assert(levGeoRef.getBoxes() == levelPhi.getBoxes());  // Overkill?

        // Add sum to running total.
        if (a_sumJPhi) {
            localSum += localMappedSum(localVol,
                                       levelPhi,
                                       levGeoRef.getDXi(),
                                       levGeoRef.getFineRefRatio(),
                                       finerGridsPtr,
                                       levGeoRef.getCCJ(),
                                       a_comp);
        } else {
            localSum += localUnmappedSum(localVol,
                                         levelPhi,
                                         levGeoRef.getDXi(),
                                         levGeoRef.getFineRefRatio(),
                                         finerGridsPtr,
                                         &(levGeoRef.getCCJ()),
                                         a_comp);
        }
    }

    Comm::reduce(localSum, MPI_SUM);
    Comm::reduce(localVol, MPI_SUM);

    a_vol = localVol;
    return localSum;
}


// -----------------------------------------------------------------------------
// Returns the integral of Jphi over the valid region.
// This version does not require a levGeo.
// This function can only handle cell-centered data.
// -----------------------------------------------------------------------------
Real
Integral::sum(const LevelData<FArrayBox>& a_Jphi,
              const RealVect&             a_dXi,
              const int                   a_comp)
{
    // Assertions are done in this function call.
    Real localVol = 0.0;
    Real localSum = localUnmappedSum(localVol,
                                     a_Jphi,
                                     a_dXi,
                                     IntVect::Unit,
                                     NULL,
                                     NULL,
                                     a_comp);

    Comm::reduce(localSum, MPI_SUM);

    return localSum;
}


// -----------------------------------------------------------------------------
// Same, but also returns the volume of the integrated region.
// If a_Jptr == nullptr, then we just use ones.
// -----------------------------------------------------------------------------
void
Integral::sum(Real&                       a_sum,
              Real&                       a_vol,
              const LevelData<FArrayBox>& a_phi,
              const LevelData<FArrayBox>* a_Jptr,
              const RealVect&             a_dXi,
              const int                   a_comp)
{
    a_vol = 0.0;

    if (a_Jptr) {
        a_sum = localMappedSum(
            a_vol, a_phi, a_dXi, IntVect::Unit, nullptr, *a_Jptr, a_comp);
    } else {
        a_sum = localUnmappedSum(
            a_vol, a_phi, a_dXi, IntVect::Unit, nullptr, nullptr, a_comp);
    }

    Comm::reduce(a_sum, MPI_SUM);
    Comm::reduce(a_vol, MPI_SUM);
}


// -----------------------------------------------------------------------------
// Returns the integral of (v.n dS) over the level's boundary,
// where n is the unit outward normal.
// It is assumed v is in the curvilinear basis and scaled by J.
// -----------------------------------------------------------------------------
Real
Integral::bdrySum(const LevelData<FluxBox>& a_v,
                  const RealVect&           a_dXi,
                  PhysBdryIter&             a_physBdryIter,
                  CFIIter&                  a_cfiIter)
{
    Real sum = 0.0;

    Lockable::Key key = a_physBdryIter.lock();
    for (a_physBdryIter.reset(key); a_physBdryIter.ok(); a_physBdryIter.next(key)) {
        const DataIndex&     di         = a_physBdryIter->di;
        const Real           rsign      = a_physBdryIter->rsign;
        const int            bdryDir    = a_physBdryIter->dir;
        const Box&           fcBdryBox  = a_physBdryIter->fcBdryBox;

        const FArrayBox& vFAB = a_v[di][bdryDir];
        const Real       dA   = a_dXi.product() / a_dXi[bdryDir];

        sum += rsign * vFAB.sum(fcBdryBox, 0, 1) * dA;
    }
    a_physBdryIter.unlock(key);

    key = a_cfiIter.lock();
    for(a_cfiIter.reset(key); a_cfiIter.ok(); a_cfiIter.next(key)) {
        const DataIndex& di      = a_cfiIter->di;
        const int        isign   = a_cfiIter->isign;
        const Real       rsign   = a_cfiIter->rsign;
        const int        bdryDir = a_cfiIter->dir;
        const CFIVS&     cfivs   = a_cfiIter->cfivs;

        const IntVect    eb   = ((isign > 0) ? IntVect::Zero : BASISV(bdryDir));
        const FArrayBox& vFAB = a_v[di][bdryDir];
        const Real       dA   = a_dXi.product() / a_dXi[bdryDir];

        if (cfivs.isPacked()) {
            Box cfiBox = cfivs.packedBox();
            cfiBox.shiftHalf(bdryDir, -isign);
            sum += rsign * vFAB.sum(cfiBox, 0, 1) * dA;

        } else {
            Real localSum = 0.0;

            const IntVectSet& ivs = cfivs.getIVS();
            IVSIterator ivsit(ivs);
            for (ivsit.reset(); ivsit.ok(); ++ivsit) {
                localSum += vFAB(ivsit() + eb);
            }

            sum += rsign * localSum * dA;
        }
    }
    a_cfiIter.unlock(key);

    Comm::reduce(sum, MPI_SUM);

    return sum;
}


// -----------------------------------------------------------------------------
// This version creates its own iterators.
// -----------------------------------------------------------------------------
Real
Integral::bdrySum(const LevelData<FluxBox>& a_v,
                  const LevelGeometry&      a_levGeo)
{
    const DisjointBoxLayout& grids    = a_levGeo.getBoxes();
    const CFRegion&          cfRegion = a_levGeo.getCFRegion();
    const RealVect&          dXi      = a_levGeo.getDXi();

    PhysBdryIter physBdryIter(grids);
    CFIIter cfiIter(grids, cfRegion);

    return Integral::bdrySum(a_v, dXi, physBdryIter, cfiIter);
}


// -----------------------------------------------------------------------------
// Performs most of the computation for the sum functions. This version does
// not perform MPI communication and only handles cell-centered data.
// This will scale both phi and dXi by J.
// -----------------------------------------------------------------------------
Real
Integral::localMappedSum (Real&                       a_vol,
                          const LevelData<FArrayBox>& a_phi,
                          const RealVect&             a_dXi,
                          const IntVect&              a_refRatio,
                          const DisjointBoxLayout*    a_finerGridsPtr,
                          const LevelData<FArrayBox>& a_ccJ,
                          const int                   a_comp)
{
    // Gather geometric structures
    const Real dXiScale = a_dXi.product();
    const DisjointBoxLayout& grids = a_phi.getBoxes();
    DataIterator dit = grids.dataIterator();

    CH_assert(0 <= a_comp);
    CH_assert(a_comp < a_phi.nComp());
    CH_assert(dXiScale > 0.0);
    CH_assert(a_ccJ.getBoxes().compatible(grids));
    CH_assert(a_ccJ.nComp() == 1);

    Real localSum = 0.0;

    // Loop over this level's grids, adding to sum.
    for (dit.reset(); dit.ok(); ++dit) {
        // Create references for convenience / allocate calucation space.
        const FArrayBox& phiFAB = a_phi[dit];
        const FArrayBox& JFAB = a_ccJ[dit];
        const Box& valid = grids[dit];

        // If needed, we can adopt other centerings later.
        CH_assert(phiFAB.box().type() == IntVect::Zero);

        // Phi should be equal to or larger than the valid region.
        CH_assert(phiFAB.box().contains(valid));

        // Create a grid mask. 0 = covered by finer grid, 1 = valid data.
        BaseFab<int> maskFAB(valid, 1);
        Masks::zeroInvalid(maskFAB, a_finerGridsPtr, a_refRatio);

        // Add to localSum
        FORT_INTEGRAL_MAPPEDSUM(
            CHF_REAL(localSum),
            CHF_REAL(a_vol),
            CHF_CONST_FRA1(phiFAB,a_comp),
            CHF_CONST_FRA1(JFAB,0),
            CHF_CONST_FIA1(maskFAB,0),
            CHF_BOX(valid),
            CHF_CONST_REAL(dXiScale));
    } // end loop over this level's grids
    return localSum;
}


// -----------------------------------------------------------------------------
// Performs most of the computation for the sum functions. This version does
// not perform MPI communication and only handles cell-centered data.
// If a_ccJPtr is not NULL:
//   a_phi will not be scaled by J.
//   a_dXi _will_ be scaled by J.      <----- The only place J is used.
// If a_ccJPtr is NULL:
//   a_phi will not be scaled by J.
//   a_dXi will not be scaled by J.
// -----------------------------------------------------------------------------
Real
Integral::localUnmappedSum (Real&                       a_vol,
                            const LevelData<FArrayBox>& a_phi,
                            const RealVect&             a_dXi,
                            const IntVect&              a_refRatio,
                            const DisjointBoxLayout*    a_finerGridsPtr,
                            const LevelData<FArrayBox>* a_ccJPtr,
                            const int                   a_comp)
{
    // Gather geometric structures
    const Real dXiScale = a_dXi.product();
    const DisjointBoxLayout& grids = a_phi.getBoxes();
    DataIterator dit = grids.dataIterator();

    CH_assert(0 <= a_comp);
    CH_assert(a_comp < a_phi.nComp());
    CH_assert(dXiScale > 0.0);
    if (a_ccJPtr) {
        CH_assert(a_ccJPtr->getBoxes().compatible(grids));
        CH_assert(a_ccJPtr->nComp() == 1);
    }

    Real localSum = 0.0;

    // Loop over this level's grids, adding to sum.
    for (dit.reset(); dit.ok(); ++dit) {
        // Create references for convenience / allocate calucation space.
        const FArrayBox& phiFAB = a_phi[dit];
        const Box& valid = grids[dit];

        // If needed, we can adopt other centerings later.
        CH_assert(phiFAB.box().type() == IntVect::Zero);

        // Phi should be equal to or larger than the valid region.
        CH_assert(phiFAB.box().contains(valid));

        // Create a grid mask. 0 = covered by finer grid, 1 = valid data.
        BaseFab<int> maskFAB(valid, 1);
        Masks::zeroInvalid(maskFAB, a_finerGridsPtr, a_refRatio);

        // Add to localSum
        if (a_ccJPtr) {
            const FArrayBox& JFAB = (*a_ccJPtr)[dit];
            FORT_INTEGRAL_UNMAPPEDSUM1(
                CHF_REAL(localSum),
                CHF_REAL(a_vol),
                CHF_CONST_FRA1(phiFAB,a_comp),
                CHF_CONST_FRA1(JFAB,0),
                CHF_CONST_FIA1(maskFAB,0),
                CHF_BOX(valid),
                CHF_CONST_REAL(dXiScale));
        } else {
            FORT_INTEGRAL_UNMAPPEDSUM2(
                CHF_REAL(localSum),
                CHF_REAL(a_vol),
                CHF_CONST_FRA1(phiFAB,a_comp),
                CHF_CONST_FIA1(maskFAB,0),
                CHF_BOX(valid),
                CHF_CONST_REAL(dXiScale));
        }
    } // end loop over this level's grids
    return localSum;
}
