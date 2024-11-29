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
 *  https://github.com/somarhub.
 ******************************************************************************/
#ifndef ___Integral_H__INCLUDED___
#define ___Integral_H__INCLUDED___

#include "LevelGeometry.H"
#include "BdryIter.H"


class Integral
{
public:
    // // Returns the integral of phi over the valid region.
    // // This version does not require a full LevelGeometry.
    // // This function can only handle cell-centered data.
    // static Real
    // sum (Real&                       a_vol,
    //      const LevelData<FArrayBox>& a_phi,
    //      const DisjointBoxLayout*    a_finerGridsPtr,
    //      const IntVect&              a_fineRefRatio,
    //      const RealVect&             a_dXi,
    //      const LevelData<FArrayBox>& a_CCJinv,
    //      const int                   a_comp);

    // Returns the integral of phi over the valid region.
    // This function can only handle cell-centered data.
    static Real
    sum(Real&                       a_vol,
        const LevelData<FArrayBox>& a_phi,
        const LevelGeometry&        a_levGeo,
        const bool                  a_sumJPhi = true,
        const int                   a_comp    = 0);

    // Returns the integral of phi over the valid region.
    // a_levGeo can be any levGeo in the hierarchy.
    // If a_sumJPhi == false, the we will sum a_phi directly, but still use J to
    // compute cell volumes.
    // This function can only handle cell-centered data.
    static Real
    sum(Real&                                a_vol,
        const Vector<LevelData<FArrayBox>*>& a_phi,
        const LevelGeometry&                 a_levGeo,
        const bool                           a_sumJPhi = true,
        const int                            a_comp    = 0,
        const int                            a_lBase   = 0);


    // These specialized versions do not require a_vol...

    // Returns the integral of phi over the valid region.
    // This function can only handle cell-centered data.
    static Real
    sum(const LevelData<FArrayBox>& a_phi,
        const LevelGeometry&        a_levGeo,
        const bool                  a_sumJPhi = true,
        const int                   a_comp    = 0)
    {
        Real vol;
        return sum(vol, a_phi, a_levGeo, a_sumJPhi, a_comp);
    }

    // Returns the integral of phi over the valid region.
    // a_levGeo can be any levGeo in the hierarchy.
    // If a_sumJPhi == false, the we will sum a_phi directly, but still use J to
    // compute cell volumes.
    // This function can only handle cell-centered data.
    static Real
    sum(const Vector<LevelData<FArrayBox>*>& a_phi,
        const LevelGeometry&                 a_levGeo,
        const bool                           a_sumJPhi = true,
        const int                            a_comp    = 0,
        const int                            a_lBase   = 0)
    {
        Real vol;
        return sum(vol, a_phi, a_levGeo, a_sumJPhi, a_comp, a_lBase);
    }

    // Returns the integral of Jphi over the valid region.
    // This version does not require a levGeo.
    // This function can only handle cell-centered data.
    static Real
    sum(const LevelData<FArrayBox>& a_Jphi,
        const RealVect&             a_dXi,
        const int                   a_comp = 0);

    // Same, but also returns the volume of the integrated region.
    // If a_Jptr == nullptr, then we just use ones.
    static void
    sum(Real&                       a_sum,
        Real&                       a_vol,
        const LevelData<FArrayBox>& a_phi,
        const LevelData<FArrayBox>* a_Jptr,
        const RealVect&             a_dXi,
        const int                   a_comp = 0);

    // Returns the integral of (v.n dS) over the level's boundary,
    // where n is the unit outward normal.
    // It is assumed v is in the curvilinear basis and scaled by J.
    static Real
    bdrySum(const LevelData<FluxBox>& a_v,
            const RealVect&           a_dXi,
            PhysBdryIter&             a_physBdryIter,
            CFIIter&                  a_cfiIter);

    // This version creates its own iterators.
    static Real
    bdrySum(const LevelData<FluxBox>& a_v,
            const LevelGeometry&      a_levGeo);

private:
    // These functions do the actual work on a single level...

    // Performs most of the computation for the sum functions. This version does
    // not perform MPI communication and only handles cell-centered data.
    // This will scale both phi and dXi by J.
    static Real
    localMappedSum (Real&                       a_vol,
                    const LevelData<FArrayBox>& a_phi,
                    const RealVect&             a_dXi,
                    const IntVect&              a_refRatio,
                    const DisjointBoxLayout*    a_finerGridsPtr,
                    const LevelData<FArrayBox>& a_ccJ,
                    const int                   a_comp);

    // Performs most of the computation for the sum functions. This version does
    // not perform MPI communication and only handles cell-centered data.
    // If a_ccJPtr is not NULL:
    //   a_phi will not be scaled by J.
    //   a_dXi _will_ be scaled by J.      <----- The only place J is used.
    // If a_ccJPtr is NULL:
    //   a_phi will not be scaled by J.
    //   a_dXi will not be scaled by J.
    static Real
    localUnmappedSum (Real&                       a_vol,
                      const LevelData<FArrayBox>& a_phi,
                      const RealVect&             a_dXi,
                      const IntVect&              a_refRatio,
                      const DisjointBoxLayout*    a_finerGridsPtr,
                      const LevelData<FArrayBox>* a_ccJPtr,
                      const int                   a_comp);
};


#endif //!___Integral_H__INCLUDED___
