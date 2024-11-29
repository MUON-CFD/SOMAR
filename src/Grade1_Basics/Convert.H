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
#ifndef ___Convert_H__INCLUDED___
#define ___Convert_H__INCLUDED___

#include "FArrayBox.H"
#include "FluxBox.H"
#include "LevelData.H"


class Convert
{
public:
    // The Swiss Army knife. Performs basic, 2nd-order averages.
    // Used to convert a CC FAB to FC FAB and so on.
    static void
    Simple (FArrayBox&       a_dest,
            const Interval&  a_destIvl,
            const Box&       a_destBox,
            const FArrayBox& a_src,
            const Interval&  a_srcIvl);

    static void
    Simple (FArrayBox&       a_dest,
            const int        a_destComp,
            const Box&       a_destBox,
            const FArrayBox& a_src,
            const int        a_srcComp);

    static void
    Simple (FArrayBox&       a_dest,
            const Box&       a_destBox,
            const FArrayBox& a_src);

    static void
    Simple (FArrayBox&       a_dest,
            const FArrayBox& a_src);

    // This version requires a_src.nComp() == SpaceDim * a_dest.nComp().
    // Comp c, FC in the d direction of a_dest is interpolated from comp
    // d + SpaceDim*c of a_dest.
    // This function is the complement of FacesToCells.
    static void
    CellsToFaces(FluxBox&         a_destFlub,
                 const FArrayBox& a_srcFAB,
                 const Box&       a_ccDestBox);

    // This version requires a_src.nComp() == SpaceDim * a_dest.nComp().
    // Comp c, FC in the d direction of a_dest is interpolated from comp
    // d + SpaceDim*c of a_dest.
    // This function is the complement of FacesToCells.
    static void
    CellsToFaces (LevelData<FluxBox>&         a_dest,
                  const LevelData<FArrayBox>& a_src,
                  const bool                  a_validOnly = false);

    // This version takes each comp of CC a_src and sends it to a comp of FC
    // a_dest. For example, if a_src has 2 comps, then a_dest must also have
    // 2 comps. Comp i of a_src will be sent to comp i of a_dest and will be
    // interpolated to each face of the original cell.
    static void
    CellsToAllFaces (FluxBox&         a_destFlub,
                     const FArrayBox& a_srcFAB);

    // This version takes each comp of CC a_src and sends it to a comp of FC
    // a_dest. For example, if a_src has 2 comps, then a_dest must also have
    // 2 comps. Comp i of a_src will be sent to comp i of a_dest and will be
    // interpolated to each face of the original cell.
    static void
    CellsToAllFaces (LevelData<FluxBox>&         a_dest,
                     const LevelData<FArrayBox>& a_src);

    // This version takes each comp and centering of FC a_src and sends it
    // to a comp of CC a_dest. Suppose a_src has 2 comps, then a_dest
    // must have 2*SpaceDim comps. If a_fcIsFastest == true, then comp c,
    // FC in the d direction of a_src, will be sent to comp d + SpaceDim*c of
    // a_dest. Otherwise, it wil be sent to comp c + a_src.nCpmp()*d.
    // This works on the overlap of the src and dest boxes.
    static void
    FacesToCells (LevelData<FArrayBox>&     a_dest,
                  const LevelData<FluxBox>& a_src,
                  const bool                a_fcIsFastest = true);
};


#endif //!___Convert_H__INCLUDED___
