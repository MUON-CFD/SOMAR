#ifdef CH_LANG_CC
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

// Provides utilities for defining and storing data on flattened grids.
#endif


#ifndef __Subspace_H__INCLUDED__
#define __Subspace_H__INCLUDED__
#ifdef CH_LANG_CC

#include "ProblemDomain.H"
#include "LevelData.H"


class Subspace
{
public:
    // Takes an N-dimensional Box and projects it to a subspace defined by
    // a_mask. a_mask should be 0 in flattened directions and 1 in unmodified
    // directions. In other words, a_mask flags the tangential directions of the
    // subspace.
    static Box
    flattenBox(const Box& a_box, const IntVect& a_mask);

    // Takes an N-dimensional Box and projects it to the (N_1)-dimensional
    // surface whose normal direction is a_normDir.
    static Box
    flattenBox(const Box& a_box, const int a_normDir);

    // Returns a box that can define holders for functions that only depend on
    // the vertical coordinate.
    static Box
    verticalDataBox(const ProblemDomain& a_domain);

    // Returns a box that can define holders for functions that only depend on
    // the vertical coordinate.
    static Box
    verticalDataBox(const Box& a_box);

    // Returns a box that can define holders for functions that only depend on
    // the horizontal coordinates.
    static Box
    horizontalDataBox(const ProblemDomain& a_domain);

    // Returns a box that can define holders for functions that only depend on
    // the horizontal coordinates.
    static Box
    horizontalDataBox(const Box& a_box);

    // Checks if a_b1 and a_b2 are the same type in all directions except a_dir.
    // This function does not check the types of the boxes in a_dir at all.
    static bool
    sameTypeTransverse(const Box& a_b1, const Box& a_b2, const int a_dir);

    // Takes a vertical line of data and injects it into a holder on a DBL.
    static void
    horizontalExtrusion(LevelData<FArrayBox>& a_dest,
                        const int             a_destComp,
                        const FArrayBox&      a_srcFAB,
                        const int             a_srcComp,
                        const int             a_numComps = 1);

    // FAB version.
    static void
    horizontalExtrusion(FArrayBox&       a_destFAB,
                        const int        a_destComp,
                        const FArrayBox& a_srcFAB,
                        const int        a_srcComp,
                        const int        a_numComps = 1);

    // Sets q = q + scale * qbar, where qbar is a vertical line of data.
    static void
    addHorizontalExtrusion(LevelData<FArrayBox>& a_q,
                           const int             a_qComp,
                           const FArrayBox&      a_qbarFAB,
                           const int             a_qbarComp,
                           const int             a_numComps,
                           const Real            a_scale);

    // FAB version.
    static void
    addHorizontalExtrusion(FArrayBox&       a_qFAB,
                           const int        a_qComp,
                           const FArrayBox& a_qbarFAB,
                           const int        a_qbarComp,
                           const int        a_numComps,
                           const Real       a_scale);
};

#else // End C++ code, begin Fortran code.

#if CH_SPACEDIM == 2
#  define FLATIX(i)       i##0     , 0
#  define FLATOFFSET(i,o) i##0 o##0, 0
#  define VERTCOMP(i)                i##1
#elif CH_SPACEDIM == 3
#  define FLATIX(i)       i##0     , i##1     , 0
#  define FLATOFFSET(i,o) i##0 o##0, i##1 o##1, 0
#  define VERTCOMP(i)                           i##2
#else
#  error Bad SpaceDim
#endif

#endif // End of fortran code
#endif //!__Subspace_H__INCLUDED__
