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
#ifndef ___Masks_H__INCLUDED___
#define ___Masks_H__INCLUDED___

#include "FArrayBox.H"
#include "FluxBox.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"


class Masks
{
public:
    // TODO: this should this be called zeroCovered?

    // Create a grid mask. 0 = covered by finer grid, 1 = valid data.
    // WARNING: a_maskFAB.box() must equal grids[dit]. No ghosts!!!
    static void
    zeroInvalid (BaseFab<int>&            a_maskFAB,
                 const DisjointBoxLayout* a_finerGridsPtr,
                 const IntVect&           a_refRatio);

    // Sets invalid data to zero.
    static void
    zeroInvalid(LevelData<FArrayBox>&    a_data,
                const DisjointBoxLayout* a_finerGridsPtr,
                const bool               a_ignoreGhosts = false,
                const Real               a_invalidVal = 0.0);

    // Sets invalid data to zero.
    static void
    zeroInvalid(Vector<LevelData<FArrayBox>*> a_data,
                const bool                    a_ignoreGhosts = false,
                const Real                    a_invalidVal = 0.0);

    // Sets invalid data to zero. FC version.
    static void
    zeroInvalid(LevelData<FluxBox>&      a_data,
                const DisjointBoxLayout* a_finerGridsPtr,
                const bool               a_ignoreGhosts = false,
                const Real               a_invalidVal = 0.0);
};


#endif //!___Masks_H__INCLUDED___
