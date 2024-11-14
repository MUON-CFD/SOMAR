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
#include "Masks.H"
#include "LayoutData.H"
#include "AnisotropicRefinementTools.H"


// -----------------------------------------------------------------------------
// Create a grid mask. 0 = covered by finer grid, 1 = valid data.
// WARNING: a_maskFAB.box() must equal grids[dit]. No ghosts!!!
// -----------------------------------------------------------------------------
void
Masks::zeroInvalid (BaseFab<int>&            a_maskFAB,
                    const DisjointBoxLayout* a_finerGridsPtr,
                    const IntVect&           a_refRatio)
{
    a_maskFAB.setVal(1);

    if (a_finerGridsPtr != NULL) {
        LayoutIterator litFine = a_finerGridsPtr->layoutIterator();
        for (litFine.reset(); litFine.ok(); ++litFine) {
            // Calculate covered region
            Box coveredBox(a_finerGridsPtr->get(litFine()));
            coveredBox.coarsen(a_refRatio);
            coveredBox.convert(a_maskFAB.box().type());
            coveredBox &= a_maskFAB.box();

            // Set mask to zero on covered region
            if (!coveredBox.isEmpty()) {
                a_maskFAB.setVal(0, coveredBox, 0, 1);
            }
        }
    } // end if there is a finer level
}


// -----------------------------------------------------------------------------
// Sets invalid data to zero.
// -----------------------------------------------------------------------------
void
Masks::zeroInvalid(LevelData<FArrayBox>&    a_data,
                   const DisjointBoxLayout* a_finerGridsPtr,
                   const bool               a_ignoreGhosts,
                   const Real               a_invalidVal)
{
    IntVect refRatio;
    if (a_finerGridsPtr) {
        const Box& crseDomBox = a_data.getBoxes().physDomain().domainBox();
        const Box& fineDomBox = a_finerGridsPtr->physDomain().domainBox();
        refRatio = calculateRefinementRatio(crseDomBox, fineDomBox);
    }

    const DisjointBoxLayout& grids = a_data.getBoxes();
    DataIterator dit = a_data.dataIterator();

    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& dataFAB = a_data[dit];
        const Box& dataBox = dataFAB.box();
        const Box& valid = grids[dit];

        FArrayBox maskFAB(dataBox, 1);
        if (a_ignoreGhosts) {
            maskFAB.setVal(1.0);
        } else {
            maskFAB.setVal(a_invalidVal);
            maskFAB.setVal(1.0, valid, 0, 1);
        }

        if (a_finerGridsPtr) {
            LayoutIterator litFine = a_finerGridsPtr->layoutIterator();
            for (litFine.reset(); litFine.ok(); ++litFine) {
                Box coveredBox(a_finerGridsPtr->get(litFine()));
                coveredBox.coarsen(refRatio);
                coveredBox.convert(dataBox.type());
                coveredBox &= maskFAB.box();

                // Set mask to zero on covered region
                if (!coveredBox.isEmpty()) {
                    maskFAB.setVal(a_invalidVal, coveredBox, 0, 1);
                }
            }
        }

        for (int comp = 0; comp < dataFAB.nComp(); ++comp) {
            dataFAB.mult(maskFAB, 0, comp, 1);
        }
    }
}


// -----------------------------------------------------------------------------
// Sets invalid data to zero.
// -----------------------------------------------------------------------------
void
Masks::zeroInvalid(Vector<LevelData<FArrayBox>*> a_data,
                   const bool                    a_ignoreGhosts,
                   const Real                    a_invalidVal)
{
    const int numLevels = a_data.size();
    for (int l = 0; l < numLevels - 1; ++l) {
        zeroInvalid(*a_data[l],
                    &a_data[l + 1]->getBoxes(),
                    a_ignoreGhosts,
                    a_invalidVal);
    }
}


// -----------------------------------------------------------------------------
// Sets invalid data to zero. FC version.
// -----------------------------------------------------------------------------
void
Masks::zeroInvalid(LevelData<FluxBox>&      a_data,
                   const DisjointBoxLayout* a_finerGridsPtr,
                   const bool               a_ignoreGhosts,
                   const Real               a_invalidVal)
{
    IntVect refRatio;
    if (a_finerGridsPtr) {
        const Box& crseDomBox = a_data.getBoxes().physDomain().domainBox();
        const Box& fineDomBox = a_finerGridsPtr->physDomain().domainBox();
        refRatio = calculateRefinementRatio(crseDomBox, fineDomBox);
    }

    const DisjointBoxLayout& grids = a_data.getBoxes();
    DataIterator dit = a_data.dataIterator();

    for (dit.reset(); dit.ok(); ++dit) {
        for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
            FArrayBox& dataFAB = a_data[dit][fcDir];
            const Box& dataBox = dataFAB.box();
            const Box  valid   = surroundingNodes(grids[dit], fcDir);

            FArrayBox maskFAB(dataBox, 1);
            if (a_ignoreGhosts) {
                maskFAB.setVal(1.0);
            } else {
                maskFAB.setVal(a_invalidVal);
                maskFAB.setVal(1.0, valid, 0, 1);
            }

            if (a_finerGridsPtr) {
                LayoutIterator litFine = a_finerGridsPtr->layoutIterator();
                for (litFine.reset(); litFine.ok(); ++litFine) {
                    Box coveredBox(a_finerGridsPtr->get(litFine()));
                    coveredBox.coarsen(refRatio);
                    coveredBox.convert(dataBox.type());
                    coveredBox &= maskFAB.box();

                    // Set mask to zero on covered region
                    if (!coveredBox.isEmpty()) {
                        maskFAB.setVal(a_invalidVal, coveredBox, 0, 1);
                    }
                }
            }

            for (int comp = 0; comp < dataFAB.nComp(); ++comp) {
                dataFAB.mult(maskFAB, 0, comp, 1);
            }
        } // fcDir
    }
}
