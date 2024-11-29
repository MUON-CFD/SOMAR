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
 *  https://github.com/MUON-CFD/somar.
 ******************************************************************************/
#include "Subspace.H"
#include "SubspaceF_F.H"


// -----------------------------------------------------------------------------
// Takes an N-dimensional Box and projects it to a subspace defined by a_mask.
// a_mask should be 0 in flattened directions and 1 in unmodified directions.
// In other words, a_mask flags the tangential directions of the subspace.
// -----------------------------------------------------------------------------
Box
Subspace::flattenBox(const Box& a_box, const IntVect& a_mask)
{
    D_TERM(CH_assert(a_mask[0] == 0 || a_mask[0] == 1);
           , CH_assert(a_mask[1] == 0 || a_mask[1] == 1);
           , CH_assert(a_mask[2] == 0 || a_mask[2] == 1);)

    const IntVect boxType = a_box.type();
    const IntVect smallIV = a_mask * a_box.smallEnd();
    const IntVect bigIV   = a_mask * a_box.bigEnd();

    Box retBox(smallIV, bigIV, boxType);
    return retBox;
}


// -----------------------------------------------------------------------------
// Takes an N-dimensional Box and projects it to the (N_1)-dimensional surface
// whose normal direction is a_normDir.
// -----------------------------------------------------------------------------
Box
Subspace::flattenBox(const Box& a_box, const int a_normDir)
{
    CH_assert(0 <= a_normDir);
    CH_assert(a_normDir < SpaceDim);

    IntVect mask = IntVect::Unit - BASISV(a_normDir);
    return Subspace::flattenBox(a_box, mask);
}


// -----------------------------------------------------------------------------
// Returns a box that can define holders for functions that only depend on the
// vertical coordinate.
// -----------------------------------------------------------------------------
Box
Subspace::verticalDataBox(const ProblemDomain& a_domain)
{
    return Subspace::flattenBox(a_domain.domainBox(), BASISV(SpaceDim - 1));
}


// -----------------------------------------------------------------------------
// Returns a box that can define holders for functions that only depend on the
// vertical coordinate.
// -----------------------------------------------------------------------------
Box
Subspace::verticalDataBox(const Box& a_box)
{
    return Subspace::flattenBox(a_box, BASISV(SpaceDim - 1));
}


// -----------------------------------------------------------------------------
// Returns a box that can define holders for functions that only depend on the
// horizontal coordinates.
// -----------------------------------------------------------------------------
Box
Subspace::horizontalDataBox(const ProblemDomain& a_domain)
{
    return Subspace::flattenBox(a_domain.domainBox(), SpaceDim - 1);
}


// -----------------------------------------------------------------------------
// Returns a box that can define holders for functions that only depend on the
// horizontal coordinates.
// -----------------------------------------------------------------------------
Box
Subspace::horizontalDataBox(const Box& a_box)
{
    return Subspace::flattenBox(a_box, SpaceDim - 1);
}


// -----------------------------------------------------------------------------
// Checks if a_b1 and a_b2 are the same type in all directions except a_dir.
// This function does not check the types of the boxes in a_dir at all.
// -----------------------------------------------------------------------------
bool
Subspace::sameTypeTransverse(const Box& a_b1, const Box& a_b2, const int a_dir)
{
#if CH_SPACEDIM == 2
    const int d1 = (a_dir + 1) % SpaceDim;
    return a_b1.type(d1) == a_b2.type(d1);
#else
    const int d1 = (a_dir + 1) % SpaceDim;
    const int d2 = (a_dir + 2) % SpaceDim;
    return (a_b1.type(d1) == a_b2.type(d1))
        && (a_b1.type(d2) == a_b2.type(d2));
#endif
}


// -----------------------------------------------------------------------------
// Static utility.
// Takes a vertical line of data and injects it into a holder on a DBL.
// -----------------------------------------------------------------------------
void
Subspace::horizontalExtrusion(LevelData<FArrayBox>& a_dest,
                              const int             a_destComp,
                              const FArrayBox&      a_srcFAB,
                              const int             a_srcComp,
                              const int             a_numComps)
{
    DataIterator dit = a_dest.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        Subspace::horizontalExtrusion(
            a_dest[dit], a_destComp, a_srcFAB, a_srcComp, a_numComps);
    }
}


// -----------------------------------------------------------------------------
// Static utility.
// Takes a vertical line of data and injects it into a holder on a DBL.
// -----------------------------------------------------------------------------
void
Subspace::horizontalExtrusion(FArrayBox&       a_destFAB,
                              const int        a_destComp,
                              const FArrayBox& a_srcFAB,
                              const int        a_srcComp,
                              const int        a_numComps)
{
    const Interval sIvl(a_srcComp, a_srcComp + a_numComps - 1);
    const Interval dIvl(a_destComp, a_destComp + a_numComps - 1);

    FArrayBox       dFAB(dIvl, a_destFAB);
    const FArrayBox sFAB(sIvl, (FArrayBox&)a_srcFAB);
    const Box&      dBox = dFAB.box();

    CH_assert(sFAB.box().type(SpaceDim - 1) == dBox.type(SpaceDim - 1));

    FORT_HORIZONTALEXTRUSION(CHF_FRA(dFAB),
                             CHF_CONST_FRA(sFAB),
                             CHF_BOX(dBox),
                             CHF_CONST_INT(a_numComps));
}


// -----------------------------------------------------------------------------
// Sets q = q + scale * qbar
// -----------------------------------------------------------------------------
void
Subspace::addHorizontalExtrusion(LevelData<FArrayBox>& a_q,
                                 const int             a_qComp,
                                 const FArrayBox&      a_qbarFAB,
                                 const int             a_qbarComp,
                                 const int             a_numComps,
                                 const Real            a_scale)
{
    DataIterator dit = a_q.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        Subspace::addHorizontalExtrusion(
            a_q[dit], a_qComp, a_qbarFAB, a_qbarComp, a_numComps, a_scale);
    }
}


// -----------------------------------------------------------------------------
// Static utility.
// Sets q = q + scale * qbar. FAB version.
// -----------------------------------------------------------------------------
void
Subspace::addHorizontalExtrusion(FArrayBox&       a_qFAB,
                                 const int        a_qComp,
                                 const FArrayBox& a_qbarFAB,
                                 const int        a_qbarComp,
                                 const int        a_numComps,
                                 const Real       a_scale)
{
    const Interval qIvl(a_qComp, a_qComp + a_numComps - 1);
    const Interval qbarIvl(a_qbarComp, a_qbarComp + a_numComps - 1);

    FArrayBox       qFAB(qIvl, a_qFAB);
    const FArrayBox qbarFAB(qbarIvl, (FArrayBox&)a_qbarFAB);
    const Box&      qBox = qFAB.box();

    CH_assert(qBox.type(SpaceDim - 1) == qbarFAB.box().type(SpaceDim - 1));

    FORT_ADDHORIZONTALEXTRUSION(CHF_FRA(a_qFAB),
                                CHF_CONST_FRA(a_qbarFAB),
                                CHF_BOX(qBox),
                                CHF_CONST_INT(a_numComps),
                                CHF_CONST_REAL(a_scale));
}
