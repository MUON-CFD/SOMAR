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
#include "BilinearInterp.H"
#include "BilinearInterpF_F.H"
#include "SOMAR_Constants.H"
#include "BCTools.H"
#include "Debug.H"


// ======================== Private helper functions ===========================
namespace {

// -----------------------------------------------------------------------------
std::vector<Real>
addGhostLayer(const std::vector<Real>& a_src) {
    std::vector<Real> growSrc(a_src.size() + 2);

    CH_verify(a_src.size() >= 2);
    growSrc[0] = 2.0 * a_src[0] - a_src[1];

    size_t gi = 1;
    for (const auto sx : a_src) {
        growSrc[gi++] = sx;
    }

    const int si = a_src.size() - 1;
    growSrc[gi] = 2.0 * a_src[si] - a_src[si - 1];

    return growSrc;
}


// -----------------------------------------------------------------------------
void
defineAndAddGhostLayer(FArrayBox&       a_growFAB,
                       const FArrayBox& a_srcFAB,
                       const Box&       a_validSrcBox)
{
    CH_assert(&a_growFAB != &a_srcFAB);
    a_growFAB.define(grow(a_validSrcBox, IntVect::Unit), a_srcFAB.nComp());

    debugInit(a_growFAB);
    a_growFAB.copy(a_srcFAB, a_validSrcBox);

    Box bx = a_validSrcBox;
    for (int dir = 0; dir < SpaceDim; ++dir) {
        constexpr int extrapOrder = 2;
        BCTools::extrap(a_growFAB, bx, dir, Side::Lo, extrapOrder);
        BCTools::extrap(a_growFAB, bx, dir, Side::Hi, extrapOrder);
        bx.grow(dir, 1);
    }
}


// -----------------------------------------------------------------------------
// Careful, there are no checks! All data holders must be correct in size.
// a_srcIV must point to the *lower* end of the line of IntVects.
void
copyLine(std::vector<Real>& a_dest,
         const FArrayBox&   a_srcFAB,
         const int          a_srcComp,
         IntVect            a_srcIV,
         const int          a_lineDir)
{
    size_t idx = 0;
    while (idx < a_dest.size()) {
        a_dest[idx] = a_srcFAB(a_srcIV, a_srcComp);
        ++idx;
        ++a_srcIV[a_lineDir];
    }
}


// -----------------------------------------------------------------------------
// Careful, there are no checks! All data holders must be correct in size.
// a_destIV must point to the *lower* end of the line of IntVects.
void
copyLine(FArrayBox&               a_destFAB,
         const int                a_destComp,
         IntVect                  a_destIV,
         const std::vector<Real>& a_src,
         const int                a_lineDir)
{
    size_t idx = 0;
    while (idx < a_src.size()) {
        a_destFAB(a_destIV, a_destComp) = a_src[idx];
        ++idx;
        ++a_destIV[a_lineDir];
    }
}


// -----------------------------------------------------------------------------
// Searches a vector for a value. Assumes vector is sorted. Returns index at or
// just below value.
size_t
binarySearch(const Real               a_val,
             const std::vector<Real>& a_vec,
             const size_t             a_startIdx = 0)
{
    size_t ilo = a_startIdx;
    size_t ihi = a_vec.size() - 1;

    while ((ihi - ilo) > 1) {
        const size_t i = (ihi + ilo) / 2;
        if (a_val < a_vec[i]) {
            ihi = i;
        } else {
            ilo = i;
        }
    }

    CH_assert(ilo < a_vec.size());
    return ilo;
};


// -----------------------------------------------------------------------------
void
LineLinearInterp(std::vector<Real>&       a_destf,
                 const std::vector<Real>& a_destx,
                 const std::vector<Real>& a_srcf,
                 const std::vector<Real>& a_srcx)
{
    CH_assert(a_destf.size() == a_destx.size());
    CH_assert(a_srcf.size() == a_srcx.size());

    size_t di = 0;
    for (size_t si = 0; si < a_srcx.size() - 1; ++si) {
        const Real xl    = a_srcx[si];
        const Real fl    = a_srcf[si];
        const Real xr    = a_srcx[si + 1];
        const Real fr    = a_srcf[si + 1];
        const Real invDx = 1.0 / (xr - xl);

        Real x = a_destx[di];

        if (RealCmp::eq(x, xr)) {
            a_destf[di] = fr;
            ++di;
            if (di == a_destx.size()) break;
            continue;
        }

        while (x < xr) {
            const Real c = (x - xl) * invDx;
            a_destf[di] = (1.0 - c) * fl + c * fr;
            ++di;
            if (di == a_destx.size()) break;
            x = a_destx[di];
        }

        if (RealCmp::eq(x, xr)) {
            a_destf[di] = fr;
            ++di;
            if (di == a_destx.size()) break;
        }
    } // si

    CH_verify(di == a_destx.size());
}


// -----------------------------------------------------------------------------
void
LineLinearInterp(FArrayBox&               a_destfFAB,
                 const std::vector<Real>& a_destx,
                 const Box&               a_destBox,
                 const FArrayBox&         a_srcfFAB,
                 const std::vector<Real>& a_srcx,
                 const Box&               a_srcBox,
                 const int                a_interpDir)
{
#ifndef NDEBUG
    CH_assert(0 <= a_interpDir && a_interpDir < SpaceDim);
    CH_assert(a_destfFAB.box().contains(a_destBox));
    CH_assert( a_srcfFAB.box().contains( a_srcBox));
    CH_assert(a_destx.size() == static_cast<size_t>(a_destBox.size(a_interpDir)));
    CH_assert( a_srcx.size() == static_cast<size_t>( a_srcBox.size(a_interpDir)));
    CH_assert(a_destfFAB.nComp() == a_srcfFAB.nComp());

    for (int d = 0; d < SpaceDim; ++d) {
        if (d != a_interpDir) {
            CH_assert(a_srcBox.smallEnd(d) == a_destBox.smallEnd(d));
            CH_assert(a_srcBox.  bigEnd(d) == a_destBox.  bigEnd(d));
        }
    }
#endif

    // Find the first and last needed src points.
    const size_t losi = binarySearch(a_destx.front(), a_srcx);
    const size_t hisi = [&]() {
        size_t idx = binarySearch(a_destx.back(), a_srcx, losi);
        if (idx < a_srcx.size() - 1) ++idx;
        return idx;
    } ();

    CH_assert(a_srcx[losi] <= a_destx.front());
    CH_assert(a_destx.back() <= a_srcx[hisi]);

    // Remove the src data that we do not need.
    const size_t croppedSrcSize = hisi - losi + 1;
    std::vector<Real> croppedSrcx(croppedSrcSize);
    for (size_t csi = 0, si = losi; si <= hisi; ++si) {
        croppedSrcx[csi++] = a_srcx[si];
    }

    // Flatten the src box to exist at the starting idx of each interp line.
    const size_t flatSrcIdx = a_srcfFAB.box().smallEnd(a_interpDir) + losi;
    const Box flatSrcBox = Box(a_srcBox).setRange(a_interpDir, flatSrcIdx, 1);
    const IntVect offset =
        BASISV(a_interpDir) * (a_destBox.smallEnd(a_interpDir) - flatSrcIdx);

    // Do it!
    // a_srcfFAB [copy to] croppedSrcf [interp to] destf [copy to] a_destfFAB.
    std::vector<Real> destf(a_destx.size());
    std::vector<Real> croppedSrcf(croppedSrcSize);
    for (BoxIterator sbit(flatSrcBox); sbit.ok(); ++sbit) {
        const IntVect& siv = sbit();
        const IntVect  div = siv + offset;

        for (int comp = 0; comp < a_srcfFAB.nComp(); ++comp) {
            copyLine(croppedSrcf, a_srcfFAB, comp, siv, a_interpDir);
            LineLinearInterp(destf, a_destx, croppedSrcf, croppedSrcx);
            copyLine(a_destfFAB, comp, div, destf, a_interpDir);
        }
    }
}

};  // namespace



// ============================= Public functions =============================

// -----------------------------------------------------------------------------
void
LinearInterp(FArrayBox&                                     a_destfFAB,
             const std::array<std::vector<Real>, SpaceDim>& a_destx,
             const Box&                                     a_destBox,
             const FArrayBox&                               a_srcfFAB,
             const std::array<std::vector<Real>, SpaceDim>& a_srcx,
             const Box&                                     a_srcBox,
             const int                                      a_srcExtrapRadius)
{
    // Sanity checks.
    // These will even be done in release mode because they are fast checks
    // and users will often need to interpolate ICs. It's better to catch
    // an error and print a helpful message before the simulation starts.

    // Check FAB comps.
    if (a_destfFAB.nComp() != a_srcfFAB.nComp()) {
        MAYDAYERROR("These two must have the number of components:"
                    << "\n\ta_srcfFAB.nComp() = " << a_srcfFAB.nComp()
                    << "\n\ta_destfFAB.nComp() = " << a_destfFAB.nComp()
                    << ".");
    }

    // Check coordinate vector sizes.
    for (int dir = 0; dir < SpaceDim; ++dir) {
        if (static_cast<size_t>(a_srcBox.size(dir)) != a_srcx[dir].size()) {
            MAYDAYERROR("These two must have the same length:"
                        << "\n\ta_srcBox.size(" << dir << ") = " << a_srcBox.size()
                        << "\n\ta_srcx[" << dir << "].size() = " << a_srcx[dir].size()
                        << ".");
        }
    }

    // Check box types.
    if (!a_srcfFAB.box().sameType(a_srcBox)) {
        MAYDAYERROR("a_srcfFAB.box() must have the same centering as a_srcBox:"
                    << "\n\ta_srcfFAB.box() = " << a_srcfFAB.box()
                    << "\n\ta_srcBox        = " << a_srcBox);
    }
    if (!a_destfFAB.box().sameType(a_destBox)) {
        MAYDAYERROR("a_destfFAB.box() must have the same centering as a_destBox:"
                    << "\n\ta_destfFAB.box() = " << a_destfFAB.box()
                    << "\n\ta_destBox        = " << a_destBox);
    }

    // Check box extents.
    if (!a_srcfFAB.box().contains(a_srcBox)) {
        MAYDAYERROR("a_srcfFAB.box() must contain a_srcBox:"
                    << "\n\ta_srcfFAB.box() = " << a_srcfFAB.box()
                    << "\n\ta_srcBox        = " << a_srcBox);
    }
    if (!a_destfFAB.box().contains(a_destBox)) {
        MAYDAYERROR("a_destfFAB.box() must contain a_destBox:"
                    << "\n\ta_destfFAB.box() = " << a_destfFAB.box()
                    << "\n\ta_destBox        = " << a_destBox);
    }

    // Initialize *srcptr to a_srcfFAB with one ghost layer added in each
    // dir we will interpolate.
    FArrayBox _f0;
    FArrayBox* srcptr = &_f0;
    if (a_srcExtrapRadius == 0) {
        srcptr->define(a_srcBox, a_srcfFAB.nComp());
        srcptr->copy(a_srcfFAB);
    } else if (a_srcExtrapRadius > 0) {
        defineAndAddGhostLayer(*srcptr, a_srcfFAB, a_srcBox);
    }

    for (int r = 1; r < a_srcExtrapRadius; ++r) {
        FArrayBox tmpFAB;
        defineAndAddGhostLayer(tmpFAB, *srcptr, srcptr->box());
        srcptr->define(tmpFAB.box(), tmpFAB.nComp());
        srcptr->copy(tmpFAB);
    }

    // Initialize destBox to srcptr's box. As we loop over the dirs,
    // this will get squashed to match the user's a_destBox in that dir.
    Box destBox = srcptr->box();
    FArrayBox _f1;
    FArrayBox* destptr = &_f1;

    for (int dir = 0; dir < SpaceDim; ++dir) {
        // Squash destBox.
        //  Set destBox to cover user's a_destBox in this dir.
        //  Dirs < this dir should already be covered.
        //  Dirs > this dir should match the grown src.
        if (a_destBox.type(dir) == IndexType::CELL) {
            destBox.setRange(dir, a_destBox.smallEnd(dir), a_destBox.size(dir));
        } else {
            destBox.setRange(dir, a_destBox.smallEnd(dir), a_destBox.size(dir) - 1);
            destBox.surroundingNodes(dir);
        }
        destptr->resize(destBox, a_srcfFAB.nComp());
        debugInit(*destptr);

        // Get src x values with ghosts included.
        // const std::vector<Real> growSrcx = addGhostLayer(a_srcx[dir]);
        // CH_assert(growSrcx.size() == a_srcx[dir].size() + 2);

        const std::vector<Real> growSrcx = [&]() {
            std::vector<Real> s = a_srcx[dir];
            for (int r = 0; r < a_srcExtrapRadius; ++r) {
                s = addGhostLayer(s);
            }
            return s;
        }();
        CH_assert(growSrcx.size() == a_srcx[dir].size() + 2 * a_srcExtrapRadius);

        CH_assert(growSrcx.front() <= a_destx[dir].front());
        CH_assert(a_destx[dir].back() <= growSrcx.back());

        // Interpolate from grown src to destBox/*destptr in this dir .
        LineLinearInterp(*destptr,
                         a_destx[dir],
                         destBox,
                         *srcptr,
                         growSrcx,
                         srcptr->box(),
                         dir);

        // The interpolated data becomes the src for the next interp dir.
        std::swap(srcptr, destptr);
    }

    // Copy the interpolated data to the user's holder.
    a_destfFAB.copy(*srcptr, a_destBox);
}


// -----------------------------------------------------------------------------
// The nodes will be labeled as A, B, C, and D in this order:   C---D
//                                                              |   |
//                                                              A---B
// Simple bilinear interp where
// f_{u,v}=fA*(1-u)*(1-v)+FB*u*(1-v)+FC*(1-v)*u+FD*u*v
// where (u,v) are the coordinate in the middle of the box
// -----------------------------------------------------------------------------
void
BilinearInterp2D(FArrayBox&               a_destf,
                 const std::vector<Real>& a_destx,
                 const std::vector<Real>& a_desty,
                 const Box&               a_destBox,
                 const FArrayBox&         a_srcf,
                 std::vector<Real>        a_srcx,
                 std::vector<Real>        a_srcy,
                 Box                      a_srcBox)
{
#ifndef NDEBUG
    {
        // Check number of comps.
        CH_assert(a_srcf .nComp() == a_destf.nComp());

        // Check centerings.
        CH_assert(a_destf.box().sameType(a_destBox));
        CH_assert(a_srcf .box().sameType(a_srcBox));

        // Check FAB regions.
        CH_assert(a_destf.box().contains(a_destBox));
        CH_assert(a_destx.size() == static_cast<size_t>(a_destBox.size(0)));
        CH_assert(a_desty.size() == static_cast<size_t>(a_destBox.size(1)));
        CH_assert(a_srcf .box().contains(a_srcBox));
        CH_assert(a_srcx.size() == static_cast<size_t>(a_srcBox.size(0)));
        CH_assert(a_srcy.size() == static_cast<size_t>(a_srcBox.size(1)));
    }
#endif

    FArrayBox srcfGrow;
    {
        // Is source region large enough?
        const bool needsLoGhostx = (a_destx.front() < a_srcx.front());
        const bool needsHiGhostx = (a_srcx.back() < a_destx.back());
        const bool needsLoGhosty = (a_desty.front() < a_srcy.front());
        const bool needsHiGhosty = (a_srcy.back() < a_desty.back());

        // Grow the src Box accordingly.
        if (needsLoGhostx) a_srcBox.growLo(0, 1);
        if (needsHiGhostx) a_srcBox.growHi(0, 1);
        if (needsLoGhosty) a_srcBox.growLo(1, 1);
        if (needsHiGhosty) a_srcBox.growHi(1, 1);

        // Extrapolate srcf.
        srcfGrow.define(a_srcBox, a_srcf.nComp());
        srcfGrow.copy(a_srcf, a_srcf.box());
        if (needsLoGhostx)
            BCTools::extrap(srcfGrow, a_srcf.box(), 0, Side::Lo, 2);
        if (needsHiGhostx)
            BCTools::extrap(srcfGrow, a_srcf.box(), 0, Side::Hi, 2);
        if (needsLoGhosty)
            BCTools::extrap(srcfGrow, a_srcf.box(), 1, Side::Lo, 2);
        if (needsHiGhosty)
            BCTools::extrap(srcfGrow, a_srcf.box(), 1, Side::Hi, 2);

        // Extrapolate coordinates.
        if (needsLoGhostx) {
            const int i = 0;
            a_srcx.insert(a_srcx.begin(),
                          BCTools::quadraticExtrap(
                              a_srcx[i], a_srcx[i + 1], a_srcx[i + 2]));
        }
        if (needsHiGhostx) {
            const int i = a_srcx.size() - 1;
            a_srcx.push_back(BCTools::quadraticExtrap(
                a_srcx[i], a_srcx[i - 1], a_srcx[i - 2]));
        }
        if (needsLoGhosty) {
            const int j = 0;
            a_srcy.insert(a_srcy.begin(),
                          BCTools::quadraticExtrap(
                              a_srcy[j], a_srcy[j + 1], a_srcy[j + 2]));
        }
        if (needsHiGhosty) {
            const int j = a_srcy.size() - 1;
            a_srcy.push_back(BCTools::quadraticExtrap(
                a_srcy[j], a_srcy[j - 1], a_srcy[j - 2]));
        }

        // Does one ghost cut it? If not, throw an error.
        bool stillTooSmall = (a_destx.front() < a_srcx.front());
        stillTooSmall |= (a_srcx.back() < a_destx.back());
        stillTooSmall |= (a_desty.front() < a_srcy.front());
        stillTooSmall |= (a_srcy.back() < a_desty.back());
        if (stillTooSmall) {
            MAYDAYERROR("a_srcBox is not large enough to accomodate stencil.");
        }
    };

    // Interpolate!
    const IntVect destBoxShift = a_destBox.smallEnd();
    const IntVect srcBoxShift = a_srcBox.smallEnd();
    FORT_BILINEARINTERP2DF (
        CHF_FRA_SHIFT(a_destf, destBoxShift),
        CHF_CONST_VR(a_destx),
        CHF_CONST_VR(a_desty),
        CHF_BOX_SHIFT(a_destBox, destBoxShift),
        CHF_CONST_FRA_SHIFT(srcfGrow, srcBoxShift),
        CHF_CONST_VR(a_srcx),
        CHF_CONST_VR(a_srcy));
}


// -----------------------------------------------------------------------------
// The nodes will be labeled as A, B, C, and D in this order:   C---D
//                                                              |   |
//                                                              A---B
// Simple bilinear interp where
// f_{u,v}=fA*(1-u)*(1-v)+FB*u*(1-v)+FC*(1-v)*u+FD*u*v
// where (u,v) are the coordinate in the middle of the box
// -----------------------------------------------------------------------------
void
BilinearInterp2D(FArrayBox&               a_fInterp,
                 const FArrayBox&         a_xInterp,
                 const FArrayBox&         a_yInterp,
                 const Box&               a_interpBox,
                 const int                a_xdir,
                 const int                a_ydir,
                 const FArrayBox&         a_f,
                 const std::vector<Real>& a_x,
                 const std::vector<Real>& a_y)
{
#ifndef NDEBUG
    {
        // Check centerings.
        CH_assert(a_fInterp.box().type() == a_interpBox.type());
        CH_assert(a_xInterp.box().type() == a_interpBox.type());
        CH_assert(a_yInterp.box().type() == a_interpBox.type());

        // Check FAB regions.
        CH_assert(a_fInterp.box().contains(a_interpBox));
        CH_assert(a_xInterp.box().contains(a_interpBox));
        CH_assert(a_yInterp.box().contains(a_interpBox));
        CH_assert(a_f.box().size(a_xdir) == static_cast<int>(a_x.size()));
        CH_assert(a_f.box().size(a_ydir) == static_cast<int>(a_y.size()));

        // Check number of comps.
        CH_assert(a_xInterp.nComp() == 1);
        CH_assert(a_yInterp.nComp() == 1);
        CH_assert(a_f      .nComp() == a_fInterp.nComp());

        // Check dirs
        CH_assert(0 <= a_xdir);
        CH_assert(a_xdir < SpaceDim);
        CH_assert(0 <= a_ydir);
        CH_assert(a_ydir < SpaceDim);
        CH_assert(a_xdir != a_ydir);
    }
#endif

    // The fortran function can only handle right-handed
    // permutations of the x, y, and z directions.
    if ((a_xdir == 0 && a_ydir == 1) ||
        (a_xdir == 1 && a_ydir == 2) ||
        (a_xdir == 2 && a_ydir == 0)) {

        // The fortran function requires the vector indices to
        // coincide with the constraint indices.
        IntVect shift = IntVect::Zero;
        shift[a_xdir] = a_f.box().smallEnd(a_xdir);
        shift[a_ydir] = a_f.box().smallEnd(a_ydir);

        // Interpolate
        FORT_BILINEARINTERP2DF_GENERAL (
            CHF_FRA(a_fInterp),
            CHF_CONST_FRA1(a_xInterp,0),
            CHF_CONST_FRA1(a_yInterp,0),
            CHF_BOX(a_interpBox),
            CHF_CONST_INT(a_xdir),
            CHF_CONST_INT(a_ydir),
            CHF_CONST_VR(a_x),
            CHF_CONST_VR(a_y),
            CHF_CONST_FRA_SHIFT(a_f,shift));
    } else {

        // The fortran function requires the vector indices to
        // coincide with the constraint indices.
        IntVect shift = IntVect::Zero;
        shift[a_ydir] = a_f.box().smallEnd(a_xdir);
        shift[a_xdir] = a_f.box().smallEnd(a_ydir);

        // Interpolate
        FORT_BILINEARINTERP2DF_GENERAL (
            CHF_FRA(a_fInterp),
            CHF_CONST_FRA1(a_yInterp,0),
            CHF_CONST_FRA1(a_xInterp,0),
            CHF_BOX(a_interpBox),
            CHF_CONST_INT(a_ydir),
            CHF_CONST_INT(a_xdir),
            CHF_CONST_VR(a_y),
            CHF_CONST_VR(a_x),
            CHF_CONST_FRA_SHIFT(a_f,shift));
    }
}

