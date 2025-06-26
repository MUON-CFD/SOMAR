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
#include "LayoutTools.H"
#include "MergeBoxesOnLines.H"
#include "AnisotropicMeshRefine.H"
#include "NeighborIterator.H"
#include "AnisotropicRefinementTools.H"
#include "Convert.H"
#include "MiscUtils.H"
#include "Debug.H"

namespace LayoutTools {

// -----------------------------------------------------------------------------
// Takes a DBL and returns a new one whose boxes are merged in a_dir.
// -----------------------------------------------------------------------------
void mergeLayout (DisjointBoxLayout&       a_newLayout,
                  const DisjointBoxLayout& a_origLayout,
                  const int                a_dir)
{
    // Sanity checks
    CH_assert(0 <= a_dir);
    CH_assert(a_dir < SpaceDim);

    // Merge the boxes
    Vector<Box> vbox = a_origLayout.boxArray();
    MergeBoxesOnLines().mergeBoxes(vbox, a_dir);

    // Create the merged layout
    a_newLayout.defineAndLoadBalance(vbox, NULL, a_origLayout.physDomain());
}


// -----------------------------------------------------------------------------
// Splits a box into a load balanced set of boxes that are not split in a_dir.
// -----------------------------------------------------------------------------
void lineLayout (DisjointBoxLayout&   a_newLayout,
                 const ProblemDomain& a_domain,
                 const int            a_dir)
{
    // Sanity checks
    CH_assert(0 <= a_dir);
    CH_assert(a_dir < SpaceDim);

    const int np = numProc();
    IntVect maxBoxSize = a_domain.size() / np;
    maxBoxSize[a_dir] = 0;

    Vector<Box> vbox;
    AnisotropicMeshRefine::domainSplit(a_domain.domainBox(), vbox, maxBoxSize, 1);

    a_newLayout.defineAndLoadBalance(vbox, NULL, a_domain);
}


// -----------------------------------------------------------------------------
// Define a dbl with just one box.
// This function must be called on ALL procs, but a_box only needs to be
// defined on a_srcProc.
// -----------------------------------------------------------------------------
void defineOneProcGrids (DisjointBoxLayout&   a_grids,
                         const ProblemDomain& a_domain,
                         Box                  a_box,
                         const int            a_srcProc)
{
    broadcast(a_box, a_srcProc);
    Vector<Box> boxArray(1, a_box);
    Vector<int> procArray(1, a_srcProc);
    a_grids.define(boxArray, procArray, a_domain);
}


// -----------------------------------------------------------------------------
void
defineImpartialCopier(Copier&              a_copier,
                      const BoxLayout&     a_srcLayout,
                      const BoxLayout&     a_destLayout,
                      const ProblemDomain& a_domain,
                      const IntVect&       a_srcGhostVect,
                      const IntVect&       a_destGhostVect,
                      const IntVect&       a_shift)
{
    BoxLayout srcLayoutGrow;
    srcLayoutGrow.deepCopy(a_srcLayout);
    srcLayoutGrow.grow(a_srcGhostVect);
    srcLayoutGrow.close();

    a_copier.define(srcLayoutGrow,
                    a_destLayout,
                    a_domain,
                    a_destGhostVect,
                    false,    // exchange copier?
                    a_shift);
}


// -----------------------------------------------------------------------------
void
defineImpartialCopier(StaggeredCopier&     a_copier,
                      const int            a_fcDir,
                      const BoxLayout&     a_srcLayout,
                      const BoxLayout&     a_destLayout,
                      const ProblemDomain& a_domain,
                      const IntVect&       a_srcGhostVect,
                      const IntVect&       a_destGhostVect,
                      const IntVect&       a_shift)
{
    BoxLayout srcLayoutGrow;
    srcLayoutGrow.deepCopy(a_srcLayout);
    srcLayoutGrow.grow(a_srcGhostVect);
    srcLayoutGrow.close();

    a_copier.define(srcLayoutGrow,
                    a_destLayout,
                    a_domain,
                    a_destGhostVect,
                    a_fcDir,
                    a_shift);
}


// -----------------------------------------------------------------------------
// This will define a copier that does not care about valid vs invalid data -
// it will just copy everything.
// -----------------------------------------------------------------------------
void
defineImpartialCopier(std::array<StaggeredCopier, CH_SPACEDIM>& a_copier,
                      const BoxLayout&                          a_srcLayout,
                      const BoxLayout&                          a_destLayout,
                      const ProblemDomain&                      a_domain,
                      const IntVect&                            a_srcGhostVect,
                      const IntVect&                            a_destGhostVect,
                      const IntVect&                            a_shift)
{
    BoxLayout srcLayoutGrow;
    srcLayoutGrow.deepCopy(a_srcLayout);
    srcLayoutGrow.grow(a_srcGhostVect);
    srcLayoutGrow.close();

    for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
        a_copier[fcDir].define(srcLayoutGrow,
                               a_destLayout,
                               a_domain,
                               a_destGhostVect,
                               fcDir,
                               a_shift);
    }
}


// -----------------------------------------------------------------------------
std::vector<Real>
domainDecompCoordVec(const std::vector<Real>& a_globalCoords,
                     const int                a_dir,
                     Box                      a_valid,
                     const Box&               a_domBox)
{
    CH_assert(0 <= a_dir && a_dir < SpaceDim);
    CH_assert(a_domBox.contains(a_valid));
    CH_verify(a_globalCoords.size() == static_cast<size_t>(a_domBox.size(a_dir)));

    // In what follows, g = global index and l = local index.

    // These casts are OK because a_valid is shifted.
    a_valid.shift(-a_domBox.smallEnd());
    const size_t gmin = static_cast<size_t>(a_valid.smallEnd(a_dir));
    const size_t gmax = static_cast<size_t>(a_valid.bigEnd(a_dir));

    std::vector<Real> localCoords(gmax - gmin + 1);
    for (size_t g = gmin, l = 0; g <= gmax; ++g, ++l) {
        localCoords[l] = a_globalCoords[g];
    }

    return localCoords;
}


// *****************************************************************************
// The GrowTransform class
//
// Grows each Box in a BoxLayout while preserving proc assignments, etc.
// *****************************************************************************

// -----------------------------------------------------------------------------
GrowTransform::GrowTransform(const IntVect& a_grow)
: m_grow(a_grow)
{
}

// -----------------------------------------------------------------------------
Box
GrowTransform::operator()(const Box& a_inputBox)
{
    return Box(a_inputBox).grow(m_grow);
}


// *****************************************************************************
// The RangeTransform class
//
// Changes the range in a single direction of all boxes in a BoxLayout
// while preserving proc assignments, etc.
// *****************************************************************************

// -----------------------------------------------------------------------------
RangeTransform::RangeTransform (const int a_dir,
                                const int a_loEnd,
                                const int a_numCells)
: m_dir(a_dir),
  m_loEnd(a_loEnd),
  m_numCells(a_numCells)
{
    CH_assert(0 <= m_dir);
    CH_assert(m_dir < SpaceDim);
    CH_assert(m_numCells >= 1);
}


// -----------------------------------------------------------------------------
RangeTransform::~RangeTransform ()
{;}


// -----------------------------------------------------------------------------
Box
RangeTransform::operator() (const Box& a_inputBox)
{
    return Box(a_inputBox).setRange(m_dir, m_loEnd, m_numCells);
}


// *****************************************************************************
// The ShiftTransform class
//
// Shifts all boxes in a BoxLayout while preserving proc assignments, etc.
// *****************************************************************************

// -----------------------------------------------------------------------------
ShiftTransform::ShiftTransform (const IntVect& a_shift)
: m_shift(a_shift)
{;}


// -----------------------------------------------------------------------------
ShiftTransform::ShiftTransform (const int a_dir,
                                const int a_shift)
: m_shift(a_shift * BASISV(a_dir))
{;}


// -----------------------------------------------------------------------------
ShiftTransform::~ShiftTransform ()
{;}


// -----------------------------------------------------------------------------
Box
ShiftTransform::operator() (const Box& a_inputBox)
{
    return Box(a_inputBox).shift(m_shift);
}


// *****************************************************************************
// The ShiftHalfTransform class
//
// Shifts all boxes in a BoxLayout while preserving proc assignments, etc.
// *****************************************************************************

// -----------------------------------------------------------------------------
ShiftHalfTransform::ShiftHalfTransform (const IntVect& a_shift)
: m_shift(a_shift)
{;}


// -----------------------------------------------------------------------------
ShiftHalfTransform::ShiftHalfTransform (const int a_dir,
                                        const int a_shift)
: m_shift(a_shift * BASISV(a_dir))
{;}


// -----------------------------------------------------------------------------
ShiftHalfTransform::~ShiftHalfTransform ()
{;}


// -----------------------------------------------------------------------------
Box
ShiftHalfTransform::operator() (const Box& a_inputBox)
{
    return Box(a_inputBox).shiftHalf(m_shift);
}


// *****************************************************************************
// The SliceTransform class
//
// Slices all boxes in a BoxLayout while preserving proc assignments, etc.
// *****************************************************************************

// -----------------------------------------------------------------------------
SliceTransform::SliceTransform (const int a_dir,
                                const int a_pos)
: m_dir(a_dir),
  m_pos(a_pos)
{;}


// -----------------------------------------------------------------------------
SliceTransform::~SliceTransform ()
{;}


// -----------------------------------------------------------------------------
Box
SliceTransform::operator() (const Box& a_inputBox)
{
    Box retBox = a_inputBox;
    retBox.setBig(m_dir, a_inputBox.smallEnd(m_dir));
    retBox.shift(m_dir, m_pos - a_inputBox.smallEnd(m_dir));
    CH_assert(!retBox.isEmpty());
    return retBox;
}


// *****************************************************************************
// The SwapTransform class
// *****************************************************************************

// -----------------------------------------------------------------------------
SwapTransform::SwapTransform(const std::vector<Box>& a_oldBoxes,
                             const std::vector<Box>& a_newBoxes)
: m_oldBoxes(a_oldBoxes), m_newBoxes(a_newBoxes)
{
    // If the a_oldBoxes was pre-sorted, this should be fast.
    const auto p = sortPermutation(m_oldBoxes);
    applyPermutation(m_oldBoxes, p);
    applyPermutation(m_newBoxes, p);
}


// -----------------------------------------------------------------------------
Box
SwapTransform::operator()(const Box& a_inputBox)
{
    const auto lower =
        std::lower_bound(m_oldBoxes.begin(), m_oldBoxes.end(), a_inputBox);

    if (lower == m_oldBoxes.end()) {
        MAYDAYERROR(
            "SwapTransform is being asked to transform an unrecognized Box. Make "
            "sure you are trying to transform a BoxLayout whose Boxes were "
            "identical to those provided in a_oldBoxes.\n\tinputBox = "
            << a_inputBox << "\n\toldBoxes = " << m_oldBoxes);
    }

    return m_newBoxes[std::distance(m_oldBoxes.begin(), lower)];
}


// *****************************************************************************
// Face overlap tools.
// *****************************************************************************

// -----------------------------------------------------------------------------
/// Identifies overlapping faces that should have identical data.
class ValidFaceOverlapCopier : public StaggeredCopier
{
public:
    ValidFaceOverlapCopier() {}

    virtual void
    define(const DisjointBoxLayout& a_grids,
           const int                a_fcDir);

private:
    // This tells C++ that our public define was meant to be an override, not
    // an overload. *These* are the overloads.
    using StaggeredCopier::define;
};


// -----------------------------------------------------------------------------
void
ValidFaceOverlapCopier::define(const DisjointBoxLayout& a_grids,
                               const int                a_fcDir)
{
    CH_assert(0 <= a_fcDir && a_fcDir < SpaceDim);

    this->clear();

    m_isDefined = true;
    m_fcDir     = a_fcDir;

    DataIterator     dit = a_grids.dataIterator();
    NeighborIterator nit(a_grids);
    const int        myRank = procID();

    for (dit.begin(); dit.ok(); ++dit) {
        // 1. Neighbor's boundary --> our boundary
        for (SideIterator sit; sit.ok(); ++sit) {
            const Box myBdry = bdryBox(a_grids[dit], a_fcDir, sit());

            for (nit.begin(dit()); nit.ok(); ++nit) {
                const int nbrRank = a_grids.procID(nit());
                const Box nbrBdry = bdryBox(nit.box(), a_fcDir, flip(sit()));

                if (nbrBdry.intersectsNotEmpty(myBdry)) {
                    const Box overlap(nbrBdry & myBdry);

                    MotionItem* item = new (s_motionItemPool.getPtr())
                        MotionItem(DataIndex(nit()),
                                    dit(),
                                    nit.unshift(overlap),
                                    overlap);

                    if (nbrRank == myRank) {  // local move
                        m_localMotionPlan.push_back(item);
                    } else {
                        item->procID = nbrRank;
                        m_toMotionPlan.push_back(item);
                    }
                }
            } // nit
        } // sit

        // 2. Our boundary --> neighbor's boundary
        for (nit.begin(dit()); nit.ok(); ++nit) {
            const int nbrRank = a_grids.procID(nit());

            // We already created the local plan.
            if (nbrRank == myRank) continue;

            for (SideIterator sit; sit.ok(); ++sit) {
                const Box nbrBdry = bdryBox(nit.box(), a_fcDir, sit());
                const Box myBdry  = bdryBox(a_grids[dit], a_fcDir, flip(sit()));

                if (nbrBdry.intersectsNotEmpty(myBdry)) {
                    const Box overlap(nbrBdry & myBdry);

                    MotionItem* item = new (s_motionItemPool.getPtr())
                        MotionItem(dit(),
                                    DataIndex(nit()),
                                    overlap,
                                    nit.unshift(overlap));

                    item->procID = nbrRank;
                    m_fromMotionPlan.push_back(item);
                }
            } // sit
        } // nit
    } // dit

    this->sort();
}


// -----------------------------------------------------------------------------
/// Checks overlapping faces to ensure data is non-NaN and identical.
class ValidFaceOverlapCheckOp : public LDOperator<FArrayBox>
{
public:
    ValidFaceOverlapCheckOp(const int    a_fcDir,
                            const char*  a_file,
                            const size_t a_line)
    : m_fcDir(a_fcDir), m_file(a_file), m_line(a_line)
    {
        CH_assert(0 <= m_fcDir);
        CH_assert(m_fcDir < SpaceDim);
    }

    /// This typically copies non-local data in a copyTo operation.
    /// Here, this is the function that does the overlap check.
    virtual void
    linearIn(FArrayBox&      arg,
             void*           buf,
             const Box&      R,
             const Interval& comps) const;

    /// This typically takes care of local copies in a copyTo operation.
    /// Here, it just calls linearIn to do the overlap check.
    virtual void
    op(FArrayBox&       dest,
       const Box&       RegionFrom,
       const Interval&  Cdest,
       const Box&       RegionTo,
       const FArrayBox& src,
       const Interval&  Csrc) const;

    /// Thrown when an error is found. Be sure to catch it!
    class overlap_error : public std::exception
    {
    public:
        overlap_error(const std::string a_msg): m_msg(a_msg) {}
        virtual const char* what() const noexcept {return m_msg.c_str();}
    private:
        const std::string m_msg;
    };

protected:
    const int    m_fcDir;
    const char*  m_file;
    const size_t m_line;
};


// -----------------------------------------------------------------------------
void
ValidFaceOverlapCheckOp::linearIn(FArrayBox&      a_arg, // This is the dest FAB.
                                  void*           a_buf,
                                  const Box&      a_R,
                                  const Interval& a_comps) const
{
    Real* buffer = (Real*)a_buf;

    if (a_arg.box().contains(a_R)) {
        for (int comp = a_comps.begin(); comp <= a_comps.end(); ++comp) {
            for (BoxIterator bit(a_R); bit.ok(); ++bit) {
                const IntVect& iv      = bit();
                const Real     srcVal  = *buffer;
                const Real     destVal = a_arg(iv, comp);

                if ((srcVal != srcVal) || (srcVal <= -1.0e10) ||
                    (1.0e10 <= srcVal)) {
                    std::ostringstream msg;
                    msg << "\n[Proc " << procID() << "]: NaN found at " << m_file << ":" << m_line
                        << "\nFluxBox[" << m_fcDir << "](" << iv << ", " << comp
                        << ") = " << Format::pushFlags << Format::scientific
                        << srcVal << Format::popFlags
                        << "\nFluxBox.box() = " << a_arg.box();
                    throw overlap_error(msg.str());
                }

                if ((destVal != destVal) || (destVal <= -1.0e10) ||
                    (1.0e10 <= destVal)) {
                    std::ostringstream msg;
                    msg << "\n[Proc " << procID() << "]: NaN found at " << m_file << ":" << m_line
                        << "\nFluxBox[" << m_fcDir << "](" << iv << ", " << comp
                        << ") = " << Format::pushFlags << Format::scientific
                        << destVal << Format::popFlags
                        << "\nFluxBox.box() = " << a_arg.box();
                    throw overlap_error(msg.str());
                }

                if (abs(srcVal - destVal) > 1.0e-12) {
                    std::ostringstream msg;
                    msg << "\n[Proc " << procID() << "]: Overlapping faces contain inconsistent data at "
                        << m_file << ":" << m_line << "\nFluxBox[" << m_fcDir
                        << "](" << iv << ", " << comp << ") = "
                        << Format::pushFlags << Format::scientific
                        << srcVal << " or " << destVal
                        << "\ndiff = " << std::abs(destVal - srcVal)
                        << Format::popFlags
                        << "\nFluxBox.box() = " << a_arg.box();
                    throw overlap_error(msg.str());
                }

                ++buffer;
            } // bit
        } // comp
    } else {
        MayDay::Error("a_arg does not contain a_R.");
    }
}


// -----------------------------------------------------------------------------
void
ValidFaceOverlapCheckOp::op(FArrayBox&       dest,
                            const Box&       RegionFrom,
                            const Interval&  Cdest,
                            const Box&       RegionTo,
                            const FArrayBox& src,
                            const Interval&  Csrc [[maybe_unused]]) const
{
    const IntVect toShift = RegionFrom.smallEnd() - RegionTo.smallEnd();
    dest.shift(toShift);
    const_cast<Box&>(RegionTo).shift(toShift);

    CH_assert(RegionFrom.type() == BASISV(m_fcDir));
    CH_assert(RegionFrom.size(m_fcDir) == 1);

    CH_assert(RegionTo.type() == BASISV(m_fcDir));
    CH_assert(RegionTo.size(m_fcDir) == 1);

    CH_assert(RegionFrom == RegionTo);

    CH_assert(dest.box().type() == BASISV(m_fcDir));
    CH_assert(dest.box().contains(RegionFrom));

    CH_assert(src.box().type() == BASISV(m_fcDir));
    CH_assert(src.box().contains(RegionTo));

    CH_assert(Csrc == Cdest);

    for (int comp = Cdest.begin(); comp <= Cdest.end(); ++comp) {
        for (BoxIterator bit(RegionTo); bit.ok(); ++bit) {
            const IntVect& iv      = bit();
            const Real     srcVal  = src(iv, comp);
            const Real     destVal = dest(iv, comp);

            if ((srcVal != srcVal) || (srcVal <= -1.0e10) ||
                (1.0e10 <= srcVal)) {

                std::ostringstream msg;
                msg << "\n[Proc " << procID() << "]: NaN found at " << this->m_file << ":" << m_line
                    << "\nFluxBox[" << m_fcDir << "](" << iv << ", " << comp
                    << ") = " << Format::pushFlags << Format::scientific
                    << srcVal << Format::popFlags
                    << "\nFluxBox.box() = " << dest.box();
                throw overlap_error(msg.str());
            }

            if ((destVal != destVal) || (destVal <= -1.0e10) ||
                (1.0e10 <= destVal)) {
                std::ostringstream msg;
                msg << "\n[Proc " << procID() << "]: NaN found at " << this->m_file << ":" << m_line
                    << "\nFluxBox[" << m_fcDir << "](" << iv << ", " << comp
                    << ") = " << Format::pushFlags << Format::scientific
                    << destVal << Format::popFlags
                    << "\nFluxBox.box() = " << dest.box();
                throw overlap_error(msg.str());
            }

            if (abs(srcVal - destVal) > 1.0e-12) {
                std::ostringstream msg;
                msg << "\n[Proc " << procID() << "]: Overlapping faces contain inconsistent data at "
                    << m_file << ":" << m_line << "\nFluxBox[" << m_fcDir
                    << "](" << iv << ", " << comp << ") = "
                    << Format::pushFlags << Format::scientific
                    << srcVal << " or " << destVal
                    << "\ndiff = " << std::abs(destVal - srcVal)
                    << Format::popFlags
                    << "\nFluxBox.box() = " << dest.box();
                throw overlap_error(msg.str());
            }
        }  // bit
    } // comp

    const_cast<Box&>(RegionTo).shift(-toShift);
    dest.shift(-toShift);
}


// -----------------------------------------------------------------------------
void
checkValidFaceOverlap(const LevelData<FluxBox>& a_data,
                      const char*               a_file,
                      const size_t              a_line)
{
    for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
        FABAliasFlBxDataFactory factory(
            const_cast<LevelData<FluxBox>*>(&a_data), a_data.interval(), fcDir);
        LevelData<FArrayBox> dataComp(
            a_data.getBoxes(), a_data.nComp(), a_data.ghostVect(), factory);

        checkValidFaceOverlap(dataComp, fcDir, a_file, a_line);
    }

    // const auto& grids = a_data.getBoxes();

    // for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
    //     LevelData<FArrayBox> ccData(grids, 1, BASISV(fcDir));

    //     for (DataIterator dit(grids); dit.ok(); ++dit) {
    //         ccData[dit].setVal(0.0);

    //         for (SideIterator sit; sit.ok(); ++sit) {
    //             const Box bx = bdryBox(grids[dit], fcDir, sit());

    //             ccData[dit].shiftHalf(fcDir, sign(sit()));
    //             ccData[dit].copy(a_data[dit][fcDir], bx);

    //             ccData[dit].shift(fcDir, -sign(sit()));
    //             ccData[dit].copy(a_data[dit][fcDir], bx);

    //             ccData[dit].shiftHalf(fcDir, sign(sit()));
    //         }
    //     }

    //     Copier cp;
    //     cp.exchangeDefine(grids, ccData.ghostVect());
    //     cp.trimEdges(grids, ccData.ghostVect());
    //     ccData.exchange(cp);

    //     for (DataIterator dit(grids); dit.ok(); ++dit) {
    //         for (SideIterator sit; sit.ok(); ++sit) {
    //             const Box fcBdryBox = bdryBox(grids[dit], fcDir, sit(), 1);
    //             const int isign = sign(sit());

    //             FArrayBox errorFAB(fcBdryBox, a_data.nComp());

    //             ccData[dit].shiftHalf(fcDir, isign);
    //             errorFAB.copy(ccData[dit], fcBdryBox);
    //             ccData[dit].shift(fcDir, -isign);
    //             errorFAB.minus(ccData[dit], fcBdryBox, 0, 0, a_data.nComp());
    //             ccData[dit].shiftHalf(fcDir, isign);

    //             const auto maxError = errorFAB.norm(0);
    //             if ((maxError != maxError) || (1.0e10 <= maxError)) {
    //                 MAYDAYERROR("\n[Proc " << procID() << "]: NaN found at "
    //                                        << a_file << ":" << a_line);
    //             }
    //             if (maxError > 1.0e-12) {
    //                 MAYDAYERROR(
    //                     "\n[Proc "
    //                     << procID()
    //                     << "]: Overlapping faces contain inconsistent data at "
    //                     << a_file << ":" << a_line << "  |diff| = " << maxError);
    //             }
    //         } // sit
    //     } // dit
    // } // fcDir
}


// -----------------------------------------------------------------------------
void
checkValidFaceOverlap(const LevelData<FArrayBox>& a_data,
                      const int                   a_fcDir,
                      const char*                 a_file,
                      const size_t                a_line)
{
    CH_assert(0 <= a_fcDir && a_fcDir < SpaceDim);

    ValidFaceOverlapCheckOp op(a_fcDir, a_file, a_line);

    ValidFaceOverlapCopier cp;
    cp.define(a_data.getBoxes(), a_fcDir);

    try {
        a_data.copyTo(const_cast<LevelData<FArrayBox>&>(a_data), cp, op);
        // const DisjointBoxLayout& grids = a_data.getBoxes();
        // LevelData<FArrayBox> tmp(grids,
        //                          a_data.nComp(),
        //                          IntVect::Zero,
        //                          NCDataFactory<FArrayBox>(BASISV(a_fcDir)));
        // for (DataIterator dit(grids); dit.ok(); ++dit) {
        //     tmp[dit].copy(a_data[dit], grids[dit].surroundingNodes(a_fcDir));
        // }
        // tmp.copyTo(tmp, cp, op);

    } catch (const ValidFaceOverlapCheckOp::overlap_error& e) {
        pout() << e.what() << std::endl;
        MayDay::Error(e.what());

    } catch (...) {
        std::ostringstream msg;
        msg << "[Proc " << procID() << "]: checkValidFaceOverlap caught an unexpected exception "
            << "at " << a_file << ": " << a_line << std::endl;
        pout() << msg.str().c_str();
        MayDay::Error(msg.str().c_str());
    }
}


// -----------------------------------------------------------------------------
void
checkValidFaceOverlap(const LevelData<FluxBox>& a_data,
                      LevelData<FluxBox>&       a_destDiscrepancies)
{
    const auto& grids = a_data.getBoxes();

    a_destDiscrepancies.define(grids, a_data.nComp());
    for (DataIterator dit(grids); dit.ok(); ++dit) {
        a_destDiscrepancies[dit].setVal(quietNAN);
    }

    for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
        LevelData<FArrayBox> ccData(grids, 1, BASISV(fcDir));

        for (DataIterator dit(grids); dit.ok(); ++dit) {
            ccData[dit].setVal(0.0);

            for (SideIterator sit; sit.ok(); ++sit) {
                const Box bx = bdryBox(grids[dit], fcDir, sit());

                ccData[dit].shiftHalf(fcDir, sign(sit()));
                ccData[dit].copy(a_data[dit][fcDir], bx);

                ccData[dit].shift(fcDir, -sign(sit()));
                ccData[dit].copy(a_data[dit][fcDir], bx);

                ccData[dit].shiftHalf(fcDir, sign(sit()));
            }
        }

        Copier cp;
        cp.exchangeDefine(grids, ccData.ghostVect());
        cp.trimEdges(grids, ccData.ghostVect());
        ccData.exchange(cp);

        for (DataIterator dit(grids); dit.ok(); ++dit) {
            for (SideIterator sit; sit.ok(); ++sit) {
                const Box fcBdryBox = bdryBox(grids[dit], fcDir, sit(), 1);
                const int isign = sign(sit());

                ccData[dit].shiftHalf(fcDir, isign);
                a_destDiscrepancies[dit][fcDir].copy(ccData[dit], fcBdryBox);
                ccData[dit].shift(fcDir, -isign);
                a_destDiscrepancies[dit][fcDir].minus(ccData[dit], fcBdryBox, 0, 0, a_data.nComp());
                ccData[dit].shiftHalf(fcDir, isign);
            }
        }
    } // fcDir
}


// -----------------------------------------------------------------------------
void
averageOverlappingValidFaces(LevelData<FluxBox>& a_data)
{
    checkForValidNAN(a_data);

    const DisjointBoxLayout& grids = a_data.getBoxes();
    DataIterator             dit   = a_data.dataIterator();

    for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
        LevelData<FArrayBox> ccData(grids, 1, BASISV(fcDir));

        for (dit.reset(); dit.ok(); ++dit) {
            ccData[dit].setVal(0.0);

            for (SideIterator sit; sit.ok(); ++sit) {
                const Box bx = bdryBox(grids[dit], fcDir, sit());

                ccData[dit].shiftHalf(fcDir, sign(sit()));
                ccData[dit].copy(a_data[dit][fcDir], bx);

                ccData[dit].shift(fcDir, -sign(sit()));
                ccData[dit].copy(a_data[dit][fcDir], bx);

                ccData[dit].shiftHalf(fcDir, sign(sit()));
            }
        }

        Copier cp;
        cp.exchangeDefine(grids, ccData.ghostVect());
        cp.trimEdges(grids, ccData.ghostVect());
        ccData.exchange(cp);

        for (dit.reset(); dit.ok(); ++dit) {
            for (SideIterator sit; sit.ok(); ++sit) {
                const Box bx = bdryBox(grids[dit], fcDir, sit());
                Convert::Simple(a_data[dit][fcDir], 0, bx, ccData[dit], 0);
            }
        }
    }

    checkForValidNAN(a_data);
    debugCheckValidFaceOverlap(a_data);
}


}; // end LayoutTools namespace
