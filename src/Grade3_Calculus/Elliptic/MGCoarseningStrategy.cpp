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
#include "MGCoarseningStrategy.H"
#include "AnisotropicRefinementTools.H"
#include "Format.H"
#include "Debug.H"

namespace Elliptic {

#if CH_SPACEDIM == 2
    const std::vector<IntVect> SemicoarseningStrategy::s_refList{
        IntVect(2, 1),
        IntVect(1, 2),
        IntVect(2, 2),
    };
#else
    const std::vector<IntVect> SemicoarseningStrategy::s_refList{
        IntVect(2, 1, 1),
        IntVect(1, 2, 1),
        IntVect(1, 1, 2),
        IntVect(1, 2, 2),
        IntVect(2, 1, 2),
        IntVect(2, 2, 1),
        IntVect(2, 2, 2)
    };
#endif

// -----------------------------------------------------------------------------
// Full constructor
// a_L is the length of the computational domain.
// -----------------------------------------------------------------------------
SemicoarseningStrategy::SemicoarseningStrategy(const RealVect& a_L,
                                               const int       a_verbosity)
: m_L(a_L)
, m_verbosity(a_verbosity)
{
    CH_assert(m_L > RealVect::Zero);
}


// -----------------------------------------------------------------------------
// Returns the MG refinement schedule.
// If m_opt.maxDepth = -1, then this function will coarsen as much as possible.
// -----------------------------------------------------------------------------
Vector<IntVect>
SemicoarseningStrategy::createMGRefSchedule(const DisjointBoxLayout& a_grids,
                                            const IntVect& a_minBoxSize,
                                            const int      a_maxDepth)
{
    DisjointBoxLayout curGrids;
    RealVect          curDx = m_L / RealVect(a_grids.physDomain().size());
    Vector<IntVect>   mgRefSchedule(0);
    Vector<RealVect>  mgDx(1, curDx);  // For debugging purposes.
    std::string       causeOfTermination;

    // We start by coarsening by minBoxSize. This ensures that the boxes at
    // all depths are at least minBoxSize.
    if (!coarsenable(a_grids, a_minBoxSize)) {
        MAYDAYERROR(
            "Received a_grids that are not coarsenable by a_minBoxSize.");
    }
    coarsen(curGrids, a_grids, a_minBoxSize);

    // Now we coarsen and coarsen...
    do {
        // Have we reached the maximum allowable depth?
        if (a_maxDepth >= 0 && mgRefSchedule.size() == (size_t)a_maxDepth) {
            causeOfTermination = std::string("maxDepth reached.");
            break;
        }

        // pout() << "\ndepth = " << mgRefSchedule.size() << Format::indent() << endl;
        const IntVect mgRefRatio = computeNextRefRatio(curDx, curGrids);
        // pout() << Format::unindent << endl;

        // Can the grids be coarsened?
        if (mgRefRatio == IntVect::Unit) {
            causeOfTermination = std::string("curGrids not coarsenable.");
            break;
        }

        // If we get here, then we can coarsen by mgRefRatio. So do it!
        DisjointBoxLayout tmpGrids;
        coarsen(tmpGrids, curGrids, mgRefRatio);
        curGrids = tmpGrids;
        curDx *= RealVect(mgRefRatio);

        // Add to schedule.
        mgRefSchedule.push_back(mgRefRatio);
        mgDx.push_back(curDx);

    } while (true);  // Muahhahah!!!

    // When we get here, the schedule has been created. In order to make some
    // looping operations simpler, I will push one more refRatio at the end
    // to serve as a termination flag.
    mgRefSchedule.push_back(IntVect::Unit);

    // Dump schedule to pout.
    if (m_verbosity > 2) {
        static std::vector<RealVect> writtenDx(0);
        if (std::find(writtenDx.begin(), writtenDx.end(), curDx) == writtenDx.end()) {
            pout() << "Semicoarsening mgRefSchedule:" << Format::indent();
            for (size_t d = 0; d < mgRefSchedule.size(); ++d) {
                pout() << "\nDepth = " << d << ",\tref = " << mgRefSchedule[d]
                    << ",\tdx = " << mgDx[d]
                    << ",\tiso = " << measureAnisotropy(mgDx[d]);
            }
            pout() << "\ncauseOfTermination = " << causeOfTermination << endl;
            pout() << Format::unindent << '\n' << std::endl;

            writtenDx.push_back(curDx);
        }
    }

    // Done!
    return mgRefSchedule;
}

// -----------------------------------------------------------------------------
IntVect
SemicoarseningStrategy::computeNextRefRatio(const RealVect&          a_dx,
                                            const DisjointBoxLayout& a_grids)
{
    constexpr Real isoLimiter = 0.75;

    size_t     isoIdx     = std::numeric_limits<size_t>::max();
    Real       isoVal     = measureAnisotropy(a_dx);
    const Real lastIsoVal = isoVal;

    for (size_t curIdx = 0; curIdx < s_refList.size(); ++curIdx) {
        const IntVect& curRef    = s_refList[curIdx];
        const Real     curIsoVal = measureAnisotropy(a_dx * curRef);
        // pout() << curRef << ": " << curIsoVal;

        // Is this amount of coarsening allowed? If not, skip.
        if (!coarsenable(a_grids, curRef)) {
            // pout() << " not coarsenable with this ref\n";
            continue;
        }

        // Skip if this reduces the iso too much AND we already found a
        // candidate ref.
        if (curIsoVal < isoLimiter * lastIsoVal && isoIdx < s_refList.size()) {
            // pout() << " skipped\n";
            continue;
        }

        // Accept if this reduces anisotropy.
        if (curIsoVal <= isoVal) {
            isoIdx = curIdx;
            isoVal = curIsoVal;
            // pout() << " Best so far";
        }
        // pout() << '\n';
    }

    // Report valid result or return (1,1,1).
    if (isoIdx < s_refList.size()) return s_refList[isoIdx];
    return IntVect::Unit;
}


// =============================================================================
#if CH_SPACEDIM == 2
    const std::vector<IntVect> HorizCoarseningStrategy::s_refList{
        IntVect(2, 1)
    };
#else
    const std::vector<IntVect> HorizCoarseningStrategy::s_refList{
        IntVect(2, 1, 1),
        IntVect(1, 2, 1),
        IntVect(2, 2, 1)
    };
#endif

// -----------------------------------------------------------------------------
HorizCoarseningStrategy::HorizCoarseningStrategy(
    const RealVect& a_L,
    const bool      a_doVertCoarsening,
    const int       a_verbosity)
: m_L(a_L)
, m_doVertCoarsening(a_doVertCoarsening)
, m_verbosity(a_verbosity)
{
    CH_assert(m_L > RealVect::Zero);
}


// -----------------------------------------------------------------------------
Vector<IntVect>
HorizCoarseningStrategy::createMGRefSchedule(
    const DisjointBoxLayout& a_grids,
    const IntVect&           a_minBoxSize,
    const int                a_maxDepth)
{
    DisjointBoxLayout curGrids;
    RealVect          curDx = m_L / RealVect(a_grids.physDomain().size());
    Vector<IntVect>   mgRefSchedule(0);
    Vector<RealVect>  mgDx(1, curDx);  // For debugging purposes.
    std::string       causeOfTermination;

    // We start by coarsening by minBoxSize. This ensures that the boxes at
    // all depths are at least minBoxSize.
    if (!coarsenable(a_grids, a_minBoxSize)) {
        MAYDAYERROR(
            "Received a_grids that are not coarsenable by a_minBoxSize.");
    }
    coarsen(curGrids, a_grids, a_minBoxSize);

    // Now we coarsen and coarsen...
    do {
        // Have we reached the maximum allowable depth?
        if (a_maxDepth >= 0 && mgRefSchedule.size() == (size_t)a_maxDepth) {
            causeOfTermination = std::string("maxDepth reached.");
            break;
        }

        // The vertical ref will be 1.
        IntVect mgRefRatio = computeNextRefRatio(curDx, curGrids);

        // Can the grids be coarsened?
        if (mgRefRatio == IntVect::Unit) {
            causeOfTermination = std::string("curGrids not coarsenable.");
            break;
        }

        // Now, decide if we should make the vertical ref 2.
        if (m_doVertCoarsening) {
            if (coarsenable(curGrids, IntVect::Unit + BASISV(SpaceDim - 1))) {
                mgRefRatio[SpaceDim - 1] = 2;
            }
        }

        // If we get here, then we can coarsen by mgRefRatio. So do it!
        DisjointBoxLayout tmpGrids;
        coarsen(tmpGrids, curGrids, mgRefRatio);
        curGrids = tmpGrids;
        curDx *= RealVect(mgRefRatio);

        // Add to schedule.
        mgRefSchedule.push_back(mgRefRatio);
        mgDx.push_back(curDx);

    } while (true);  // Muahhahah!!!

    // When we get here, the schedule has been created. In order to make some
    // looping operations simpler, I will push one more refRatio at the end
    // to serve as a termination flag.
    mgRefSchedule.push_back(IntVect::Unit);

    // Dump schedule to pout.
    if (m_verbosity > 4) {
        pout() << "HorizCoarseningStrategy mgRefSchedule:" << Format::indent();
        for (size_t d = 0; d < mgRefSchedule.size(); ++d) {
            pout() << "\nDepth = " << d << ",\tref = " << mgRefSchedule[d]
                   << ",\tdx = " << mgDx[d]
                   << ",\tiso = " << measureAnisotropy(mgDx[d]);
        }
        pout() << "\ncauseOfTermination = " << causeOfTermination << endl;
        pout() << Format::unindent << "\n" << std::endl;
    }

    // Done!
    return mgRefSchedule;
}

// -----------------------------------------------------------------------------
IntVect
HorizCoarseningStrategy::computeNextRefRatio(const RealVect&          a_dx,
                                             const DisjointBoxLayout& a_grids)
{
    constexpr Real isoLimiter = 0.75;

    size_t     isoIdx     = std::numeric_limits<size_t>::max();
    Real       isoVal     = measureAnisotropy(a_dx);
    const Real lastIsoVal = isoVal;

    for (size_t curIdx = 0; curIdx < s_refList.size(); ++curIdx) {
        const IntVect& curRef    = s_refList[curIdx];
        const Real     curIsoVal = measureAnisotropy(a_dx * curRef);

        // Is this amount of coarsening allowed? If not, skip.
        if (!coarsenable(a_grids, curRef)) continue;

        // Skip if this reduces the iso too much AND we already found a candidate ref.
        if (curIsoVal < isoLimiter * lastIsoVal && isoIdx < s_refList.size())
            continue;

        // Accept if this reduces anisotropy.
        if (curIsoVal <= isoVal) {
            isoIdx = curIdx;
            isoVal = curIsoVal;
        }
    }

    // Report valid result or return (1,1,1).
    if (isoIdx < s_refList.size()) return s_refList[isoIdx];
    return IntVect::Unit;
}





// -----------------------------------------------------------------------------
// Full constructor
// a_L is the length of the computational domain.
// -----------------------------------------------------------------------------
MiniVCycleStrategy::MiniVCycleStrategy(const RealVect& a_L,
                                       const IntVect&  a_amrRefRatio,
                                       const int       a_verbosity)
: m_L(a_L)
, m_amrRefRatio(a_amrRefRatio)
, m_verbosity(a_verbosity)
{
    CH_assert(m_L > RealVect::Zero);
    CH_assert(m_amrRefRatio >= IntVect::Unit);
    CH_assert(m_amrRefRatio.product() >= 2);
}


// -----------------------------------------------------------------------------
// In AMRMGSolvers, you might need to coarsen from l to l-1 where
// the refinement ratio is > 2 in some or all of the directions.
// this is bad for MG. For this purpose, you must perform a mini
// V-cycle before coarsening to l-1. The mini V-Cycle will coarsen
// in such a way that it tries to go to isotropy as fast as possible
// and it terminates at a depth that is one coarsening (by 2) away
// from l-1. In short, we need to eliminate all of the error modes
// not captured by the grid at l-1.

// a_maxDepth will be ignored. We coarsen as much as is needed.
// -----------------------------------------------------------------------------
Vector<IntVect>
MiniVCycleStrategy::createMGRefSchedule(const DisjointBoxLayout& a_grids,
                                        const IntVect&           a_minBoxSize,
                                        const int                /*a_maxDepth*/)
{
    DisjointBoxLayout curGrids;
    RealVect          curDx = m_L / RealVect(a_grids.physDomain().size());
    IntVect           totalCoarsening = IntVect::Unit;
    Vector<IntVect>   mgRefSchedule(0);
    Vector<RealVect>  mgDx(1, curDx);  // For debugging purposes.
    std::string       causeOfTermination;


    // We start by coarsening by minBoxSize. This ensures that the boxes at
    // all depths are at least minBoxSize.
    if (!coarsenable(a_grids, a_minBoxSize)) {
        MAYDAYERROR(
            "Received a_grids that are not coarsenable by a_minBoxSize.");
    }
    coarsen(curGrids, a_grids, a_minBoxSize);

    // Now we coarsen and coarsen...
    do {
        // // Have we reached the maximum allowable depth?
        // if (a_maxDepth >= 0 && mgRefSchedule.size() == (size_t)a_maxDepth) {
        //     causeOfTermination = std::string("maxDepth reached.");
        //     break;
        // }

        // Figure out which directions will benefit from coarsening.
        Real maxDx = 0.0;
        for (int dir = 0; dir < SpaceDim; ++dir) {
            if (totalCoarsening[dir] == m_amrRefRatio[dir]) continue;
            maxDx = max(maxDx, curDx[dir]);
        }

        IntVect mgRefRatio = IntVect::Unit;
        for (int dir = 0; dir < SpaceDim; ++dir) {
            if (totalCoarsening[dir] == m_amrRefRatio[dir]) continue;
            if (curDx[dir] <= maxDx / 2.0) mgRefRatio[dir] = 2;
        }

        // Is anisotropic coarsening worth it?
        // If not, just use the standard mgRefRatio = (2,2,2).
        if (mgRefRatio.product() == 1) {
            for (int dir = 0; dir < SpaceDim; ++dir) {
                if (totalCoarsening[dir] == m_amrRefRatio[dir]) continue;
                mgRefRatio[dir] = 2;
            }
        }

        // Would the mgRefRatio coarsen beyond the coarser AMR level?
        // If so, don't do it!
        {
            IntVect tc = totalCoarsening * mgRefRatio;
            if (tc == m_amrRefRatio) {
                causeOfTermination =
                    std::string("Coarsening would reach crse AMR level.");
                break;
            } else if (!(tc <= m_amrRefRatio)) {
                causeOfTermination =
                    std::string("Coarsening would go past crse AMR level.");
                break;
            }
        }

        // Is this amount of coarsening allowed?
        if (!coarsenable(curGrids, mgRefRatio)) {
            // Nope. Put backup plan here.
            // For now, just cut off MG.
            causeOfTermination = std::string("curGrids not coarsenable.");
            break;
        }

        // If we get here, then we can coarsen by mgRefRatio. So do it!
        DisjointBoxLayout tmpGrids;
        coarsen(tmpGrids, curGrids, mgRefRatio);
        curGrids = tmpGrids;
        curDx *= RealVect(mgRefRatio);
        totalCoarsening *= mgRefRatio;

        // Add to schedule.
        mgRefSchedule.push_back(mgRefRatio);
        mgDx.push_back(curDx);

    } while (true);  // Muahhahah!!!

    // When we get here, the schedule has been created. In order to make some
    // looping operations simpler, I will push one more refRatio at the end
    // to serve as a termination flag.
    mgRefSchedule.push_back(IntVect::Unit);

    // Dump schedule to pout.
    if (m_verbosity > 4) {
        pout() << "Mini V-Cycle with m_amrRefRatio = " << m_amrRefRatio
               << " mgRefSchedule:" << Format::indent();
        for (size_t d = 0; d < mgRefSchedule.size(); ++d) {
            pout() << "\nDepth = " << d << ",\tref = " << mgRefSchedule[d]
                   << ",\tdx = " << mgDx[d];
        }
        pout() << "\ncauseOfTermination = " << causeOfTermination << endl;
        pout() << Format::unindent << "\n" << std::endl;
    }

    // Done!
    return mgRefSchedule;
}


}; // namespace Elliptic
