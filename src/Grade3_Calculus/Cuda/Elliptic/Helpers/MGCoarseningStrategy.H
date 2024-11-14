#pragma once
#include "EllipticSolver.H"
#include "LevelGeometry.H"
#include "MGEllipticOp.H"

#include "AnisotropicCoarseAverage.H"
#include "AnisotropicCoarseAverageF_F.H"
#include "AnisotropicRefinementTools.H"
#include "Comm.H"
#include "Debug.H"
#include "Integral.H"
#include "MGSolver.H"
#include "MGSolverF_F.H"
namespace Elliptic {
namespace Helpers {
// =========================== MGCoarseningStrategy ============================

// -----------------------------------------------------------------------------
// The interface for all coarsening strategies. These choose the amount of
// refinement at each MG depth (the MG refinement schedule).
// -----------------------------------------------------------------------------
class MGCoarseningStrategy {
public:
  // Returns the MG refinement schedule.
  // If maxDepth = -1, then this function will coarsen as much as possible.
  virtual Vector<IntVect>
  createMGRefSchedule(const DisjointBoxLayout &a_grids,
                      const IntVect &a_minBoxSize = IntVect::Unit,
                      const int a_maxDepth = -1) = 0;
};

// -----------------------------------------------------------------------------
// This will create a set of refRatios that push the grids towards
// isotropy (dx = dy = dz) as quickly as possible.
// -----------------------------------------------------------------------------
class SemicoarseningStrategy : public MGCoarseningStrategy {
public:
  // Full constructor
  // a_L is the length of the computational domain.
  SemicoarseningStrategy(const RealVect &a_L) : m_L(a_L) {
  CH_assert(m_L > RealVect::Zero);
}

  // Returns the MG refinement schedule.
  // If maxDepth = -1, then this function will coarsen as much as possible.
  virtual Vector<IntVect>
  createMGRefSchedule(const DisjointBoxLayout &a_grids,
                      const IntVect &a_minBoxSize = IntVect::Unit,
                      const int a_maxDepth = -1)
                      {
  DisjointBoxLayout curGrids;
  RealVect curDx = m_L / RealVect(a_grids.physDomain().size());
  Vector<IntVect> mgRefSchedule(0);
  Vector<RealVect> mgDx(1, curDx); // For debugging purposes.
  std::string causeOfTermination;

  // We start by coarsening by minBoxSize. This ensures that the boxes at
  // all depths are at least minBoxSize.
  if (!coarsenable(a_grids, a_minBoxSize)) {
    MAYDAYERROR("Received a_grids that are not coarsenable by a_minBoxSize.");
  }
  coarsen(curGrids, a_grids, a_minBoxSize);

  do {
    // Have we reached the maximum allowable depth?
    if (a_maxDepth >= 0 && mgRefSchedule.size() == a_maxDepth) {
      causeOfTermination = std::string("maxDepth reached.");
      break;
    }

    // Figure out which directions will benefit from coarsening.
    Real maxDx = 0.0;
    for (int dir = 0; dir < SpaceDim; ++dir) {
      maxDx = max(maxDx, curDx[dir]);
    }

    IntVect mgRefRatio = IntVect::Unit;
    for (int dir = 0; dir < SpaceDim; ++dir) {
      if (curDx[dir] <= maxDx / 2.0)
        mgRefRatio[dir] = 2;
    }

    // Is anisotropic coarsening worth it?
    if (mgRefRatio.product() == 1) {
      mgRefRatio = IntVect(D_DECL(2, 2, 2));
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

    // Add to schedule.
    mgRefSchedule.push_back(mgRefRatio);
    mgDx.push_back(curDx);

  } while (true); // Muahhahah!!!

  // When we get here, the schedule has been created. In order to make some
  // looping operations simpler, I will push one more refRatio at the end
  // to serve as a termination flag.
  mgRefSchedule.push_back(IntVect::Unit);

  // Dump schedule to pout.
  if (false) {
    const int a_maxDepth = mgRefSchedule.size() - 1;
    pout() << "Semicoarsening mgRefSchedule:" << Format::indent();
    for (int d = 0; d <= a_maxDepth; ++d) {
      pout() << "\nDepth = " << d << ",\tref = " << mgRefSchedule[d]
             << ",\tdx = " << mgDx[d];
    }
    pout() << "\ncauseOfTermination = " << causeOfTermination << endl;
    pout() << Format::unindent << "\n" << std::endl;
  }

  // Done!
  return mgRefSchedule;
}

protected:
  RealVect m_L;
};
} // namespace Helpers
} // namespace Elliptic