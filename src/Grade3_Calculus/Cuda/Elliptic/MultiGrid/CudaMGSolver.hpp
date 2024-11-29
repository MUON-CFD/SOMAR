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
#pragma once
#include "FunctionTimer.hpp"
#include "IntVect.H"
#include "isLinearOperatorOverHilbertMGCompatible.hpp"
#include <vector>

namespace Elliptic {
namespace MultiGrid {
// =============================== MGSolver ====================================

/**
 * MGSolver class, generic. It requires a LinearOperator over a Hilbert class,
 * a BottomSolver class to handle the solver at the coarsest level, and a BC
 * class to provide the boundary conditions.
 * See Elliptic/traits/isLinearOperatorOverHilbertMGCompatible.hpp for a
 * list of required members.
 */
template <typename Hilbert, typename LinearOperator, typename BottomSolver,
          typename BC>
class MGSolver {
  static_assert(Elliptic::traits::isLinearOperatorOverHilbertMGCompatible<
                    LinearOperator, Hilbert>::value,
                "The linear operator failed \
to provide the functions required.\
See isLinearOperatorOverHilbertMGCompatible.hpp for more details");

public:
  /// Knobs for you to turn.
  struct Options {
    Options() {} // A bug in the C++11 standard requires this.

    Real absTol = 1.0e-16; //!< Stop if |res| < absTol
    Real relTol = 1.0e-14; //!< Stop if |res| / |init res| < relTol
    int numSmoothDown = 4;
    int numSmoothUp = 4;
    int numSmoothBottom = 2;
    int numSmoothPrecond = 2; // Is this needed?
    int prolongOrderFMG = 3;
    int numSmoothUpFMG = 0;
    int maxDepth = -1; // -1 = coarsen as much as possible
    int numCycles = 1;
    int maxIters = 10;
    int normType = 2; // 0 = inf-norm, etc...
    int verbosity = 0;

    // AMRNSLevel::ProjectorBC phiBC;
  };

  /// Default constructor. Leaves object unusable.
  MGSolver();

  /// Destructor
  ~MGSolver();

  //! a_opPtr must be well-defined at depth = 0.
  //! If you have something up your sleeve, you can provide your own refRatio
  //! at each depth. Otherwise, SemicoarseningStrategy will be used.
  void define(std::vector<std::shared_ptr<LinearOperator>> &a_opPtr,
              std::vector<IntVect> a_refSchedule, Options a_opt = Options());

  /// Sinmple check.
  bool isDefined() const;

  /// Free's memory. Leaves object unusable.
  void clear();

  //! returns pointer to operator at level d
  //! defaults to level 0 if d not provided
  LinearOperator *getOpPtr(int d = 0) const { return m_opPtrs[d].get(); }

  //! returns hang to Bottom solver
  BottomSolver *getBottomSolver() const { return m_bottomSolverPtr; }

  //! Solves L[phi] = rhs using the most appropriate methods.
  //! returns a vector storing the convergence pattern.
  auto solve(Hilbert &a_phi, const Hilbert &a_rhs, const Real a_time,
             const bool a_useHomogBCs = false,
             const bool a_setPhiToZero = false) -> std::vector<Real>;

  //! Solves L[phi] = rhs using V-cycles.
  std::vector<Real> vCycle(Hilbert &a_phi, const Hilbert &a_rhs, const Real a_time,
              const bool a_useHomogBCs = false,
              const bool a_setPhiToZero = false);

  //! Solves L[phi] = rhs using Full Multigrid.
  std::vector<Real> fmg(Hilbert &a_phi, const Hilbert &a_rhs, const Real a_time,
           const bool a_useHomogBCs = false, const bool a_setPhiToZero = false);

protected:
  //! pout formatting utility.
  void indent(const int a_depth = -1) const;

  //! pout formatting utility.
  void unindent(const int a_depth = -1) const;

  //! Solves L[cor] = res with homogBCs at a specified depth.
  void vCycle_noInit(Hilbert &a_cor, const Hilbert &a_res, const Real a_time,
                     const int a_depth);

  //! Solves L[cor] = res with homogBCs at a specified depth.
  void fmg_noInit(Hilbert &a_cor, const Hilbert &a_res, const Real a_time,
                  const int a_depth);
//! Returns a pointer to the operator at level, defaults to level 0 if not specified.
  LinearOperator *getOp(int level = 0) const { return getOpPtr(level); }
  // Member variables...
  Vector<IntVect> m_refSchedule;
  bool m_isDefined;
  std::vector<std::shared_ptr<LinearOperator>> m_opPtrs;
  BottomSolver *m_bottomSolverPtr;

  Options m_opt;

private:
  MGSolver(const MGSolver &source) {}
  MGSolver &operator=(const MGSolver &that) {}
};

#include "CudaMGSolverImpl.hpp"

} // namespace MultiGrid
} // namespace Elliptic
