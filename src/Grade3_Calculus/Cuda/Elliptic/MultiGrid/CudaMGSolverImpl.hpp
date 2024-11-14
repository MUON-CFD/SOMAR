
/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2019
 *    Jefferson University and
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
// -----------------------------------------------------------------------------
// Default constructor. Leaves object unusable.
// -----------------------------------------------------------------------------
template <typename Hilbert, typename LinearOperator, typename BottomSolver,
          typename BC>
MGSolver<Hilbert, LinearOperator, BottomSolver, BC>::MGSolver()
    : m_refSchedule(std::vector<IntVect>()),
      m_opPtrs(std::vector<std::shared_ptr<LinearOperator>>()),
      m_bottomSolverPtr(nullptr),
      m_opt(),
      m_isDefined(false) {}

// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
template <typename Hilbert, typename LinearOperator, typename BottomSolver,
          typename BC>
MGSolver<Hilbert, LinearOperator, BottomSolver, BC>::~MGSolver() {
  this->clear();
}

// -----------------------------------------------------------------------------
// a_opPtr must be well-defined at depth = 0.
// If you have something up your sleeve, you can provide your own refRatio
// at each depth. Otherwise, SemicoarseningStrategy will be used.
// -----------------------------------------------------------------------------
template <typename Hilbert, typename LinearOperator, typename BottomSolver,
          typename BC>
void MGSolver<Hilbert, LinearOperator, BottomSolver, BC>::define(
    std::vector<std::shared_ptr<LinearOperator>> &a_opPtrs,
    std::vector<IntVect> a_refSchedule, Options a_opt) {
  CH_assert(a_opPtrs.size() > 0 && a_opPtrs[0] != nullptr);

  m_opt = a_opt;

  m_refSchedule = a_refSchedule;

  CH_assert(m_refSchedule.size() > 0);
  CH_assert(m_refSchedule.back() == IntVect::Unit);

  // m_refSchedule is well defined. Set MGOptions to describe the schedule.
  m_opt.maxDepth = m_refSchedule.size() - 1;

  m_opPtrs = a_opPtrs;

  // The bottom solver.
  CH_assert(m_opPtrs[m_opt.maxDepth] != nullptr);
  m_bottomSolverPtr = new BottomSolver;
  m_bottomSolverPtr->define(m_opPtrs[m_opt.maxDepth]);
  m_bottomSolverPtr->options().numSmoothPrecond = m_opt.numSmoothPrecond;
  m_bottomSolverPtr->options().verbosity = 0;
  m_bottomSolverPtr->options().normType = m_opt.normType;
  m_isDefined = true;
}

// -----------------------------------------------------------------------------
// Sinmple check.
// -----------------------------------------------------------------------------
template <typename Hilbert, typename LinearOperator, typename BottomSolver,
          typename BC>
bool MGSolver<Hilbert, LinearOperator, BottomSolver, BC>::isDefined()
    const {
  return m_isDefined;
}

// -----------------------------------------------------------------------------
// Free's Hilbert. Leaves object unusable.
// -----------------------------------------------------------------------------
template <typename Hilbert, typename LinearOperator, typename BottomSolver,
          typename BC>
void MGSolver<Hilbert, LinearOperator, BottomSolver, BC>::clear() {
  // Don't reset the options. Doing so may cause frustration. Consider this:
  //   MGSovler mySolver;
  //   mySolver.verbosity = 10
  //   mySolver.define(&myOp);
  //   mySolver.fmg(...);
  // If define() calls clear() and clear() erases the options, the user
  // will be wondering why the solver isn't as verbose as requested.

  delete m_bottomSolverPtr;
  m_bottomSolverPtr = nullptr;

  m_opPtrs.resize(0);

  m_refSchedule.resize(0);
}

// -----------------------------------------------------------------------------
// Solves L[phi] = rhs using the most appropriate methods.
// -----------------------------------------------------------------------------
template <typename Hilbert, typename LinearOperator, typename BottomSolver,
          typename BC>
auto MGSolver<Hilbert, LinearOperator, BottomSolver, BC>::solve(
    Hilbert &a_phi, const Hilbert &a_rhs, const Real a_time,
    const bool a_useHomogBCs, const bool a_setPhiToZero) -> std::vector<Real> {
  return this->vCycle(a_phi, a_rhs, a_time, a_useHomogBCs, a_setPhiToZero);
  // return this->fmg(a_phi, a_rhs, a_time, a_useHomogBCs, a_setPhiToZero);
}

// -----------------------------------------------------------------------------
// Solves L[phi] = rhs using V-cycles.
// -----------------------------------------------------------------------------
template <typename Hilbert, typename LinearOperator, typename BottomSolver,
          typename BC>
std::vector<Real> MGSolver<Hilbert, LinearOperator, BottomSolver, BC>::vCycle(
    Hilbert &a_phi, const Hilbert &a_rhs, const Real a_time,
    const bool a_useHomogBCs, const bool a_setPhiToZero) {
  // Sanity checks
  CH_assert(a_phi.nComp() == a_rhs.nComp());  // This may be unnecessary.

  std::vector<Real> absResNorms(0);
  std::vector<Real> relResNorms(0);

  // Create workspace.
  auto op = m_opPtrs[0].get();
  Hilbert res;
  res.create(a_rhs);
  Hilbert cor;
  cor.create(a_phi);

  // Initialize to zero, if requested.
  if (a_setPhiToZero) a_phi.setToZero();

  // Set up residual equation at depth 0.

  op->residual(res, a_phi, a_rhs, a_time, a_useHomogBCs);

  // Compute initial diagnostics.
  absResNorms.push_back(res.norm());
  relResNorms.push_back(1.0);
  if (m_opt.verbosity > 3) {
    pout() << "Initial absolute |res| = " << absResNorms.back() << endl;
  }

  //Is initial guess good enough?
  if (absResNorms.back() < m_opt.absTol)
  {
      pout() << "Initial guess solves problem. Exiting solver." << endl;
      relResNorms.back();
  }

  int iter;
  for (iter = 1; iter <= m_opt.maxIters; ++iter) {
    // Set initial guess. Downward relaxation will happen in vCycle.

    op->preCond(cor, res, a_time, 0);

    // Solve residual eq. L[cor] = res.
    this->indent(0);

    this->vCycle_noInit(cor, res, a_time, 0);
    this->unindent(0);

    // Add correction to phi at depth 0.
    a_phi.incr(cor, 1.0);

    // Compute new residual at depth 0.
    op->residual(res, a_phi, a_rhs, a_time, a_useHomogBCs);

    // Compute diagnostics.
    absResNorms.push_back(res.norm());
    relResNorms.push_back(absResNorms.back() / absResNorms[0]);
    if (m_opt.verbosity >= 3) {
      pout() << "Relative |res| at iter " << iter << " = " << relResNorms.back()
             << endl;
    }

    // Did we converge?
    if (absResNorms.back() < m_opt.absTol) {
      pout() << "Absolute convergence achieved." << endl;
      break;
    }

    if (relResNorms.back() < m_opt.relTol) {
      pout() << "Relative convergence achieved." << endl;
     break;
    }

    // Are we diverging? If so, undo last correction and exit.
    if (relResNorms[iter] > relResNorms[iter - 1]) {
      pout() << "Diverging." << endl;
      a_phi.incr(cor, -1.0);
      absResNorms.pop_back();
      relResNorms.pop_back();
      --iter;
      break;
    }
  }

  // If m_opt.maxIters was reached, we need to reset iter.
  iter = min(iter, m_opt.maxIters);

  if (m_opt.verbosity >= 5) {
    pout() << "\nAbsolute convergence pattern:" << Format::indent() << endl;
    pout() << "Iter " << 0 << ": " << absResNorms[0] << "\n";
    for (int i = 1; i <= iter; ++i) {
      pout() << "Iter " << i << ": " << absResNorms[i]
             << "\t ratio = " << absResNorms[i] / absResNorms[i - 1] << "\n";
    }
    pout() << Format::unindent;

    pout() << "\nRelative convergence pattern:" << Format::indent() << endl;
    pout() << "Iter " << 0 << ": " << relResNorms[0] << "\n";
    for (int i = 1; i <= iter; ++i) {
      pout() << "Iter " << i << ": " << relResNorms[i]
             << "\t ratio = " << relResNorms[i] / relResNorms[i - 1] << "\n";
    }
    pout() << Format::unindent << endl;
  }
  return relResNorms;
}

// -----------------------------------------------------------------------------
// Solves L[phi] = rhs using Full Multigrid.
// -----------------------------------------------------------------------------
template <typename Hilbert, typename LinearOperator, typename BottomSolver,
          typename BC>
std::vector<Real> MGSolver<Hilbert, LinearOperator, BottomSolver, BC>::fmg(
    Hilbert &a_phi, const Hilbert &a_rhs, const Real a_time,
    const bool a_useHomogBCs, const bool a_setPhiToZero) {
  // Sanity checks

  std::vector<Real> absResNorms(0);
  std::vector<Real> relResNorms(0);

  // Create workspace.
  auto op = m_opPtrs[0].get();
  Hilbert res;
  res.create(a_rhs);
  Hilbert cor;
  cor.create(a_phi);

  // Initialize to zero, if requested.
  if (a_setPhiToZero) a_phi.setToZero();

  // Set up residual equation at depth 0.
  op->residual(res, a_phi, a_rhs, a_time, a_useHomogBCs);

  // Compute initial diagnostics.
  absResNorms.push_back(res.norm());
  relResNorms.push_back(1.0);
  if (m_opt.verbosity > 3) {
    pout() << "Initial absolute |res| = " << absResNorms.back() << endl;
  }

  // Is initial guess good enough?
  if (absResNorms.back() < m_opt.absTol) {
    pout() << "Initial guess solves problem. Exiting solver." << endl;
    return relResNorms;
  }

  int iter;
  for (iter = 1; iter <= m_opt.maxIters; ++iter) {
    // Solve residual eq. L[cor] = res.
    this->indent(0);
    this->fmg_noInit(cor, res, a_time, 0);
    this->unindent(0);

    // Add correction to phi at depth 0.
    a_phi.incr(cor, 1.0);

    // Compute new residual at depth 0.
    op->residual(res, a_phi, a_rhs, a_time, a_useHomogBCs);

    // Compute diagnostics.
    absResNorms.push_back(res.norm());
    relResNorms.push_back(absResNorms.back() / absResNorms[0]);
    if (m_opt.verbosity >= 3) {
      pout() << "Relative |res| at iter " << iter << " = " << relResNorms.back()
             << endl;
    }

    // Did we converge?
    if (absResNorms.back() < m_opt.absTol) {
      pout() << "Absolute convergence achieved." << endl;
      break;
    }

    if (relResNorms.back() < m_opt.relTol) {
      pout() << "Relative convergence achieved." << endl;
      break;
    }

    // Are we diverging? If so, undo last correction and exit.
    if (relResNorms[iter] > relResNorms[iter - 1]) {
      pout() << "Diverging." << endl;
      a_phi.incr(cor, -1.0);
      absResNorms.pop_back();
      relResNorms.pop_back();
      --iter;
      break;
    }
  }

  // If m_opt.maxIters was reached, we need to reset iter.
  iter = min(iter, m_opt.maxIters);

  if (m_opt.verbosity >= 5) {
    pout() << "\nAbsolute convergence pattern:" << Format::indent() << endl;
    pout() << "Iter " << 0 << ": " << absResNorms[0] << "\n";
    for (int i = 1; i <= iter; ++i) {
      pout() << "Iter " << i << ": " << absResNorms[i]
             << "\t ratio = " << absResNorms[i] / absResNorms[i - 1] << "\n";
    }
    pout() << Format::unindent;

    pout() << "\nRelative convergence pattern:" << Format::indent() << endl;
    pout() << "Iter " << 0 << ": " << relResNorms[0] << "\n";
    for (int i = 1; i <= iter; ++i) {
      pout() << "Iter " << i << ": " << relResNorms[i]
             << "\t ratio = " << relResNorms[i] / relResNorms[i - 1] << "\n";
    }
    pout() << Format::unindent << endl;
  }
  return relResNorms;
}

// -----------------------------------------------------------------------------
// pout formatting utility.
// -----------------------------------------------------------------------------
template <typename Hilbert, typename LinearOperator, typename BottomSolver,
          typename BC>
void MGSolver<Hilbert, LinearOperator, BottomSolver, BC>::indent(
    const int a_depth) const {
  if (m_opt.verbosity >= 5) {
    if (a_depth < 0) {
      pout() << Format::indent();
    } else {
      std::string bullet = std::to_string(a_depth) + std::string(": ");
      pout() << Format::indent(4, bullet.c_str());
    }
  }
}

// -----------------------------------------------------------------------------
// pout formatting utility.
// -----------------------------------------------------------------------------
template <typename Hilbert, typename LinearOperator, typename BottomSolver,
          typename BC>
void MGSolver<Hilbert, LinearOperator, BottomSolver, BC>::unindent(
    const int a_depth) const {
  if (m_opt.verbosity >= 5) {
    pout() << Format::unindent;
  }
}

// -----------------------------------------------------------------------------
// Solves L[phi] = Jrhs with homogBCs at a specified depth.
// -----------------------------------------------------------------------------
template <typename Hilbert, typename LinearOperator, typename BottomSolver,
          typename BC>
void MGSolver<Hilbert, LinearOperator, BottomSolver, BC>::vCycle_noInit(
    Hilbert &a_cor, const Hilbert &a_res, const Real a_time,
    const int a_depth) {
  if (m_opt.verbosity >= 5) {
    pout() << "MG depth = " << a_depth << endl;
  }

  // Create workspace.
  auto op = m_opPtrs[a_depth].get();
  // auto lg = m_levGeoPtrs[a_depth];
  Hilbert tmpRes;
  tmpRes.create(a_res);

  const Real downOmega = 1.0;
  const Real upOmega = 1.0;
  const Real bottomOmega = 1.0;
  const int prolongOrder = 0;

  // Compute rhs integral.
  if (m_opt.verbosity >= 8) {
    Real norm = a_res.norm();
    // Integral::sum is not defined for device_memory objects.
    // Real sumRes = Integral::sum(a_res, lg->getDXi(), 0);
    // Real sumCor = Integral::sum(a_cor, lg->getDXi(), 0);
    pout() << "|rhs| = "
           << norm
           //<< ", sum res = " << sumRes
           //<< ", sum cor = " << sumCor
           << endl;
  }

  if (a_depth == m_opt.maxDepth) {
    // Use bottom solver..

    // --- Bottom relaxation ---
    if (m_opt.verbosity >= 7) {
      pout() << "Bottom relax" << endl;
    }

    op->relax(a_cor, a_res, a_time, m_opt.numSmoothBottom, bottomOmega);

    // Diagnostics
    if (m_opt.verbosity >= 8) {
      op->residual(tmpRes, a_cor, a_res, a_time, true);
      Real norm = tmpRes.norm();
      pout() << "|rhs| = " << norm << endl;
    }

    // --- Bottom solver ---
    if (m_opt.verbosity >= 7) {
      pout() << "Bottom solver" << endl;
    }
    this->indent();
    // TODO: Do something with exitStatus?
    m_bottomSolverPtr->solve(a_cor, a_res, a_time, true, false);
    this->unindent();

    // Diagnostics
    if (m_opt.verbosity >= 8) {
      op->residual(tmpRes, a_cor, a_res, a_time, true);
      Real norm = tmpRes.norm();
      pout() << "|rhs| = " << norm << endl;
    }
  } else {
    // V-Cycle...

    // Create needed structures.

    Hilbert crseCor, crseRes;

    op->createCoarsened(crseCor, a_cor, crseRes, a_res);
    // --- Downward relaxation ---
    if (m_opt.verbosity >= 7) {
      pout() << "Smooth down" << endl;
    }

    op->relax(a_cor, a_res, a_time, m_opt.numSmoothDown, downOmega);

    op->applyOp(tmpRes, a_cor, 0.0, true);

    // Compute residual
    op->residual(tmpRes, a_cor, a_res, a_time, true);

    // Diagnostics
    if (m_opt.verbosity >= 8) {
      Real norm = tmpRes.norm();
      // Real sumRes = Integral::sum(tmpRes, lg->getDXi(), 0);
      pout() << "|rhs| = "
             << norm
             //<< " sum tmpRes = " << sumRes
             << endl;
    }

    // --- Restriction ---
    if (m_opt.verbosity >= 7) {
      pout() << "Restrict residual" << endl;
    }
    op->restrictResidual(crseRes, tmpRes, a_time, m_refSchedule[a_depth]);

    // --- Coarse level solve ---

    // Set initial guess using operator at coarser level
    m_opPtrs[a_depth + 1]->preCond(crseCor, crseRes, a_time, 0);

    // Solve
    this->indent(a_depth + 1);
    for (int i = 0; i < m_opt.numCycles; ++i) {
      this->vCycle_noInit(crseCor, crseRes, a_time, a_depth + 1);
    }
    this->unindent(a_depth + 1);

    // --- Prolong correction ---
    if (m_opt.verbosity >= 7) {
      pout() << "Prolong and add correction" << endl;
    }

    op->prolongIncrement(a_cor, crseCor, a_time, m_refSchedule[a_depth],
                         prolongOrder);

    // Diagnostics

    if (m_opt.verbosity >= 8) {
      op->residual(tmpRes, a_cor, a_res, a_time, true);
      Real norm = tmpRes.norm();
      pout() << "|rhs| = " << norm << endl;
    }

    // --- Upward relaxation ---
    if (m_opt.verbosity >= 7) {
      pout() << "Smooth up" << endl;
    }

    op->relax(a_cor, a_res, a_time, m_opt.numSmoothUp, upOmega);

    // Diagnostics
    if (m_opt.verbosity >= 8) {
      op->residual(tmpRes, a_cor, a_res, a_time, true);
      Real norm = tmpRes.norm();
      // Real sum = Integral::sum(a_cor, lg->getDXi(), 0);
      pout() << "|rhs| = "
             << norm
             //<< ", sum cor = " << sum
             << endl;
    }
  }
}

// -----------------------------------------------------------------------------
// Solves L[phi] = Jrhs with homogBCs at a specified depth.
// -----------------------------------------------------------------------------
template <typename Hilbert, typename LinearOperator, typename BottomSolver,
          typename BC>
void MGSolver<Hilbert, LinearOperator, BottomSolver, BC>::fmg_noInit(
    Hilbert &a_cor, const Hilbert &a_res, const Real a_time,
    const int a_depth) {
  auto op = m_opPtrs[a_depth].get();
  const int upOmega = 1.0;

  // Initialize solution at this depth.
  a_cor.setToZero();

  // Go to the coarser MG level if possible.
  if (a_depth < m_opt.maxDepth) {
    // Create needed structures.
    auto crseOp = m_opPtrs[a_depth + 1].get();

    const auto &crseRef = m_refSchedule[a_depth];
    Hilbert crseCor, crseRes;
    op->createCoarsened(crseCor, a_cor, crseRes, a_res);

    // --- Restrict residual ---
    if (m_opt.verbosity >= 7) {
      pout() << "Restrict residual" << endl;
    }
    op->restrictResidual(crseRes, a_res, a_time, crseRef);

    // --- Solve ---
    this->indent(a_depth + 1);
    this->fmg_noInit(crseCor, crseRes, a_time, a_depth + 1);
    this->unindent(a_depth + 1);

    // --- Prolong correction ---
    if (m_opt.verbosity >= 7) {
      pout() << "Prolong and add correction" << endl;
    }
    op->prolongIncrement(a_cor, crseCor, a_time, crseRef,
                         m_opt.prolongOrderFMG);

    // --- Additional relaxation for FMG prolong ---
    if (m_opt.verbosity >= 7) {
      pout() << "Smooth up for FMG prolong" << endl;
    }
    op->relax(a_cor, a_res, a_time, m_opt.numSmoothUpFMG, upOmega);

    // Diagnostics
    if (m_opt.verbosity >= 8) {
      Hilbert tmpRes;
      tmpRes.create(a_res);
      op->residual(tmpRes, a_cor, a_res, a_time, true);
      Real norm = tmpRes.norm();
      pout() << "|rhs| = " << norm << endl;
    }
  }

  // Solve with a V-cycle.
  for (int i = 0; i < m_opt.numCycles; ++i) {
    this->vCycle_noInit(a_cor, a_res, a_time, a_depth);
  }
}
