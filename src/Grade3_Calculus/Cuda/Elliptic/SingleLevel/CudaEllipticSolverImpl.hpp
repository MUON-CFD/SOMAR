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
#ifndef H9929fb0cc38edce433e3b0842ce1c397
#error CudaEllipticSolverImpl.hpp should only be used by EllipticSolver.H
#endif

// ============================ EllipticSolver =================================

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template <class T, class L>
EllipticSolver<T, L>::EllipticSolver() : m_opPtr(nullptr) {}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template <class T, class L> EllipticSolver<T, L>::~EllipticSolver() {
  m_opPtr = nullptr;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template <class T, class L>
void EllipticSolver<T, L>::define(const std::shared_ptr<L> a_opPtr) {

  CH_assert(a_opPtr.get());
  m_opPtr = a_opPtr;
}

// -----------------------------------------------------------------------------
// Convert the exit status to a readable message.
// -----------------------------------------------------------------------------
template <class T, class L>
const char *
EllipticSolver<T, L>::getExitMessage(const int a_exitStatus) const {
  switch (a_exitStatus) {
  case UNDEFINED:
    return "undefined";
  case CONVERGED:
    return "converged";
  case DIVERGED:
    return "diverged";
  case HANG:
    return "hanging";
  case MAXITERS:
    return "max iters reached";
  default:
    return "unknown";
  }
}

// ============================ BiCGStabSolver =================================

// -----------------------------------------------------------------------------
// Constructor. You must still call define() to supply an operator.
// -----------------------------------------------------------------------------
template <class T, class L>
BiCGStabSolver<T, L>::BiCGStabSolver(Options a_opt)
    : EllipticSolver<T, L>(), m_opt(a_opt) {}

// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
template <class T, class L> BiCGStabSolver<T, L>::~BiCGStabSolver() {}

// -----------------------------------------------------------------------------
// Solve!
// -----------------------------------------------------------------------------
template <class T, class L>
int BiCGStabSolver<T, L>::solve(T &a_phi, const T &a_rhs, const Real a_time,
                                    const bool a_useHomogBCs,
                                    const bool a_setPhiToZero) const {
  this->setExitStatus(EllipticSolver<T, L>::UNDEFINED);

  const L &op = this->getOp();

  if (a_setPhiToZero) {
    a_phi.setToZero();
  }

  T r, r_tilde, e, p, p_tilde, s_tilde, t, v;

  r.create(a_rhs);
  r_tilde.create(a_rhs);
  e.create(a_phi);
  p.create(a_rhs);
  p_tilde.create(a_phi);
  s_tilde.create(a_phi);
  t.create(a_rhs);
  v.create(a_rhs);

  int recount = 0;

  op.residual(r, a_phi, a_rhs, a_time, a_useHomogBCs);

  r_tilde.assign(r);
  e.setToZero();
  // (DFM 2/1/07) these next two need to be set to zero to prevent
  // problems in the multilevel case
  p_tilde.setToZero();
  s_tilde.setToZero();

  int i = 0;

  // rho[0] = r_i , rho[1] = r_(i-1), etc.
  Real rho[4] = {0, 0, 0, 0};
  Real norm[2];
  norm[0] = r.norm();
  Real initial_norm = norm[0];
  Real initial_rnorm = norm[0];
  norm[1] = norm[0];

  Real alpha[2] = {0, 0};
  Real beta[2] = {0, 0};
  Real omega[2] = {0, 0};

  bool init = true;
  int restarts = 0;

  if (m_opt.verbosity >= 3) {
    pout() << "BiCGStab:: initial Residual norm = " << initial_norm << "\n";
  }

  // if a convergence metric has been supplied, replace initial residual
  // with the supplied convergence metric...
  if (m_opt.convergenceMetric > 0.0) {
    initial_norm = m_opt.convergenceMetric;
  }

  while ((i < m_opt.maxIters && norm[0] > m_opt.absTol * norm[1]) &&
         (norm[1] > 0)) {
    i++;

    norm[1] = norm[0];
    alpha[1] = alpha[0];
    beta[1] = beta[0];
    omega[1] = omega[0];

    if (m_opt.verbosity >= 5) {
      pout() << "BiCGStab::       norm[0]  = " << norm[0] << ", "
             << "norm[1]  = " << norm[1] << "\n";
      pout() << "BiCGStab::       alpha[0] = " << alpha[0] << ", "
             << "alpha[1] = " << alpha[1] << "\n";
      pout() << "BiCGStab::       beta[0]  = " << beta[0] << ", "
             << "beta[1]  = " << beta[1] << "\n";
      pout() << "BiCGStab::       omega[0] = " << omega[0] << ", "
             << "omega[1] = " << omega[1] << "\n";
    }

    rho[3] = rho[2];
    rho[2] = rho[1];
    rho[1] = r_tilde.dotProduct(r);

    if (m_opt.verbosity >= 5) {
      pout() << "BiCGStab::       rho[1] = " << rho[1] << ", "
             << "rho[2] = " << rho[2] << ", "
             << "rho[3] = " << rho[3] << "\n";
    }

    if (rho[1] == 0.0) {
      // we are finished, we will not converge anymore
      a_phi.incr(e, 1.0);

      if (m_opt.verbosity >= 3) {
        pout() << "BiCGStab:: rho = 0, returning"
               << " -- Residual norm = " << norm[0] << "\n";
      }

      // need to call clear just in case
      // this is not ideal -- maybe change the return
      // to exit
      r.clear();
      r_tilde.clear();
      e.clear();
      p.clear();
      p_tilde.clear();
      s_tilde.clear();
      t.clear();
      v.clear();

      this->setExitStatus(EllipticSolver<T, L>::SINGULAR);
      return this->getExitStatus();
    }

    if (init) {
      p.assign(r);
      init = false;
    } else {
      beta[1] = (rho[1] / rho[2]) * (alpha[1] / omega[1]);
      p.scale(beta[1]);
      p.incr(v, -beta[1] * omega[1]);
      p.incr(r, 1.0);
    }

    if (m_opt.verbosity >= 5) {
      pout() << "BiCGStab::       beta[1]  = " << beta[1] << "\n";
    }

    op.preCond(p_tilde, p, a_time, m_opt.numSmoothPrecond);
    op.applyOp(v, p_tilde, a_time, true);
    Real m = r_tilde.dotProduct(v);
    alpha[0] = rho[1] / m;

    if (m_opt.verbosity >= 5) {
      pout() << "BiCGStab::       rho[1] = " << rho[1] << ", "
             << "m = " << m << ", "
             << "alpha[0] = " << alpha[0] << "\n";
    }

    if (Abs(m) > m_opt.small * Abs(rho[1])) {
      r.incr(v, -alpha[0]);
      norm[0] = r.norm();
      e.incr(p_tilde, alpha[0]);
    } else {
      r.setToZero();
      norm[0] = 0.0;
    }

    if (norm[0] > m_opt.absTol * initial_norm &&
        norm[0] > m_opt.relTol * initial_rnorm) {
      op.preCond(s_tilde, r, a_time, m_opt.numSmoothPrecond);

      // Util::FuncTimer("BiCGStab applyOp ", [&]() -> void {
      op.applyOp(t, s_tilde, a_time, true);
      // });
      omega[0] = t.dotProduct(r) / t.dotProduct(t);
      e.incr(s_tilde, omega[0]);
      r.incr(t, -omega[0]);
      norm[0] = r.norm();
    }

    if (m_opt.verbosity >= 4) {
      pout() << "BiCGStab::     iteration = " << i
             << ", error norm = " << norm[0] << ", rate = " << norm[1] / norm[0]
             << "\n";
    }

    if (norm[0] <= m_opt.absTol * initial_norm ||
        norm[0] <= m_opt.relTol * initial_rnorm) {
      // converged to tolerance
      this->setExitStatus(EllipticSolver<T, L>::CONVERGED);
      break;
    }

    if (omega[0] == 0.0 || norm[0] > (1.0 - m_opt.hang) * norm[1]) {
      if (recount == 0) {
        recount = 1;
      } else {
        recount = 0;
        a_phi.incr(e, 1.0);

        if (restarts == m_opt.maxRestarts) {
          if (m_opt.verbosity >= 3) {
            pout() << "BiCGStab: max restarts reached" << endl;
            pout() << "init  norm = " << initial_norm << endl;
            pout() << "final norm = " << norm[0] << endl;
          }

          // need to call clear just in case
          // this is not ideal -- maybe change the return
          // to exit
          r.clear();
          r_tilde.clear();
          e.clear();
          p.clear();
          p_tilde.clear();
          s_tilde.clear();
          t.clear();
          v.clear();

          this->setExitStatus(EllipticSolver<T, L>::MAXITERS);
          return this->getExitStatus();
        }

        {
          op.residual(r, a_phi, a_rhs, a_time, a_useHomogBCs);
          norm[0] = r.norm();
          rho[1] = 0.0;
          rho[1] = 0.0;
          rho[2] = 0.0;
          rho[3] = 0.0;
          alpha[0] = 0;
          beta[0] = 0;
          omega[0] = 0;
          r_tilde.assign(r);
          r.setToZero();

          restarts++;
        }

        if (m_opt.verbosity >= 4) {
          pout() << "BiCGStab::   restart =  " << restarts << "\n";
        }

        init = true;
      }
    }
  }

  if (m_opt.verbosity >= 3) {
    pout() << "BiCGStab:: " << i
           << " iterations, final Residual norm = " << norm[0] << "\n";
  }
  if (m_opt.verbosity >= 2) {
    pout() << "BiCGStab: " << i
           << " iterations, relative residual = " << norm[0] / initial_norm
           << endl;
  }

  a_phi.incr(e, 1.0);

  r.clear();
  r_tilde.clear();
  e.clear();
  p.clear();
  p_tilde.clear();
  s_tilde.clear();
  t.clear();
  v.clear();

  return this->getExitStatus();
}
