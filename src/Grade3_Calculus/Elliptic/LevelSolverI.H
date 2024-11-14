#ifndef H9c9778ae61b21da3b93db6cca54d9de6
#error LevelSolverI.H should only be used by LevelSolver.H
#endif

#include "LevelSolver.H" // To help VSCode find symbols.
#include "SOMAR_Constants.H"
#include "Format.H"

namespace Elliptic {

// ============================== LevelSolver ==================================

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template <class StateType>
LevelSolver<StateType>::LevelSolver()
{
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template <class StateType>
LevelSolver<StateType>::~LevelSolver()
{
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template <class StateType>
void
LevelSolver<StateType>::define(
    std::shared_ptr<const LevelOperator<StateType>> a_opPtr)
{
    CH_assert(a_opPtr);
    m_opPtr = a_opPtr;
    m_solverStatus.clear();
}


// ============================== RelaxSolver ==================================

// -----------------------------------------------------------------------------
// Constructor. You must still call define() to supply an operator.
// -----------------------------------------------------------------------------
template <class StateType>
RelaxSolver<StateType>::RelaxSolver(Options a_opt)
: m_opt(a_opt)
{
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
template <class StateType>
RelaxSolver<StateType>::~RelaxSolver()
{
}


// -----------------------------------------------------------------------------
// Solve!
// If a_convergenceMetric is provided, it will take precedence over
// Options::convergenceMetric.
// -----------------------------------------------------------------------------
template <class StateType>
SolverStatus
RelaxSolver<StateType>::solve(StateType&       a_phi,
                              const StateType* a_crsePhiPtr,
                              const StateType& a_rhs,
                              const Real       a_time,
                              const bool       a_useHomogBCs,
                              const bool       a_setPhiToZero,
                              const Real       a_convergenceMetric) const
{
    SolverStatus& solverStatus = this->getSolverStatusRef();
    const LevelOperator<StateType>& op = this->getOp();

    pout() << Format::pushFlags;
    solverStatus.clear();

    if (a_setPhiToZero) {
        op.setToZero(a_phi);
    }

    StateType r, e;

    op.create(r, a_rhs);
    op.create(e, a_phi);

    // Set up residual equation and compute initial norms.
    op.residual(
        r, a_phi, a_crsePhiPtr, a_rhs, a_time, a_useHomogBCs, a_useHomogBCs);

    Vector<Real> absResNorms(1);
    absResNorms.reserve(m_opt.maxIters + 1);
    absResNorms[0] = op.norm(r, m_opt.normType);
    this->getSolverStatus().setInitResNorm(absResNorms[0]);

    Vector<Real> relResNorms(1);
    relResNorms.reserve(m_opt.maxIters + 1);
    relResNorms[0] = 1.0;

    if (m_opt.verbosity >= 3) {
        pout() << "RelaxSolver: initial residual norm = "
               << Format::number(absResNorms[0]) << "\n";
    }

    // if a convergence metric has been supplied, replace initial residual
    // with the supplied convergence metric...
    if (m_opt.convergenceMetric > 0.0) {
        absResNorms[0] = m_opt.convergenceMetric;
    }
    if (a_convergenceMetric > 0.0) {
        absResNorms[0] = a_convergenceMetric;
    }

    // Are we already at absTol?
    if (absResNorms[0] < m_opt.absTol && m_opt.minIters <= 0) {
        if (m_opt.verbosity >= 3) {
            pout() << "RelaxSolver: Initial guess solves problem. Exiting.\n";
        }
        pout() << Format::popFlags;

        solverStatus.setFinalResNorm(absResNorms.back());
        solverStatus.setSolverStatus(SolverStatus::CONVERGED);
        return solverStatus;
    }

    // Main loop
    int i = 1;
    do {
        // Relax
        op.setToZero(e);
        op.relax(e, r, a_time, m_opt.itersPerNormCheck);

        // Recompute residual and norms.
        op.incr(a_phi, e, 1.0);
        op.residual(r,
                    a_phi,
                    a_crsePhiPtr,
                    a_rhs,
                    a_time,
                    a_useHomogBCs,
                    a_useHomogBCs);
        const Real thisResNorm = op.norm(r, m_opt.normType);
        absResNorms.push_back(thisResNorm);
        relResNorms.push_back(thisResNorm / absResNorms[0]);

        // Report residual norm, if requested.
        if (m_opt.verbosity >= 4) {
            pout() << "RelaxSolver: iter " << i * m_opt.itersPerNormCheck
                   << " rel res norm = " << Format::scientific
                   << relResNorms.back() << '\n';
        }

        if (i * m_opt.itersPerNormCheck >= m_opt.minIters) {
            // Did we reach absTol?
            if (absResNorms[i] < m_opt.absTol) {
                if (m_opt.verbosity >= 3) {
                    pout() << "RelaxSolver: Converged to absTol.\n";
                }
                solverStatus.setSolverStatus(SolverStatus::CONVERGED);
                break;
            }

            // Did we reach relTol?
            if (relResNorms[i] < m_opt.relTol) {
                if (m_opt.verbosity >= 3) {
                    pout() << "RelaxSolver: Converged to relTol.\n";
                }
                solverStatus.setSolverStatus(SolverStatus::CONVERGED);
                break;
            }

            // Did we hang?
            if (abs(absResNorms[i] - absResNorms[i - 1]) < m_opt.hang) {
                if (m_opt.verbosity >= 3) {
                    pout() << "RelaxSolver: Hanging.\n";
                }
                solverStatus.setSolverStatus(SolverStatus::HANG);
                break;
            }

            // Are we diverging? If so, remove the last correction.
            if (absResNorms[i] > absResNorms[i - 1]) {
                op.incr(a_phi, e, -1.0);

                if (m_opt.verbosity >= 3) {
                    pout() << "RelaxSolver: Diverging.\n";
                }
                solverStatus.setSolverStatus(SolverStatus::DIVERGED);
                break;
            }
        }

        // Did we reach maxIters?
        if (i * m_opt.itersPerNormCheck >= m_opt.maxIters) {
            if (m_opt.verbosity >= 3) {
                pout() << "RelaxSolver: maxIters reached.\n";
            }
            solverStatus.setSolverStatus(SolverStatus::MAXITERS);
            break;
        }

        ++i;
    } while (true);


    if (1 <= m_opt.verbosity && m_opt.verbosity <= 3) {
        pout() << "RelaxSolver: " << i * m_opt.itersPerNormCheck
               << " iters. Final relative res norm = " << Format::scientific
               << relResNorms.back() << '\n';
    }

    pout() << Format::popFlags;
    solverStatus.setFinalResNorm(absResNorms.back());
    return solverStatus;
}


// ============================ BiCGStabSolver =================================

// -----------------------------------------------------------------------------
// Constructor. You must still call define() to supply an operator.
// -----------------------------------------------------------------------------
template <class StateType>
BiCGStabSolver<StateType>::BiCGStabSolver(Options a_opt)
: m_opt(a_opt)
{
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
template <class StateType>
BiCGStabSolver<StateType>::~BiCGStabSolver()
{
}


// -----------------------------------------------------------------------------
// Solve!
// If a_convergenceMetric is provided, it will take precedence over
// Options::convergenceMetric.
// -----------------------------------------------------------------------------
template <class StateType>
SolverStatus
BiCGStabSolver<StateType>::solve(StateType&       a_phi,
                                 const StateType* a_crsePhiPtr,
                                 const StateType& a_rhs,
                                 const Real       a_time,
                                 const bool       a_useHomogBCs,
                                 const bool       a_setPhiToZero,
                                 const Real       a_convergenceMetric) const
{
    SolverStatus& solverStatus = this->getSolverStatusRef();

    pout() << Format::pushFlags;
    solverStatus.setSolverStatus(SolverStatus::UNDEFINED);

    const LevelOperator<StateType>& op = this->getOp();

    if (a_setPhiToZero) {
        op.setToZero(a_phi);
    }

    StateType r, r_tilde, e, p, p_tilde, s_tilde, t, v;

    op.create(r, a_rhs);
    op.create(r_tilde, a_rhs);
    op.create(e, a_phi);
    op.create(p, a_rhs);
    op.create(p_tilde, a_phi);
    op.create(s_tilde, a_phi);
    op.create(t, a_rhs);
    op.create(v, a_rhs);

    int recount = 0;

    op.residual(
        r, a_phi, a_crsePhiPtr, a_rhs, a_time, a_useHomogBCs, a_useHomogBCs);

    op.assignLocal(r_tilde, r);
    op.setToZero(e);
    // (DFM 2/1/07) these next two need to be set to zero to prevent
    // problems in the multilevel case
    op.setToZero(p_tilde);
    op.setToZero(s_tilde);

    int i = 0;

    // rho[0] = r_i , rho[1] = r_(i-1), etc.
    Real rho[4] = { 0, 0, 0, 0 };
    Real norm[2];
    norm[0]            = op.norm(r, m_opt.normType);
    Real initial_norm  = norm[0];
    Real initial_rnorm = norm[0];
    norm[1]            = norm[0];
    solverStatus.setInitResNorm(initial_norm);

    Real alpha[2] = { 0, 0 };
    Real beta[2]  = { 0, 0 };
    Real omega[2] = { 0, 0 };

    bool init     = true;
    int  restarts = 0;

    // if a convergence metric has been supplied, replace initial residual
    // with the supplied convergence metric...
    if (m_opt.convergenceMetric > 0.0) {
        initial_norm = m_opt.convergenceMetric;
    }
    if (a_convergenceMetric > 0.0) {
        initial_norm = a_convergenceMetric;
    }

    if (m_opt.verbosity >= 3) {
        pout() << "BiCGStab: initial Residual norm = " << Format::number(initial_norm) << "\n";
    }


    while ((i < m_opt.maxIters && norm[0] > m_opt.absTol * norm[1]) &&
           (norm[1] > 0)) {
        i++;

        norm[1]  = norm[0];
        alpha[1] = alpha[0];
        beta[1]  = beta[0];
        omega[1] = omega[0];

        if (m_opt.verbosity >= 5) {
            pout() << "BiCGStab:       norm[0]  = " << Format::number(norm[0]) << ", "
                   << "norm[1]  = " << norm[1] << "\n";
            pout() << "BiCGStab:       alpha[0] = " << Format::number(alpha[0]) << ", "
                   << "alpha[1] = " << alpha[1] << "\n";
            pout() << "BiCGStab:       beta[0]  = " << Format::number(beta[0]) << ", "
                   << "beta[1]  = " << beta[1] << "\n";
            pout() << "BiCGStab:       omega[0] = " << Format::number(omega[0]) << ", "
                   << "omega[1] = " << omega[1] << "\n";
        }

        rho[3] = rho[2];
        rho[2] = rho[1];
        rho[1] = op.dotProduct(r_tilde, r);

        if (m_opt.verbosity >= 5) {
            pout() << "BiCGStab:       "
                   << "rho[1] = " << Format::number(rho[1]) << ", "
                   << "rho[2] = " << Format::number(rho[2]) << ", "
                   << "rho[3] = " << Format::number(rho[3]) << "\n";
        }

        if (RealCmp::isZero(rho[1])) {
            // we are finished, we will not converge anymore
            op.incr(a_phi, e, 1.0);

            if (m_opt.verbosity >= 3) {
                pout() << "BiCGStab: rho = " << Format::number(rho[1]) << ", returning"
                       << " -- Residual norm = " << Format::number(norm[0])
                       << "\n";
            }

            // need to call clear just in case
            // this is not ideal -- maybe change the return
            // to exit
            op.clear(r);
            op.clear(r_tilde);
            op.clear(e);
            op.clear(p);
            op.clear(p_tilde);
            op.clear(s_tilde);
            op.clear(t);
            op.clear(v);

            pout() << Format::popFlags;
            solverStatus.setFinalResNorm(initial_norm);
            solverStatus.setSolverStatus(SolverStatus::SINGULAR);
            return solverStatus;
        }

        if (init) {
            op.assignLocal(p, r);
            init = false;
        } else {
            beta[1] = (rho[1] / rho[2]) * (alpha[1] / omega[1]);
            op.scale(p, beta[1]);
            op.incr(p, v, -beta[1] * omega[1]);
            op.incr(p, r, 1.0);
        }

        if (m_opt.verbosity >= 5) {
            pout() << "BiCGStab:       beta[1]  = " << Format::number(beta[1])
                   << "\n";
        }

        op.preCond(p_tilde, p, a_time, m_opt.numSmoothPrecond);
        op.applyOp(v, p_tilde, nullptr, a_time, true, true);
        Real m   = op.dotProduct(r_tilde, v);
        alpha[0] = rho[1] / m;

        if (m_opt.verbosity >= 5) {
            pout() << "BiCGStab:       "
                   << "rho[1] = " << Format::number(rho[1]) << ", "
                   << "m = " << Format::number(m) << ", "
                   << "alpha[0] = " << Format::number(alpha[0]) << "\n";
        }

        if (Abs(m) > m_opt.small * Abs(rho[1])) {
            op.incr(r, v, -alpha[0]);
            norm[0] = op.norm(r, m_opt.normType);
            op.incr(e, p_tilde, alpha[0]);
        } else {
            op.setToZero(r);
            norm[0] = 0.0;
        }

        if (norm[0] > m_opt.absTol * initial_norm &&
            norm[0] > m_opt.relTol * initial_rnorm) {
            op.preCond(s_tilde, r, a_time, m_opt.numSmoothPrecond);
            op.applyOp(t, s_tilde, nullptr, a_time, true, true);
            omega[0] = op.dotProduct(t, r) / op.dotProduct(t, t);
            op.incr(e, s_tilde, omega[0]);
            op.incr(r, t, -omega[0]);
            norm[0] = op.norm(r, m_opt.normType);
        }

        if (m_opt.verbosity >= 4) {
            pout() << "BiCGStab:     iteration = " << Format::fixed << i
                   << ", error norm = " << Format::number(norm[0])
                   << ", rate = " << Format::number(norm[1] / norm[0]) << "\n";
        }

        if (norm[0] <= m_opt.absTol * initial_norm ||
            norm[0] <= m_opt.relTol * initial_rnorm) {
            // converged to tolerance
            solverStatus.setFinalResNorm(norm[0]);
            solverStatus.setSolverStatus(SolverStatus::CONVERGED);
            break;
        }

        if (omega[0] == 0.0 || norm[0] > (1.0 - m_opt.hang) * norm[1]) {
            if (recount == 0) {
                recount = 1;
            } else {
                recount = 0;
                op.incr(a_phi, e, 1.0);

                if (restarts == m_opt.maxRestarts) {
                    if (m_opt.verbosity >= 3) {
                        pout() << "BiCGStab: max restarts reached"
                               << "\ninit  norm = " << Format::number(initial_norm)
                               << "\nfinal norm = " << Format::number(norm[0])
                               << '\n';
                    }


                    // need to call clear just in case
                    // this is not ideal -- maybe change the return
                    // to exit
                    op.clear(r);
                    op.clear(r_tilde);
                    op.clear(e);
                    op.clear(p);
                    op.clear(p_tilde);
                    op.clear(s_tilde);
                    op.clear(t);
                    op.clear(v);

                    pout() << Format::popFlags;
                    solverStatus.setFinalResNorm(norm[0]);
                    solverStatus.setSolverStatus(SolverStatus::MAXITERS);
                    return solverStatus;
                }

                {
                    op.residual(r,
                                a_phi,
                                a_crsePhiPtr,
                                a_rhs,
                                a_time,
                                a_useHomogBCs,
                                a_useHomogBCs);
                    norm[0]  = op.norm(r, m_opt.normType);
                    rho[1]   = 0.0;
                    rho[1]   = 0.0;
                    rho[2]   = 0.0;
                    rho[3]   = 0.0;
                    alpha[0] = 0;
                    beta[0]  = 0;
                    omega[0] = 0;
                    op.assignLocal(r_tilde, r);
                    op.setToZero(e);

                    restarts++;
                }

                if (m_opt.verbosity >= 4) {
                    pout() << "BiCGStab:   restart =  "
                           << Format::fixed << restarts << "\n";
                }

                init = true;
            }
        }
    }

    if (m_opt.verbosity >= 3) {
        pout() << "BiCGStab: " << Format::fixed << i
               << " iterations, final Residual norm = "
               << Format::number(norm[0]) << "\n";
    }
    if (m_opt.verbosity >= 2) {
        pout() << "BiCGStab: " << Format::fixed << i
               << " iterations, relative residual = "
               << Format::number(norm[0] / initial_norm)
               << '\n';
    }

    op.incr(a_phi, e, 1.0);

    op.clear(r);
    op.clear(r_tilde);
    op.clear(e);
    op.clear(p);
    op.clear(p_tilde);
    op.clear(s_tilde);
    op.clear(t);
    op.clear(v);

    pout() << Format::popFlags;
    solverStatus.setFinalResNorm(norm[0]);
    return solverStatus;
}

}; // namespace Elliptic
