#ifndef H955d3c5f3254f48107821d4c7c238980
#error MGSolverI.H should only be used by MGSolver.H
#endif

#include "MGSolver.H" // This helps VSCode.
#include "MGCoarseningStrategy.H"
#include "Integral.H" // Assumes StateType = LevelData<FArrayBox>!
#include "ProblemContext.H"

namespace Elliptic {


// ======================== Option setters / getters ===========================

// -----------------------------------------------------------------------------
template <class StateType>
typename MGSolver<StateType>::Options
MGSolver<StateType>::getDefaultOptions()
{
    Options opt;
    const auto& proj = ProblemContext::getInstance()->proj;

    opt.absTol            = proj.absTol;
    opt.relTol            = proj.relTol;
    opt.convergenceMetric = -1.0;
    opt.numSmoothDown     = proj.numSmoothDown;
    opt.numSmoothUp       = proj.numSmoothUp;
    opt.numSmoothBottom   = proj.numSmoothBottom;
    opt.numSmoothPrecond  = proj.numSmoothPrecond;
    opt.prolongOrder      = proj.prolongOrder;
    opt.prolongOrderFMG   = proj.prolongOrderFMG;
    opt.numSmoothUpFMG    = proj.numSmoothUpFMG;
    opt.maxDepth          = proj.maxDepth;
    opt.numCycles         = proj.numCycles;
    opt.maxIters          = proj.maxIters;
    opt.hang              = proj.hang;
    opt.normType          = proj.normType;
    opt.verbosity         = proj.verbosity;

    opt.bottomOptions.absTol           = proj.bottom_absTol;
    opt.bottomOptions.relTol           = proj.bottom_relTol;
    opt.bottomOptions.small            = proj.bottom_small;
    opt.bottomOptions.hang             = proj.bottom_hang;
    opt.bottomOptions.maxIters         = proj.bottom_maxIters;
    opt.bottomOptions.maxRestarts      = proj.bottom_maxRestarts;
    opt.bottomOptions.normType         = proj.bottom_normType;
    opt.bottomOptions.verbosity        = proj.bottom_verbosity;
    opt.bottomOptions.numSmoothPrecond = proj.bottom_numSmoothPrecond;

    return opt;
}


// -----------------------------------------------------------------------------
template <class StateType>
typename MGSolver<StateType>::Options
MGSolver<StateType>::getQuickAndDirtyOptions()
{
    Options opt = getDefaultOptions();

    opt.absTol            = 1.0e-300;
    opt.relTol            = 1.0e-300;
    opt.numCycles         = -1;
    opt.maxIters          = 1;
    opt.verbosity         = 0;

    opt.bottomOptions.absTol           = 1.0e-300;
    opt.bottomOptions.relTol           = 1.0e-300;
    // opt.bottomOptions.maxIters         = 80;
    // opt.bottomOptions.maxRestarts      = 5;
    opt.bottomOptions.verbosity        = 0;

    return opt;
}


// -----------------------------------------------------------------------------
template <class StateType>
const typename MGSolver<StateType>::Options&
MGSolver<StateType>::getOptions() const
{
    return m_opt;
}


// -----------------------------------------------------------------------------
template <class StateType>
void
MGSolver<StateType>::modifyOptionsExceptMaxDepth(
    const MGSolver<StateType>::Options& a_opt)
{
    if (!this->isDefined()) {
        MAYDAYERROR("This can only be called AFTER MGSolver is defined.");
    }

    const auto oldMaxDepth = m_opt.maxDepth;
    m_opt = a_opt;
    m_opt.maxDepth = oldMaxDepth;

    if (m_bottomSolverPtr) {
        m_bottomSolverPtr->setOptions(a_opt.bottomOptions);
    }
}


// ======================== Constructors / destructors =========================

// -----------------------------------------------------------------------------
template <class StateType>
MGSolver<StateType>::MGSolver()
{
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
template <class StateType>
MGSolver<StateType>::~MGSolver()
{
    this->clear();
}


// -----------------------------------------------------------------------------
template <class StateType>
void
MGSolver<StateType>::define(std::shared_ptr<const MGOpType> a_topOpPtr,
                            Options                         a_opt,
                            Vector<IntVect>                 a_mgRefSchedule,
                            const bool                      a_useBottomSolver)
{
    this->define(*a_topOpPtr, a_opt, a_mgRefSchedule, a_useBottomSolver);
}


// -----------------------------------------------------------------------------
template <class StateType>
void
MGSolver<StateType>::define(const MGOpType& a_topOp,
                            Options         a_opt,
                            Vector<IntVect> a_refSchedule,
                            bool            a_useBottomSolver)
{
    m_opt = a_opt;

    // Use our own MGCoarseningStrategy? If so, use SemicoarseningStrategy.
    if (a_refSchedule.size() == 0) {
        const RealVect&  topDXi       = a_topOp.getDXi();
        const GridsType& topGrids     = a_topOp.getBoxes();
        const IntVect&   topNx        = topGrids.physDomain().size();
        const RealVect   domainLength = RealVect(topNx) * topDXi;

        const ProblemContext* ctx = ProblemContext::getInstance();
        if (ctx->proj.relaxMethod ==
            ProjectorParameters::RelaxMethod::VERTLINE) {
            constexpr bool doVertCoarsening = true;
            m_refSchedule =
                HorizCoarseningStrategy(
                    domainLength, doVertCoarsening, m_opt.verbosity)
                    .createMGRefSchedule(
                        topGrids, a_topOp.minBoxSize(), m_opt.maxDepth);
        } else {
            m_refSchedule =
                SemicoarseningStrategy(domainLength, m_opt.verbosity)
                    .createMGRefSchedule(
                        topGrids, a_topOp.minBoxSize(), m_opt.maxDepth);
        }

    } else {
        m_refSchedule = a_refSchedule;
    }
    CH_assert(m_refSchedule.size() > 0);
    CH_assert(m_refSchedule.back() == IntVect::Unit);

    // m_refSchedule is well defined. Set MGOptions to describe the schedule.
    m_opt.maxDepth = m_refSchedule.size() - 1;

    // Set the dXis.
    m_dXi.resize(m_refSchedule.size());
    m_dXi[0] = a_topOp.getDXi();
    for (size_t d = 1; d < m_refSchedule.size(); ++d) {
        m_dXi[d] = m_dXi[d - 1] * RealVect(m_refSchedule[d - 1]);
    }

    // Create the op hierarchy. We will own the ops at all depths.
    m_opPtrs.resize(m_opt.maxDepth + 1);
    m_opPtrs[0].reset(a_topOp.newMGOperator(IntVect::Unit));
    for (int d = 1; d <= m_opt.maxDepth; ++d) {
        m_opPtrs[d].reset(m_opPtrs[d - 1]->newMGOperator(m_refSchedule[d - 1]));
    }

    // The bottom solver.
    if (a_useBottomSolver) {
        auto ptr = new BiCGStabSolver<StateType>;
        ptr->define(m_opPtrs[m_opt.maxDepth]);
        ptr->setOptions(m_opt.bottomOptions);
        m_bottomSolverPtr.reset(ptr);
    }

    m_solverStatus.setSolverStatus(SolverStatus::UNDEFINED);
}


// -----------------------------------------------------------------------------
template <class StateType>
bool
MGSolver<StateType>::isDefined() const
{
    return (m_opPtrs.size() > 0);
}


// -----------------------------------------------------------------------------
template <class StateType>
void
MGSolver<StateType>::clear()
{
    m_bottomSolverPtr.reset();
    m_opPtrs.resize(0);
    m_refSchedule.resize(0);
    m_solverStatus.clear();
    m_opt = Options();
}


// ================================== Solvers ==================================

// -----------------------------------------------------------------------------
template <class StateType>
SolverStatus
MGSolver<StateType>::solve(StateType&       a_phi,
                           const StateType* a_crsePhiPtr,
                           const StateType& a_rhs,
                           const Real       a_time,
                           const bool       a_useHomogBCs,
                           const bool       a_setPhiToZero,
                           const Real       a_convergenceMetric)
{
    if (m_opt.numCycles < 0) {
        this->fmg(a_phi,
                  a_crsePhiPtr,
                  a_rhs,
                  a_time,
                  a_useHomogBCs,
                  a_setPhiToZero,
                  a_convergenceMetric);
    } else {
        this->vCycle(a_phi,
                     a_crsePhiPtr,
                     a_rhs,
                     a_time,
                     a_useHomogBCs,
                     a_setPhiToZero,
                     a_convergenceMetric);
    }

    return this->getSolverStatus();
}


// -----------------------------------------------------------------------------
template <class StateType>
SolverStatus
MGSolver<StateType>::vCycle(StateType&       a_phi,
                            const StateType* a_crsePhiPtr,
                            const StateType& a_rhs,
                            const Real       a_time,
                            const bool       a_useHomogBCs,
                            const bool       a_setPhiToZero,
                            const Real       a_convergenceMetric)
{
    SolverStatus& solverStatus = this->getSolverStatusRef();
    solverStatus.clear();

    Vector<Real> absResNorms(0);
    Vector<Real> relResNorms(0);

    // Create workspace.
    auto& op = m_opPtrs[0];
    StateType res, cor;
    op->create(res, a_rhs);
    op->create(cor, a_phi);

    // Initialize to zero, if requested.
    if (a_setPhiToZero) op->setToZero(a_phi);

    // If provided, use convergenceMetric as the initial residual norm.
    if (a_convergenceMetric > 0.0) {
        absResNorms.push_back(a_convergenceMetric);
        relResNorms.push_back(1.0);
    } else if (m_opt.convergenceMetric > 0.0) {
        absResNorms.push_back(m_opt.convergenceMetric);
        relResNorms.push_back(1.0);
    }

    pout() << Format::scientific;

    // Set up residual equation at depth 0.
    op->residual(
        res, a_phi, a_crsePhiPtr, a_rhs, a_time, a_useHomogBCs, a_useHomogBCs);

    // Compute initial diagnostics.
    absResNorms.push_back(op->norm(res, m_opt.normType));
    if (m_opt.verbosity > 3) {
        pout() << "Initial absolute |res| = " << absResNorms.back() << endl;
    }
    solverStatus.setInitResNorm(absResNorms[0]);

    // // Is initial guess good enough?
    // if (absResNorms.back() < m_opt.absTol) {
    //     if (m_opt.verbosity >= 1) {
    //         pout() << "Absolute convergence with initial guess. Exiting solver."
    //                << endl;
    //     }
    //     solverStatus.setFinalResNorm(absResNorms.back());
    //     solverStatus.setSolverStatus(SolverStatus::CONVERGED);
    //     return solverStatus;
    // }

    // This comes after the absTol check to avoid division by zero.
    relResNorms.push_back(absResNorms.back() / absResNorms[0]);

    // Is initial guess good enough?
    if (relResNorms.back() < m_opt.relTol) {
        if (m_opt.verbosity >= 1) {
            pout() << "Relative convergence with initial guess. Exiting solver."
                   << endl;
        }
        solverStatus.setFinalResNorm(absResNorms.back());
        solverStatus.setSolverStatus(SolverStatus::CONVERGED);
        return solverStatus;
    }

    int iter;
    for (iter = 1; iter <= m_opt.maxIters; ++iter) {
        // Set initial guess. Downward relaxation will happen in vCycle.
        op->preCond(cor, res, a_time, 0);

        // Solve residual eq. L[cor] = res.
        this->indent(0);
        this->vCycle_residualEq(cor, res, a_time, 0);
        this->unindent(0);

        // Add correction to phi at depth 0.
        op->incr(a_phi, cor, 1.0);

        // Compute new residual at depth 0.
        op->residual(res,
                     a_phi,
                     a_crsePhiPtr,
                     a_rhs,
                     a_time,
                     a_useHomogBCs,
                     a_useHomogBCs);

        // Compute diagnostics.
        absResNorms.push_back(op->norm(res, m_opt.normType));
        relResNorms.push_back(absResNorms.back() / absResNorms[0]);
        if (m_opt.verbosity >= 3) {
            pout() << "Relative |res| at iter " << iter << " = "
                   << relResNorms.back() << endl;
        }

        // Did we converge?
        if (absResNorms.back() < m_opt.absTol) {
            if (m_opt.verbosity >= 1) {
                pout() << "Absolute convergence achieved." << endl;
            }
            solverStatus.setSolverStatus(SolverStatus::CONVERGED);
            break;
        }

        if (relResNorms.back() < m_opt.relTol) {
            if (m_opt.verbosity >= 1) {
                pout() << "Relative convergence achieved." << endl;
            }
            solverStatus.setSolverStatus(SolverStatus::CONVERGED);
            break;
        }

        // Are we diverging? If so, undo last correction and exit.
        if (relResNorms[iter] > relResNorms[iter - 1]) {
            if (m_opt.verbosity >= 1) {
                pout() << "Diverging." << endl;
            }
            op->incr(a_phi, cor, -1.0);
            absResNorms.pop_back();
            relResNorms.pop_back();
            --iter;
            solverStatus.setSolverStatus(SolverStatus::DIVERGED);
            break;
        }

        // Are we hanging?
        if (relResNorms[iter] > (1.0 - m_opt.hang) * relResNorms[iter - 1]) {
            if (m_opt.verbosity >= 1) {
                pout() << "Hanging." << endl;
            }
            solverStatus.setSolverStatus(SolverStatus::HANG);
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
                   << "\t ratio = " << absResNorms[i - 1] / absResNorms[i] << "\n";
        }
        pout() << Format::unindent;

        pout() << "\nRelative convergence pattern:" << Format::indent() << endl;
        pout() << "Iter " << 0 << ": " << relResNorms[0] << "\n";
        for (int i = 1; i <= iter; ++i) {
            pout() << "Iter " << i << ": " << relResNorms[i]
                   << "\t ratio = " << relResNorms[i - 1] / relResNorms[i] << "\n";
        }
        pout() << Format::unindent << endl;
    }

    solverStatus.setFinalResNorm(absResNorms.back());
    return solverStatus;
}


// -----------------------------------------------------------------------------
template <class StateType>
SolverStatus
MGSolver<StateType>::fmg(StateType&       a_phi,
                         const StateType* a_crsePhiPtr,
                         const StateType& a_rhs,
                         const Real       a_time,
                         const bool       a_useHomogBCs,
                         const bool       a_setPhiToZero,
                         const Real       a_convergenceMetric)
{
    // If true, this will run consecutive FMG cycles.
    // If false, this will run one FMG cycle, followed by consecutive V-Cycles.
    constexpr bool allFMG = true;

    SolverStatus& solverStatus = this->getSolverStatusRef();
    solverStatus.clear();

    Vector<Real> absResNorms(0);
    Vector<Real> relResNorms(0);

    // Create workspace.
    auto& op = m_opPtrs[0];
    StateType res, cor;
    op->create(res, a_rhs);
    op->create(cor, a_phi);

    // Initialize to zero, if requested.
    if (a_setPhiToZero) op->setToZero(a_phi);

    checkForValidNAN(a_phi);
    checkForValidNAN(a_rhs);

    // If provided, use convergenceMetric as the initial residual norm.
    if (a_convergenceMetric > 0.0) {
        absResNorms.push_back(a_convergenceMetric);
        relResNorms.push_back(1.0);
    } else if (m_opt.convergenceMetric > 0.0) {
        absResNorms.push_back(m_opt.convergenceMetric);
        relResNorms.push_back(1.0);
    }

    pout() << Format::scientific;

    // Set up residual equation at depth 0.
    op->residual(
        res, a_phi, a_crsePhiPtr, a_rhs, a_time, a_useHomogBCs, a_useHomogBCs);

    // Compute initial diagnostics.
    absResNorms.push_back(op->norm(res, m_opt.normType));
    if (m_opt.verbosity > 3) {
        pout() << "Initial absolute |res| = " << absResNorms.back() << endl;
    }
    solverStatus.setInitResNorm(absResNorms[0]);

    // // Is initial guess good enough?
    // if (absResNorms.back() < m_opt.absTol) {
    //     if (m_opt.verbosity >= 1) {
    //         pout() << "Absolute convergence with initial guess. Exiting solver."
    //                << endl;
    //     }
    //     solverStatus.setFinalResNorm(absResNorms.back());
    //     solverStatus.setSolverStatus(SolverStatus::CONVERGED);
    //     return solverStatus;
    // }

    // This comes after the absTol check to avoid division by zero.
    relResNorms.push_back(absResNorms.back() / absResNorms[0]);

    // Is initial guess good enough?
    if (relResNorms.back() < m_opt.relTol) {
        if (m_opt.verbosity >= 1) {
            pout() << "Relative convergence with initial guess. Exiting solver."
                   << endl;
        }
        solverStatus.setFinalResNorm(absResNorms.back());
        solverStatus.setSolverStatus(SolverStatus::CONVERGED);
        return solverStatus;
    }

    int iter;
    for (iter = 1; iter <= m_opt.maxIters; ++iter) {
        // Solve residual eq. L[cor] = res.
        this->indent(0);
        checkForValidNAN(res);
        if (iter == 1 || allFMG) {
            this->fmg_residualEq(cor, res, a_time, 0);
        } else {
            op->preCond(cor, res, a_time, 0);
            checkForValidNAN(cor);
            this->vCycle_residualEq(cor, res, a_time, 0);
        }
        checkForValidNAN(cor);
        this->unindent(0);

        // Add correction to phi at depth 0.
        op->incr(a_phi, cor, 1.0);
        checkForValidNAN(a_phi);

        // Compute new residual at depth 0.
        op->residual(res,
                     a_phi,
                     a_crsePhiPtr,
                     a_rhs,
                     a_time,
                     a_useHomogBCs,
                     a_useHomogBCs);

        // Compute diagnostics.
        absResNorms.push_back(op->norm(res, m_opt.normType));
        relResNorms.push_back(absResNorms.back() / absResNorms[0]);
        if (m_opt.verbosity >= 3) {
            pout() << "Relative |res| at iter " << iter << " = "
                    << relResNorms.back() << endl;
        }

        // Did we converge?
        if (absResNorms.back() < m_opt.absTol) {
            if (m_opt.verbosity >= 1) {
                pout() << "Absolute convergence achieved." << endl;
            }
            solverStatus.setSolverStatus(SolverStatus::CONVERGED);
            break;
        }

        if (relResNorms.back() < m_opt.relTol) {
            if (m_opt.verbosity >= 1) {
                pout() << "Relative convergence achieved." << endl;
            }
            solverStatus.setSolverStatus(SolverStatus::CONVERGED);
            break;
        }

        // Are we diverging? If so, undo last correction and exit.
        if (relResNorms[iter] > relResNorms[iter - 1]) {
            if (m_opt.verbosity >= 1) {
                pout() << "Diverging." << endl;
            }
            op->incr(a_phi, cor, -1.0);
            absResNorms.pop_back();
            relResNorms.pop_back();
            --iter;
            solverStatus.setSolverStatus(SolverStatus::DIVERGED);
            break;
        }

        // Are we hanging?
        if (relResNorms[iter] > (1.0 - m_opt.hang) * relResNorms[iter - 1]) {
            if (m_opt.verbosity >= 1) {
                pout() << "Hanging." << endl;
            }
            solverStatus.setSolverStatus(SolverStatus::HANG);
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
                   << "\t ratio = " << absResNorms[i - 1] / absResNorms[i] << "\n";
        }
        pout() << Format::unindent;

        pout() << "\nRelative convergence pattern:" << Format::indent() << endl;
        pout() << "Iter " << 0 << ": " << relResNorms[0] << "\n";
        for (int i = 1; i <= iter; ++i) {
            pout() << "Iter " << i << ": " << relResNorms[i]
                   << "\t ratio = " << relResNorms[i - 1] / relResNorms[i] << "\n";
        }
        pout() << Format::unindent << endl;
    }

    solverStatus.setFinalResNorm(absResNorms.back());
    return solverStatus;
}


// -----------------------------------------------------------------------------
template <class StateType>
void
MGSolver<StateType>::vCycle_residualEq(StateType&       a_cor,
                                       const StateType& a_res,
                                       const Real       a_time,
                                       const int        a_depth)
{
    if (m_opt.verbosity >= 6) {
        pout() << "MG depth = " << a_depth << endl;
    }

    // Create workspace.
    auto&     op = m_opPtrs[a_depth];
    StateType tmpRes;
    op->create(tmpRes, a_res);

    // Compute rhs integral.
    if (m_opt.verbosity >= 8) {
        Real norm   = op->norm(a_res, m_opt.normType);
        // Real sumRes = Integral::sum(a_res, lg->getDXi(), 0);
        // Real sumCor = Integral::sum(a_cor, lg->getDXi(), 0);
        pout() << "|rhs| = " << norm
            //    << ", sum res = " << sumRes
            //    << ", sum cor = " << sumCor
               << endl;
    }

    if (a_depth == m_opt.maxDepth) {
        // Use bottom solver..

        // --- Bottom relaxation ---
        if (m_opt.verbosity >= 7) {
            pout() << "Bottom relax" << endl;
        }
        CH_assert(a_cor.getBoxes().compatible(a_res.getBoxes()));
        op->relax(a_cor, a_res, a_time, m_opt.numSmoothBottom);

        // Diagnostics
        if (m_opt.verbosity >= 8) {
            op->residual(tmpRes, a_cor, nullptr, a_res, a_time, true, true);
            Real norm = op->norm(tmpRes, m_opt.normType);
            pout() << "|rhs| = " << norm << endl;
        }

        // --- Bottom solver ---
        if (m_bottomSolverPtr) {
            if (m_opt.verbosity >= 7) {
                pout() << "Bottom solver" << endl;
            }
            this->indent();
            // TODO: Do something with exitStatus?
            m_bottomSolverPtr->solve(
                a_cor, nullptr, a_res, a_time, true, false);
            this->unindent();
        }

        // Diagnostics
        if (m_opt.verbosity >= 8) {
            op->residual(tmpRes, a_cor, nullptr, a_res, a_time, true, true);
            Real norm = op->norm(tmpRes, m_opt.normType);
            pout() << "|rhs| = " << norm << endl;
        }

    } else {
        // V-Cycle...

        // Create needed structures.
        auto&          crseOp  = m_opPtrs[a_depth + 1];
        const IntVect& crseRef = m_refSchedule[a_depth];
        StateType      crseCor, crseRes;
        op->createCoarsened(crseCor, a_cor, crseRef);
        op->createCoarsened(crseRes, a_res, crseRef);

        // --- Downward relaxation ---
        if (m_opt.verbosity >= 7) {
            pout() << "Smooth down" << endl;
        }
        op->relax(a_cor, a_res, a_time, m_opt.numSmoothDown);

        // Compute residual
        op->residual(tmpRes, a_cor, nullptr, a_res, a_time, true, true);

        // Diagnostics
        if (m_opt.verbosity >= 8) {
            Real norm = op->norm(tmpRes, m_opt.normType);
            pout() << "|rhs| = " << norm << endl;
        }

        // --- Restriction ---
        if (m_opt.verbosity >= 7) {
            pout() << "Restrict residual" << endl;
        }
        op->MGRestrict(crseRes, tmpRes, a_time, crseRef, *crseOp);

        // --- Coarse level solve ---

        // Set initial guess
        crseOp->preCond(crseCor, crseRes, a_time, 0);

        // Solve
        this->indent(a_depth + 1);
        const int numCycles = abs(m_opt.numCycles);
        for (int i = 0; i < numCycles; ++i) {
            this->vCycle_residualEq(crseCor, crseRes, a_time, a_depth + 1);
        }
        this->unindent(a_depth + 1);

        // --- Prolong correction ---
        if (m_opt.verbosity >= 7) {
            pout() << "Prolong and add correction" << endl;
        }
        op->MGProlong(
            a_cor, crseCor, a_time, crseRef, *crseOp, m_opt.prolongOrder);

        // Diagnostics
        if (m_opt.verbosity >= 8) {
            op->residual(tmpRes, a_cor, nullptr, a_res, a_time, true, true);
            Real norm = op->norm(tmpRes, m_opt.normType);
            pout() << "|rhs| = " << norm << endl;
        }

        // --- Upward relaxation ---
        if (m_opt.verbosity >= 7) {
            pout() << "Smooth up" << endl;
        }
        op->relax(a_cor, a_res, a_time, m_opt.numSmoothUp);

        // Diagnostics
        if (m_opt.verbosity >= 8) {
            op->residual(tmpRes, a_cor, nullptr, a_res, a_time, true, true);
            Real norm = op->norm(tmpRes, m_opt.normType);
            // Real sum  = Integral::sum(a_cor, lg->getDXi(), 0);
            pout() << "|rhs| = " << norm
                //    << ", sum cor = " << sum
                   << endl;
        }
    }
}


// -----------------------------------------------------------------------------
template <class StateType>
void
MGSolver<StateType>::fmg_residualEq(StateType&       a_cor,
                                    const StateType& a_res,
                                    const Real       a_time,
                                    const int        a_depth)
{
    CH_assert(a_cor.getBoxes().compatible(a_res.getBoxes()));

    auto& op = m_opPtrs[a_depth];

    // Initialize solution at this depth.
    op->setToZero(a_cor);

    // Go to the coarser MG level if possible.
    if (a_depth < m_opt.maxDepth) {
        // Create needed structures.
        auto&          crseOp  = m_opPtrs[a_depth + 1];
        const IntVect& crseRef = m_refSchedule[a_depth];
        StateType      crseCor, crseRes;
        op->createCoarsened(crseCor, a_cor, crseRef);
        op->createCoarsened(crseRes, a_res, crseRef);

        // --- Restrict residual ---
        if (m_opt.verbosity >= 7) {
            pout() << "Restrict residual" << endl;
        }
        op->MGRestrict(crseRes, a_res, a_time, crseRef, *crseOp);

        // --- Solve ---
        this->indent(a_depth + 1);
        this->fmg_residualEq(crseCor, crseRes, a_time, a_depth + 1);
        this->unindent(a_depth + 1);

        // --- Prolong correction ---
        if (m_opt.verbosity >= 7) {
            pout() << "Prolong and add correction" << endl;
        }
        op->MGProlong(
            a_cor, crseCor, a_time, crseRef, *crseOp, m_opt.prolongOrderFMG);

        // --- Additional relaxation for FMG prolong ---
        if (m_opt.verbosity >= 7) {
            pout() << "Smooth up for FMG prolong" << endl;
        }
        op->relax(a_cor, a_res, a_time, m_opt.numSmoothUpFMG);

        // Diagnostics
        if (m_opt.verbosity >= 8) {
            StateType tmpRes;
            op->create(tmpRes, a_res);
            op->residual(tmpRes, a_cor, nullptr, a_res, a_time, true, true);
            Real norm = op->norm(tmpRes, m_opt.normType);
            pout() << "|rhs| = " << norm << endl;
        }
    }

    // Solve with a V-cycle.
    const int numCycles = abs(m_opt.numCycles);
    for (int i = 0; i < numCycles; ++i) {
        this->vCycle_residualEq(a_cor, a_res, a_time, a_depth);
    }
}


// -----------------------------------------------------------------------------
template <class StateType>
void
MGSolver<StateType>::indent(const int a_depth) const
{
    if (m_opt.verbosity >= 6) {
        if (a_depth < 0) {
            pout() << Format::indent();
        } else {
            std::string bullet = std::to_string(a_depth) + std::string(": ");
            pout() << Format::indent(4, bullet.c_str());
        }
    }
}


// -----------------------------------------------------------------------------
template <class StateType>
void
MGSolver<StateType>::unindent(const int /*a_depth*/) const
{
    if (m_opt.verbosity >= 6) {
        pout() << Format::unindent;
    }
}


}; // namespace Elliptic
