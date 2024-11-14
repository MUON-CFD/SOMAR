#ifndef Hb85c5095bf0dc1fbde78ea67277de9cb
#error AMRMGSolverI.H should only be used by AMRMGSolver.H
#endif

#include "AnisotropicRefinementTools.H"


namespace Elliptic {

// ======================= AMRMGSolver Implementation ==========================

// -----------------------------------------------------------------------------
// Default constructor. Leaves object unusable.
// -----------------------------------------------------------------------------
template <class StateType>
AMRMGSolver<StateType>::AMRMGSolver()
: m_isDefined(false)
, m_opt()
, m_solverStatus()
, m_lbase(size_t(-1))
, m_lmin(size_t(-1))
, m_lmax(size_t(-1))
{
    this->setDefaultOptions();
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
template <class StateType>
AMRMGSolver<StateType>::~AMRMGSolver()
{
    this->clear();
}


// -----------------------------------------------------------------------------
// Set up this solver to solve L[phi] = rhs on levels a_lmin to a_lmax.
// If a_lmin > 0, then a_vAMRMGOps must also be defined on a_lmin - 1.
// a_vAMRMGOps can be defined or undefined on the levels we don't use.
// -----------------------------------------------------------------------------
template <class StateType>
void
AMRMGSolver<StateType>::define(
    Vector<std::shared_ptr<const AMRMGOpType>> a_vAMRMGOps,
    const size_t                               a_lmin,
    const size_t                               a_lmax,
    Options                                    a_opt)
{
    m_opt = a_opt;
    this->define(a_vAMRMGOps, a_lmin, a_lmax);
}


// -----------------------------------------------------------------------------
// Set up this solver to solve L[phi] = rhs on levels a_lmin to a_lmax.
// If a_lmin > 0, then a_vAMRMGOps must also be defined on a_lmin - 1.
// a_vAMRMGOps can be defined or undefined on the levels we don't use.
// This version does not set Options, just in case you already did that.
// -----------------------------------------------------------------------------
template <class StateType>
void
AMRMGSolver<StateType>::define(
    Vector<std::shared_ptr<const AMRMGOpType>> a_vAMRMGOps,
    const size_t                               a_lmin,
    const size_t                               a_lmax)
{
    CH_assert(a_lmax >= a_lmin);

    m_lbase = ((a_lmin > 0) ? (a_lmin - 1) : a_lmin);
    m_lmin  = a_lmin;
    m_lmax  = a_lmax;

    // Collect the AMRMG ops.
    m_vAMRMGOps.resize(m_lmax + 1);
    for (size_t l = m_lbase; l <= m_lmax; ++l) {
        CH_assert(a_vAMRMGOps[l]);
        m_vAMRMGOps[l] = a_vAMRMGOps[l];
    }

    // Compute the ref ratios.
    m_vAMRCrseRefRatios = Vector<IntVect>(m_lmax + 1, IntVect::Zero);
    m_vAMRFineRefRatios = Vector<IntVect>(m_lmax + 1, IntVect::Zero);
    for (size_t l = m_lbase + 1; l <= m_lmax; ++l) {
        Box crseBox = a_vAMRMGOps[l - 1]->getBoxes().physDomain().domainBox();
        Box fineBox = a_vAMRMGOps[l]->getBoxes().physDomain().domainBox();

        m_vAMRCrseRefRatios[l] = calculateRefinementRatio(crseBox, fineBox);
        m_vAMRFineRefRatios[l - 1] = m_vAMRCrseRefRatios[l];

        CH_assert(m_vAMRCrseRefRatios[l] >= IntVect::Unit);
        CH_assert(m_vAMRCrseRefRatios[l].product() > 1);
    }
    m_vAMRFineRefRatios[m_lmax] = IntVect::Unit;

    // Define the mini VCycle solvers, where needed.
    m_vMiniVCycleSolvers.resize(m_lmax + 1);
    for (size_t l = m_lmin + 1; l <= m_lmax; ++l) {
        const IntVect&           crseRefRatio = m_vAMRCrseRefRatios[l];
        const IntVect            minBoxSize   = m_vAMRMGOps[l]->minBoxSize();
        const DisjointBoxLayout& grids        = m_vAMRMGOps[l]->getBoxes();
        const RealVect&          dXi          = m_vAMRMGOps[l]->getDXi();

        const IntVect nx = grids.physDomain().size();
        const RealVect L = RealVect(nx) * dXi;

        if (crseRefRatio <= 2*IntVect::Unit) {
            // No mini V-cycle.
        } else {
            // Mini V-cycle.
            MiniVCycleStrategy mgStrategy(
                L, crseRefRatio, m_opt.mgOptions.verbosity);
            Vector<IntVect> mgSchedule =
                mgStrategy.createMGRefSchedule(grids, minBoxSize);

            typename MGSolverType::Options miniVOpts = m_opt.mgOptions;
            miniVOpts.verbosity = 0;
            miniVOpts.maxDepth = mgSchedule.size() - 1;
            miniVOpts.maxIters = 1;
            miniVOpts.absTol   = 1.0e-30;
            miniVOpts.relTol   = 1.0e-30;

            const bool useBottomSolver = false;
            m_vMiniVCycleSolvers[l].reset(new MGSolverType);
            m_vMiniVCycleSolvers[l]->define(
                *m_vAMRMGOps[l], miniVOpts, mgSchedule, useBottomSolver);
        }
    } // l

    // Bottom solver
    typename MGSolverType::Options bottomMGOptions = m_opt.mgOptions;
    // Be as loud as any bottom solver.
    bottomMGOptions.verbosity = m_opt.mgOptions.bottomOptions.verbosity;
    m_bottomMGSolver.reset(new MGSolverType);
    m_bottomMGSolver->define(*m_vAMRMGOps[m_lmin], bottomMGOptions);

    // This object is ready for use.
    m_isDefined = true;
}


// -----------------------------------------------------------------------------
// Free's memory. Leaves object unusable.
// -----------------------------------------------------------------------------
template <class StateType>
void
AMRMGSolver<StateType>::clear()
{
    m_bottomMGSolver.reset();
    m_vMiniVCycleSolvers.resize(0);
    m_vAMRMGOps.resize(0);

    m_lmax  = size_t(-1);
    m_lmin  = size_t(-1);
    m_lbase = size_t(-1);

    m_solverStatus.clear();
    // Do not reset options.

    m_isDefined = false;
}


// -----------------------------------------------------------------------------
template <class StateType>
void
AMRMGSolver<StateType>::setDefaultOptions()
{
    const auto& proj = ProblemContext::getInstance()->proj;

    m_opt.absTol            = proj.absTol;
    m_opt.relTol            = proj.relTol;
    m_opt.convergenceMetric = -1.0;
    m_opt.numCycles         = proj.numCycles;
    m_opt.maxIters          = proj.maxIters;
    m_opt.hang              = proj.hang;

    m_opt.mgOptions.absTol            = proj.absTol;
    m_opt.mgOptions.relTol            = proj.relTol;
    m_opt.mgOptions.convergenceMetric = -1.0;
    m_opt.mgOptions.numSmoothDown     = proj.numSmoothDown;
    m_opt.mgOptions.numSmoothUp       = proj.numSmoothUp;
    m_opt.mgOptions.numSmoothBottom   = proj.numSmoothBottom;
    m_opt.mgOptions.numSmoothPrecond  = proj.numSmoothPrecond;
    m_opt.mgOptions.prolongOrder      = proj.prolongOrder;
    m_opt.mgOptions.prolongOrderFMG   = proj.prolongOrderFMG;
    m_opt.mgOptions.numSmoothUpFMG    = proj.numSmoothUpFMG;
    m_opt.mgOptions.maxDepth          = proj.maxDepth;
    m_opt.mgOptions.numCycles         = proj.numCycles;
    m_opt.mgOptions.maxIters          = proj.maxIters;
    m_opt.mgOptions.hang              = proj.hang;
    m_opt.mgOptions.normType          = proj.normType;
    m_opt.mgOptions.verbosity         = proj.verbosity;

    m_opt.mgOptions.bottomOptions.absTol           = proj.bottom_absTol;
    m_opt.mgOptions.bottomOptions.relTol           = proj.bottom_relTol;
    m_opt.mgOptions.bottomOptions.small            = proj.bottom_small;
    m_opt.mgOptions.bottomOptions.hang             = proj.bottom_hang;
    m_opt.mgOptions.bottomOptions.maxIters         = proj.bottom_maxIters;
    m_opt.mgOptions.bottomOptions.maxRestarts      = proj.bottom_maxRestarts;
    m_opt.mgOptions.bottomOptions.normType         = proj.bottom_normType;
    m_opt.mgOptions.bottomOptions.verbosity        = proj.bottom_verbosity;
    m_opt.mgOptions.bottomOptions.numSmoothPrecond = proj.bottom_numSmoothPrecond;
}


// -----------------------------------------------------------------------------
// @brief      Solves L[phi] = rhs using the most appropriate methods.

// @param[inout] a_phi                On entry, the initial guess.
//                                    On exit, the solution.
// @param[in]    a_rhs
// @param[in]    a_time               Used to set BCs.
// @param[in]    a_useHomogBCs
// @param[in]    a_setPhiToZero       Should we set initial guess to zero?
// @param[in]    a_convergenceMetric  If supplied, this will be used as the
//                                    initial residual when computing
//                                    relative residuals. This setting takes
//                                    precedence over
//                                    Options::convergenceMetric.

// @return     The SolverStatus object.
//             This tells the caller if we converged.
// -----------------------------------------------------------------------------
template <class StateType>
SolverStatus
AMRMGSolver<StateType>::solve(Vector<StateType*>&             a_vphi,
                              const Vector<const StateType*>& a_vrhs,
                              const Real                      a_time,
                              const bool                      a_useHomogBCs,
                              const bool                      a_setPhiToZero,
                              const Real a_convergenceMetric) const
{
#ifndef NDEBUG
    CH_assert(m_isDefined);
    CH_assert(a_vphi[m_lbase]);
    for (size_t l = m_lmin; l <= m_lmax; ++l) {
        CH_assert(a_vphi[l]);
        CH_assert(a_vrhs[l]);
    }
#endif

    pout() << Format::pushFlags << Format::scientific;

    SolverStatus& solverStatus = this->getSolverStatusRef();
    solverStatus.clear();

    // Allocation
    this->initializeSolve(a_vphi, a_vrhs);

    // ................... Set up initial residual eq ..........................

    // Initialize phi, if requested.
    if (a_setPhiToZero) {
        for (size_t l = m_lbase; l <= m_lmax; ++l) {
            m_vAMRMGOps[l]->setToZero(*a_vphi[l]);
        }
    }

    // Initialize residual.
    Real resNorm = this->computeAMRResidual(
        m_vr, a_vphi, a_vrhs, a_time, a_useHomogBCs, true);

    // If provided, use convergenceMetric as the initial residual norm.
    Real initResNorm = resNorm;
    if (a_convergenceMetric > 0.0) {
        initResNorm = a_convergenceMetric;
    } else if (m_opt.convergenceMetric > 0.0) {
        initResNorm = m_opt.convergenceMetric;
    }

    // Report initial diagnostics.
    Vector<Real> absResNorm(1, resNorm);
    Vector<Real> relResNorm(1, 1.0);
    if (m_opt.verbosity >= 3) {
        pout() << "Initial absolute |res| = " << absResNorm.back() <<
        endl;
    }
    solverStatus.setInitResNorm(absResNorm[0]);

    // Is initial guess good enough? Check abs norm.
    if (absResNorm.back() < m_opt.absTol) {
        if (m_opt.verbosity >= 1) {
            pout() << "Absolute convergence with initial guess. Exiting solver."
                   << endl;
        }

        solverStatus.setFinalResNorm(absResNorm.back());
        solverStatus.setSolverStatus(SolverStatus::CONVERGED);
        this->finalizeSolve();
        pout() << Format::popFlags;
        return solverStatus;
    }

    // ........................ Solve / main loop ..............................

    // Call VCycle on L[e] = r and update vphi += e, vres = vrhs - L[vphi].
    int iter;
    for (iter = 1; iter <= m_opt.maxIters; ++iter) {
        // pout() << "iter = " << iter << endl;

        // Set initial guess.
        for (size_t l = m_lbase; l <= m_lmax; ++l) {
            m_vAMRMGOps[l]->setToZero(*m_ve[l]);
        }

        // Relax L[e] = r.
        this->indent(m_lmax);
        this->amrVCycle_residualEq(m_ve, m_vr, a_time, m_lmax);
        this->unindent(m_lmax);

        // vphi += e
        for (size_t l = m_lmin; l <= m_lmax; ++l) {
            m_vAMRMGOps[l]->incr(*a_vphi[l], *m_ve[l], 1.0);
        }

        // CF consistency sweep for ops with orderOfAccuracy() > 2.
        // Apply correction from lmax - 1 -> lmin.
        if (m_vAMRMGOps[m_lmin]->orderOfAccuracy() > 2) {
            for (size_t l = m_lmax; l > m_lmin; --l) {
                m_vAMRMGOps[l]->enforceCFConsistency(*a_vphi[l - 1],
                                                     *a_vphi[l]);
            }
        }

        // Compute new AMR residual.
        resNorm = this->computeAMRResidual(
            m_vr, a_vphi, a_vrhs, a_time, a_useHomogBCs, true);

        // Report diagnostics.
        absResNorm.push_back(resNorm);
        relResNorm.push_back(absResNorm.back() / initResNorm);
        if (m_opt.verbosity >= 3) {
            pout() << "Relative |res| after iter " << iter << " = "
                   << relResNorm.back() << endl;
        }

        // Did we converge?
        if (absResNorm.back() < m_opt.absTol) {
            if (m_opt.verbosity >= 1) {
                pout() << "Absolute convergence achieved." << endl;
            }
            solverStatus.setSolverStatus(SolverStatus::CONVERGED);
            break;
        }

        if (relResNorm.back() < m_opt.relTol) {
            if (m_opt.verbosity >= 1) {
                pout() << "Relative convergence achieved." << endl;
            }
            solverStatus.setSolverStatus(SolverStatus::CONVERGED);
            break;
        }

        // Are we diverging? If so, undo last correction and exit.
        if (relResNorm[iter] > relResNorm[iter - 1]) {
            if (m_opt.verbosity >= 1) {
                pout() << "Diverging." << endl;
            }
            if (m_vAMRMGOps[m_lmin]->orderOfAccuracy() <= 2) {
                for (size_t l = m_lmin; l <= m_lmax; ++l) {
                    m_vAMRMGOps[l]->incr(*a_vphi[l], *m_ve[l], -1.0);
                }
            } else {
                // We need to remove the CF consistency update before
                // removing the faulty increment.
                UNDEFINED_FUNCTION();
            }
            absResNorm.pop_back();
            relResNorm.pop_back();
            --iter;
            solverStatus.setSolverStatus(SolverStatus::DIVERGED);
            break;
        }

        // Are we hanging?
        if (relResNorm[iter] > (1.0 - m_opt.hang) * relResNorm[iter - 1]) {
            if (m_opt.verbosity >= 1) {
                pout() << "Hanging." << endl;
            }
            solverStatus.setSolverStatus(SolverStatus::HANG);
            break;
        }
    }  // iter

    // ....................... Report final results ............................

    // If m_opt.maxIters was reached, tell the caler and reset iter.
    if (iter == m_opt.maxIters) {
        if (m_opt.verbosity >= 1) {
            pout() << "Max iters reached." << endl;
        }
        solverStatus.setSolverStatus(SolverStatus::MAXITERS);
    }
    iter = min(iter, m_opt.maxIters);

    if (m_opt.verbosity >= 5) {
        pout() << "\nAbsolute convergence pattern:" << Format::indent() << endl;
        pout() << "Iter " << 0 << ": " << absResNorm[0] << "\n";
        for (int i = 1; i <= iter; ++i) {
            pout() << "Iter " << i << ": " << absResNorm[i]
                   << "\t ratio = " << absResNorm[i - 1] / absResNorm[i]
                   << "\n";
        }
        pout() << Format::unindent;

        pout() << "\nRelative convergence pattern:" << Format::indent() << endl;
        pout() << "Iter " << 0 << ": " << relResNorm[0] << "\n";
        for (int i = 1; i <= iter; ++i) {
            pout() << "Iter " << i << ": " << relResNorm[i]
                   << "\t ratio = " << relResNorm[i - 1] / relResNorm[i]
                   << "\n";
        }
        pout() << Format::unindent << endl;
    }

    this->finalizeSolve();
    pout() << Format::popFlags;
    return solverStatus;
}


// -----------------------------------------------------------------------------
// Relaxes L[phi] = rhs over the AMR hierarchy with a V-Cycle.
// It is assumed the problem is homogeneous, so send in the residual eq.
// a_vphi needs to be set to zero before entry.
// -----------------------------------------------------------------------------
template <class StateType>
void
AMRMGSolver<StateType>::amrVCycle_residualEq(
    Vector<unique_ptr<StateType>>& a_vphi,
    Vector<unique_ptr<StateType>>& a_vrhs,
    const Real                     a_time,
    const size_t                   a_lev) const
{
    const int prolongOrderAMR = 0; // TODO
    const bool homogPhysBCs = true;

    // Some shortcuts.
    auto&      op  = m_vAMRMGOps[a_lev];
    StateType& phi = *a_vphi[a_lev];
    StateType& rhs = *a_vrhs[a_lev];

    if (a_lev > m_lmin) {
        // Some more shortcuts.
        auto&          crseOp      = m_vAMRMGOps[a_lev - 1];
        StateType&     crsePhi     = *a_vphi[a_lev - 1];
        StateType&     crseRhs     = *a_vrhs[a_lev - 1];
        StateType&     scratch     = *m_vScratchPhi[a_lev];
        StateType&     crseScratch = *m_vScratchPhi[a_lev - 1];
        StateType&     resC        = *m_vResC[a_lev];  // on coarsened res grids
        const IntVect& refRatio    = m_vAMRCrseRefRatios[a_lev];

        const StateType* reallyCrsePhiPtr = nullptr;
        if (a_lev >= 2) reallyCrsePhiPtr = &*a_vphi[a_lev - 2];

        // --- Smooth down ---
        this->smoothDown(phi, rhs, a_time, a_lev);

        // --- Restrict residual ---
        {
            // Recompute residual on entire coarser level.
            // res[l-1] = rhs[l-1] - Lcomp[l-1](phi[l], phi[l-1], phi[l-2]) on G[l-1],
            // where G[l-1] are the grids on level l-1.
            crseOp->assignLocal(crseScratch, crseRhs);
            this->computeAMRResidualLevel(crseRhs,
                                          &phi,
                                          crsePhi,
                                          reallyCrsePhiPtr,
                                          crseScratch,
                                          a_lev - 1,
                                          a_time,
                                          homogPhysBCs);

            // pout() << "* Crse level residual before restriction:\n";
            // crseOp->levelEquationIsConsistent(crsePhi,
            //                                   reallyCrsePhiPtr,
            //                                   crseRhs,
            //                                   a_time,
            //                                   homogPhysBCs,
            //                                   false);

            // Restrict residual where fine grids exist.
            // Res[l-1] = Restrict( Res[l] - L[l](Cor[l], Cor[l-1]) ) on PG[l],
            // where PG[l] is the projection of the level l grids onto level l-1.
            //
            // We do this in 3 steps:
            //  1. Compute tmpRes[l] = Res[l] - L[l](Cor[l], Cor[l-1]) on G[l].
            //  2. resC[l] = MGRestrict( tmpRes[l] ), sending tmpRes[l] to PG[l].
            //  3. Copy resC[l] -> Res[l] only on PG[l].
            op->AMRResidualNF(
                scratch, phi, crsePhi, rhs, refRatio, a_time, homogPhysBCs);

            // pout() << "* Fine level residual being restricted:\n";
            // op->levelEquationIsConsistent(
            //     phi, &crsePhi, scratch, a_time, homogPhysBCs, false);

            op->MGRestrict(resC, scratch, a_time, refRatio, *crseOp);

            crseOp->assign(crseRhs, resC, &*m_vResCopier[a_lev]);

            // pout() << "* Crse level residual after restriction:\n";
            // crseOp->levelEquationIsConsistent(crsePhi,
            //                                   reallyCrsePhiPtr,
            //                                   crseRhs,
            //                                   a_time,
            //                                   homogPhysBCs,
            //                                   false);
        }

        // --- Coarse solve ---
        const int numCycles = abs(m_opt.numCycles);
        for (int cycle = 0; cycle < numCycles; ++cycle) {
            this->indent(a_lev - 1);
            this->amrVCycle_residualEq(a_vphi, a_vrhs, a_time, a_lev - 1);
            this->unindent(a_lev - 1);
        }

        // --- Prolong increment and update residual ---
        {
            // Prolong
            crseOp->assign(resC, crsePhi, &*m_vReverseCopier[a_lev]);
            op->MGProlong(
                phi, resC, a_time, refRatio, *crseOp, prolongOrderAMR);

            // Update residual.
            // Res[l] = Res[l] - L[l](Cor[l], Cor[l-1])
            op->AMRUpdateResidualNF(
                rhs, phi, crsePhi, scratch, refRatio, a_time, homogPhysBCs);
        }

        // --- Smooth up ---
        // Do this in incremental form so CFBCs are homog during relaxation.
        op->setToZero(scratch);
        this->smoothUp(scratch, rhs, a_time, a_lev);
        op->incr(phi, scratch, 1.0);

    } else {
        // --- Bottom solve ---

        // // Diagnostics
        // if (doDiagnostics) {
        //     diagnostics(op, tmpRes, cor, crseCorPtr, res);
        //     pout() << "AMR Bottom solver:" << endl;
        // }
        // pout() << "* Residual on coarse level:\n";
        // op->levelEquationIsConsistent(
        //     phi, nullptr, rhs, a_time, homogPhysBCs, true);

        // op->removeKernel(rhs);  // TODO: shouldn't be needed.

        // m_bottomMGSolver->vCycle_residualEq(phi, rhs, a_time, 0); // Just one iter!
        m_bottomMGSolver->solve(phi, nullptr, rhs, a_time, homogPhysBCs, false);

        // op->preCond(phi, rhs, a_time, 10);
        // BiCGStabSolver<StateType> solver;
        // solver.define(op);
        // solver.options().absTol    = 1e-12;
        // solver.options().relTol    = 1e-10;
        // solver.options().verbosity = 3;
        // solver.solve(phi, nullptr, rhs, a_time, true, false);

        // if (doDiagnostics) {
        //     diagnostics(op, tmpRes, cor, crseCorPtr, res);
        // }
    }
}


// -----------------------------------------------------------------------------
// A solve requires a number of temporaries over the entire AMR hierarchy.
// Allocation + definition is potentially expensive, so I shoved it all
// into this function. If the expense is high, we'll see it in our
// profilers.
// -----------------------------------------------------------------------------
template <class StateType>
void
AMRMGSolver<StateType>::initializeSolve(
    const Vector<StateType*>&       a_vphi,
    const Vector<const StateType*>& a_vrhs) const
{
    CH_assert(m_isDefined); // For m_lbase, m_lmin, and m_lmax.

    m_ve.resize(m_lmax + 1);
    m_vr.resize(m_lmax + 1);
    m_vScratchPhi.resize(m_lmax + 1);
    m_vResC.resize(m_lmax + 1);
    m_vResCopier.resize(m_lmax + 1);
    m_vReverseCopier.resize(m_lmax + 1);

    for (size_t l = m_lbase; l <= m_lmax; ++l) {
        m_ve[l].reset(new StateType);
        m_vAMRMGOps[l]->create(*m_ve[l], *a_vphi[l]);
    }

    for (size_t l = m_lmin; l <= m_lmax; ++l) {
        const auto     op = m_vAMRMGOps[l];
        const IntVect& r  = m_vAMRCrseRefRatios[l];

        m_vr[l].reset(new StateType);
        op->create(*m_vr[l], *a_vrhs[l]);

        m_vScratchPhi[l].reset(new StateType);
        op->create(*m_vScratchPhi[l], *a_vphi[l]);

        if (l > m_lmin) {
            m_vResC[l].reset(new StateType);
            op->createCoarsened(*m_vResC[l], *a_vrhs[l], r);

            m_vResCopier[l].reset(new CopierType);
            op->buildCopier(*m_vResCopier[l], *a_vrhs[l-1], *m_vResC[l]);

            m_vReverseCopier[l].reset(new CopierType);
            op->buildReverseCopier(*m_vReverseCopier[l],
                                   *m_vResCopier[l],
                                   *a_vrhs[l - 1],
                                   *m_vResC[l]);
        }
    }
}


// -----------------------------------------------------------------------------
// Deallocation of temporaries.
// -----------------------------------------------------------------------------
template <class StateType>
void
AMRMGSolver<StateType>::finalizeSolve() const
{
    for (size_t l = m_lbase; l <= m_lmax; ++l) {
        m_vAMRMGOps[l]->clear(*m_ve[l]);
    }

    for (size_t l = m_lmin; l <= m_lmax; ++l) {
        m_vAMRMGOps[l]->clear(*m_vr[l]);
        m_vAMRMGOps[l]->clear(*m_vScratchPhi[l]);
        if (l > m_lmin) {
            m_vAMRMGOps[l]->clear(*m_vResC[l]);
        }
    }

    m_vReverseCopier.resize(0);
    m_vResCopier.resize(0);
    m_vResC.resize(0);
    m_vScratchPhi.resize(0);
    m_vr.resize(0);
    m_ve.resize(0);
}


// -----------------------------------------------------------------------------
// Computes lhs = rhs - L[phi] from m_lmin to m_lmax.
// This calls computeAMRResidualLevel on each level and computes the
// inf-norm.
// -----------------------------------------------------------------------------
template <class StateType>
Real
AMRMGSolver<StateType>::computeAMRResidual(
    Vector<unique_ptr<StateType>>&  a_vres,
    Vector<StateType*>&             a_vphi,
    const Vector<const StateType*>& a_vrhs,
    const Real                      a_time,
    const bool                      a_useHomogBCs,
    const bool                      a_computeNorm) const
{
    Real      rnorm     = 0.0;
    Real      localNorm = 0.0;
    const int p         = m_opt.normType;

    for (size_t l = m_lmin; l <= m_lmax; ++l) {
        StateType* finePhiPtr = nullptr;
        if (l < m_lmax) finePhiPtr = a_vphi[l + 1];

        const StateType* crsePhiPtr = nullptr;
        if (l > 0) crsePhiPtr = a_vphi[l - 1];

        this->computeAMRResidualLevel(*a_vres[l],
                                      finePhiPtr,
                                      *a_vphi[l],
                                      crsePhiPtr,
                                      *a_vrhs[l],
                                      l,
                                      a_time,
                                      a_useHomogBCs);

        if (a_computeNorm) {
            if (l == m_lmax) {
                localNorm = m_vAMRMGOps[l]->AMRNormLevel(
                    *a_vres[l], nullptr, IntVect::Unit, p);
            } else {
                localNorm = m_vAMRMGOps[l]->AMRNormLevel(
                    *a_vres[l], &*a_vres[l + 1], m_vAMRFineRefRatios[l], p);
            }

            if (p == 0) {
                rnorm = max(rnorm, localNorm);
            } else {
                rnorm = rnorm + pow(localNorm, p);
            }
        }
    }

    if (p != 0) {
        rnorm = pow(rnorm, 1.0 / Real(p));
    }

    return rnorm;
}


// -----------------------------------------------------------------------------
// Computes the residual at a single level using the appropriate version
// of AMROperator.
// a_finePhiPtr ghosts may be reset.
// -----------------------------------------------------------------------------
template <class StateType>
void
AMRMGSolver<StateType>::computeAMRResidualLevel(StateType&       a_res,
                                                StateType*       a_finePhiPtr,
                                                StateType&       a_phi,
                                                const StateType* a_crsePhiPtr,
                                                const StateType& a_rhs,
                                                const size_t     a_lev,
                                                const Real       a_time,
                                                const bool a_useHomogBCs) const
{
    const IntVect& crseRefRatio = m_vAMRCrseRefRatios[a_lev];
    const IntVect& fineRefRatio = m_vAMRFineRefRatios[a_lev];

    // AMR operator, L[phi].
    if (m_lmax != m_lmin) {
        if (a_lev == m_lmax) {
            // Nothing above, something below.
            m_vAMRMGOps[a_lev]->AMRResidualNF(a_res,
                                              a_phi,
                                              *a_crsePhiPtr,
                                              a_rhs,
                                              crseRefRatio,
                                              a_time,
                                              a_useHomogBCs);
        } else if (a_lev == m_lmin) {
            if (a_lev == 0) {
                // Something above, nothing below.
                m_vAMRMGOps[a_lev]->AMRResidualNC(a_res,
                                                  *a_finePhiPtr,
                                                  a_phi,
                                                  a_rhs,
                                                  fineRefRatio,
                                                  a_time,
                                                  a_useHomogBCs,
                                                  *m_vAMRMGOps[a_lev + 1]);
            } else {
                // Something above, something below.
                m_vAMRMGOps[a_lev]->AMRResidual(a_res,
                                                *a_finePhiPtr,
                                                a_phi,
                                                *a_crsePhiPtr,
                                                a_rhs,
                                                fineRefRatio,
                                                crseRefRatio,
                                                a_time,
                                                a_useHomogBCs,
                                                *m_vAMRMGOps[a_lev + 1]);
            }
        } else {
            // Something above, something below.
            m_vAMRMGOps[a_lev]->AMRResidual(a_res,
                                            *a_finePhiPtr,
                                            a_phi,
                                            *a_crsePhiPtr,
                                            a_rhs,
                                            fineRefRatio,
                                            crseRefRatio,
                                            a_time,
                                            a_useHomogBCs,
                                            *m_vAMRMGOps[a_lev + 1]);
        }
    } else {
        CH_assert(a_lev == m_lmin);
        if (a_lev == 0) {
            // Nothing above, nothing below.
            m_vAMRMGOps[a_lev]->residual(a_res,
                                         a_phi,
                                         nullptr,
                                         a_rhs,
                                         a_time,
                                         a_useHomogBCs,
                                         a_useHomogBCs);
        } else {
            // Nothing above, something below.
            m_vAMRMGOps[a_lev]->AMRResidualNF(a_res,
                                              a_phi,
                                              *a_crsePhiPtr,
                                              a_rhs,
                                              crseRefRatio,
                                              a_time,
                                              a_useHomogBCs);
        }
    }
}


// -----------------------------------------------------------------------------
// Eliminate high wavenumber errors on a_lev that will not be captured by
// the grids on a_lev - 1. This either calls the MG relax function or
// performs a mini V-cycle.
// -----------------------------------------------------------------------------
template <class StateType>
void
AMRMGSolver<StateType>::smoothDown(StateType&   a_cor,
                                   StateType&   a_res,
                                   const Real   a_time,
                                   const size_t a_lev) const
{
    // Diagnostics
    if (m_opt.verbosity >= 7) {
        auto& op = m_vAMRMGOps[a_lev];

        StateType tmpCor, tmpRes;
        op->create(tmpCor, a_cor);
        op->create(tmpRes, a_res);

        op->assignLocal(tmpCor, a_cor);
        op->residual(tmpRes, tmpCor, nullptr, a_res, a_time, true, true);
        Real norm = op->norm(tmpRes, m_opt.normType);
        pout() << "Pre  smooth |rhs| = " << norm << endl;

        op->clear(tmpRes);
        op->clear(tmpCor);
    }

    if (m_vMiniVCycleSolvers[a_lev]) {
        if (m_opt.verbosity >= 7) {
            pout() << "Begin mini V-Cycle:" << endl;
        }
        m_vMiniVCycleSolvers[a_lev]->vCycle_residualEq(a_cor, a_res, a_time, 0);
        if (m_opt.verbosity >= 7) {
            pout() << "End mini V-Cycle." << endl;
        }
    } else {
        m_vAMRMGOps[a_lev]->relax(
            a_cor, a_res, a_time, m_opt.mgOptions.numSmoothDown);
    }

    // Diagnostics
    if (m_opt.verbosity >= 7) {
        auto& op = m_vAMRMGOps[a_lev];

        StateType tmpCor, tmpRes;
        op->create(tmpCor, a_cor);
        op->create(tmpRes, a_res);

        op->assignLocal(tmpCor, a_cor);
        op->residual(tmpRes, tmpCor, nullptr, a_res, a_time, true, true);
        Real norm = op->norm(tmpRes, m_opt.normType);
        pout() << "Post smooth |rhs| = " << norm << endl;

        op->clear(tmpRes);
        op->clear(tmpCor);
    }
}


// -----------------------------------------------------------------------------
// Eliminate high wavenumber errors on a_lev that were not captured by
// the grids on a_lev - 1. This either calls the MG relax function or
// performs a mini V-cycle.
// -----------------------------------------------------------------------------
template <class StateType>
void
AMRMGSolver<StateType>::smoothUp(StateType&   a_cor,
                                 StateType&   a_res,
                                 const Real   a_time,
                                 const size_t a_lev) const
{
    // Diagnostics
    if (m_opt.verbosity >= 7) {
        auto& op = m_vAMRMGOps[a_lev];

        StateType tmpCor, tmpRes;
        op->create(tmpCor, a_cor);
        op->create(tmpRes, a_res);

        op->assignLocal(tmpCor, a_cor);
        op->residual(tmpRes, tmpCor, nullptr, a_res, a_time, true, true);
        Real norm = op->norm(tmpRes, m_opt.normType);
        pout() << "Pre  smooth |rhs| = " << norm << endl;

        op->clear(tmpRes);
        op->clear(tmpCor);
    }

    if (m_vMiniVCycleSolvers[a_lev]) {
        if (m_opt.verbosity >= 7) {
            pout() << "Begin mini V-Cycle:" << endl;
        }
        m_vMiniVCycleSolvers[a_lev]->vCycle_residualEq(a_cor, a_res, a_time, 0);
        if (m_opt.verbosity >= 7) {
            pout() << "End mini V-Cycle." << endl;
        }
    } else {
        m_vAMRMGOps[a_lev]->relax(
            a_cor, a_res, a_time, m_opt.mgOptions.numSmoothUp);
    }

    // Diagnostics
    if (m_opt.verbosity >= 7) {
        auto& op = m_vAMRMGOps[a_lev];

        StateType tmpCor, tmpRes;
        op->create(tmpCor, a_cor);
        op->create(tmpRes, a_res);

        op->assignLocal(tmpCor, a_cor);
        op->residual(tmpRes, tmpCor, nullptr, a_res, a_time, true, true);
        Real norm = op->norm(tmpRes, m_opt.normType);
        pout() << "Post smooth |rhs| = " << norm << endl;

        op->clear(tmpRes);
        op->clear(tmpCor);
    }
}


// -----------------------------------------------------------------------------
// pout formatting utility.
// -----------------------------------------------------------------------------
template <class StateType>
void
AMRMGSolver<StateType>::indent(const int a_lev) const
{
    using std::string;
    using std::to_string;

    if (m_opt.verbosity >= 7) {
        if (a_lev < 0) {
            pout() << Format::indent();
        } else {
            string bullet = string("(") + to_string(a_lev) + string(",0): ");
            pout() << Format::indent(8, bullet.c_str());
        }
    }
}


// -----------------------------------------------------------------------------
// pout formatting utility.
// -----------------------------------------------------------------------------
template <class StateType>
void
AMRMGSolver<StateType>::unindent(const int /*a_lev*/) const
{
    if (m_opt.verbosity >= 7) {
        pout() << Format::unindent;
    }
}


}; // namespace Elliptic
