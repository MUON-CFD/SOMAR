#include "LevelHybridSolver.H"

namespace Elliptic {

// ======================== Option setters / getters ===========================

// -----------------------------------------------------------------------------
LevelHybridSolver::Options
LevelHybridSolver::getDefaultOptions()
{
    // BUG: A lot of this is hard-coded.

    Options opt;
    const auto& proj = ProblemContext::getInstance()->proj;

    opt.absTol         = proj.absTol;
    opt.relTol         = proj.relTol;
    opt.maxSolverSwaps = 10;
    opt.normType       = proj.normType;
    opt.verbosity      = proj.verbosity;
    opt.hang           = proj.hang;

    opt.mgOptions = MGSolver<StateType>::getDefaultOptions();
    opt.mgOptions.normType  = opt.normType;
    opt.mgOptions.verbosity = 0;
    // Use default smoothing options.

    opt.mgOptions.bottomOptions.normType = opt.normType;
    opt.mgOptions.bottomOptions.verbosity = 0;
    // Use default convergence options.

    opt.lepticOptions = LevelLepticSolver::getDefaultOptions();
    opt.lepticOptions.normType = opt.normType;
    // opt.lepticOptions.verbosity = 0;

    opt.lepticOptions.horizOptions.normType = opt.normType;
    opt.lepticOptions.horizOptions.verbosity = 0;

    opt.lepticOptions.horizOptions.bottomOptions.normType = opt.normType;
    opt.lepticOptions.horizOptions.bottomOptions.verbosity = 0;

    return opt;
}


// -----------------------------------------------------------------------------
LevelHybridSolver::Options
LevelHybridSolver::getQuickAndDirtyOptions()
{
    Options opt = getDefaultOptions();

    opt.absTol         = 1.0e-300;
    opt.relTol         = 1.0e-2;
    opt.maxSolverSwaps = 1;

    opt.mgOptions = MGSolver<StateType>::getQuickAndDirtyOptions();
    opt.lepticOptions = LevelLepticSolver::getQuickAndDirtyOptions();

    return opt;
}


// -----------------------------------------------------------------------------
void
LevelHybridSolver::modifyOptionsExceptMaxDepth(const Options& a_options)
{
    if (!this->isDefined()) {
        MAYDAYERROR("This can only be called AFTER LevelHybridSolver is defined.");
    }

    const auto oldLepticVerb = m_options.lepticOptions.verbosity;
    m_options = a_options;
    m_options.lepticOptions.verbosity = oldLepticVerb;

    if (m_mgSolverPtr) {
        m_mgSolverPtr->modifyOptionsExceptMaxDepth(m_options.mgOptions);
    }
    if (m_lepticSolverPtr) {
        m_lepticSolverPtr->modifyOptionsExceptMaxDepth(m_options.lepticOptions);
    }
}



// ======================== Constructors / destructors =========================

// -----------------------------------------------------------------------------
LevelHybridSolver::LevelHybridSolver()
: m_isDefined(false)
, m_options()
, m_resNorms()
, m_solverStatus()
, m_mgOpPtr()
, m_solveMode(SolveMode::Undefined)
, m_lepticSolverPtr()
, m_mgSolverPtr()
{
}


// -----------------------------------------------------------------------------
LevelHybridSolver::~LevelHybridSolver()
{
    this->clear();
}


// -----------------------------------------------------------------------------
void
LevelHybridSolver::clear()
{
    m_options = Options();

    m_mgSolverPtr.reset();
    m_lepticSolverPtr.reset();
    m_solveMode = SolveMode::Undefined;
    m_mgOpPtr.reset();
    m_solverStatus.clear();

    m_resNorms.clear();
    m_isDefined = false;
}


// -----------------------------------------------------------------------------
void
LevelHybridSolver::define(
    std::shared_ptr<const MGOperator<StateType>> a_mgOpPtr,
    const Options&                               a_opts)
{
    if (m_isDefined) {
        this->clear();
    }

    m_options = a_opts;
    m_solverStatus.clear();

    CH_assert(a_mgOpPtr);
    m_mgOpPtr = a_mgOpPtr;

    m_solveMode = LevelHybridSolver::computeSolveMode(m_mgOpPtr, m_options);
    bool useLeptic = false;
    if (m_solveMode == SolveMode::Leptic || m_solveMode == SolveMode::Leptic_MG) {
        useLeptic = true;
    }
    bool useMG = false;
    if (m_solveMode == SolveMode::MG || m_solveMode == SolveMode::Leptic_MG) {
        useMG = true;
    }

    if (useLeptic && useMG) {
        // Silence the individual solvers. We will summarize at the end.
        m_options.lepticOptions.verbosity                            = 0;
        m_options.lepticOptions.horizOptions.verbosity               = 0;
        m_options.lepticOptions.horizOptions.bottomOptions.verbosity = 0;
        m_options.mgOptions.verbosity               = 0;
        m_options.mgOptions.bottomOptions.verbosity = 0;

        m_lepticSolverPtr.reset(new LevelLepticSolver);
        m_lepticSolverPtr->define(m_mgOpPtr, m_options.lepticOptions);
        m_options.lepticOptions = m_lepticSolverPtr->getOptions();

        m_mgSolverPtr.reset(new MGSolver<StateType>);
        m_mgSolverPtr->define(*m_mgOpPtr, m_options.mgOptions);
        m_options.mgOptions = m_mgSolverPtr->getOptions();

    } else if (useLeptic) {
        m_options.lepticOptions.verbosity = m_options.verbosity;

        m_lepticSolverPtr.reset(new LevelLepticSolver);
        m_lepticSolverPtr->define(m_mgOpPtr, m_options.lepticOptions);
        m_options.lepticOptions = m_lepticSolverPtr->getOptions();

    } else if (useMG) {
        m_options.mgOptions.verbosity = m_options.verbosity;

        m_mgSolverPtr.reset(new MGSolver<StateType>);
        m_mgSolverPtr->define(*m_mgOpPtr, m_options.mgOptions);
        m_options.mgOptions = m_mgSolverPtr->getOptions();

    } else {
        // ??? Just use MG
        useMG = true;
        m_mgSolverPtr.reset(new MGSolver<StateType>);
        m_mgSolverPtr->define(*m_mgOpPtr, m_options.mgOptions);
        m_options.mgOptions = m_mgSolverPtr->getOptions();
    }

    m_isDefined = true;
}


// -----------------------------------------------------------------------------
void
LevelHybridSolver::define(
    std::shared_ptr<LevelLepticSolver>           a_lepticSolverPtr,
    std::shared_ptr<MGSolver<StateType>>         a_mgSolverPtr,
    std::shared_ptr<const MGOperator<StateType>> a_mgOpPtr,
    const Options&                               a_opts)
{
    if (m_isDefined) {
        this->clear();
    }

    m_options = a_opts;
    m_solverStatus.clear();

    CH_assert(a_mgOpPtr);
    m_mgOpPtr = a_mgOpPtr;

    m_solveMode = LevelHybridSolver::computeSolveMode(m_mgOpPtr, m_options);
    bool useLeptic = false;
    if (m_solveMode == SolveMode::Leptic || m_solveMode == SolveMode::Leptic_MG) {
        useLeptic = true;
    }
    bool useMG = false;
    if (m_solveMode == SolveMode::MG || m_solveMode == SolveMode::Leptic_MG) {
        useMG = true;
    }

    if (useLeptic && useMG) {
        // Silence the individual solvers. We will summarize at the end.
        m_options.lepticOptions.verbosity                            = 0;
        m_options.lepticOptions.horizOptions.verbosity               = 0;
        m_options.lepticOptions.horizOptions.bottomOptions.verbosity = 0;
        m_options.mgOptions.verbosity               = 0;
        m_options.mgOptions.bottomOptions.verbosity = 0;

        CH_verify(a_lepticSolverPtr);
        m_lepticSolverPtr       = a_lepticSolverPtr;
        m_options.lepticOptions = a_lepticSolverPtr->getOptions();

        CH_verify(a_mgSolverPtr);
        CH_verify(a_mgSolverPtr->isDefined());
        m_mgSolverPtr       = a_mgSolverPtr;
        m_options.mgOptions = a_mgSolverPtr->getOptions();

    } else if (useLeptic) {
        CH_verify(a_lepticSolverPtr);
        m_lepticSolverPtr       = a_lepticSolverPtr;
        m_options.lepticOptions = a_lepticSolverPtr->getOptions();

    } else if (useMG) {
        CH_verify(a_mgSolverPtr);
        CH_verify(a_mgSolverPtr->isDefined());
        m_mgSolverPtr       = a_mgSolverPtr;
        m_options.mgOptions = a_mgSolverPtr->getOptions();

    } else {
        // ??? Just use MG
        CH_verify(a_mgSolverPtr);
        CH_verify(a_mgSolverPtr->isDefined());
        useMG               = true;
        m_mgSolverPtr       = a_mgSolverPtr;
        m_options.mgOptions = a_mgSolverPtr->getOptions();
    }

    m_isDefined = true;
}


// ================================== Solvers ==================================

// -----------------------------------------------------------------------------
SolverStatus
LevelHybridSolver::solve(StateType&       a_phi,
                         const StateType* a_crsePhiPtr,
                         const StateType& a_rhs,
                         const Real       a_time,
                         const bool       a_useHomogBCs,
                         const bool       a_setPhiToZero,
                         const Real       a_convergenceMetric) const
{
    // Sanity checks
    CH_verify(m_isDefined);
    CH_assert(a_phi.getBoxes() == m_mgOpPtr->getBoxes());
    CH_assert(a_rhs.getBoxes() == m_mgOpPtr->getBoxes());

    // Allocation
    StateType cor, res;
    m_mgOpPtr->create(cor, a_phi);
    m_mgOpPtr->create(res, a_rhs);

    // Solve
    if (a_setPhiToZero) m_mgOpPtr->setToZero(a_phi);
    m_mgOpPtr->residual(res, a_phi, a_crsePhiPtr, a_rhs, a_time, a_useHomogBCs, a_useHomogBCs);
    const auto solverStatus = this->solveResidualEq(cor, res, a_time, a_convergenceMetric);
    m_mgOpPtr->incr(a_phi, cor, 1.0);

    return solverStatus;
}


// -----------------------------------------------------------------------------
SolverStatus
LevelHybridSolver::solveResidualEq(StateType&       a_cor,
                                   const StateType& a_res,
                                   const Real       a_time,
                                   const Real       a_convergenceMetric) const
{
    // Sanity checks
    CH_verify(m_isDefined);
    CH_assert(a_cor.getBoxes() == m_mgOpPtr->getBoxes());
    CH_assert(a_res.getBoxes() == m_mgOpPtr->getBoxes());

    m_solverStatus.setSolverStatus(SolverStatus::UNDEFINED);
    pout() << Format::pushFlags << Format::scientific;

    StateType localRes;
    m_mgOpPtr->create(localRes, a_res);

    // Initialize correction
    m_mgOpPtr->setToZero(a_cor);

    // Initialize convergence metrics
    m_resNorms.clear();
    m_resNorms.push_back(m_mgOpPtr->norm(a_res, m_options.normType));
    if (a_convergenceMetric > 0.0) {
        m_resNorms[0] = a_convergenceMetric;
    }
    std::vector<std::string> solverTypes(1, "");
    m_solverStatus.setInitResNorm(m_resNorms.back());

    // Select solver.
    if (m_solveMode == SolveMode::Leptic) {
        constexpr bool homogBCs     = true;
        constexpr bool setCorToZero = false;
        m_solverStatus = m_lepticSolverPtr->solve(
            a_cor, nullptr, a_res, a_time, homogBCs, setCorToZero);

        m_resNorms.insert(m_resNorms.end(),
                          m_lepticSolverPtr->getResNorms().begin() + 1,
                          m_lepticSolverPtr->getResNorms().end());

        const std::vector<std::string> vl(m_lepticSolverPtr->getResNorms().size() - 1, "Leptic");
        solverTypes.insert(solverTypes.end(), vl.begin(), vl.end());

    } else if (m_solveMode == SolveMode::Leptic_MG) {
        for (int swaps = 0; swaps < m_options.maxSolverSwaps; ++swaps) {
            const Real preLepticResNorm = m_resNorms.back();

            // Leptic method as a preconditioner until it diverges.
            constexpr bool homogBCs     = true;
            constexpr bool setCorToZero = false;

            m_lepticSolverPtr->solve(a_cor, nullptr, a_res, a_time, homogBCs, setCorToZero);

            m_resNorms.insert(m_resNorms.end(),
                              m_lepticSolverPtr->getResNorms().begin() + 1,
                              m_lepticSolverPtr->getResNorms().end());
            const std::vector<std::string> vl(m_lepticSolverPtr->getResNorms().size() - 1, "Leptic");
            solverTypes.insert(solverTypes.end(), vl.begin(), vl.end());

            if (m_resNorms.back() <= m_options.absTol) {
                m_solverStatus.setSolverStatus(SolverStatus::CONVERGED);
                break;
            }
            if (m_resNorms.back() <= m_options.relTol * m_resNorms[0]) {
                m_solverStatus.setSolverStatus(SolverStatus::CONVERGED);
                break;
            }


            // MG to solve until it stalls or diverges.
            m_mgSolverPtr->vCycle_residualEq(a_cor, a_res, a_time, 0);
            m_mgOpPtr->residual(localRes, a_cor, nullptr, a_res, a_time, true, true);
            m_resNorms.push_back(m_mgOpPtr->norm(localRes, m_options.normType));
            solverTypes.push_back("V-cycle");

            if (m_resNorms.back() <= m_options.absTol) {
                m_solverStatus.setSolverStatus(SolverStatus::CONVERGED);
                break;
            }
            if (m_resNorms.back() <= m_options.relTol * m_resNorms[0]) {
                m_solverStatus.setSolverStatus(SolverStatus::CONVERGED);
                break;
            }
            if (m_resNorms.back() > preLepticResNorm) {
                m_solverStatus.setSolverStatus(SolverStatus::DIVERGED);
                break;
            }

            if (swaps == m_options.maxSolverSwaps - 1) {
                m_solverStatus.setSolverStatus(SolverStatus::MAXITERS);
                break;
            }
        }
    } else if (m_solveMode == SolveMode::MG) {
        // Fully leptic problem. Just use leptic solver.
        constexpr bool homogBCs     = true;
        constexpr bool setCorToZero = false;
        m_solverStatus = m_mgSolverPtr->solve(
            a_cor, nullptr, a_res, a_time, homogBCs, setCorToZero);

        // m_resNorms.insert(m_resNorms.end(),
        //                   m_mgSolverPtr->getResNorms().begin() + 1,
        //                   m_mgSolverPtr->getResNorms().end());
        // const std::vector<std::string> vl(m_lepticSolverPtr->getResNorms().size() - 1, "MG");
        // solverTypes.insert(solverTypes.end(), vl.begin(), vl.end());

        m_resNorms.push_back(m_solverStatus.getFinalResNorm());
        solverTypes.push_back("MG");
    } else {
        MAYDAYERROR("Solve mode not recognized.");
    }

    CH_assert(m_resNorms.size() == solverTypes.size());

    // Summarize convergence pattern.
    if (m_solveMode == SolveMode::Leptic_MG && m_options.verbosity >= 1) {
        for (size_t iter = 0; iter < m_resNorms.size(); ++iter) {
            pout() << "iter " << iter
                   << ": |rel res| = " << m_resNorms[iter] / m_resNorms[0]
                   << '\t' << solverTypes[iter] << '\n';
        }
        pout() << "LevelHybridSolver " << m_solverStatus.getExitMessage() << ".\n";
    }

    m_solverStatus.setFinalResNorm(m_resNorms.back());
    return m_solverStatus;
}


// // -----------------------------------------------------------------------------
// void
// LevelHybridSolver::setDefaultOptions()
// {
//     const auto& proj = ProblemContext::getInstance()->proj;

//     m_options.absTol         = proj.absTol;
//     m_options.relTol         = proj.relTol;
//     m_options.maxSolverSwaps = 10; // BUG: Hard-coded.
//     m_options.normType       = proj.normType;
//     m_options.verbosity      = proj.verbosity;
//     m_options.hang           = proj.hang;

//     m_options.mgOptions.absTol    = m_options.absTol;
//     m_options.mgOptions.relTol    = m_options.relTol;
//     m_options.mgOptions.hang      = m_options.hang;
//     m_options.mgOptions.maxDepth  = -1;
//     m_options.mgOptions.maxIters  = 10;
//     m_options.mgOptions.normType  = m_options.normType;
//     m_options.mgOptions.numCycles = -1;
//     m_options.mgOptions.verbosity = 0;
//     // Use default smoothing options.

//     m_options.mgOptions.bottomOptions.verbosity = 0;
//     // Use default convergence options.

//     m_options.lepticOptions.horizOptions.verbosity = 0;
//     m_options.lepticOptions.horizOptions.bottomOptions.verbosity = 0;
// }


// -----------------------------------------------------------------------------
LevelHybridSolver::SolveMode
LevelHybridSolver::computeSolveMode(
    std::shared_ptr<const MGOperator<StateType>>& a_mgOpPtr,
    const Options&                                a_options)
{
    const auto& grids  = a_mgOpPtr->getBoxes();
    const auto& dXi    = a_mgOpPtr->getDXi();
    // const Real  minDXi = std::min(std::min(dXi[0], dXi[1]), dXi[SpaceDim - 1]);
    const auto  Nx     = grids.physDomain().domainBox().size();
    const auto  L      = Nx * dXi;

    const Real lepticity  = std::min(dXi[0], dXi[SpaceDim - 2]) / L[SpaceDim - 1];
    // const Real anisotropy = (dXi / minDXi).product();

    // const Real coarsestAnisotropy = [&] () {
    //     SemicoarseningStrategy mgStrategy(L, a_options.verbosity);
    //     const auto refSchedule = mgStrategy.createMGRefSchedule(
    //         grids, a_mgOpPtr->minBoxSize(), a_options.mgOptions.maxDepth);

    //     IntVect totalRef = IntVect::Unit;
    //     for (const auto& r : refSchedule) {
    //         totalRef *= r;
    //     }
    //     const RealVect coarsestDXi = dXi * totalRef;
    //     const Real minCoarsestDXi =
    //         std::min(std::min(coarsestDXi[0], coarsestDXi[1]), coarsestDXi[SpaceDim - 1]);
    //     return (coarsestDXi / minCoarsestDXi).product();
    // } ();
    // const bool useLineRelax = (coarsestAnisotropy > 2.0);

    // Set main solve mode.
    SolveMode solveMode = SolveMode::Undefined;
    if (lepticity > 1.0) {
        solveMode = SolveMode::Leptic;
    } else if (lepticity > 0.2) { // eps in [1, ~30].
        solveMode = SolveMode::Leptic_MG;
    } else {
        solveMode = SolveMode::MG;
    }

    return solveMode;
}


}; // end namespace Elliptic
