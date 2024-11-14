#include "LevelHybridSolver.H"

namespace Elliptic {

// -----------------------------------------------------------------------------
LevelHybridSolver::LevelHybridSolver()
: m_isDefined(false)
, m_options()
, m_resNorms()
, m_solverStatus()
, m_mgOpPtr()
, m_solveModeOrder()
, m_lepticSolverPtr()
, m_mgSolverPtr()
{
    this->setDefaultOptions();
}


// -----------------------------------------------------------------------------
LevelHybridSolver::~LevelHybridSolver()
{
    this->undefine();
}


// -----------------------------------------------------------------------------
void
LevelHybridSolver::undefine()
{
    m_mgSolverPtr.reset();
    m_lepticSolverPtr.reset();
    m_solveModeOrder.clear();
    m_mgOpPtr.reset();
    m_solverStatus.clear();

    m_resNorms.clear();
    m_isDefined = false;
}


// -----------------------------------------------------------------------------
void
LevelHybridSolver::define(std::shared_ptr<const MGOperator<StateType>> a_mgOpPtr)
{
    if (m_isDefined) {
        this->undefine();
    }

    m_solverStatus.clear();

    CH_assert(a_mgOpPtr);
    m_mgOpPtr = a_mgOpPtr;

    m_solveModeOrder = LevelHybridSolver::computeSolveModes(m_mgOpPtr, m_options);

    bool useLeptic = false;
    bool useMG = false;
    for (const auto& mode : m_solveModeOrder) {
        if (mode == SolveMode::Leptic || mode == SolveMode::Leptic_MG) {
            useLeptic = true;
        }
        if (mode == SolveMode::MG || mode == SolveMode::Leptic_MG) {
            useMG = true;
        }
    }

    if (useLeptic) {
        m_lepticSolverPtr.reset(new LevelLepticSolver);
        m_lepticSolverPtr->define(m_mgOpPtr, m_options.lepticOptions);
    }

    if (useMG) {
        m_mgSolverPtr.reset(new MGSolver<StateType>);
        m_mgSolverPtr->define(*m_mgOpPtr, m_options.mgOptions);
    }

    m_isDefined = true;
}


// -----------------------------------------------------------------------------
void
LevelHybridSolver::define(
    std::shared_ptr<const MGOperator<StateType>> a_mgOpPtr,
    const Options&                               a_opts)
{
    if (m_isDefined) {
        this->undefine();
    }

    m_options = a_opts;
    m_solverStatus.clear();

    CH_assert(a_mgOpPtr);
    m_mgOpPtr = a_mgOpPtr;

    m_solveModeOrder = LevelHybridSolver::computeSolveModes(m_mgOpPtr, m_options);

    bool useLeptic = false;
    bool useMG = false;
    for (const auto& mode : m_solveModeOrder) {
        if (mode == SolveMode::Leptic || mode == SolveMode::Leptic_MG) {
            useLeptic = true;
        }
        if (mode == SolveMode::MG || mode == SolveMode::Leptic_MG) {
            useMG = true;
        }
    }

    if (useLeptic) {
        m_lepticSolverPtr.reset(new LevelLepticSolver);
        m_lepticSolverPtr->define(m_mgOpPtr, m_options.lepticOptions);
        m_options.lepticOptions = m_lepticSolverPtr->getOptions();
    }

    if (useMG) {
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
        this->undefine();
    }

    m_options = a_opts;
    m_solverStatus.clear();

    CH_assert(a_mgOpPtr);
    m_mgOpPtr = a_mgOpPtr;

    m_solveModeOrder = LevelHybridSolver::computeSolveModes(m_mgOpPtr, m_options);

    bool useLeptic = false;
    bool useMG = false;
    for (const auto& mode : m_solveModeOrder) {
        if (mode == SolveMode::Leptic || mode == SolveMode::Leptic_MG) {
            useLeptic = true;
        }
        if (mode == SolveMode::MG || mode == SolveMode::Leptic_MG) {
            useMG = true;
        }
    }

    if (useLeptic) {
        CH_assert(a_lepticSolverPtr);
        m_lepticSolverPtr = a_lepticSolverPtr;
        m_options.lepticOptions = a_lepticSolverPtr->getOptions();
    }

    if (useMG) {
        CH_assert(a_mgSolverPtr);
        CH_assert(a_mgSolverPtr->isDefined());
        m_mgSolverPtr = a_mgSolverPtr;
        m_options.mgOptions = a_mgSolverPtr->getOptions();
    }

    m_isDefined = true;
}


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
    const auto solveMode = m_solveModeOrder[0];
    if (solveMode == SolveMode::Leptic) {
        // Fully leptic problem. Just use leptic solver.
        constexpr bool homogBCs     = true;
        constexpr bool setCorToZero = false;
        m_solverStatus = m_lepticSolverPtr->solve(
            a_cor, nullptr, a_res, a_time, homogBCs, setCorToZero);

        m_resNorms.insert(m_resNorms.end(),
                          m_lepticSolverPtr->getResNorms().begin() + 1,
                          m_lepticSolverPtr->getResNorms().end());

        const std::vector<std::string> vl(m_lepticSolverPtr->getResNorms().size() - 1, "Leptic");
        solverTypes.insert(solverTypes.end(), vl.begin(), vl.end());

    } else if (solveMode == SolveMode::Leptic_MG) {
        for (int swaps = 0; swaps < m_options.maxSolverSwaps; ++swaps) {
            const Real preLepticResNorm = m_resNorms.back();

            // Leptic method as a preconditioner until it diverges.
            constexpr bool homogBCs     = true;
            constexpr bool setCorToZero = false;

            const auto oldOptions = m_lepticSolverPtr->getOptions();
            auto newOptions = oldOptions;
            newOptions.verbosity = 0;
            m_lepticSolverPtr->setOptions(newOptions);
            m_lepticSolverPtr->solve(a_cor, nullptr, a_res, a_time, homogBCs, setCorToZero);
            m_lepticSolverPtr->setOptions(oldOptions);

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
    } else if (solveMode == SolveMode::MG) {
        // Fully leptic problem. Just use leptic solver.
        constexpr bool homogBCs     = true;
        constexpr bool setCorToZero = false;
        const auto oldOptions = m_mgSolverPtr->getOptions();
        auto newOptions = oldOptions;
        newOptions.verbosity = 4;
        m_mgSolverPtr->setOptions(newOptions);
        m_solverStatus = m_mgSolverPtr->solve(
            a_cor, nullptr, a_res, a_time, homogBCs, setCorToZero);
        m_mgSolverPtr->setOptions(oldOptions);

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

    if (m_options.verbosity >= 1) {
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


// -----------------------------------------------------------------------------
void
LevelHybridSolver::setDefaultOptions()
{
    const auto& proj = ProblemContext::getInstance()->proj;

    m_options.absTol         = proj.absTol;
    m_options.relTol         = proj.relTol;
    m_options.maxSolverSwaps = 10; // BUG: Hard-coded.
    m_options.normType       = proj.normType;
    m_options.verbosity      = proj.verbosity;
    m_options.hang           = proj.hang;

    m_options.mgOptions.absTol    = m_options.absTol;
    m_options.mgOptions.relTol    = m_options.relTol;
    m_options.mgOptions.hang      = m_options.hang;
    m_options.mgOptions.maxDepth  = -1;
    m_options.mgOptions.maxIters  = 10;
    m_options.mgOptions.normType  = m_options.normType;
    m_options.mgOptions.numCycles = -1;
    m_options.mgOptions.verbosity = 0;
    // Use default smoothing options.

    m_options.mgOptions.bottomOptions.verbosity = 0;
    // Use default convergence options.
}


// -----------------------------------------------------------------------------
std::vector<LevelHybridSolver::SolveMode>
LevelHybridSolver::computeSolveModes(
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
    constexpr size_t numSolveModeAttempts = 1;
    std::vector<SolveMode> solveMode;
    solveMode.reserve(numSolveModeAttempts);
    if (lepticity > 1.0) {
        solveMode.push_back(SolveMode::Leptic);
    } else if (lepticity > 0.2) { // eps in [1, ~30].
        solveMode.push_back(SolveMode::Leptic_MG);
    } else {
        solveMode.push_back(SolveMode::MG);
    }

    return solveMode;
}


}; // end namespace Elliptic
