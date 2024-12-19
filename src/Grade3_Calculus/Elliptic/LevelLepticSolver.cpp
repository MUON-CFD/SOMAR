#include "LevelLepticSolver.H"
#include "AMRMGOperator.H"
#include "AnisotropicRefinementTools.H"
#include "Debug.H"
#include "LayoutTools.H"
#include "LepticBoxTools.H"
#include "PoissonOpF_F.H"  // TODO: Don't reference PoissonOp.
#include "ProjectorParameters.H"
#include "SetValLevel.H"  // May not be needed
#include "Subspace.H"
#include "SubspaceF_F.H"

namespace Elliptic
{


// ======================== Option setters / getters ===========================

// -----------------------------------------------------------------------------
LevelLepticSolver::Options
LevelLepticSolver::getDefaultOptions()
{
    Options opt;
    const auto& proj = ProblemContext::getInstance()->proj;

    // Leptic solver parameters
    opt.absTol             = proj.absTol;
    opt.relTol             = proj.relTol;
    opt.maxOrder           = proj.maxIters;
    opt.normType           = proj.normType;
    opt.verbosity          = proj.verbosity;
    opt.hang               = proj.hang;
    opt.maxDivergingOrders = 2; // BUG: Hard-coded.

    // Horizontal solver parameters
    opt.horizOptions.absTol           = 1.0e-15;
    opt.horizOptions.relTol           = 1.0e-15;
    opt.horizOptions.numSmoothDown    = 4;
    opt.horizOptions.numSmoothUp      = 4;
    opt.horizOptions.numSmoothBottom  = 2;
    opt.horizOptions.numSmoothPrecond = 2;
    opt.horizOptions.prolongOrder     = 1;
    opt.horizOptions.prolongOrderFMG  = 3;
    opt.horizOptions.numSmoothUpFMG   = 0;
    opt.horizOptions.maxDepth         = -1;
    opt.horizOptions.numCycles        = 1;
    opt.horizOptions.maxIters         = 20;
    opt.horizOptions.hang             = 0.01;
    opt.horizOptions.normType         = opt.normType;
    opt.horizOptions.verbosity        = 0;

    // Horizontal bottom solver parameters
    opt.horizOptions.bottomOptions.absTol           = 1.0e-15;
    opt.horizOptions.bottomOptions.relTol           = 1.0e-15;
    opt.horizOptions.bottomOptions.small            = 1.0e-30;
    opt.horizOptions.bottomOptions.hang             = 0.01;
    opt.horizOptions.bottomOptions.maxIters         = 80;
    opt.horizOptions.bottomOptions.maxRestarts      = 5;
    opt.horizOptions.bottomOptions.normType         = opt.normType;
    opt.horizOptions.bottomOptions.verbosity        = 0;
    opt.horizOptions.bottomOptions.numSmoothPrecond = 2;

    return opt;
}


// -----------------------------------------------------------------------------
LevelLepticSolver::Options
LevelLepticSolver::getQuickAndDirtyOptions()
{
    Options opt = getDefaultOptions();

    // Forget this...the leptic method is already cheap!

    // // Leptic solver parameters
    // opt.absTol    = 1.0e-300;
    // opt.relTol    = 1.0e-2;
    // opt.maxOrder  = 2;

    // // Horizontal solver parameters
    // opt.horizOptions.verbosity = 0;

    // // Horizontal bottom solver parameters
    // opt.horizOptions.bottomOptions.verbosity = 0;

    return opt;
}


// -----------------------------------------------------------------------------
const LevelLepticSolver::Options&
LevelLepticSolver::getOptions() const
{
    return m_options;
}


// -----------------------------------------------------------------------------
void
LevelLepticSolver::modifyOptionsExceptMaxDepth(
    const LevelLepticSolver::Options& a_opt)
{
    if (!this->isDefined()) {
        MAYDAYERROR("This can only be called AFTER LevelLepticSolver is defined.");
    }

    const auto oldMaxDepth = m_options.horizOptions.maxDepth;
    m_options = a_opt;
    m_options.horizOptions.maxDepth = oldMaxDepth;
}


// ======================== Constructors / destructors =========================

// -----------------------------------------------------------------------------
LevelLepticSolver::LevelLepticSolver()
: LevelSolver()  // Initializes m_opPtr and m_solverStatus.
, m_isDefined(false)
, m_options()
, m_resNorms()
//
, m_L(D_DECL(quietNAN, quietNAN, quietNAN))
, m_dXi(D_DECL(quietNAN, quietNAN, quietNAN))
, m_dXiCrse(D_DECL(quietNAN, quietNAN, quietNAN))
, m_domain()
, m_origGrids()
//
, m_grids()
, m_JgupPtr()
, m_origToVertCopier()
, m_vertToOrigCopier()
, m_vertOpPtr()
//
, m_flatGrids()
, m_flatDI()
, m_flatDIComplement()
, m_shiftedFlatGrids()
//
, m_doHorizSolve(true)
, m_horizRemoveAvg(true)
, m_horizDomain()
, m_horizGrids()
, m_shiftedFlatToHorizCopier()
, m_horizToShiftedFlatCopier()
, m_horizMGOpPtr()
, m_horizSolverPtr()
{
}


// -----------------------------------------------------------------------------
LevelLepticSolver::~LevelLepticSolver()
{
    this->clear();
}


// -----------------------------------------------------------------------------
void
LevelLepticSolver::define(
    std::shared_ptr<const LevelOperator<StateType>> a_linearOpPtr)
{
    this->define(a_linearOpPtr, getDefaultOptions());
}


// -----------------------------------------------------------------------------
void
LevelLepticSolver::define(
    std::shared_ptr<const LevelOperator<StateType>> a_linearOpPtr,
    const Options&                                  a_opts)
{
    if (m_isDefined) {
        this->clear();
    }

    // This defines m_opPtr and m_solverStatus.
    m_options = a_opts;
    LevelSolver<StateType>::define(a_linearOpPtr);

    // Make sure we were handed an appropriate operator.
    CH_verify(a_linearOpPtr);
    auto lepticOpPtr =
        std::dynamic_pointer_cast<const LepticOperator>(a_linearOpPtr);
    if (!lepticOpPtr) {
        MAYDAYERROR(
            "LevelLepticSolver::define received a LinearOperator that cannot "
            "be cast into a LepticOperator.");
    }


    // Original grid stuff
    m_L                  = lepticOpPtr->getDomainLength();
    m_dXi                = lepticOpPtr->getDXi();
    m_dXiCrse            = lepticOpPtr->getDXiCrse();
    const auto& origJgup = lepticOpPtr->getFCJgup();
    m_origGrids          = origJgup.getBoxes();
    m_domain             = m_origGrids.physDomain();

    // Get domain and pack periodicity info into a C array.
    m_domain = m_origGrids.physDomain();
    bool isPeriodic[SpaceDim] = {D_DECL(
        m_domain.isPeriodic(0),
        m_domain.isPeriodic(1),
        m_domain.isPeriodic(2)
    )};
    CH_verify(!isPeriodic[SpaceDim - 1]);


    // Vertical grid stuff...
    {
        // Create boxes suitable for the vertical solver.
        Vector<Box> vertBoxArray;
        LepticBoxTools::createVerticalSolverGrids(
            vertBoxArray, m_origGrids.boxArray(), m_domain.domainBox());
        m_grids.defineAndLoadBalance(vertBoxArray, nullptr, m_domain);

        // Copy metric data.
        m_JgupPtr.reset(new LevelData<FluxBox>(m_grids, 1));
        CH_assert(m_JgupPtr);
        debugInitLevel(*m_JgupPtr);
        origJgup.copyTo(*m_JgupPtr);

        // This helps copy the rhs.
        m_origToVertCopier.define(m_origGrids, m_grids, m_domain, IntVect::Zero, false);
        m_vertToOrigCopier = m_origToVertCopier;
        m_vertToOrigCopier.reverse();

        // Create an operator that can calulate residuals and such.
        m_vertOpPtr = lepticOpPtr->createLevelOperator(m_grids);

    }  // end of vertical grid stuff.


    // Should be true. I only set this to false when debugging.
    m_doHorizSolve = true;

    // Set up flat and horizontal stuctures, if needed.
    if (m_doHorizSolve) {
        const Box& domBox = m_domain.domainBox();
        const int loIdx = domBox.smallEnd(SpaceDim - 1);

        // 1. We need to create a flat set of grids that are compatible with
        //    m_grids. The result will NOT be disjoint.
        m_flatGrids.deepCopy(m_grids);
        LepticBoxTools::FlattenTransform stampy(loIdx);
        m_flatGrids.transform(stampy);
        m_flatGrids.close();
        CH_assert(m_grids.compatible(m_flatGrids));

        // 2. Create an array of data indices that point to grids that
        //    vertically span the domain. This subset of m_flatGrids should be
        //    disjoint.
        m_flatDI.clear();
        m_flatDIComplement.clear();

        for (DataIterator dit(m_grids); dit.ok(); ++dit) {
            const DataIndex& di    = dit();
            const Box&       valid = m_grids[di];

            if (LepticBoxTools::vertSpanCheck(valid, domBox)) {
                m_flatDI.push_back(di);
            } else {
                m_flatDIComplement.push_back(di);
            }
        }

        // 3. Create flat grids at vertical index = 0 so that we can copy data
        //    to/from m_horizGrids without the need for a shifted copier.
        //    NOTE: This would not be needed if we could create a copier that
        //    supports shifting and periodicity.
        m_shiftedFlatGrids.deepCopy(m_flatGrids);
        LayoutTools::ShiftTransform shifty(-loIdx * s_vmask);
        m_shiftedFlatGrids.transform(shifty);
        m_shiftedFlatGrids.close();
        CH_assert(m_shiftedFlatGrids.compatible(m_flatGrids));

        // 4. We need to make a version of m_flatGrids that are suitable for
        //    horizontal solves. That is, we want to remove the unused grids and
        //    perform load balancing.
        const auto& horizBlockFactor =
            lepticOpPtr->getBlockFactor() * s_hmask + s_vmask;

        // 4a. Calculate horizontal domain.
        const Box flatDomBox = Subspace::flattenBox(domBox, SpaceDim-1);
        m_horizDomain.define(flatDomBox, isPeriodic);

        // 4b. Create the load-balanced horizontal grids.
        Vector<Box> horizBoxArray;
        LepticBoxTools::createHorizontalSolverGrids(horizBoxArray,
                                                    m_grids.boxArray(),
                                                    domBox,
                                                    horizBlockFactor);
        CH_assert(horizBoxArray.size() > 0);
        m_horizGrids.defineAndLoadBalance(horizBoxArray, nullptr, m_horizDomain);

        // 5. Create the flat <--> horiz copiers.
        m_shiftedFlatToHorizCopier.define(m_shiftedFlatGrids,
                                          m_horizGrids,
                                          m_horizDomain,
                                          IntVect::Zero,
                                          false);
        m_horizToShiftedFlatCopier = m_shiftedFlatToHorizCopier;
        m_horizToShiftedFlatCopier.reverse();

        // Create the horiz ops and solvers.
        m_horizMGOpPtr = lepticOpPtr->createHorizontalMGOperator(m_horizGrids);

        constexpr bool doVertCoarsening = false; // Can't coarsen a flat domain.
        const int      schedVerbosity   = m_options.horizOptions.verbosity;
        const IntVect  horizMinBoxSize  = m_horizMGOpPtr->minBoxSize();
        const int      horizMaxDepth    = m_options.horizOptions.maxDepth;

        const auto mgRefSchedule =
            HorizCoarseningStrategy(m_L, doVertCoarsening, schedVerbosity)
                .createMGRefSchedule(
                    m_horizGrids, horizMinBoxSize, horizMaxDepth);

        m_horizSolverPtr.reset(new MGSolver<LevelData<FArrayBox>>);
        m_horizSolverPtr->define(
            *m_horizMGOpPtr, m_options.horizOptions, mgRefSchedule);
        m_options.horizOptions = m_horizSolverPtr->getOptions();

        // 6. Do we need to remove the average from the horizontal solution?
        do {
            // I figure that we will need to remove the average if the
            // horizontal domain is completely spanned by m_horizGrids.
            // To check, we can just count the number of cells in each.
            // NOTE: This assumes Neumann BCs all around the domain!

            // Get the number of cells in the horizontal domain.
            const Box& horizDomBox = m_horizDomain.domainBox();
            if (!horizDomBox.numPtsOK()) {
                MAYDAYWARNING("LevelLepticSolver: horizDomBox.numPtsOK() failed");
                break;
            }
            const long domNumPts = horizDomBox.numPts();

            // Add up the number of cells in each grid.
            long numPts = 0;
            LayoutIterator lit = m_horizGrids.layoutIterator();
            for (lit.reset(); lit.ok(); ++lit) {
                numPts += m_horizGrids[lit].numPts();
            }

            // Compare.
            m_horizRemoveAvg = (numPts == domNumPts);
        } while(0);
    } // end flat and horiz grid stuff

    m_isDefined = true;
}


// -----------------------------------------------------------------------------
void
LevelLepticSolver::clear()
{
    // LinearSolver variables.
    m_solverStatus.clear();
    m_opPtr.reset();

    // The rest are LevelLepticSolver variables.
    m_horizSolverPtr.reset();
    m_horizMGOpPtr.reset();
    m_horizToShiftedFlatCopier.clear();
    m_shiftedFlatToHorizCopier.clear();
    m_horizGrids     = DisjointBoxLayout();
    m_horizDomain    = ProblemDomain();
    m_horizRemoveAvg = true;
    m_doHorizSolve   = true;

    m_shiftedFlatGrids = DisjointBoxLayout();
    m_flatDIComplement.clear();
    m_flatDI.clear();
    m_flatGrids = DisjointBoxLayout();

    m_vertOpPtr.reset();
    m_vertToOrigCopier.clear();
    m_origToVertCopier.clear();
    m_JgupPtr.reset();
    m_grids = DisjointBoxLayout();

    m_origGrids = DisjointBoxLayout();
    m_domain    = ProblemDomain();
    m_dXiCrse   = RealVect(D_DECL(quietNAN, quietNAN, quietNAN));
    m_dXi       = RealVect(D_DECL(quietNAN, quietNAN, quietNAN));
    m_L         = RealVect(D_DECL(quietNAN, quietNAN, quietNAN));

    m_resNorms.clear();
    m_options   = Options();
    m_isDefined = false;
}


// ================================== Solvers ==================================

// -----------------------------------------------------------------------------
SolverStatus
LevelLepticSolver::solve(StateType&       a_phi,
                         const StateType* a_crsePhiPtr,
                         const StateType& a_rhs,
                         const Real       a_time,
                         const bool       a_useHomogBCs,
                         const bool       a_setPhiToZero,
                         const Real       a_convergenceMetric) const
{
    SolverStatus& solverStatus = this->getSolverStatusRef();
    solverStatus.setSolverStatus(SolverStatus::UNDEFINED);
    pout() << Format::pushFlags << Format::scientific;

    // Sanity checks
    CH_verify(m_isDefined);
    CH_assert(a_phi.nComp() == 1);
    CH_assert(a_rhs.nComp() == 1);
    CH_assert(a_phi.getBoxes().physDomain() == m_domain);
    CH_assert(a_rhs.getBoxes().physDomain() == m_domain);
    CH_assert(a_phi.getBoxes() == m_origGrids);
    CH_assert(a_rhs.getBoxes() == m_origGrids);
    CH_assert(m_dXi[SpaceDim - 1] * m_domain.size(SpaceDim - 1) == m_L[SpaceDim - 1]);

    // Settings
    const int maxOrder           = m_options.maxOrder;
    bool      useExcess          = m_doHorizSolve;
    bool      useHorizPhi        = m_doHorizSolve;
    int       numDivergingOrders = 0;


    // Initialization ----------------------------------------------------------
    using std::shared_ptr;

    // Allocation.
    StateType             corTotal;                 // The cummulative solution.
    StateType             cor;                      // The current vertical solution.
    shared_ptr<StateType> rhsPtr(new StateType);    // The current vertical rhs.
    shared_ptr<StateType> tmpRhsPtr(new StateType); // Temp storage
    corTotal.define(m_grids, a_phi.nComp(), a_phi.ghostVect());
    cor.define(m_grids, a_phi.nComp(), a_phi.ghostVect());
    rhsPtr->define(m_grids, a_rhs.nComp(), a_rhs.ghostVect());
    tmpRhsPtr->define(m_grids, a_rhs.nComp(), a_rhs.ghostVect());

    BoxLayoutData<FArrayBox> excess;    // The excess function.
    BoxLayoutData<FArrayBox> flatRhs;   // Used to compute the horizontal rhs.
    StateType                horizPhi;  // The horizontal problem's solution.
    StateType                horizRhs;  // The horizontal problem's rhs.
    if (m_doHorizSolve) {
        excess  .define(m_flatGrids, 1);
        flatRhs .define(m_flatGrids, 1);
        horizPhi.define(m_horizGrids, 1, s_hmask);
        horizRhs.define(m_horizGrids, 1);
    }

    // Setup residual equation.
    {
        if (a_setPhiToZero) {
            m_opPtr->setToZero(a_phi);
        }

        LevelData<FArrayBox> res(m_origGrids, 1);
        m_opPtr->residual(res,
                          a_phi,
                          a_crsePhiPtr,
                          a_rhs,
                          a_time,
                          a_useHomogBCs,
                          a_useHomogBCs);
        res.copyTo(*rhsPtr, m_origToVertCopier);

        m_vertOpPtr->setToZero(corTotal);
    }

    // Initialize convergence metrics
    m_resNorms.clear();
    m_resNorms.reserve(maxOrder + 1);
    m_resNorms.push_back(m_vertOpPtr->norm(*rhsPtr, m_options.normType));
    if (m_options.verbosity >= 4) {
        pout() << "Absolute initial residual norm = "
               << std::scientific << m_resNorms[0] << '\n'
               << "Relative residual norm = 1.0\n";
    }

    // Initialize BC values.
    BoundaryData<Real> bdryData(m_grids, m_domain);
    CH_assert(!bdryData.isFlat());
    bdryData.setVal(0.0);
#ifndef NDEBUG
    // Check for consistency. Integral[rhs]-Integral[fluxes] should be zero.
    // This only needs to be done if we have Neum or periodic BCs all around.
    if (m_doHorizSolve && m_horizRemoveAvg && m_options.verbosity >= 4) {
        if (m_options.verbosity >= 4) {
            pout() << SpaceDim << "D consistency check (should be zero) = "
                   << bdryData.consistencyCheck(*rhsPtr, m_dXi)
                   << '\n';
        }
    }
#endif


    // Solve -------------------------------------------------------------------
    for (int order = 0; order <= maxOrder; ++order) {
        // Set up the BCs for all-Neum problem...
        // Do this even if we have abandoned the horizontal solves!
        if (m_doHorizSolve) {
            // Add previous excess to upper BC
            if (order >= 1 && useExcess) {
                bdryData.vertPlus(excess, 1.0, Side::Hi);
            }

            // Compute new excess or delete if no longer needed.
            if (useExcess) {
                this->computeVerticalExcess(excess, *rhsPtr, bdryData);

                // Report size of excess
                if (m_options.verbosity >= 5) {
                    const Real excessNorm = norm(
                        excess, excess.interval(), m_options.normType);
                    pout() << "Absolute excess norm = " << excessNorm << '\n';
                }

                // The algo does not require the excess function after O(1).
                if (order == 1) {
                    useExcess = false;
                }
            }

            // Remove excess from upper vertical BC.
            if (order == 0 && useExcess) {
                bdryData.vertPlus(excess, -1.0, Side::Hi);
            }
            // At this point, BCs and rhs should agree.
            // The vertical line solver will check for us.

        } // end if m_doHorizSolve

        // Solve the vertical problem and accumulate correction...
        this->verticalLineSolver(cor, *rhsPtr, bdryData);

        // Set up the horizontal problem if needed...
        if (useHorizPhi) {
            setValLevel(flatRhs, 0.0);

            // Subtract excess/H from horizontal rhs if needed.
            if (useExcess) {
                for (size_t idx = 0; idx < m_flatDI.size(); ++idx) {
                    const DataIndex& di = m_flatDI[idx];
                    flatRhs[di].plus(excess[di], -1.0 / m_L[SpaceDim - 1]);
                }
            }

            Real flatRhsNorm = norm(flatRhs, flatRhs.interval(), m_options.normType);
            if (m_options.verbosity >= 5) {
                pout() << "Flat relative rhs norm = "
                       << flatRhsNorm / m_resNorms[0] << '\n';
            }

            // Send flatRhs to m_horizGrids.
            setValLevel(horizRhs, 0.0);
            this->addFlatToHoriz(horizRhs, flatRhs);

            const Real horizRhsNorm =
                norm(horizRhs, horizRhs.interval(), m_options.normType);

            if (m_options.verbosity >= 5) {
                pout() << "Horizontal relative rhs norm = "
                       << horizRhsNorm / m_resNorms[0] << '\n';
            }

            // if (m_options.horizOptions.rhsRelTol * m_resNorms[0] > horizRhsNorm) {
            //     if (m_options.verbosity >= 4) {
            //         pout() << "\tAbandoning horizontal solves.\n";
            //     }
            //     useHorizPhi = false;
            // }
        } // end if useHorizPhi

        // Solve the horizontal problem and accumulate the correction...
        if (useHorizPhi) {
            const auto hStatus =
                this->horizontalSolver(horizPhi, horizRhs, a_time);

            if (hStatus.getSolverStatus() == SolverStatus::CONVERGED) {
                if (m_options.verbosity >= 3) {
                    pout() << "O(0h) relative |res| = "
                           << hStatus.getFinalResNorm() / hStatus.getInitResNorm()
                           << ' ' << hStatus.getExitMessage() << '\n';
                } else if (m_options.verbosity >= 2) {
                    pout() << "O(0h): " << hStatus.getExitMessage() << '\n';
                }
            } else {
                if (m_options.verbosity >= 3) {
                    pout() << "O(0h) relative |res| = "
                           << hStatus.getFinalResNorm() / hStatus.getInitResNorm()
                           << ' ' << hStatus.getExitMessage() << '\n';
                } else if (m_options.verbosity >= 1) {
                    pout() << "O(0h): " << hStatus.getExitMessage() << '\n';
                }
            }

            this->addHorizontalCorrection(cor, horizPhi);
            // At this point, cor holds phi_p^v + phi_p^h.
        }

        // Finalize order...
        {
            // Apply the correction.
            m_vertOpPtr->incr(corTotal, cor, 1.0);

            // Compute the new residual.
            {
                constexpr LevelData<FArrayBox>* crseVertPhiPtr = nullptr;
                constexpr bool                  homogPhysBCs   = true;
                constexpr bool                  homogCFBCs     = true;
                m_vertOpPtr->residual(*tmpRhsPtr,
                                      cor,
                                      crseVertPhiPtr,
                                      *rhsPtr,
                                      a_time,
                                      homogPhysBCs,
                                      homogCFBCs);
                m_resNorms.push_back(
                    m_vertOpPtr->norm(*tmpRhsPtr, m_options.normType));

                if (m_options.verbosity >= 2) {
                    pout() << "O(" << order << ") Relative |res| = "
                           << m_resNorms.back() / m_resNorms[0] << '\n';
                }
            }

            // Did we absolutely converge?
            if (m_resNorms.back() <= m_options.absTol) {
                if (m_options.verbosity >= 2) {
                    pout() << "Absolute convergence achieved.\n";
                }
                solverStatus.setSolverStatus(SolverStatus::CONVERGED);
                break;

            // Did we relatively converge?
            } else if (m_resNorms.back() <= m_options.relTol * m_resNorms[0]) {
                if (m_options.verbosity >= 1) {
                    pout() << "Relative convergence achieved.\n";
                }
                solverStatus.setSolverStatus(SolverStatus::CONVERGED);
                break;

            // Are we diverging?
            // Unlike AMRMGSolver, we keep the last diverging correction so
            // it can serve as a preconditioned initial guess for a standard
            // solver.
            } else  if (m_resNorms[order + 1] > m_resNorms[order]) {
                if (numDivergingOrders < m_options.maxDivergingOrders) {
                    ++numDivergingOrders;
                } else {
                    if (m_options.verbosity >= 1) {
                        pout() << "Diverging.\n";
                    }
                    solverStatus.setSolverStatus(SolverStatus::DIVERGED);
                    break;
                }

            // Are we hanging?
            } else if (m_resNorms[order + 1] > (1.0 - m_options.hang) * m_resNorms[order]) {
                if (m_options.verbosity >= 1) {
                    pout() << "Hanging.\n";
                }
                solverStatus.setSolverStatus(SolverStatus::HANG);
                break;

            // We are converging, but not finished yet.
            } else {
                numDivergingOrders = 0;

                if (order == maxOrder) {
                    solverStatus.setSolverStatus(SolverStatus::MAXITERS);
                }
            }

            // Prepare for next order.
            std::swap(rhsPtr, tmpRhsPtr);
            useHorizPhi = false;
        }
    } // end loop over orders


    // Finale ------------------------------------------------------------------

    // Send result back to original holder.
    corTotal.addTo(corTotal.interval(),
                   a_phi,
                   a_phi.interval(),
                   m_domain,
                   m_vertToOrigCopier);

    // Inform user of the result if we didn't already do so.
    if (m_options.verbosity == 3) {
        pout() << "Final relative residual norm = "
               << m_resNorms.back() / m_resNorms[0] << '\n';
    }

    pout() << Format::popFlags;
    m_solverStatus.setFinalResNorm(m_resNorms.back());
    return m_solverStatus;

}


// -----------------------------------------------------------------------------
// Computes the excess function, excess = hiNeumBC - loNeumBC - Integral[rhs].
// ERROR FOUND: If there are no grids on this proc, this function will crash.
// -----------------------------------------------------------------------------
void
LevelLepticSolver::computeVerticalExcess(
    BoxLayoutData<FArrayBox>&   a_excess,
    const LevelData<FArrayBox>& a_rhs,
    const BoundaryData<Real>&   a_bdryData) const
{
    // Sanity checks
    CH_assert(m_isDefined);

    CH_assert(m_doHorizSolve);
    CH_assert(a_excess.isDefined());

    CH_assert(a_excess.boxLayout() == m_flatGrids);
    CH_assert(a_rhs.getBoxes() == m_grids);

    CH_assert(a_excess.nComp() == 1);
    CH_assert(a_rhs.nComp() == 1);

    // Initialize phi with a bogus value
    debugInitLevel(a_excess);

    // Set the excess to zero in regions that don't need it.
    for (size_t idx = 0; idx < m_flatDIComplement.size(); ++idx) {
        const DataIndex& di = m_flatDIComplement[idx];
        a_excess[di].setVal(0.0);
    }

    // Loop over grids that require an excess
    for (size_t idx = 0; idx < m_flatDI.size(); ++idx) {
        const DataIndex& di        = m_flatDI[idx];
        FArrayBox&       excessFAB = a_excess[di];

        {  // Copy hi BC
            const FArrayBox& hiBCFAB = a_bdryData.getData(di, SpaceDim - 1, Side::Hi);
            CH_assert(!hiBCFAB.box().isEmpty());

            LayoutTools::VertShifter<Real> shifty(excessFAB, hiBCFAB);
            excessFAB.copy(hiBCFAB);
            shifty.restore();
        }
        {  // Subtract lo BC
            const FArrayBox& loBCFAB = a_bdryData.getData(di, SpaceDim - 1, Side::Lo);
            CH_assert(!loBCFAB.box().isEmpty());

            LayoutTools::VertShifter<Real> shifty(excessFAB, loBCFAB);
            excessFAB.plus(loBCFAB, -1.0);
            shifty.restore();
        }
        {  // Subtract integral
            const FArrayBox& rhsFAB  = a_rhs[di];
            const Box&       valid   = m_grids[di];
            const Real       dzScale = -m_dXi[SpaceDim - 1];
            const IntVect    shift   = excessFAB.box().smallEnd() * s_vmask;

            FORT_UNMAPPEDVERTINTEGRAL(
                CHF_FRA1_SHIFT(excessFAB, 0, shift),
                CHF_CONST_FRA1(rhsFAB, 0),
                CHF_BOX(valid),
                CHF_CONST_REAL(dzScale));
        }
    }
}


// -----------------------------------------------------------------------------
void
LevelLepticSolver::verticalLineSolver(
    LevelData<FArrayBox>&     a_vertPhi,
    LevelData<FArrayBox>&     a_vertRhs,
    const BoundaryData<Real>& a_bdryData) const
{
    // Sanity checks
    CH_assert(m_isDefined);
    CH_assert(m_JgupPtr);

    CH_assert(a_vertPhi.ghostVect()[SpaceDim - 1] >= 1);

    CH_assert(a_vertPhi.getBoxes() == m_grids);
    CH_assert(a_vertRhs.getBoxes() == m_grids);
    // CH_assert(a_vertBCTypes.boxLayout() == m_grids);
    CH_assert(m_JgupPtr->getBoxes().compatible(m_grids));

    CH_assert(a_vertPhi.nComp() == 1);
    CH_assert(a_vertRhs.nComp() == 1);

#ifndef NDEBUG
    // If we have Neum-Neum BCs, check that each Poisson problem is consistent.
    if (m_doHorizSolve && m_options.verbosity >= 5) {
        BoxLayoutData<FArrayBox> consistency(m_flatGrids, 1);
        this->computeVerticalExcess(consistency, a_vertRhs, a_bdryData);

        const Real consistencyNorm =
            norm(consistency, consistency.interval(), m_options.normType);

        pout() << "\tVertical problem consistency check (should be zero) = "
               << consistencyNorm << '\n';
    }
#endif

    // Initialize phi with a bogus value
    debugInitLevel(a_vertPhi);

    // Loop over grids and solve.
    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        FArrayBox&       phiFAB  = a_vertPhi[dit];
        FArrayBox&       rhsFAB  = a_vertRhs[dit];
        const FArrayBox& JgzzFAB = (*m_JgupPtr)[dit][SpaceDim - 1];
        const Box&       valid   = m_grids[dit];
        const int        Nz      = valid.size(SpaceDim - 1);
        const Real       dz      = m_dXi[SpaceDim - 1];

        const FArrayBox& upperBCFAB =
            a_bdryData.getData(dit(), SpaceDim - 1, Side::Hi);

        // Sanity checks
        CH_assert(rhsFAB.box().size(SpaceDim - 1) == Nz);
        CH_assert(phiFAB.box().size(SpaceDim - 1) == Nz + 2);
        CH_assert(phiFAB.box().contains(valid));
        CH_assert(rhsFAB.box().contains(valid));
        CH_assert(enclosedCells(JgzzFAB.box()).contains(valid));
        CH_assert(!m_domain.isPeriodic(SpaceDim - 1));

        BUG("Can only handle Neum-Neum BCs for now.");

        // Use the modified homogeneous Neumann tridiagonal solver.
        const int vertDir   = SpaceDim - 1;
        Box       bottomBox = adjCellLo(valid, vertDir, 1);
        bottomBox.shift(vertDir, 1);

        FORT_TRIDIAGPOISSONNN1DFAB(CHF_FRA(phiFAB),
                                   CHF_CONST_FRA(rhsFAB),
                                   CHF_CONST_FRA(upperBCFAB),
                                   CHF_CONST_FRA1(JgzzFAB, 0),
                                   CHF_BOX(bottomBox),
                                   CHF_CONST_INT(Nz),
                                   CHF_CONST_REAL(dz),
                                   CHF_CONST_INT(vertDir));
    }  // end loop over grids (dit)

    // writeLevelHDF5(a_vertRhs, 0.0, false);
    // writeLevelHDF5(a_vertPhi, 0.0, false);
}


// -----------------------------------------------------------------------------
// Computes the horizontal solutions.
// -----------------------------------------------------------------------------
SolverStatus
LevelLepticSolver::horizontalSolver(LevelData<FArrayBox>&       a_phi,
                                    const LevelData<FArrayBox>& a_rhs,
                                    const Real                  a_time) const
{
    // Sanity checks
    CH_assert(m_isDefined);
    CH_assert(m_doHorizSolve);
    CH_assert(a_phi.getBoxes() == m_horizGrids);
    CH_assert(a_rhs.getBoxes() == m_horizGrids);

#ifndef NDEBUG
    if (m_horizRemoveAvg && m_options.verbosity >= 5) {
        // Check for consistency. Integral[rhs]-Integral[fluxes] should be zero.
        BoundaryData<Real> horizBdryData(m_horizGrids, m_horizDomain);
        horizBdryData.setVal(0.0);
        Real consistency = horizBdryData.consistencyCheck(a_rhs, m_dXi);
        pout() << "\tFlat, " << SpaceDim - 1
               << "D consistency check (should be zero) = " << consistency
               << '\n';
    }
#endif

    constexpr LevelData<FArrayBox>* crsePhiPtr   = nullptr;
    constexpr bool                  useHomogBCs  = true;
    constexpr bool                  setPhiToZero = true;

    const SolverStatus status = m_horizSolverPtr->solve(
        a_phi, crsePhiPtr, a_rhs, a_time, useHomogBCs, setPhiToZero);

    if (m_horizRemoveAvg) {
        this->setZeroAvg(a_phi);
    }

    return status;
}


// -----------------------------------------------------------------------------
// Adds a horizontal correciton to the vertical solution.
// -----------------------------------------------------------------------------
void
LevelLepticSolver::addHorizontalCorrection(
    LevelData<FArrayBox>&       a_vertPhi,
    const LevelData<FArrayBox>& a_horizCor) const
{
    CH_assert(m_isDefined);
    CH_assert(a_vertPhi.nComp() == 1);
    CH_assert(a_horizCor.nComp() == 1);
    CH_assert(a_vertPhi.getBoxes() == m_grids);
    CH_assert(a_horizCor.getBoxes() == m_horizGrids);

    // Send correction to m_flatGrids.
    BoxLayoutData<FArrayBox> flatData(m_flatGrids, 1);
    debugInitLevel(flatData);
    this->copyHorizToFlat(flatData, a_horizCor);

    // Add the extrusion to the vertical solution.
    for (size_t idx = 0; idx < m_flatDI.size(); ++idx) {
        const DataIndex& di      = m_flatDI[idx];
        FArrayBox&       phiFAB  = a_vertPhi[di];
        const FArrayBox& flatFAB = flatData[di];
        const Box&       fullRegion =
            m_grids[di];  // This was phiFAB.box() and caused an invalid memory
                          // access error.

        FORT_ADDVERTICALEXTRUSION(CHF_FRA1(phiFAB, 0),
                                  CHF_CONST_FRA1(flatFAB, 0),
                                  CHF_BOX(fullRegion));
    }
}


// //
// -----------------------------------------------------------------------------
// // Collects each box's vertical BC type.
// // If all BCTypes are Neum or periodic, a_doHorizSolve will be set to true.
// //
// -----------------------------------------------------------------------------
// void
// LevelLepticSolver::gatherVerticalBCTypes(LevelData<FArrayBox>& a_loBCCoeffs,
//                                          LevelData<FArrayBox>& a_hiBCCoeffs,
//                                          bool& a_doHorizSolve, const
//                                          DisjointBoxLayout& a_grids, const
//                                          CFRegion&          a_CFRegion, const
//                                          LepticOperator&    a_lepticOp)
// {
//     // Gather grid info.
//     const auto& domain     = a_grids.physDomain();
//     const auto& domBox     = domain.domainBox();
//     const int   isPeriodic = domain.isPeriodic(SpaceDim - 1);
//     const auto  bcFunction = a_lepticOp.getBCFunction();

//     constexpr int maxSrcLayers = 2;
//     {
//         DisjointBoxLayout loBdryGrids;
//         loBdryGrids.deepCopy(a_grids);
//         loBdryGrids.adjCellSide(SpaceDim - 1, -maxSrcLayers, Side::Lo);
//         loBdryGrids.close();
//         a_loBCCoeffs.define(loBdryGrids, 1);

//         PhysBdryIter physBdryIter(a_grids, LoZSide);
//         for (physBdryIter.reset(); physBdryIter.ok(); ++physBdryIter) {
//             const DataIndex& di       = physBdryIter->di;
//             const Box&       valid    = physBdryIter->ccValidBox;
//             const Box&       ghostBox = physBdryIter->ccAdjGhostBox;
//             const int        isign    = physBdryIter->isign;
//             FArrayBox&       loFAB    = a_loBCCoeffs[di];

//         //     FArrayBox stateFAB(valid, 1);
//         //     stateFAB.setVal(0.0);

//         //     Box oneValidBox = ghostBox;
//         //     oneValidBox.shift(SpaceDim - 1, -isign);
//         //     stateFAB.setVal(1.0, oneValidBox, 0);

//         //     // Get important regions.
//         //     Box stateBdry, stateGhosts;
//         //     BCTools::getBoxesForApplyBC(
//         //         stateBdry, stateGhosts, stateFAB, ccValid, bdryDir, side);

//         // // Get physical coordinates at boundary.
//         // FArrayBox xFAB(stateBdry, SpaceDim);
//         // a_geoSrc.fill_physCoor(xFAB, a_dXi);

//         // // Get grid spacing at boundary, dx = dx/dXi * dXi.
//         // FArrayBox dxFAB(stateBdry, 1);
//         // a_geoSrc.fill_dxdXi(dxFAB, 0, bdryDir, a_dXi, a_dXi[bdryDir]);

//         }
//     }

//     {
//         DisjointBoxLayout hiBdryGrids;
//         hiBdryGrids.deepCopy(a_grids);
//         hiBdryGrids.adjCellSide(SpaceDim - 1, -maxSrcLayers, Side::Hi);
//         hiBdryGrids.close();
//         a_hiBCCoeffs.define(hiBdryGrids, 1);
//     }


//     // // Gather BCTypes from the BCDescriptors.
//     // // NOTE: Is is assumed that if we do not have Neum BCs, then the
//     fluxBCType is irrelevant.
//     // std::array<int, 2> physBCType;
//     // {
//     //     std::array<int, 2> ghostBCType;
//     //     ghostBCType[Side::Lo] =
//     a_bc.getGhostDescriptor()[SpaceDim-1][Side::Lo];
//     //     ghostBCType[Side::Hi] =
//     a_bc.getGhostDescriptor()[SpaceDim-1][Side::Hi];

//     //     std::array<int, 2> fluxBCType;
//     //     fluxBCType[Side::Lo] =
//     a_bc.getFluxDescriptor()[SpaceDim-1][Side::Lo];
//     //     fluxBCType[Side::Hi] =
//     a_bc.getFluxDescriptor()[SpaceDim-1][Side::Hi];

//     //     physBCType = ghostBCType;
//     //     if (ghostBCType[Side::Lo] == BCType::Neum || fluxBCType[Side::Lo]
//     == BCType::Neum) {
//     //         physBCType[Side::Lo] = BCType::Neum;
//     //     }
//     //     if (ghostBCType[Side::Hi] == BCType::Neum || fluxBCType[Side::Hi]
//     == BCType::Neum) {
//     //         physBCType[Side::Hi] = BCType::Neum;
//     //     }
//     // }

//     // for (SideIterator sit; sit.ok(); ++sit) {
//     //     const Side::LoHiSide iside = sit();
//     //     CH_assert(iside == 0 || iside == 1);
//     //     const int domSideIdx = domBox.sideEnd(iside)[SpaceDim-1];

//     //     for (dit.reset(); dit.ok(); ++dit) {
//     //         int& thisBCType = a_BCTypes[dit][iside];
//     //         const Box& valid = grids[dit];
//     //         const int validSideIdx = valid.sideEnd(iside)[SpaceDim-1];

//     //         // Initialize to error BCs.
//     //         thisBCType = BCType::Undefined;

//     //         if (validSideIdx == domSideIdx) {
//     //             // We are at a physical boundary...

//     //             if (isPeriodic) {
//     //                 thisBCType = BCType::Periodic;
//     //             } else {
//     //                 thisBCType = physBCType[iside];
//     //                 if (thisBCType != BCType::Neum) {
//     //                     a_doHorizSolve = false;
//     //                 }
//     //             }

//     //         } else {
//     //             // We are NOT at a physical boundary...

//     //             const Box ghostBox = adjCellBox(valid, SpaceDim-1, iside,
//     1);
//     //             const CFIVS& cfivs = (iside == Side::Lo?
//     //                                   a_CFRegion.loCFIVS(dit(),
//     SpaceDim-1):
//     //                                   a_CFRegion.hiCFIVS(dit(),
//     SpaceDim-1));

//     //             // If the cfivs is empty, we are not at a CF interface.
//     //             // We throw an error because this implies decomposition in
//     the vertical.
//     //             if (cfivs.isEmpty()) {
//     //                 MayDay::Error("Vertical grids are ill-formed");

//     //             } else  if (cfivs.isPacked()) {
//     //                 const Box& cfBox = cfivs.packedBox();
//     //                 if (cfBox.contains(ghostBox)) {
//     //                     thisBCType = BCType::CF;
//     //                     a_doHorizSolve = false;
//     //                 } else {
//     //                     MayDay::Error("Vertical grids are ill-formed");
//     //                 }

//     //             } else {
//     //                 const IntVectSet& ivs = cfivs.getIVS();
//     //                 if (ivs.contains(ghostBox)) {
//     //                     thisBCType = BCType::CF;
//     //                     a_doHorizSolve = false;
//     //                 } else {
//     //                     MayDay::Error("Vertical grids are ill-formed");
//     //                 }
//     //             }
//     //         } // end if region abuts physical boundary or not
//     //     } // end loop over grids (dit)
//     // } // end loop over sides (sit)

// //     // Search for Neum-Neum BCs, which demand horizontal solves.
// //     a_doHorizSolve = false;
// //     for (dit.reset(); dit.ok(); ++dit) {
// //         const int loBCType = a_BCTypes[dit][0];
// //         const int hiBCType = a_BCTypes[dit][1];

// //         const bool loNotDiri = (loBCType == BCType::Neum || loBCType ==
// BCType::Periodic);
// //         const bool hiNotDiri = (hiBCType == BCType::Neum || hiBCType ==
// BCType::Periodic);

// //         if (loNotDiri && hiNotDiri) {
// //             a_doHorizSolve = true;
// //             break;
// //         }
// //     }

// // #if CH_MPI
// //     // If any of the processors reported true, then set all to true.
// //     int globalDoSolve;
// //     int localDoSolve = (a_doHorizSolve? 1: 0);

// //     int ierr = MPI_Allreduce(&localDoSolve, &globalDoSolve, 1, MPI_INT,
// MPI_SUM, Chombo_MPI::comm);
// //     if (ierr != MPI_SUCCESS) {
// //         std::ostringstream errmsg;
// //         errmsg << "MPI_Allreduce failed. Error " << ierr << std::endl;
// //         MayDay::Error(errmsg.str().c_str());
// //     }

// //     a_doHorizSolve = (globalDoSolve > 0);
// // #endif
// }


// -----------------------------------------------------------------------------
void
LevelLepticSolver::addFlatToHoriz(LevelData<FArrayBox>&           a_horiz,
                                  const BoxLayoutData<FArrayBox>& a_flat) const
{
    // Sanity checks
    CH_assert(m_isDefined);
    CH_assert(a_horiz.getBoxes() == m_horizGrids);
    CH_assert(a_flat.boxLayout() == m_flatGrids);

    // Shift flat data vertically to coincide with a_horiz's boxes.
    BoxLayoutData<FArrayBox> shiftedFlat(m_shiftedFlatGrids, a_flat.nComp());
    for (DataIterator dit(m_flatGrids); dit.ok(); ++dit) {
        const IntVect shiftIV = m_domain.domainBox().smallEnd() * s_vmask;
        shiftedFlat[dit].shift(shiftIV);
        CH_assert(shiftedFlat[dit].box().smallEnd() == a_flat[dit].box().smallEnd());
        shiftedFlat[dit].copy(a_flat[dit]);
        shiftedFlat[dit].shift(-shiftIV);
    }

    // Add the data into a_horiz.
    shiftedFlat.addTo(shiftedFlat.interval(),
                      a_horiz,
                      a_horiz.interval(),
                      m_horizDomain,
                      m_shiftedFlatToHorizCopier);
}


// -----------------------------------------------------------------------------
void
LevelLepticSolver::copyHorizToFlat(BoxLayoutData<FArrayBox>&   a_flat,
                                   const LevelData<FArrayBox>& a_horiz) const
{
    // Sanity checks
    CH_assert(m_isDefined);
    CH_assert(a_horiz.getBoxes() == m_horizGrids);
    CH_assert(a_flat.boxLayout() == m_flatGrids);

    // Copy the data from a_horiz.
    BoxLayoutData<FArrayBox> shiftedFlat(m_shiftedFlatGrids, a_flat.nComp());
    a_horiz.copyTo(shiftedFlat, m_horizToShiftedFlatCopier);

    // Shift the new data vertically to coincide with a_flat's boxes.
    for (DataIterator dit(m_flatGrids); dit.ok(); ++dit) {
        const IntVect shiftIV = m_domain.domainBox().smallEnd() * s_vmask;
        shiftedFlat[dit].shift(shiftIV);
        CH_assert(shiftedFlat[dit].box().smallEnd() == a_flat[dit].box().smallEnd());
        a_flat[dit].copy(shiftedFlat[dit]);
        shiftedFlat[dit].shift(-shiftIV);
    }
}


// -----------------------------------------------------------------------------
// Static utility
// Brings average[a_phi] to zero.
// -----------------------------------------------------------------------------
void
LevelLepticSolver::setZeroAvg(LevelData<FArrayBox>& a_phi)
{
    const DisjointBoxLayout& grids = a_phi.getBoxes();
    DataIterator             dit   = a_phi.dataIterator();

    const int numproc  = numProc();
    const int thisproc = procID();
    const int srcproc  = uniqueProc(SerialTask::compute);

    Vector<Real> sumVect(numProc());
    Real         localSum  = 0.0;
    Real         globalSum = 0.0;

    Vector<long> volVect(numProc());
    long         localVol  = 0;
    long         globalVol = 0;

    Real globalAvg = 0.0;

    // Collect sums on each processor.
    for (dit.reset(); dit.ok(); ++dit) {
        const FArrayBox& phiFAB = a_phi[dit];
        const Box&       valid  = grids[dit];

        localSum += phiFAB.sum(valid, 0, 1);
        localVol += valid.numPts();
    }

    // Gather / broadcast total average.
    gather(sumVect, localSum, srcproc);
    gather(volVect, localVol, srcproc);
    if (thisproc == srcproc) {
        for (int idx = 0; idx < numproc; ++idx) {
            globalSum += sumVect[idx];
            globalVol += volVect[idx];
        }
        globalAvg = globalSum / Real(globalVol);
    }
    broadcast(globalAvg, srcproc);

    // Remove average from a_phi.
    for (dit.reset(); dit.ok(); ++dit) {
        a_phi[dit] -= globalAvg;
    }
}


};  // namespace Elliptic
