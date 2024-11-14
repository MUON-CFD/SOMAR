/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2018
 *    Jefferson (Philadelphia University + Thomas Jefferson University) and
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
#include "AMRNSLevel.H"
#include "AMRNSLevelF_F.H"
#include "CFInterp.H"
#include "ProblemContext.H"
#include "SetValLevel.H"

#include "BoxIterator.H"
#include "Debug.H"
#include "IO.H"
#include "LayoutTools.H"
#include "MGSolver.H"
#include "lapack.H"
// #include "FETools.H" // TODO: IB

#include <chrono> // for high_resolution_clock
#include "AnisotropicRefinementTools.H"
#include "Subspace.H"
#include "BCTools.H"
#include "PoissonOp.H"

#include "FluxRegisterFace.H"


// -----------------------------------------------------------------------------
// Weak constructor. Leaves object unusable.
// Partial constructor. Rest of construction happens when AnisotropicAMR
// calls define.
// -----------------------------------------------------------------------------
AMRNSLevel::AMRNSLevel()
: m_isFactoryDefined(false),
  m_isActivated(false),
  m_levelOpsValid(false),
  m_levelSolversValid(false),
  m_amrOpsValid(false),
  m_amrSolversValid(false),
  m_geoSrcPtr(nullptr),
  m_levGeoPtr(nullptr),
  m_finiteDiffPtr(nullptr),
  m_velPtr(nullptr),
  m_pPtr(nullptr),
  m_qPtr(nullptr),
  m_statePtr(nullptr),
  m_cfInterpPtr(nullptr),
  m_velFluxRegPtr(nullptr),
  m_qFluxRegPtr(nullptr),
  m_c1(quietNAN),
  m_parkPtr(nullptr)
{
    // Set by base constructor:
    // m_coarser_level_ptr = nullptr;
    // m_finer_level_ptr = nullptr;
    // m_level = 0;
    // m_time = 0;
    // m_dt = 0;
    // m_initial_dt_multiplier = 0.1;
    // m_isDefined = false;
}


// -----------------------------------------------------------------------------
// Off limits. This will throw an error.
// -----------------------------------------------------------------------------
void
AMRNSLevel::define(AnisotropicAMRLevel* /*a_coarse_level_ptr*/,
                   const Box&           /*a_problem_domain*/,
                   int                  /*a_level*/,
                   const IntVect&       /*a_ref_ratio*/)
{
    MayDay::Error(
        "AMRNSLevel::define that takes a Box instead of a "
        "ProblemDomain is off limits");
}


// -----------------------------------------------------------------------------
// Full virtual constructor
// The factory will create the new AMRNSLevel object, then
// AnisotropicAMR will call this define function. So, you should expect
// all NavierStokes-specific / level-independent parameters to be set.
// -----------------------------------------------------------------------------
void
AMRNSLevel::define(AnisotropicAMRLevel* a_coarse_level_ptr,
                   const ProblemDomain& a_problem_domain,
                   int                  a_level,
                   const IntVect&       a_ref_ratio)
{
    // Only the factory can create an AMRNSLevel object.
    // If someone else is trying to do so, catch the bug here.
    CH_assert(m_isFactoryDefined);

    // Define the generic stuff. This will set:
    // m_coarser_level_ptr = a_coarser_level_ptr;
    // m_problem_domain = a_problem_domain;
    // m_level = a_level;
    // m_fineRefRatio = a_ref_ratio;
    // m_finer_level_ptr = nullptr;
    // m_isDefined = true;
    AnisotropicAMRLevel::define(
        a_coarse_level_ptr, a_problem_domain, a_level, a_ref_ratio);

    // Initialize State. This does not allocate or attach, just sets
    // some static info.
    m_statePtr = new State(this->numScalars());
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
AMRNSLevel::~AMRNSLevel()
{
    this->deactivateLevel();  // This detaches m_statePtr from m_*Ptr.
    delete m_statePtr;        // deactivateLevel() does not do this.
}


// -----------------------------------------------------------------------------
// Initializes this level to have the specified domain a_new_grids.
// -----------------------------------------------------------------------------
void
AMRNSLevel::initialGrid(const Vector<Box>& a_new_grids)
{
    BEGIN_FLOWCHART();

    // Chack if level is inactive.
    if (a_new_grids.size() == 0) {
        this->deactivateLevel();
    } else {
        this->activateLevel(a_new_grids);
    }
}


// -----------------------------------------------------------------------------
// Performs operations required after the grid has been defined but before
// data initialization. This will also be called after readCheckpointLevel
// during a restart procedure with argument a_restart set to true.
// Called from top -> bottom.
// -----------------------------------------------------------------------------
void
AMRNSLevel::postInitialGrid(const bool a_restart)
{
    BEGIN_FLOWCHART();

    // How much initialization can be moved into here and postRegrid()?

    if (a_restart) {
        // Remember, postInitialize will NOT be called after a restart,
        // so don't rely on it to validate your ops and solvers.

        // This must be done before anything else, but it also needs to be called
        // on the coarsest level, so...
        if (this->isFinestLevel()) {
            this->coarsestNSPtr()->validateOpsAndSolvers();
        }
    }
}


// -----------------------------------------------------------------------------
// Set initial conditions.
// -----------------------------------------------------------------------------
void
AMRNSLevel::initialData()
{
    BEGIN_FLOWCHART();

    if (this->isEmpty()) return;

    // Set everything to zero except T and S, which are set to the background
    // profile. The user can add perturbations
    // or completely overwrite these these.
    setValLevel(m_statePtr->vel, 0.0);
    setValLevel(m_statePtr->p, 0.0);
    setValLevel(m_statePtr->q, 0.0);

    DataIterator dit = m_statePtr->dit;
    for (dit.reset(); dit.ok(); ++dit) {
        Subspace::horizontalExtrusion(m_statePtr->T[dit], 0, *m_TbarPtr, 0, 1);
        Subspace::horizontalExtrusion(m_statePtr->S[dit], 0, *m_SbarPtr, 0, 1);
    }

    // Call the workhorse.
    this->setICs(*m_statePtr);

    LevelData<FArrayBox> eddyNu;
    aliasLevelData(eddyNu, &(m_statePtr->q), m_statePtr->eddyNuInterval);
    this->computeEddyNu(eddyNu, m_statePtr->vel, 1.0);
}


// -----------------------------------------------------------------------------
// Tag cells during initialization.
// -----------------------------------------------------------------------------
void
AMRNSLevel::tagCellsInit(IntVectSet& a_tags)
{
    BEGIN_FLOWCHART();
    this->tagCells(a_tags);
}


// -----------------------------------------------------------------------------
// Perform initial projection and set up the initial pressure.
// This must be done over the entire hierarchy at once, so we wait until
// this is called on level 0 and do everything there.
// Called from top -> bottom.
// -----------------------------------------------------------------------------
void
AMRNSLevel::postInitialize()
{
    BEGIN_FLOWCHART();

    if (this->isEmpty()) return;

    checkForValidNAN(m_statePtr->vel);
    debugCheckValidFaceOverlap(m_statePtr->vel);
    checkForValidNAN(m_statePtr->p);
    checkForValidNAN(m_statePtr->q);

    const ProblemContext* ctx = ProblemContext::getInstance();

    // This must be done before anything else, but it also needs to be called
    // on the coarsest level, so...
    if (this->isFinestLevel()) {
        this->coarsestNSPtr()->validateOpsAndSolvers();
    }

    // Average down and project.
    if (m_level == 0) {
        // Set all the BCs on every level. I do this before averaging down just
        // in case the interpolation schemes someday require ghosts.
        AMRNSLevel* levPtr = this;
        while (levPtr) {
            levPtr->setBC(*levPtr->m_statePtr, 0.0);
            levPtr = levPtr->fineNSPtr();
        };

        // Make sure coarse data is average of fine data.
        this->averageDownToThis();

        // Now that data is averaged down (and altered on most levels),
        // go through and reset all the BCs on every level.
        this->setBCsDownToThis();

        // Initial projection
        if (ctx->proj.doInitProj) {
            pout() << "\nInitialization projection:" << endl;
            this->projectDownToThis();
            pout() << endl;
        }

        // Make sure coarse data is average of fine data.
        this->averageDownToThis();

        // Now that data is averaged down (and altered on most levels),
        // go through and reset all the BCs on every level.
        this->setBCsDownToThis();
    }

    checkForValidNAN(m_statePtr->vel);
    debugCheckValidFaceOverlap(m_statePtr->vel);
    checkForValidNAN(m_statePtr->p);
    checkForValidNAN(m_statePtr->q);

    // Initial pressure estimate
    if (m_level == 0 && ctx->rhs.computeInitPressure && m_levelProjSolverPtr) {
        TODONOTE("Limit dt as if we are not using implicit diffusion.");
        this->initAMRPressure(0.1);
    }

    checkForValidNAN(m_statePtr->vel);
    debugCheckValidFaceOverlap(m_statePtr->vel);
    checkForValidNAN(m_statePtr->p);
    checkForValidNAN(m_statePtr->q);

    // Write initial diagnostic info to terminal.
    if (m_level == 0) {
        this->printDiagnostics(0, true);
    }}


// -----------------------------------------------------------------------------
// Returns maximum stable time step for this level with initial data.
// -----------------------------------------------------------------------------
Real
AMRNSLevel::computeInitialDt()
{
    BEGIN_FLOWCHART();
    return this->computeDt() * m_initial_dt_multiplier;
}


// -----------------------------------------------------------------------------
// Frees all memory on this level and clears grids. Called during
// construction, destruction, and regridding when a level vanishes. Note
// that this also deletes the solvers. To completely re-activate the level,
// you need to call activateLevel() AND allocAndDefine*Solvers().
//
// Since it doesn't make sense for level l to be deactivated and level
// l+1 to be active, this function will ensure all levels above this level
// are properly deactivated as well.
// -----------------------------------------------------------------------------
void
AMRNSLevel::deactivateLevel()
{
    // Do this before printing the flowchart stuff so that the user sees
    // the order that levels are actually being deactivated, not the order
    // that the functions are called.
    {
        // Do we have work to do?
        if (!m_isActivated) return;

        // Level's must be deactivated from top to bottom.
        AMRNSLevel* finePtr = this->fineNSPtr();
        if (finePtr) {
            finePtr->deactivateLevel();
        }
    }

    BEGIN_FLOWCHART();
    if (s_verbosity >= 3) {
        pout() << "Deactivating level " << m_level << "." << endl;
    }

    m_isActivated = false;

    m_IBPtr.reset(nullptr); // TODO: IB

    // This level gets completely invalidated.
    m_levelOpsValid     = false;
    m_levelSolversValid = false;
    m_amrOpsValid       = false;
    m_amrSolversValid   = false;

    if (m_level > 0) {
        // The next coarser level's AMR ops and solvers depend on this level
        // for refluxing, so they get invalidated.
        AMRNSLevel* levPtr = this->crseNSPtr();
        CH_assert(levPtr);
        levPtr->m_amrOpsValid       = false;
        levPtr->m_amrSolversValid   = false;

        // The AMR solvers on all levels become invalidated.
        levPtr = levPtr->crseNSPtr();
        while (levPtr) {
            levPtr->m_amrSolversValid = false;
            levPtr = levPtr->crseNSPtr();
        }
    }

    // Notes and warnings:
    //
    // 1. No need to detach from coarser level. The fineNSPtr* functions
    // will return nullptr when it sees the level is empty.
    //
    // 2. Do not deallocate m_statePtr until destructor.
    //
    // 3. m_geoSrcPtr is maintained by the factory. Do not delete it.

    delete m_parkPtr;
    m_parkPtr = nullptr;

    m_statePtr->detachFromQ();
    delete m_velPtr;
    delete m_pPtr;
    delete m_qPtr;
    m_velPtr = nullptr;
    m_pPtr   = nullptr;
    m_qPtr   = nullptr;

    delete m_cfInterpPtr;
    m_cfInterpPtr = nullptr;

    // Each level takes care of the coarser level's flux register.
    if (m_level > 0) {
        AMRNSLevel* crsePtr = this->crseNSPtr();

        delete crsePtr->m_velFluxRegPtr;
        crsePtr->m_velFluxRegPtr = nullptr;

        delete crsePtr->m_qFluxRegPtr;
        crsePtr->m_qFluxRegPtr = nullptr;
    }

    m_amrProjSolverPtr.reset();
    m_levelProjSolverPtr.reset();
    m_projOpPtr.reset();

    delete m_finiteDiffPtr;
    m_finiteDiffPtr = nullptr;

    for (int d = 0; d < SpaceDim; ++d) {
        m_singlePhysBdryIter[d][0].clear();
        m_singlePhysBdryIter[d][1].clear();
    }
    m_allPhysBdryIter.clear();

    delete m_levGeoPtr;
    m_levGeoPtr = nullptr;

    m_level_grids.clear();
}


// -----------------------------------------------------------------------------
// Allocates and defines structures when this level comes into existence.
// This does NOT define the solvers. Do that after level activation
// by calling validateOpsAndSolvers().
//
// The solvers are initialized separately because it may be a very
// expensive operation, and activatLevel may be called many times during
// construction.
// -----------------------------------------------------------------------------
void
AMRNSLevel::activateLevel(const Vector<Box>& a_new_grids,
                          const Vector<int>& a_procMap)
{
    BEGIN_FLOWCHART();

    CH_assert(a_procMap.size() == 0 || a_procMap.size() == a_new_grids.size());

    // Prevent a double new.
    this->deactivateLevel();

    m_isActivated = true;

    if (s_verbosity >= 3) {
        pout() << "Activating level " << m_level << "." << endl;
    }

    // Copy this level's new grids, sort them, then balance the load.
    // From Wikipedia: Morton ordering maps multidimensional data to one
    // dimension while preserving locality of the data points.

    m_level_grids = a_new_grids;

    mortonOrdering(m_level_grids);

    const DisjointBoxLayout grids =
        a_procMap.size() == 0
            ? this->loadBalance(m_level_grids)
            : DisjointBoxLayout(a_new_grids, a_procMap, m_problem_domain);


    // Grids are now well-defined on this level.

    // Set up the geometry.
    if (m_levGeoPtr == nullptr) {
        const ProblemContext* ctx = ProblemContext::getInstance();

        AMRNSLevel* crsePtr = this->crseNSPtr();

        LevelGeometry* crseLevGeoPtr =
            (crsePtr ? crsePtr->m_levGeoPtr : nullptr);

        m_levGeoPtr = new LevelGeometry(
            m_problem_domain, ctx->base.L, crseLevGeoPtr, m_geoSrcPtr);
    }
    m_levGeoPtr->createMetricCache(grids);

    // Set up boundary iterators.
    m_allPhysBdryIter.define(grids);
    for (int d = 0; d < SpaceDim; ++d) {
        SideArray sideArray = NoSides;
        sideArray[d][0] = 1;
        m_singlePhysBdryIter[d][0].define(grids, sideArray);

        sideArray[d][0] = 0;
        sideArray[d][1] = 1;
        m_singlePhysBdryIter[d][1].define(grids, sideArray);
    }

    // Set up finite difference ops.
    delete m_finiteDiffPtr;
    m_finiteDiffPtr = new FiniteDiff(*m_levGeoPtr);

    // Stratification, structure functions, etc...
    this->setStratificationMembers();

    // Allocate and define BCFunctions.
    {
        using std::array;
        using std::shared_ptr;
        using namespace BCTools;

        array<array<shared_ptr<BCFunction>, 2>, SpaceDim> velBCPtrs;
        array<array<shared_ptr<BCFunction>, 2>, SpaceDim> pBCPtrs;
        array<array<shared_ptr<BCFunction>, 2>, SpaceDim> TBCPtrs;
        array<array<shared_ptr<BCFunction>, 2>, SpaceDim> SBCPtrs;
        array<array<shared_ptr<BCFunction>, 2>, SpaceDim> scalarsBCPtrs;

        for (int bdryDir = 0; bdryDir < SpaceDim; ++bdryDir) {
            for (SideIterator sit; sit.ok(); ++sit) {
                velBCPtrs[bdryDir][int(sit())].reset(
                    this->createVelPhysBC(bdryDir, sit()));
                pBCPtrs[bdryDir][int(sit())].reset(
                    this->createPressurePhysBC(bdryDir, sit()));
                TBCPtrs[bdryDir][int(sit())].reset(
                    this->createTemperaturePhysBC(bdryDir, sit()));
                SBCPtrs[bdryDir][int(sit())].reset(
                    this->createSalinityPhysBC(bdryDir, sit()));
                scalarsBCPtrs[bdryDir][int(sit())].reset(
                    this->createScalarsPhysBC(bdryDir, sit()));
            }
        }
        m_velBCPtr.reset(new UberBCFunction(velBCPtrs));
        m_pBCPtr.reset(new UberBCFunction(pBCPtrs));
        m_TBCPtr.reset(new UberBCFunction(TBCPtrs));
        m_SBCPtr.reset(new UberBCFunction(SBCPtrs));
        m_scalarsBCPtr.reset(new UberBCFunction(scalarsBCPtrs));
    }

    // Allocate + define the state variables.
    // They will be filled in initialData().
    m_statePtr->detachFromQ();
    delete m_velPtr;
    delete m_pPtr;
    delete m_qPtr;
    m_velPtr = new LevelData<FluxBox>;
    m_pPtr   = new LevelData<FArrayBox>;
    m_qPtr   = new LevelData<FArrayBox>;
    m_statePtr->defineAndAttachToQ(*m_velPtr,
                                   *m_pPtr,
                                   *m_qPtr,
                                   grids,
                                   IntVect::Unit,   // vel ghosts
                                   IntVect::Unit,   // p ghosts
                                   IntVect::Unit);  // q ghosts

    // Set up CFI interpolator.
    delete m_cfInterpPtr;
    m_cfInterpPtr = nullptr;
    if (m_level > 0) {
        const DisjointBoxLayout& crseGrids = this->crseNSPtr()->getBoxes();
        m_cfInterpPtr = new CFInterp(*m_levGeoPtr, crseGrids);
    }

    // Each level takes care of the coarser level's flux register.
    if (m_level > 0) {
        crseNSPtr()->m_velFluxRegPtr =
            new FluxRegisterFace(m_levGeoPtr,
                                 crseNSPtr()->m_levGeoPtr,
                                 1,
                                 false,  // doNormalFluxes
                                 true); // doTransverseFluxes

        crseNSPtr()->m_qFluxRegPtr =
            new AnisotropicFluxRegister(this->getBoxes(),
                                        *this->getCrseGridsPtr(),
                                        this->getDomain(),
                                        this->getCrseRefRatio(),
                                        m_qPtr->nComp());
    }

    // Set up the time integrator.
    m_parkPtr = new PARKType(grids,
                             m_velPtr->nComp(),
                             m_velPtr->ghostVect(),
                             m_pPtr->nComp(),
                             m_pPtr->ghostVect(),
                             m_qPtr->nComp(),
                             m_qPtr->ghostVect());

    // Set up IB on this level.
    {
        const auto start = std::chrono::high_resolution_clock::now();
        this->constructIB();
        const auto finish = std::chrono::high_resolution_clock::now();

        const std::chrono::duration<double> elapsed = finish - start;
        pout() << "IB setup time on level " << m_level << ": "
               << Format::textTime(elapsed.count()) << '\n';
    }
}


// -----------------------------------------------------------------------------
// This defines all level solvers on this level.
// This is called AFTER activateLevel() is called during initialization
// and regridding.
// Called from 0 -> lmax, even if just one level changed.
// -----------------------------------------------------------------------------
void
AMRNSLevel::validateOpsAndSolvers()
{
    BEGIN_FLOWCHART();
    CH_assert(m_level == 0);

    constexpr int verbThresh = 4;
    if (s_verbosity >= verbThresh) {
        pout() << "Validating ops and solvers:" << Format::indent() << endl;
    }

    const ProblemContext* ctx = ProblemContext::getInstance();

    bool usingAnyProjector = ctx->proj.doLevelProj;
    usingAnyProjector |= ctx->proj.doInitProj;
    usingAnyProjector |= ctx->proj.doSyncProj;

    // Just in case you need these...
    const size_t lmax  = this->finestNSPtr()->getLevel();
    // const size_t lmin  = m_level;
    // const size_t lbase = (lmin == 0) ? 0 : lmin - 1;

    // Projection stuff.
    if (usingAnyProjector) {
        // Operators. (These act as both level and amr ops.)
        using OpType = Elliptic::AMRMGOperator<LevelData<FArrayBox>>;
        Vector<std::shared_ptr<const OpType>> vAMRMGOps(lmax + 1);

        AMRNSLevel* levPtr = this->coarsestNSPtr();
        while (levPtr) {
            if (!levPtr->m_amrOpsValid) {
                // Redefining this op invalidates its associated level solver.
                levPtr->m_levelSolversValid = false;

                DisjointBoxLayout fineGrids;
                {
                    const AMRNSLevel* finePtr = levPtr->fineNSPtr();
                    if (finePtr) fineGrids = finePtr->getBoxes();
                }

                DisjointBoxLayout crseGrids;
                {
                    const AMRNSLevel* crsePtr = levPtr->crseNSPtr();
                    if (crsePtr) crseGrids = crsePtr->getBoxes();
                }

                levPtr->m_projOpPtr.reset(new PoissonOp(*levPtr->m_levGeoPtr,
                                                        fineGrids,
                                                        crseGrids,
                                                        1,
                                                        this->pressurePhysBC()));

                if (s_verbosity >= verbThresh) {
                    pout() << "m_projOpPtr redefined on level "
                           << levPtr->m_level << endl;
                }
            }

            vAMRMGOps[levPtr->m_level] = levPtr->m_projOpPtr;
            levPtr = levPtr->fineNSPtr();
        }

        // Level solvers.
        levPtr = this->coarsestNSPtr();
        while (levPtr) {
            if (!levPtr->m_levelSolversValid) {
                levPtr->m_levelProjSolverPtr.reset(new LevelProjSolver);
                levPtr->m_levelProjSolverPtr->define(*levPtr->m_projOpPtr);

                if (s_verbosity >= verbThresh) {
                    pout() << "m_levelProjSolverPtr redefined on level "
                           << levPtr->m_level << endl;
                }
            }
            levPtr = levPtr->fineNSPtr();
        }

        // AMR solvers.
        levPtr = this->coarsestNSPtr();
        while (levPtr) {
            if (!levPtr->m_amrSolversValid) {
                levPtr->m_amrProjSolverPtr.reset(new AMRProjSolver);
                levPtr->m_amrProjSolverPtr->define(
                    vAMRMGOps, levPtr->m_level, lmax);

                if (s_verbosity >= verbThresh) {
                    pout() << "m_amrProjSolverPtr redefined on level "
                           << levPtr->m_level << endl;
                }
            }
            levPtr = levPtr->fineNSPtr();
        }
    }

    m_levelOpsValid     = true;
    m_levelSolversValid = true;
    m_amrOpsValid       = true;
    m_amrSolversValid   = true;

    if (s_verbosity >= verbThresh) {
        pout() << Format::unindent << endl;
    }
}
