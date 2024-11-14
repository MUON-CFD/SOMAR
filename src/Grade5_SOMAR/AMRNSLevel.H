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
#ifndef ___AMRNSLevel_H__INCLUDED___
#define ___AMRNSLevel_H__INCLUDED___

#ifndef CH_USE_CUDA_PROJECTOR
#   define CH_USE_CUDA_PROJECTOR 0
#else
#   define CH_USE_CUDA_PROJECTOR 1
#endif

// Projector stuff
#include "PoissonOp.H"
#include "AMRMGSolver.H"
#include "AMRLepticSolver.H"
#include "LevelHybridSolver.H"


#include "AnisotropicAMRLevel.H"
#include "FluxRegisterFace.H"
#include "AnisotropicFluxRegister.H"
#include "CFInterp.H"
#include "SOMAR_Constants.H"
#include "Debug.H"
#include "FiniteDiff.H"
#include "LevelGeometry.H"
#include "State.H"
#include "GNUC_Extensions.H"
#include "SideArray.H"
#include "CubicSpline.H"

// Time stepping
#include "PARK.H"
#include "ForwardEuler.H"
#include "SDC.H"
#include "RKW3CN.H"

// ViscousOp stuff
#include "BCTools.H"
#include "BoxIterator.H"
#include "FABAlgebra.H"
#include "IO.H"
#include "NeighborIterator.H"
#include "ProblemContext.H"
#include "LoadBalance.H"



/**
 * @interface AMRNSLevel
 * @brief     Defines the Navier-Stokes solver at a single AMR level.
 * @author    Edward.Santilli@jefferson.edu
 * @details
 *   This inherits from AnisotropicAMRLevel and PARKRHS, so it can be used
 *   to solve the Navier-Stokes equation via PARK methods on a level.
 *
 *   There are a few functions that are pure virtual. Your problem-specific
 *   physics class should inherit from this and override those functions.
 */
class AMRNSLevel : public AnisotropicAMRLevel, public PARKRHS
{
protected:
    // -------------------------------------------------------------------------
    /// \name AMRNSLevelInit.cpp
    /// These are listed in the order they will be called during a subcycled run.
    /// \{

    /// The constructors are protected because you shouldn't need to make a
    /// AMRNSLevel directly. Use the factory!
    friend class AMRNSLevelFactory;

    /// Partial constructor.
    /// Rest of construction happens when AnisotropicAMR calls define.
    AMRNSLevel();

    /// Off limits. This will throw an error.
    virtual void
    define(AnisotropicAMRLevel* a_coarser_level_ptr,
           const Box&           a_problem_domain,
           int                  a_level,
           const IntVect&       a_ref_ratio) final;

    /// Full virtual constructor
    /// The factory will create the new AMRNSLevel object, then
    /// AnisotropicAMR will call this define function. So, you should expect
    /// all NavierStokes-specific / level-independent parameters to be set.
    virtual void
    define(AnisotropicAMRLevel* a_coarse_level_ptr,
           const ProblemDomain& a_problem_domain,
           int                  a_level,
           const IntVect&       a_ref_ratio);

public:
    /// Virtual destructor
    virtual ~AMRNSLevel();

    /// Initializes this level to have the specified domain a_new_grids.
    virtual void
    initialGrid(const Vector<Box> &a_new_grids);

    /// Performs operations required after the grid has been defined but before
    /// data initialization. This will also be called after readCheckpointLevel
    /// during a restart procedure with argument a_restart set to true.
    /// Called from top -> bottom.
    virtual void
    postInitialGrid(const bool a_restart);

    /// Set initial conditions.
    /// This does some basic setup, but the real work is done by setICs.
    virtual void
    initialData();

    /// Tag cells during initialization.
    virtual void
    tagCellsInit(IntVectSet &a_tags);

    /// Perform initial projection and set up the initial pressure.
    /// This must be done over the entire hierarchy at once, so we wait until
    /// this is called on level 0 and do everything there.
    /// Called from top -> bottom.
    virtual void
    postInitialize();

    /// Returns maximum stable time step for this level with initial data.
    virtual Real
    computeInitialDt();

    /// Frees all memory on this level and clears grids. Called during
    /// construction, destruction, and regridding when a level vanishes. Note
    /// that this also deletes the solvers. To completely re-activate the level,
    /// you need to call activateLevel() AND validateOpsAndSolvers().
    ///
    /// Since it doesn't make sense for level l to be deactivated and level
    /// l+1 to be active, this function will ensure all levels above this level
    /// are properly deactivated as well.
    virtual void
    deactivateLevel();

    /// Allocates and defines structures when this level comes into existence.
    /// This does NOT define the solvers. Do that after level activation
    /// by calling validateOpsAndSolvers().
    ///
    /// The solvers are initialized separately because it may be a very
    /// expensive operation, and activatLevel may be called many times during
    /// construction. The second optional parameter allows the specification of
    /// a processor map, which is otherwise calculated based on load balancing
    virtual void
    activateLevel(const Vector<Box>& a_new_grids,
                  const Vector<int>& a_procMap = Vector<int>());

    /// This defines all level solvers on this level.
    /// This is called AFTER activateLevel() is called during initialization
    /// and regridding.
    /// Called from 0 -> lmax, even if just one level changed.
    virtual void
    validateOpsAndSolvers();
    /// \}

    // -------------------------------------------------------------------------
    /// \name AMRNSLevelRegrid.cpp
    /// \{

    /// Tag cells during runtime.
    virtual void
    tagCells(IntVectSet &a_tags);

    /// Perform any pre-regridding ops -- lBase is the finest unchanged level
    /// This function will be called from the finest level downward to a_lBase.
    /// \note If a new level is about to be created, this function will be run
    /// on that level as well...even if it doesn't exist yet!
    virtual void
    preRegrid(int a_lBase, const Vector<Vector<Box>> &a_newGrids);

    /// Set up data on this level after regridding.
    /// Called from coarsest regridded level --> finest regridded level
    /// Load balancing occurs here.
    virtual void
    regrid(const Vector<Box> &a_new_grids);

    /// perform any post-regridding ops -- lBase is the finest unchanged level
    /// Called from finest regridded level --> a_lBase (coarsest regridded
    /// level - 1)
    virtual void
    postRegrid(int a_lBase);

    /// Assigns procs to the new set of boxes.
    virtual DisjointBoxLayout
    loadBalance(const Vector<Box> &a_grids);

    /// Recomputes the pressure on this level.
    /// This will update m_statePtr->p and set its BCs.
    /// None of the other state variables will be touched.
    /// If this is level 0, we assume p contains an initial guess or zeros.
    virtual void
    initLevelPressure(const Real a_dt);

    /// Recomputes the pressure from this level up to the finest level.
    /// This will only update m_statePtr->p and set its BCs. None of the other
    ///   state variables will be touched.
    /// This function will not check if ctx->rhs.computeInitPressure is true.
    /// This function will take a small, non-subcycled, composite timestep.
    ///   the dt used will be the finest computed dt among all level times
    ///   the multiplier a_dtMult. Typically, a_dtMult ~ 0.5.
    virtual void
    initAMRPressure(const Real a_dtMult);

    /// \}

    // -------------------------------------------------------------------------
    /// \name AMRNSLevelAdvance.cpp
    /// \{

    /// Computes the maximum stable timestep using the current state.
    virtual Real
    computeDt();

    /// advance by one timestep (returns maximum allowable timestep)
    virtual Real
    advance();

    /// Things to do after a timestep.
    /// a_step is the step number that is being completed.
    /// Unlike in AnisotropicAMR.cpp, we start counting from 1.
    /// \todo Maybe this should be merged with postStep. Or we should rename
    /// them to be more descriptive.
    virtual void
    postTimeStep(int a_step);

    /// Refluxes via Q = Q + dQ where
    ///  dQ are the contents of the flux registers.
    /// You must call this from lbase = the coarsest refluxed level.
    virtual void
    simpleRefluxing(Vector<LevelData<FluxBox>*>&   a_amrVel,
                    Vector<LevelData<FArrayBox>*>& a_amrQ);

    /// Refluxes via Q = Q + [1 + dt*D](dQ) where
    ///  D is the diffusive operator and
    ///  dQ are the contents of the flux registers.
    /// You must call this from lbase = the coarsest refluxed level.
    virtual void
    explicitRefluxing(Vector<LevelData<FluxBox>*>&   a_amrVel,
                      Vector<LevelData<FArrayBox>*>& a_amrQ);

    /// Writes basic simulation info to the terminal.
    virtual void
    printDiagnostics(const int a_step, const bool a_writeBanner) const;

    /// \}

    // -------------------------------------------------------------------------
    /// \name AMRNSLevelBC.cpp
    /// \{

    // These functions only need to set the physical BCs. do not worry about
    // CF-interface BCs or periodic BCs. They will be taken care of elsewhere.

    virtual BCTools::BCFunction*
    createVelPhysBC(const int a_dir, const Side::LoHiSide& a_side) const;

    virtual BCTools::BCFunction*
    createPressurePhysBC(const int a_dir, const Side::LoHiSide& a_side) const;

    virtual BCTools::BCFunction*
    createTemperaturePhysBC(const int a_dir, const Side::LoHiSide& a_side) const;

    virtual BCTools::BCFunction*
    createSalinityPhysBC(const int a_dir, const Side::LoHiSide& a_side) const;

    virtual BCTools::BCFunction*
    createScalarsPhysBC(const int a_dir, const Side::LoHiSide& a_side) const;


    // These functions retrieve the physical BC functions that are created and
    // stored by AMRNSLevel.

    inline virtual std::shared_ptr<BCTools::BCFunction>
    velPhysBC() const final
    {
        return m_velBCPtr;
    }

    inline virtual std::shared_ptr<BCTools::BCFunction>
    pressurePhysBC() const final
    {
        return m_pBCPtr;
    }

    inline virtual std::shared_ptr<BCTools::BCFunction>
    temperaturePhysBC() const final
    {
        return m_TBCPtr;
    }

    inline virtual std::shared_ptr<BCTools::BCFunction>
    salinityPhysBC() const final
    {
        return m_SBCPtr;
    }

    inline virtual std::shared_ptr<BCTools::BCFunction>
    scalarsPhysBC() const final
    {
        return m_scalarsBCPtr;
    }


    /// Set physical BCs on the velocity.
    inline virtual void
    setVelPhysBC(LevelData<FluxBox>& a_vel,
                 const Real          a_time,
                 const bool          a_homogBCs) const
    {
        AMRNSLevel::setVelPhysBC(a_vel, a_time, a_homogBCs, AllSides);
    }

    /// Set physical BCs on the pressure.
    inline virtual void
    setPressurePhysBC(LevelData<FArrayBox>& a_p,
                      const Real            a_time,
                      const bool            a_homogBCs) const
    {
        AMRNSLevel::setPressurePhysBC(a_p, a_time, a_homogBCs, AllSides);
    }

    /// Set physical BCs on temperature.
    inline virtual void
    setTemperaturePhysBC(LevelData<FArrayBox>& a_T,
                         const Real            a_time,
                         const bool            a_homogBCs) const
    {
        AMRNSLevel::setTemperaturePhysBC(a_T, a_time, a_homogBCs, AllSides);
    }

    /// Set physical BCs on salinity.
    inline virtual void
    setSalinityPhysBC(LevelData<FArrayBox>& a_S,
                      const Real            a_time,
                      const bool            a_homogBCs) const
    {
        AMRNSLevel::setSalinityPhysBC(a_S, a_time, a_homogBCs, AllSides);
    }

    /// Set physical BCs on user-defined scalars.
    inline virtual void
    setScalarsPhysBC(LevelData<FArrayBox>& a_scalars,
                     const Real            a_time,
                     const bool            a_homogBCs) const
    {
        AMRNSLevel::setScalarsPhysBC(a_scalars, a_time, a_homogBCs, AllSides);
    }


    /// Set physical BCs on velocity.
    /// This assumes a_vel is in the Cartesian (not stretched) basis.
    ///
    /// This is probably NOT the version you want to override in your custom
    /// physics class. Override the non-SideArray version and call this one as
    /// needed.
    virtual void
    setVelPhysBC(LevelData<FluxBox>& a_vel,
                 const Real          a_time,
                 const bool          a_homogBCs,
                 const SideArray&    a_doSides) const final;


    /// This version is useful when you want to call this function to impose the
    /// default BCs on some sides and create custom BCs on other sides.
    ///
    /// This is probably NOT the version you want to override in your custom
    /// physics class. Override the non-SideArray version and call this one as
    /// needed.
    virtual void
    setPressurePhysBC(LevelData<FArrayBox>& a_p,
                      const Real            a_time,
                      const bool            a_homogBCs,
                      const SideArray&      a_doSides) const final;

    /// This version is useful when you want to call this function to impose the
    /// default BCs on some sides and create custom BCs on other sides.
    ///
    /// This is probably NOT the version you want to override in your custom
    /// physics class. Override the non-SideArray version and call this one as
    /// needed.
    virtual void
    setTemperaturePhysBC(LevelData<FArrayBox>& a_T,
                         const Real            a_time,
                         const bool            a_homogBCs,
                         const SideArray&      a_doSides) const final;

    /// This version is useful when you want to call this function to impose the
    /// default BCs on some sides and create custom BCs on other sides.
    ///
    /// This is probably NOT the version you want to override in your custom
    /// physics class. Override the non-SideArray version and call this one as
    /// needed.
    virtual void
    setSalinityPhysBC(LevelData<FArrayBox>& a_S,
                      const Real            a_time,
                      const bool            a_homogBCs,
                      const SideArray&      a_doSides) const final;

    /// Set physical BCs on all scalars.
    virtual void
    setScalarsPhysBC(LevelData<FArrayBox>& a_s,
                     const Real            a_time,
                     const bool            a_homogBCs,
                     const SideArray&      a_doSides) const final;

    // These functions set ALL BCs -- physical, CFI, and exchange.
    // These are not the physical BCs that the user should override. They
    // are simply convenience functions so that we don't need to keep defining
    // interpolators, etc.

    /// Resets all BCs on all state variables. This will involve several
    /// exchanges and CFI interpolations. Use sparingly.
    /// NOTE: Unlike the other setBC functions, this assumes vel is in the
    ///  Cartesian basis!
    virtual void
    setBC(State& a_state, const Real a_time) const;

    /// Set all BCs on velocity.
    /// This assumes a_vel is in the curvilinear (mapped) basis.
    virtual void
    setVelBC(LevelData<FluxBox>&       a_vel,
             const Real                a_time,
             const bool                a_homogBCs,
             const LevelData<FluxBox>* a_crseVelPtr = nullptr) const;

    /// Set all BCs on the pressure.
    virtual void
    setPressureBC(LevelData<FArrayBox>&       a_p,
                  const Real                  a_time,
                  const bool                  a_homogBCs,
                  const LevelData<FArrayBox>* a_crsePPtr = nullptr) const;

    /// Set all BCs on temperature.
    virtual void
    setTemperatureBC(LevelData<FArrayBox>&       a_state,
                     const Real                  a_time,
                     const bool                  a_homogBCs,
                     const LevelData<FArrayBox>* a_crseBPtr = nullptr) const;

    /// Set all BCs on salinity.
    virtual void
    setSalinityBC(LevelData<FArrayBox>&       a_state,
                  const Real                  a_time,
                  const bool                  a_homogBCs,
                  const LevelData<FArrayBox>* a_crseBPtr = nullptr) const;

    /// Set all BCs on all scalars.
    virtual void
    setScalarBC(LevelData<FArrayBox>&       a_s,
                const Real                  a_time,
                const bool                  a_homogBCs,
                const LevelData<FArrayBox>* a_crseSPtr = nullptr) const;

    // -------------------------------------------------------------------------
    /// \name AMRNSLevelFill.cpp
    /// \{

    // Fill with time-interpolated state...

    /// Fills a data holder with the Cartesian-based velocity.
    /// All ghosts will be filled, even at edges and vertices.
    virtual void
    fillVelocity(LevelData<FluxBox>& a_vel,
                 const Real          a_time,
                 const bool          a_setBCs = true) const;

    /// Fills a data holder with the pressure. All ghosts will be filled,
    /// even at edges and vertices.
    virtual void
    fillPressure(LevelData<FArrayBox>& a_p,
                 const Real            a_time,
                 const bool            a_setBCs = true) const;

    /// Fills a data holder with the temperature. All ghosts will be filled,
    /// even at edges and vertices.
    virtual void
    fillTemperature(LevelData<FArrayBox>& a_state,
                    const Real            a_time,
                    const bool            a_setBCs           = true,
                    const bool            a_removeBackground = false) const;

    /// Fills a data holder with the salinity. All ghosts will be filled,
    /// even at edges and vertices.
    virtual void
    fillSalinity(LevelData<FArrayBox>& a_state,
                 const Real            a_time,
                 const bool            a_setBCs           = true,
                 const bool            a_removeBackground = false) const;

    /// Fills a data holder with the scalars. All ghosts will be filled,
    /// even at edges and vertices.
    virtual void
    fillScalar(LevelData<FArrayBox>& a_s,
               const Real            a_time,
               const bool            a_setBCs = true) const;

    /// Computes advVel^i = J * dXi^i/dx^i * cartVel^i (No sum)
    /// a_cartVel must have exactly 1 component.
    /// This function can work in-place.
    virtual void
    sendToAdvectingVelocity(LevelData<FluxBox>&       a_advVel,
                            const LevelData<FluxBox>& a_cartVel) const;

    /// Undo the sendToAdvectingVelocity operation. Result will be in the
    /// Cartesian basis. This function can work in-place.
    virtual void
    sendToCartesianVelocity(LevelData<FluxBox>&       a_cartVel,
                            const LevelData<FluxBox>& a_advVel) const;

    /// \brief Fills a data holder with the buoyancy given temperature and
    /// salinity.
    ///
    /// \details
    /// The FABs will be at this level's index space, but may be defined over a
    /// patch that is not in this level's DisjointBoxLayout. Therefore, you
    /// can compare the FAB's boxes to this level's ProblemDomain to figure out
    /// where the ghosts are, but you should not loop over the grids to find
    /// an associated DataIndex. For any reasonable EoS, this should not be
    /// a problem.
    ///
    /// T and S must be the total temperature and salinity, not just the
    /// deviation from the background. Likewise, this function returns the total
    /// buoyancy.
    ///
    /// By default, this uses a simple linear equation of state:
    ///  b = -\alpha*g*T + \beta*g*S
    /// with \alpha = 2.0e-4 K^-1, \beta = 8.0e-4 psu^-1, and g = 9.81 m/s^2.
    virtual void
    equationOfState(FArrayBox&       a_bFAB,
                    const FArrayBox& a_TFAB,
                    const FArrayBox& a_SFAB,
                    const FArrayBox& a_zFAB) const;

    // Misc fill functions...

    /// Fills a data holder with the total energy on a single level.
    /// Don't expect conservation on the upper levels as the grids change.
    virtual void
    computeTotalEnergy(LevelData<FArrayBox>& a_energy, const Real a_time) const;

    /// Adds the kinetic energy to a data holder on a single level.
    virtual void
    addKE(LevelData<FArrayBox>&     a_energy,
          const LevelData<FluxBox>& a_cartVel) const;

    /// Adds the potential energy to a data holder on a single level.
    /// Use this when there is no background stratification.
    virtual void
    addGPE(LevelData<FArrayBox>&       a_energy,
           const LevelData<FArrayBox>& a_btot,
           const LevelData<FArrayBox>& a_z) const;

    /// Adds the available potential energy to a data holder on a single level.
    /// Use this when Nsq != 0.0.
    virtual void
    addAPE(LevelData<FArrayBox>&       a_energy,
           const LevelData<FArrayBox>& a_bpert) const;

    /// In SOMAR, the "pressure" is actually whatever it needs to be in order to
    /// make the velocity divergenceless. It is usually an okay estimate of the
    /// physical pressure in the domain's interior, but not at the boundaries
    /// where homogeneous Neumann BCs have been imposed. This function computes
    /// the physical pressure properly.
    ///
    /// This must be called from level 0 once all levels are synchronized.
    ///
    /// \warning This is a work in progress! I need to reflux in AMR and I
    /// haven't checked if it produces more accurate results than the projector
    /// yet.
    virtual void
    computePhysicalPressure(Vector<LevelData<FArrayBox>*>& a_amrPressure,
                            const Real                     a_oldTime,
                            const Real                     a_curTime) const;
    /// \}


    // -------------------------------------------------------------------------
    /// \name AMRNSLevelStrat.cpp
    /// \{
    /// \brief Sets the T and S background stratification at specified heights.
    /// \details Provide as many elements as you need to resolve the
    /// stratification. This function will be called just once on level 0.
    /// a_zmin and a_zmax are provided so that you know where the extrapolation
    /// will take place. Keep in mind this range may include some space for
    /// ghost values.
    virtual void
    setStratification(Vector<Real>& a_vTbar,
                      Vector<Real>& a_vSbar,
                      Vector<Real>& a_vz,
                      const Real    a_zmin,
                      const Real    a_zmax) const;

    /// \brief Sets the stratification in b, T, and S and computes the vertical
    /// structure.
    /// \details To do it's job, this function calls setStratification(...)
    /// and equationOfState(...). Those are the functions you want to override
    /// in your physics class, not this...unless you know what you are doing
    /// and need to do something exotic.
    ///
    /// This will only set the members once per level. If the members are not
    /// nullptr, this function just returns.
    virtual void
    setStratificationMembers();

    /// Solve for the internal wave speed and structure functions due to Nsq.
    ///
    /// When called on level l, you can be sure that the stratification on
    /// levels [0,l-1] are well defined.
    virtual void
    computeStructureFunctions(Real&            a_c1,
                              FArrayBox&       a_phi1FAB,
                              const FArrayBox& a_NsqFAB) const;
    /// \}


    // -------------------------------------------------------------------------
    /// \name AMRNSLevelSponge.cpp
    /// \{

    /// Static utility.
    /// Splits the domain into its sponge layers and the interior. The locations
    /// of the splitting are given as face indices. If the domain is not split,
    /// these indices will lie outside of the domain.
    static void
    computeSpongeRegions(Tuple<Box, 2>&       a_spongeBox,
                         Tuple<int, 2>&       a_splitFaceIndex,
                         Box&                 a_interior,
                         const int            a_dir,
                         const LevelGeometry& a_levGeo);

    /// By default, this does nothing. You can override this in your physics
    /// class if you need to perform global computations before filling targets.
    virtual void
    prepareSpongeTargets(const LevelData<FluxBox>&   a_vel,
                         const LevelData<FArrayBox>& a_p,
                         const LevelData<FArrayBox>& a_q,
                         const Real                  a_time);

    /// By default, this throws an error for all BCs except TIDAL.
    /// You can override this to fit your own needs.
    virtual void
    fillVelSpongeTarget(FArrayBox&            a_targetFAB,
                        const int             a_targetFABComp,
                        const int             a_velComp,
                        const FArrayBox&      a_physCoorFAB,
                        const Box&            a_spongeBox,
                        const DataIndex&      a_di,
                        const Real            a_time,
                        const int             a_bdryDir,
                        const Side::LoHiSide& a_side) const;

    /// By default, this sets T and S to zero and does nothing
    /// to the other q comps. You can override this to fit your own needs.
    virtual void
    fillQSpongeTarget(FArrayBox&            a_targetFAB,
                      const int             a_targetFABComp,
                      const int             a_qComp,
                      const FArrayBox&      a_physCoorFAB,
                      const Box&            a_spongeBox,
                      const DataIndex&      a_di,
                      const Real            a_time,
                      const int             a_bdryDir,
                      const Side::LoHiSide& a_side) const;

    /// Adds ramp * (target - state) * invTimeScale to the forcing, a_k*.
    /// invTimeScale is typically 1.0 / (c * dt), where c ~ 10 or so.
    /// By default, this only operates on vel, T, and S. Override if you like.
    virtual void
    addSpongeForcing(LevelData<FluxBox>&         a_kvel,
                     LevelData<FArrayBox>&       a_kq,
                     const LevelData<FluxBox>&   a_vel,
                     const LevelData<FArrayBox>& a_p,
                     const LevelData<FArrayBox>& a_q,
                     const Real                  a_time,
                     const Real                  a_invTimeScale);

    /// a_ratio is the fractional distance through the sponge, a measure of the
    /// distance to the target boundary. You can use this function to make a
    /// linear ramp (a_ratio), a quadratic ramp (a_ratio * a_ratio), and so on.
    /// Default = linear.
    inline Real
    spongeRampFunction(Real a_ratio) const
    {
        return a_ratio * a_ratio;
    }

    /// \}

    // -------------------------------------------------------------------------
    /// \name AMRNSLevelSGS.cpp
    /// \{

    // /// @brief Computes JS^{ia}, where the last index is Cartesian-based.
    // /// @param a_JSia            Note that by default, this is scaled by J!
    // /// @param a_cartVel
    // /// @param a_multByJ         True by default!
    // /// @param a_makeCovariant   False by default. If true, JS_{ia} is returned.
    // /// @param a_primaryScale    The coeff of the Grad[cartVel] term.
    // /// @param a_transposeScale  The coeff of the Grad[cartVel]^T term.
    // /// @note Careful, you may need to run the result through
    // ///       LayoutTools::averageOverlappingValidFaces.
    // virtual void
    // rateOfStrain(StaggeredFluxLD&          a_JSia,
    //              const LevelData<FluxBox>& a_cartVel,
    //              const bool                a_multByJ        = true,
    //              const bool                a_makeCovariant  = false,
    //              const Real                a_primaryScale   = 0.5,
    //              const Real                a_transposeScale = 0.5) const;

    // /// @brief Computes T^{ia} = 2 * J * (nu + eddyNu) * S^{ia},
    // ///        where the last index is Cartesian-based.
    // /// @param a_Tia    Can be the same holder as a_Sia for an in-place calc.
    // /// @param a_JSia   The results of rateOfStrain(...) with default options.
    // /// @param a_nu     The kinematic viscosity in each direction.
    // /// @param a_eddyNu The CC eddy viscosity. Zero if not using an LES.
    // virtual void
    // viscousStress(StaggeredFluxLD&            a_Tia,
    //               const StaggeredFluxLD&      a_JSia,
    //               const RealVect&             a_nu,
    //               const LevelData<FArrayBox>& a_eddyNu) const;

    // This file only contains LES-specific routines. This does not include the
    // compute*Diffusion(...) functions, which are needed even without the LES.
    // Those are in AMRNSLevelRHS.cpp.

    /// @brief Computes eddyNu given a velocity field.
    /// @param a_eddyNu
    /// @param a_cartVel
    /// @param a_time
    virtual void
    computeEddyNu(LevelData<FArrayBox>&     a_eddyNu,
                  const LevelData<FluxBox>& a_cartVel,
                  const Real                a_time) const;

    /// CC Laplacian filter used by Ducros, et. al.
    /// DOI: http://dx.doi.org/10.1017/S0022112096008221
    virtual void
    LaplacianFilter(LevelData<FArrayBox>& a_ccCartVel,
                    const int             a_numFilterSweeps,
                    const RealVect&       a_dirScale) const;

    /// Computes the eddy viscosity using the method of Ducros, et. al.,
    /// Journal of Fluid Mechanics / Volume 326 / November 1996, pp 1- 36
    /// DOI: http://dx.doi.org/10.1017/S0022112096008221
    ///
    /// The filter is the Laplacian filter described in the paper.
    ///
    /// a_dirScale = RealVect(1.0, 1.0, 1.0) for the 6-neighbor formulation and
    /// a_dirScale = RealVect(1.0, 1.0, 0.0) for the 4-neighbor formulation.
    /// You are not restricted to these two cases, but they are most common.
    virtual void
    SGSModel_Ducros(LevelData<FArrayBox>&     a_nuT,
                    const LevelData<FluxBox>& a_cartVel,
                    const Real                a_time,
                    const int                 a_numFilterSweeps,
                    const RealVect&           a_dirScale) const;
    /// \}


    // -------------------------------------------------------------------------
    /// \name AMRNSLevelRHS.cpp
    /// \{

    /// ARKRHS override.
    /// The explicit Navier-Stokes forcing. Be sure to set the comps with no
    /// forcing to zero. Q is non-const so that you can set BCs.
    ///
    /// Both a_vel and a_kvel are in the Cartesian basis. If you need
    /// a_vel in the mapped basis, use LevelGeometry::sendToMappedBasis on
    /// a local copy.
    virtual void
    setExplicitRHS(LevelData<FluxBox>&   a_kvel,
                   LevelData<FArrayBox>& a_kq,
                   LevelData<FluxBox>&   a_vel,
                   LevelData<FArrayBox>& a_p,
                   LevelData<FArrayBox>& a_q,
                   const Real            a_time,
                   const Real            a_refluxDt);

    /// \brief
    ///  Override this to modify the Navier-Stokes forces in a_kvel and a_kq.
    /// \details
    ///  This is called by setExplicitRHS *after* a_kvel and a_kq are filled
    ///  with the typical Navier-Stokes forces, but before the Immersed Boundary
    ///  forces are computed/applied. Override this function to add your own
    ///  custom forces or modify what is already there.
    ///
    ///  Both a_kvel and a_vel are Cartesian-based. Unlike setExplicitRHS,
    ///  this function also receives a_advVel, the contravariant, J-scaled
    ///  velocity that was used by setExplicitRHS to advect a_vel and a_q.
    ///  Feel free to use it as-is, or use a_advVel as workspace. You can
    ///  clobber its data if you like, but you should NOT clobber a_vel, a_p,
    ///  or a_q! Those are non-const do that you can change the BCs if needed.
    ///
    ///  By default, this function does nothing.
    virtual void
    addExplicitRHS(LevelData<FluxBox>&   a_kvel     __nowarn_unused,
                   LevelData<FArrayBox>& a_kq       __nowarn_unused,
                   LevelData<FluxBox>&   a_vel      __nowarn_unused,
                   LevelData<FArrayBox>& a_p        __nowarn_unused,
                   LevelData<FArrayBox>& a_q        __nowarn_unused,
                   LevelData<FluxBox>&   a_advVel   __nowarn_unused,
                   const Real            a_time     __nowarn_unused,
                   const Real            a_refluxDt __nowarn_unused)
    {
    }


    /// ARKRHS override.
    /// The implicit Navier-Stokes forcing. Be sure to set the comps with no
    /// forcing to zero. Q is non-const so that you can set BCs.
    ///
    /// Both a_vel and a_kvel are in the Cartesian basis. If you need
    /// a_vel in the mapped basis, use LevelGeometry::sendToMappedBasis on
    /// a local copy.
    virtual void
    setImplicitRHS(LevelData<FluxBox>&   a_kvel,
                   LevelData<FArrayBox>& a_kq,
                   LevelData<FluxBox>&   a_vel,
                   LevelData<FArrayBox>& a_p,
                   LevelData<FArrayBox>& a_q,
                   const Real            a_time,
                   const Real            a_refluxDt);

    /// Solves [1 - gammaDt * D] Q^{n+1} = Q^{n} in place.
    ///
    /// a_vel is in the Cartesian basis.
    virtual void
    solveImplicit(LevelData<FluxBox>&   a_vel,
                  LevelData<FArrayBox>& a_q,
                  const Real            a_gammaDt,
                  const Real            a_time,
                  const Real            a_refluxDt);

    /// Computes viscous force and fluxes for the FC momentum (aka velocity).
    /// a_kvel and a_momentumFlux will be overwritten.
    /// Flux registers will not be updated. That is left to the caller.
    /// Result will not be scaled by 1/J. That is left to the caller.
    virtual void
    computeMomentumDiffusion(LevelData<FluxBox>&         a_kvel,
                             StaggeredFluxLD&            a_momentumFlux,
                             const LevelData<FluxBox>&   a_cartVel,
                             const RealVect&             a_nu,
                             const LevelData<FArrayBox>& a_eddyNu,
                             const Real                  a_primaryScale   = 1.0,
                             const Real                  a_transposeScale = 1.0) const;

    /// Computes diffusive force and fluxes for any CC scalar.
    /// a_kq and a_qFlux will be overwritten.
    /// Flux registers will not be updated. That is left to the caller.
    /// Result will not be scaled by 1/J. That is left to the caller.
    virtual void
    computeScalarDiffusion(LevelData<FArrayBox>&       a_kq,
                           LevelData<FluxBox>&         a_qFlux,
                           const LevelData<FArrayBox>& a_q,
                           const Real                  a_kappa,
                           const LevelData<FArrayBox>& a_eddyNu,
                           const Real                  a_eddyPrandtl) const;

    /// Solves [1 - gammaDt * D] cartVel^{n+1} = cartVel^{n} in place.
    // a_totalNu = nu + eddyNu.
    virtual void
    solveMomentumDiffusion(LevelData<FluxBox>&         a_cartVel,
                           const LevelData<FArrayBox>& a_totalNu,
                           const Real                  a_gammaDt,
                           const Real                  a_time) const;

    /// Solves [1 - gammaDt * D.(kappa + eddyKappa).G] phi^{n+1} = phi^{n}
    /// in place, where eddyKappa = eddyNu / eddyPrandtl.
    virtual void
    solveScalarDiffusion(
        LevelData<FArrayBox>&                       a_phi,
        const LevelData<FArrayBox>*                 a_crsePhiPtr,
        const std::vector<Real>&                    a_vKappa,
        const LevelData<FArrayBox>&                 a_eddyNu,
        const std::vector<Real>&                    a_vEddyPrandtl,
        const Real                                  a_gammaDt,
        const Real                                  a_time,
        const std::shared_ptr<BCTools::BCFunction>& a_bcFuncPtr,
        const std::string                           a_scalarName) const;

    /// Increments both the registers on this and the coarser level as needed.
    /// This function can be used to reflux the FC momentum, but not the
    /// CC scalars.
    virtual void
    incrementFluxRegisters(const StaggeredFluxLD& a_flux,
                           const Real             a_refluxDt);

    /// Increments both the registers on this and the coarser level as needed.
    /// This function can be used to reflux any of the CC scalars, but not the
    /// FC momentum.
    virtual void
    incrementFluxRegisters(const LevelData<FluxBox>& a_flux,
                           const Interval&           a_regInterval,
                           const Real                a_refluxDt);

    /// ARKRHS override.
    /// This function is called at the end of RKW3CN::advance(...) so that the
    /// user can perform postprocessing. If not overridden in a physics class, it
    /// does nothing.
    virtual void
    postStep(const LevelData<FluxBox>&   a_vel  __nowarn_unused,
             const LevelData<FArrayBox>& a_p    __nowarn_unused,
             const LevelData<FArrayBox>& a_q    __nowarn_unused,
             const Real                  a_time __nowarn_unused) const
    {
    }
    /// \}


    // -------------------------------------------------------------------------
    /// \name AMRNSLevelAdvect.cpp
    /// \{

    /// Computes advective force and fluxes for the FC momentum (aka velocity).
    /// a_kvel will be ***incremented***, not overwritten.
    /// a_momentumFlux will be overwritten.
    /// Flux registers will not be updated. That is left to the caller.
    /// Result will not be scaled by 1/J. That is left to the caller.
    virtual void
    computeMomentumAdvection(LevelData<FluxBox>&       a_kvel,
                             StaggeredFluxLD&          a_momentumFlux,
                             const LevelData<FluxBox>& a_cartVel,
                             const LevelData<FluxBox>& a_advVel) const;

    ///
    virtual void
    momAdvFlux_Cons_Form_SecondOrder(
        StaggeredFluxLD&          a_momentumFlux,
        const LevelData<FluxBox>& a_cartVel,
        const LevelData<FluxBox>& a_advVel) const;

    ///
    virtual void
    momAdvFlux_Cons_Form_FourthOrder(
        StaggeredFluxLD&          a_momentumFlux,
        const LevelData<FluxBox>& a_cartVel,
        const LevelData<FluxBox>& a_advVel) const;

    ///
    virtual void
    momAdvForce_AdvForm_SecondOrder(
        LevelData<FluxBox>&       a_kvel,
        const LevelData<FluxBox>& a_cartVel,
        const LevelData<FluxBox>& a_advVel,
        const Real                a_scale = 1.0) const;

    /// Computes advective force and fluxes for any CC scalar.
    /// a_kq and a_qFlux will be overwritten.
    /// Flux registers will not be updated. That is left to the caller.
    /// Result will not be scaled by 1/J. That is left to the caller.
    virtual void
    computeScalarAdvection(LevelData<FArrayBox>&       a_kq,
                           LevelData<FluxBox>&         a_qFlux,
                           const LevelData<FArrayBox>& a_q,
                           const LevelData<FluxBox>&   a_advVel) const;

    ///
    virtual void
    scalAdvFlux_Cons_Form_SecondOrder(
        LevelData<FluxBox>&         a_qFlux,
        const LevelData<FArrayBox>& a_q,
        const LevelData<FluxBox>&   a_advVel) const;

    ///
    virtual void
    scalAdvFlux_Cons_Form_FourthOrder(
        LevelData<FluxBox>&         a_qFlux,
        const LevelData<FArrayBox>& a_q,
        const LevelData<FluxBox>&   a_advVel) const;
    /// \}


    // -------------------------------------------------------------------------
    /// \name AMRNSLevelProject.cpp
    /// \{

    /// ARKRHS override.
    /// Use an approximation of p to project. Do not solve for p.
    /// The p provided is the lagged pressure computed in the previous RK stage.
    ///
    /// This is a single-level projector.
    virtual void
    projectPredict(LevelData<FluxBox>&   a_vel,
                   LevelData<FArrayBox>& a_p,
                   const Real            a_time,
                   const Real            a_projDt);

    /// ARKRHS override.
    /// Solve for p and project. This will be an incremental solve
    /// since we already called projectPredict using the lagged p.
    ///
    /// To project the velocity without affecting p, set projDt = 0.
    ///
    /// This is a single-level projector.
    virtual void
    projectCorrect(LevelData<FluxBox>&   a_vel,
                   LevelData<FArrayBox>& a_p,
                   const Real            a_time,
                   const Real            a_projDt,
                   const Real            a_gammaDt = 0.0);

    /// Project m_statePtr->vel the top level to *this.
    /// If a_sync is true, m_statePtr->p will be updated. Otherwise, the phi
    /// generated in this projection will be thrown away.
    ///
    /// If you call this function, vel will be projected no matter what. This
    /// function will not check the problem parameters. If you don't want to use
    /// the projector, wrap your call to this function in an if statement.
    virtual void
    projectDownToThis(const bool a_sync = false);

    /// \}

    // -------------------------------------------------------------------------
    /// \name AMRNSLevelIO.cpp
    /// \{

    /// write checkpoint header
    virtual void
    writeCheckpointHeader(HeaderData&        a_handle,
                          const std::string& a_fileName) const;

    /// write checkpoint data for this level
    virtual void
    writeCheckpointLevel(const std::string& a_fileName, int level) const;

    /// read checkpoint header
    virtual void
    readCheckpointHeader(const std::string& a_filename);

    /// read checkpoint data for this level
    virtual void
    readCheckpointLevel(const std::string& a_filename);

    /// Open plotfile. Only called on level 0.
    virtual void
    openFile(const std::string& a_filename, const bool checkpoint) const;

    /// Close plotfile. Only called on level 0.
    virtual void
    closeFile(const std::string& a_filename ) const;

    /// Write plotfile header. Only called on level 0.
    virtual void
    writePlotHeader(HeaderData& a_handle, const std::string& a_fileName) const;

    /// Write plotfile data for this level. Called from bottom->top.
    virtual void
    writePlotLevel(const std::string& a_fileName, int level) const;

    /// \}

    // -------------------------------------------------------------------------
    /// \name AMRNSLevelUtil.cpp
    /// \{

    /// What level is *this at?
    inline virtual int
    getLevel() const
    {
        return m_level;
    }

    /// Returns a pointer to a specified level or NULL if level DNE.
    virtual AMRNSLevel*
    getLevel(const int a_level);

    /// Returns a pointer to a specified level or NULL if level DNE.
    /// (const version)
    virtual const AMRNSLevel*
    getLevel(const int a_level) const;

    /// Returns a pointer to a specified level or NULL if level DNE.
    virtual AMRNSLevel*
    getLevel(const ProblemDomain& a_domain);

    /// Returns a pointer to a specified level or NULL if level DNE.
    /// (const version)
    virtual const AMRNSLevel*
    getLevel(const ProblemDomain& a_domain) const;

    /// Casts the fine level ptr to an AMRNSLevel object.
    virtual AMRNSLevel*
    fineNSPtr() const;

    /// Casts the coarse level ptr to an AMRNSLevel object.
    virtual AMRNSLevel*
    crseNSPtr() const;

    /// Returns a pointer to the finest level
    virtual AMRNSLevel*
    finestNSPtr();

    /// Returns a pointer to the finest level (const version)
    virtual const AMRNSLevel*
    finestNSPtr() const;

    /// Returns a pointer to the coarsest level
    virtual AMRNSLevel*
    coarsestNSPtr();

    /// Returns a pointer to the coarsest level (const version)
    virtual const AMRNSLevel*
    coarsestNSPtr() const;

    /// Returns the smallest grid size in the entire AMR hierarchy.
    /// This searches entire domains, not just extant grids.
    /// This returns a physical dx, not the mapped dXi.
    virtual Real
    getMinAMRDx(const int a_coorDir) const;

    /// Is this this finest level?
    virtual bool
    isFinestLevel() const;

    /// Is this an empty level?
    virtual bool
    isEmpty() const;

    /// Returns this level's DisjointBoxLayout.
    virtual const DisjointBoxLayout&
    getBoxes() const;

    /// Returns this level's ProblemDomain.
    virtual const ProblemDomain&
    getDomain() const;

    /// Returns this level's ProblemDomain as a simple Box.
    virtual const Box&
    getDomainBox() const;

    /// Returns a pointer to the coarser grids. nullptr if this is level 0.
    virtual const DisjointBoxLayout*
    getCrseGridsPtr() const;

    /// Returns a pointer to the finer grids. nullptr if this is the finest level.
    virtual const DisjointBoxLayout*
    getFineGridsPtr() const;

    /// Retrieves the State.
    virtual State&
    getState();

    /// Retrieves the State.
    /// (const version)
    virtual const State&
    getState() const;

    /// Retrieves the LevelGeometry.
    virtual const LevelGeometry&
    getLevGeo() const;

    /// Packages all LevGeos in the AMR hierarchy into a vector.
    virtual Vector<const LevelGeometry*>
    getAMRLevGeos () const;

    /// \brief Allocates and defines m_amrData on the indicated levels.
    /// \details
    /// If a_lmax = -1, then we allocate up to the top level.
    /// All elements outside of the range [a_lbase, a_lmax] will be nullptr.
    ///
    /// The resulting vector's indices will correspond to the levels and the
    /// vector will have as many elements as there are extant levels.
    ///
    /// If the levels in the range [a_lbase, a_lmax] are not at a_time, an
    /// error will be thrown.
    ///
    /// It is the caller's responsibility to delete the allocated memory!
    virtual void
    allocateAndDefine(Vector<LevelData<FArrayBox>*>& a_amrData,
                      const int                      a_numComps,
                      const IntVect& a_ghostVect = IntVect::Zero,
                      int            a_lbase     = 0,
                      int            a_lmax      = -1) const;

    /// FC version of allocateAndDefine. See it's comments for details.
    virtual void
    allocateAndDefine(Vector<LevelData<FluxBox>*>& a_amrData,
                      const int                    a_numComps,
                      const IntVect&               a_ghostVect = IntVect::Zero,
                      int                          a_lbase     = 0,
                      int                          a_lmax      = -1) const;

    /// \brief Allocates and aliases components of m_pPtr on the indicated
    /// levels.
    /// \details
    /// If a_ivl = Interval(), then we collect all comps.
    /// If a_lmax = -1, then we collect up to the top level.
    /// All elements outside of the range [a_lbase, a_lmax] will be nullptr.
    ///
    /// The resulting vector's indices will correspond to the levels and the
    /// vector will have as many elements as there are extant levels.
    ///
    /// If the levels in the range [a_lbase, a_lmax] are not at a_time, an
    /// error will be thrown.
    ///
    /// It is the caller's responsibility to delete the allocated memory!
    virtual void
    allocateAndAliasPressure(Vector<LevelData<FArrayBox>*>& a_amrData,
                             const Real                     a_time,
                             Interval                       a_ivl = Interval(),
                             int                            a_lbase = 0,
                             int                            a_lmax  = -1);

    /// m_qPtr version of allocateAndAlias. See it's comments for details.
    virtual void
    allocateAndAliasScalars(Vector<LevelData<FArrayBox>*>& a_amrData,
                            const Real                     a_time,
                            Interval                       a_ivl   = Interval(),
                            int                            a_lbase = 0,
                            int                            a_lmax  = -1);

    /// m_velPtr version of allocateAndAlias. See it's comments for details.
    virtual void
    allocateAndAliasVel(Vector<LevelData<FluxBox>*>& a_amrData,
                        const Real                   a_time,
                        Interval                     a_ivl   = Interval(),
                        int                          a_lbase = 0,
                        int                          a_lmax  = -1);

    /// Deallocates Vectors created by allocateAndAlias. CC version.
    /// This would be better as a template that works on any Vector<T*>.
    virtual void deallocate(Vector<LevelData<FArrayBox>*> &a_amrData) const;

    /// Deallocates Vectors created by allocateAndAlias. FC version.
    /// This would be better as a template that works on any Vector<T*>.
    virtual void deallocate(Vector<LevelData<FluxBox>*> &a_amrData) const;

    /// \brief Average data down. CC version.
    /// \details This function starts at the highest, non-null entry and ends at
    /// a_lmin or the coarsest non-null level. It is okay to let a_lmin be < 0.
    /// The indices of a_amrData MUST correspond to the level numbers.
    virtual void
    averageDown(Vector<LevelData<FArrayBox>*>& a_amrData,
                int                            a_lmin,
                const bool                     a_useJWeighting) const;

    /// Average data down. FC version.
    /// See the CC version's comments for details.
    virtual void
    averageDown(Vector<LevelData<FluxBox>*>& a_amrData,
                int                          a_lmin,
                const bool                   a_useJWeighting) const;

    /// \brief Average *m_statePtr down from a_lmax to *this.
    /// \details This function looks at the input file parameters and averages
    /// down eddyNu and p as instructed.
    ///
    /// If a_lmax < 0, then the finest extant level will be used.
    ///
    /// This will not reset BCs or exchange after performing the average.
    /// that is up to you.
    virtual void
    averageDownToThis(const int a_lmax = -1);

    /// Set BCs on *m_statePtr on all levels from a_lmax to *this.
    /// See averageDownToThis for details on a_lmax.
    virtual void
    setBCsDownToThis(const int a_lmax = -1);
    /// \}


    // -------------------------------------------------------------------------
    /// \name Must be overridden in a physics class
    /// \{

    // Custom scalars...

    // Custom scalars will automatically be advected and, if requested,
    // diffused. If your needs are more specialized than that, then it
    // is up to you to define a forcing function (see below).

    /// This tells the Navier-Stokes class how many additional cell-centered
    /// scalars to create. Do not include T, S, or eddyNu.
    virtual int numScalars() const = 0;

    /// Returns the names of the CC scalars.
    /// This will be used to put scalar names in the plot files.
    virtual std::string getScalarName(const int a_comp) const = 0;

    // ICs...

    /// Set the ICs on all fields at once. This includes the velocity, buoyancy,
    ///  eddy viscosity and diffusivity, and all scalars. Note that I did not
    ///  say pressure. Upon entry, the pressure will already be set to it's
    ///  initial guess of 0. If you have a better estimate, feel free to provide
    ///  it. We will use it as an initial guess for the projector.
    /// By the time this function is called, *m_levGeoPtr will be ready for use.
    virtual void setICs(State &a_state) = 0;

    /// \}

protected:
    /// If true and m_isDefined is false, then we are in the state between
    /// weak constructor and the full virtual constructor calls.
    /// m_isDefined is defined and set in AnisotropicAMRLevel.
    bool m_isFactoryDefined;

    /// Remember, this being true does not mean the ops or solvers are valid.
    bool m_isActivated;
    /// Set to false if any level has changed from [lev-1, lev].
    bool m_levelOpsValid;
    /// Set to false if any level has changed from [lev-1, lev].
    bool m_levelSolversValid;
    /// Set to false if any level has changed from [lev-1, lev+1].
    bool m_amrOpsValid;
    /// Set to false if any level has changed from [lev-1, lmax].
    bool m_amrSolversValid;

    /// This level's geometry/coordinate system.
    GeoSourceInterface* m_geoSrcPtr;
    LevelGeometry*      m_levGeoPtr;

    /// \name BC setting objects
    /// \{
    std::shared_ptr<BCTools::BCFunction> m_velBCPtr;
    std::shared_ptr<BCTools::BCFunction> m_pBCPtr;
    std::shared_ptr<BCTools::BCFunction> m_TBCPtr;
    std::shared_ptr<BCTools::BCFunction> m_SBCPtr;
    std::shared_ptr<BCTools::BCFunction> m_scalarsBCPtr;

    /// m_physBdryIter[dir][side] = PhysBdryIter for just one domain side.
    mutable PhysBdryIter                                      m_allPhysBdryIter;
    mutable std::array<std::array<PhysBdryIter, 2>, SpaceDim> m_singlePhysBdryIter;
    /// \}

    FiniteDiff* m_finiteDiffPtr;

    /// The current state and its descriptor.
    LevelData<FluxBox>*   m_velPtr;
    LevelData<FArrayBox>* m_pPtr;
    LevelData<FArrayBox>* m_qPtr;
    State*                m_statePtr;

    CFInterp*                m_cfInterpPtr;
    FluxRegisterFace*        m_velFluxRegPtr;
    AnisotropicFluxRegister* m_qFluxRegPtr;

    /// Stratification...
    bool m_hasStrat;

    /// All levels use these to interpolate to thier levels.
    /// We use a separate spline for the derivs to avoid introducing extrema.
    std::shared_ptr<CubicSpline> m_bbarSplinePtr;
    std::shared_ptr<CubicSpline> m_TbarSplinePtr;
    std::shared_ptr<CubicSpline> m_SbarSplinePtr;

    /// This level's stratification.
    std::shared_ptr<FArrayBox> m_bbarPtr;
    std::shared_ptr<FArrayBox> m_TbarPtr;
    std::shared_ptr<FArrayBox> m_SbarPtr;

    /// This level's mode-1 structure functions.
    Real                       m_c1;
    std::shared_ptr<FArrayBox> m_phi1Ptr;
    std::shared_ptr<FArrayBox> m_NsqPtr;  // = -dbdz


    /// Projector stuff...
    /// Use smart ptrs because we use validation functions, and I don't want
    /// to worry about invalidation functions.
    using PoissonOp       = Elliptic::PoissonOp;
    // using LevelProjSolver = Elliptic::LevelHybridSolver;
    // using AMRProjSolver   = Elliptic::AMRLepticSolver;
    using LevelProjSolver = Elliptic::MGSolver<LevelData<FArrayBox>>;
    using AMRProjSolver   = Elliptic::AMRMGSolver<LevelData<FArrayBox>>;

    std::shared_ptr<PoissonOp>       m_projOpPtr;
    std::shared_ptr<LevelProjSolver> m_levelProjSolverPtr;
    std::shared_ptr<AMRProjSolver>   m_amrProjSolverPtr;

    /// The RK method.
    // #/#/# = explicit/implicit/implicit stage order.
    // typedef PARK<ForwardEuler_Coeffs>   PARKType;  // 1/-/1
    // typedef PARK<Midpoint_Coeffs>   PARKType;
    // typedef PARK<CN_RKW3_Coeffs>        PARKType;  // 3/2A/1
    // typedef PARK<IMEXRKCB2_Coeffs>      PARKType;  // 2SSP/2L/1
    // typedef PARK<IMEXRKCB3c_Coeffs>     PARKType;  // 3SSP/3L/1 large ERK stability extent
    // typedef PARK<IMEXRKCB3f_Coeffs>     PARKType;  // 3/3L/2
    typedef PARK<ARK3_2_4L_2_SA_Coeffs> PARKType;  // 3/3L/2
    // typedef PARK<ARK4_3_6L_2_SA_Coeffs> PARKType;  // 4/3L/2  (Best for IB)
    // typedef PARK<RK4_Coeffs> PARKType;

    // typedef ForwardEuler PARKType;
    // typedef SDC PARKType;
    // typedef RKW3CN PARKType;

    PARKType *m_parkPtr;


    // These are used during regridding to save data when calling deactivate.
    std::shared_ptr<LevelData<FluxBox>>   m_oldVelPtr;
    std::shared_ptr<LevelData<FArrayBox>> m_oldPPtr;
    std::shared_ptr<LevelData<FArrayBox>> m_oldQPtr;
};


// -----------------------------------------------------------------------------
#ifndef NDEBUG

/// Print a "flowchart" to pout.
/// Basically just tells you the order major functions are called.
class __flowchart_auto
{
public:
    __flowchart_auto(bool a_doFlowchart, const char* a_fn, int a_level)
    : m_doFlowchart(a_doFlowchart)
    {
        if (m_doFlowchart) {
            pout() << "FLOWCHART: "
                   << a_fn << " [" << a_level << "]" << std::endl
                   << Format::indent();
        }
    }
    ~__flowchart_auto()
    {
        if (m_doFlowchart) {
            pout() << Format::unindent;
        }
    }
private:
    bool m_doFlowchart;
};

#define BEGIN_FLOWCHART()                                               \
    __flowchart_auto(ProblemContext::getInstance()->output.doFlowchart, \
                     __PRETTY_FUNCTION__,                               \
                     m_level)
#else
#define BEGIN_FLOWCHART()
#endif


#endif //!___AMRNSLevel_H__INCLUDED___
