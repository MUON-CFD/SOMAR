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
#ifndef ___AMRHybridSolver_H__INCLUDED___
#define ___AMRHybridSolver_H__INCLUDED___

#include "LevelData.H"
#include "FArrayBox.H"
#include "LevelHybridSolver.H"
#include "AMRMGOperator.H"


namespace Elliptic {


/// @brief The AMR version of the leptic solver.
/// @details
/// You are looking at the prototype of the AMR leptic method. This solver has
/// been working well for us, but the class itself is not 100% user-friendly.
/// There are a lot of parameters to set -- many of which are hard coded. In the
/// next version of this code, we hope to create a more versatile solver that is
/// a bit more automatic. That is, if you just need to solve equations without
/// getting under the hood, we should be able to provide a leptic solver with
/// less knobs to turn.
class AMRHybridSolver
{
public:
    typedef LevelData<FArrayBox>                        StateType;
    typedef AMRMGOperator<StateType>                    AMRMGOpType;
    typedef typename StateTraits<StateType>::CopierType CopierType;

    // -------------------------------------------------------------------------
    /// \name Option setters / getters
    /// \{

    /// Knobs for you to turn.
    struct Options {
        /// A bug in the C++11 standard requires this,
        /// even though it does nothing.
        Options() {}

        /// Stop if |res| < absTol
        Real absTol            = -1.0;
        /// Stop if |res| / |init res| < relTol
        Real relTol            = -1.0;
        Real convergenceMetric = -1.0;
        int  numCycles         = 0;
        int  maxIters          = -1;
        Real hang              = -1;

        /// @brief 0 = inf-norm, 1 = 1-norm, etc...
        /// @details
        /// I suggest using the 2-norm (default) so that norms across different
        /// MG depths are comparable.
        int  normType          = -1;
        int  verbosity         = -1;
    };

    /// @brief Creates an Options struct set to defaults as defined by the input
    /// file.
    static Options
    getDefaultOptions();

    /// Retrns the current solver's options.
    inline const Options&
    getOptions() const {
        return m_opt;
    }

    /// @brief Gets the current {absTol, relTol}.
    inline std::array<Real, 2>
    getTols() const {
        return {m_opt.absTol, m_opt.relTol};
    }

    /// @brief Sets the new {absTol, relTol}.
    /// @param a_tols
    inline void
    setTols(const std::array<Real, 2>& a_tols) {
        m_opt.absTol = a_tols[0];
        m_opt.relTol = a_tols[1];
    }
    /// \}


    // -------------------------------------------------------------------------
    /// \name Constructors / destructors
    /// \{

    /// Default constructor. Leaves object unusable.
    AMRHybridSolver ();

    /// Destructor
    ~AMRHybridSolver ();

    /// Set up this solver to solve L[phi] = rhs on levels a_lmin to a_lmax.
    /// If a_lmin > 0, then a_vAMRMGOps must also be defined on a_lmin - 1.
    /// a_vAMRMGOps can be defined or undefined on the levels we don't use.
    void
    define(Vector<std::shared_ptr<const AMRMGOpType>> a_vAMRMGOps,
           const size_t                               a_lmin,
           const size_t                               a_lmax,
           Options                                    a_opt);

    // Is this solver ready for use?
    inline bool
    isDefined() const
    {
        return m_isDefined;
    }

    /// Free's memory. Leaves object unusable.
    void
    clear();
    /// \}


    // -------------------------------------------------------------------------
    /// \name Solvers
    /// \{

    /**
     * @brief      Solves L[phi] = rhs using the most appropriate methods.
     *
     * @param[inout] a_phi                On entry, the initial guess.
     *                                    On exit, the solution.
     * @param[in]    a_rhs
     * @param[in]    a_time               Used to set BCs.
     * @param[in]    a_useHomogBCs
     * @param[in]    a_setPhiToZero       Should we set initial guess to zero?
     * @param[in]    a_convergenceMetric  If supplied, this will be used as the
     *                                    initial residual when computing
     *                                    relative residuals. This setting takes
     *                                    precedence over
     *                                    Options::convergenceMetric.
     *
     * @return     The SolverStatus object.
     *             This tells the caller if we converged.
     */
    SolverStatus
    solve(Vector<StateType*>&             a_vphi,
          const Vector<const StateType*>& a_vrhs,
          const Real                      a_time,
          const bool                      a_useHomogBCs       = false,
          const bool                      a_setPhiToZero      = false,
          const Real                      a_convergenceMetric = -1.0) const;
    /// \}

    /// Returns the exit status.
    inline SolverStatus
    getSolverStatus() const
    {
        return m_solverStatus;
    }

protected:
    /// Relaxes L[e] = r over the AMR hierarchy with a V-Cycle.
    /// It is assumed the problem is homogeneous, so send in the residual eq.
    /// a_vphi needs to be set to zero before entry.
    void
    amrVCycle_residualEq(Vector<unique_ptr<StateType>>& a_e,
                         Vector<unique_ptr<StateType>>& a_r,
                         const Real                     a_time,
                         const size_t                   a_lev) const;

    /// A solve requires a number of temporaries over the entire AMR hierarchy.
    /// Allocation + definition is potentially expensive, so I shoved it all
    /// into this function. If the expense is high, we'll see it in our
    /// profilers.
    void
    initializeSolve(const Vector<StateType*>&       a_vphi,
                    const Vector<const StateType*>& a_vrhs) const;

    /// Deallocation of temporaries.
    void
    finalizeSolve() const;

    /// Computes lhs = rhs - L[phi] from m_minSolveLevel to m_maxSolveLevel.
    /// This calls computeAMRResidualLevel on each level and computes the
    /// inf-norm.
    Real
    computeAMRResidual(Vector<unique_ptr<StateType>>&  a_vres,
                       Vector<StateType*>&             a_vphi,
                       const Vector<const StateType*>& a_vrhs,
                       const Real                      a_time,
                       const bool                      a_useHomogBCs,
                       const bool                      a_computeNorm) const;

    /// Computes the residual at a single level using the appropriate version
    /// of AMROperator.
    /// a_finePhiPtr ghosts may be reset.
    void
    computeAMRResidualLevel(StateType&       a_res,
                            StateType*       a_finePhiPtr,
                            StateType&       a_phi,
                            const StateType* a_crsePhiPtr,
                            const StateType& a_rhs,
                            const size_t     a_lev,
                            const Real       a_time,
                            const bool       a_useHomogBCs) const;

    /// Eliminate high wavenumber errors on a_lev that will not be captured
    /// by the grids on a_lev - 1. This either calls the MG relax function
    /// or performs a mini V-cycle.
    void
    smoothDown(StateType&   a_cor,
               StateType&   a_res,
               const Real   a_time,
               const size_t a_lev) const;

    /// Eliminate high wavenumber errors on a_lev that were not captured by
    /// the grids on a_lev - 1. This either calls the MG relax function or
    /// performs a mini V-cycle.
    void
    smoothUp(StateType&   a_cor,
             StateType&   a_res,
             const Real   a_time,
             const size_t a_lev) const;

    /// pout formatting utility.
    void
    indent(const int a_lev = -1) const;

    /// pout formatting utility.
    void
    unindent(const int a_lev = -1) const;

    /// References our exit status.
    SolverStatus&
    getSolverStatusRef() const
    {
        return m_solverStatus;
    }

    bool                   m_isDefined;
    Options                m_opt;
    mutable SolverStatus   m_solverStatus;

    // All AMR solves must happen between these levels.
    size_t m_lbase;
    size_t m_lmin;
    size_t m_lmax;

    // These Vector indices are level numbers.
    Vector<std::shared_ptr<const AMRMGOpType>> m_vAMRMGOps;
    Vector<IntVect>                            m_vAMRCrseRefRatios;
    Vector<IntVect>                            m_vAMRFineRefRatios;
    Vector<std::unique_ptr<LevelHybridSolver>> m_vHybridSolver;

    mutable Vector<unique_ptr<StateType>>  m_ve;    // Defined on lbase too!
    mutable Vector<unique_ptr<StateType>>  m_vr;
    mutable Vector<unique_ptr<StateType>>  m_vScratchPhi;
    mutable Vector<unique_ptr<StateType>>  m_vResC;          // coarsened rhs[l]
    mutable Vector<unique_ptr<CopierType>> m_vResCopier;     // from ResC[l] to Res[l-1]
    mutable Vector<unique_ptr<CopierType>> m_vReverseCopier; // from Res[l-1] to ResC[l]

private:
    AMRHybridSolver(const AMRHybridSolver&) = delete;
    AMRHybridSolver& operator=(const AMRHybridSolver&) = delete;
};



}; // namespace Elliptic
#endif //!___AMRHybridSolver_H__INCLUDED___
