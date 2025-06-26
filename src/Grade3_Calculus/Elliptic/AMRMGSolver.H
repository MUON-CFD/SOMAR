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
#ifndef ___EllipticAMRMGSolver_H__INCLUDED___
#define ___EllipticAMRMGSolver_H__INCLUDED___

#include "MGSolver.H"
#include "AMRMGOperator.H"


namespace Elliptic {


/**
 * @interface AMRMGSolver
 * @brief     The AMR multi-grid elliptic solver.
 * @date      2020/04/30
 * @author    Edward.Santilli@jefferson.edu
 * @note
 *  Currently, this only works when solving L[phi]=rhs where phi and rhs are
 *  StateTypees.
 */
template <class T>
class AMRMGSolver
{
public:
    typedef T                                   StateType;
    typedef MGSolver<T>                         MGSolverType;
    typedef AMRMGOperator<T>                    AMRMGOpType;
    typedef typename StateTraits<T>::CopierType CopierType;

    /// Knobs for you to turn.
    struct Options {
        /// A bug in the C++11 standard requires this,
        /// even though it does nothing.
        Options() {}

        /// Stop if |res| < absTol
        Real absTol            = 1.0e-12;
        /// Stop if |res| / |init res| < relTol
        Real relTol            = 1.0e-12;
        Real convergenceMetric = -1.0;
        int  numCycles         = 1;
        int  maxIters          = 10;
        Real hang              = 0.1;
        // /**
        //  * @brief 0 = inf-norm, 1 = 1-norm, etc...
        //  * @detail
        //  *  I suggest using the 2-norm (defualt) so that norms across different
        //  *  MG depths are comparable.
        //  */
        int  normType          = 2;
        int  verbosity         = 4;

        typename MGSolverType::Options mgOptions;
    };

    /// Default constructor. Leaves object unusable.
    AMRMGSolver ();

    /// Destructor
    virtual ~AMRMGSolver ();

    /// Set up this solver to solve L[phi] = rhs on levels a_lmin to a_lmax.
    /// If a_lmin > 0, then a_vAMRMGOps must also be defined on a_lmin - 1.
    /// a_vAMRMGOps can be defined or undefined on the levels we don't use.
    virtual void
    define(Vector<std::shared_ptr<const AMRMGOpType>> a_vAMRMGOps,
           const size_t                               a_lmin,
           const size_t                               a_lmax,
           Options                                    a_opt);

    /// Set up this solver to solve L[phi] = rhs on levels a_lmin to a_lmax.
    /// If a_lmin > 0, then a_vAMRMGOps must also be defined on a_lmin - 1.
    /// a_vAMRMGOps can be defined or undefined on the levels we don't use.
    /// This version does not set Options, just in case you already did that.
    virtual void
    define(Vector<std::shared_ptr<const AMRMGOpType>> a_vAMRMGOps,
           const size_t                               a_lmin,
           const size_t                               a_lmax);

    // Is this solver ready for use?
    virtual bool
    isDefined() const
    {
        return m_isDefined;
    }

    /// Free's memory. Leaves object unusable.
    virtual void
    clear();

    /// Returns the exit status.
    virtual SolverStatus
    getSolverStatus() const
    {
        return m_solverStatus;
    }

    /// \brief Initializes all parameters with default values.
    virtual void
    setDefaultOptions();

    /// Returns the current solver's options.
    virtual const Options&
    getOptions() const
    {
        return m_opt;
    }

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
    virtual SolverStatus
    solve(Vector<StateType*>&             a_vphi,
          const Vector<const StateType*>& a_vrhs,
          const Real                      a_time,
          const bool                      a_useHomogBCs       = false,
          const bool                      a_setPhiToZero      = false,
          const Real                      a_convergenceMetric = -1.0) const;

protected:
    /// Relaxes L[e] = r over the AMR hierarchy with a V-Cycle.
    /// It is assumed the problem is homogeneous, so send in the residual eq.
    /// a_vphi needs to be set to zero before entry.
    virtual void
    amrVCycle_residualEq(Vector<unique_ptr<StateType>>& a_e,
                         Vector<unique_ptr<StateType>>& a_r,
                         const Real                     a_time,
                         const size_t                   a_lev) const;

    /// A solve requires a number of temporaries over the entire AMR hierarchy.
    /// Allocation + definition is potentially expensive, so I shoved it all
    /// into this function. If the expense is high, we'll see it in our
    /// profilers.
    virtual void
    initializeSolve(const Vector<StateType*>&       a_vphi,
                    const Vector<const StateType*>& a_vrhs) const;

    /// Deallocation of temporaries.
    virtual void
    finalizeSolve() const;

    /// Computes lhs = rhs - L[phi] from m_minSolveLevel to m_maxSolveLevel.
    /// This calls computeAMRResidualLevel on each level and computes the
    /// inf-norm.
    virtual Real
    computeAMRResidual(Vector<unique_ptr<StateType>>&  a_vres,
                       Vector<StateType*>&             a_vphi,
                       const Vector<const StateType*>& a_vrhs,
                       const Real                      a_time,
                       const bool                      a_useHomogBCs,
                       const bool                      a_computeNorm) const;

    /// Computes the residual at a single level using the appropriate version
    /// of AMROperator.
    /// a_finePhiPtr ghosts may be reset.
    virtual void
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
    virtual void
    smoothDown(StateType&   a_cor,
               StateType&   a_res,
               const Real   a_time,
               const size_t a_lev) const;

    /// Eliminate high wavenumber errors on a_lev that were not captured by
    /// the grids on a_lev - 1. This either calls the MG relax function or
    /// performs a mini V-cycle.
    virtual void
    smoothUp(StateType&   a_cor,
             StateType&   a_res,
             const Real   a_time,
             const size_t a_lev) const;

    /// pout formatting utility.
    virtual void
    indent(const int a_lev = -1) const;

    /// pout formatting utility.
    virtual void
    unindent(const int a_lev = -1) const;

    /// References our exit status.
    virtual SolverStatus&
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
    Vector<IntVect>                            m_vAMRCrseRefRatios;
    Vector<IntVect>                            m_vAMRFineRefRatios;
    Vector<std::shared_ptr<const AMRMGOpType>> m_vAMRMGOps;
    Vector<std::unique_ptr<MGSolverType>>      m_vMiniVCycleSolvers;
    std::unique_ptr<MGSolverType>              m_bottomMGSolver;

    mutable Vector<unique_ptr<StateType>>  m_ve;    // Defined on lbase too!
    mutable Vector<unique_ptr<StateType>>  m_vr;
    mutable Vector<unique_ptr<StateType>>  m_vScratchPhi;
    mutable Vector<unique_ptr<StateType>>  m_vResC;          // coarsened rhs[l]
    mutable Vector<unique_ptr<CopierType>> m_vResCopier;     // from ResC[l] to Res[l-1]
    mutable Vector<unique_ptr<CopierType>> m_vReverseCopier; // from Res[l-1] to ResC[l]

private:
    AMRMGSolver(const AMRMGSolver&) = delete;
    AMRMGSolver& operator=(const AMRMGSolver&) = delete;
};


}; // namespace Elliptic

#define Hb85c5095bf0dc1fbde78ea67277de9cb
#include "AMRMGSolverI.H"
#undef Hb85c5095bf0dc1fbde78ea67277de9cb

#endif //!___EllipticAMRMGSolver_H__INCLUDED___
