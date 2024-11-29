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
 *  https://github.com/MUON-CFD/somar.
 ******************************************************************************/
#ifndef ___EllipticLevelSolver_H__INCLUDED___
#define ___EllipticLevelSolver_H__INCLUDED___

#include "LevelOperator.H"
#include "SolverStatus.H"
#include "CH_assert.H"
#include <memory>

namespace Elliptic {


/**
 * @interface LevelSolver
 * @brief     The generic LevelSolver interface.
 * @date      2019/11/26
 * @author    Edward.Santilli@jefferson.edu
 * @tparam    StateType     The type of the state variable given to this solver.
 */
template <class StateType>
class LevelSolver
{
public:
    /// Default constructor
    LevelSolver();

    /// Virtual destructor
    virtual ~LevelSolver();

    /// Full virtual constructor
    virtual void
    define (std::shared_ptr<const LevelOperator<StateType>> a_opPtr);

    /// Retrieves the operator used during the solve.
    virtual const LevelOperator<StateType>&
    getOp() const
    {
        CH_assert(m_opPtr);
        return *m_opPtr;
    }

    /// Returns a copy of our exit status.
    virtual SolverStatus
    getSolverStatus() const
    {
        return m_solverStatus;
    }

    /**
     * @brief        Solves L[a_phi] = a_rhs, where L is given by *m_opPtr.
     *
     * @param[inout] a_phi                The initial gues on entry.
     *                                    The solution on exit.
     * @param[in]    a_crsePhiPtr         Used to set CFI BCs, if needed.
     *                                    Can be null if not needed.
     * @param[in]    a_rhs
     * @param[in]    a_time               Potentially needed by the BCs.
     * @param[in]    a_useHomogBCs
     * @param[in]    a_setPhiToZero       Should we initialize the initial
     *                                    guess to zero?
     * @param[in]    a_convergenceMetric  If suppplied, all residual norms will
     *                                    be divided by this.
     *
     * @return     The exit status. (Did we converge?)
     *
     * @detail
     *  If a_convergenceMetric is provided, it will take precedence over
     *  Options::convergenceMetric.
     */
    virtual SolverStatus
    solve(StateType&       a_phi,
          const StateType* a_crsePhiPtr,
          const StateType& a_rhs,
          const Real       a_time,
          const bool       a_useHomogBCs,
          const bool       a_setPhiToZero,
          const Real       a_convergenceMetric = -1.0) const = 0;

protected:
    /// References our exit status.
    virtual SolverStatus&
    getSolverStatusRef() const
    {
        return m_solverStatus;
    }

    std::shared_ptr<const LevelOperator<StateType>> m_opPtr;
    mutable SolverStatus                            m_solverStatus;
};


/**
 * @interface NullSolver
 * @brief     This does nothing! It's useful for debugging V-Cycles.
 * @date      2019/11/26
 * @author    Edward.Santilli@jefferson.edu
 * @tparam    StateType     The type of the state variable given to this solver.
 */
template<class StateType>
class NullSolver: public LevelSolver<StateType>
{
public:
    ///
    NullSolver()
    {
    }

    ///
    virtual ~NullSolver()
    {
    }

    ///
    virtual SolverStatus
    solve(StateType&       a_phi,
          const StateType* a_crsePhiPtr,
          const StateType& a_rhs,
          const Real       a_time,
          const bool       a_useHomogBCs,
          const bool       a_setPhiToZero,
          const Real       a_convergenceMetric = -1.0) const
    {
        return this->getSolverStatus();
    }
};


/**
 * @interface BiCGStabSolver
 * @brief     Stabilized Bi-Conjugate Gradient solver.
 *            Copied almost directly from Chombo.
 * @date      2019/11/26
 * @author    Edward.Santilli@jefferson.edu
 *
 * @tparam    StateType The type of the state variable given to this solver.
 */
template<class StateType>
class RelaxSolver: public LevelSolver<StateType>
{
public:
    /// Knobs for you to turn.
    struct Options {
        Options() {}  // A bug in the C++11 standard requires this.

        Real absTol = 1.0e-6;  // Stop if |res| < absTol.
        Real relTol = 1.0e-6;  // Stop if |res| / |init res| < relTol.
        Real hang   = 1.0e-7;
        Real convergenceMetric = -1.0;  // If < 0, use initial residual norm.
        int  minIters          = 1;
        int  maxIters          = 100;
        int  itersPerNormCheck = 10;
        int  normType          = 2;  // 0 = inf-norm, etc...
        int  verbosity         = 4;
    };

    /// Constructor. You must still call define() to supply an operator.
    RelaxSolver(Options a_opt = Options());

    /// Destructor
    virtual ~RelaxSolver();

    /// Retrieves the solver's options.
    virtual const Options&
    getOptions() const
    {
        return m_opt;
    }

    /// Replace this solver's options.
    virtual void
    setOptions(const Options& a_opt)
    {
        m_opt = a_opt;
    }

    /// This allows you to manipulate individual options.
    virtual Options&
    options()
    {
        return m_opt;
    }

    /**
     * @brief        Solves L[a_phi] = a_rhs, where L is given by *m_opPtr.
     *
     * @param[inout] a_phi                The initial gues on entry.
     *                                    The solution on exit.
     * @param[in]    a_rhs
     * @param[in]    a_time               Potentially needed by the BCs.
     * @param[in]    a_useHomogBCs
     * @param[in]    a_setPhiToZero       Should we initialize the initial
     *                                    guess to zero?
     * @param[in]    a_convergenceMetric  If suppplied, all residual norms will
     *                                    be divided by this.
     *
     * @return     The exit status. (Did we converge?)
     *
     * @detail
     *  If a_convergenceMetric is provided, it will take precedence over
     *  Options::convergenceMetric.
     */
    virtual SolverStatus
    solve(StateType&       a_phi,
          const StateType* a_crsePhiPtr,
          const StateType& a_rhs,
          const Real       a_time,
          const bool       a_useHomogBCs,
          const bool       a_setPhiToZero,
          const Real       a_convergenceMetric = -1.0) const;

protected:
    Options m_opt;
};


/**
 * @interface BiCGStabSolver
 * @brief     Stabilized Bi-Conjugate Gradient solver.
 *            Copied almost directly from Chombo.
 * @date      2019/11/26
 * @author    Edward.Santilli@jefferson.edu
 *
 * @tparam    StateType The type of the state variable given to this solver.
 */
template<class StateType>
class BiCGStabSolver: public LevelSolver<StateType>
{
public:
    /// Knobs for you to turn.
    struct Options {
        Options() {}  // A bug in the C++11 standard requires this.

        Real absTol = 1.0e-6;  // Stop if |res| < absTol.
        Real relTol = 1.0e-6;  // Stop if |res| / |init res| < relTol.
        Real small  = 1.0e-30;
        Real hang   = 1.0e-7;
        Real convergenceMetric = -1.0;  // If < 0, use initial residual norm.
        int  maxIters          = 80;
        int  maxRestarts       = 5;
        int  normType          = 2;  // 0 = inf-norm, etc...
        int  verbosity         = 3;
        int  numSmoothPrecond  = 2;
    };

    /// Constructor. You must still call define() to supply an operator.
    BiCGStabSolver(Options a_opt = Options());

    /// Destructor
    virtual ~BiCGStabSolver();

    /// Retrieves the solver's options.
    virtual const Options&
    getOptions() const
    {
        return m_opt;
    }

    /// Replace this solver's options.
    virtual void
    setOptions(const Options& a_opt)
    {
        m_opt = a_opt;
    }

    /// This allows you to manipulate individual options.
    virtual Options&
    options()
    {
        return m_opt;
    }

    /**
     * @brief        Solves L[a_phi] = a_rhs, where L is given by *m_opPtr.
     *
     * @param[inout] a_phi                The initial gues on entry.
     *                                    The solution on exit.
     * @param[in]    a_rhs
     * @param[in]    a_time               Potentially needed by the BCs.
     * @param[in]    a_useHomogBCs
     * @param[in]    a_setPhiToZero       Should we initialize the initial
     *                                    guess to zero?
     * @param[in]    a_convergenceMetric  If suppplied, all residual norms will
     *                                    be divided by this.
     *
     * @return     The exit status. (Did we converge?)
     *
     * @detail
     *  If a_convergenceMetric is provided, it will take precedence over
     *  Options::convergenceMetric.
     */
    virtual SolverStatus
    solve(StateType&       a_phi,
          const StateType* a_crsePhiPtr,
          const StateType& a_rhs,
          const Real       a_time,
          const bool       a_useHomogBCs,
          const bool       a_setPhiToZero,
          const Real       a_convergenceMetric = -1.0) const;

protected:
    Options m_opt;
};


}; // namespace Elliptic

#define H9c9778ae61b21da3b93db6cca54d9de6
#include "LevelSolverI.H"
#undef H9c9778ae61b21da3b93db6cca54d9de6

#endif //!___EllipticLevelSolver_H__INCLUDED___
