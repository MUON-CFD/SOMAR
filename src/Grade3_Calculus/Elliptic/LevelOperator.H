/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2019
 *    Jefferson University and
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
#ifndef ___EllipticLevelOperator_H__INCLUDED___
#define ___EllipticLevelOperator_H__INCLUDED___

#include "StateOps.H"
#include "REAL.H"

namespace Elliptic {


/**
 * @interface LevelOperator
 * @brief     The generic operator with no MG or AMR capabilities.
 * @date      2019/11/26
 * @author    Edward.Santilli@jefferson.edu
 *
 * @details
 *  The generic LevelOperator. This can be used in LevelSolvers, but not
 *  MG or AMR solvers.
 */
template <class T>
class LevelOperator: public StateOps<T>
{
public:
    typedef T StateType;

    virtual
    ~LevelOperator() {}

    /// Pure virtual.
    /// Sets BCs and fills all ghosts of a_phi.
    /// This is needed by the MG and AMR prolongation functions, but we include
    /// it here because it is likely to be used in applyOp. If you are not using
    /// ghosts or if you are not creating an MG or AMR operator, this can be a
    /// simple no-op. a_crsePhiPtr can be nullptr if not needed.
    virtual void
    applyBCs(StateType&       a_phi,
             const StateType* a_crsePhiPtr,
             const Real       a_time,
             const bool       a_homogPhysBCs,
             const bool       a_homogCFIBCs) const = 0;

    /// Pure virtual.
    /// Computes lhs = L[phi].
    virtual void
    applyOp(StateType&       a_lhs,
            StateType&       a_phi,
            const StateType* a_crsePhiPtr,
            const Real       a_time,
            const bool       a_homogPhysBCs,
            const bool       a_homogCFIBCs) const = 0;

    /// Pure virtual
    /// If this operator has a null space, a_phi needs to be projected.
    /// In other words, if L[e^1] = 0, and phi = c_1*e^1 + c_2*e^2 + ...,
    /// then this function should remove c_1*e^1.
    /// (Typically, you would just remove the average.)
    virtual void
    removeKernel(StateType& a_phi) const = 0;

    /// Pure virtual.
    /// Sets phi to a preconditioned value.
    /// Typically, this will be phi = rhs / invDiags, then relaxed, but you
    /// can ignore a_relaxIters if you have a better idea.
    virtual void
    preCond(StateType&       a_phi,
            const StateType& a_rhs,
            const Real       a_time,
            const int        a_relaxIters) const = 0;

    /// Pure virtual
    /// Applies relaxation to the residual equation.
    /// If you don't plan to use this in preCond, or if you aren't creating
    /// an MG or AMR op, then this function can be a simple no-op.
    virtual void
    relax(StateType&       a_cor,
          const StateType& a_res,
          const Real       a_time,
          const int        a_relaxIters) const = 0;

    /// Computes lhs = rhs - L[phi].
    /// You can override the default functionality if you have a better idea.
    virtual void
    residual(StateType&       a_lhs,
             StateType&       a_phi,
             const StateType* a_crsePhiPtr,
             const StateType& a_rhs,
             const Real       a_time,
             const bool       a_homogPhysBCs,
             const bool       a_homogCFIBCs) const
    {
        this->applyOp(
            a_lhs, a_phi, a_crsePhiPtr, a_time, a_homogPhysBCs, a_homogCFIBCs);
        this->scale(a_lhs, -1.0);
        this->incr(a_lhs, a_rhs, 1.0);
    }

    /// Checks if L[phi] = rhs is a solvable problem.
    /// By default, this returns true. Override if you have a better plan.
    virtual bool
    levelEquationIsConsistent(const StateType& /*a_phi*/,
                              const StateType* /*a_crsePhiPtr*/,
                              const StateType& /*a_rhs*/,
                              const Real       /*a_time*/,
                              const bool       /*a_homogPhysBCs*/,
                              const bool       /*a_homogCFIBCs*/) const
    {
        return true;
    }

protected:
    LevelOperator() = default;
    LevelOperator(const LevelOperator&) = delete;
    LevelOperator& operator=(const LevelOperator&) = delete;
};


}; // namespace Elliptic
#endif //!___EllipticLevelOperator_H__INCLUDED___
