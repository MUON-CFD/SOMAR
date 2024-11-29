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
 *  https://github.com/somarhub.
 ******************************************************************************/
#ifndef ___EllipticAMRMGOperator_H__INCLUDED___
#define ___EllipticAMRMGOperator_H__INCLUDED___

#include "MGOperator.H"

namespace Elliptic {


/**
 * @interface AMRMGOperator
 * @brief     The generic AMRMGOperator for use in all AMRMGSolvers.
 * @date      2020/04/30
 * @author    Edward.Santilli@jefferson.edu
 *
 * @details
 *  The generic AMRMGOperator for use in all AMRMGSolvers. This will act as
 *  its own factory. You create the op at (level,depth) = (0,0), and AMRMGSolver
 *  will do the rest by calling newAMRMGOperator().
 */
template <class T>
class AMRMGOperator: public MGOperator<T>
{
public:
    typedef T                                   StateType;
    typedef typename StateTraits<T>::CopierType CopierType;

    /// Apply the AMR operator, including coarse-fine matching.
    /// The CFI BCs are always inhomogeneous in this function.
    virtual void
    AMROperator(StateType&       a_LofPhi,
                StateType&       a_phiFine,  // You may set ghosts.
                StateType&       a_phi,      // You may set ghosts.
                const StateType& a_phiCoarse,
                const IntVect&   a_fineRefRatio,
                const IntVect&   a_crseRefRatio,
                const Real       a_time,
                const bool       a_homogPhysBCs,
                const AMRMGOperator<StateType>& a_finerOp) const = 0;

    /// Apply the AMR operator, including coarse-fine matching.
    /// The CFI BCs are always inhomogeneous in this function.
    /// Assume no finer AMR level.
    virtual void
    AMROperatorNF(StateType&       a_LofPhi,
                  StateType&       a_phi,  // You may set ghosts.
                  const StateType& a_phiCoarse,
                  const IntVect&   a_crseRefRatio,
                  const Real       a_time,
                  const bool       a_homogPhysBCs) const = 0;

    /// Apply the AMR operator, including coarse-fine matching.
    /// Assume no coarser AMR level.
    virtual void
    AMROperatorNC(StateType&     a_LofPhi,
                  StateType&     a_phiFine,  // You may set ghosts.
                  StateType&     a_phi,      // You may set ghosts.
                  const IntVect& a_fineRefRatio,
                  const Real     a_time,
                  const bool     a_homogPhysBCs,
                  const AMRMGOperator<StateType>& a_finerOp) const = 0;

    /// \brief
    ///  Compute the p-Norm of a_x only on where not covered by a_fineRes.
    ///  If a_p = 0, then we compute an inf-norm.
    /// \details
    ///  This function must perform all MPI communication.
    ///  We want the norm on this level to be comparable to the norms on
    ///  other levels, so we choose powScale = dx*dy*dz and compute the norm as
    ///   |a_res|_p = ( Sum |a_res|^p * powScale ) ^(1/p).
    ///
    ///  Sorry if the name is confusing. It means we want the norm of a single
    ///  level, but that level may be covered by another level in an AMR solve.
    ///  If a_fineResPtr is nullptr, then we just need plain ol' norm().
    virtual Real
    AMRNormLevel(const StateType& a_res,
                 const StateType* a_fineResPtr,
                 const IntVect&   a_refRatio,
                 const int        a_p) const = 0;

    // -------------------------------------------------------------------------
    // Optional optimizations for an AMRMGOperator.  These are not pure virtual
    // functions The AMRMGSolver algorithm calls these functions, which a smart
    // operator might be able to perform in faster ways.
    // -------------------------------------------------------------------------

    /// Uses AMROperator to compute a_lhs = a_rhs - L[a_phi].
    virtual void
    AMRResidual(StateType&       a_lhs,
                StateType&       a_phiFine,  // You may set ghosts.
                StateType&       a_phi,      // You may set ghosts.
                const StateType& a_phiCoarse,
                const StateType& a_rhs,
                const IntVect&   a_fineRefRatio,
                const IntVect&   a_crseRefRatio,
                const Real       a_time,
                const bool       a_homogPhysBCs,
                const AMRMGOperator<StateType>& a_finerOp) const

    {
        this->AMROperator(a_lhs,
                          a_phiFine,
                          a_phi,
                          a_phiCoarse,
                          a_fineRefRatio,
                          a_crseRefRatio,
                          a_time,
                          a_homogPhysBCs,
                          a_finerOp);
        this->scale(a_lhs, -1.0);
        this->incr(a_lhs, a_rhs, 1.0);
    }

    /// Uses AMROperatorNF to compute a_lhs = a_rhs - L[a_phi].
    virtual void
    AMRResidualNF(StateType&       a_lhs,
                  StateType&       a_phi,  // You may set ghosts.
                  const StateType& a_phiCoarse,
                  const StateType& a_rhs,
                  const IntVect&   a_crseRefRatio,
                  const Real       a_time,
                  const bool       a_homogPhysBCs) const
    {
        this->AMROperatorNF(a_lhs,
                            a_phi,
                            a_phiCoarse,
                            a_crseRefRatio,
                            a_time,
                            a_homogPhysBCs);
        this->scale(a_lhs, -1.0);
        this->incr(a_lhs, a_rhs, 1.0);
    }

    /// Uses AMROperatorNC to compute a_lhs = a_rhs - L[a_phi].
    virtual void
    AMRResidualNC(StateType&       a_lhs,
                  StateType&       a_phiFine,  // You may set ghosts.
                  StateType&       a_phi,      // You may set ghosts.
                  const StateType& a_rhs,
                  const IntVect&   a_fineRefRatio,
                  const Real       a_time,
                  const bool       a_homogPhysBCs,
                  const AMRMGOperator<StateType>& a_finerOp) const
    {
        this->AMROperatorNC(a_lhs,
                            a_phiFine,
                            a_phi,
                            a_fineRefRatio,
                            a_time,
                            a_homogPhysBCs,
                            a_finerOp);
        this->scale(a_lhs, -1.0);
        this->incr(a_lhs, a_rhs, 1.0);
    }

    /// Uses AMROperatorNF to compute a_res = a_res - L[a_phi].
    /// a_scratchRes is scratch space defined to hold a_res.
    virtual void
    AMRUpdateResidualNF(StateType&       a_res,
                        StateType&       a_phi,  // You may set ghosts.
                        const StateType& a_phiCoarse,
                        StateType&       a_scratchRes,
                        const IntVect&   a_crseRefRatio,
                        const Real       a_time,
                        const bool       a_homogPhysBCs) const
    {
        this->assignLocal(a_scratchRes, a_res);
        this->AMRResidualNF(a_res,
                            a_phi,
                            a_phiCoarse,
                            a_scratchRes,
                            a_crseRefRatio,
                            a_time,
                            a_homogPhysBCs);
    }

    // -------------------------------------------------------------------------
    // Functions for ops with high-order of spatial accuracy.
    // -------------------------------------------------------------------------

    /// Order of spatial accuracy of this op.
    /// Default = 2.
    virtual unsigned int
    orderOfAccuracy() const
    {
        return 2;
    };

    /// This is needed vor operators with orderOfAccuracy() > 2.
    /// Default is no-op.
    virtual void
    enforceCFConsistency(T&       /*a_coarseCorrection*/,
                         const T& /*a_correction*/) const
    {
        // Do nothing.
    }
};


}; // namespace Elliptic
#endif //!___EllipticAMRMGOperator_H__INCLUDED___
