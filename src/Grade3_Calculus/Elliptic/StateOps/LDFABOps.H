/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2020
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
#ifndef ___EllipticLDFABOps_H__INCLUDED___
#define ___EllipticLDFABOps_H__INCLUDED___

#include "StateOps.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "Copier.H"

namespace Elliptic {


// Traits
template<>
struct StateTraits<LevelData<FArrayBox>>
{
    typedef DisjointBoxLayout GridsType;
    typedef Copier            CopierType;
};


// Policies
template<>
class StateOps< LevelData<FArrayBox>, StateTraits<LevelData<FArrayBox>> >
{
public:
    typedef LevelData<FArrayBox>                                   StateType;
    typedef typename StateTraits<LevelData<FArrayBox>>::CopierType CopierType;
    typedef typename StateTraits<LevelData<FArrayBox>>::GridsType  GridsType;

    virtual ~StateOps() {}

    // -------------------------------------------------------------------------
    /// \{
    /// \name {Needed by LevelOperators}

    /// Create data holder a_lhs that mirrors a_rhs. You do not need to copy
    /// the data of a_rhs, just make a holder the same size.
    virtual void
    create(StateType& a_lhs, const StateType& a_rhs) const;

    /// Opposite of create -- perform any operations required before lhs goes
    /// out of scope. In general, the only time this needs to be defined in
    /// a derived class is if the create() function called new. Otherwise, the
    /// default no-op function is sufficient.
    virtual void
    clear(StateType& a_lhs) const;

    /// Set a_lhs equal to a_rhs.
    /// If StateType does not have copiers, just send in a_copierPtr = nullptr.
    virtual void
    assign(StateType&        a_lhs,
           const StateType&  a_rhs,
           const CopierType* a_copierPtr = nullptr) const;

    /// Set a_lhs equal to a_rhs. Non-blocking version.
    /// a_lhs and a_rhs MUST be defined on compatible grids.
    virtual void
    assignLocal(StateType& a_lhs, const StateType& a_rhs) const;

    /// Compute and return the dot product of a_1 and a_2. In most contexts,
    /// this means return the sum over all data points of a_1*a_2.
    virtual Real
    dotProduct(const StateType& a_1, const StateType& a_2) const;

    /// \brief
    ///  Compute the p-Norm of a_x. If a_p = 0, then we compute an inf-norm.
    /// \details
    ///  This function must perform all MPI communication.
    ///  If you want the norm on this level to be comparable to the norms on
    ///  other levels, then a_powScale must be chosen carefully. For example,
    ///  with cell-averaged data, a_powScale = dx * dy * dz would work.
    ///  It is used like this: ( Sum[|a_x|^p] * a_powScale )^(1/p).
    virtual Real
    norm(const StateType& a_x,
         const int        a_p,
         const Real       a_powScale = 1.0) const;

    /// Increment by scaled amount (a_lhs += a_scale*a_x).
    virtual void
    incr(StateType& a_lhs, const StateType& a_x, Real a_scale) const;

    /// Set input to a scaled sum (a_lhs = a_a*a_x + a_b*a_y).
    virtual void
    axby(StateType&       a_lhs,
         const StateType& a_x,
         const StateType& a_y,
         Real             a_a,
         Real             a_b) const;

    /// Multiply the input by a given scale (a_lhs *= a_scale).
    virtual void
    scale(StateType& a_lhs, const Real& a_scale) const;

    /// Set a_lhs to zero.
    virtual void
    setToZero(StateType& a_lhs) const;
    /// \}

    // -------------------------------------------------------------------------
    /// \{
    /// \name {Needed by MGOperators}

    /// Create a coarsened version of a_fine.
    /// You do not need to fill a_crse with data, just define it properly.
    virtual void
    createCoarsened(StateType&       a_crse,
                    const StateType& a_fine,
                    const IntVect&   a_refRatio) const;
    /// \}


    // -------------------------------------------------------------------------
    /// \{
    /// \name {Needed by AMRMGOperators}

    /// Create a copier for a_lhs <-- a_rhs
    /// If your StateType does not support copiers, leave the function empty.
    virtual void
    buildCopier(CopierType&      a_copier,
                const StateType& a_lhs,
                const StateType& a_rhs) const;

    /// Create a copier for a_lhs --> a_rhs.
    /// If your StateType does not support copiers, leave the function empty.
    virtual void
    buildReverseCopier(CopierType&       a_reverseCopier,
                       const CopierType& a_copier,
                       const StateType&  a_lhs,
                       const StateType&  a_rhs) const;

    // /// Set lhs = 0 where covered by rhs.
    // /// rhs can be clobbered, it is temporary storage.
    // /// If you don't have a copier for any reason, just send in nullptr.
    // virtual void
    // zeroCovered(StateType&        a_lhs,
    //             StateType&        a_rhs,
    //             const CopierType* a_copierPtr = nullptr) const;

    // /// Compute the max norm. Don't worry about MPI communication, that is left
    // /// to the caller. This function should be non-blocking.
    // virtual Real
    // localMaxNorm(const StateType& a_x) const;
    /// \}
};



}; // namespace Elliptic
#endif //___EllipticLDFABOps_H__INCLUDED___
