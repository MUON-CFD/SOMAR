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
#ifndef ___EllipticMGOperator_H__INCLUDED___
#define ___EllipticMGOperator_H__INCLUDED___

#include "LevelOperator.H"
#include "IntVect.H"
#include "RealVect.H"
#include "GNUC_Extensions.H"

namespace Elliptic {


/**
 * @interface MGOperator
 * @brief     The generic MGOperator for use in all MGSolvers.
 * @date      2019/11/26
 * @author    Edward.Santilli@jefferson.edu
 *
 * @details
 *  The generic MGOperator for use in all MGSolvers. This will act as
 *  its own factory. You create the op at depth 0, and MGSolver will do
 *  the rest by calling newMGOperator().
 */
template <class T>
class MGOperator: public LevelOperator<T>
{
public:
    typedef T                                  StateType;
    typedef typename StateTraits<T>::GridsType GridsType;

    virtual
    ~MGOperator() {}

    /// ------------------------------------------------------------------------
    /// \{
    /// \name {You must override these}

    /// Pure virtual. Factory method.
    /// Allocate + define a coarsened version of *this.
    /// Deletion is left to the caller.
    ///
    /// a_refRatio specifies how much *coarsening* will be needed to get from
    /// this to the new operator.
    ///
    /// The new operator will only be used with homogeneous BCs.
    virtual MGOperator<StateType>*
    newMGOperator(const IntVect& a_refRatio) const __malloc = 0;

    /// What is the smallest Box we can use in MG?
    virtual IntVect
    minBoxSize() const = 0;

    /// What grids is this op defined over?
    virtual const GridsType&
    getBoxes() const = 0;

    /// What is the resolution of this op's grids?
    virtual const RealVect&
    getDXi() const = 0;

    /// Pure virtual.
    /// Restrict to coarser MG depth: a_crseRes = I[h->2h](a_fineRes).
    /// This op is at the fine level.
    /// a_crseRes and a_fineRes should be defined over compatible grids.
    virtual void
    MGRestrict(StateType&                   a_crseRes,
               const StateType&             a_fineRes,
               const Real                   a_time,
               const IntVect&               a_refRatio,
               const MGOperator<StateType>& a_crseOp) const = 0;

    /// Pure virtual.
    /// Prolong (interpolate) a_crse to finer MG depth and add to a_fine:
    ///   a_finePhi += I[2h->h](a_crseCor).
    /// This op is at the fine level.
    /// a_crse and a_fine should be defined over compatible grids.
    ///
    /// We do not overwrite a_finePhi! We ADD the correction.
    ///
    /// If a_interpOrder > 1, then all ghosts of a_crseCor must be extrapolated,
    /// including edge and vertex ghosts.
    virtual void
    MGProlong(StateType&                   a_finePhi,
              StateType&                   a_crseCor,
              const Real                   a_time,
              const IntVect&               a_refRatio,
              const MGOperator<StateType>& a_crseOp,
              const int                    a_interpOrder) const = 0;
    /// \}

protected:
    MGOperator() = default;
    MGOperator(const MGOperator&) = delete;
    MGOperator& operator=(const MGOperator&) = delete;
};


}; // namespace Elliptic
#endif //!___EllipticMGOperator_H__INCLUDED___
