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
#ifndef ___LepticOperator_H__INCLUDED___
#define ___LepticOperator_H__INCLUDED___

#include <memory>
#include "LDFABOps.H"
#include "RealVect.H"
#include "BCTools.H"
#include "MGOperator.H"

namespace Elliptic {


/// @brief Provides addl interface methods for use by LepticSolver.
/// For an op to work in LepticSolver, it must be a LevelOperator AND a LepticOp.
class LepticOperator
{
public:
    typedef LevelData<FArrayBox>                        StateType;
    typedef typename StateTraits<StateType>::GridsType  GridsType;
    typedef typename StateTraits<StateType>::CopierType CopierType;

    virtual ~LepticOperator() {}

    virtual const RealVect&           getDomainLength() const = 0;
    virtual const RealVect&           getDXi()          const = 0;
    virtual const RealVect&           getDXiCrse()      const = 0;
    virtual const LevelData<FluxBox>& getFCJgup()       const = 0;
    virtual IntVect                   getBlockFactor()  const = 0;

    virtual std::shared_ptr<LevelOperator<LevelData<FArrayBox>>>
    createLevelOperator(const DisjointBoxLayout& a_grids) const = 0;

    virtual std::shared_ptr<MGOperator<LevelData<FArrayBox>>>
    createHorizontalMGOperator(const DisjointBoxLayout& a_hgrids) const = 0;
};


}; // namespace Elliptic
#endif //!___LepticOperator_H__INCLUDED___
