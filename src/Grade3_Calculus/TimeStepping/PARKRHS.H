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
#ifndef ___PARKRHS_H__INCLUDED___
#define ___PARKRHS_H__INCLUDED___

#include "LevelData.H"
#include "FluxBox.H"
#include "GNUC_Extensions.H"


/**
 * \class   PARKRHS
 * \brief   Specifies the RHS of dQ/dt = kqE(q(t),t) + kqI(q(t),t) - Gp.
 * \date    2019/12/06
 * \author  Edward.Santilli@jefferson.edu
 *
 * \details
 *  We want to use PARK to solve
 *    dQ/dt = kQE(Q(t),t) + kQI(Q(t),t) - Gp,
 *  where
 *    Q   = the state vector including vel, q, and p
 *    kQE = non-stiff forcing
 *    kQI = stiff forcing
 *    vel = the face-centered velocity
 *    p   = cell-centered pressure
 *    q   = cell-centered scalars not including pressure
 *
 *  To do this, PARK needs a few problem-specific things that this class must
 *  be overridden to provide.
 */
class PARKRHS
{
public:
    virtual ~PARKRHS() {}

    /// The explicit forcing. Be sure to set the comps with no
    /// forcing to zero. Q is non-const so that you can set BCs.
    virtual void
    setExplicitRHS(LevelData<FluxBox>&   a_kvel,
                   LevelData<FArrayBox>& a_kq,
                   LevelData<FluxBox>&   a_vel,
                   LevelData<FArrayBox>& a_p,
                   LevelData<FArrayBox>& a_q,
                   const Real            a_time,
                   const Real            a_refluxDt) = 0;

    /// The implicit forcing. Be sure to set the comps with no
    /// forcing to zero. Q is non-const so that you can set BCs.
    virtual void
    setImplicitRHS(LevelData<FluxBox>&   a_kvel,
                   LevelData<FArrayBox>& a_kq,
                   LevelData<FluxBox>&   a_vel,
                   LevelData<FArrayBox>& a_p,
                   LevelData<FArrayBox>& a_q,
                   const Real            a_time,
                   const Real            a_refluxDt) = 0;

    /// Solves [1 - gammaDt * D] Q^{n+1} = Q^{n} in place.
    virtual void
    solveImplicit(LevelData<FluxBox>&   a_vel,
                  LevelData<FArrayBox>& a_q,
                  const Real            a_gammaDt,
                  const Real            a_time,
                  const Real            a_refluxDt) = 0;

    /// Use an approximation of p to project. Do not solve for p.
    /// The p provided is the lagged pressure computed in the previous RK stage.
    virtual void
    projectPredict(LevelData<FluxBox>&   a_vel,
                   LevelData<FArrayBox>& a_p,
                   const Real            a_time,
                   const Real            a_projDt) = 0;

    /// Solve for p and project. This will be an incremental solve
    /// since we already called projectPredict using the lagged p.
    virtual void
    projectCorrect(LevelData<FluxBox>&   a_vel,
                   LevelData<FArrayBox>& a_p,
                   const Real            a_time,
                   const Real            a_projDt,
                   const Real            a_gammaDt = 0.0) = 0;
};


#endif //!___PARKRHS_H__INCLUDED___
