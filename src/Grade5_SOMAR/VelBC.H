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
#ifndef ___VelBC_H__INCLUDED___
#define ___VelBC_H__INCLUDED___

#include "BCTools.H"

namespace VelBC {


// -----------------------------------------------------------------------------
/// No flow into wall, no flow across wall.
// -----------------------------------------------------------------------------
struct NoSlip: public BCTools::BCFunction
{
    virtual void
    operator()(FArrayBox&            a_alpha,
               FArrayBox&            a_beta,
               FArrayBox&            a_bcFAB,
               const FArrayBox&      a_stateFAB,
               const FArrayBox&      a_xFAB,
               const DataIndex&      a_di,
               const int             a_bdryDir,
               const Side::LoHiSide& a_side,
               const Real            a_time,
               const bool            a_homogBCs) const override;
};


// -----------------------------------------------------------------------------
/// No flow into wall, no shear stress perpendicular to wall.
// -----------------------------------------------------------------------------
struct FreeSlip: public BCTools::BCFunction
{
    virtual void
    operator()(FArrayBox&            a_alpha,
               FArrayBox&            a_beta,
               FArrayBox&            a_bcFAB,
               const FArrayBox&      a_stateFAB,
               const FArrayBox&      a_xFAB,
               const DataIndex&      a_di,
               const int             a_bdryDir,
               const Side::LoHiSide& a_side,
               const Real            a_time,
               const bool            a_homogBCs) const override;
};


// -----------------------------------------------------------------------------
// Sinusoidal in each direction.
// -----------------------------------------------------------------------------
struct Tidal: public BCTools::BCFunction
{
    Tidal();

    virtual void
    operator()(FArrayBox&            a_alpha,
               FArrayBox&            a_beta,
               FArrayBox&            a_bcFAB,
               const FArrayBox&      a_stateFAB,
               const FArrayBox&      a_xFAB,
               const DataIndex&      a_di,
               const int             a_bdryDir,
               const Side::LoHiSide& a_side,
               const Real            a_time,
               const bool            a_homogBCs) const override;

protected:
    RealVect m_tidalU0;
    RealVect m_tidalOmega;
    RealVect m_tidalInitPhase;
};


// -----------------------------------------------------------------------------
struct Inflow: public BCTools::BCFunction
{
    Inflow(const RealVect& a_vel);
    Inflow(std::function<RealVect(const RealVect&)> a_func);

    virtual ~Inflow() {}

    virtual void
    operator()(FArrayBox&            a_alpha,
               FArrayBox&            a_beta,
               FArrayBox&            a_bcFAB,
               const FArrayBox&      a_stateFAB,
               const FArrayBox&      a_xFAB,
               const DataIndex&      a_di,
               const int             a_bdryDir,
               const Side::LoHiSide& a_side,
               const Real            a_time,
               const bool            a_homogBCs) const override;

protected:
    const enum {UNIFORM_VALUE, FUNCTION}         m_type;
    union {
        RealVect                                 m_vel;
        std::function<RealVect(const RealVect&)> m_func;
    };
};


// -----------------------------------------------------------------------------
/// Homogeneous Neumann.
// -----------------------------------------------------------------------------
struct Outflow: public BCTools::BCFunction
{
    virtual void
    operator()(FArrayBox&            a_alpha,
               FArrayBox&            a_beta,
               FArrayBox&            a_bcFAB,
               const FArrayBox&      a_stateFAB,
               const FArrayBox&      a_xFAB,
               const DataIndex&      a_di,
               const int             a_bdryDir,
               const Side::LoHiSide& a_side,
               const Real            a_time,
               const bool            a_homogBCs) const override;
};


}; // end namespace VelBC
#endif //!___VelBC_H__INCLUDED___
