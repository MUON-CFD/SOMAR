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
#ifndef ___ScalarBC_H__INCLUDED___
#define ___ScalarBC_H__INCLUDED___

#include "BCTools.H"
#include "LevelGeometry.H"

namespace ScalarBC {


// -----------------------------------------------------------------------------
/// Sets Grad[state].n = Grad[bgState].n, where
///   bgState = a static vertical background stratification.
///
/// bgState can only be a function of z. It must have FAB indices of 0
///   in the x and y directions.
///
/// If bgState is not provided, it will be set to zero.
///
/// \warning This will only set BCs on a 1-component, CC a_stateFAB.
// -----------------------------------------------------------------------------
struct NormalGradient: public BCTools::BCFunction
{
    /// Constructor -- no background stratification
    NormalGradient();

    /// Constructor -- with a background stratification
    /// a_bgStateFABPtr must have 1 componentm, be CC, and have horizontal
    /// indices equal to zero.
    NormalGradient(const FArrayBox&     a_bgStateFAB,
                   const LevelGeometry& a_levGeo);

    /// The BCFunction override.
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
    std::array<Real, 2> m_bcVal;
};


// -----------------------------------------------------------------------------
struct Tidal: public BCTools::BCFunction
{
    /// Constructor
    Tidal();

    /// The BCFunction override.
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
    RealVect m_tidalOmega;
    RealVect m_tidalInitPhase;
};


// -----------------------------------------------------------------------------
//  Typically, we set ghosts on phi via
//      alpha*phi + beta*dphi/dn = BC.
//
//  If phi is the perturbation to a vertical background stratification,
//  that is, phi = phiTotal - phiBar, then we must alter the rhs via
//      alpha*phi + beta*dphi/dn = BC - (alpha*phiBar + beta*dphiBar/dn).

//  This class modifies a typical BCFunction and sets
//      BC = BC - (alpha*phiBar + beta*dphiBar/dn).
// -----------------------------------------------------------------------------
struct BackgroundScalarWrapper: public BCTools::BCFunction
{
    BackgroundScalarWrapper() = delete;
    BackgroundScalarWrapper(const BackgroundScalarWrapper&) = delete;
    BackgroundScalarWrapper& operator=(const BackgroundScalarWrapper&) = delete;

    virtual ~BackgroundScalarWrapper() {}

    /// Constructor

    /// @brief Constructor
    /// @param a_baseBCPtr     The BCFunction for phiTotal.
    /// @param a_bgStateFABPtr The vertical background stratification.
    /// @param a_dzFABPtr      The dz values along a vertical line of faces that
    ///                        span the domain.
    BackgroundScalarWrapper(
        const std::shared_ptr<BCTools::BCFunction>& a_baseBCPtr,
        const std::shared_ptr<FArrayBox>&           a_bgStateFABPtr,
        const std::shared_ptr<FArrayBox>&           a_dzFABPtr);

    /// The BCFunction override.
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
    const std::shared_ptr<BCTools::BCFunction> m_baseBCPtr;
    const std::shared_ptr<FArrayBox>           m_bgStateFABPtr;
    const std::shared_ptr<FArrayBox>           m_dzFABPtr;
};


}; // end namespace ScalarBC
#endif //!___ScalarBC_H__INCLUDED___
