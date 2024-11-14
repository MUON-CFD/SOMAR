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
#ifndef ___BCTools_H__INCLUDED___
#define ___BCTools_H__INCLUDED___

#include "FiniteDiff.H"
#include "CornerCopier.H"
#include "SideArray.H"
#include "FluxBox.H"
#include "BdryIter.H"
#include "GNUC_Extensions.H"


namespace BCTools {


// -----------------------------------------------------------------------------
// A suite of helper functions to set BCs, extrapolate to ghosts, etc.
//
// WARNING:
//  These tools are designed to work at any boundary. It does not check if we
//  are at a physical boundary, exchange boundary, or CFI. That job is up to
//  the caller.
// -----------------------------------------------------------------------------


// =============================================================================
/// \name Inline interpolation functions
/// \{

/// Linear extrapolation formula. Assumes all data points are equidistant.
inline Real
linearExtrap(const Real& a_nearVal, const Real& a_farVal)
{
    return 2.0 * a_nearVal - a_farVal;
}


/// Quadratic extrapolation formula. Assumes all data points are equidistant.
inline Real
quadraticExtrap(const Real& a_nearVal,
                const Real& a_midVal,
                const Real& a_farVal)
{
    return 3.0 * (a_nearVal - a_midVal) + a_farVal;
}


/// Cubic extrapolation formula. Assumes all data points are equidistant.
inline Real
cubicExtrap(const Real& a_val1,  // closest to extraped value
            const Real& a_val2,
            const Real& a_val3,
            const Real& a_val4)  // farthest extraped value
{
    return 4.0 * (a_val1 + a_val3) - 6.0 * a_val2 - a_val4;
}


/// Quartic extrapolation formula. Assumes all data points are equidistant.
inline Real
quarticExtrap(const Real& a_val1,  // closest to extraped value
              const Real& a_val2,
              const Real& a_val3,
              const Real& a_val4,
              const Real& a_val5)  // farthest from extraped value
{
    return 5.0 * (a_val1 - a_val4) + 10.0 * (a_val3 - a_val2) + a_val5;
}
/// \}


// =============================================================================
/// \name BC-setting functions on a single FAB
/// \note Unlike the LevelData versions, these functions that set BCs on
///       FABS do not check if we are at a computational boundary. They blindly
///       fill the ghosts you feed them.
/// \{


/// \brief Extrapolation boundary conditions for a side.
/// \details
/// This can be used to extrapolate ghosts or the boundary face of a
/// FC FArrayBox if you choose your a_valid Box correctly.
/// a_state and a_valid must have the same centering.
/// If a_state does not contain ghosts, we just exit quietly.
/// \note This function does not check if we are at a computational boundary!
void
extrap(FArrayBox&           a_state,
       const Box&           a_valid,
       const int            a_dir,
       const Side::LoHiSide a_side,
       const int            a_order);

/// Assumes all ghosts except edges and vertices are filled and uses that
/// data to extrapolate to edges and vertices.
/// \note This only fills ghosts in the first layer.
void
extrapCorners(FArrayBox& a_state, const Box& a_valid, const int a_order);


/// a_xFAB must have SpaceDim comps filled with physical coordinates.
/// a_alpha, a_beta, a_bcFAB, and a_stateFAB must have the same # of comps.
/// Where the BC will be set according to alpha * s + beta * ds/dn = B.
/// Note that ds/dn is an outward normal derivative!
struct BCFunction {
    virtual ~BCFunction() {}

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
               const bool            a_homogBCs) const = 0;
};


/// Collects 2*SpaceDim BCFunctions, one for each bdryDir/side, and turns them
/// into a single BCFunction.
struct UberBCFunction : public BCFunction
{
    UberBCFunction(std::array<std::array<std::shared_ptr<BCFunction>, 2>,
                              SpaceDim>& a_bcFuncPtrs)
    : m_bcFuncPtrs(a_bcFuncPtrs)
    {
    }

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
               const bool            a_homogBCs) const override
    {
        const auto bcptr = m_bcFuncPtrs[a_bdryDir][int(a_side)];
        CH_assert(bcptr);
        bcptr->operator()(a_alpha,
                          a_beta,
                          a_bcFAB,
                          a_stateFAB,
                          a_xFAB,
                          a_di,
                          a_bdryDir,
                          a_side,
                          a_time,
                          a_homogBCs);
    }

protected:
    std::array<std::array<std::shared_ptr<BCFunction>, 2>, SpaceDim>
        m_bcFuncPtrs;
};


/// Collects 2*SpaceDim BCFunctions, one for each bdryDir/side, and turns them
/// into a single BCFunction pointer.
inline std::shared_ptr<BCFunction>
makeUberBCFunctionPtr(std::array<std::array<
                      std::shared_ptr<BCFunction>, 2>, SpaceDim>& a_bcFuncPtrs)
{
    return std::shared_ptr<BCFunction>(new UberBCFunction(a_bcFuncPtrs));
}


/// Very simple BCFunction - Homogeneous Dirichlet.
struct HomogDiriBC : public BCTools::BCFunction {
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
               const bool            a_homogBCs) const override
    {
        a_alpha.setVal(1.0);
        a_beta.setVal(0.0);
        a_bcFAB.setVal(0.0);
    }
};


/// Very simple BCFunction - Homogeneous Neumann.
struct HomogNeumBC : public BCTools::BCFunction {
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
               const bool            a_homogBCs) const override
    {
        a_alpha.setVal(0.0);
        a_beta.setVal(1.0);
        a_bcFAB.setVal(0.0);
    }
};


/// Inhomogeneous BCFunction for a single side.
/// Sets BCs according to alpha*phi + beta*d(phi)/dn = bcVals,
/// where d/dn is the *outward* normal derivative.
/// \warning This will cause issues if you define it on one level and use it on
/// another!
struct RobinBC_CCSide : public BCTools::BCFunction {
    RobinBC_CCSide(const int             a_bdryDir,
                   const Side::LoHiSide& a_side,
                   const ProblemDomain&  a_domain,
                   const Real            a_alpha,
                   const Real            a_beta,
                   FArrayBox&            a_bcValsFAB);

    virtual void
    operator()(FArrayBox&            a_alphaFAB,
               FArrayBox&            a_betaFAB,
               FArrayBox&            a_bcFAB,
               const FArrayBox&      a_stateFAB,
               const FArrayBox&      a_xFAB,
               const DataIndex&      a_di,
               const int             a_bdryDir,
               const Side::LoHiSide& a_side,
               const Real            a_time,
               const bool            a_homogBCs) const override;

protected:
    int            m_bdryDir;
    Side::LoHiSide m_side;
    Real           m_alpha;
    Real           m_beta;
    FArrayBox      m_bcValsFAB;
};


/// Inhomogeneous BCFunction for a single side.
/// Sets BCs according to alpha*phi + beta*d(phi)/dn = bcVals,
/// where d/dn is the *outward* normal derivative.
/// \warning This will cause issues if you define it on one level and use it on
/// another!
struct RobinBC_FCSide : public BCTools::BCFunction {
    RobinBC_FCSide(const int                        a_bdryDir,
                   const Side::LoHiSide&            a_side,
                   const ProblemDomain&             a_domain,
                   const RealVect&                  a_alpha,
                   const RealVect&                  a_beta,
                   std::array<FArrayBox, SpaceDim>& a_bcValsFAB);

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
    int                             m_bdryDir;
    Side::LoHiSide                  m_side;
    RealVect                        m_alpha;
    RealVect                        m_beta;
    std::array<FArrayBox, SpaceDim> m_bcValsFAB;
};


/// Defines two boxes:
///  a_stateBdry   = where the BCs need to be specified.
///          For a CC a_stateFAB, this will be FC at the boundary.
///          For a_stateFAB FC in a_bdryDir, this will be FC at the boundary.
///          For a_stateFAB FC in any other direction, this will be
///             edge-centered in 3D / node-centered in 2D at the boundary.
///  a_stateGhosts = a_state's ghost region that needs to be filled.
/// The return value will be the IndexType of a_stateFAB in a_bdryDir.
int
getBoxesForApplyBC(Box&                  a_stateBdry,
                   Box&                  a_stateGhosts,
                   const FArrayBox&      a_stateFAB,
                   const Box&            a_ccValid,
                   const int             a_bdryDir,
                   const Side::LoHiSide& a_side);


/// Sets the BC according to alpha * s + beta * ds/dn = B.
/// a_bcFuncPtr is a custom BCFunction that sets alpha, beta, and B.
/// a_bcFuncPtr must be able to set BCs on all of a_stateFAB's comps.
/// Note that ds/dn is an outward normal derivative!
/// a_xFAB and a_dxFAB must be defined at the physical boundary face,
///   and have the same centering as a_stateFAB in all other directions.
/// a_xFAB must have SpaceDim comps. a_dxFAB just needs 1 comp.
/// a_dxFAB must contain dx = dx/dXi * dXi.
/// a_ccValid must be CC. It might not coincide with a_stateFAB's centering.
void
applyBC(FArrayBox&            a_stateFAB,
        const FArrayBox&      a_xFAB,
        const FArrayBox&      a_dxFAB,
        const Box&            a_ccValid,
        const DataIndex&      a_di,
        const int             a_bdryDir,
        const Side::LoHiSide& a_side,
        const bool            a_homogBCs,
        const BCFunction&     a_bcFunc,
        const Real            a_time);
/// \}


// =============================================================================
/// \name BC-setting functions on an entire level
/// \{

/// \brief Sets the BC according to alpha * s + beta * ds/dn = B.
/// \details
/// a_bcFuncPtr is a custom BCFunction that sets alpha, beta, and B.
/// a_bcFuncPtr must be able to set BCs on all of a_state's comps.
/// Note that ds/dn is an outward normal derivative!
void
applyBC(LevelData<FArrayBox>&     a_state,
        const bool                a_homogBCs,
        const BCFunction&         a_bcFunc,
        const Real                a_time,
        const GeoSourceInterface& a_geoSrc,
        const RealVect&           a_dXi,
        PhysBdryIter&             a_physBdryIter);


/// \brief Sets the BC according to alpha * s + beta * ds/dn = B.
/// \details
/// a_bcFuncPtr is a custom BCFunction that sets alpha, beta, and B.
/// a_bcFuncPtr must be able to set BCs on all of a_state's comps.
/// Note that ds/dn is an outward normal derivative!
inline void
applyBC(LevelData<FArrayBox>& a_state,
        const bool            a_homogBCs,
        const BCFunction&     a_bcFunc,
        const Real            a_time,
        const LevelGeometry&  a_levGeo,
        PhysBdryIter&         a_physBdryIter)
{
    BCTools::applyBC(a_state,
                     a_homogBCs,
                     a_bcFunc,
                     a_time,
                     a_levGeo.getGeoSource(),
                     a_levGeo.getDXi(),
                     a_physBdryIter);
}


/// \brief Sets the BC according to alpha * s + beta * ds/dn = B.
/// \details
/// a_bcFuncPtr is a custom BCFunction that sets alpha, beta, and B.
/// a_bcFuncPtr must be able to set BCs on all of a_state's comps.
/// Note that ds/dn is an outward normal derivative!
///
/// \note Constructing the PhysBdryIter is a bottleneck! Consider using the
/// version that takes in a cached PhysBdryIter.
inline void
applyBC(LevelData<FArrayBox>& a_state,
        const bool            a_homogBCs,
        const BCFunction&     a_bcFunc,
        const Real            a_time,
        const LevelGeometry&  a_levGeo,
        const SideArray&      a_sideArray = AllSides)
{
    PhysBdryIter it(a_state.getBoxes(), a_sideArray);
    BCTools::applyBC(a_state,
                     a_homogBCs,
                     a_bcFunc,
                     a_time,
                     a_levGeo.getGeoSource(),
                     a_levGeo.getDXi(),
                     it);
}


/// This is the FC version of applyBC.
/// This will only set BCs on one FC component of each FluxBox.
/// So, for example, if you want to set velocity BCs in 3D, you'll need to
/// call this function 3 times - once for each vector component.
void
applyBC(LevelData<FluxBox>&  a_state,
        const int            a_fcDir,
        const bool           a_homogBCs,
        const BCFunction&    a_bcFunc,
        const Real           a_time,
        const LevelGeometry& a_levGeo,
        PhysBdryIter&        a_physBdryIter);


/// This is the FC version of applyBC.
/// This will only set BCs on one FC component of each FluxBox.
/// So, for example, if you want to set velocity BCs in 3D, you'll need to
/// call this function 3 times - once for each vector component.
///
/// \note Constructing the PhysBdryIter is a bottleneck! Consider using the
/// version that takes in a cached PhysBdryIter.
inline void
applyBC(LevelData<FluxBox>&  a_state,
        const int            a_fcDir,
        const bool           a_homogBCs,
        const BCFunction&    a_bcFunc,
        const Real           a_time,
        const LevelGeometry& a_levGeo,
        const SideArray&     a_sideArray = AllSides)
{
    PhysBdryIter it(a_state.getBoxes(), a_sideArray);
    BCTools::applyBC(
        a_state, a_fcDir, a_homogBCs, a_bcFunc, a_time, a_levGeo, it);
}


/// Extrapolates one layer of ghosts. This does not dellineate between
/// exchange, physical boundary, or CFI ghosts. It just fills 'em all.
///
/// Order matters! If one or more dirs are periodic:
///  If you exchange before calling extrap, corner ghosts will be valid.
///  If you exchange after, you'll need to use a CornerCopier too.
///  If all dirs are periodic, then you must use a CornerCopier.
///
/// If you want to be safe, exchange ALL ghosts before or after extrap.
void
extrapAllGhosts(LevelData<FArrayBox>& a_state,
                const int             a_order,
                const IntVect&        a_skipGhosts = IntVect::Zero);


/// FC version.
void
extrapAllGhosts(LevelData<FluxBox>& a_state,
                const int           a_order,
                const IntVect&      a_skipGhosts = IntVect::Zero);


/// Fill ghosts at edges and vertices of the domain.
/// This only fills ghosts in the first layer.
void
extrapDomainCorners(LevelData<FArrayBox>& a_state, const int a_order);



/// Fill ghosts at edges and vertices of the domain.
/// This only fills ghosts in the first layer.
void
extrapDomainCorners(LevelData<FluxBox>& a_state, const int a_order);


///
void
setValAtPhysBdry(LevelData<FluxBox>& a_flux,
                 const Real          a_val,
                 const SideArray&    a_sideArray = AllSides);


///
void
setValAtCFI(LevelData<FluxBox>& a_flux,
            const Real          a_val,
            const CFRegion&     a_cfRegion,
            const SideArray&    a_sideArray = AllSides);
/// \}


/// A simple tool for finding the FC dir of a box.
/// Returns -1 if the box is not FC.
int
getNodalDir(const Box& a_box);



}; // end namespace BCTools
#endif //!___BCTools_H__INCLUDED___
