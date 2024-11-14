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
#include "AMRNSLevel.H"

// #include "EllipticOp.H"
#include "ProblemContext.H"
#include "BCTools.H"
#include "FABAlgebra.H"

#include "BoxIterator.H"


// // ======================= BCs for the velocity projector ======================

// // -----------------------------------------------------------------------------
// // Full constructor.
// // -----------------------------------------------------------------------------
// AMRNSLevel::HomogProjBC::HomogProjBC(const LevelGeometry* a_levGeoPtr)
// : m_cfInterp(a_levGeoPtr)
// , m_finiteDiff(*a_levGeoPtr)
// {
//     const IntVect            ghostVect     = IntVect::Unit;
//     const DisjointBoxLayout& grids         = a_levGeoPtr->getBoxes();
//     const ProblemDomain&     domain        = a_levGeoPtr->getDomain();
//     const CFRegion&          cfRegion      = a_levGeoPtr->getCFRegion();
//     const LevelGeometry*     crseLevGeoPtr = a_levGeoPtr->getCoarserPtr();

//     if (crseLevGeoPtr) {
//         m_crseAMRDXi = crseLevGeoPtr->getDXi();
//     } else {
//         m_crseAMRDXi = RealVect(D_DECL(quietNAN, quietNAN, quietNAN));
//     }
//     m_cfiIter.define(grids, cfRegion);

//     m_exCopier.exchangeDefine(grids, ghostVect);
//     m_exCornerCopier.define(grids, grids, domain, ghostVect, true);

//     m_physBdryIter.define(grids);
// }


// // -----------------------------------------------------------------------------
// // Full constructor.
// // -----------------------------------------------------------------------------
// AMRNSLevel::HomogProjBC::HomogProjBC(const LevelGeometry* a_levGeoPtr,
//                                      const RealVect&      a_crseAMRDXi)
//   : m_crseAMRDXi(a_crseAMRDXi)
//   , m_cfInterp(a_levGeoPtr)
//   , m_finiteDiff(*a_levGeoPtr)
// {
//     const IntVect            ghostVect     = IntVect::Unit;
//     const DisjointBoxLayout& grids         = a_levGeoPtr->getBoxes();
//     const ProblemDomain&     domain        = a_levGeoPtr->getDomain();
//     const CFRegion&          cfRegion      = a_levGeoPtr->getCFRegion();
//     const LevelGeometry* __attribute__((unused))     crseLevGeoPtr = a_levGeoPtr->getCoarserPtr();

//     m_cfiIter.define(grids, cfRegion);

//     m_exCopier.exchangeDefine(grids, ghostVect);
//     m_exCornerCopier.define(grids, grids, domain, ghostVect, true);

//     m_physBdryIter.define(grids);
// }


// // -----------------------------------------------------------------------------
// // Destructor.
// // -----------------------------------------------------------------------------
// AMRNSLevel::HomogProjBC::~HomogProjBC()
// {
//     m_physBdryIter.clear();
//     m_cfInterp.undefine();
//     m_finiteDiff.clear();
// }


// // -----------------------------------------------------------------------------
// // Factory method.
// // Be careful! The levGeoPtr may be stored by the new object for use.
// // This version allows a custom a_crseAMRDXi. (Used by AMR solvers).
// // -----------------------------------------------------------------------------
// EllipticOpBC<LevelData<FArrayBox>>*
// AMRNSLevel::HomogProjBC::newEllipticOpBC(const LevelGeometry* a_levGeoPtr,
//                                          const RealVect&      a_crseAMRDXi) const
// {
//     return new AMRNSLevel::HomogProjBC(a_levGeoPtr, a_crseAMRDXi);
// }


// // -----------------------------------------------------------------------------
// // If a_crseStatePtr == NULL, then it is assumed that this function will
// // not need the coarse data or will be able to access the coarse data
// // through some other means.
// // -----------------------------------------------------------------------------
// void
// AMRNSLevel::HomogProjBC::fillGhosts(
//     LevelData<FArrayBox>&       a_state,
//     const Real                  a_time,
//     const LevelGeometry&        a_levGeo,
//     const bool                  a_homogBCs,
//     const LevelData<FArrayBox>* a_crseStatePtr) const
// {
//     CH_assert(m_exCornerCopier.isDefined());

//     const DisjointBoxLayout& grids  = a_levGeo.getBoxes();
//     DataIterator             dit    = grids.dataIterator();
//     FArrayBox                dummyFAB;

//     // Initialize with bogus values, just to make sure we don't use any
//     // uninitialized data.
//     debugInitLevelGhosts(a_state);
//     // BCTools::extrapAllGhosts(a_state, 2);

//     // Set physical BCs.
//     for (m_physBdryIter.reset(); m_physBdryIter.ok(); ++m_physBdryIter) {
//         const DataIndex&     di       = m_physBdryIter->di;
//         const Box&           valid    = m_physBdryIter->ccValidBox;
//         const int            bdryDir  = m_physBdryIter->dir;
//         const Side::LoHiSide side     = m_physBdryIter->side;
//         FArrayBox&           stateFAB = a_state[di];

//         BCTools::neum(stateFAB,
//                       0,
//                       dummyFAB,
//                       -1,
//                       m_finiteDiff,
//                       valid,
//                       di,
//                       bdryDir,
//                       side,
//                       true); // homog?
//     }

//     // Fill CFI ghosts.
//     if (!m_cfiIter.isEmpty()) {
//         if (a_homogBCs) {
//             CFInterp::homogInterpAtCFI(
//                 a_state, a_levGeo.getDXi(), m_crseAMRDXi, m_cfiIter);
//         } else {
//             UNDEFINED_FUNCTION();
//         }
//     }
//     checkForValidNAN(a_state);

//     // Do exchanges.
//     a_state.exchange(m_exCopier);
//     a_state.exchange(m_exCornerCopier);

//     // Fill ghosts at corners of domain. This must happen last!
//     BCTools::extrapDomainCorners(a_state, 2);
// }


// // -----------------------------------------------------------------------------
// // Fills the ghosts of a_state using non-homog BCs at the CFI.
// // Used by AMRMGEllipticOp::AMROperator and AMRMGEllipticOp::AMROperatorNF.
// // -----------------------------------------------------------------------------
// void
// AMRNSLevel::HomogProjBC::fillGhostsNonhomogCFI(
//     LevelData<FArrayBox>&       a_state,
//     const LevelData<FArrayBox>& a_crseState,
//     const Real                  a_time,
//     const LevelGeometry&        a_levGeo,
//     const bool                  a_homogBCs) const
// {
//     CH_assert(m_exCornerCopier.isDefined());

//     const DisjointBoxLayout& grids  = a_levGeo.getBoxes();
//     DataIterator             dit    = grids.dataIterator();
//     FArrayBox                dummyFAB;

//     // Initialize with bogus values, just to make sure we don't use any
//     // uninitialized data.
//     debugInitLevelGhosts(a_state);
//     // BCTools::extrapAllGhosts(a_state, 2);

//     // Set physical BCs.
//     for (m_physBdryIter.reset(); m_physBdryIter.ok(); ++m_physBdryIter) {
//         const DataIndex&     di       = m_physBdryIter->di;
//         const Box&           valid    = m_physBdryIter->ccValidBox;
//         const int            bdryDir  = m_physBdryIter->dir;
//         const Side::LoHiSide side     = m_physBdryIter->side;
//         FArrayBox&           stateFAB = a_state[di];

//         BCTools::neum(stateFAB,
//                       0,
//                       dummyFAB,
//                       -1,
//                       m_finiteDiff,
//                       valid,
//                       di,
//                       bdryDir,
//                       side,
//                       true); // homog?
//     }
// UNDEFINED_FUNCTION();
//     // Fill CFI ghosts.
//     m_cfInterp.interpAtCFI(a_state, a_crseState);
//     checkForValidNAN(a_state);

//     // Do exchanges.
//     a_state.exchange(m_exCopier);
//     a_state.exchange(m_exCornerCopier);

//     // Fill ghosts at corners of domain. This must happen last!
//     BCTools::extrapDomainCorners(a_state, 2);
// }


// // -----------------------------------------------------------------------------
// void
// AMRNSLevel::HomogProjBC::setFluxes(
//     LevelData<FluxBox>&       a_flux,
//     const Real                a_time,
//     const LevelGeometry&      a_levGeo,
//     const bool                a_homogBCs,
//     const LevelData<FluxBox>* a_crseFluxPtr) const
// {
//     BCTools::setValAtPhysBdry(a_flux, 0.0);
// }



