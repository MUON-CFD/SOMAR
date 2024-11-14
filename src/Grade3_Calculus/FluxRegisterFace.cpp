#include "FluxRegisterFace.H"
#include "AnisotropicRefinementTools.H"
#include "CFInterp.H"
#include "Masks.H"
#include "FABAlgebra.H"
#include "SOMAR_Constants.H"
// #include "AnisotropicFluxRegisterF_F.H"
#include "FluxRegisterFaceF_F.H"
#include "Debug.H"

#include "IO.H" // TEMPORARY!!!
#include "SetValLevel.H"


// -----------------------------------------------------------------------------
// Default constructor. Leaves object undefined.
// -----------------------------------------------------------------------------
FluxRegisterFace::FluxRegisterFace()
: m_isDefined(false)
, m_isNoOp(false)
, m_doNormalFluxes(false)
, m_doTransverseFluxes(false)
, m_refRatio(D_DECL(-1, -1, -1))
, m_numComps(-1)
{
}


// -----------------------------------------------------------------------------
// Full constructor. Leaves object in a usable state.
// -----------------------------------------------------------------------------
FluxRegisterFace::FluxRegisterFace(const LevelGeometry* a_fineLevGeoPtr,
                                   const LevelGeometry* a_crseLevGeoPtr,
                                   const int            a_numComps,
                                   const bool           a_doNormalFluxes,
                                   const bool           a_doTransverseFluxes)
: m_isDefined(false)
, m_isNoOp(false)
, m_doNormalFluxes(false)
, m_doTransverseFluxes(false)
, m_refRatio(D_DECL(-1, -1, -1))
, m_numComps(-1)
{
    this->define(a_fineLevGeoPtr,
                 a_crseLevGeoPtr,
                 a_numComps,
                 a_doNormalFluxes,
                 a_doTransverseFluxes);
}


// -----------------------------------------------------------------------------
// Full virtual constructor. Leaves object in a usable state.
// -----------------------------------------------------------------------------
void
FluxRegisterFace::define(const LevelGeometry* a_fineLevGeoPtr,
                         const LevelGeometry* a_crseLevGeoPtr,
                         const int            a_numComps,
                         const bool           a_doNormalFluxes,
                         const bool           a_doTransverseFluxes)
{
    this->clear();

    // Is there a coarser level?
    m_isNoOp = false;
    if (a_crseLevGeoPtr == nullptr) {
        m_isNoOp = true;
        m_isDefined = true;
        return;
    }

    // These members must be set before calling define's helper functions.
    {
        m_doNormalFluxes     = a_doNormalFluxes;
        m_doTransverseFluxes = a_doTransverseFluxes;

        if (a_doNormalFluxes) {
            MAYDAYWARNING("I think normal refluxing is useless. Not coded up.");
        }

        m_refRatio = calculateRefinementRatio(a_crseLevGeoPtr->getDomainBox(),
                                              a_fineLevGeoPtr->getDomainBox());
        CH_assert(m_refRatio >= IntVect::Unit);

        m_fineGrids = a_fineLevGeoPtr->getBoxes();
        ::coarsen(m_coarsenedFineGrids, m_fineGrids, m_refRatio);
        m_crseGrids = a_crseLevGeoPtr->getBoxes();

        m_numComps = a_numComps;

        for (size_t fcDir = 0; fcDir < SpaceDim; ++fcDir) {
            m_reverseCopier[fcDir].staggeredGhostDefine(
                m_coarsenedFineGrids,
                m_crseGrids,
                m_crseGrids.physDomain(),
                IntVect::Unit - BASISV(fcDir),
                fcDir);
        }
    }

    // Does the fine level cover the entire domain?
// #define DISABLE_TEMPORARY_FLUX_REGISTER_OPTIMIZATION
#ifndef DISABLE_TEMPORARY_FLUX_REGISTER_OPTIMIZATION
    {
        LayoutIterator lit;
        int numPts = 0;

        lit = m_crseGrids.layoutIterator();
        for (lit.reset(); lit.ok(); ++lit) {
            numPts += m_crseGrids[lit].numPts();
        }

        lit = m_coarsenedFineGrids.layoutIterator();
        for (lit.reset(); lit.ok(); ++lit) {
            numPts -= m_coarsenedFineGrids[lit].numPts();
        }

        if (numPts == 0) {
            m_isNoOp = true;
            m_isDefined = true;
            return;
        }
    }
#endif

    // Identify FC coarse velocity locations to be updated.
    for (size_t n = 0; n < SpaceDim; ++n) {
        for (size_t f = 0; f < SpaceDim; ++f) {
            for (SideIterator sit; sit.ok(); ++sit) {
                const size_t iside = static_cast<size_t>(sit()); // cfi side

                if (f == n) {
                    FluxRegisterFace::normalLocations(
                        m_crseLocations[n][f][iside],
                        m_coarsenedFineGrids,
                        m_crseGrids,
                        n,
                        sit());
                } else {
                    FluxRegisterFace::transverseLocations(
                        m_crseLocations[n][f][iside],
                        m_coarsenedFineGrids,
                        m_crseGrids,
                        n,
                        f,
                        sit());
                }
            } // sit
        } // fcDir
    } // bdryDir

    // Define flux registers.
    m_coarsenedFineFlux.define(m_coarsenedFineGrids, a_numComps, IntVect::Unit);
    m_crseFlux.define(m_crseGrids, a_numComps);

    // This object is ready for use.
    m_isDefined = true;

    this->setToZero();
}


// -----------------------------------------------------------------------------
// Virtual destructor.
// -----------------------------------------------------------------------------
FluxRegisterFace::~FluxRegisterFace()
{
    this->clear();
}


// -----------------------------------------------------------------------------
// Free memory. Leaves object undefined.
// -----------------------------------------------------------------------------
void
FluxRegisterFace::clear()
{
    for (size_t fcDir = 0; fcDir < SpaceDim; ++fcDir) {
        m_reverseCopier[fcDir].clear();
    }

    m_crseFlux.clear();
    m_coarsenedFineFlux.clear();
    m_numComps = -1;

    for (int n = 0; n < SpaceDim; ++n) {     // = cfiDir / derivDir
        for (int f = 0; f < SpaceDim; ++f) { // = fcDir
            for (SideIterator sit; sit.ok(); ++sit) {
                const size_t iside = static_cast<size_t>(sit());
                auto& loc = m_crseLocations[n][f][iside];

                for (DataIterator dit = loc.dataIterator(); dit.ok(); ++dit) {
                    loc[dit].clear();
                }
            } // sit
        } // f
    } // n

    m_crseGrids          = DisjointBoxLayout();
    m_coarsenedFineGrids = DisjointBoxLayout();
    m_fineGrids          = DisjointBoxLayout();

    m_refRatio = IntVect(D_DECL(-1, -1, -1));

    m_doTransverseFluxes = false;
    m_doNormalFluxes     = false;

    m_isNoOp    = false;
    m_isDefined = false;
}


// -----------------------------------------------------------------------------
// Initialize the flux register.
// -----------------------------------------------------------------------------
void
FluxRegisterFace::setToZero()
{
    CH_assert(m_isDefined);
    if (m_isNoOp) return;

    setValLevel(m_coarsenedFineFlux, 0.0);
    setValLevel(m_crseFlux, 0.0);
}


// -----------------------------------------------------------------------------
void
FluxRegisterFace::incrementCoarse(const FArrayBox&      a_crseFluxFAB,
                                  const Real            a_scale,
                                  const DataIndex&      a_crseDI,
                                  const Interval&       a_srcInterval,
                                  const Interval&       a_destInterval,
                                  const int             a_fcDir,
                                  const int             a_divDir,
                                  const Side::LoHiSide& a_cfiSide)
{
    CH_assert(m_isDefined);

    if (m_isNoOp) return;
    if (RealCmp::isZero(a_scale)) return;

    const int n     = a_divDir;  // advectING dir
    const int f     = a_fcDir;   // advectED dir
    const int s     = int(a_cfiSide);
    const int isign = sign(a_cfiSide);

    const int numComps = a_srcInterval.size();
    const int srcComp  = a_srcInterval.begin();
    const int destComp = a_destInterval.begin();

    FArrayBox&         regFAB       = m_crseFlux[a_crseDI][f];
    const Vector<Box>& crseVelBoxes = m_crseLocations[n][f][s][a_crseDI];
    const Real         localScale   = Real(isign) * a_scale;

    // Subtract data from the existing register.
    if ((n == f) && m_doNormalFluxes) {
        MAYDAYERROR("Normal refluxing is not complete.");

    } else if ((n != f) && m_doTransverseFluxes) {
        // Shift register to flux location.
        regFAB.shiftHalf(n, -isign);

        // Create copy mask.
        FArrayBox maskFAB(regFAB.box(), 1);
        maskFAB.setVal(0.0);
        for (Box bx : crseVelBoxes) {
            bx.shiftHalf(n, -isign);
            maskFAB.setVal(localScale, bx, 0);
        }

        // Add to register.
        for (int c = 0; c < numComps; ++c) {
            FABAlgebra::AddProd2(regFAB,
                                 destComp + c,
                                 a_crseFluxFAB,
                                 srcComp + c,
                                 maskFAB,
                                 0,
                                 regFAB.box());
        }

        // Shift register back to vel location.
        regFAB.shiftHalf(n, isign);

    }  // if / if not n == f
}


// -----------------------------------------------------------------------------
void
FluxRegisterFace::incrementFine(const FArrayBox&      a_fineFluxFAB,
                                const Real            a_scale,
                                const DataIndex&      a_fineDI,
                                const Interval&       a_srcInterval,
                                const Interval&       a_destInterval,
                                const int             a_fcDir,
                                const int             a_divDir,
                                const Side::LoHiSide& a_cfiSide)
{
    // UNDEFINED_FUNCTION();

    CH_assert(m_isDefined);

    if (m_isNoOp) return;
    if (RealCmp::isZero(a_scale)) return;

    const int  n          = a_divDir;
    const int  f          = a_fcDir;
    const int  isign      = sign(a_cfiSide);

    const int  numComps   = a_srcInterval.size();
    const int  srcComp    = a_srcInterval.begin();
    const int  destComp   = a_destInterval.begin();

    FArrayBox& regFAB     = m_coarsenedFineFlux[a_fineDI][f];
    const Real localScale = -Real(isign) * a_scale;

    // const bool       doHarmonicAvg = false;
    const FArrayBox* fineJFABPtr   = nullptr;  // Fluxes are already J-scaled.

    // Add data to the existing register.
    if ((n == f) & m_doNormalFluxes) {
        MAYDAYERROR("Normal refluxing is not complete.");

    } else if ((n != f) && m_doTransverseFluxes) {
        // Locate fluxes to be updated.
        const Box edgeFluxBox =
            bdryBox(m_coarsenedFineGrids[a_fineDI], n, a_cfiSide)
                .surroundingNodes(f);

        // Shift register to flux location.
        regFAB.shiftHalf(n, -isign);

        // Average fluxes to coarse resolution.
        FArrayBox coarsenedFluxesFAB(edgeFluxBox, regFAB.nComp());
        if (SpaceDim == 2) {
            CFInterp::localCoarsenNode(coarsenedFluxesFAB,
                                       a_fineFluxFAB,
                                       edgeFluxBox,
                                       m_refRatio);
        } else {
            const int edgeDir = SpaceDim - n - f;
            CFInterp::localCoarsenEdge(coarsenedFluxesFAB,
                                       a_fineFluxFAB,
                                       edgeFluxBox,
                                       edgeDir,
                                       m_refRatio,
                                       fineJFABPtr);
        }

        // Add to register.
        regFAB.plus(coarsenedFluxesFAB,
                    edgeFluxBox,
                    edgeFluxBox,
                    localScale,
                    srcComp,
                    destComp,
                    numComps);

        // Shift register back to vel location.
        regFAB.shiftHalf(n, isign);

    } // if / if not n == f
}


// -----------------------------------------------------------------------------
// Reflux does:
//  deltaF = fine flux - crse flux
//  face divergence -= iside * scale * deltaF
// -----------------------------------------------------------------------------
void
FluxRegisterFace::reflux(LevelData<FluxBox>& a_crseDiv,
                         const Real          a_scale)
{
    CH_assert(m_isDefined);
    if (m_isNoOp) return;
    if (RealCmp::isZero(a_scale)) return;

    CH_assert(a_crseDiv.getBoxes() == m_crseGrids);
    CH_assert(a_crseDiv.nComp() == m_numComps);

    // Do fine side.
    LevelData<FluxBox> dF(m_crseGrids, m_numComps);
    setValLevel(dF, 0.0);
    m_coarsenedFineFlux.exchange(); // TODO: Is this needed?
    m_coarsenedFineFlux.copyTo(dF, m_reverseCopier);

    // Do coarse side.
    for (DataIterator dit(m_crseGrids); dit.ok(); ++dit) {
        for (int f = 0; f < SpaceDim; ++f) {
            FArrayBox&       dFFAB        = dF[dit][f];
            const FArrayBox& crseFluxFAB  = m_crseFlux[dit][f];
            const Box&       bx           = dFFAB.box();

            dFFAB.plus(crseFluxFAB, bx, 0, 0, m_numComps);
        } // f
    } // dit

    // Just in case the user averaged down before refluxing, let's set the
    // invalid dF to zero.
    Masks::zeroInvalid(dF, &m_fineGrids);

    // At this point, dF = the reflux divergence.
    // Apply it to a_crseDiv with the user's scaling.
    for (DataIterator dit(m_crseGrids); dit.ok(); ++dit) {
        for (int f = 0; f < SpaceDim; ++f) {
            FArrayBox&       crseDivFAB = a_crseDiv[dit][f];
            const FArrayBox& dFFAB      = dF[dit][f];
            const Box        fcValid    = m_crseGrids[dit].surroundingNodes(f);

            FArrayBox maskFAB(fcValid, 1);
            maskFAB.setVal(0.0);

            for (int n = 0; n < SpaceDim; ++n) {
                if (f == n && !m_doNormalFluxes) continue;
                if (f != n && !m_doTransverseFluxes) continue;

                for (SideIterator sit; sit.ok(); ++sit) {
                    const size_t iside = size_t(sit());
                    const auto&  boxes = m_crseLocations[n][f][iside][dit()];

                    for (const auto& bx : boxes) {
                        maskFAB.setVal(a_scale, bx, 0);
                    }
                } // sit
            } // n

            for (int c = 0; c < m_numComps; ++c) {
                FABAlgebra::AddProd2(
                    crseDivFAB, c, dFFAB, c, maskFAB, 0, fcValid);
            }
        }  // f
    } // dit

    // CFInterp::validateAtCoarseCFI(a_crseDiv, m_fineGrids, m_refRatio);
}


// -----------------------------------------------------------------------------
void
FluxRegisterFace::normalLocations(
    LayoutData<Vector<Box>>& a_locations,
    const DisjointBoxLayout& a_coarsenedFineGrids,
    const DisjointBoxLayout& a_crseGrids,
    const size_t             a_bdryDir,
    const Side::LoHiSide     a_cfiSide)
{
    a_locations.define(a_crseGrids);

//     // for (DataIterator dit(a_crseGrids); dit.ok(); ++dit) {
//     //     const Box crseTestBox = a_crseGrids[dit]; // Keep CC

//     //     a_locations[dit].clear();

//     //     LayoutIterator lit = a_coarsenedFineGrids.layoutIterator();
//     //     for (lit.reset(); lit.ok(); ++lit) {
//     //         const Box fineFluxBox =
//     //             adjCellBox(a_coarsenedFineGrids[lit], a_bdryDir, a_side, 1);

//     //         if (crseTestBox.intersectsNotEmpty(fineFluxBox)) {
//     //             a_locations[dit].push_back(crseTestBox & fineFluxBox);
//     //         }
//     //     } // lit over coarsenedFineGrids
//     // } // dit over crseGrids

//     // if (a_crseGrids.physDomain().isPeriodic()) {
//     //     UNDEFINED_FUNCTION();
//     // } // end if periodic
}


// -----------------------------------------------------------------------------
void
FluxRegisterFace::transverseLocations(
    LayoutData<Vector<Box>>& a_locations,
    const DisjointBoxLayout& a_coarsenedFineGrids,
    const DisjointBoxLayout& a_crseGrids,
    const size_t             a_bdryDir,
    const size_t             a_fcDir,
    const Side::LoHiSide     a_cfiSide)
{
    a_locations.define(a_crseGrids);

    // NOTE: Currently a no-op until I work out the periodic code.
    // Non-periodic code works.

    // const ProblemDomain& crseDomain = a_crseGrids.physDomain();

    // for (DataIterator dit(a_crseGrids); dit.ok(); ++dit) {
    //     const Box crseTestBox = a_crseGrids[dit]
    //                            .surroundingNodes(a_fcDir);  // Make FC

    //     a_locations[dit].clear();

    //     LayoutIterator lit = a_coarsenedFineGrids.layoutIterator();
    //     for (lit.reset(); lit.ok(); ++lit) {
    //         Box fineFluxBox = a_coarsenedFineGrids[lit];
    //         fineFluxBox = adjCellBox(fineFluxBox, a_bdryDir, a_cfiSide, 1);
    //         fineFluxBox.surroundingNodes(a_fcDir);

    //         if (crseTestBox.intersectsNotEmpty(fineFluxBox)) {
    //             a_locations[dit].push_back(crseTestBox & fineFluxBox);
    //         }
    //     } // lit over coarsenedFineGrids
    // } // dit over crseGrids

    // if (crseDomain.isPeriodic()) {
    //     UNDEFINED_FUNCTION();
    // }  // end if periodic
}
