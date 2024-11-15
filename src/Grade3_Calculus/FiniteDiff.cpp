#include "FiniteDiff.H"
#include "FiniteDiffF_F.H"
#include "FABAlgebra.H"
#include "SOMAR_Constants.H"
#include "Debug.H"
#include "IO.H"
#include "Convert.H"
#include "AnisotropicFluxRegister.H"
#include "AnisotropicRefinementTools.H"
#include "BCTools.H"


// -----------------------------------------------------------------------------
// Default constructor
// -----------------------------------------------------------------------------
FiniteDiff::FiniteDiff ()
:m_levGeoPtr(nullptr)
{
}


// -----------------------------------------------------------------------------
// Full constructor
// -----------------------------------------------------------------------------
FiniteDiff::FiniteDiff (const LevelGeometry& a_levGeo)
:m_levGeoPtr(nullptr)
{
    this->define(a_levGeo);
}


// -----------------------------------------------------------------------------
// Full virtual constructor
// -----------------------------------------------------------------------------
void
FiniteDiff::define (const LevelGeometry& a_levGeo)
{
    CH_assert(a_levGeo.hasGrids());
    m_levGeoPtr = &a_levGeo;
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
FiniteDiff::~FiniteDiff ()
{
    this->clear();
}


// -----------------------------------------------------------------------------
// "Destructor"
// -----------------------------------------------------------------------------
void
FiniteDiff::clear()
{
}


// -----------------------------------------------------------------------------
// Static utility
// Computes the staggered derivative of a scalar onto a_pdBox.
// This funciton does not fill ghosts.
//
// a_pd      : The dest holder.
// a_pdComp  : The comp of a_pd that will be filled.
// a_pdBox   : Where to compute the derivs.
//             Must have same centering as a_pd.
// a_phi     : The source data. Ghosts must be filled if needed.
//             Must have the opposite centering of a_pd in a_derivDir.
// a_phiComp : The comp of a_phi that will be differentiated.
// a_derivDir: The differencing direction.
// a_dXi     : The grid spacing in a_derivDir.
// -----------------------------------------------------------------------------
void
FiniteDiff::partialD(FArrayBox&       a_pd,
                     const int        a_pdComp,
                     const Box&       a_pdBox,
                     const FArrayBox& a_phi,
                     const int        a_phiComp,
                     const int        a_derivDir,
                     const Real       a_dXi,
                     const bool       a_addToPD)
{
    // Sanity checks.
    CH_assert(0 <= a_derivDir && a_derivDir < SpaceDim);

    CH_assert(a_pd.box().type() == a_pdBox.type());
    CH_assert(a_pd.box().contains(a_pdBox));
    CH_assert(0 <= a_pdComp && a_pdComp < a_pd.nComp());

    CH_assert(a_phi.box().type(a_derivDir) == 1 - a_pdBox.type(a_derivDir));
    CH_assert(0 <= a_phiComp && a_phiComp < a_phi.nComp());

    if (a_addToPD) {
        if (a_pdBox.type(a_derivDir) == 1) {
            FORT_FINITEDIFF_ADDPARTIALD_CC2NC (
                CHF_FRA1(a_pd, a_pdComp),
                CHF_CONST_FRA1(a_phi, a_phiComp),
                CHF_CONST_INT(a_derivDir),
                CHF_CONST_REAL(a_dXi),
                CHF_BOX(a_pdBox));
        } else {
            FORT_FINITEDIFF_ADDPARTIALD_NC2CC (
                CHF_FRA1(a_pd, a_pdComp),
                CHF_CONST_FRA1(a_phi, a_phiComp),
                CHF_CONST_INT(a_derivDir),
                CHF_CONST_REAL(a_dXi),
                CHF_BOX(a_pdBox));
        }
    } else {
        if (a_pdBox.type(a_derivDir) == 1) {
            FORT_FINITEDIFF_PARTIALD_CC2NC (
                CHF_FRA1(a_pd, a_pdComp),
                CHF_CONST_FRA1(a_phi, a_phiComp),
                CHF_CONST_INT(a_derivDir),
                CHF_CONST_REAL(a_dXi),
                CHF_BOX(a_pdBox));
        } else {
            FORT_FINITEDIFF_PARTIALD_NC2CC (
                CHF_FRA1(a_pd, a_pdComp),
                CHF_CONST_FRA1(a_phi, a_phiComp),
                CHF_CONST_INT(a_derivDir),
                CHF_CONST_REAL(a_dXi),
                CHF_BOX(a_pdBox));
        }
    }
}


// -----------------------------------------------------------------------------
void
FiniteDiff::computeSlopes(FArrayBox&       a_pd,
                          const int        a_pdComp,
                          const FArrayBox& a_phi,
                          const int        a_phiComp,
                          const Box&       a_validPhiBox,
                          const int        a_derivDir,
                          const Real       a_dXi)
{
    // Sanity checks.
    CH_assert(0 <= a_derivDir && a_derivDir < SpaceDim);

    CH_assert(0 <= a_pdComp && a_pdComp < a_pd.nComp());
    CH_assert(0 <= a_phiComp && a_phiComp < a_phi.nComp());

    CH_assert(a_pd.box().type() == a_phi.box().type());
    CH_assert(a_phi.box().type() == a_validPhiBox.type());

    CH_assert(a_validPhiBox.contains(a_pd.box()));

    Box loBox, centerBox, hiBox;
    if (a_validPhiBox.type(a_derivDir) == IndexType::CELL) {
        loBox     = adjCellLo(a_validPhiBox, a_derivDir, -1) & a_pd.box();
        centerBox = grow(a_validPhiBox, -BASISV(a_derivDir)) & a_pd.box();
        hiBox     = adjCellHi(a_validPhiBox, a_derivDir, -1) & a_pd.box();
    } else {
        loBox     = bdryLo(a_validPhiBox, a_derivDir, 1)     & a_pd.box();
        centerBox = grow(a_validPhiBox, -BASISV(a_derivDir)) & a_pd.box();
        hiBox     = bdryHi(a_validPhiBox, a_derivDir, 1)     & a_pd.box();
    }

    CH_assert(loBox.isEmpty() || a_phi.box().contains(loBox));
    CH_assert(a_phi.box().contains(centerBox));
    CH_assert(hiBox.isEmpty() || a_phi.box().contains(hiBox));

    FORT_FINITEDIFF_SLOPES(
        CHF_FRA1(a_pd, a_pdComp),
        CHF_CONST_FRA1(a_phi, a_phiComp),
        CHF_CONST_INT(a_derivDir),
        CHF_CONST_REAL(a_dXi),
        CHF_BOX(loBox),
        CHF_BOX(centerBox),
        CHF_BOX(hiBox)
    );
}


// -----------------------------------------------------------------------------
void
FiniteDiff::addSecondDifference(FArrayBox&       a_pd,
                                const int        a_pdComp,
                                const FArrayBox& a_phi,
                                const int        a_phiComp,
                                const Box&       a_validPhiBox,
                                const int        a_derivDir,
                                const Real       a_scale)
{
    // Sanity checks.
    CH_assert(0 <= a_derivDir && a_derivDir < SpaceDim);

    CH_assert(0 <= a_pdComp && a_pdComp < a_pd.nComp());
    CH_assert(0 <= a_phiComp && a_phiComp < a_phi.nComp());

    CH_assert(a_pd.box().type() == a_phi.box().type());
    CH_assert(a_phi.box().type() == a_validPhiBox.type());

    CH_assert(a_validPhiBox.contains(a_pd.box()));

    Box loBox, centerBox, hiBox;
    if (a_validPhiBox.type(a_derivDir) == IndexType::CELL) {
        loBox     = adjCellLo(a_validPhiBox, a_derivDir, -1) & a_pd.box();
        centerBox = grow(a_validPhiBox, -BASISV(a_derivDir)) & a_pd.box();
        hiBox     = adjCellHi(a_validPhiBox, a_derivDir, -1) & a_pd.box();
    } else {
        loBox     = bdryLo(a_validPhiBox, a_derivDir, 1)     & a_pd.box();
        centerBox = grow(a_validPhiBox, -BASISV(a_derivDir)) & a_pd.box();
        hiBox     = bdryHi(a_validPhiBox, a_derivDir, 1)     & a_pd.box();
    }

    CH_assert(loBox.isEmpty() || a_phi.box().contains(loBox));
    CH_assert(a_phi.box().contains(centerBox));
    CH_assert(hiBox.isEmpty() || a_phi.box().contains(hiBox));

    FORT_FINITEDIFF_ADDSECONDDIFFERENCE(
        CHF_FRA1(a_pd, a_pdComp),
        CHF_CONST_FRA1(a_phi, a_phiComp),
        CHF_CONST_INT(a_derivDir),
        CHF_BOX(loBox),
        CHF_BOX(centerBox),
        CHF_BOX(hiBox),
        CHF_CONST_REAL(a_scale)
    );
}


// -----------------------------------------------------------------------------
void
FiniteDiff::addLaplacian(FArrayBox&       a_pd,
                         const int        a_pdComp,
                         const FArrayBox& a_phi,
                         const int        a_phiComp,
                         const Box&       a_validPhiBox,
                         const RealVect&  a_scale)
{
    for (int derivDir = 0; derivDir < SpaceDim; ++derivDir) {
        if (RealCmp::isZero(a_scale[derivDir])) continue;
        FiniteDiff::addSecondDifference(a_pd,
                                        a_pdComp,
                                        a_phi,
                                        a_phiComp,
                                        a_validPhiBox,
                                        derivDir,
                                        a_scale[derivDir]);
    }
}


// -----------------------------------------------------------------------------
// Computes the FC, single-level gradient of a CC scalar.
//
// This function assumes:
//   1. All BCs have been set.
//   2. *this, a_gradPhi, and a_phi are defined over the same grids.
//   3. a_gradPhi and a_phi have the same number of components.
// -----------------------------------------------------------------------------
void
FiniteDiff::levelGradientMAC (LevelData<FluxBox>&         a_gradPhi,
                              const LevelData<FArrayBox>& a_phi) const
{
    // Sanity checks
    CH_assert(this->isDefined());
    CH_assert(m_levGeoPtr->getBoxes().compatible(a_phi.getBoxes()));
    CH_assert(m_levGeoPtr->getBoxes().compatible(a_gradPhi.getBoxes()));
    CH_assert(a_phi.nComp() == a_gradPhi.nComp());

    // Collect references
    const int                 numPhiComps = a_phi.nComp();
    const LevelData<FluxBox>& Jgup        = m_levGeoPtr->getFCJgup();
    const RealVect&           dXi         = m_levGeoPtr->getDXi();
    const DisjointBoxLayout&  grids       = m_levGeoPtr->getBoxes();
    DataIterator              dit         = grids.dataIterator();

    for (dit.reset(); dit.ok(); ++dit) {
        const FArrayBox& phiFAB = a_phi[dit];

        for (int gradDir = 0; gradDir < SpaceDim; ++gradDir) {
            // Gather references
            FArrayBox&       gradFAB = a_gradPhi[dit][gradDir];
            const FArrayBox& JgupFAB = Jgup[dit][gradDir];
            const Real       dXiDir  = dXi[gradDir];

            // Compute destination box.
            Box fcValid = grids[dit];
            fcValid.surroundingNodes(gradDir);
            fcValid &= gradFAB.box();

            // Compute derivs in each direction.
            for (int phiComp = 0; phiComp < numPhiComps; ++phiComp) {
                FiniteDiff::partialD(gradFAB,
                                     phiComp,
                                     fcValid,
                                     phiFAB,
                                     phiComp,
                                     gradDir,
                                     dXiDir);
                gradFAB.mult(JgupFAB, 0, phiComp, 1);
            } // phiComp
        } // gradDir
    } // dit
}


// -----------------------------------------------------------------------------
// Computes the CC, single-level divergence of FC fluxes.
//
// This function assumes:
//   1. all BCs have been set.
//   2. The fluxes are scaled by J.
//
// a_scaleByJinv: A true divergence will scale the result by 1/J, but
//    there are some algorithms that avoid this for efficiency reasons.
//
// a_scale: These will be multiplied into each partial derivative.
//    For example, if you wish to compute a horizontal divergence in 3D,
//    set this to (1.0, 1.0, 0.0).
// -----------------------------------------------------------------------------
void
FiniteDiff::levelDivergenceMAC(LevelData<FArrayBox>&     a_div,
                               const LevelData<FluxBox>& a_flux,
                               const bool                a_scaleByJinv,
                               const RealVect&           a_scale) const
{
    // Sanity checks
    CH_assert(this->isDefined());
    CH_assert(a_div.getBoxes().compatible(a_flux.getBoxes()));
    CH_assert(a_div.nComp() == a_flux.nComp());

    // Collect references
    const int                divComps = a_div.nComp();
    const RealVect&          dXi      = m_levGeoPtr->getDXi();
    const DisjointBoxLayout& grids    = a_div.getBoxes();
    DataIterator             dit      = a_div.dataIterator();

    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox&     divFAB = a_div[dit];
        const FluxBox& fluxFB = a_flux[dit];
        const Box&     valid  = grids[dit];
        CH_assert(fluxFB.box().contains(valid));

        if (SpaceDim == 2) {
            FORT_FINITEDIFF_DIV2D (
                CHF_FRA(divFAB),
                CHF_CONST_FRA(fluxFB[0]),
                CHF_CONST_FRA(fluxFB[1]),
                CHF_BOX(valid),
                CHF_CONST_REALVECT(dXi),
                CHF_CONST_REALVECT(a_scale));
        } else {
            FORT_FINITEDIFF_DIV3D (
                CHF_FRA(divFAB),
                CHF_CONST_FRA(fluxFB[0]),
                CHF_CONST_FRA(fluxFB[1]),
                CHF_CONST_FRA(fluxFB[2]),
                CHF_BOX(valid),
                CHF_CONST_REALVECT(dXi),
                CHF_CONST_REALVECT(a_scale));
        }

        // Scale by Jinv, if requested.
        if (a_scaleByJinv) {
            CH_assert(m_levGeoPtr->hasGrids());
            CH_assert(m_levGeoPtr->getBoxes().compatible(grids));

            const FArrayBox& JinvFAB = m_levGeoPtr->getCCJinv()[dit];
            CH_assert(JinvFAB.box().contains(valid));

            for (int comp = 0; comp < divComps; ++comp) {
                divFAB.mult(JinvFAB, valid, 0, comp, 1);
            }
        } // end if a_scaleByJinv
    } // end loop over grids (dit)
}


// -----------------------------------------------------------------------------
void
FiniteDiff::levelDivergenceMAC_4thOrder(LevelData<FArrayBox>&     a_div,
                                        const LevelData<FluxBox>& a_flux,
                                        const bool                a_scaleByJinv,
                                        const RealVect&           a_scale) const
{
    // Sanity checks
    CH_assert(this->isDefined());
    CH_assert(a_div.getBoxes().compatible(a_flux.getBoxes()));
    CH_assert(a_div.nComp() == a_flux.nComp());

    // Collect references
    const RealVect&          dXi   = m_levGeoPtr->getDXi();
    const DisjointBoxLayout& grids = a_div.getBoxes();

    LevelData<FArrayBox> deriv2(grids, SpaceDim, IntVect::Unit);

    for (int divComp = 0; divComp < a_div.nComp(); ++divComp) {
        for (DataIterator dit(grids); dit.ok(); ++dit) {
            for (int derivDir = 0; derivDir < SpaceDim; ++derivDir) {
                FArrayBox&       derivFAB = deriv2[dit];
                const FArrayBox& fluxFAB  = a_flux[dit][derivDir];
                const Box        destBox  = grids[dit];
                const Real       dXiDir   = dXi[derivDir] / a_scale[derivDir];
                constexpr bool   accum    = false;

                FiniteDiff::partialD(derivFAB,
                                     derivDir,
                                     destBox,
                                     fluxFAB,
                                     divComp,
                                     derivDir,
                                     dXiDir,
                                     accum);
            }  // derivDir
        } // end loop over grids (dit)

        // Fill ghosts.
        BCTools::extrapAllGhosts(deriv2, 1);
        deriv2.exchange();

        // Compute divergence
        for (DataIterator dit(grids); dit.ok(); ++dit) {
            FArrayBox&       divFAB   = a_div[dit];
            const FArrayBox& derivFAB = deriv2[dit];
            const Box        destBox  = grids[dit];

            divFAB.setVal(0.0, divComp);
            for (int derivDir = 0; derivDir < SpaceDim; ++derivDir) {
                const IntVect e = BASISV(derivDir);
                for (BoxIterator bit(destBox); bit.ok(); ++bit) {
                    const IntVect& cc = bit();
                    divFAB(cc, divComp) += (-        derivFAB(cc + e, derivDir)
                                            + 26.0 * derivFAB(cc    , derivDir)
                                            -        derivFAB(cc - e, derivDir))
                                         / 24.0;
                    // divFAB(cc, divComp) += derivFAB(cc, derivDir);
                }  // bit
            } // derivDir
        } // dit

        // Scale by Jinv, if requested.
        if (a_scaleByJinv) {
            for (DataIterator dit(grids); dit.ok(); ++dit) {
                FArrayBox& divFAB  = a_div[dit];
                const Box  destBox = grids[dit];

                CH_assert(m_levGeoPtr->hasGrids());
                CH_assert(m_levGeoPtr->getBoxes().compatible(grids));

                const FArrayBox& JinvFAB = m_levGeoPtr->getCCJinv()[dit];
                CH_assert(JinvFAB.box().contains(destBox));

                divFAB.mult(JinvFAB, destBox, 0, divComp, 1);
            } // dit
        } // end if a_scaleByJinv
    } // divComp
}


// -----------------------------------------------------------------------------
void
FiniteDiff::levelVectorGradient(StaggeredFluxLD&          a_grad,
                                const LevelData<FluxBox>& a_cartVect) const
{
    CH_assert(this->isDefined());
    CH_assert(a_cartVect.nComp() == 1);
    CH_assert(a_cartVect.ghostVect() >= IntVect::Unit);
    CH_assert(m_levGeoPtr->getBoxes() == a_cartVect.getBoxes());
    CH_assert(m_levGeoPtr->getBoxes() == a_grad.getBoxes());

    debugCheckValidFaceOverlap(a_cartVect);
    debugCheckValidFaceOverlap(m_levGeoPtr->getFCJgup());

    const RealVect&          dXi      = m_levGeoPtr->getDXi();
    const DisjointBoxLayout& grids    = m_levGeoPtr->getBoxes();
    DataIterator             dit      = grids.dataIterator();
    int                      i, j;

    for (dit.reset(); dit.ok(); ++dit) {
        const Box& ccValid = grids[dit];

        for (i = 0; i < SpaceDim; ++i) {   // i = derivDir
            const FArrayBox& fcJgiiFAB = m_levGeoPtr->getFCJgup()[dit][i];

            // Diagonal element is CC.
            {
                j = i;

                FArrayBox&       gradFAB = a_grad[i][j][dit];
                const FArrayBox& vFAB    = a_cartVect[dit][j];
                const Box        pdBox   = grids[dit].grow(j, 1);

                CH_assert(gradFAB.box().contains(pdBox));
                debugInit(gradFAB);

                FiniteDiff::partialD(gradFAB, 0, pdBox, vFAB, 0, i, dXi[i]);

                // Multiply by metric.
                FArrayBox ccJgiiFAB(pdBox, 1);
                Convert::Simple(ccJgiiFAB, fcJgiiFAB);
                gradFAB.mult(ccJgiiFAB);
            }

            // Off diagonal element is nodal in both i and j dirs.
            {
                j = (i + 1) % SpaceDim;

                FArrayBox&       gradFAB = a_grad[i][j][dit];
                const FArrayBox& vFAB    = a_cartVect[dit][j];
                const Box pdBox = grids[dit].surroundingNodes(i)
                                            .surroundingNodes(j);

                CH_assert(gradFAB.box().contains(pdBox));
                debugInit(gradFAB);

                FiniteDiff::partialD(gradFAB, 0, pdBox, vFAB, 0, i, dXi[i]);

                // Multiply by metric.
                FArrayBox ecJgiiFAB(pdBox, 1);
                Convert::Simple(ecJgiiFAB, fcJgiiFAB);
                gradFAB.mult(ecJgiiFAB);
            }

            // Off diagonal element is nodal in both i and j dirs.
            if (SpaceDim > 2) {
                j = (i + 2) % SpaceDim;

                FArrayBox&       gradFAB = a_grad[i][j][dit];
                const FArrayBox& vFAB    = a_cartVect[dit][j];
                const Box pdBox = grids[dit].surroundingNodes(i)
                                            .surroundingNodes(j);

                CH_assert(gradFAB.box().contains(pdBox));
                debugInit(gradFAB);

                FiniteDiff::partialD(gradFAB, 0, pdBox, vFAB, 0, i, dXi[i]);

                // Multiply by metric.
                FArrayBox ecJgiiFAB(pdBox, 1);
                Convert::Simple(ecJgiiFAB, fcJgiiFAB);
                gradFAB.mult(ecJgiiFAB);
            }
        }
    } // end loop over grids (dit)

    // At this point, a_grad has valid data everywhere it is defined.

    // If a_cartVect was exchanged, we shouldn't need to exchange a_grad, but for now..
    // a_grad.exchange();
}


// -----------------------------------------------------------------------------
void
FiniteDiff::levelVectorDivergence(LevelData<FluxBox>&    a_div,
                                  const StaggeredFluxLD& a_flux,
                                  const bool             a_scaleByJinv,
                                  const bool             a_addToDiv,
                                  const Real             a_scale) const
{
    CH_assert(this->isDefined());
    CH_assert(a_div.nComp() == 1);
    CH_assert(m_levGeoPtr->getBoxes() == a_div.getBoxes());
    CH_assert(m_levGeoPtr->getBoxes() == a_flux.getBoxes());
    CH_assert(!RealCmp::isZero(a_scale));

    const RealVect&          dXi      = m_levGeoPtr->getDXi();
    const DisjointBoxLayout& grids    = m_levGeoPtr->getBoxes();

    for (DataIterator dit(grids); dit.ok(); ++dit) {
        for (int divComp = 0; divComp < SpaceDim; ++divComp) {
            FArrayBox& divFAB = a_div[dit][divComp];
            const Box destBox = grids[dit].surroundingNodes(divComp);

            if (!a_addToDiv) {
                divFAB.setVal(0.0);
            }

            for (int derivDir = 0; derivDir < SpaceDim; ++derivDir) {
                const FArrayBox& fluxFAB = a_flux[derivDir][divComp][dit];
                const bool       accum   = true;
                const Real       dXiDir  = dXi[derivDir] / a_scale;

                FiniteDiff::partialD(
                    divFAB, 0, destBox, fluxFAB, 0, derivDir, dXiDir, accum);
            }
        } // divComp
    } // end loop over grids (dit)

    // Scale by Jinv, if requested.
    if (a_scaleByJinv) {
        m_levGeoPtr->divByJ(a_div);
    }
}


// -----------------------------------------------------------------------------
void
FiniteDiff::levelVectorDivergence_4thOrder(LevelData<FluxBox>&    a_div,
                                           const StaggeredFluxLD& a_flux,
                                           const bool             a_scaleByJinv,
                                           const bool             a_addToDiv,
                                           const Real             a_scale) const
{
    CH_assert(this->isDefined());
    CH_assert(a_div.nComp() == 1);
    CH_assert(m_levGeoPtr->getBoxes() == a_div.getBoxes());
    CH_assert(m_levGeoPtr->getBoxes() == a_flux.getBoxes());
    CH_assert(!RealCmp::isZero(a_scale));

    const RealVect&          dXi   = m_levGeoPtr->getDXi();
    const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();

    // Compute 2nd-order, FC derivatives of CC and EC fluxes.
    // Format: deriv2[dit][velComp](fciv, derivDir)
    LevelData<FluxBox> deriv2(grids, SpaceDim, IntVect::Unit);
    for (DataIterator dit(grids); dit.ok(); ++dit) {
        for (int velComp = 0; velComp < SpaceDim; ++velComp) { // = fc dir
            for (int derivDir = 0; derivDir < SpaceDim; ++derivDir) {
                FArrayBox&       derivFAB = deriv2[dit][velComp];
                const FArrayBox& fluxFAB  = a_flux[derivDir][velComp][dit];
                const Real       dXiDir   = dXi[derivDir] / a_scale;
                const Box        destBox  = grids[dit].surroundingNodes(velComp);
                constexpr bool   accum    = false;

                FiniteDiff::partialD(derivFAB,
                                     derivDir,
                                     destBox,
                                     fluxFAB,
                                     0,
                                     derivDir,
                                     dXiDir / a_scale,
                                     accum);
            }  // derivDir
        } // velComp
    } // dit

    // Fill ghosts.
    // Linear extrap leads to a 2nd-order derivative at the boundaries.
    BCTools::extrapAllGhosts(deriv2, 1);
    deriv2.exchange();

    // Upgrade to 4th order derivatives.
    LevelData<FluxBox> deriv4(grids, SpaceDim);  // Don't need this. Just add to div.
    for (DataIterator dit(grids); dit.ok(); ++dit) {
        for (int velComp = 0; velComp < SpaceDim; ++velComp) { // = fc dir
            FArrayBox&       d4FAB   = deriv4[dit][velComp];
            const FArrayBox& d2FAB   = deriv2[dit][velComp];
            const Box        destBox = grids[dit].surroundingNodes(velComp);

            for (int derivDir = 0; derivDir < SpaceDim; ++derivDir) {
                const IntVect e = BASISV(derivDir);
                for (BoxIterator bit(destBox); bit.ok(); ++bit) {
                    const IntVect& fc = bit();
                    d4FAB(fc, derivDir) = (-        d2FAB(fc + e, derivDir)
                                           + 26.0 * d2FAB(fc    , derivDir)
                                           -        d2FAB(fc - e, derivDir))
                                        / 24.0;
                }  // bit
            } // derivDir
        } // velComp
    } // dit

    // Compute divergence.
    for (DataIterator dit(grids); dit.ok(); ++dit) {
        for (int velComp = 0; velComp < SpaceDim; ++velComp) { // = fc dir
            FArrayBox&       divFAB   = a_div[dit][velComp];
            const FArrayBox& derivFAB = deriv4[dit][velComp];
            const Box        destBox = grids[dit].surroundingNodes(velComp);

            if (!a_addToDiv) divFAB.setVal(0.0);

            D_TERM(
            divFAB.plus(derivFAB, destBox, 0, 0, 1);,
            divFAB.plus(derivFAB, destBox, 1, 0, 1);,
            divFAB.plus(derivFAB, destBox, 2, 0, 1);)
        } // velComp
    } // dit

    // Scale by Jinv, if requested.
    if (a_scaleByJinv) {
        m_levGeoPtr->divByJ(a_div);
    }
}


// -----------------------------------------------------------------------------
void
FiniteDiff::compGradientMAC(LevelData<FluxBox>&         a_gradPhi,
                            const LevelData<FArrayBox>& a_phi) const
{
    FiniteDiff::levelGradientMAC(a_gradPhi, a_phi);
}


// -----------------------------------------------------------------------------
void
FiniteDiff::compDivergenceMAC(LevelData<FArrayBox>&     a_div,
                              const LevelData<FluxBox>& a_flux,
                              const LevelData<FluxBox>& a_fineFlux,
                              const bool                a_scaleByJinv,
                              const RealVect&           a_scale) const
{
    FiniteDiff::levelDivergenceMAC(a_div, a_flux, a_scaleByJinv, a_scale);

    const DisjointBoxLayout& grids  = m_levGeoPtr->getBoxes();
    const ProblemDomain&     domain = m_levGeoPtr->getDomain();
    const Box&               domBox = m_levGeoPtr->getDomainBox();

    const DisjointBoxLayout& fineGrids  = a_fineFlux.getBoxes();
    const Box&               fineDomBox = fineGrids.physDomain().domainBox();

    const IntVect   refRatio = calculateRefinementRatio(domBox, fineDomBox);
    const RealVect& dXi      = m_levGeoPtr->getDXi();
    const int       numComps = a_div.nComp();
    const Interval  ivl(0, numComps - 1);

    // 1. Initialize flux register.
    AnisotropicFluxRegister reg(fineGrids, grids, domain, refRatio, numComps);
    reg.setToZero();

    // 2. Increment the coarse side.
    for (DataIterator ditc(grids); ditc.ok(); ++ditc) {
        for (int idir = 0; idir < SpaceDim; ++idir) {
            const FArrayBox& crseFlux = a_flux[ditc][idir];
            const Real       scale    = 1.0 / dXi[idir];
            reg.incrementCoarse(crseFlux, scale, ditc(), ivl, ivl, idir);
        }
    }

    // 3. Increment the fine side.
    for (DataIterator ditf(fineGrids); ditf.ok(); ++ditf) {
        for (int idir = 0; idir < SpaceDim; ++idir) {
            const FArrayBox& fineFlux = a_fineFlux[ditf][idir];
            const Real       scale = 1.0 / dXi[idir];  // Yes, the coarse dXi!
            reg.incrementFine(fineFlux, scale, ditf(), ivl, ivl, idir);
        }  // idir
    } // ditf

    // 4. Reflux.
    reg.reflux(a_div);
}
