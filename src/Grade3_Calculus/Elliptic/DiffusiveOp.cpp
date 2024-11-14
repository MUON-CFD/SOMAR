#include "DiffusiveOp.H"
#include "DiffusiveOpF_F.H"
#include "FABAlgebra.H"
#include "AnisotropicRefinementTools.H"

#include "PoissonOpF_F.H" // TODO: For prolongation functions. Move these to a neutral location.


// NOTE: The nanChecks in applyOp(...) take the most time by far!
#ifndef NDEBUG
#   define nanCheck(x) checkForValidNAN(x)
#else
#   define nanCheck(x)
#endif


namespace Elliptic {


// -----------------------------------------------------------------------------
DiffusiveOp::DiffusiveOp(
    const LevelGeometry&                        a_levGeo,
    const DisjointBoxLayout*                    a_crseGridsPtr,
    const Real                                  a_alpha,
    const Real                                  a_beta,
    const std::vector<Real>&                    a_vKappa,
    const LevelData<FArrayBox>&                 a_eddyNu,
    const std::vector<Real>&                    a_vEddyPrandtl,
    const std::shared_ptr<BCTools::BCFunction>& a_bcFuncPtr,
    const int                                   a_numComps)
: m_alpha(a_alpha)
, m_beta(a_beta)
, m_domain(a_levGeo.getDomain())
, m_grids(a_levGeo.getBoxes())
, m_dXi(a_levGeo.getDXi())
, m_geoSrc(a_levGeo.getGeoSource())
, m_crseAMRGridsPtr(a_crseGridsPtr)
, m_amrCrseDXi(D_DECL(quietNAN, quietNAN, quietNAN))
, m_Jptr(new LevelData<FArrayBox>)
, m_JgupPtr(new LevelData<FluxBox>(a_levGeo.getBoxes(), a_numComps))
, m_invDiagsPtr(new LevelData<FArrayBox>(a_levGeo.getBoxes(), a_numComps))
, m_bcFuncPtr(a_bcFuncPtr)
, m_physBdryIter(a_levGeo.getBoxes())
, m_cfiIter(a_levGeo.getBoxes(), a_levGeo.getCFRegion())
, m_cfInterpPtr()
, m_exCopier()
, m_resPtr(new LevelData<FArrayBox>(a_levGeo.getBoxes(), a_numComps))
{
    // amrCrseDXi
    if (m_crseAMRGridsPtr) {
        const Box&    domBox     = m_domain.domainBox();
        const Box&    crseDomBox = m_crseAMRGridsPtr->physDomain().domainBox();
        const IntVect ref        = calculateRefinementRatio(crseDomBox, domBox);
        CH_assert(ref >= IntVect::Unit);
        CH_assert(ref.product() > 1);
        m_amrCrseDXi = m_dXi * RealVect(ref);
    }

    // m_cfInterpPtr
    if (m_crseAMRGridsPtr) {
        m_cfInterpPtr.reset(new CFInterp);
        m_cfInterpPtr->define(m_grids, m_dXi, *m_crseAMRGridsPtr);
    }

    // m_exCopier
    m_exCopier.exchangeDefine(m_grids, IntVect::Unit);
    m_exCopier.trimEdges(m_grids, IntVect::Unit);

    // J
    aliasLevelData(*m_Jptr,
                   const_cast<LevelData<FArrayBox>*>(&a_levGeo.getCCJ()),
                   Interval(0, 0));

    // Jgup <- beta * (kappa + eddyNu / eddyPrandtl) * Jgup
    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
            FArrayBox&       destFAB   = (*m_JgupPtr)[dit][fcDir];
            const FArrayBox& JgupFAB   = a_levGeo.getFCJgup()[dit][fcDir];
            const FArrayBox& eddyNuFAB = a_eddyNu[dit];

            for (int comp = 0; comp < a_numComps; ++comp) {
                destFAB.copy(JgupFAB, 0, comp, 1);

                if (RealCmp::isZero(a_vEddyPrandtl[comp])) {
                    FABAlgebra::FCmultCC(destFAB,
                                         comp,
                                         destFAB.box(),
                                         eddyNuFAB,
                                         0,
                                         a_beta * a_vKappa[comp],
                                         0.0);  // set eddy diffusivity to zero
                } else {
                    FABAlgebra::FCmultCC(destFAB,
                                         comp,
                                         destFAB.box(),
                                         eddyNuFAB,
                                         0,
                                         a_beta * a_vKappa[comp],
                                         a_beta / a_vEddyPrandtl[comp]);
                }
            } // comp
        } // fcDir
    } // dit

    // invDiags
    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        FArrayBox&       invDiagsFAB = (*m_invDiagsPtr)[dit];
        const Box&       destBox     = invDiagsFAB.box();
        const FArrayBox& JgxxFAB     = (*m_JgupPtr)[dit][0];
        const FArrayBox& JgyyFAB     = (*m_JgupPtr)[dit][1];
        const FArrayBox& JgzzFAB     = (*m_JgupPtr)[dit][SpaceDim - 1];  // Not used in 2D
        const FArrayBox& JFAB        = (*m_Jptr)[dit];

        debugInit(invDiagsFAB);
        FORT_DIFFUSIVEOP_COMPUTEDIAGS(
            CHF_FRA(invDiagsFAB),
            CHF_CONST_FRA(JgxxFAB),
            CHF_CONST_FRA(JgyyFAB),
            CHF_CONST_FRA(JgzzFAB),
            CHF_CONST_FRA1(JFAB, 0),
            CHF_CONST_REALVECT(m_dXi),
            CHF_BOX(destBox),
            CHF_CONST_REAL(m_alpha));

        invDiagsFAB.invert(1.0);
    }
    nanCheck(*m_invDiagsPtr);
}


// -----------------------------------------------------------------------------
void
DiffusiveOp::setAlphaAndBeta(const Real a_alpha, const Real a_beta)
{
    m_alpha = a_alpha;
    m_beta  = a_beta;

    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        FArrayBox&       invDiagsFAB = (*m_invDiagsPtr)[dit];
        const Box&       destBox     = invDiagsFAB.box();
        const FArrayBox& JgxxFAB     = (*m_JgupPtr)[dit][0];
        const FArrayBox& JgyyFAB     = (*m_JgupPtr)[dit][1];
        const FArrayBox& JgzzFAB     = (*m_JgupPtr)[dit][SpaceDim - 1];  // Not used in 2D
        const FArrayBox& JFAB        = (*m_Jptr)[dit];

        debugInit(invDiagsFAB);
        FORT_DIFFUSIVEOP_COMPUTEDIAGS(
            CHF_FRA(invDiagsFAB),
            CHF_CONST_FRA(JgxxFAB),
            CHF_CONST_FRA(JgyyFAB),
            CHF_CONST_FRA(JgzzFAB),
            CHF_CONST_FRA1(JFAB, 0),
            CHF_CONST_REALVECT(m_dXi),
            CHF_BOX(destBox),
            CHF_CONST_REAL(m_alpha));

        invDiagsFAB.invert(1.0);
    }
    nanCheck(*m_invDiagsPtr);
}


// -----------------------------------------------------------------------------
void
DiffusiveOp::applyBCs(StateType&       a_phi,
                      const StateType* a_crsePhiPtr,
                      const Real       a_time,
                      const bool       a_homogPhysBCs,
                      const bool       a_homogCFIBCs) const
{
    CH_assert(a_phi.nComp() == this->numComps());
    CH_assert(a_phi.getBoxes().physDomain() == m_domain);
    CH_assert(a_phi.getBoxes().compatible(m_grids));
    CH_assert(a_phi.ghostVect() == IntVect::Unit);

    // Begin exchange
    a_phi.exchangeBegin(m_exCopier);

    // CFI BCs
    if (!m_cfiIter.isEmpty()) {
        if (a_homogCFIBCs) {
            CFInterp::homogInterpAtCFI(a_phi, m_dXi, m_amrCrseDXi, m_cfiIter);
        } else {
            // if MG depth > 0, we should be calling homogInterpAtCFI!
            CH_assert(a_crsePhiPtr);
            CH_assert(a_crsePhiPtr->getBoxes() == *m_crseAMRGridsPtr);

            nanCheck(*a_crsePhiPtr);
            m_cfInterpPtr->interpAtCFI(a_phi, *a_crsePhiPtr);
        }
    }

    // Phys BCs
    BCTools::applyBC(a_phi,
                     a_homogPhysBCs,
                     *m_bcFuncPtr,
                     a_time,
                     m_geoSrc,
                     m_dXi,
                     m_physBdryIter);

    // Finish exchange
    a_phi.exchangeEnd();  // Bottleneck!

    nanCheck(a_phi);
}


// -----------------------------------------------------------------------------
void
DiffusiveOp::applyOp(StateType&       a_lhs,
                     StateType&       a_phi,
                     const StateType* a_crsePhiPtr,
                     const Real       a_time,
                     const bool       a_homogPhysBCs,
                     const bool       a_homogCFIBCs) const
{
    CH_assert(a_lhs.nComp() == this->numComps());
    CH_assert(a_phi.nComp() == this->numComps());

    this->applyBCs(a_phi, a_crsePhiPtr, a_time, a_homogPhysBCs, a_homogCFIBCs);

    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        FArrayBox&       lhsFAB      = a_lhs[dit];
        const FArrayBox& phiFAB      = a_phi[dit];
        const FArrayBox& JgxxFAB     = (*m_JgupPtr)[dit][0];
        const FArrayBox& JgyyFAB     = (*m_JgupPtr)[dit][1];
        const FArrayBox& JgzzFAB     = (*m_JgupPtr)[dit][SpaceDim - 1];  // Not used in 2D
        const FArrayBox& invDiagsFAB = (*m_invDiagsPtr)[dit];
        const Box        valid       = m_grids[dit];

        debugInit(lhsFAB); // Should do nothing.

        FORT_DIFFUSIVEOP_APPLYOP(
            CHF_FRA(lhsFAB),
            CHF_CONST_FRA(phiFAB),
            CHF_CONST_FRA(JgxxFAB),
            CHF_CONST_FRA(JgyyFAB),
            CHF_CONST_FRA(JgzzFAB),
            CHF_CONST_FRA(invDiagsFAB),
            CHF_CONST_REALVECT(m_dXi),
            CHF_BOX(valid));
    }

    nanCheck(a_lhs);
}


// -----------------------------------------------------------------------------
void
DiffusiveOp::preCond(StateType&       a_phi,
                     const StateType& a_rhs,
                     const Real       a_time,
                     const int        a_relaxIters) const
{
    CH_assert(a_phi.nComp() == this->numComps());
    CH_assert(a_rhs.nComp() == this->numComps());

    this->assignLocal(a_phi, a_rhs);
    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        for (int comp = 0; comp < a_phi.nComp(); ++comp) {
            a_phi[dit].mult((*m_invDiagsPtr)[dit], m_grids[dit], 0, comp, 1);
        }
    }
    this->relax(a_phi, a_rhs, a_time, a_relaxIters);
}


// -----------------------------------------------------------------------------
void
DiffusiveOp::relax(StateType&       a_cor,
                   const StateType& a_res,
                   const Real       a_time,
                   const int        a_relaxIters) const
{
    // this->jacobi_relax(a_cor, a_res, a_time, a_relaxIters);
    this->gsrb_relax(a_cor, a_res, a_time, a_relaxIters);
}


// =========================== MGOperator methods ==============================

// -----------------------------------------------------------------------------
MGOperator<LevelData<FArrayBox>>*
DiffusiveOp::newMGOperator(const IntVect& a_refRatio) const
{
    DiffusiveOp* newOpPtr = nullptr;

    if (a_refRatio == IntVect::Unit) {
        newOpPtr = new DiffusiveOp(*this);
    } else {
        DisjointBoxLayout crseGrids;
        ::coarsen(crseGrids, m_grids, a_refRatio);

        newOpPtr = new DiffusiveOp(*this, crseGrids, a_refRatio);
    }

    return newOpPtr;
}


// -----------------------------------------------------------------------------
void
DiffusiveOp::MGRestrict(StateType&       a_crseRes,
                        const StateType& a_fineRes,
                        const Real       /*a_time*/,
                        const IntVect&   a_refRatio,
                        const MGOpType&  /*a_crseOp*/) const
{
    const DisjointBoxLayout& crseGrids = a_crseRes.getBoxes();
    DataIterator             dit       = a_crseRes.dataIterator();

    nanCheck(a_fineRes);
    for (dit.reset(); dit.ok(); ++dit) {
        // Because the residual is already scaled by J.
        const FArrayBox* JFABPtr = nullptr;

        CFInterp::localCoarsen(a_crseRes[dit],
                               a_fineRes[dit],
                               crseGrids[dit],
                               a_refRatio,
                               false,
                               JFABPtr);
    }
    nanCheck(a_crseRes);
}


// -----------------------------------------------------------------------------
void
DiffusiveOp::MGProlong(StateType&      a_finePhi,
                       StateType&      a_crseCor,
                       const Real      a_time,
                       const IntVect&  a_refRatio,
                       const MGOpType& a_crseOp,
                       const int       a_interpOrder) const
{
    const DisjointBoxLayout& crseGrids = a_crseCor.getBoxes();
    const DisjointBoxLayout& fineGrids = a_finePhi.getBoxes();
    DataIterator             dit       = a_finePhi.dataIterator();

    CH_assert(a_refRatio ==
              fineGrids.physDomain().size() / crseGrids.physDomain().size());

    // This is a local operation.
    CH_assert(crseGrids.compatible(fineGrids));

    const Box refBox(IntVect::Zero, a_refRatio - IntVect::Unit);

    // Do constant interp.
    nanCheck(a_finePhi);
    nanCheck(a_crseCor);
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox&       fineFAB   = a_finePhi[dit];
        const FArrayBox& crseFAB   = a_crseCor[dit];
        const Box&       fineValid = fineGrids[dit];

        const IntVect& fiv = fineValid.smallEnd();
        const IntVect  civ = coarsen(fiv, a_refRatio);

        FORT_MGOPERATOR_PROLONG_CONSTANT(
            CHF_FRA_SHIFT(fineFAB, fiv),
            CHF_CONST_FRA_SHIFT(crseFAB, civ),
            CHF_BOX_SHIFT(fineValid, fiv),
            CHF_CONST_INTVECT(a_refRatio));
    }
    nanCheck(a_finePhi);

    // Upgrade to linear interp, if requested.
    if (a_interpOrder >= 1) {
        nanCheck(a_crseCor);
        a_crseOp.applyBCs(a_crseCor, nullptr, a_time, true, true);
        nanCheck(a_crseCor);

        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox&       fineFAB   = a_finePhi[dit];
            const FArrayBox& crseFAB   = a_crseCor[dit];
            const Box&       crseValid = crseGrids[dit];

            FORT_MGOPERATOR_PROLONG_LINEARUPGRADE(
                CHF_FRA(fineFAB),
                CHF_CONST_FRA(crseFAB),
                CHF_BOX(crseValid),
                CHF_BOX(refBox),
                CHF_CONST_INTVECT(a_refRatio));
        }
        nanCheck(a_finePhi);
    }

    // Upgrade to quadratic interp, if requested.
    while (a_interpOrder >= 2) {
        // Is there room for the stencil?
        bool noRoom = false;
        LayoutIterator lit = crseGrids.layoutIterator();
        for (lit.reset(); lit.ok(); ++lit) {
            const IntVect& size = crseGrids[lit].size();
            if ( !( size >= IntVect(D_DECL(4,4,4)) ) ) {
                noRoom = true;
                break;
            }
        }
        if (noRoom) break;

        // Update with d^/dx^2, d^2/dy^2, and d^2/dz^2 terms.
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox&       fineFAB   = a_finePhi[dit];
            const FArrayBox& crseFAB   = a_crseCor[dit];
            const Box&       crseValid = crseGrids[dit];

            FORT_MGOPERATOR_PROLONG_QUADUPGRADE1(
                CHF_FRA(fineFAB),
                CHF_CONST_FRA(crseFAB),
                CHF_BOX(crseValid),
                CHF_BOX(refBox),
                CHF_CONST_INTVECT(a_refRatio));
        }

        // Fill *ALL* ghosts. Some are already filled, the corners are not.
        //
        // Although this is operating on a coarse MG level, the CornerCopier
        // definition might be expensive. If so, you can cache it.
        for (dit.reset(); dit.ok(); ++dit) {
            BCTools::extrapCorners(a_crseCor[dit], crseGrids[dit], 2);
        }
        CornerCopier ccp(a_crseCor.getBoxes(),
                        a_crseCor.getBoxes(),
                        a_crseCor.getBoxes().physDomain(),
                        a_crseCor.ghostVect(),
                        true);
        a_crseCor.exchange(ccp);

        // Update with mixed second derivative terms.
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox&       fineFAB   = a_finePhi[dit];
            const FArrayBox& crseFAB   = a_crseCor[dit];
            const Box&       crseValid = crseGrids[dit];

            FORT_MGOPERATOR_PROLONG_QUADUPGRADE2(
                CHF_FRA(fineFAB),
                CHF_CONST_FRA(crseFAB),
                CHF_BOX(crseValid),
                CHF_BOX(refBox),
                CHF_CONST_INTVECT(a_refRatio));
        }

        break;
    }

    // Remove the average, if necessary.
    nanCheck(a_finePhi);
    this->removeKernel(a_finePhi);
    nanCheck(a_finePhi);
}


// ============================ Private methods ================================

// -----------------------------------------------------------------------------
DiffusiveOp::DiffusiveOp(const DiffusiveOp& a_srcOp)
: m_alpha(a_srcOp.m_alpha)
, m_beta(a_srcOp.m_beta)
, m_domain(a_srcOp.m_domain)
, m_grids(a_srcOp.m_grids)
, m_dXi(a_srcOp.m_dXi)
, m_geoSrc(a_srcOp.m_geoSrc)
, m_crseAMRGridsPtr(a_srcOp.m_crseAMRGridsPtr)
, m_amrCrseDXi(a_srcOp.m_amrCrseDXi)
, m_Jptr(a_srcOp.m_Jptr)
, m_JgupPtr(a_srcOp.m_JgupPtr)
, m_invDiagsPtr(a_srcOp.m_invDiagsPtr)
, m_bcFuncPtr(a_srcOp.m_bcFuncPtr)
, m_physBdryIter(a_srcOp.m_physBdryIter)
, m_cfiIter(a_srcOp.m_cfiIter)
, m_cfInterpPtr(a_srcOp.m_cfInterpPtr)
, m_exCopier(a_srcOp.m_exCopier)
, m_resPtr(a_srcOp.m_resPtr)
{
}


// -----------------------------------------------------------------------------
DiffusiveOp::DiffusiveOp(const DiffusiveOp&       a_srcOp,
                         const DisjointBoxLayout& a_crseGrids,
                         const IntVect&           a_refRatio)
: m_alpha(a_srcOp.m_alpha)
, m_beta(a_srcOp.m_beta)
, m_domain(a_crseGrids.physDomain())
, m_grids(a_crseGrids)
, m_dXi(a_srcOp.m_dXi * a_refRatio)
, m_geoSrc(a_srcOp.m_geoSrc)
, m_crseAMRGridsPtr(a_srcOp.m_crseAMRGridsPtr)
, m_amrCrseDXi(a_srcOp.m_amrCrseDXi)
, m_Jptr(new LevelData<FArrayBox>(a_crseGrids, 1))
, m_JgupPtr(new LevelData<FluxBox>(a_crseGrids, a_srcOp.numComps()))
, m_invDiagsPtr(new LevelData<FArrayBox>(a_crseGrids, a_srcOp.numComps()))
, m_bcFuncPtr(a_srcOp.m_bcFuncPtr)
, m_physBdryIter(a_crseGrids)
, m_cfiIter(a_crseGrids, CFRegion(a_crseGrids, a_crseGrids.physDomain()))
, m_cfInterpPtr()
, m_exCopier()
, m_resPtr(new LevelData<FArrayBox>(a_crseGrids, a_srcOp.numComps()))
{
    CH_assert(a_srcOp.m_grids.compatible(a_crseGrids));
    CH_assert(a_refRatio >= IntVect::Unit);
    CH_assert(a_refRatio != IntVect::Unit);

    // m_cfInterpPtr remains undefined

    // m_exCopier
    m_exCopier = a_srcOp.m_exCopier;
    ::coarsen(m_exCopier, a_refRatio);

    // J, Jgup
    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        {
            FArrayBox&       destFAB = (*m_Jptr)[dit];
            const FArrayBox& srcFAB  = (*a_srcOp.m_Jptr)[dit];
            CFInterp::localCoarsen(
                destFAB, srcFAB, destFAB.box(), a_refRatio, false, nullptr);
        }
        {
            FluxBox&       destFlub = (*m_JgupPtr)[dit];
            const FluxBox& srcFlub  = (*a_srcOp.m_JgupPtr)[dit];
            CFInterp::localCoarsen(
                destFlub, srcFlub, destFlub.box(), a_refRatio, nullptr);
        }

    }  // dit
    nanCheck(*m_Jptr);
    nanCheck(*m_JgupPtr);

    // invDiags
    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        FArrayBox&       invDiagsFAB = (*m_invDiagsPtr)[dit];
        const Box&       destBox     = invDiagsFAB.box();
        const FArrayBox& JgxxFAB     = (*m_JgupPtr)[dit][0];
        const FArrayBox& JgyyFAB     = (*m_JgupPtr)[dit][1];
        const FArrayBox& JgzzFAB     = (*m_JgupPtr)[dit][SpaceDim - 1];  // Not used in 2D
        const FArrayBox& JFAB        = (*m_Jptr)[dit];

        debugInit(invDiagsFAB);
        FORT_DIFFUSIVEOP_COMPUTEDIAGS(
            CHF_FRA(invDiagsFAB),
            CHF_CONST_FRA(JgxxFAB),
            CHF_CONST_FRA(JgyyFAB),
            CHF_CONST_FRA(JgzzFAB),
            CHF_CONST_FRA1(JFAB, 0),
            CHF_CONST_REALVECT(m_dXi),
            CHF_BOX(destBox),
            CHF_CONST_REAL(m_alpha));

        invDiagsFAB.invert(1.0);
    }
    nanCheck(*m_invDiagsPtr);
}


// -----------------------------------------------------------------------------
void
DiffusiveOp::jacobi_relax(StateType&       a_phi,
                          const StateType& a_rhs,
                          const Real       a_time,
                          const int        a_iters) const
{
    for (int iter = 0; iter < a_iters; ++iter) {
        this->residual(*m_resPtr, a_phi, nullptr, a_rhs, a_time, true, true);

        for (DataIterator dit(m_grids); dit.ok(); ++dit) {
            const Box valid = m_grids[dit];

            FORT_DIFFUSIVEOP_JACOBI(
                CHF_FRA(a_phi[dit]),
                CHF_CONST_FRA((*m_resPtr)[dit]),
                CHF_CONST_FRA((*m_invDiagsPtr)[dit]),
                CHF_BOX(valid));
        }
    }
}


// -----------------------------------------------------------------------------
void
DiffusiveOp::gsrb_relax(StateType&       a_phi,
                        const StateType& a_rhs,
                        const Real       a_time,
                        const int        a_iters) const
{
    for (int iter = 0; iter < a_iters; ++iter) {
        for (int whichPass = 0; whichPass < 2; ++whichPass) {
            // Bottleneck!
            if (whichPass == 0) {
                this->applyBCs(a_phi, nullptr, a_time, true, true);
            } else {
                a_phi.exchange(m_exCopier);
            }

            for (DataIterator dit(m_grids); dit.ok(); ++dit) {
                FArrayBox&       phiFAB  = a_phi[dit];
                const FArrayBox& rhsFAB  = a_rhs[dit];
                const FArrayBox& JgxxFAB = (*m_JgupPtr)[dit][0];
                const FArrayBox& JgyyFAB = (*m_JgupPtr)[dit][1];
                const FArrayBox& JgzzFAB = (*m_JgupPtr)[dit][SpaceDim - 1]; // Not used in 2D
                const FArrayBox& invDiagsFAB = (*m_invDiagsPtr)[dit];
                const Box        destBox     = m_grids[dit];

                FORT_DIFFUSIVEOP_GSRB(
                    CHF_FRA(phiFAB),
                    CHF_CONST_FRA(rhsFAB),
                    CHF_CONST_FRA(JgxxFAB),
                    CHF_CONST_FRA(JgyyFAB),
                    CHF_CONST_FRA(JgzzFAB),
                    CHF_CONST_FRA(invDiagsFAB),
                    CHF_CONST_REALVECT(m_dXi),
                    CHF_BOX(destBox),
                    CHF_CONST_INT(whichPass));
            }  // dit
        } // whichPass (red or black)
    } // iter
}


}; // namespace Elliptic

