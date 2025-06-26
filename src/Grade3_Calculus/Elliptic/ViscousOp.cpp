#include "ViscousOp.H"
#include "ViscousOpF_F.H"
#include "Subspace.H"
#include "Convert.H"
#include "Debug.H"
#include "SOMAR_Constants.H"
#include "MiscUtils.H"
// #include"Comm.H"
#include "CFInterp.H"
#include "CFInterpF_F.H"
#include "AnisotropicRefinementTools.H"
#include "FABAlgebra.H"


namespace Elliptic {

// NOTE: The nanChecks in applyOp(...) take the most time by far!
#ifndef NDEBUG
#   define nanCheck(x) checkForValidNAN(x)
#else
#   define nanCheck(x)
#endif


// -----------------------------------------------------------------------------
ViscousOp::ViscousOp(const LevelGeometry*                        a_levGeoPtr,
                     const Real                                  a_alpha,
                     const Real                                  a_beta,
                     const LevelData<FArrayBox>&                 a_ccNu,
                     const std::shared_ptr<BCTools::BCFunction>& a_bcFuncPtr)
: m_alpha(a_alpha)
, m_beta(a_beta)
, m_domain(a_levGeoPtr->getDomain())
, m_grids(a_levGeoPtr->getBoxes())
, m_dXi(a_levGeoPtr->getDXi())
, m_amrCrseDXi(D_DECL(quietNAN, quietNAN, quietNAN))
, m_geoSrc(a_levGeoPtr->getGeoSource())
, m_Jptr(new LevelData<FluxBox>(a_levGeoPtr->getBoxes(), 1, IntVect::Unit))
, m_nuJgupPtr(new StaggeredFluxLD(a_levGeoPtr->getBoxes()))
, m_invDiagsPtr(new LevelData<FluxBox>(a_levGeoPtr->getBoxes(), 1))
, m_bcFuncPtr(a_bcFuncPtr)
, m_physBdryIter(a_levGeoPtr->getBoxes())
// m_cfInterp
// exchange copiers
// m_interiorBox
// m_bdryBox
, m_resPtr(new LevelData<FluxBox>(a_levGeoPtr->getBoxes(), 1))
{
    // m_amrCrseDXi
    if (a_levGeoPtr->getCrseGridsPtr()) {
        m_amrCrseDXi = m_dXi * a_levGeoPtr->getCrseRefRatio();
    }

    // m_cfRegion
    for (int velComp = 0; velComp < SpaceDim; ++velComp) {
        m_cfRegion[velComp].define(m_grids, m_domain, BASISV(velComp));
    }

    // Exchange copiers
    for (int d = 0; d < SpaceDim; ++d) {
        constexpr IntVect ghostVect(D_DECL(1, 1, 1));
        constexpr bool doValidCorners = true;

        m_exCopier[d].defineValidExchange(m_grids, ghostVect, d, doValidCorners);
        m_exCornerCopier1[d].defineInvalidCornerExchange1(m_grids, ghostVect, d);
        m_exCornerCopier2[d].defineInvalidCornerExchange2(m_grids, ghostVect, d);
    }

    // Where L[u] needs to be computed.
    for (int velComp = 0; velComp < SpaceDim; ++velComp) {
        m_interiorBox[velComp] = m_domain.domainBox();
        m_interiorBox[velComp].surroundingNodes(velComp);
        if (!m_domain.isPeriodic(velComp)) {
            m_bdryBox[velComp][0] = bdryLo(m_interiorBox[velComp], velComp);
            m_bdryBox[velComp][1] = bdryHi(m_interiorBox[velComp], velComp);
            m_interiorBox[velComp].grow(velComp, -1);
        }
    }

    // Metric
    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        const FArrayBox& ccNuFAB = a_ccNu[dit];
        checkForNAN(ccNuFAB, ccNuFAB.box());

        for (int velComp = 0; velComp < SpaceDim; ++velComp) {
            // FC J
            {
                auto& destFAB = (*m_Jptr)[dit][velComp];
                Convert::Simple(destFAB, a_levGeoPtr->getCCJ()[dit]);
            }

            // CC nu*Jgup
            {
                auto& destFAB = (*m_nuJgupPtr)[velComp][velComp][dit];
                Convert::Simple(destFAB, a_levGeoPtr->getFCJgup()[dit][velComp]);
                destFAB.mult(ccNuFAB);
            }

            // EC nu*Jgup
            for (int derivDir = 0; derivDir < SpaceDim; ++derivDir) {
                if (derivDir == velComp) continue;

                const Box ecRegion = m_grids[dit].surroundingNodes(velComp)
                                                 .surroundingNodes(derivDir);

                auto& destFAB = (*m_nuJgupPtr)[velComp][derivDir][dit];
                Convert::Simple(destFAB, a_levGeoPtr->getFCJgup()[dit][derivDir]);
                FABAlgebra::ECmultCC(destFAB, 0, ecRegion, ccNuFAB, 0);
            }
        }
    }

    BCTools::extrapAllGhosts(*m_Jptr, 2);
    m_Jptr->exchange();
    debugCheckValidFaceOverlap(*m_Jptr);

    for (int velComp = 0; velComp < SpaceDim; ++velComp) {
        for (int derivDir = 0; derivDir < SpaceDim; ++derivDir) {
            BCTools::extrapAllGhosts((*m_nuJgupPtr)[velComp][derivDir], 2);
        }
    }
    m_nuJgupPtr->exchange();


    // Compute diagonals, etc.
    this->cacheMatrixElements();

    // CFInterp
    if (a_levGeoPtr->getCrseGridsPtr()) {
        m_cfInterp.define(*a_levGeoPtr, *a_levGeoPtr->getCrseGridsPtr());
    }
}


// -----------------------------------------------------------------------------
void
ViscousOp::setAlphaAndBeta(const Real a_alpha, const Real a_beta)
{
    m_alpha = a_alpha;
    m_beta  = a_beta;
    this->cacheMatrixElements();
}


// -----------------------------------------------------------------------------
void
ViscousOp::applyBCs(StateType&       a_vel,
                    const StateType* a_crseVelPtr,
                    const Real       a_time,
                    const bool       a_homogPhysBCs,
                    const bool       a_homogCFIBCs) const
{
    checkForValidNAN(a_vel);
    // debugCheckValidFaceOverlap(a_vel);
    LayoutTools::averageOverlappingValidFaces(a_vel);

    BCTools::extrapAllGhosts(a_vel, 2);

    // // Begin valid exchange (See StaggeredCopier::defineValidExchange for details.)
    // std::array<LevelData<FArrayBox>, SpaceDim> alias;
    // for (int d = 0; d < SpaceDim; ++d) {
    //     FABAliasFlBxDataFactory factory(&a_vel, Interval(0,0), d);
    //     alias[d].define(m_grids, 1, a_vel.ghostVect(), factory);
    //     alias[d].exchangeBegin(m_exCopier[d]);
    // }

    // CFI BCs
    const bool cfRegionIsEmpty = D_TERM(m_cfRegion[0].isEmpty(),
                                        &&m_cfRegion[1].isEmpty(),
                                        &&m_cfRegion[2].isEmpty());
    if (!cfRegionIsEmpty) {
        if (a_homogCFIBCs) {
            CFInterp::homogInterpAtCFI(a_vel, m_dXi, m_amrCrseDXi, m_cfRegion);

        } else {
            StateType crseVel;
            CH_assert(a_crseVelPtr);
            this->create(crseVel, *a_crseVelPtr);
            this->assignLocal(crseVel, *a_crseVelPtr);
            this->sendToAdvectingVelocity(crseVel, m_dXi * m_cfInterp.getRefRatio());

            // This is inefficient in theory. Will revisit if necessary.
            this->sendToAdvectingVelocity(a_vel, m_dXi);
            m_cfInterp.interpGhostsAtCFI(a_vel, crseVel, a_homogCFIBCs);
            this->sendToCartesianVelocity(a_vel, m_dXi);
        }
    }

    // Phys BCs
    const auto key = m_physBdryIter.lock();
    for (m_physBdryIter.reset(key); m_physBdryIter.ok(); m_physBdryIter.next(key)) {
        const DataIndex&     di      = m_physBdryIter->di;
        const Box&           ccValid = m_physBdryIter->ccValidBox;
        const int            bdryDir = m_physBdryIter->dir;
        const Side::LoHiSide side    = m_physBdryIter->side;

        for (int velComp = 0; velComp < SpaceDim; ++velComp) {
            FArrayBox& stateFAB = a_vel[di][velComp];

            // Get important regions.
            Box stateBdry, stateGhosts;
            BCTools::getBoxesForApplyBC(
                stateBdry, stateGhosts, stateFAB, ccValid, bdryDir, side);

            // Get physical coordinates at boundary.
            FArrayBox xFAB(stateBdry, SpaceDim);
            m_geoSrc.fill_physCoor(xFAB, m_dXi);

            // Get grid spacing at boundary, dx = dx/dXi * dXi.
            FArrayBox dxFAB(stateBdry, 1);
            m_geoSrc.fill_dxdXi(dxFAB, 0, bdryDir, m_dXi, m_dXi[bdryDir]);

            // Set BCs, fill ghosts.
            BCTools::applyBC(stateFAB,
                             xFAB,
                             dxFAB,
                             ccValid,
                             di,
                             bdryDir,
                             side,
                             a_homogPhysBCs,
                             *m_bcFuncPtr,
                             a_time);
        }
    }
    m_physBdryIter.unlock(key);

    // Begin valid exchange (See StaggeredCopier::defineValidExchange for details.)
    std::array<LevelData<FArrayBox>, SpaceDim> alias;
    for (int d = 0; d < SpaceDim; ++d) {
        FABAliasFlBxDataFactory factory(&a_vel, Interval(0,0), d);
        alias[d].define(m_grids, 1, a_vel.ghostVect(), factory);
        alias[d].exchangeBegin(m_exCopier[d]);
    }
    // End valid exchange and exchange corner ghosts filled by BC-setting
    // functions(See StaggeredCopier::defineInvalidCornerExchange1 for details.)
    for (int d = 0; d < SpaceDim; ++d) {
        alias[d].exchangeEnd();
        // alias[d].exchangeBegin(m_exCornerCopier1[d]);
    }
    // for (int d = 0; d < SpaceDim; ++d) {
    //     alias[d].exchangeEnd();
    //     // alias[d].exchangeBegin(m_exCornerCopier2[d]);
    // }
    // for (int d = 0; d < SpaceDim; ++d) {
    //     alias[d].exchangeEnd();
    // }

    nanCheck(a_vel);
    debugCheckValidFaceOverlap(a_vel);
}


// -----------------------------------------------------------------------------
void
ViscousOp::applyOp(StateType&       a_lhs,
                   StateType&       a_vel,
                   const StateType* a_crseVelPtr,
                   const Real       a_time,
                   const bool       a_homogPhysBCs,
                   const bool       a_homogCFIBCs) const
{
    this->applyBCs(a_vel, a_crseVelPtr, a_time, a_homogPhysBCs, a_homogCFIBCs);

    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        for (int velComp = 0; velComp < SpaceDim; ++velComp) {
            const int d0 = velComp;
            const int d1 = (velComp + 1) % SpaceDim;
            const int d2 = (velComp + 2) % SpaceDim;

            FArrayBox&       lhsFAB  = a_lhs[dit][d0];
            const FArrayBox& velFAB  = a_vel[dit][d0];
            const FArrayBox& JFAB    = (*m_Jptr)[dit][d0];
            const FArrayBox& nuJgup0 = (*m_nuJgupPtr)[d0][d0][dit];
            const FArrayBox& nuJgup1 = (*m_nuJgupPtr)[d0][d1][dit];
            const FArrayBox& nuJgup2 = (*m_nuJgupPtr)[d0][d2][dit]; // Not used in 2D
            const FArrayBox& invDiagsFAB = (*m_invDiagsPtr)[dit][d0];

            debugInit(lhsFAB); // Should do nothing.

            const Box fcValid = m_grids[dit].surroundingNodes(d0);

            if (!m_bdryBox[d0][0].isEmpty()) {
                const Box destBox = fcValid & m_bdryBox[d0][0];
                if (!destBox.isEmpty()) {
                    lhsFAB.copy(velFAB, destBox);
                    lhsFAB.mult(JFAB, destBox, 0, 0, 1);
                    lhsFAB.mult(m_alpha);
                }
            }

            if (!m_bdryBox[d0][1].isEmpty()) {
                const Box destBox = fcValid & m_bdryBox[d0][1];
                if (!destBox.isEmpty()) {
                    lhsFAB.copy(velFAB, destBox);
                    lhsFAB.mult(JFAB, destBox, 0, 0, 1);
                    lhsFAB.mult(m_alpha);
                }
            }

            const Box destBox = fcValid & m_interiorBox[d0];
            CH_assert(!destBox.isEmpty());
            FORT_VISCOUSOP_APPLYOP(
                CHF_FRA1(lhsFAB, 0),
                CHF_CONST_FRA1(velFAB, 0),
                CHF_CONST_INT(d0),
                CHF_CONST_FRA1(nuJgup0, 0),
                CHF_CONST_FRA1(nuJgup1, 0),
                CHF_CONST_FRA1(nuJgup2, 0),
                CHF_CONST_REALVECT(m_dXi),
                CHF_BOX(destBox),
                CHF_CONST_FRA1(invDiagsFAB, 0),
                CHF_CONST_REAL(m_beta));
        } // velComp
    } // end loop over grids (dit)

    nanCheck(a_lhs);
    LayoutTools::averageOverlappingValidFaces(a_lhs);
}


// -----------------------------------------------------------------------------
void
ViscousOp::preCond(StateType&       a_vel,
                   const StateType& a_rhs,
                   const Real       a_time,
                   const int        a_relaxIters) const
{
    this->assignLocal(a_vel, a_rhs);
    {
        for (DataIterator dit(m_grids); dit.ok(); ++dit) {
            for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                const Box fcValid = m_grids[dit].surroundingNodes(velComp);
                a_vel[dit][velComp].mult(
                    (*m_invDiagsPtr)[dit][velComp], fcValid, 0, 0, 1);
            }
        }
    }
    this->relax(a_vel, a_rhs, a_time, a_relaxIters);
}


// -----------------------------------------------------------------------------
void
ViscousOp::relax(StateType&       a_cor,
                 const StateType& a_res,
                 const Real       a_time,
                 const int        a_relaxIters) const
{
    // this->jacobi_relax(a_cor, a_res, a_time, a_relaxIters);
    this->gsrb_relax(a_cor, a_res, a_time, a_relaxIters);
}


// =========================== MGOperator methods ==============================

// -----------------------------------------------------------------------------
MGOperator<LevelData<FluxBox>>*
ViscousOp::newMGOperator(const IntVect& a_refRatio) const
{
    ViscousOp* newOpPtr = nullptr;

    if (a_refRatio == IntVect::Unit) {
        newOpPtr = new ViscousOp(*this);
    } else {
        DisjointBoxLayout crseGrids;
        ::coarsen(crseGrids, m_grids, a_refRatio);

        newOpPtr = new ViscousOp(*this, crseGrids, a_refRatio);
    }

    return newOpPtr;
}


// -----------------------------------------------------------------------------
void
ViscousOp::MGRestrict(StateType&       a_crseRes,
                      const StateType& a_fineRes,
                      const Real       a_time,
                      const IntVect&   a_refRatio,
                      const MGOpType&  a_crseOp) const
{
    nanCheck(a_fineRes);

    constexpr bool           useFullWeighting = true;
    const DisjointBoxLayout& crseGrids        = a_crseRes.getBoxes();
    DataIterator             dit              = a_crseRes.dataIterator();
    const ViscousOp& crseViscousOp = dynamic_cast<const ViscousOp&>(a_crseOp);

    if (useFullWeighting) {
        for (dit.reset(); dit.ok(); ++dit) {
            for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                const Box fcCrseValid = crseGrids[dit].surroundingNodes(velComp);

                Box refBox(IntVect::Zero, a_refRatio - IntVect::Unit);
                refBox.setBig(velComp, 0);

                const IntVect& civ = fcCrseValid.smallEnd();
                const IntVect  fiv = civ * a_refRatio;

                FORT_VISCOUSOP_MGRESTRICT_PARENTAVERAGE(
                    CHF_FRA1_SHIFT(a_crseRes[dit][velComp], 0, civ),
                    CHF_FRA1_SHIFT(a_fineRes[dit][velComp], 0, fiv),
                    CHF_BOX_SHIFT(fcCrseValid, civ),
                    CHF_CONST_INT(velComp),
                    CHF_CONST_INTVECT(a_refRatio),
                    CHF_BOX(refBox));
            }
        }

    } else {
        for (dit.reset(); dit.ok(); ++dit) {
            for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                const Box fcCrseValid = crseGrids[dit].surroundingNodes(velComp)
                                      & crseViscousOp.m_interiorBox[velComp];

                Box refBox(IntVect::Zero, a_refRatio - IntVect::Unit);
                refBox.setBig(velComp, 0);

                const IntVect& civ = fcCrseValid.smallEnd();
                const IntVect  fiv = civ * a_refRatio;

                FORT_VISCOUSOP_MGRESTRICT_FULLWEIGHTING(
                    CHF_FRA1_SHIFT(a_crseRes[dit][velComp], 0, civ),
                    CHF_FRA1_SHIFT(a_fineRes[dit][velComp], 0, fiv),
                    CHF_BOX_SHIFT(fcCrseValid, civ),
                    CHF_CONST_INT(velComp),
                    CHF_CONST_INTVECT(a_refRatio),
                    CHF_BOX(refBox));
            }
        }
    }

    crseViscousOp.applyBCs(a_crseRes, nullptr, a_time, true, true);
    nanCheck(a_crseRes);
}


// -----------------------------------------------------------------------------
void
ViscousOp::MGProlong(StateType&      a_finePhi,
                     StateType&      a_crseCor,
                     const Real      /* a_time */,
                     const IntVect&  a_refRatio,
                     const MGOpType& /* a_crseOp */,
                     const int       a_interpOrder) const
{
    const DisjointBoxLayout& crseGrids = a_crseCor.getBoxes();
    DataIterator             dit       = a_finePhi.dataIterator();

    CH_assert(a_refRatio == a_finePhi.getBoxes().physDomain().size() /
                                crseGrids.physDomain().size());

    // This is a local operation.
    CH_assert(crseGrids.compatible(a_finePhi.getBoxes()));

    nanCheck(a_finePhi);
    debugCheckValidFaceOverlap(a_finePhi);

    nanCheck(a_crseCor);
    debugCheckValidFaceOverlap(a_crseCor);

    // Constant parent interp
    for (dit.reset(); dit.ok(); ++dit) {
        for (int velComp = 0; velComp < SpaceDim; ++velComp) {
            FArrayBox&       fineFAB   = a_finePhi[dit][velComp];
            const FArrayBox& crseFAB   = a_crseCor[dit][velComp];
            const Box        crseValid = crseGrids[dit].surroundingNodes(velComp);

            const IntVect& civ = crseValid.smallEnd();
            const IntVect  fiv = civ * a_refRatio;
            CH_assert(m_grids[dit].surroundingNodes(velComp).smallEnd() == fiv);

            // Interpolation to parent faces.
            Box refBox(IntVect::Zero, a_refRatio - IntVect::Unit, BASISV(velComp));
            refBox.setBig(velComp, 0);

            FORT_VISCOUSOP_MGPROLONG_PARENTINJECTION (
                CHF_FRA1_SHIFT(fineFAB, 0, fiv),
                CHF_CONST_FRA1_SHIFT(crseFAB, 0, civ),
                CHF_BOX_SHIFT(crseValid, civ),
                CHF_CONST_INTVECT(a_refRatio),
                CHF_BOX(refBox),
                CHF_CONST_INT(velComp));
        }
    }

    // Linear parent upgrade
    if (a_interpOrder >= 1) {
        for (dit.reset(); dit.ok(); ++dit) {
            for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                FArrayBox&       fineFAB   = a_finePhi[dit][velComp];
                const FArrayBox& crseFAB   = a_crseCor[dit][velComp];
                const Box        crseValid = crseGrids[dit].surroundingNodes(velComp);

                const IntVect& civ = crseValid.smallEnd();
                const IntVect  fiv = civ * a_refRatio;
                CH_assert(m_grids[dit].surroundingNodes(velComp).smallEnd() == fiv);

                // Interpolation to parent faces.
                Box refBox(IntVect::Zero, a_refRatio - IntVect::Unit, BASISV(velComp));
                refBox.setBig(velComp, 0);

                FORT_VISCOUSOP_MGPROLONG_PARENTLINEARUPGRADE (
                    CHF_FRA1_SHIFT(fineFAB, 0, fiv),
                    CHF_CONST_FRA1_SHIFT(crseFAB, 0, civ),
                    CHF_BOX_SHIFT(crseValid, civ),
                    CHF_CONST_INTVECT(a_refRatio),
                    CHF_BOX(refBox),
                    CHF_CONST_INT(velComp));
            }
        }
    }

    // Quadratic parent upgrade.
    // Is there room for the stencil?
    bool noRoom = false;
    if (a_interpOrder >= 2) {
        LayoutIterator lit = crseGrids.layoutIterator();
        for (lit.reset(); lit.ok(); ++lit) {
            const IntVect& size = crseGrids[lit].size();
            if ( !( size >= IntVect(D_DECL(4,4,4)) ) ) {
                noRoom = true;
                break;
            }
        }
    }
    if (a_interpOrder >= 2 && !noRoom) {
        // Fill *ALL* ghosts. Some are already filled, the corners are not.
        for (dit.reset(); dit.ok(); ++dit) {
            for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                const Box fcValid = crseGrids[dit].surroundingNodes(velComp);
                BCTools::extrapCorners(a_crseCor[dit][velComp], fcValid, 2);
            }
        }
        // const auto& crseViscousOp = static_cast<const ViscousOp&>(a_crseOp);
        // CH_verify(crseViscousOp.m_HOProlongCornerCopier.isDefined());
        // a_crseCor.exchange(crseViscousOp.m_HOProlongCornerCopier);
        a_crseCor.exchange();

        for (dit.reset(); dit.ok(); ++dit) {
            for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                FArrayBox&       fineFAB   = a_finePhi[dit][velComp];
                const FArrayBox& crseFAB   = a_crseCor[dit][velComp];
                const Box        crseValid = crseGrids[dit].surroundingNodes(velComp);

                const IntVect& civ = crseValid.smallEnd();
                const IntVect  fiv = civ * a_refRatio;
                CH_assert(m_grids[dit].surroundingNodes(velComp).smallEnd() == fiv);

                Box refBox(IntVect::Zero, a_refRatio - IntVect::Unit, BASISV(velComp));
                refBox.setBig(velComp, 0);

                FORT_VISCOUSOP_MGPROLONG_PARENTQUADUPGRADE1 (
                    CHF_FRA_SHIFT(fineFAB, fiv),
                    CHF_CONST_FRA_SHIFT(crseFAB, civ),
                    CHF_BOX_SHIFT(crseValid, civ),
                    CHF_CONST_INTVECT(a_refRatio),
                    CHF_BOX(refBox),
                    CHF_CONST_INT(velComp));

                FORT_VISCOUSOP_MGPROLONG_PARENTQUADUPGRADE2 (
                    CHF_FRA_SHIFT(fineFAB, fiv),
                    CHF_CONST_FRA_SHIFT(crseFAB, civ),
                    CHF_BOX_SHIFT(crseValid, civ),
                    CHF_CONST_INTVECT(a_refRatio),
                    CHF_BOX(refBox),
                    CHF_CONST_INT(velComp));
            }
        }
    }


    // Linear interior faces
    for (dit.reset(); dit.ok(); ++dit) {
        for (int velComp = 0; velComp < SpaceDim; ++velComp) {
            FArrayBox&       fineFAB   = a_finePhi[dit][velComp];
            // const FArrayBox& crseFAB   = a_crseCor[dit][velComp];
            const Box        crseValid = crseGrids[dit].surroundingNodes(velComp);

            const IntVect& civ = crseValid.smallEnd();
            const IntVect  fiv = civ * a_refRatio;
            CH_assert(m_grids[dit].surroundingNodes(velComp).smallEnd() == fiv);

            // Offsets from the lower parent face to the interior faces.
            Box interiorRefBox(IntVect::Zero, a_refRatio - IntVect::Unit);
            interiorRefBox.surroundingNodes(velComp);
            interiorRefBox.grow(velComp, -1);

            Box ccCrseInterpBox = crseValid;
            ccCrseInterpBox.enclosedCells(velComp);

            FORT_VISCOUSOP_MGPROLONG_LINEARINTERIORFACE (
                CHF_FRA1_SHIFT(fineFAB, 0, fiv),
                CHF_BOX_SHIFT(ccCrseInterpBox, civ),
                CHF_CONST_INTVECT(a_refRatio),
                CHF_CONST_INT(velComp),
                CHF_BOX(interiorRefBox));
        }
    }

    // Quadratic interior upgrade
    if (a_interpOrder >= 2 && !noRoom) {
        UNDEFINED_FUNCTION();

        // for (dit.reset(); dit.ok(); ++dit) {
        //     for (int velComp = 0; velComp < SpaceDim; ++velComp) {
        //         if (a_refRatio[velComp] == 1) continue;

        //         FArrayBox&       fineFAB   = a_finePhi[dit][velComp];
        //         const FArrayBox& crseFAB   = a_crseCor[dit][velComp];
        //         const Box        crseValid = crseGrids[dit].surroundingNodes(velComp);

        //         const IntVect& civ = crseValid.smallEnd();
        //         const IntVect  fiv = civ * a_refRatio;
        //         CH_assert(m_grids[dit].surroundingNodes(velComp).smallEnd() == fiv);

        //         // Offsets from the lower parent face to the interior faces.
        //         Box interiorRefBox(IntVect::Zero, a_refRatio - IntVect::Unit);
        //         interiorRefBox.surroundingNodes(velComp);
        //         interiorRefBox.grow(velComp, -1);

        //         Box ccCrseInterpBox = crseValid;
        //         ccCrseInterpBox.enclosedCells(velComp);

        //         FORT_VISCOUSOP_MGPROLONG_LINEARINTERIORFACE (
        //             CHF_FRA1_SHIFT(fineFAB, 0, fiv),
        //             CHF_BOX_SHIFT(ccCrseInterpBox, civ),
        //             CHF_CONST_INTVECT(a_refRatio),
        //             CHF_CONST_INT(velComp),
        //             CHF_BOX(interiorRefBox));
        //     }
        // }
    }

        // const ViscousOp& crseViscousOp =
        //     dynamic_cast<const ViscousOp&>(a_crseOp);

        // for (dit.reset(); dit.ok(); ++dit) {
        //     for (int velComp = 0; velComp < SpaceDim; ++velComp) {
        //         FArrayBox&       fineFAB   = a_finePhi[dit][velComp];
        //         const FArrayBox& crseFAB   = a_crseCor[dit][velComp];
        //         const Box fcCrseInterpBox  = crseGrids[dit].surroundingNodes(velComp);

        //         CFInterp::localRefineFace(
        //             fineFAB, crseFAB, fcCrseInterpBox, velComp, a_refRatio);
        //     }
        // }
    // }

    nanCheck(a_finePhi);
}



// ============================ Private methods ================================

// -----------------------------------------------------------------------------
ViscousOp::ViscousOp(const ViscousOp& a_srcOp)
: m_alpha(a_srcOp.m_alpha)
, m_beta(a_srcOp.m_beta)
, m_domain(a_srcOp.m_domain)
, m_grids(a_srcOp.m_grids)
, m_dXi(a_srcOp.m_dXi)
, m_amrCrseDXi(a_srcOp.m_amrCrseDXi)
, m_geoSrc(a_srcOp.m_geoSrc)
, m_Jptr(a_srcOp.m_Jptr)
, m_nuJgupPtr(a_srcOp.m_nuJgupPtr)
, m_invDiagsPtr(a_srcOp.m_invDiagsPtr)
, m_bcFuncPtr(a_srcOp.m_bcFuncPtr)
, m_physBdryIter(a_srcOp.m_physBdryIter)
, m_cfRegion(a_srcOp.m_cfRegion)
// m_cfInterp
, m_exCopier(a_srcOp.m_exCopier)
, m_exCornerCopier1(a_srcOp.m_exCornerCopier1)
, m_exCornerCopier2(a_srcOp.m_exCornerCopier2)
, m_interiorBox(a_srcOp.m_interiorBox)
, m_bdryBox(a_srcOp.m_bdryBox)
, m_resPtr(a_srcOp.m_resPtr)
{
    // CFInterp
    if (a_srcOp.m_cfInterp.isDefined()) {
        m_cfInterp.define(m_grids, m_dXi, a_srcOp.m_cfInterp.getCrseBoxes());
    }
}


// -----------------------------------------------------------------------------
ViscousOp::ViscousOp(const ViscousOp&         a_srcOp,
                     const DisjointBoxLayout& a_crseGrids,
                     const IntVect&           a_refRatio)
: m_alpha(a_srcOp.m_alpha)
, m_beta(a_srcOp.m_beta)
, m_domain(a_crseGrids.physDomain())
, m_grids(a_crseGrids)
, m_dXi(a_srcOp.m_dXi * a_refRatio)
, m_amrCrseDXi(a_srcOp.m_amrCrseDXi)
, m_geoSrc(a_srcOp.m_geoSrc)
, m_Jptr(new LevelData<FluxBox>(a_crseGrids, 1, IntVect::Unit))
, m_nuJgupPtr(new StaggeredFluxLD(a_crseGrids))
, m_invDiagsPtr(new LevelData<FluxBox>(a_crseGrids, 1))
, m_bcFuncPtr(a_srcOp.m_bcFuncPtr)
, m_physBdryIter(a_crseGrids)
// , m_cfRegion(a_srcOp.m_cfRegion) // We will coarsen in the body!
// m_cfInterp (Leave undefined)
// exchange copiers
// m_interiorBox
// m_bdryBox
, m_resPtr(new LevelData<FluxBox>(a_crseGrids, 1))
{
    // m_cfRegion
    for (int velComp = 0; velComp < SpaceDim; ++velComp) {
        m_cfRegion[velComp].define(m_grids, m_domain, BASISV(velComp));
    }

    // Exchange copiers
    for (int d = 0; d < SpaceDim; ++d) {
        constexpr IntVect ghostVect(D_DECL(1, 1, 1));
        constexpr bool doValidCorners = true;

        m_exCopier[d].defineValidExchange(m_grids, ghostVect, d, doValidCorners);
        m_exCornerCopier1[d].defineInvalidCornerExchange1(m_grids, ghostVect, d);
        m_exCornerCopier2[d].defineInvalidCornerExchange2(m_grids, ghostVect, d);
    }

    // Where L[u] needs to be computed.
    for (int velComp = 0; velComp < SpaceDim; ++velComp) {
        m_interiorBox[velComp] = m_domain.domainBox();
        m_interiorBox[velComp].surroundingNodes(velComp);
        if (!m_domain.isPeriodic(velComp)) {
            m_bdryBox[velComp][0] = bdryLo(m_interiorBox[velComp], velComp);
            m_bdryBox[velComp][1] = bdryHi(m_interiorBox[velComp], velComp);
            m_interiorBox[velComp].grow(velComp, -1);
        }
    }

    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        for (int velComp = 0; velComp < SpaceDim; ++velComp) {
            CFInterp::localCoarsenFace((*m_Jptr)[dit][velComp],
                                       (*a_srcOp.m_Jptr)[dit][velComp],
                                       m_grids[dit].surroundingNodes(velComp),
                                       velComp,
                                       a_refRatio,
                                       nullptr);
            CFInterp::localCoarsen(
                (*m_nuJgupPtr)[velComp][velComp][dit],
                (*a_srcOp.m_nuJgupPtr)[velComp][velComp][dit],
                m_grids[dit],
                a_refRatio,
                false,
                nullptr);

            for (int derivDir = 0; derivDir < SpaceDim; ++derivDir) {
                if (velComp == derivDir) continue;

                const Box crseBox =
                    m_grids[dit].surroundingNodes(velComp).surroundingNodes(
                        derivDir);
                if (SpaceDim == 2) {
                    CFInterp::localCoarsenNode(
                        (*m_nuJgupPtr)[velComp][derivDir][dit],
                        (*a_srcOp.m_nuJgupPtr)[velComp][derivDir][dit],
                        crseBox,
                        a_refRatio);
                } else {
                    CFInterp::localCoarsenEdge(
                        (*m_nuJgupPtr)[velComp][derivDir][dit],
                        (*a_srcOp.m_nuJgupPtr)[velComp][derivDir][dit],
                        crseBox,
                        3 - velComp - derivDir,
                        a_refRatio,
                        nullptr);
                }
            }
        }
    }  // dit

    BCTools::extrapAllGhosts(*m_Jptr, 2);
    m_Jptr->exchange();

    for (int velComp = 0; velComp < SpaceDim; ++velComp) {
        for (int derivDir = 0; derivDir < SpaceDim; ++derivDir) {
            BCTools::extrapAllGhosts((*m_nuJgupPtr)[velComp][derivDir], 2);
        }
    }
    m_nuJgupPtr->exchange();

    // Compute diagonals, etc.
    this->cacheMatrixElements();
}


// -----------------------------------------------------------------------------
void
ViscousOp::cacheMatrixElements()
{
    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        for (int velComp = 0; velComp < SpaceDim; ++velComp) {
            FArrayBox& invDiagsFAB = (*m_invDiagsPtr)[dit][velComp];
            const Box  destBox     = invDiagsFAB.box();

            const int d0 = velComp;
            const int d1 = (velComp + 1) % SpaceDim;
            const int d2 = (velComp + 2) % SpaceDim;

            const FArrayBox& nuJgup0 = (*m_nuJgupPtr)[d0][d0][dit];
            const FArrayBox& nuJgup1 = (*m_nuJgupPtr)[d0][d1][dit];
            const FArrayBox& nuJgup2 = (*m_nuJgupPtr)[d0][d2][dit]; // Not used in 2D
            const FArrayBox& JFAB    = (*m_Jptr)[dit][d0];

            debugInit(invDiagsFAB);
            FORT_VISCOUSOP_COMPUTEDIAGS(
                CHF_FRA1(invDiagsFAB, 0),
                CHF_CONST_INT(velComp),
                CHF_CONST_FRA1(nuJgup0, 0),
                CHF_CONST_FRA1(nuJgup1, 0),
                CHF_CONST_FRA1(nuJgup2, 0),
                CHF_CONST_FRA1(JFAB, 0),
                CHF_CONST_REALVECT(m_dXi),
                CHF_BOX(destBox),
                CHF_CONST_REAL(m_alpha),
                CHF_CONST_REAL(m_beta));

            invDiagsFAB.invert(1.0);
        }
    }

    nanCheck(*m_invDiagsPtr);
    LayoutTools::averageOverlappingValidFaces(*m_invDiagsPtr);
}


// -----------------------------------------------------------------------------
void
ViscousOp::sendToAdvectingVelocity(StateType&      a_vel,
                                   const RealVect& a_dXi) const
{
    const DisjointBoxLayout& grids = a_vel.getBoxes();

    for (DataIterator dit(grids); dit.ok(); ++dit) {
        for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
            FArrayBox& velFAB = a_vel[dit][fcDir];
            FArrayBox  dxdXiFAB(velFAB.box(), 1);

            for (int offset = 1; offset < SpaceDim; ++offset) {
                const int mu = (fcDir + offset) % SpaceDim;
                m_geoSrc.fill_dxdXi(dxdXiFAB, 0, mu, a_dXi);
                velFAB.mult(dxdXiFAB, 0, 0, 1);
            }
        } // fcDir
    } // dit
}


// -----------------------------------------------------------------------------
void
ViscousOp::sendToCartesianVelocity(StateType&      a_vel,
                                   const RealVect& a_dXi) const
{
    const DisjointBoxLayout& grids = a_vel.getBoxes();

    for (DataIterator dit(grids); dit.ok(); ++dit) {
        for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
            FArrayBox& velFAB = a_vel[dit][fcDir];
            FArrayBox  dxdXiFAB(velFAB.box(), 1);

            for (int offset = 1; offset < SpaceDim; ++offset) {
                const int mu = (fcDir + offset) % SpaceDim;
                m_geoSrc.fill_dxdXi(dxdXiFAB, 0, mu, a_dXi);
                velFAB.divide(dxdXiFAB, 0, 0, 1);
            }
        } // fcDir
    } // dit
}


// -----------------------------------------------------------------------------
void
ViscousOp::jacobi_relax(StateType&       a_vel,
                        const StateType& a_rhs,
                        const Real       a_time,
                        const int        a_iters) const
{
    for (int iter = 0; iter < a_iters; ++iter) {
        this->residual(*m_resPtr, a_vel, nullptr, a_rhs, a_time, true, true);

        for (DataIterator dit(m_grids); dit.ok(); ++dit) {
            for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                const Box fcValid = m_grids[dit].surroundingNodes(velComp)
                                  & m_interiorBox[velComp];

                FORT_VISCOUSOP_JACOBI(
                    CHF_FRA1(a_vel[dit][velComp], 0),
                    CHF_CONST_FRA1((*m_resPtr)[dit][velComp], 0),
                    CHF_CONST_FRA1((*m_invDiagsPtr)[dit][velComp], 0),
                    CHF_BOX(fcValid));
            }
        }
    }
}


// -----------------------------------------------------------------------------
void
ViscousOp::gsrb_relax(StateType&       a_vel,
                      const StateType& a_rhs,
                      const Real       a_time,
                      const int        a_iters) const
{
    for (int iter = 0; iter < a_iters; ++iter) {
        for (int whichPass = 0; whichPass < 2; ++whichPass) {
            // Bottleneck!
            if (whichPass == 0) {
                this->applyBCs(a_vel, nullptr, a_time, true, true);
            } else {
                a_vel.exchange();
            }

            for (DataIterator dit(m_grids); dit.ok(); ++dit) {
                for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                    const int d0 = velComp;
                    const int d1 = (velComp + 1) % SpaceDim;
                    const int d2 = (velComp + 2) % SpaceDim;

                    FArrayBox&       velFAB  = a_vel[dit][d0];
                    const FArrayBox& rhsFAB  = a_rhs[dit][d0];
                    const FArrayBox& nuJgup0 = (*m_nuJgupPtr)[d0][d0][dit];
                    const FArrayBox& nuJgup1 = (*m_nuJgupPtr)[d0][d1][dit];
                    const FArrayBox& nuJgup2 = (*m_nuJgupPtr)[d0][d2][dit]; // Not used in 2D
                    const FArrayBox& invDiagsFAB = (*m_invDiagsPtr)[dit][d0];

                    const Box destBox = m_grids[dit].surroundingNodes(d0)
                                      & m_interiorBox[d0];

                    //checkForNAN(velFAB, velFAB.box());
                    checkForNAN(rhsFAB, destBox);
                    checkForNAN(nuJgup0, nuJgup0.box());
                    checkForNAN(nuJgup1, nuJgup1.box());
                    checkForNAN(nuJgup2, nuJgup2.box());
                    checkForNAN(invDiagsFAB, invDiagsFAB.box());

                    FORT_VISCOUSOP_GSRB(
                        CHF_FRA1(velFAB, 0),
                        CHF_CONST_FRA1(rhsFAB, 0),
                        CHF_CONST_INT(d0),
                        CHF_CONST_FRA1(nuJgup0, 0),
                        CHF_CONST_FRA1(nuJgup1, 0),
                        CHF_CONST_FRA1(nuJgup2, 0),
                        CHF_CONST_REALVECT(m_dXi),
                        CHF_BOX(destBox),
                        CHF_CONST_FRA1(invDiagsFAB, 0),
                        CHF_CONST_REAL(m_beta),
                        CHF_CONST_INT(whichPass));
                } // velComp
            }  // dit
        } // whichPass (red or black)
    } // iter
}



}; // end namespace Elliptic
