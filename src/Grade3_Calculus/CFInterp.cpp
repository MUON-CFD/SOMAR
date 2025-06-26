#include "CFInterp.H"
#include "CFInterpF_F.H"
#include "BoxIterator.H"
#include "LayoutTools.H"
#include "SOMAR_Constants.H"
#include "AnisotropicRefinementTools.H"
#include "Convert.H"
#include "BCTools.H"
#include "Debug.H"

#include "AnisotropicLinearCFInterp.H" // TEMPORARY!!!
#include "MappedQuadCFInterp.H"        // TEMPORARY!!!


// ======================= Construction / destruction ==========================

// -----------------------------------------------------------------------------
CFInterp::CFInterp()
: m_isDefined(false)
{
}


// -----------------------------------------------------------------------------
CFInterp::CFInterp(const LevelGeometry&     a_levGeo,
                   const DisjointBoxLayout& a_userCrseGrids)
: m_isDefined(false)
{
    this->define(a_levGeo.getBoxes(), a_levGeo.getDXi(), a_userCrseGrids);
}


// -----------------------------------------------------------------------------
CFInterp::~CFInterp()
{
    this->undefine();
}


// -----------------------------------------------------------------------------
void
CFInterp::define(const LevelGeometry&     a_levGeo,
                 const DisjointBoxLayout& a_userCrseGrids)
{
    this->define(a_levGeo.getBoxes(), a_levGeo.getDXi(), a_userCrseGrids);
}


// -----------------------------------------------------------------------------
void
CFInterp::define(const DisjointBoxLayout& a_grids,
                 const RealVect&          a_dXi,
                 const DisjointBoxLayout& a_userCrseGrids)
{
    // Ensure we are starting clean.
    if (m_isDefined) this->undefine();

    // User's fine-level grids.
    m_domain    = a_grids.physDomain();
    m_grids     = a_grids;
    m_refRatio  = IntVect::Unit;
    m_fineDXi   = a_dXi;
    m_cfiIter.define(m_grids, CFRegion(m_grids, m_domain));

    m_userCrseGrids = a_userCrseGrids;
    {
        // Fix the ref ratio.
        const ProblemDomain& crseDomain = m_userCrseGrids.physDomain();
        m_refRatio = calculateRefinementRatio(crseDomain.domainBox(),
                                              m_domain.domainBox());

        // Coarsened fine-level stuff.
        // These grids are compatible with m_grids.
        ::coarsen(m_crseGrids, m_grids, m_refRatio);
        CH_assert(m_crseGrids.physDomain() == crseDomain);
        m_crseGhosts = IntVect(D_DECL(1, 1, 1));
        m_crseExCopier.exchangeDefine(m_crseGrids, m_crseGhosts);
        m_crseExCornerCopier.define(m_crseGrids,
                                    m_crseGrids,
                                    m_crseGrids.physDomain(),
                                    m_crseGhosts,
                                    true);

        m_crseDXi = m_fineDXi * RealVect(m_refRatio);

        // Copiers for m_userCrseGrids <-> m_crseGrids operations.
        m_crseToUserCopier.define(m_crseGrids, m_userCrseGrids);

#ifdef ALLOW_DIVFREEINTERP
        // Matrix inversion.
        this->defineRowIdx(m_refRatio);

        if (s_mode == 0) {
            this->defineColIdx_Original(m_refRatio);
            this->invertMatrix_Original(m_refRatio);
        } else if (s_mode == 1) {
            this->defineColIdx_StdLS(m_refRatio);
            this->invertMatrix_StdLS(m_refRatio);
        } else {
            MAYDAYERROR("CFInterp does not recognize s_mode = " << s_mode);
        }

        // Clear unnecessary holders.
        if (m_refRatio[0] == 1) m_rowIdx[0].clear();
        if (m_refRatio[1] == 1) m_rowIdx[1].clear();
        if (m_refRatio[0] == 1 || m_refRatio[1] == 1) m_colIdx[1].clear();
#endif // ALLOW_DIVFREEINTERP

    } // end coarser level stuff.

    m_isDefined = true;
}


#ifdef ALLOW_DIVFREEINTERP
// -----------------------------------------------------------------------------
void
CFInterp::defineRowIdx(const IntVect& a_refRatio)
{
    const int rx = a_refRatio[0];
    const int ry = a_refRatio[1];

    if (rx > 1 || ry > 1) {
        // const int rank = 2 * (rx * ry) - (rx + ry);
        const Box ccRefBox(IntVect::Zero, IntVect(D_DECL(rx - 1, ry - 1, 0)));
        Box       region;
        int       idx = 0;

        // Set the u indices.
        if (rx > 1) {
            region = ccRefBox;
            region.surroundingNodes(0);
            region.grow(0, -1);
            CH_assert(region.numPts() == (rx - 1) * ry);

            m_rowIdx[0].define(region, 1);

            for (BoxIterator bit(region); bit.ok(); ++bit) {
                m_rowIdx[0](bit()) = idx++;
            }
        } else {
            // m_rowIdx will not be used, but I need to define it
            // in order to pass it to Fortran. We will undefine it
            // after the Fortran call.
            region.define(IntVect::Zero, IntVect::Zero);
            m_rowIdx[0].define(region, 1);
        }

        // Set the v indices.
        if (ry > 1) {
            region = ccRefBox;
            region.surroundingNodes(1);
            region.grow(1, -1);
            CH_assert(region.numPts() == rx * (ry - 1));

            m_rowIdx[1].define(region, 1);

            for (BoxIterator bit(region); bit.ok(); ++bit) {
                m_rowIdx[1](bit()) = idx++;
            }
        } else {
            // See m_rowIdx[0] comment.
            region.define(IntVect::Zero, IntVect::Zero);
            m_rowIdx[1].define(region, 1);
        }
        CH_assert(idx == 2 * (rx * ry) - (rx + ry));
    }
}


// -----------------------------------------------------------------------------
void
CFInterp::defineColIdx_Original(const IntVect& a_refRatio)
{
    const int rx = a_refRatio[0];
    const int ry = a_refRatio[1];

    if (rx > 1 || ry > 1) {
        // const int rank = 2 * (rx * ry) - (rx + ry);
        const Box ccRefBox(IntVect::Zero, IntVect(D_DECL(rx - 1, ry - 1, 0)));
        Box       region;
        int       idx = 0;

        // Set the div[u] indices
        {
            region = ccRefBox;
            CH_assert(region.numPts() == rx * ry);

            m_colIdx[0].define(region, 1);

            BoxIterator bit(region);
            m_colIdx[0](bit()) = -1; // Skip first one. It's lin. dep.
            ++bit;
            for (; bit.ok(); ++bit) {
                m_colIdx[0](bit()) = idx++;
            }
        }

        // Set the curl indices.
        if (rx > 1 && ry > 1) {
            region = ccRefBox;
            region.surroundingNodes();
            region.grow(0, -1);
            region.grow(1, -1);
            if (SpaceDim == 3) {
                region.setBig(2, 0);
            }
            CH_assert(region.numPts() == (rx - 1) * (ry - 1));

            m_colIdx[1].define(region, 1);

            for (BoxIterator bit(region); bit.ok(); ++bit) {
                m_colIdx[1](bit()) = idx++;
            }
        } else {
            // See m_rowIdx[0] comment.
            region.define(IntVect::Zero, IntVect::Zero);
            m_colIdx[1].define(region, 1);
        }
        CH_assert(idx == 2 * (rx * ry) - (rx + ry));
    }
}


// -----------------------------------------------------------------------------
void
CFInterp::invertMatrix_Original(const IntVect& a_refRatio)
{
    const int rx = a_refRatio[0];
    const int ry = a_refRatio[1];

    if (rx > 1 || ry > 1) {
        const int rank = 2 * (rx * ry) - (rx + ry);
        const Box ccRefBox(IntVect::Zero, IntVect(D_DECL(rx - 1, ry - 1, 0)));
        Box       region;

        // Define V matrix and workspace.
        {
            const Box MBox(IntVect::Zero,
                            IntVect(D_DECL(rank - 1, rank - 1, 0)));
            m_Minv.define(MBox, 1);

            // Invert M in-place.
            FORT_CFINTERP_MATRIXINVERT_ORIGINAL(
                CHF_FRA1(m_Minv, 0),
                CHF_CONST_FIA1(m_rowIdx[0], 0),
                CHF_CONST_FIA1(m_rowIdx[1], 0),
                CHF_CONST_FIA1(m_colIdx[0], 0),
                CHF_CONST_FIA1(m_colIdx[1], 0),
                CHF_CONST_REALVECT(m_fineDXi),
                CHF_CONST_INT(rx),
                CHF_CONST_INT(ry));
        }
    }
}


// -----------------------------------------------------------------------------
void
CFInterp::defineColIdx_StdLS(const IntVect& a_refRatio)
{
    const int rx = a_refRatio[0];
    const int ry = a_refRatio[1];

    if (rx > 1 || ry > 1) {
        // const int rank = 2 * (rx * ry) - (rx + ry);
        const Box ccRefBox(IntVect::Zero, IntVect(D_DECL(rx - 1, ry - 1, 0)));
        Box       region;
        int       idx = 0;

        // Set the div[u] indices
        {
            region = ccRefBox;
            CH_assert(region.numPts() == rx * ry);

            m_colIdx[0].define(region, 1);

            BoxIterator bit(region);
            for (; bit.ok(); ++bit) {
                m_colIdx[0](bit()) = idx++;
            }
        }

        // Set the curl indices.
        if (rx > 1 && ry > 1) {
            region = ccRefBox;
            region.surroundingNodes();
            region.grow(0, -1);
            region.grow(1, -1);
            if (SpaceDim == 3) {
                region.setBig(2, 0);
            }
            CH_assert(region.numPts() == (rx - 1) * (ry - 1));

            m_colIdx[1].define(region, 1);

            for (BoxIterator bit(region); bit.ok(); ++bit) {
                m_colIdx[1](bit()) = idx++;
            }
        } else {
            // See m_rowIdx[0] comment.
            region.define(IntVect::Zero, IntVect::Zero);
            m_colIdx[1].define(region, 1);
        }
        CH_assert(idx == 2 * (rx * ry) - (rx + ry) + 1);
    }
}


// -----------------------------------------------------------------------------
void
CFInterp::invertMatrix_StdLS(const IntVect& a_refRatio)
{
    const int rx = a_refRatio[0];
    const int ry = a_refRatio[1];

    if (rx > 1 || ry > 1) {
        const int rank = 2 * (rx * ry) - (rx + ry);
        const Box ccRefBox(IntVect::Zero, IntVect(D_DECL(rx - 1, ry - 1, 0)));
        Box       region;

        // Define V matrix and workspace.
        {
            Box MBox(IntVect::Zero, IntVect(D_DECL(rank - 1, rank - 1, 0)));
            m_Minv.define(MBox, 1);

            MBox.define(IntVect::Zero, IntVect(D_DECL(rank, rank - 1, 0)));
            BaseFab<Real> DC(MBox, 1);

            MBox.define(IntVect::Zero, IntVect(D_DECL(rank - 1, rank, 0)));
            m_DCt.define(MBox, 1);

            // Invert M in-place.
            FORT_CFINTERP_MATRIXINVERT_STDLS(
                CHF_FRA1(m_Minv, 0),
                CHF_FRA1(DC, 0),
                CHF_FRA1(m_DCt, 0),
                CHF_CONST_FIA1(m_rowIdx[0], 0),
                CHF_CONST_FIA1(m_rowIdx[1], 0),
                CHF_CONST_FIA1(m_colIdx[0], 0),
                CHF_CONST_FIA1(m_colIdx[1], 0),
                CHF_CONST_REALVECT(m_fineDXi),
                CHF_CONST_INT(rx),
                CHF_CONST_INT(ry));
        }
    }
}
#endif // ALLOW_DIVFREEINTERP


// -----------------------------------------------------------------------------
void
CFInterp::undefine()
{
    if (!m_isDefined) return;

    m_rowIdx[0].clear();
    m_rowIdx[1].clear();
    m_colIdx[0].clear();
    m_colIdx[1].clear();
    m_Minv.clear();

    m_crseToUserCopier.clear();

    m_crseGrids  = DisjointBoxLayout();
    m_crseGhosts = IntVect(D_DECL(-1, -1, -1));
    m_crseExCopier.clear();
    m_crseExCornerCopier.clear();

    m_userCrseGrids = DisjointBoxLayout();

    m_cfiIter.clear();
    m_domain   = ProblemDomain();
    m_grids    = DisjointBoxLayout();
    m_refRatio = IntVect::Zero;
    m_fineDXi  = RealVect::Zero;
    m_crseDXi  = RealVect::Zero;

    m_isDefined = false;
}


// ============================ CC interpolators ===============================

// -----------------------------------------------------------------------------
void
CFInterp::interpAtCFI(LevelData<FArrayBox>&       a_fine,
                      const LevelData<FArrayBox>& a_crse,
                      const bool                  a_isHomog) const
{
    // Is there anything to do?
    if (!this->hasCFI()) return;

    CH_assert(m_isDefined);
    CH_assert(a_fine.getBoxes() == m_grids);

    if (a_isHomog) {
        CFInterp::homogInterpAtCFI(a_fine, m_fineDXi, m_crseDXi, m_cfiIter);
    } else {
        // AnisotropicLinearCFInterp interpObj(a_fine.getBoxes(),
        //                                     a_crse.getBoxes(),
        //                                     a_fine.nComp(),
        //                                     a_crse.getBoxes().physDomain(),
        //                                     m_refRatio,
        //                                     a_fine.ghostVect());
        // interpObj.fillInterp(a_fine, a_crse, a_crse, 1.0, 0, 0, a_fine.nComp());

        MappedQuadCFInterp interpObj(a_fine.getBoxes(),
                                     &a_crse.getBoxes(),
                                     m_fineDXi,
                                     m_refRatio,
                                     a_fine.nComp(),
                                     a_fine.getBoxes().physDomain());
        interpObj.coarseFineInterp(a_fine, a_crse);
    }
}


// -----------------------------------------------------------------------------
void
CFInterp::homogInterpAtCFI(LevelData<FArrayBox>& a_fine) const
{
    CFInterp::homogInterpAtCFI(a_fine, m_fineDXi, m_crseDXi, m_cfiIter);
}


// -----------------------------------------------------------------------------
void
CFInterp::homogInterpAtCFI(LevelData<FArrayBox>& a_fine,
                           const RealVect&       a_fineDXi,
                           const RealVect&       a_crseDXi,
                           CFIIter&              a_cfiIter)
{
    // Is there anything to do?
    if (a_cfiIter.isEmpty()) return;

    for (a_cfiIter.reset(); a_cfiIter.ok(); ++a_cfiIter) {
        const DataIndex& di      = a_cfiIter->di;
        const int        isign   = a_cfiIter->isign;
        const int        bdryDir = a_cfiIter->dir;
        const CFIVS&     cfivs   = a_cfiIter->cfivs;
        const Box&       valid   = a_cfiIter->ccValidBox;

        FArrayBox& fineFAB  = a_fine[di];
        const int  fineSize = valid.size(bdryDir);
        const int  ncomp    = fineFAB.nComp();

        const Real dxf = a_fineDXi[bdryDir];
        const Real dxc = a_crseDXi[bdryDir];
        CH_assert(dxc > 0.0);

        if (cfivs.isPacked()) {
            const Box& ghostAdjCFI = cfivs.packedBox();

            if (fineSize > 1) {
                FORT_CFINTERP_HOMOGINTERPATCFI_QUAD (
                    CHF_FRA(fineFAB),
                    CHF_BOX(ghostAdjCFI),
                    CHF_CONST_REAL(dxf),
                    CHF_CONST_REAL(dxc),
                    CHF_CONST_INT(bdryDir),
                    CHF_CONST_INT(isign));
            } else {
                FORT_CFINTERP_HOMOGINTERPATCFI_LINEAR (
                    CHF_FRA(fineFAB),
                    CHF_BOX(ghostAdjCFI),
                    CHF_CONST_REAL(dxf),
                    CHF_CONST_REAL(dxc),
                    CHF_CONST_INT(bdryDir),
                    CHF_CONST_INT(isign));
            }

        } else {
            IVSIterator ivsit(cfivs.getIVS());

            if (fineSize > 1) {
                const Real c1 = 2.0 * (dxc - dxf) / (dxc + dxf);
                const Real c2 = -(dxc - dxf) / (dxc + 3.0 * dxf);
                Real       p1, p2;
                IntVect    ivf;

                for (int comp = 0; comp < ncomp; ++comp) {
                    for (ivsit.reset(); ivsit.ok(); ++ivsit) {
                        ivf = ivsit();

                        ivf[bdryDir] -= 2*isign;
                        p2 = fineFAB(ivf, comp);

                        ivf[bdryDir] += isign;
                        p1 = fineFAB(ivf, comp);

                        ivf[bdryDir] += isign;
                        fineFAB(ivf, comp) = c1 * p1 + c2 * p2;
                    }
                }

            } else {
                const Real factor = 1.0 - 2.0 * dxf / (dxf + dxc);
                Real       pa;
                IntVect    ivf;

                for (int comp = 0; comp < ncomp; ++comp) {
                    for (ivsit.reset(); ivsit.ok(); ++ivsit) {
                        ivf = ivsit();

                        ivf[bdryDir] -= isign;
                        pa = fineFAB(ivf, comp);

                        ivf[bdryDir] += isign;
                        fineFAB(ivf, comp) = factor * pa;
                    }
                }
            }
        } // if / if not packed
    }
}


// -----------------------------------------------------------------------------
void
CFInterp::refine(LevelData<FArrayBox>&       a_fine,
                 const LevelData<FArrayBox>& a_crse,
                 const int                   a_interpOrder,
                 const bool                  a_doSlopeLimiter) const
{
    // Sanity checks
    CH_assert(m_isDefined);
    CH_assert(a_fine.getBoxes() == m_grids);
    CH_assert(a_crse.getBoxes() == m_userCrseGrids);
    CH_assert(a_fine.nComp() == a_crse.nComp());

    // Localize the coarse data.
    LevelData<FArrayBox> localCrse(m_crseGrids, a_crse.nComp(), IntVect::Unit);
    this->localizeCrseData(localCrse, a_crse);

    // Call the workhorse.
    this->localRefine(a_fine, localCrse, a_interpOrder, a_doSlopeLimiter);
}


// -----------------------------------------------------------------------------
void
CFInterp::localRefine(LevelData<FArrayBox>&       a_fine,
                      const LevelData<FArrayBox>& a_crse,
                      const int                   a_interpOrder,
                      const bool                  a_doSlopeLimiter) const
{

    // Sanity checks
    CH_assert(m_isDefined);
    CH_assert(a_fine.getBoxes() == m_grids);
    CH_assert(a_crse.getBoxes() == m_crseGrids);
    CH_assert(a_fine.nComp() == a_crse.nComp());

    DataIterator dit = a_fine.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        this->localRefine(a_fine[dit],
                          a_crse[dit],
                          m_crseGrids[dit],
                          a_interpOrder,
                          a_doSlopeLimiter);
    }
}


// -----------------------------------------------------------------------------
void
CFInterp::localRefine(FArrayBox&       a_fineFAB,
                      const FArrayBox& a_crseFAB,
                      const Box&       a_crseInterpBox,
                      const int        a_interpOrder,
                      const bool       a_doSlopeLimiter) const
{
    // Sanity checks
    CH_assert(m_isDefined);
    CH_assert(a_fineFAB.nComp() == a_crseFAB.nComp());
    CH_assert(a_crseFAB.box().contains(a_crseInterpBox));
    CH_assert(a_fineFAB.box().contains(::refine(a_crseInterpBox, m_refRatio)));

    // This is used to effectively remove the non-contributing directions
    // from the interpolation. The non-contributors are the directions
    // in which the refRatio is 1.
    //
    // NOTE: To bring this back to the original, full dimensional interpolation,
    // just set refMask to IntVect::Unit;
    // const IntVect& refMask = IntVect::Unit;
    const IntVect refMask(D_DECL((m_refRatio[0] > 1) ? 1 : 0,
                                 (m_refRatio[1] > 1) ? 1 : 0,
                                 (m_refRatio[2] > 1) ? 1 : 0));

    const Box refbox(IntVect::Zero, m_refRatio - IntVect::Unit);
    const ProblemDomain& crseDomain = m_crseGrids.physDomain();
    const int numComps = a_crseFAB.nComp();

    // Begin with piecewise constant interpolation (order = 0).
    FORT_CFINTERP_UNMAPPEDINTERPCONSTANT(CHF_FRA(a_fineFAB),
                                         CHF_CONST_FRA(a_crseFAB),
                                         CHF_BOX(a_crseInterpBox),
                                         CHF_CONST_INTVECT(m_refRatio),
                                         CHF_BOX(refbox));

    // Upgrade to linear interpolation (order = 1).
    if (a_interpOrder >= 1) {
        // hardwired to 3 due to lack of variable number of arguments in chfpp
        BaseFab<Real> slopes[3];
        for (int dir = 0; dir < 3; ++dir) {
            BaseFab<Real>& dir_slope = slopes[dir];
            dir_slope.resize(a_crseInterpBox, numComps);
        }

        for (int dir = 0; dir < SpaceDim; ++dir) {
            BaseFab<Real>& dir_slope = slopes[dir];

            if (refMask[dir] > 0) {
                Box interiorBox = crseDomain.domainBox();
                if (!crseDomain.isPeriodic(dir)) {
                    interiorBox.grow(-BASISV(dir));
                }

                const Box bcenter = interiorBox & a_crseInterpBox;
                if (!bcenter.isEmpty()) {
                    FORT_CFINTERP_INTERPCENTRALSLOPE(CHF_FRA(dir_slope),
                                                     CHF_CONST_FRA(a_crseFAB),
                                                     CHF_BOX(bcenter),
                                                     CHF_CONST_INT(dir));
                }

                const Box blo =
                    a_crseInterpBox & adjCellLo(interiorBox, dir, 1);
                if (!blo.isEmpty()) {
                    FORT_CFINTERP_INTERPHISIDESLOPE(CHF_FRA(dir_slope),
                                                    CHF_CONST_FRA(a_crseFAB),
                                                    CHF_BOX(blo),
                                                    CHF_CONST_INT(dir));
                }

                const Box bhi =
                    a_crseInterpBox & adjCellHi(interiorBox, dir, 1);
                if (!bhi.isEmpty()) {
                    FORT_CFINTERP_INTERPLOSIDESLOPE(CHF_FRA(dir_slope),
                                                    CHF_CONST_FRA(a_crseFAB),
                                                    CHF_BOX(bhi),
                                                    CHF_CONST_INT(dir));
                }
            } else {
                // Unrefined directions should not contribute. The entire
                // refinement process should reduce in dimensionality.
                dir_slope.setVal(0.0);
            }
        } // dir

        if (a_doSlopeLimiter) {
            // To do limits, we need to have a box which includes
            // the neighbors of a given point (to check for the
            // local maximum...
            Box neighborBox(-1 * refMask, refMask);

            // GHM 7/12/01
            // interplimit iterates over box b_mod (was a_crseInterpBox), but cells
            // within 1 of the physical boundary never enter result (and this
            // wasted calculation may call upon uninitialized memory).
            // DFM 10/8/01
            // note that this turns off slope limiting for cells adjacent to the
            // boundary -- may want to revisit this in the future
            Box b_mod(a_crseInterpBox);
            b_mod.grow(refMask);
            b_mod = crseDomain.domainBox() & b_mod;
            b_mod.grow(-refMask);

            // create a box grown big enough to remove periodic BCs from domain
            Box domBox = grow(a_crseInterpBox, 2 * refMask);
            domBox &= crseDomain.domainBox();

            FORT_CFINTERP_INTERPLIMIT(CHF_FRA(slopes[0]),
                                      CHF_FRA(slopes[1]),
                                      CHF_FRA(slopes[2]),
                                      CHF_CONST_FRA(a_crseFAB),
                                      CHF_BOX(b_mod),
                                      CHF_BOX(neighborBox),
                                      CHF_BOX(domBox));
        }

        // Update to a linear interpolation (the correction).
        for (int dir = 0; dir < SpaceDim; ++dir) {
            // Don't do unrefined directions since it should add nothing anyway.
            if (m_refRatio[dir] == 1) continue;

            BaseFab<Real>& dir_slope = slopes[dir];

            FORT_CFINTERP_UNMAPPEDINTERPLINEAR(CHF_FRA(a_fineFAB),
                                               CHF_CONST_FRA(dir_slope),
                                               CHF_BOX(a_crseInterpBox),
                                               CHF_CONST_INT(dir),
                                               CHF_CONST_INTVECT(m_refRatio),
                                               CHF_BOX(refbox));
        } // dir
    } // if linear upgrade

    if (a_interpOrder >= 2) {
        MAYDAYERROR("Cannot handle a_interpOrder = " << a_interpOrder
                                                     << "yet.");
    }
}


// -----------------------------------------------------------------------------
// Static utility
void
CFInterp::localRefine(FArrayBox&       a_fineFAB,
                      const FArrayBox& a_crseFAB,
                      const Box&       a_crseInterpBox,
                      const Box&       a_crseCentralBox,
                      const IntVect&   a_refRatio,
                      const int        a_interpOrder,
                      const bool       a_doSlopeLimiter)
{
    // Sanity checks
    CH_assert(a_refRatio >= IntVect::Unit);
    CH_assert(a_fineFAB.nComp() == a_crseFAB.nComp());
    CH_assert(a_crseFAB.box().contains(a_crseInterpBox));
    CH_assert(a_fineFAB.box().contains(::refine(a_crseInterpBox, a_refRatio)));

    // This is used to effectively remove the non-contributing directions
    // from the interpolation. The non-contributors are the directions
    // in which the refRatio is 1.
    //
    // NOTE: To bring this back to the original, full dimensional interpolation,
    // just set refMask to IntVect::Unit;
    // const IntVect& refMask = IntVect::Unit;
    const IntVect refMask(D_DECL((a_refRatio[0] > 1) ? 1 : 0,
                                 (a_refRatio[1] > 1) ? 1 : 0,
                                 (a_refRatio[2] > 1) ? 1 : 0));

    const Box refbox(IntVect::Zero, a_refRatio - IntVect::Unit);
    // const ProblemDomain& crseDomain = m_crseGrids.physDomain();
    const int numComps = a_crseFAB.nComp();

    // Begin with piecewise constant interpolation (order = 0).
    FORT_CFINTERP_UNMAPPEDINTERPCONSTANT(CHF_FRA(a_fineFAB),
                                         CHF_CONST_FRA(a_crseFAB),
                                         CHF_BOX(a_crseInterpBox),
                                         CHF_CONST_INTVECT(a_refRatio),
                                         CHF_BOX(refbox));

    // Upgrade to linear interpolation (order = 1).
    if (a_interpOrder >= 1) {
        // hardwired to 3 due to lack of variable number of arguments in chfpp
        BaseFab<Real> slopes[3];
        for (int dir = 0; dir < 3; ++dir) {
            BaseFab<Real>& dir_slope = slopes[dir];
            dir_slope.resize(a_crseInterpBox, numComps);
        }

        for (int dir = 0; dir < SpaceDim; ++dir) {
            BaseFab<Real>& dir_slope = slopes[dir];

            if (refMask[dir] > 0) {
                const Box bcenter = a_crseInterpBox & a_crseCentralBox;
                if (!bcenter.isEmpty()) {
                    FORT_CFINTERP_INTERPCENTRALSLOPE(CHF_FRA(dir_slope),
                                                     CHF_CONST_FRA(a_crseFAB),
                                                     CHF_BOX(bcenter),
                                                     CHF_CONST_INT(dir));
                }

                const Box blo =
                    a_crseInterpBox & adjCellLo(a_crseCentralBox, dir, 1);
                if (!blo.isEmpty()) {
                    FORT_CFINTERP_INTERPHISIDESLOPE(CHF_FRA(dir_slope),
                                                    CHF_CONST_FRA(a_crseFAB),
                                                    CHF_BOX(blo),
                                                    CHF_CONST_INT(dir));
                }

                const Box bhi =
                    a_crseInterpBox & adjCellHi(a_crseCentralBox, dir, 1);
                if (!bhi.isEmpty()) {
                    FORT_CFINTERP_INTERPLOSIDESLOPE(CHF_FRA(dir_slope),
                                                    CHF_CONST_FRA(a_crseFAB),
                                                    CHF_BOX(bhi),
                                                    CHF_CONST_INT(dir));
                }
            } else {
                // Unrefined directions should not contribute. The entire
                // refinement process should reduce in dimensionality.
                dir_slope.setVal(0.0);
            }
        } // dir

        if (a_doSlopeLimiter) {
            // To do limits, we need to have a box which includes
            // the neighbors of a given point (to check for the
            // local maximum...
            Box neighborBox(-1 * refMask, refMask);

            // GHM 7/12/01
            // interplimit iterates over box b_mod (was a_crseInterpBox), but cells
            // within 1 of the physical boundary never enter result (and this
            // wasted calculation may call upon uninitialized memory).
            // DFM 10/8/01
            // note that this turns off slope limiting for cells adjacent to the
            // boundary -- may want to revisit this in the future
            Box b_mod(a_crseInterpBox);
            b_mod.grow(refMask);
            b_mod = a_crseCentralBox & b_mod;
            b_mod.grow(-refMask);

            // create a box grown big enough to remove periodic BCs from domain
            Box domBox = grow(a_crseInterpBox, 2 * refMask);
            domBox &= a_crseCentralBox;

            FORT_CFINTERP_INTERPLIMIT(CHF_FRA(slopes[0]),
                                      CHF_FRA(slopes[1]),
                                      CHF_FRA(slopes[2]),
                                      CHF_CONST_FRA(a_crseFAB),
                                      CHF_BOX(b_mod),
                                      CHF_BOX(neighborBox),
                                      CHF_BOX(domBox));
        }

        // Update to a linear interpolation (the correction).
        for (int dir = 0; dir < SpaceDim; ++dir) {
            // Don't do unrefined directions since it should add nothing anyway.
            if (a_refRatio[dir] == 1) continue;

            BaseFab<Real>& dir_slope = slopes[dir];

            FORT_CFINTERP_UNMAPPEDINTERPLINEAR(CHF_FRA(a_fineFAB),
                                               CHF_CONST_FRA(dir_slope),
                                               CHF_BOX(a_crseInterpBox),
                                               CHF_CONST_INT(dir),
                                               CHF_CONST_INTVECT(a_refRatio),
                                               CHF_BOX(refbox));
        } // dir
    } // if linear upgrade

    if (a_interpOrder >= 2) {
        MAYDAYERROR("Cannot handle a_interpOrder = " << a_interpOrder
                                                     << "yet.");
    }
}


// -----------------------------------------------------------------------------
void
CFInterp::coarsen(LevelData<FArrayBox>&       a_crse,
                  const LevelData<FArrayBox>& a_fine,
                  const bool                  a_doHarmonicAvg,
                  const LevelData<FArrayBox>* a_fineJPtr) const
{
    // Sanity checks
    CH_assert(m_isDefined);
    CH_assert(a_crse.getBoxes() == m_userCrseGrids);
    CH_assert(a_fine.getBoxes() == m_grids);
    CH_assert(a_fine.nComp() == a_crse.nComp());

    // Coarsen locally.
    LevelData<FArrayBox> localCrse(m_crseGrids, a_crse.nComp(), IntVect::Zero);
    this->localCoarsen(localCrse, a_fine, a_doHarmonicAvg, a_fineJPtr);

    // Copy coarsened data to user's grids.
    checkForValidNAN(a_crse);
    localCrse.copyTo(a_crse, m_crseToUserCopier);
    checkForValidNAN(a_crse);
}


// -----------------------------------------------------------------------------
void
CFInterp::localCoarsen(LevelData<FArrayBox>&       a_localCrse,
                       const LevelData<FArrayBox>& a_fine,
                       const bool                  a_doHarmonicAvg,
                       const LevelData<FArrayBox>* a_fineJPtr)
{
    const DisjointBoxLayout& fineGrids  = a_fine.getBoxes();
    const DisjointBoxLayout& crseGrids  = a_localCrse.getBoxes();
    const Box&               crseDomBox = crseGrids.physDomain().domainBox();
    const Box&               fineDomBox = fineGrids.physDomain().domainBox();
    const IntVect refRatio = calculateRefinementRatio(crseDomBox, fineDomBox);

    CH_assert(crseGrids.compatible(fineGrids));
    CH_assert(a_fine.nComp() == a_localCrse.nComp());
    CH_assert(!(refRatio <= IntVect::Unit));

    debugInitLevel(a_localCrse);

    DataIterator dit = a_localCrse.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        const FArrayBox* fineJFABPtr = nullptr;
        if (a_fineJPtr) {
            CH_assert(a_fineJPtr->getBoxes() == fineGrids);
            fineJFABPtr = &((*a_fineJPtr)[dit]);
        }

        CFInterp::localCoarsen(a_localCrse[dit],
                               a_fine[dit],
                               crseGrids[dit],
                               refRatio,
                               a_doHarmonicAvg,
                               fineJFABPtr);
    }

    checkForValidNAN(a_localCrse);
}


// -----------------------------------------------------------------------------
void
CFInterp::localCoarsen(FArrayBox&       a_crseFAB,
                       const FArrayBox& a_fineFAB,
                       const Box&       a_crseInterpBox,
                       const IntVect&   a_refRatio,
                       const bool       a_doHarmonicAvg,
                       const FArrayBox* a_fineJFABPtr)
{
    // Sanity checks
    CH_assert(a_fineFAB.nComp() == a_crseFAB.nComp());
    CH_assert(a_crseFAB.box().contains(a_crseInterpBox));
    CH_assert(a_fineFAB.box().contains(::refine(a_crseInterpBox, a_refRatio)));

    // refBox contains the offsets used to average each cell.
    Box refBox(IntVect::Zero, (a_refRatio - IntVect::Unit));

    // This is the region of fine data required to fill a_crseFAB.
    Box srcRegion;
    {
        IntVect srcRegionSmallEnd = a_crseInterpBox.smallEnd() * a_refRatio + refBox.smallEnd();
        IntVect srcRegionBigEnd = a_crseInterpBox.bigEnd() * a_refRatio + refBox.bigEnd();
        srcRegion.define(srcRegionSmallEnd, srcRegionBigEnd);
    }

    if (a_fineJFABPtr) {
        // Use volume weighting.

        // Make sure our source data is defined where we need it.
        const FArrayBox& fineJFAB = *a_fineJFABPtr;
        CH_assert(fineJFAB.box().contains(srcRegion));
        CH_assert(a_fineFAB.box().contains(srcRegion));

        if (!a_doHarmonicAvg) {
            // Arithemtic average.
            FORT_CFINTERP_MAPPEDAVERAGE (
                CHF_FRA(a_crseFAB),
                CHF_CONST_FRA(a_fineFAB),
                CHF_CONST_FRA1(fineJFAB, 0),
                CHF_BOX(a_crseInterpBox),
                CHF_CONST_INTVECT(a_refRatio),
                CHF_BOX(refBox));

        } else {
            // Harmonic average.
            MAYDAYERROR(
                "AnisotropicCoarseAverage::averageGridData -- Mapped harmonic "
                "averaging not complete");
        }
    } else {
        // No not use volume weighting. (Maybe the FAB is already
        // scaled by J.)

        // Make sure our source data is defined where we need it.
        CH_assert(a_fineFAB.box().contains(srcRegion));

        if (!a_doHarmonicAvg) {
            // Arithemtic average.
            FORT_CFINTERP_UNMAPPEDAVERAGE (
                CHF_FRA(a_crseFAB),
                CHF_CONST_FRA(a_fineFAB),
                CHF_BOX(a_crseInterpBox),
                CHF_CONST_INTVECT(a_refRatio),
                CHF_BOX(refBox));

        } else {
            // Harmonic average.
            FORT_CFINTERP_UNMAPPEDAVERAGEHARMONIC (
                CHF_FRA(a_crseFAB),
                CHF_CONST_FRA(a_fineFAB),
                CHF_BOX(a_crseInterpBox),
                CHF_CONST_INTVECT(a_refRatio),
                CHF_BOX(refBox));

        }
    }
}


// ============================ FC interpolators ===============================

#if 1
// -----------------------------------------------------------------------------
void
CFInterp::interpGhostsAtCFI(LevelData<FluxBox>&       a_fineAdvVel,
                            const LevelData<FluxBox>& a_crseAdvVel,
                            const bool                a_isHomog) const
{
    debugCheckValidFaceOverlap(a_fineAdvVel);

    // Is there anything to do?
    if (!this->hasCFI()) return;

    // Sanity checks.
    CH_assert(m_isDefined);
    CH_assert(a_fineAdvVel.getBoxes() == m_grids);
    // CH_assert(a_fineAdvVel.ghostVect() == IntVect::Unit); // For now. Modify if needed.
    CH_assert(a_isHomog || a_crseAdvVel.nComp() == a_fineAdvVel.nComp());
    CH_assert(a_isHomog || a_crseAdvVel.getBoxes() == m_userCrseGrids);

    const int      numComps  = a_fineAdvVel.nComp();
    const IntVect& numGhosts = a_fineAdvVel.ghostVect();

    // TEMPORARY!!!
    // Maintain a copy of the valid faces. We will use this later to ensure
    // they weren't overwritten.
    //
    // The only reason I need to do this is because there is no FC CFIVS!
    //
    // And yes, this produces slightly different results than a simple
    // a_fineAdvVel.copyTo(validSave). Sigh.
    debugCheckValidFaceOverlap(a_fineAdvVel);
    LevelData<FluxBox> validSave(m_grids, numComps);
    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
            FArrayBox& saveFAB = validSave[dit][fcDir];
            FArrayBox& userFAB = a_fineAdvVel[dit][fcDir];
            const Box fcValid = surroundingNodes(m_grids[dit], fcDir);

            saveFAB.copy(userFAB, fcValid);

#ifndef NDEBUG
            // Since we have validSave lying around, we might as well use
            // it to debug. This will initialize all non-exchange ghosts to NAN.
            debugInit(userFAB);
            userFAB.copy(saveFAB);
#endif
        }
    }
    // a_fineAdvVel.copyTo(validSave);
    debugCheckValidFaceOverlap(validSave);
    debugCheckValidFaceOverlap(a_fineAdvVel);

    // Send crse data to compatible grids and fill ghosts.
    LevelData<FluxBox> localCrse;
    {
        // We need enough crse ghosts to fit the stencil.
        // // If fineAdvVel has [1,ref] ghosts, then localCrse needs 2 ghosts.
        // // If fineAdvVel has [ref+1, 2*ref] ghosts, then localCrse needs 3 ghosts.
        // // etc...
        // IntVect ghostVect;
        // for (int d = 0; d < SpaceDim; ++d) {
        //     const int g = max(1, a_fineAdvVel.ghostVect()[d]);
        //     ghostVect[d] = (g - 1) / m_refRatio[d] + 2;
        // }


        // If fineAdvVel has [0,ref) ghosts, then localCrse needs 2 ghosts.
        // If fineAdvVel has [ref, 2*ref) ghosts, then localCrse needs 3 ghosts.
        // etc...
        const IntVect ghostVect = a_fineAdvVel.ghostVect() / m_refRatio + 2;
        localCrse.define(m_crseGrids, numComps, ghostVect);

        if (a_isHomog) {
            for (DataIterator dit(m_crseGrids); dit.ok(); ++dit) {
                localCrse[dit].setVal(0.0);
            }
        } else {
            debugCheckValidFaceOverlap(a_crseAdvVel);
            this->localizeCrseData(localCrse, a_crseAdvVel);
        }
    }

    // Average data down to coarse level. Do NOT assume data is div free.
    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        CFInterp::localCoarsen(localCrse[dit],
                               a_fineAdvVel[dit],
                               m_crseGrids[dit],
                               m_refRatio,
                               nullptr); // advVel is already scaled by J.
    }
    debugCheckValidFaceOverlap(localCrse); // This shouldn't trip. We'll see.

    // Loop through CFI.
    const auto key = m_cfiIter.lock();
    for (m_cfiIter.reset(key); m_cfiIter.ok(); m_cfiIter.next(key)) {
        const DataIndex&      di            = m_cfiIter->di;
        const Box&            fcBdryBox     = m_cfiIter->fcBdryBox;
        const Box&            ccAdjGhostBox = m_cfiIter->ccAdjGhostBox;
        const int             bdryDir       = m_cfiIter->dir;
        const Side::LoHiSide& side          = m_cfiIter->side;
        const int             isign         = m_cfiIter->isign;
        const CFIVS&          cfivs         = m_cfiIter->cfivs;

        const IntVect eb = BASISV(bdryDir);
        const IntVect et = IntVect::Unit - eb;

        // The ghosts.
        Box fineGhostBox = ccAdjGhostBox;
        fineGhostBox.grow(et * numGhosts); // Comment out to eliminate corners.
        fineGhostBox.growDir(bdryDir, side, numGhosts[bdryDir] - 1);
        fineGhostBox &= a_fineAdvVel[di].box();

        // The coarse-level shadow of the fine-level ghosts.
        Box crseInterpBox = fineGhostBox;
        crseInterpBox.coarsen(m_refRatio);
        CH_assert(!crseInterpBox.isEmpty());

        // This fine patch, with refRatio ghosts to cover crseInterpBox.
        Box fineInterpBox = crseInterpBox;
        fineInterpBox.refine(m_refRatio);
        CH_assert(!fineInterpBox.isEmpty());

        for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
            FArrayBox&       fineFAB  = a_fineAdvVel[di][fcDir];
            const FArrayBox& crseFAB  = localCrse[di][fcDir];
            const Box fcCrseInterpBox = surroundingNodes(crseInterpBox, fcDir);
            const Box fcFineInterpBox = surroundingNodes(fineInterpBox, fcDir);

            // Step 1: Interp to parent faces.
            FArrayBox fineInterpFAB(fcFineInterpBox, numComps);
            if (SpaceDim == 2) {
                this->interpAlongParentEdgesIn2D(fineInterpFAB,
                                                 crseFAB,
                                                 fcCrseInterpBox,
                                                 1 - fcDir,  // interpDir
                                                 m_refRatio);
            } else {
                this->interpAlongParentFacesIn3D(
                    fineInterpFAB, crseFAB, fcCrseInterpBox, fcDir, m_refRatio);
            }

            // Step 1.5: Replace with more accurate, valid fine data.
            if (fcDir == bdryDir) {
                CH_assert(fineInterpFAB.box().contains(fcBdryBox));
                CH_assert(fineFAB.box().contains(fcBdryBox));
                fineInterpFAB.copy(
                    fineFAB, fcBdryBox, 0, fcBdryBox, 0, numComps);
            }

            // Step 2: Interpolate from parent faces to interior faces.
            {
                // Offsets from the lower parent face to the interior faces.
                Box interiorRefBox(IntVect::Zero, m_refRatio - IntVect::Unit);
                interiorRefBox.surroundingNodes(fcDir);
                interiorRefBox.grow(fcDir, -1);

                // Do it!
                FORT_CFINTERP_UNMAPPEDINTERPLINEARINTERIORFACE(
                    CHF_FRA(fineInterpFAB),
                    CHF_BOX(crseInterpBox),  // Yes, the CC version!
                    CHF_CONST_INTVECT(m_refRatio),
                    CHF_CONST_INT(fcDir),
                    CHF_BOX(interiorRefBox)
                );
            }

            // Copy ghost faces to user's fine holder.
            if (cfivs.isPacked()) {
                Box copyBox = cfivs.packedBox();
                copyBox.grow(et * numGhosts);
                if (side == Side::Lo) {
                    copyBox.growLo(bdryDir, numGhosts[bdryDir]);
                } else {
                    copyBox.growHi(bdryDir, numGhosts[bdryDir]);
                }
                if (fcDir == bdryDir) {
                    copyBox.shiftHalf(fcDir, isign);
                } else {
                    copyBox.surroundingNodes(fcDir);
                }
                copyBox &= fineFAB.box();
                fineFAB.copy(fineInterpFAB, copyBox);
            } else {
                // TEMPORARY!!!
                // If I ever create a proper FC CFIVS, I will fix up the
                // commented code below. Until then, just clobber the valid
                // data, then restore it later with validSave.

                // cout << "!isPacked" << endl;

                Box copyBox = cfivs.minBox();
                copyBox.grow(et * numGhosts);
                if (side == Side::Lo) {
                    copyBox.growLo(bdryDir, numGhosts[bdryDir]);
                } else {
                    copyBox.growHi(bdryDir, numGhosts[bdryDir]);
                }
                if (fcDir == bdryDir) {
                    copyBox.shiftHalf(fcDir, isign);
                } else {
                    copyBox.surroundingNodes(fcDir);
                }
                copyBox &= fineFAB.box();
                fineFAB.copy(fineInterpFAB, copyBox);

                // // This is the weakest link because it can only handle one
                // // ghost and cannot handle corners. In theory, if you ever need
                // // multiple ghost layers or corners, you should just need to
                // // modify this bit of code.
                // IVSIterator ivsit = cfivs.getIVS();
                // if (fcDir == bdryDir) {
                //     if (side == Side::Lo) {
                //         for (int comp = 0; comp < numComps; ++comp) {
                //             for (ivsit.reset(); ivsit.ok(); ++ivsit) {
                //                 const IntVect& iv = ivsit();
                //                 fineFAB(iv, comp) = fineInterpFAB(iv, comp);
                //             }
                //         }
                //     } else {
                //         IntVect iv;
                //         for (int comp = 0; comp < numComps; ++comp) {
                //             for (ivsit.reset(); ivsit.ok(); ++ivsit) {
                //                 iv = ivsit() + eb;
                //                 fineFAB(iv, comp) = fineInterpFAB(iv, comp);
                //             }
                //         }
                //     }
                // } else {
                //     BUG("I need to think about this...");
                // }
            } // if / if not packed
        } // fcDir
    } // m_cfiIter
    m_cfiIter.unlock(key);

    // TEMPORARY!!! Restore valid data.
    // And yes, this produces slightly different results than a simple
    // validSave.copyTo(a_fineAdvVel). Sigh.
    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
            FArrayBox&       userFAB = a_fineAdvVel[dit][fcDir];
            const FArrayBox& saveFAB = validSave[dit][fcDir];
            const Box        fcValid = surroundingNodes(m_grids[dit], fcDir);

            userFAB.copy(saveFAB, fcValid);
        }
    }
    // validSave.copyTo(a_fineAdvVel);

    debugCheckValidFaceOverlap(a_fineAdvVel);
}

#else

// -----------------------------------------------------------------------------
void
CFInterp::interpGhostsAtCFI(LevelData<FluxBox>&       a_fineAdvVel,
                            const LevelData<FluxBox>& a_crseAdvVel,
                            const bool                a_isHomog) const
{
    debugCheckValidFaceOverlap(a_fineAdvVel);

    // Is there anything to do?
    if (!this->hasCFI()) return;

    // Sanity checks.
    CH_assert(m_isDefined);
    CH_assert(a_fineAdvVel.nComp() == 1);
    CH_assert(a_fineAdvVel.getBoxes() == m_grids);
    CH_assert(a_fineAdvVel.ghostVect() == IntVect::Unit); // For now. Modify if needed.
    CH_assert(a_isHomog || a_crseAdvVel.nComp() == 1);
    CH_assert(a_isHomog || a_crseAdvVel.getBoxes() == m_userCrseGrids);

    const int      numComps  = 1;
    const IntVect& numGhosts = a_fineAdvVel.ghostVect();

    // Send crse data to compatible grids and fill ghosts.
    LevelData<FluxBox> localCrse;
    {
        const IntVect ghostVect = 3 * IntVect::Unit;
        localCrse.define(m_crseGrids, numComps, ghostVect);

        if (a_isHomog) {
            for (DataIterator dit(m_crseGrids); dit.ok(); ++dit) {
                localCrse[dit].setVal(0.0);
            }
        } else {
            debugCheckValidFaceOverlap(a_crseAdvVel);
            this->localizeCrseData(localCrse, a_crseAdvVel);
            debugCheckValidFaceOverlap(localCrse);
        }
    }

    // Average data down to coarse level. Do NOT assume data is div free.
    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        CFInterp::localCoarsen(localCrse[dit],
                               a_fineAdvVel[dit],
                               m_crseGrids[dit],
                               m_refRatio,
                               nullptr); // advVel is already scaled by J.
    }
    debugCheckValidFaceOverlap(localCrse); // This shouldn't trip. We'll see.

    // This version interpolates over a grown valid box, which makes it
    // wasteful but simple to code up.
    // Honestly, I don't think this function is called enough for efficiency
    // to be an issue. Still, I coded up the more efficient version because
    // it will help when I re-instate the div-free interp code.

    LevelData<FluxBox> fineInterp(m_grids, numComps, m_refRatio);
    debugInitLevel(fineInterp);

    // Interp to parent faces.
    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        FluxBox&       fineFB = fineInterp[dit];
        const FluxBox& crseFB = localCrse[dit];

        // The shadow of this fine patch, including the ghost.
        Box crseInterpBox = fineFB.box();
        crseInterpBox.coarsen(m_refRatio);

        this->localDivFreeInterp_Step1(fineFB, crseFB, crseInterpBox);
    }

    // Copy whatever valid data the user wishes to keep around.
    if (true) { //a_fineDataPtr) {
        debugCheckValidFaceOverlap(a_fineAdvVel);

        // Copy.
        a_fineAdvVel.copyTo(fineInterp);

        // Repair the vertical interpolation of w.
        if (SpaceDim > 2) {
            for (DataIterator dit(m_grids); dit.ok(); ++dit) {
                FluxBox&       fineFB = fineInterp[dit];
                const FluxBox& crseFB = localCrse[dit];

                // The shadow of this fine patch, including the ghost.
                Box crseInterpBox = fineFB.box();
                crseInterpBox.coarsen(m_refRatio);

                Box crseCCValid = m_grids[dit];
                crseCCValid.coarsen(m_refRatio);
                this->vertVelInterp(
                    fineFB[2], crseFB[2], crseCCValid);
            }
        }

        fineInterp.exchange(); // TODO: Is this needed? If so, use a copier.
    }

    // Interp to parent interiors.
    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        FluxBox&       fineFB = fineInterp[dit];
        const FluxBox& crseFB = localCrse[dit];

        // The shadow of this fine patch, including the ghost.
        Box crseInterpBox = fineFB.box();
        crseInterpBox.coarsen(m_refRatio);

        this->localDivFreeInterp_Step2(fineFB, crseFB, crseInterpBox);
    }

    // Repair the interiors that were overwritten in step 2.
    if (true) { //a_fineDataPtr) {
        a_fineAdvVel.copyTo(fineInterp);
        // if (a_fineDataCopierPtr) {
        //     a_fineDataPtr->copyTo(a_fine, *a_fineDataCopierPtr);
        // } else {
        //     parentCopier.define(a_fineDataPtr->getBoxes(), m_grids, m_domain);
        //     a_fineDataPtr->copyTo(a_fine, parentCopier);
        // }
    }

    fineInterp.exchange();
    debugCheckValidFaceOverlap(fineInterp);

    // Copy ghost faces to user's fine holder.
    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
            Box validFCBox = m_grids[dit];
            validFCBox.surroundingNodes(fcDir);

            fineInterp[dit][fcDir].copy(a_fineAdvVel[dit][fcDir], validFCBox);
            a_fineAdvVel[dit][fcDir].copy(fineInterp[dit][fcDir]);
        }
    } // dit
    a_fineAdvVel.exchange();

    debugCheckValidFaceOverlap(a_fineAdvVel);
}
#endif


#ifdef ALLOW_DIVFREEINTERP
// -----------------------------------------------------------------------------
void
CFInterp::divFreeInterp(LevelData<FluxBox>&       a_fine,
                        const LevelData<FluxBox>& a_crse,
                        const LevelData<FluxBox>* a_fineDataPtr,
                        const Copier*             a_fineDataCopierPtr) const
{
    CH_assert(m_isDefined);

    MAYDAYERROR("divFreeInterp needs to take in J ptrs.");

    // Localize the coarse data, then call the workhorse.
    LevelData<FluxBox> localCrse(m_crseGrids, a_crse.nComp(), IntVect::Unit);
    this->localizeCrseData(localCrse, a_crse);
    this->localDivFreeInterp(
        a_fine, localCrse, a_fineDataPtr, a_fineDataCopierPtr);

    debugCheckValidFaceOverlap(a_fine);
}
#endif // ALLOW_DIVFREEINTERP


// -----------------------------------------------------------------------------
void
CFInterp::homogInterpAtCFI(LevelData<FluxBox>&                   a_fine,
                           const RealVect&                       a_fineDXi,
                           const RealVect&                       a_crseDXi,
                           const std::array<CFRegion, SpaceDim>& a_cfRegion)
{
    const DisjointBoxLayout& grids = a_fine.getBoxes();
    const RealVect refRatio = a_crseDXi / a_fineDXi;

    for (int velComp = 0; velComp < SpaceDim; ++velComp) {
        if (a_cfRegion[velComp].isEmpty()) continue;

        for (DataIterator dit(grids); dit.ok(); ++dit) {
            FArrayBox& uFAB = a_fine[dit][velComp];

            for (int bdryDir = 0; bdryDir < SpaceDim; ++bdryDir) {
                const Real r = refRatio[bdryDir];

                for (SideIterator sit; sit.ok(); ++sit) {
                    const CFIVS& cfivs =
                        a_cfRegion[velComp].getCFIVS(dit(), bdryDir, sit());
                    if (cfivs.isEmpty()) continue;

                    const int isign = sign(sit());
                    IVSIterator ivsit(cfivs.getIVS());
                    const Real  cb = 1.0 - 1.0 / r;

                    // const Real  qb = 2.0 * cb;
                    // const Real  qn = -(r - 1.0) / (r + 1.0);

                    for (int comp = 0; comp < a_fine.nComp(); ++comp) {
                        // Linear
                        for (ivsit.reset(); ivsit.ok(); ++ivsit) {
                            IntVect ivf = ivsit(); // ghost

                            ivf[bdryDir] -= isign; // first valid face
                            const Real ub = uFAB(ivf, comp);

                            ivf[bdryDir] += isign; // Back to ghost
                            uFAB(ivf, comp) = cb * ub;
                        } // ivsit

                        // // Quadratic
                        // for (ivsit.reset(); ivsit.ok(); ++ivsit) {
                        //     IntVect ivf = ivsit(); // ghost

                        //     ivf[bdryDir] -= isign; // first valid face
                        //     const Real ub = uFAB(ivf, comp);

                        //     ivf[bdryDir] -= isign; // second valid face
                        //     const Real un = uFAB(ivf, comp);

                        //     ivf[bdryDir] += 2*isign; // Back to ghost
                        //     uFAB(ivf, comp) = qb*ub + qn*un;
                        // } // ivsit
                    } // comp
                } // sit
            } // bdryDir
        } // dit
    } // velComp (= fcDir)
}


// -----------------------------------------------------------------------------
void
CFInterp::refine(LevelData<FluxBox>&       a_fineAdvVel,
                 const LevelData<FluxBox>& a_crseAdvVel) const
{
    CH_assert(m_isDefined);
    CH_assert(a_fineAdvVel.getBoxes() == m_grids);
    CH_assert(a_crseAdvVel.getBoxes() == m_userCrseGrids);
    CH_assert(a_crseAdvVel.nComp() == a_fineAdvVel.nComp());

    // Send crse data to compatible grids and fill ghosts.
    LevelData<FluxBox> localCrse;
    {
        const IntVect ghostVect = IntVect::Unit;
        localCrse.define(m_crseGrids, a_crseAdvVel.nComp(), ghostVect);
        debugCheckValidFaceOverlap(a_crseAdvVel);
        this->localizeCrseData(localCrse, a_crseAdvVel);
    }

    // Call workhorse.
    this->localRefine(a_fineAdvVel, localCrse);
    debugCheckValidFaceOverlap(a_fineAdvVel);
}


// -----------------------------------------------------------------------------
void
CFInterp::localRefine(LevelData<FluxBox>&       a_fineAdvVel,
                      const LevelData<FluxBox>& a_crseAdvVel) const
{
    CH_assert(m_isDefined);
    CH_assert(a_fineAdvVel.getBoxes() == m_grids);
    CH_assert(a_crseAdvVel.getBoxes() == m_crseGrids);
    CH_assert(a_crseAdvVel.nComp() == a_fineAdvVel.nComp());

    DataIterator dit = m_grids.dataIterator();

    for (dit.reset(); dit.ok(); ++dit) {
        Box ccCrseInterpBox = m_grids[dit];
        ccCrseInterpBox.coarsen(m_refRatio);

        this->localRefine(a_fineAdvVel[dit],
                          a_crseAdvVel[dit],
                          ccCrseInterpBox);
    }
}


// -----------------------------------------------------------------------------
void
CFInterp::localRefine(FluxBox&       a_fineAdvVelFB,
                      const FluxBox& a_crseAdvVelFB,
                      const Box&     a_ccCrseInterpBox) const
{
    CH_assert(m_isDefined);
    CH_assert(a_fineAdvVelFB.nComp() == a_crseAdvVelFB.nComp());
    CH_assert(a_crseAdvVelFB.box().contains(a_ccCrseInterpBox));

    // Box fineInterpBox = crseInterpBox;
    // fineInterpBox.refine(m_refRatio);

    for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
        FArrayBox&       fineFAB  = a_fineAdvVelFB[fcDir];
        const FArrayBox& crseFAB  = a_crseAdvVelFB[fcDir];
        const Box fcCrseInterpBox = surroundingNodes(a_ccCrseInterpBox, fcDir);

        CFInterp::localRefineFace(
            fineFAB, crseFAB, fcCrseInterpBox, fcDir, m_refRatio);
    }
}


// -----------------------------------------------------------------------------
void
CFInterp::localRefineFace(FArrayBox&       a_fineAdvVelFAB,
                          const FArrayBox& a_crseAdvVelFAB,
                          const Box&       a_fcCrseInterpBox,
                          const int        a_fcDir,
                          const IntVect&   a_refRatio)
{
    CH_assert(a_fineAdvVelFAB.box().type() == BASISV(a_fcDir));
    CH_assert(a_crseAdvVelFAB.box().type() == BASISV(a_fcDir));
    CH_assert(a_fcCrseInterpBox.type() == BASISV(a_fcDir));
    CH_assert(a_crseAdvVelFAB.box().contains(a_fcCrseInterpBox));
    CH_assert(a_fineAdvVelFAB.nComp() == a_crseAdvVelFAB.nComp());

    // Conservative interpolation to parent faces.
    if (SpaceDim == 2) {
        CFInterp::interpAlongParentEdgesIn2D(a_fineAdvVelFAB,
                                             a_crseAdvVelFAB,
                                             a_fcCrseInterpBox,
                                             1 - a_fcDir,  // interpDir
                                             a_refRatio);
    } else {
        CFInterp::interpAlongParentFacesIn3D(a_fineAdvVelFAB,
                                             a_crseAdvVelFAB,
                                             a_fcCrseInterpBox,
                                             a_fcDir,
                                             a_refRatio);
    }

    // Interpolate from parent faces to interior faces...

    // Offsets from the lower parent face to the interior faces.
    Box interiorRefBox(IntVect::Zero, a_refRatio - IntVect::Unit);
    interiorRefBox.surroundingNodes(a_fcDir);
    interiorRefBox.grow(a_fcDir, -1);

    Box ccCrseInterpBox = a_fcCrseInterpBox;
    ccCrseInterpBox.enclosedCells(a_fcDir);

    // Do it!
    FORT_CFINTERP_UNMAPPEDINTERPLINEARINTERIORFACE(
        CHF_FRA(a_fineAdvVelFAB),
        CHF_BOX(ccCrseInterpBox),
        CHF_CONST_INTVECT(a_refRatio),
        CHF_CONST_INT(a_fcDir),
        CHF_BOX(interiorRefBox));
}


// -----------------------------------------------------------------------------
void
CFInterp::coarsen(LevelData<FluxBox>&       a_crse,
                  const LevelData<FluxBox>& a_fine,
                  const LevelData<FluxBox>* a_fineJPtr) const
{
    CH_assert(m_isDefined);
    CH_assert(a_fine.getBoxes() == m_grids);
    CH_assert(a_crse.getBoxes() == m_userCrseGrids);
    CH_assert(a_fine.nComp() == a_crse.nComp());

    // Send crse data to compatible grids.
    LevelData<FluxBox> localCrse(m_crseGrids, a_crse.nComp());
    CFInterp::localCoarsen(localCrse, a_fine, a_fineJPtr);

    // Send localCrse to user's holder.
    localCrse.copyTo(a_crse);

    debugCheckValidFaceOverlap(a_crse);
}


// -----------------------------------------------------------------------------
void
CFInterp::localCoarsen(LevelData<FluxBox>&       a_localCrse,
                       const LevelData<FluxBox>& a_fine,
                       const LevelData<FluxBox>* a_fineJPtr)
{
    const DisjointBoxLayout& fineGrids  = a_fine.getBoxes();
    const DisjointBoxLayout& crseGrids  = a_localCrse.getBoxes();
    const Box&               crseDomBox = crseGrids.physDomain().domainBox();
    const Box&               fineDomBox = fineGrids.physDomain().domainBox();
    const IntVect refRatio = calculateRefinementRatio(crseDomBox, fineDomBox);

    CH_assert(crseGrids.compatible(fineGrids));
    CH_assert(a_fine.nComp() == a_localCrse.nComp());
    CH_assert(refRatio >= IntVect::Unit);

    debugInitLevel(a_localCrse);
    debugCheckValidFaceOverlap(a_fine);

    DataIterator dit = a_fine.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        const FluxBox* fineJFlubPtr = nullptr;
        if (a_fineJPtr) {
            CH_assert(a_fineJPtr->getBoxes() == fineGrids);
            fineJFlubPtr = &((*a_fineJPtr)[dit]);
        }

        CFInterp::localCoarsen(a_localCrse[dit],
                               a_fine[dit],
                               crseGrids[dit],
                               refRatio,
                               fineJFlubPtr);
    }

    checkForValidNAN(a_localCrse);
    debugCheckValidFaceOverlap(a_localCrse);
}


// -----------------------------------------------------------------------------
void
CFInterp::localCoarsen(FluxBox&       a_crse,
                       const FluxBox& a_fine,
                       const Box&     a_crseCCValid,
                       const IntVect& a_refRatio,
                       const FluxBox* a_fineJPtr)
{
    CH_assert(a_fine.nComp() == a_crse.nComp());
    CH_assert(a_crse.box().contains(a_crseCCValid));
    CH_assert(a_fine.box().contains( Box(a_crseCCValid).refine(a_refRatio) ));

    for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
        FArrayBox&       crseFAB = a_crse[fcDir];
        const FArrayBox& fineFAB = a_fine[fcDir];
        const int        ncomp   = crseFAB.nComp();

        Box crseFCValid = surroundingNodes(a_crseCCValid, fcDir);
        Box refBox;
        {
            IntVect hiVect(D_DECL(a_refRatio[0]-1,
                                  a_refRatio[1]-1,
                                  a_refRatio[2]-1));
            hiVect[fcDir] = 0;
            refBox.define(IntVect::Zero, hiVect, BASISV(fcDir));
        }

        // We must shift to avoid negative numbers.
        const IntVect crseShift = crseFAB.box().smallEnd();
        const IntVect fineShift = crseShift * a_refRatio;

        // Set equal to averaged fine data.
        crseFAB.setVal(0.0, crseFCValid, 0, ncomp);
        if (a_fineJPtr) {
            const FArrayBox& fineJFAB = (*a_fineJPtr)[fcDir];

            FORT_CFINTERP_MAPPEDAVERAGEFACE(
                CHF_FRA_SHIFT(crseFAB, crseShift),
                CHF_CONST_FRA_SHIFT(fineFAB, fineShift),
                CHF_CONST_FRA1_SHIFT(fineJFAB, 0, fineShift),
                CHF_BOX_SHIFT(crseFCValid, crseShift),
                CHF_CONST_INT(fcDir),
                CHF_CONST_INTVECT(a_refRatio),
                CHF_BOX(refBox));

        } else {
            FORT_CFINTERP_UNMAPPEDAVERAGEFACE(
                CHF_FRA_SHIFT(crseFAB, crseShift),
                CHF_CONST_FRA_SHIFT(fineFAB, fineShift),
                CHF_BOX_SHIFT(crseFCValid, crseShift),
                CHF_CONST_INT(fcDir),
                CHF_CONST_INTVECT(a_refRatio),
                CHF_BOX(refBox));
        }
    }
}


// -----------------------------------------------------------------------------
void
CFInterp::localCoarsenFace(FArrayBox&       a_crse,
                           const FArrayBox& a_fine,
                           const Box&       a_crseFCValid,
                           const int        a_fcDir,
                           const IntVect&   a_refRatio,
                           const FArrayBox* a_fineJPtr)
{
    CH_assert(0 <= a_fcDir && a_fcDir < SpaceDim);
    CH_assert(a_crseFCValid.type() == BASISV(a_fcDir));
    CH_assert(a_fine.nComp() == a_crse.nComp());
    CH_assert(a_crse.box().contains(a_crseFCValid));
    CH_assert(a_fine.box().contains( Box(a_crseFCValid).refine(a_refRatio) ));

    Box refBox;
    {
        IntVect hiVect(D_DECL(a_refRatio[0]-1,
                              a_refRatio[1]-1,
                              a_refRatio[2]-1));
        hiVect[a_fcDir] = 0;
        refBox.define(IntVect::Zero, hiVect, BASISV(a_fcDir));
    }

    // We must shift to avoid negative numbers.
    const IntVect crseShift = a_crse.box().smallEnd();
    const IntVect fineShift = crseShift * a_refRatio;

    // Set equal to averaged fine data.
    a_crse.setVal(0.0, a_crseFCValid, 0, a_crse.nComp());
    if (a_fineJPtr) {
        const FArrayBox& fineJFAB = *a_fineJPtr;

        FORT_CFINTERP_MAPPEDAVERAGEFACE(
            CHF_FRA_SHIFT(a_crse, crseShift),
            CHF_CONST_FRA_SHIFT(a_fine, fineShift),
            CHF_CONST_FRA1_SHIFT(fineJFAB, 0, fineShift),
            CHF_BOX_SHIFT(a_crseFCValid, crseShift),
            CHF_CONST_INT(a_fcDir),
            CHF_CONST_INTVECT(a_refRatio),
            CHF_BOX(refBox));

    } else {
        FORT_CFINTERP_UNMAPPEDAVERAGEFACE(
            CHF_FRA_SHIFT(a_crse, crseShift),
            CHF_CONST_FRA_SHIFT(a_fine, fineShift),
            CHF_BOX_SHIFT(a_crseFCValid, crseShift),
            CHF_CONST_INT(a_fcDir),
            CHF_CONST_INTVECT(a_refRatio),
            CHF_BOX(refBox));
    }
}


// -----------------------------------------------------------------------------
void
CFInterp::localCoarsenEdge(FArrayBox&       a_crse,
                           const FArrayBox& a_fine,
                           const Box&       a_crseECValid,
                           const int        a_ecDir,
                           const IntVect&   a_refRatio,
                           const FArrayBox* a_fineJPtr)
{
    CH_assert(0 <= a_ecDir && a_ecDir < SpaceDim);
    CH_assert(a_crseECValid.type() == IntVect::Unit - BASISV(a_ecDir));
    CH_assert(a_fine.nComp() == a_crse.nComp());
    CH_assert(a_crse.box().contains(a_crseECValid));
    CH_assert(a_fine.box().contains( Box(a_crseECValid).refine(a_refRatio) ));

    Box refBox;
    {
        IntVect hiVect(D_DECL(0, 0, 0));
        hiVect[a_ecDir] = a_refRatio[a_ecDir] - 1;

        IntVect boxType = IntVect::Unit;
        boxType[a_ecDir] = 0;

        refBox.define(IntVect::Zero, hiVect, boxType);
    }

    // We must shift to avoid negative numbers.
    const IntVect crseShift = a_crse.box().smallEnd();
    const IntVect fineShift = crseShift * a_refRatio;

    // Set equal to averaged fine data.
    a_crse.setVal(0.0, a_crseECValid, 0, a_crse.nComp());
    if (a_fineJPtr) {
        const FArrayBox& fineJFAB = *a_fineJPtr;

        FORT_CFINTERP_MAPPEDAVERAGEEDGE(
            CHF_FRA_SHIFT(a_crse, crseShift),
            CHF_CONST_FRA_SHIFT(a_fine, fineShift),
            CHF_CONST_FRA1_SHIFT(fineJFAB, 0, fineShift),
            CHF_BOX_SHIFT(a_crseECValid, crseShift),
            CHF_CONST_INT(a_ecDir),
            CHF_CONST_INTVECT(a_refRatio),
            CHF_BOX(refBox));

    } else {
        FORT_CFINTERP_UNMAPPEDAVERAGEEDGE(
            CHF_FRA_SHIFT(a_crse, crseShift),
            CHF_CONST_FRA_SHIFT(a_fine, fineShift),
            CHF_BOX_SHIFT(a_crseECValid, crseShift),
            CHF_CONST_INT(a_ecDir),
            CHF_CONST_INTVECT(a_refRatio),
            CHF_BOX(refBox));
    }
}


// -----------------------------------------------------------------------------
void
CFInterp::localCoarsenNode(FArrayBox&       a_crse,
                           const FArrayBox& a_fine,
                           const Box&       a_crseNCValid,
                           const IntVect&   a_refRatio)
{
    CH_assert(a_crseNCValid.type() == IntVect::Unit);
    CH_assert(a_fine.nComp() == a_crse.nComp());
    CH_assert(a_crse.box().contains(a_crseNCValid));
    CH_assert(a_fine.box().contains( Box(a_crseNCValid).refine(a_refRatio) ));

    // We must shift to avoid negative numbers.
    const IntVect crseShift = a_crse.box().smallEnd();
    const IntVect fineShift = crseShift * a_refRatio;

    // Set equal to averaged fine data.
    a_crse.setVal(0.0, a_crseNCValid, 0, a_crse.nComp());
    FORT_CFINTERP_INTERPTOCOARSENODES(
        CHF_FRA_SHIFT(a_crse, crseShift),
        CHF_CONST_FRA_SHIFT(a_fine, fineShift),
        CHF_BOX_SHIFT(a_crseNCValid, crseShift),
        CHF_CONST_INTVECT(a_refRatio));
}


// ============================ Helper functions ===============================

// -----------------------------------------------------------------------------
void
CFInterp::validateAtCoarseCFI(LevelData<FluxBox>& a_crse,
                              const BoxLayout&    a_fineLayout,
                              const IntVect&      a_refRatio)
{
    // This function has been deprecated since StaggeredCopier was created.
    // UNDEFINED_FUNCTION();

    // CH_assert(a_refRatio >= IntVect::Unit);
    // CH_assert(a_refRatio.product() > 1);
    CH_assert(a_refRatio >= IntVect::Unit);

    const int                numComps = a_crse.nComp();
    const DisjointBoxLayout& grids    = a_crse.getBoxes();
    DataIterator             dit      = grids.dataIterator();

    for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
        // Suppose leftFlub and rightFlub share an exchange boundary.
        //             leftFlub           rightFlub
        //       |-----|-----|-----U
        //                         L-----|-----|-----|
        // where L and U represents data that is potentially different.
        // Also, suppose rightFlub is covered by finer data and leftFlub is not.
        // Then L is valid data since it came from the finer level and needs to
        // override U.
        // This situation is identified and corrected in the following steps.
        //
        // 1. We send FC crse data a half-cell toward the valid interior.
        //       leftFlub and leftFAB
        //       ?-----|-----|-----U
        // |.....|--?--|-----|--U--|.....|
        //
        //                        rightFlub and rightFAB
        //                         L-----|-----|-----?
        //                   |.....|--L--|-----|--?--|.....|
        // The ?s represent data not at an exchange boundary.
        // Also, the FABs have a ghost layer drawn as |.....|
        //
        // 2. We perform an exchange.
        // Since leftFAB and rightFAB share an exchange boundary, this
        //             leftFAB           rightFAB
        // |.....|--?--|-----|--U--|.....|
        //                   |.....|--L--|-----|--?--|.....|
        //
        // will become this
        //             leftFAB           rightFAB
        // |.....|--?--|-----|--U--|..L..|
        //                   |..U..|--L--|-----|--?--|.....|
        //
        // 3. Create a mask. 1 at cells covered by finer grids, 0 elsewhere.
        // The mask does not have a ghost layer.
        // Graphically, if rightFAB is covered and leftFAB is not,
        // then the mask will be
        //                             fine grids
        //                         |--|--|--|--|--|--|
        //        leftFAB and mask
        // |..0..|--0--|--0--|--0--|..1..|
        // |.....|--?--|-----|--U--|..L..|
        //                          rightFAB and mask
        //                   |..0..|--1--|--1--|--1--|..0..|
        //                   |..U..|--L--|-----|--?--|.....|
        // Note that at this point, L is valid and masked with a 1
        // while U is invalid and masked with a 0. This is true on both FABs.
        //
        // 4. On each FAB, search for adjacent cells with different mask values,
        // aka adjacent valid/invalid data. Copy the valid data in the CC FAB
        // to the boundary face in the original FC holder.
        //
        // On the left side...
        // Before:
        //   left maskFAB  |..0..|--0--|--0--|--0--|..1..|
        //   leftFAB       |.....|--?--|-----|--U--|..L..|
        //
        //   leftFlub            ?-----|-----|-----U
        //
        // After:
        //   left maskFAB  |..0..|--0--|--0--|--0--|..1..|
        //   leftFAB       |.....|--?--|-----|--U--|..L..|
        //                                          <─┘ copy this way
        //   leftFlub            ?-----|-----|-----L
        //
        // On the right side...
        // Before:
        //   right maskFAB  |..0..|--1--|--1--|--1--|..0..|
        //   rightFAB       |..U..|--L--|-----|--?--|.....|
        //
        //   rightFlub            |--L--|-----|-----?
        //
        // After:
        //   right maskFAB  |..0..|--1--|--1--|--1--|..0..|
        //   rightFAB       |..U..|--L--|-----|--?--|.....|
        //           copy this way <─┘           └─>
        //   rightFlub            L-----|-----|-----?

        LevelData<FArrayBox> exData(grids, numComps, BASISV(fcDir));
        debugInitLevel(exData);

        // 1. Copy bdry faces to adjacent interior cells.
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox&       exFAB   = exData[dit];
            const FArrayBox& crseFAB = a_crse[dit][fcDir];
            const Box&       ccValid = grids[dit];

            for (SideIterator sit; sit.ok(); ++sit) {
                const Side::LoHiSide& side  = sit();
                const int             isign = sign(side);

                const Box copyBox = bdryBox(ccValid, fcDir, side);

                exFAB.shiftHalf(fcDir, isign);
                exFAB.copy(crseFAB, copyBox, 0, copyBox, 0, numComps);
                exFAB.shiftHalf(fcDir, -isign);
            } // sit
        } // dit

        // 2. Exchange. exData will have data on each side of the bdry faces.
        // The interior cells will contain this patch's face values.
        // The exterior cells will contain the neighbor's face values.
        {
            Copier cp;
            cp.exchangeDefine(grids, exData.ghostVect());
            exData.exchange(cp);
        }

        // 3. Create a mask.
        // 0 at uncovered cells. 1 under the fine level.
        LevelData<BaseFab<int>> ccMask(grids, 1, exData.ghostVect());
        for (dit.reset(); dit.ok(); ++dit) {
            BaseFab<int>& maskFAB = ccMask[dit];
            const Box&    maskBox = maskFAB.box();

            maskFAB.setVal(0);

            LayoutIterator fineLit = a_fineLayout.layoutIterator();
            for (fineLit.reset(); fineLit.ok(); ++fineLit) {
                Box coveredBox = a_fineLayout[fineLit];
                coveredBox.coarsen(a_refRatio);
                coveredBox &= maskBox;

                if (!coveredBox.isEmpty()) {
                    maskFAB.setVal(1, coveredBox, 0);
                }
            } // fineLit
        } // dit

        // 4. At each bdry face, copy the valid data.
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox&    crseFAB = a_crse[dit][fcDir];
            FArrayBox&    exFAB   = exData[dit];
            BaseFab<int>& maskFAB = ccMask[dit];
            const Box&    ccValid = grids[dit];

            for (SideIterator sit; sit.ok(); ++sit) {
                const Side::LoHiSide& side  = sit();
                const int             isign = sign(side);
                const Box ghostBox = adjCellBox(ccValid, fcDir, side, 1);

                // Make these overlie the ghost.
                crseFAB.shiftHalf(fcDir, isign);

                FORT_CFINTERP_COPYMASKEDDATATOBDRYFACES(
                    CHF_FRA(crseFAB),
                    CHF_CONST_FRA(exFAB),
                    CHF_CONST_FIA1(maskFAB, 0),
                    CHF_BOX(ghostBox),
                    CHF_CONST_INT(fcDir),
                    CHF_CONST_INT(isign));

                // Restore centering.
                crseFAB.shiftHalf(fcDir, -isign);
            } // sit
        } // dit
    } // fcDir

    debugCheckValidFaceOverlap(a_crse);
}


// -----------------------------------------------------------------------------
void
CFInterp::localizeCrseData(LevelData<FArrayBox>&       a_localCrseData,
                           const LevelData<FArrayBox>& a_userCrseData) const
{
    // Create references.
    LevelData<FArrayBox>&       dest       = a_localCrseData;
    const LevelData<FArrayBox>& src        = a_userCrseData;
    const BoxLayout&            destLayout = dest.boxLayout();
    const BoxLayout&            srcLayout  = src.boxLayout();
    const ProblemDomain&        domain     = src.getBoxes().physDomain();

    // Copy valid data, including one ghost layer.
    {
        CH_assert(dest.ghostVect() >= IntVect::Unit);
        CH_assert(src.ghostVect() >= IntVect::Unit);

        // TODO: It would be better to copy all ghosts wherever the data is
        // available rather than extrapolate it.

        // Should we cache this?
        Copier copier;
        LayoutTools::defineImpartialCopier(copier,
                                           srcLayout,
                                           destLayout,
                                           domain,
                                           IntVect::Unit,  // srcGhostVect
                                           IntVect::Unit); // destGhostVect

        debugInitLevel(dest);
        src.copyTo(dest, copier);
    }

    // Fill all ghosts beyond the first layer.
    BCTools::extrapAllGhosts(dest, 2, IntVect::Unit);
    dest.exchange(m_crseExCopier);
    dest.exchange(m_crseExCornerCopier);
}


// -----------------------------------------------------------------------------
void
CFInterp::localizeCrseData(LevelData<FluxBox>&       a_localCrseData,
                           const LevelData<FluxBox>& a_userCrseData) const
{
    // Create references
    LevelData<FluxBox>&      dest          = a_localCrseData;
    // const DisjointBoxLayout& destGrids     = dest.getBoxes();
    const BoxLayout&         destLayout    = dest.boxLayout();
    // const IntVect&           destGhostVect = dest.ghostVect();
    // const Interval&          destIvl       = dest.interval();

    const LevelData<FluxBox>& src          = a_userCrseData;
    const DisjointBoxLayout&  srcGrids     = src.getBoxes();
    const BoxLayout&          srcLayout    = src.boxLayout();
    // const IntVect&            srcGhostVect = src.ghostVect();
    // const Interval&           srcIvl       = src.interval();

    const ProblemDomain& domain = srcGrids.physDomain();

    // Copy valid data.
    // Copier creation is not a bottleneck, so don't bother creating a cache.
    {
        CH_assert(src.ghostVect() >= IntVect::Unit);
        CH_assert(dest.ghostVect() >= IntVect::Unit);

        std::array<StaggeredCopier, CH_SPACEDIM> copier;
        LayoutTools::defineImpartialCopier(copier,
                                           srcLayout,
                                           destLayout,
                                           domain,
                                           IntVect::Unit,  // srcGhostVect
                                           IntVect::Unit); // destGhostVect

        debugCheckValidFaceOverlap(src);
        debugInitLevel(dest);
        src.copyTo(dest, copier);
    }

    // // I wrote this version to reduce the MPI_Wait time. Needs testing.
    // {
    //     std::array<StaggeredCopier, SpaceDim> copier;
    //     std::array<LevelData<FArrayBox>, SpaceDim> srcAlias;
    //     std::array<LevelData<FArrayBox>, SpaceDim> destAlias;

    //     // Should we cache this?
    //     LayoutTools::defineImpartialCopier(copier,
    //                                        srcLayout,
    //                                        destLayout,
    //                                        domain,
    //                                        IntVect::Unit,   // srcGhostVect
    //                                        IntVect::Unit);  // destGhostVect

    //     for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
    //         FABAliasFlBxDataFactory destFactory(&dest, destIvl, fcDir);
    //         destAlias[fcDir].define(destGrids, dest.nComp(), destGhostVect, destFactory);

    //         LevelData<FluxBox>* srcPtr = const_cast<LevelData<FluxBox>*>(&src);
    //         FABAliasFlBxDataFactory srcFactory(srcPtr, srcIvl, fcDir);
    //         srcAlias[fcDir].define(srcGrids, src.nComp(), srcGhostVect, srcFactory);

    //         srcAlias[fcDir].copyToBegin(srcIvl, destAlias[fcDir], destIvl, copier[fcDir]);
    //     }

    //     for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
    //         srcAlias[fcDir].copyToEnd(destAlias[fcDir], destIvl, copier[fcDir]);
    //     }
    // }

    // Fill all ghosts beyond the first layer.
    BCTools::extrapAllGhosts(dest, 2, IntVect::Unit);
    dest.exchange();  // Does corners too!

    debugCheckValidFaceOverlap(dest);
}


// -----------------------------------------------------------------------------
void
CFInterp::interpAlongParentEdgesIn2D(FArrayBox&       a_fineFAB,
                                     const FArrayBox& a_crseFAB,
                                     const Box&       a_crseBox,
                                     const int        a_interpDir,
                                     const IntVect&   a_refRatio)
{
    // Sanity checks
    CH_assert(a_interpDir >= 0);
    CH_assert(a_interpDir < SpaceDim);

    CH_assert(a_refRatio >= IntVect::Unit);

    CH_assert(a_crseBox.type(a_interpDir) == IndexType::CELL);
    CH_assert(a_fineFAB.box().type() == a_crseBox.type());
    CH_assert(a_crseFAB.box().type() == a_crseBox.type());
    CH_assert(a_crseFAB.box().contains(grow(a_crseBox, BASISV(a_interpDir))));

    CH_assert(a_crseFAB.nComp() == a_fineFAB.nComp());

    const int     ncomp = a_crseFAB.nComp();
    const Real    r     = Real(a_refRatio[a_interpDir]);
    const IntVect ei    = BASISV(a_interpDir);
    IntVect       fineIV;

    const Box refBox(
        IntVect::Zero, (a_refRatio - IntVect::Unit) * ei, a_crseBox.type());
    BoxIterator refBit(refBox);
    BoxIterator crseBit(a_crseBox);

    // Loop over the FAB comps.
    for (int comp = 0; comp < ncomp; ++comp) {
        // Loop over the coarse faces.
        for (crseBit.reset(); crseBit.ok(); ++crseBit) {
            const IntVect& crseIV = crseBit();

            // Get coarse values
            const Real uc = a_crseFAB(crseIV     , comp);
            const Real uw = a_crseFAB(crseIV - ei, comp);
            const Real ue = a_crseFAB(crseIV + ei, comp);

            // Compute the interpolation coefficients.
            const Real d = 0.5 * (ue - 2. * uc + uw);
            const Real b = 0.5 * (ue - uw);
            const Real a = uc - d * (r - 1.0) * (r + 1.0) / (12.0 * r * r);

            for (refBit.reset(); refBit.ok(); ++refBit) {
                const IntVect& refIV = refBit();
                fineIV               = a_refRatio * crseIV + refIV;
                const Real i         = Real(refBit()[a_interpDir]);

                // Quadratic order (fixes uc, 1st deriv, and 2nd deriv)
                // This is the 3D stencil with y = 0.
                const Real x = (2.0 * i + 1. - r) / (2.0 * r);
                a_fineFAB(fineIV, comp) = a + x*(b + d*x);

            }  // end loop over fine faces within parent face (refBit, refIV)
        } // end loop over coarse (parent) face (crseBit, crseIV)
    } // end loop over FAB comps (comp)
}


// -----------------------------------------------------------------------------
void
CFInterp::interpAlongParentFacesIn3D(FArrayBox&       a_fineFAB,
                                     const FArrayBox& a_crseFAB,
                                     const Box&       a_crseBox,
                                     const int        a_fcDir,
                                     const IntVect&   a_refRatio)
{
    // Sanity checks
    CH_assert(a_fcDir >= 0);
    CH_assert(a_fcDir < SpaceDim);

    CH_assert(a_refRatio >= IntVect::Unit);

    CH_assert(a_crseBox.type() == BASISV(a_fcDir));
    CH_assert(a_fineFAB.box().type() == a_crseBox.type());
    CH_assert(a_crseFAB.box().type() == a_crseBox.type());
    CH_assert(a_crseFAB.box().contains(
        grow(a_crseBox, IntVect::Unit - BASISV(a_fcDir))));

    CH_assert(a_crseFAB.nComp() == a_fineFAB.nComp());

    const int     d1 = (a_fcDir + 1) % SpaceDim;  // interp dir 1
    const IntVect e1 = BASISV(d1);                // unit vector along d1
    const Real    r1 = Real(a_refRatio[d1]);      // refRatio along d1

    const int     d2 = (a_fcDir + 2) % SpaceDim;  // ditto in
    const IntVect e2 = BASISV(d2);                // interp dir 2
    const Real    r2 = Real(a_refRatio[d2]);

    IntVect fineIV;
    const int ncomp = a_crseFAB.nComp();

    const Box refBox(IntVect::Zero,
                     (a_refRatio - IntVect::Unit) * (e1 + e2),
                     a_crseBox.type());
    BoxIterator refBit(refBox);
    BoxIterator crseBit(a_crseBox);

    // Loop over the FAB comps.
    for (int comp = 0; comp < ncomp; ++comp) {
        // Loop over the coarse faces.
        for (crseBit.reset(); crseBit.ok(); ++crseBit) {
            const IntVect& crseIV = crseBit();

            // Get coarse values
            const Real uc  = a_crseFAB(crseIV, comp);
            const Real ue  = a_crseFAB(crseIV + e1     , comp);
            const Real une = a_crseFAB(crseIV + e1 + e2, comp);
            const Real un  = a_crseFAB(crseIV      + e2, comp);
            const Real unw = a_crseFAB(crseIV - e1 + e2, comp);
            const Real uw  = a_crseFAB(crseIV - e1     , comp);
            const Real usw = a_crseFAB(crseIV - e1 - e2, comp);
            const Real us  = a_crseFAB(crseIV      - e2, comp);
            const Real use = a_crseFAB(crseIV + e1 - e2, comp);


            // Compute the interpolation coefficients.
            const Real f = 0.5 * (un - 2. * uc + us);
            const Real e = 0.25 * (une - use - unw + usw);
            const Real d = 0.5 * (ue - 2. * uc + uw);
            const Real c = 0.5 * (un - us);
            const Real b = 0.5 * (ue - uw);
            const Real a = uc - d * (r1 - 1.0) * (r1 + 1.0) / (12.0 * r1 * r1)
                              - f * (r2 - 1.0) * (r2 + 1.0) / (12.0 * r2 * r2);

            for (refBit.reset(); refBit.ok(); ++refBit) {
                const IntVect& refIV = refBit();
                fineIV               = a_refRatio * crseIV + refIV;
                const Real i         = Real(refBit()[d1]);
                const Real j         = Real(refBit()[d2]);

                // Quadratic order (fixes uc, two 1st derivs, and three 2nd derivs)
                const Real x = (2.0 * i + 1. - r1) / (2.0 * r1);
                const Real y = (2.0 * j + 1. - r2) / (2.0 * r2);
                a_fineFAB(fineIV, comp) = a + x*(b + d*x + e*y) + y*(c + f*y);

            }  // end loop over fine faces within parent face (refBit, refIV)
        } // end loop over coarse (parent) face (crseBit, crseIV)
    } // end loop over FAB comps (comp)
}


#ifdef ALLOW_DIVFREEINTERP
// -----------------------------------------------------------------------------
// * \brief      Divergence-free CF interpolator.
// * \see        divFreeInterp
// * \details    a_crse's grids must be a coarsened version of a_fine's grids.
// *             a_fineParentPtr can be defined on any set of grids.
// -----------------------------------------------------------------------------
void
CFInterp::localDivFreeInterp(LevelData<FluxBox>&       a_fine,
                             const LevelData<FluxBox>& a_crse,
                             const LevelData<FluxBox>* a_fineDataPtr,
                             const Copier*             a_fineDataCopierPtr) const
{
    // UNDEFINED_FUNCTION();
    CH_assert(m_isDefined);
    CH_assert(a_fine.getBoxes() == m_grids);
    CH_assert(a_crse.getBoxes() == m_crseGrids);
    CH_assert(a_crse.nComp() == a_fine.nComp());

    DataIterator dit = m_grids.dataIterator();
    Copier parentCopier;

    // Interp to parent faces.
    for (dit.reset(); dit.ok(); ++dit) {
        this->localDivFreeInterp_Step1(
            a_fine[dit], a_crse[dit], m_crseGrids[dit]);
    }

    // Copy whatever valid data the user wishes to keep around.
    if (a_fineDataPtr) {
        debugCheckValidFaceOverlap(*a_fineDataPtr);

        // Copy.
        if (a_fineDataPtr) {
            a_fineDataPtr->copyTo(a_fine);
            // if (a_fineDataCopierPtr) {
            //     a_fineDataPtr->copyTo(a_fine, *a_fineDataCopierPtr);
            // } else {
            //     parentCopier.define(
            //         a_fineDataPtr->getBoxes(), m_grids, m_domain);
            //     a_fineDataPtr->copyTo(a_fine, parentCopier);
            // }
        }

        // Repair the vertical interpolation of w.
        if (SpaceDim > 2) {
            for (dit.reset(); dit.ok(); ++dit) {
                Box crseCCValid = m_grids[dit];
                crseCCValid.coarsen(m_refRatio);
                this->vertVelInterp(
                    a_fine[dit][2], a_crse[dit][2], crseCCValid);
            }
        }

        a_fine.exchange(); // TODO: Is this needed? If so, use a copier.
    }

    // Interp to parent interiors.
    for (dit.reset(); dit.ok(); ++dit) {
        this->localDivFreeInterp_Step2(
            a_fine[dit], a_crse[dit], m_crseGrids[dit]);
    }

    // Repair the interiors that were overwritten in step 2.
    if (a_fineDataPtr) {
        a_fineDataPtr->copyTo(a_fine);
        // if (a_fineDataCopierPtr) {
        //     a_fineDataPtr->copyTo(a_fine, *a_fineDataCopierPtr);
        // } else {
        //     parentCopier.define(a_fineDataPtr->getBoxes(), m_grids, m_domain);
        //     a_fineDataPtr->copyTo(a_fine, parentCopier);
        // }
    }

    debugCheckValidFaceOverlap(a_fine);
}


// -----------------------------------------------------------------------------
void
CFInterp::localDivFreeInterp_Step1(FluxBox&       a_fineAdvVelFB,
                                   const FluxBox& a_crseAdvVelFB,
                                   const Box&     a_ccCrseInterpBox) const
{
    CH_assert(m_isDefined);
    CH_assert(a_fineAdvVelFB.nComp() == 1);
    CH_assert(a_crseAdvVelFB.nComp() == 1);
    CH_assert(a_crseAdvVelFB.box().contains(a_ccCrseInterpBox));

    debugInit(a_fineAdvVelFB);

    if (SpaceDim == 2) {
        // Interpolate at parent faces.
        for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
            this->interpAlongParentEdgesIn2D(
                a_fineAdvVelFB[fcDir],
                a_crseAdvVelFB[fcDir],
                surroundingNodes(a_ccCrseInterpBox, fcDir),
                1 - fcDir,  // interpDir
                m_refRatio);
        }
    } else {
        // TODO: Why refine then coarsen?
        const Box fineCCValid = ::refine(a_ccCrseInterpBox, m_refRatio);

        // Interpolate u along parent faces.
        Box crseBox = surroundingNodes(fineCCValid, 0);
        crseBox.coarsen(m_refRatio);
        this->interpAlongParentFacesIn3D(
            a_fineAdvVelFB[0], a_crseAdvVelFB[0], crseBox, 0, m_refRatio);

        // Interpolate v along parent faces.
        crseBox = surroundingNodes(fineCCValid, 1);
        crseBox.coarsen(m_refRatio);
        this->interpAlongParentFacesIn3D(
            a_fineAdvVelFB[1], a_crseAdvVelFB[1], crseBox, 1, m_refRatio);

        // Interpolate w to all faces, not just parents.
        this->vertVelInterp(
            a_fineAdvVelFB[2], a_crseAdvVelFB[2], a_ccCrseInterpBox);
    }
}


// -----------------------------------------------------------------------------
void
CFInterp::localDivFreeInterp_Step2(FluxBox&       a_fineAdvVelFB,
                                   const FluxBox& a_crseAdvVelFB,
                                   const Box&     a_ccCrseInterpBox/*,
                                   const FluxBox& a_crseJgupFB*/) const
{
    CH_assert(m_isDefined);
    CH_assert(a_fineAdvVelFB.nComp() == 1);
    CH_assert(a_crseAdvVelFB.nComp() == 1);
    CH_assert(a_crseAdvVelFB.box().contains(a_ccCrseInterpBox));

    // Gather references
    const int     rx    = m_refRatio[0];
    const int     ry    = m_refRatio[1];
    const int     rz    = m_refRatio[SpaceDim - 1];
    const int     rank  = 2 * (rx * ry) - (rx + ry);
    const IntVect ex    = BASISV(0);
    const IntVect ey    = BASISV(1);
    const IntVect ez    = BASISV(SpaceDim - 1);
    const Real    invDx = 1.0 / m_fineDXi[0];
    const Real    invDy = 1.0 / m_fineDXi[1];
    const Real    invDz = 1.0 / m_fineDXi[SpaceDim - 1];
    const IntVect horizRefRatio(D_DECL(rx, ry, 1));
    const IntVect vertRefRatio(D_DECL(1, 1, rz));

    Box       crseCCValid = a_ccCrseInterpBox;
    const Box fineCCValid = ::refine(crseCCValid, m_refRatio);

    // No horiz refinement means no furter interp needed.
    if (rx == 1 && ry == 1) return;

    // Interpolate crse vel horizontal divergence to fine CCs.
    FArrayBox fineDivFAB;
    if (SpaceDim == 2) {
        // Assume crse vel is projected.
        fineDivFAB.define(fineCCValid, 1);
        fineDivFAB.setVal(0.0);
    } else {
        if (rx > 1 || ry > 1) {
            // Compute horiz div = -dw/dz.
            fineDivFAB.define(fineCCValid, 1);
            Real wt, wb;
            for (BoxIterator bit(fineCCValid); bit.ok(); ++bit) {
                const IntVect& cc = bit();
                wt = a_fineAdvVelFB[2](cc + ez);
                wb = a_fineAdvVelFB[2](cc);
                fineDivFAB(cc) = -(wt - wb) * invDz;
            }
        }
    }

    // Interpolate Curl to fine nodes.
    FArrayBox fineCurlFAB;
    if (rx > 1 && ry > 1) {
        // // Compute the covariant velocity, u_i = (Ju^i) / (Jg^{ii}).
        // Box crseOneGhost = crseCCValid;
        // crseOneGhost.grow(0, 1);
        // crseOneGhost.grow(1, 1);

        // FluxBox crseCovelFB(crseOneGhost, 1);
        // debugInit(crseCovelFB);
        // for (int d = 0; d < SpaceDim - 1; ++d) {
        //     Box thisCompBox = crseCCValid;
        //     thisCompBox.surroundingNodes(d);
        //     thisCompBox.grow(1 - d, 1);

        //     CH_assert(a_crseAdvVelFB[d].box().contains(thisCompBox));
        //     CH_assert(a_crseJgupFB[d].box().contains(thisCompBox));

        //     crseCovelFB[d].copy(a_crseAdvVelFB[d], thisCompBox);
        //     crseCovelFB[d].divide(a_crseJgupFB[d], thisCompBox, 0, 0, 1);
        // }

        // Compute the covariant velocity, u_i = (Ju^i) / (Jg^{ii}).
        Box crseOneGhost = crseCCValid;
        crseOneGhost.grow(0, 1);
        crseOneGhost.grow(1, 1);

        FluxBox crseCovelFB(a_crseAdvVelFB.box(), 1);
        debugInit(crseCovelFB);
        crseCovelFB.copy(a_crseAdvVelFB);
        // crseCovelFB.divide(a_crseJgupFB, a_crseAdvVelFB.box(), 0, 0, 1);

        if (SpaceDim == 2) {
            // Compute Curl at coarse nodes (2D).
            Box crseCurlBox = crseCCValid;
            crseCurlBox.surroundingNodes(0);
            crseCurlBox.surroundingNodes(1);

            FArrayBox crseCurlFAB(crseCurlBox, 1);
            this->horizNodeCurl(
                crseCurlFAB, crseCurlBox, crseCovelFB, m_crseDXi);

            // Interpolate Curl to the fine nodes (2D).
            Box fineCurlBox = fineCCValid;
            fineCurlBox.surroundingNodes(0);
            fineCurlBox.surroundingNodes(1);

            fineCurlFAB.define(fineCurlBox, 1);
            this->horizNodeInterp(
                fineCurlFAB, crseCurlFAB, crseCurlBox, horizRefRatio);

        } else {
            // Compute Curl at coarse edges (3D).
            Box crseCurlBox = crseCCValid;
            crseCurlBox.surroundingNodes(0);
            crseCurlBox.surroundingNodes(1);
            crseCurlBox.grow(2, 1);

            FArrayBox crseCurlFAB(crseCurlBox, 1);
            this->horizNodeCurl(
                crseCurlFAB, crseCurlBox, crseCovelFB, m_crseDXi);

            // Interpolate horizontally to fine edges (3D).
            Box semiCurlBox = fineCCValid;
            semiCurlBox.coarsen(vertRefRatio);
            semiCurlBox.surroundingNodes(0);
            semiCurlBox.surroundingNodes(1);
            semiCurlBox.grow(2, 1);

            FArrayBox semiCurlFAB(semiCurlBox, 1);
            this->horizNodeInterp(
                semiCurlFAB, crseCurlFAB, crseCurlBox, horizRefRatio);

            // Interpolate vertically.
            Box fineCurlBox = fineCCValid;
            fineCurlBox.surroundingNodes(0);
            fineCurlBox.surroundingNodes(1);

            fineCurlFAB.define(fineCurlBox, 1);

            TODONOTE("Use conservative vertical interpolation for curl?");
            Box refBox(IntVect::Zero, IntVect(D_DECL(0, 0, rz - 1)));
            IntVect civ, fiv;
            Real ct, cm, cb, ctm, cmb, z;

            semiCurlBox.grow(2, -1);
            for (BoxIterator cbit(semiCurlBox); cbit.ok(); ++cbit) {
                const IntVect& civ = cbit();
                for (BoxIterator rbit(refBox); rbit.ok(); ++rbit) {
                    const IntVect& riv = rbit();

                    fiv[0] = civ[0];
                    fiv[1] = civ[1];
                    fiv[2] = civ[2] * rz + riv[2];

                    ct  = semiCurlFAB(civ + ez);
                    cm  = semiCurlFAB(civ);
                    cb  = semiCurlFAB(civ - ez);

                    ctm = 0.5 * (ct + cm);
                    cmb = 0.5 * (cm + cb);

                    z = (Real(riv[2]) + 0.5) / Real(rz);
                    fineCurlFAB(fiv) = z * ctm + (1.0 - z) * cmb;
                } // rbit
            } // cbit
        }
    }  // end computing fineCurlFAB

    // Loop over each coarse-level cell and interpolate the velocity to
    // non-parent faces.
    crseCCValid.refine(vertRefRatio); // No-op in 2D.
    BoxIterator cbit(crseCCValid);

    if (s_mode == 0) {
        // Original interp method
        Vector<Real> rhs(rank, 0.0);

        for (cbit.reset(); cbit.ok(); ++cbit) {
            const IntVect& crseIV = cbit();
            const Box      parentCellBox(crseIV, crseIV);
            const Box      fineCellBox = ::refine(parentCellBox, horizRefRatio);
            int            row, col;

            // cols <--> [div, curl]
            // rows <--> [u,v]

            // Setup rhs vector - Div comps
            {
                const IntVect& smallEnd = fineCellBox.smallEnd();
                BoxIterator    rbit(m_colIdx[0].box());
                IntVect        fineIV;

                for (rbit.reset(); rbit.ok(); ++rbit) {
                    const IntVect& refIV = rbit();  // smallEnd = (0,0)
                    fineIV = smallEnd + refIV;      // smallEnd is correct.
                    col    = m_colIdx[0](refIV);

                    // Skip linearly dependent eq.
                    if (col < 0) continue;

                    // Set div value
                    rhs[col] = fineDivFAB(fineIV);

                    // Add lhs terms
                    if (refIV[0] == 0) {
                        rhs[col] += a_fineAdvVelFB[0](fineIV) * invDx;
                    }

                    if (refIV[0] == rx - 1) {
                        rhs[col] -= a_fineAdvVelFB[0](fineIV + ex) * invDx;
                    }

                    if (refIV[1] == 0) {
                        rhs[col] += a_fineAdvVelFB[1](fineIV) * invDy;
                    }

                    if (refIV[1] == ry - 1) {
                        rhs[col] -= a_fineAdvVelFB[1](fineIV + ey) * invDy;
                    }
                }
            }

            // Setup rhs vector - Curl comps
            {
                Box fineCurlBox = fineCellBox;
                fineCurlBox.surroundingNodes(0);
                fineCurlBox.surroundingNodes(1);
                const IntVect& smallEnd = fineCurlBox.smallEnd();

                BoxIterator bit(m_colIdx[1].box());
                IntVect     fineIV;

                for (bit.reset(); bit.ok(); ++bit) {
                    const IntVect& refIV = bit();
                    fineIV = smallEnd + refIV;
                    col    = m_colIdx[1](refIV);

                    // Set curl value
                    rhs[col] = fineCurlFAB(fineIV);
                }
            }

            // Compute unknowns as Minv*rhs -- uFAB
            for (int fcDir = 0; fcDir <= 1; ++fcDir) {
                const IntVect smallEnd =
                    surroundingNodes(fineCellBox, fcDir).smallEnd();

                BoxIterator bit(m_rowIdx[fcDir].box());
                IntVect     mIV = IntVect::Zero;
                IntVect     fineIV;

                for (bit.reset(); bit.ok(); ++bit) {
                    const IntVect& refIV = bit();
                    fineIV = smallEnd + refIV;
                    row    = m_rowIdx[fcDir](refIV);
                    mIV[0] = row;

                    a_fineAdvVelFB[fcDir](fineIV) = 0.0;
                    for (int col = 0; col < rank; ++col) {
                        mIV[1] = col;
                        CH_assert(a_fineAdvVelFB[fcDir].box().contains(fineIV));
                        CH_assert(m_Minv.box().contains(mIV));
                        a_fineAdvVelFB[fcDir](fineIV) += m_Minv(mIV) * rhs[col];
                    }
                }
            }
        }  // end loop over parent cells (cbit, parentCellBox)

    } else if (s_mode == 1) {
        Vector<Real> rhs(rank + 1, 0.0);
        Vector<Real> DCt_rhs(rank, 0.0);

        for (cbit.reset(); cbit.ok(); ++cbit) {
            const IntVect& crseIV = cbit();
            const Box      parentCellBox(crseIV, crseIV);
            const Box      fineCellBox = ::refine(parentCellBox, horizRefRatio);
            int            row, col;

            // cols <--> [div, curl]
            // rows <--> [u,v]

            // Setup rhs vector - Div comps
            {
                const IntVect& smallEnd = fineCellBox.smallEnd();
                BoxIterator    rbit(m_colIdx[0].box());
                IntVect        fineIV;

                for (rbit.reset(); rbit.ok(); ++rbit) {
                    const IntVect& refIV = rbit();  // smallEnd = (0,0)
                    fineIV = smallEnd + refIV;      // smallEnd is correct.
                    col    = m_colIdx[0](refIV);

                    // Set div value
                    rhs[col] = fineDivFAB(fineIV);

                    // Add lhs terms
                    if (refIV[0] == 0) {
                        rhs[col] += a_fineAdvVelFB[0](fineIV) * invDx;
                    }

                    if (refIV[0] == rx - 1) {
                        rhs[col] -= a_fineAdvVelFB[0](fineIV + ex) * invDx;
                    }

                    if (refIV[1] == 0) {
                        rhs[col] += a_fineAdvVelFB[1](fineIV) * invDy;
                    }

                    if (refIV[1] == ry - 1) {
                        rhs[col] -= a_fineAdvVelFB[1](fineIV + ey) * invDy;
                    }
                }
            }

            // Setup rhs vector - Curl comps
            {
                Box fineCurlBox = fineCellBox;
                fineCurlBox.surroundingNodes(0);
                fineCurlBox.surroundingNodes(1);
                const IntVect& smallEnd = fineCurlBox.smallEnd();

                BoxIterator bit(m_colIdx[1].box());
                IntVect     fineIV;

                for (bit.reset(); bit.ok(); ++bit) {
                    const IntVect& refIV = bit();  // smallEnd = (0,0)
                    fineIV = smallEnd + refIV;     // smallEnd is correct.
                    col    = m_colIdx[1](refIV);

                    // Set curl value
                    rhs[col] = fineCurlFAB(fineIV);
                }
            }

            // Compute DCt * rhs.
            IntVect iv(D_DECL(0,0,0));
            for (int i = 0; i < rank; ++i) {
                DCt_rhs[i] = 0.0;
                iv[0] = i;
                for (int j = 0; j < rank + 1; ++j) {
                    iv[1] = j;
                    DCt_rhs[i] += m_DCt(iv) * rhs[j];
                }
            }

            // Compute unknowns as Minv*rhs -- uFAB
            for (int fcDir = 0; fcDir <= 1; ++fcDir) {
                const IntVect smallEnd =
                    surroundingNodes(fineCellBox, fcDir).smallEnd();

                BoxIterator bit(m_rowIdx[fcDir].box());
                IntVect     mIV = IntVect::Zero;
                IntVect     fineIV;

                for (bit.reset(); bit.ok(); ++bit) {
                    const IntVect& refIV = bit();  // .smallEnd = (1,0)
                    fineIV = smallEnd + refIV;     // .smallEnd is correct.
                    row    = m_rowIdx[fcDir](refIV);
                    mIV[0] = row;

                    a_fineAdvVelFB[fcDir](fineIV) = 0.0;
                    for (int col = 0; col < rank; ++col) {
                        mIV[1] = col;

                        CH_assert(a_fineAdvVelFB[fcDir].box().contains(fineIV));
                        CH_assert(m_Minv.box().contains(mIV));

                        a_fineAdvVelFB[fcDir](fineIV) +=
                            m_Minv(mIV) * DCt_rhs[col];
                    }
                }
            }
        }  // end loop over parent cells (cbit, parentCellBox)

    } else {
        MAYDAYERROR("CFInterp does not recognize s_mode = " << s_mode);
    }
}


// -----------------------------------------------------------------------------
// \brief      Computes the curl of velocity in xy-plane.
//
// \param[out] a_curlFAB   NC in 2D, EC in 3D.
// \param[in]  a_curlBox   Where do you want the curl? Same NC/EC in 2D/3D.
// \param[in]  a_velFB     The FC velocity.
// \param[in]  a_dXi       The grid spacing.
// -----------------------------------------------------------------------------
void
CFInterp::horizNodeCurl(FArrayBox&      a_curlFAB,
                        const Box&      a_curlBox,
                        const FluxBox&  a_vel,
                        const RealVect& a_dXi) const
{
    const IntVect ex    = BASISV(0);
    const IntVect ey    = BASISV(1);
    const Real    invDx = 1.0 / a_dXi[0];
    const Real    invDy = 1.0 / a_dXi[1];
    BoxIterator   ncBit(a_curlBox);

#ifndef NDEBUG
    CH_assert(a_curlBox.type(0) == IndexType::NODE);
    CH_assert(a_curlBox.type(1) == IndexType::NODE);
    CH_assert(a_curlFAB.box().type() == a_curlBox.type());
    CH_assert(a_curlFAB.box().contains(a_curlBox));
    CH_assert(a_curlFAB.nComp() == 1);
    CH_assert(a_vel.nComp() == 1);

    {
        Box srcBox = a_curlBox;
        srcBox.grow(ey);
        srcBox.enclosedCells(1);
        CH_assert(a_vel[0].box().contains(srcBox));

        srcBox = a_curlBox;
        srcBox.grow(ex);
        srcBox.enclosedCells(0);
        CH_assert(a_vel[1].box().contains(srcBox));
    }

    debugInit(a_curlFAB);
#endif

    if (a_vel.nComp() == 1) {
        for (ncBit.reset(); ncBit.ok(); ++ncBit) {
            const IntVect& nc = ncBit();

            a_curlFAB(nc) = (a_vel[1](nc) - a_vel[1](nc - ex)) * invDx -
                            (a_vel[0](nc) - a_vel[0](nc - ey)) * invDy;
        }
    } else {
        CH_verify(false); // Is this used?
        for (ncBit.reset(); ncBit.ok(); ++ncBit) {
            const IntVect& nc = ncBit();

            a_curlFAB(nc) = (a_vel[1](nc, 1) - a_vel[1](nc - ex, 1)) * invDx -
                            (a_vel[0](nc, 0) - a_vel[0](nc - ey, 0)) * invDy;
        }
    }
}


// -----------------------------------------------------------------------------
// Bi-linear interpolation used on horiz curl.
// In 3D, this does NOT interpolate in the z dir, it only fills a_fineFAB
// at the bottom z-centered edge of each parent cell. It is up to you to
// interpolate to the remainiing z-centered edges.
// -----------------------------------------------------------------------------
void
CFInterp::horizNodeInterp(FArrayBox&       a_fineFAB,
                          const FArrayBox& a_crseFAB,
                          Box              a_crseBox,
                          const IntVect&   a_refRatio) const
{
    CH_assert(a_fineFAB.nComp() == 1);
    CH_assert(a_crseFAB.nComp() == 1);

    debugInit(a_fineFAB);

    const IntVect ex = BASISV(0);
    const IntVect ey = BASISV(1);

    IntVect href = a_refRatio;
    if (SpaceDim > 2) href[2] = 0;

    a_crseBox.growHi(0, -1);  // We only want cbit to loop over the
    a_crseBox.growHi(1, -1);  // lower node of each parent cell.
    BoxIterator cbit(a_crseBox);

    Box         refBox(IntVect::Zero, href, ex + ey);
    BoxIterator rbit(refBox);

    IntVect fb;
    Real    x0, y0, x1, y1;
    Real    c00, c01, c10, c11;

    for (cbit.reset(); cbit.ok(); ++cbit) {
        const IntVect& cb = cbit();

        // Get coarse-level curl values at corner nodes.
        c00 = a_crseFAB(cb);
        c01 = a_crseFAB(cb + ey);
        c10 = a_crseFAB(cb + ex);
        c11 = a_crseFAB(cb + ex + ey);

        for (rbit.reset(); rbit.ok(); ++rbit) {
            // refBox indices
            const IntVect& rb = rbit();

            // Fine-level node indices
            D_TERM(
            fb[0] = a_refRatio[0] * cb[0] + rb[0];,
            fb[1] = a_refRatio[1] * cb[1] + rb[1];,
            fb[2] = a_refRatio[2] * cb[2] + rb[2];)

            // Right/top value multiplier
            x1 = Real(rb[0]) / Real(m_refRatio[0]);
            y1 = Real(rb[1]) / Real(m_refRatio[1]);
            // Left/bottom value multiplier
            x0 = 1.0 - x1;
            y0 = 1.0 - y1;

            // Interpolate
            a_fineFAB(fb) =
                x0 * y0 * c00 + x1 * y0 * c10 + x0 * y1 * c01 + x1 * y1 * c11;
        }
    }
}


// -----------------------------------------------------------------------------
//  \brief      Vertical, explicit, cubic interpolation for 3D only.
//              FABs must be FC in vertical dir.
//  \details
//   Used to interpolate the vertical velocity.
//   This will call interpAlongParentFacesIn3D first, then interpolate
//   to faces between.
// -----------------------------------------------------------------------------
void
CFInterp::vertVelInterp(FArrayBox&       a_fineFAB,
                        const FArrayBox& a_crseFAB,
                        const Box&       a_ccCrseInterpBox) const
{
    BUG("Try using the convervative interp because the FABs are scaled by J.");

    CH_assert(SpaceDim == 3);
    CH_assert(m_isDefined);
    CH_assert(a_fineFAB.nComp() == 1);
    CH_assert(a_crseFAB.nComp() == 1);

    // Gather references
    const int     rx = m_refRatio[0];
    const int     ry = m_refRatio[1];
    const int     rz = m_refRatio[2];
    const IntVect ez = BASISV(2);

    const Box& crseCCValid = a_ccCrseInterpBox;
    const Box& fineCCValid = ::refine(crseCCValid, m_refRatio);

    // Interpolate w horizontally along parent faces.
    // Include rz fine ghosts (= one crse ghost) at top and bottom.
    FArrayBox localFineWFAB;
    {
        Box fineBox = surroundingNodes(fineCCValid, 2);
        fineBox.grow(2, rz);
        localFineWFAB.define(fineBox, 1);

        Box crseBox = surroundingNodes(fineCCValid, 2);
        crseBox.coarsen(m_refRatio);
        crseBox.grow(2, 1);
        this->interpAlongParentFacesIn3D(
            localFineWFAB, a_crseFAB, crseBox, 2, m_refRatio);
    }

    // Interpolate vertically to all fine faces.
    Real    a, b, c, d;
    Real    wtt, wt, wb, wbb, z;
    IntVect basefiv, fiv;

    Box crseBox = surroundingNodes(fineCCValid, 2);
    crseBox.coarsen(m_refRatio);
    crseBox.growHi(2, -1);

    Box refBox(IntVect::Zero, IntVect(D_DECL(rx - 1, ry - 1, rz)));

    for (BoxIterator cbit(crseBox); cbit.ok(); ++cbit) {
        const IntVect& civ = cbit();

        for (BoxIterator rbit(refBox); rbit.ok(); ++rbit) {
            const IntVect& riv = rbit();
            z = Real(riv[2]) / Real(rz);
            fiv = m_refRatio * civ + riv;

            basefiv = m_refRatio * civ;
            basefiv[0] += riv[0];
            basefiv[1] += riv[1];

            wtt = localFineWFAB(basefiv + 2 * rz * ez);
            wt  = localFineWFAB(basefiv + rz * ez);
            wb  = localFineWFAB(basefiv);
            wbb = localFineWFAB(basefiv - rz * ez);

            // // Linear interp (O(dz^2) error)
            // a_fineFAB(fiv) = (1.0 - z) * wb + z * wt;

            // Cubic interp (Gives a smoother result in all vel comps.)
            a = (         6.*wb               ) / 6.;
            b = (-2.*wbb -3.*wb +6.*wt -1.*wtt) / 6.;
            c = ( 3.*wbb -6.*wb +3.*wt        ) / 6.;
            d = (-1.*wbb +3.*wb -3.*wt +1.*wtt) / 6.;
            a_fineFAB(fiv) = ((d*z + c)*z + b)*z + a;
        } // rbit
    } // cbit
}
#endif // ALLOW_DIVFREEINTERP
