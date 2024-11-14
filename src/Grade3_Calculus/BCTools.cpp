#include "BCTools.H"
#include "BCToolsF_F.H"
#include "Convert.H"
#include "Debug.H"
#include "Subspace.H"

namespace BCTools
{
// ================= BC-setting functions on a single FAB ======================

// NOTE: Unlike the LevelData versions, these functions that set BCs on
// FABS do not check if we are at a computational boundary. They blindly
// fill the ghosts you feed them.


// -----------------------------------------------------------------------------
// Static utility
// Extrapolation boundary conditions for a side.
// This can be used to extrapolate ghosts or the boundary face of a
// FC FArrayBox if you choose your a_valid Box correctly.
// a_state and a_valid must have the same centering.
// If a_state does not contain ghosts, we just exit quietly.
// NOTE: This function does not check if we are at a computational boundary!
// -----------------------------------------------------------------------------
void
extrap(FArrayBox&           a_state,
       const Box&           a_valid,
       const int            a_dir,
       const Side::LoHiSide a_side,
       const int            a_order)
{
    // Sanity checks
    CH_assert(0 <= a_order && a_order <= 4);
    CH_assert(a_state.box().type() == a_valid.type());

    const int isign = sign(a_side);

    // Find the ghost region.
    Box ghostBox = a_valid;
    if (a_valid.type(a_dir) == IndexType::CELL) {
        ghostBox = adjCellBox(ghostBox, a_dir, a_side, 1);
    } else {
        ghostBox = bdryBox(ghostBox, a_dir, a_side);
        ghostBox.shift(a_dir, isign);
    }

    // Does the state have ghosts?
    ghostBox &= a_state.box();
    if (ghostBox.isEmpty()) return;

    // Do the extrapolation.
    FORT_BCTOOLS_EXTRAPSIDE(CHF_FRA(a_state),
                            CHF_BOX(ghostBox),
                            CHF_CONST_INT(a_dir),
                            CHF_CONST_INT(isign),
                            CHF_CONST_INT(a_order));
}


// -----------------------------------------------------------------------------
// Assumes all ghosts except edges and vertices are filled and uses that
// data to extrapolate to edges and vertices.
// NOTE: This only fills ghosts in the first layer.
// -----------------------------------------------------------------------------
void
extrapCorners(FArrayBox& a_state, const Box& a_valid, const int a_order)
{
    auto adjBdryBox =
        [](const Box& bx, const int dir, const Side::LoHiSide& sd) {
            Box retBx(bx);
            if (sd == Side::Lo) {
                retBx.setRange(dir, retBx.smallEnd(dir) - 1, 1);
            } else {
                retBx.setRange(dir, retBx.bigEnd(dir) + 1, 1);
            }
            return retBx;
        };

    Box overlap;

    if (SpaceDim == 2) {
        const int adir = 0;
        const int bdir = 1;

        for (SideIterator asit; asit.ok(); ++asit) {
            const Side::LoHiSide& aside = asit();
            const int             asign = sign(aside);
            const Box             aBox  = adjBdryBox(a_valid, adir, aside);

            for (SideIterator bsit; bsit.ok(); ++bsit) {
                const Side::LoHiSide& bside = bsit();
                const int             bsign = sign(bside);
                const Box             abBox = adjBdryBox(aBox, bdir, bside);

                overlap = abBox & a_state.box();
                if (!overlap.isEmpty()) {
                    FORT_BCTOOLS_EXTRAPEDGE(CHF_FRA(a_state),
                                            CHF_BOX(overlap),
                                            CHF_CONST_INT(adir),
                                            CHF_CONST_INT(asign),
                                            CHF_CONST_INT(bdir),
                                            CHF_CONST_INT(bsign),
                                            CHF_CONST_INT(a_order));
                }
            }  // bsit
        }      // asit

    } else {
        for (int cdir = 0; cdir < SpaceDim; ++cdir) {
            const int adir = (cdir + 1) % SpaceDim;
            const int bdir = (cdir + 2) % SpaceDim;

            for (SideIterator asit; asit.ok(); ++asit) {
                const Side::LoHiSide& aside = asit();
                const int             asign = sign(aside);
                const Box             aBox  = adjBdryBox(a_valid, adir, aside);

                for (SideIterator bsit; bsit.ok(); ++bsit) {
                    const Side::LoHiSide& bside = bsit();
                    const int             bsign = sign(bside);
                    const Box             abBox = adjBdryBox(aBox, bdir, bside);

                    overlap = abBox & a_state.box();
                    if (!overlap.isEmpty()) {
                        FORT_BCTOOLS_EXTRAPEDGE(CHF_FRA(a_state),
                                                CHF_BOX(overlap),
                                                CHF_CONST_INT(adir),
                                                CHF_CONST_INT(asign),
                                                CHF_CONST_INT(bdir),
                                                CHF_CONST_INT(bsign),
                                                CHF_CONST_INT(a_order));
                    }

                    for (SideIterator csit; csit.ok(); ++csit) {
                        const Side::LoHiSide& cside = csit();
                        const int             csign = sign(cside);
                        const Box abcBox = adjBdryBox(abBox, cdir, cside);

                        overlap = abcBox & a_state.box();
                        if (!overlap.isEmpty()) {
                            CH_assert(abcBox.numPts() == 1);
                            const IntVect iv = abcBox.smallEnd();
                            const IntVect ea = -asign * BASISV(adir);
                            const IntVect eb = -bsign * BASISV(bdir);
                            const IntVect ec = -csign * BASISV(cdir);
                            Real          aval, bval, cval;

                            for (int comp = 0; comp < a_state.nComp(); ++comp) {
                                if (a_order == 0) {
                                    aval = a_state(iv + ea, comp);
                                    bval = a_state(iv + eb, comp);
                                    cval = a_state(iv + ec, comp);
                                } else if (a_order == 1) {
                                    aval = 2.0 * a_state(iv + ea, comp) -
                                           a_state(iv + 2 * ea, comp);
                                    bval = 2.0 * a_state(iv + eb, comp) -
                                           a_state(iv + 2 * eb, comp);
                                    cval = 2.0 * a_state(iv + ec, comp) -
                                           a_state(iv + 2 * ec, comp);
                                } else {
                                    aval = 3.0 * a_state(iv + ea, comp) -
                                           3.0 * a_state(iv + 2 * ea, comp) +
                                           a_state(iv + 3 * ea, comp);
                                    bval = 3.0 * a_state(iv + eb, comp) -
                                           3.0 * a_state(iv + 2 * eb, comp) +
                                           a_state(iv + 3 * eb, comp);
                                    cval = 3.0 * a_state(iv + ec, comp) -
                                           3.0 * a_state(iv + 2 * ec, comp) +
                                           a_state(iv + 3 * ec, comp);
                                }
                                a_state(iv, comp) = (aval + bval + cval) / 3.0;
                            }  // comp
                        }
                    }  // csit
                }      // bsit
            }          // asit
        }              // cdir
    }                  // if 2D / 3D
}


// -----------------------------------------------------------------------------
RobinBC_CCSide::RobinBC_CCSide(const int             a_bdryDir,
                               const Side::LoHiSide& a_side,
                               const ProblemDomain&  a_domain,
                               const Real            a_alpha,
                               const Real            a_beta,
                               FArrayBox&            a_bcValsFAB)
: m_bdryDir(a_bdryDir)
, m_side(a_side)
, m_alpha(a_alpha)
, m_beta(a_beta)
{
    const Box targetBox = bdryBox(a_domain.domainBox(), a_bdryDir, a_side);
    CH_verify(targetBox.type() == a_bcValsFAB.box().type());
    CH_verify(a_bcValsFAB.box().contains(targetBox));

    m_bcValsFAB.define(targetBox, a_bcValsFAB.nComp());
    m_bcValsFAB.copy(a_bcValsFAB);
}


// -----------------------------------------------------------------------------
void
RobinBC_CCSide::operator()(FArrayBox&            a_alphaFAB,
                           FArrayBox&            a_betaFAB,
                           FArrayBox&            a_bcFAB,
                           const FArrayBox&      a_stateFAB,
                           const FArrayBox&      a_xFAB,
                           const DataIndex&      a_di,
                           const int             a_bdryDir,
                           const Side::LoHiSide& a_side,
                           const Real            a_time,
                           const bool            a_homogBCs) const
{
    a_alphaFAB.setVal(m_alpha);
    a_betaFAB.setVal(m_beta);

    if (!a_homogBCs) {
        CH_assert(m_bcValsFAB.box().type() == a_bcFAB.box().type());
        CH_assert(m_bcValsFAB.box().contains(a_bcFAB.box()));
        a_bcFAB.copy(m_bcValsFAB);
    }
}


// -----------------------------------------------------------------------------
RobinBC_FCSide::RobinBC_FCSide(const int                        a_bdryDir,
                               const Side::LoHiSide&            a_side,
                               const ProblemDomain&             a_domain,
                               const RealVect&                  a_alpha,
                               const RealVect&                  a_beta,
                               std::array<FArrayBox, SpaceDim>& a_bcValsFAB)
: m_bdryDir(a_bdryDir)
, m_side(a_side)
, m_alpha(a_alpha)
, m_beta(a_beta)
{
    const Box fcTargetBox = bdryBox(a_domain.domainBox(), a_bdryDir, a_side);
    CH_verify(fcTargetBox.type() == a_bcValsFAB[a_bdryDir].box().type());
    CH_verify(a_bcValsFAB[a_bdryDir].box().contains(fcTargetBox));

    m_bcValsFAB[a_bdryDir].define(fcTargetBox, a_bcValsFAB[a_bdryDir].nComp());
    m_bcValsFAB[a_bdryDir].copy(a_bcValsFAB[a_bdryDir]);

    for (int flubComp = 0; flubComp < SpaceDim; ++flubComp) {
        if (flubComp == a_bdryDir) continue;

        const Box ecTargetBox = surroundingNodes(fcTargetBox, flubComp);
        CH_verify(ecTargetBox.type() == a_bcValsFAB[flubComp].box().type());
        CH_verify(a_bcValsFAB[flubComp].box().contains(ecTargetBox));

        CH_verify(a_bcValsFAB[flubComp].nComp() == a_bcValsFAB[a_bdryDir].nComp());

        m_bcValsFAB[flubComp].define(ecTargetBox, a_bcValsFAB[flubComp].nComp());
        m_bcValsFAB[flubComp].copy(a_bcValsFAB[flubComp]);
    }
}


// -----------------------------------------------------------------------------
void
RobinBC_FCSide::operator()(FArrayBox&            a_alpha,
                           FArrayBox&            a_beta,
                           FArrayBox&            a_bcFAB,
                           const FArrayBox&      a_stateFAB,
                           const FArrayBox&      a_xFAB,
                           const DataIndex&      a_di,
                           const int             a_bdryDir,
                           const Side::LoHiSide& a_side,
                           const Real            a_time,
                           const bool            a_homogBCs) const
{
    CH_verify(m_bdryDir == a_bdryDir);
    CH_verify(m_side == a_side);

    const int velComp = getNodalDir(a_stateFAB.box());
    a_alpha.setVal(m_alpha[velComp]);
    a_beta.setVal(m_beta[velComp]);

    if (!a_homogBCs) {
        CH_verify(m_bcValsFAB[velComp].box().type() == a_bcFAB.box().type());
        CH_verify(m_bcValsFAB[velComp].box().contains(a_bcFAB.box()));
        a_bcFAB.copy(m_bcValsFAB[velComp]);
    }
}


// -----------------------------------------------------------------------------
int
getBoxesForApplyBC(Box&                  a_stateBdry,
                   Box&                  a_stateGhosts,
                   const FArrayBox&      a_stateFAB,
                   const Box&            a_ccValid,
                   const int             a_bdryDir,
                   const Side::LoHiSide& a_side)
{
    CH_assert(a_ccValid.type() == IntVect::Zero);

    const int      isign        = sign(a_side);
    const Box&     stateBox     = a_stateFAB.box();
    const IntVect& stateBoxType = stateBox.type();

    // Get state's valid region.
    Box stateValid = a_ccValid;
    stateValid.convert(stateBoxType);
    stateValid &= stateBox;

    // Find the boundary faces.
    a_stateBdry = bdryBox(stateValid, a_bdryDir, a_side, 1);
    if (stateBoxType[a_bdryDir] == IndexType::NODE) {
        a_stateBdry &= stateBox;
    }

    // Find the ghosts.
    a_stateGhosts = a_stateBdry;
    if (stateBoxType[a_bdryDir] == IndexType::CELL) {
        a_stateGhosts.shiftHalf(a_bdryDir, isign);
        a_stateGhosts &= stateBox;
    } else {
        a_stateGhosts.shift(a_bdryDir, isign);
        a_stateGhosts &= stateBox;
    }

    // Return the boundary-normal centering.
    return stateBoxType[a_bdryDir];
}


// -----------------------------------------------------------------------------
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
        const Real            a_time)
{
    const int iHomogBCs = int(a_homogBCs);
    const int isign     = sign(a_side);
    const int numComps  = a_stateFAB.nComp();

    // Get important regions.
    Box stateBdry, stateGhosts;
    const int normIXType = BCTools::getBoxesForApplyBC(
        stateBdry, stateGhosts, a_stateFAB, a_ccValid, a_bdryDir, a_side);

    // Check sanity of input FABs.
    CH_assert(a_xFAB.box().type() == stateBdry.type());
    CH_assert(a_xFAB.box().contains(stateBdry));
    CH_assert(a_xFAB.nComp() == SpaceDim);

    CH_assert(a_dxFAB.box().type() == stateBdry.type());
    CH_assert(a_dxFAB.box().contains(stateBdry));
    CH_assert(a_dxFAB.nComp() == 1);

    // Gather BC data.
    FArrayBox alphaFAB(stateBdry, numComps);
    FArrayBox betaFAB(stateBdry, numComps);
    FArrayBox bcFAB(stateBdry, numComps);
    a_bcFunc(alphaFAB,
             betaFAB,
             bcFAB,
             a_stateFAB,
             a_xFAB,
             a_di,
             a_bdryDir,
             a_side,
             a_time,
             a_homogBCs);
    if (a_homogBCs) bcFAB.setVal(0.0);


    // Determine how many cells are available for the stencil.
    // const int numValidCells = std::min(a_ccValid.size(a_bdryDir),
    //                                    a_stateFAB.box().size(a_bdryDir));
    constexpr int numValidCells = 1;
    BUG("IBs close to the domain boundary don't work with second order BC extrap.");

    // Set BCs on state. This happens differently for different centerings.
    if (normIXType == IndexType::CELL) {
        IntVect shiftIV = IntVect::Zero;
        if (a_side == Side::Lo) shiftIV[a_bdryDir] = 1;

        if (numValidCells >= 2) {
            FORT_BCTOOLS_FILLGHOSTCELLS_MIXEDBC_2CELLS(
                CHF_FRA(a_stateFAB),
                CHF_CONST_FRA_SHIFT(alphaFAB, shiftIV),
                CHF_CONST_FRA_SHIFT(betaFAB, shiftIV),
                CHF_CONST_FRA_SHIFT(bcFAB, shiftIV),
                CHF_CONST_FRA1_SHIFT(a_dxFAB, 0, shiftIV),
                CHF_BOX(stateGhosts),
                CHF_CONST_INT(a_bdryDir),
                CHF_CONST_INT(isign),
                CHF_CONST_INT(iHomogBCs));

        } else {
            FORT_BCTOOLS_FILLGHOSTCELLS_MIXEDBC_1CELL(
                CHF_FRA(a_stateFAB),
                CHF_CONST_FRA_SHIFT(alphaFAB, shiftIV),
                CHF_CONST_FRA_SHIFT(betaFAB, shiftIV),
                CHF_CONST_FRA_SHIFT(bcFAB, shiftIV),
                CHF_CONST_FRA1_SHIFT(a_dxFAB, 0, shiftIV),
                CHF_BOX(stateGhosts),
                CHF_CONST_INT(a_bdryDir),
                CHF_CONST_INT(isign),
                CHF_CONST_INT(iHomogBCs));
        }

    } else {
        if (numValidCells >= 2) {
            // Set faces at boundary.
            FORT_BCTOOLS_FILLBDRYFACES_MIXEDBC_2CELLS(
                CHF_FRA(a_stateFAB),
                CHF_CONST_FRA(alphaFAB),
                CHF_CONST_FRA(betaFAB),
                CHF_CONST_FRA(bcFAB),
                CHF_CONST_FRA1(a_dxFAB, 0),
                CHF_BOX(stateBdry),
                CHF_CONST_INT(a_bdryDir),
                CHF_CONST_INT(isign),
                CHF_CONST_INT(iHomogBCs));

            // Then extrapolate to ghost
            constexpr int extrapOrder = 2;
            FORT_BCTOOLS_EXTRAPSIDE(CHF_FRA(a_stateFAB),
                                    CHF_BOX(stateGhosts),
                                    CHF_CONST_INT(a_bdryDir),
                                    CHF_CONST_INT(isign),
                                    CHF_CONST_INT(extrapOrder));
        } else {
            // Set faces at boundary.
            FORT_BCTOOLS_FILLBDRYFACES_MIXEDBC_1CELL(
                CHF_FRA(a_stateFAB),
                CHF_CONST_FRA(alphaFAB),
                CHF_CONST_FRA(betaFAB),
                CHF_CONST_FRA(bcFAB),
                CHF_CONST_FRA1(a_dxFAB, 0),
                CHF_BOX(stateBdry),
                CHF_CONST_INT(a_bdryDir),
                CHF_CONST_INT(isign),
                CHF_CONST_INT(iHomogBCs));

            // Then extrapolate to ghost
            constexpr int extrapOrder = 1;
            FORT_BCTOOLS_EXTRAPSIDE(CHF_FRA(a_stateFAB),
                                    CHF_BOX(stateGhosts),
                                    CHF_CONST_INT(a_bdryDir),
                                    CHF_CONST_INT(isign),
                                    CHF_CONST_INT(extrapOrder));
        }
    }
}


// =============== BC-setting functions on an entire level =====================

// -----------------------------------------------------------------------------
void
applyBC(LevelData<FArrayBox>&     a_state,
        const bool                a_homogBCs,
        const BCFunction&         a_bcFunc,
        const Real                a_time,
        const GeoSourceInterface& a_geoSrc,
        const RealVect&           a_dXi,
        PhysBdryIter&             a_physBdryIter)
{
    const auto key = a_physBdryIter.lock();
    for (a_physBdryIter.reset(key); a_physBdryIter.ok(); a_physBdryIter.next(key)) {
        const DataIndex&     di       = a_physBdryIter->di;
        const Box&           ccValid  = a_physBdryIter->ccValidBox;
        const int            bdryDir  = a_physBdryIter->dir;
        const Side::LoHiSide side     = a_physBdryIter->side;
        FArrayBox&           stateFAB = a_state[di];

        // Get important regions.
        Box stateBdry, stateGhosts;
        BCTools::getBoxesForApplyBC(
            stateBdry, stateGhosts, stateFAB, ccValid, bdryDir, side);

        // Get physical coordinates at boundary.
        FArrayBox xFAB(stateBdry, SpaceDim);
        a_geoSrc.fill_physCoor(xFAB, a_dXi);

        // Get grid spacing at boundary, dx = dx/dXi * dXi.
        FArrayBox dxFAB(stateBdry, 1);
        a_geoSrc.fill_dxdXi(dxFAB, 0, bdryDir, a_dXi, a_dXi[bdryDir]);

        // Set BCs, fill ghosts.
        BCTools::applyBC(stateFAB,
                         xFAB,
                         dxFAB,
                         ccValid,
                         di,
                         bdryDir,
                         side,
                         a_homogBCs,
                         a_bcFunc,
                         a_time);
    }
    a_physBdryIter.unlock(key);
}


// -----------------------------------------------------------------------------
void
applyBC(LevelData<FluxBox>&  a_state,
        const int            a_fcDir,
        const bool           a_homogBCs,
        const BCFunction&    a_bcFunc,
        const Real           a_time,
        const LevelGeometry& a_levGeo,
        PhysBdryIter&        a_physBdryIter)
{
    const GeoSourceInterface& geoSrc = a_levGeo.getGeoSource();
    const RealVect&           dXi    = a_levGeo.getDXi();

    const auto key = a_physBdryIter.lock();
    for (a_physBdryIter.reset(key); a_physBdryIter.ok(); a_physBdryIter.next(key)) {
        const DataIndex&     di        = a_physBdryIter->di;
        const int            bdryDir   = a_physBdryIter->dir;
        const Side::LoHiSide side      = a_physBdryIter->side;
        FArrayBox&           stateFAB  = a_state[di][a_fcDir];
        const IntVect&       stateType = stateFAB.box().type();

        // Compute regions, centerings, etc.
        IntVect bdryRegionType  = stateType;
        bdryRegionType[bdryDir] = 1;

        Box bdryRegion = a_physBdryIter->fcBdryBox;
        bdryRegion.convert(bdryRegionType);

        Box validRegion = a_physBdryIter->ccValidBox;
        validRegion &= Box(stateFAB.box()).enclosedCells();

        // Get physical coordinates at boundary.
        FArrayBox xFAB(bdryRegion, SpaceDim);
        a_levGeo.fill_physCoor(xFAB);

        // Get grid spacing at boundary, dx = dx/dXi * dXi.
        FArrayBox dxFAB(bdryRegion, 1);
        geoSrc.fill_dxdXi(dxFAB, 0, bdryDir, dXi, dXi[bdryDir]);

        // Set BCs, fill ghosts.
        BCTools::applyBC(stateFAB,
                         xFAB,
                         dxFAB,
                         validRegion,
                         di,
                         bdryDir,
                         side,
                         a_homogBCs,
                         a_bcFunc,
                         a_time);
    }
    a_physBdryIter.unlock(key);
}


// -----------------------------------------------------------------------------
// Static utility
// Extrapolates one layer of ghosts. This does not dellineate between
// exchange, physical boundary, or CFI ghosts. It just fills 'em all.
//
// Order matters! If one or more dirs are periodic:
//  If you exchange before calling extrap, corner ghosts will be valid.
//  If you exchange after, you'll need to use a CornerCopier too.
//  If all dirs are periodic, then you must use a CornerCopier.
//
// If you want to be safe, exchange ALL ghosts before or after extrap.
// -----------------------------------------------------------------------------
void
extrapAllGhosts(LevelData<FArrayBox>& a_state,
                const int             a_order,
                const IntVect&        a_skipGhosts)
{
    const DisjointBoxLayout& grids  = a_state.getBoxes();
    DataIterator             dit    = grids.dataIterator();

    for (dit.reset(); dit.ok(); ++dit) {
        const Box& stateBox  = a_state[dit].box();
        IntVect    ghostVect = a_state.ghostVect();
        int        sumGhosts = ghostVect.sum() - a_skipGhosts.sum();

        Box validDomain = grids[dit];
        validDomain.convert(stateBox.type());
        validDomain.grow(a_skipGhosts);

        while (sumGhosts > 0) {
            for (int dir = 0; dir < SpaceDim; ++dir) {
                const Box valid = stateBox & validDomain;

                BCTools::extrap(a_state[dit], valid, dir, Side::Lo, a_order);
                BCTools::extrap(a_state[dit], valid, dir, Side::Hi, a_order);

                validDomain.grow(dir, 1);
                --ghostVect[dir];
                --sumGhosts;
            }  // end loop over directions (dir)
        }
    }  // end loop over grids (dit)
}


// -----------------------------------------------------------------------------
// Static utility
// FC version.
// -----------------------------------------------------------------------------
void
extrapAllGhosts(LevelData<FluxBox>& a_state,
                const int           a_order,
                const IntVect&      a_skipGhosts)
{
    const DisjointBoxLayout& grids  = a_state.getBoxes();
    DataIterator             dit    = grids.dataIterator();

    for (dit.reset(); dit.ok(); ++dit) {
        for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
            FArrayBox& stateFAB  = a_state[dit][fcDir];
            const Box& stateBox  = stateFAB.box();
            IntVect    ghostVect = a_state.ghostVect();
            int        sumGhosts = ghostVect.sum() - a_skipGhosts.sum();

            Box validDomain = grids[dit];
            validDomain.grow(a_skipGhosts);
            validDomain.surroundingNodes(fcDir);

            while (sumGhosts > 0) {
                for (int i = 0; i < SpaceDim; ++i) {
                    const int bdryDir = (fcDir + i) % SpaceDim;
                    const Box valid   = stateBox & validDomain;

                    BCTools::extrap(
                        stateFAB, valid, bdryDir, Side::Lo, a_order);
                    BCTools::extrap(
                        stateFAB, valid, bdryDir, Side::Hi, a_order);

                    validDomain.grow(bdryDir, 1);
                    --ghostVect[bdryDir];
                    --sumGhosts;
                }  // end loop over directions (i, bdryDir)
            }
        }  // end loop over FC directions (fcDir)
    } // end loop over grids (dit)
}


// -----------------------------------------------------------------------------
// Static utility
// Fill ghosts at edges and vertices of the domain.
// This only fills ghosts in the first layer.
// -----------------------------------------------------------------------------
void
extrapDomainCorners(LevelData<FArrayBox>& a_state, const int a_order)
{
    const IntVect&           ghostVect = a_state.ghostVect();
    const DisjointBoxLayout& grids     = a_state.getBoxes();
    const ProblemDomain&     domain    = grids.physDomain();
    DataIterator             dit       = grids.dataIterator();

    Box domBox = domain.domainBox();
    for (int dir = 0; dir < SpaceDim; ++dir) {
        if (domain.isPeriodic(dir)) {
            domBox.grow(BASISV(dir) * (ghostVect[dir] + 1));
        }
    }

    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& stateFAB = a_state[dit];

        Box ccValid = domBox;
        ccValid &= stateFAB.box();

        BCTools::extrapCorners(stateFAB, ccValid, a_order);
    }
}


// -----------------------------------------------------------------------------
// Static utility
// Fill ghosts at edges and vertices of the domain.
// This only fills ghosts in the first layer.
// -----------------------------------------------------------------------------
void
extrapDomainCorners(LevelData<FluxBox>& a_state, const int a_order)
{
    const IntVect&           ghostVect = a_state.ghostVect();
    const DisjointBoxLayout& grids     = a_state.getBoxes();
    const ProblemDomain&     domain    = grids.physDomain();
    DataIterator             dit       = grids.dataIterator();

    Box domBox = domain.domainBox();
    for (int dir = 0; dir < SpaceDim; ++dir) {
        if (domain.isPeriodic(dir)) {
            domBox.grow(BASISV(dir) * (ghostVect[dir] + 1));
        }
    }

    for (dit.reset(); dit.ok(); ++dit) {
        for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
            FArrayBox& stateFAB = a_state[dit][fcDir];

            Box fcValid = domBox;
            fcValid.surroundingNodes(fcDir);
            fcValid &= stateFAB.box();

            BCTools::extrapCorners(stateFAB, fcValid, a_order);
        }
    }
}


// -----------------------------------------------------------------------------
// Static utility
// -----------------------------------------------------------------------------
void
setValAtPhysBdry(LevelData<FluxBox>& a_flux,
                 const Real          a_val,
                 const SideArray&    a_sideArray)
{
    const int                numComps = a_flux.nComp();
    const DisjointBoxLayout& grids    = a_flux.getBoxes();
    PhysBdryIter             physBdryIter(grids, a_sideArray);

    for (physBdryIter.reset(); physBdryIter.ok(); ++physBdryIter) {
        const DataIndex& di                     = physBdryIter->di;
        const int __attribute__((unused)) isign = physBdryIter->isign;
        const Side::LoHiSide& __attribute__((unused)) side = physBdryIter->side;
        const int  bdryDir                                 = physBdryIter->dir;
        const Box& fcBdryBox = physBdryIter->fcBdryBox;
        FArrayBox& fluxFAB   = a_flux[di][bdryDir];

        fluxFAB.setVal(a_val, fcBdryBox, 0, numComps);
    }
}


// -----------------------------------------------------------------------------
// Static utility
// -----------------------------------------------------------------------------
void
setValAtCFI(LevelData<FluxBox>& a_flux,
            const Real          a_val,
            const CFRegion&     a_cfRegion,
            const SideArray&    a_sideArray)
{
    const int                numComps = a_flux.nComp();
    const DisjointBoxLayout& grids    = a_flux.getBoxes();
    CFIIter                  cfiIter(grids, a_cfRegion, a_sideArray);

    for (cfiIter.reset(); cfiIter.ok(); ++cfiIter) {
        const DataIndex&      di      = cfiIter->di;
        const int             isign   = cfiIter->isign;
        const Side::LoHiSide& side    = cfiIter->side;
        const int             bdryDir = cfiIter->dir;
        const CFIVS&          cfivs   = cfiIter->cfivs;

        FArrayBox& fluxFAB = a_flux[di][bdryDir];

        if (cfivs.isPacked()) {
            Box region = cfivs.packedBox();
            region.shiftHalf(bdryDir, -isign);

            fluxFAB.setVal(a_val, region, 0, numComps);
        } else {
            IVSIterator ivsit = cfivs.getIVS();
            if (side == Side::Lo) {
                const IntVect e = BASISV(bdryDir);
                for (int comp = 0; comp < numComps; ++comp) {
                    for (ivsit.reset(); ivsit.ok(); ++ivsit) {
                        fluxFAB(ivsit() + e, comp) = a_val;
                    }
                }
            } else {
                for (int comp = 0; comp < numComps; ++comp) {
                    for (ivsit.reset(); ivsit.ok(); ++ivsit) {
                        fluxFAB(ivsit(), comp) = a_val;
                    }
                }
            }
        }
    }
}


// -----------------------------------------------------------------------------
// Static utility
// A simple tool for finding the FC dir of a box.
// Returns -1 if the box is not FC.
// -----------------------------------------------------------------------------
int
getNodalDir(const Box& a_box)
{
    // Find the first nodal dir.
    int fcDir  = -1;
    int curDir = 0;
    for (; curDir < SpaceDim; ++curDir) {
        if (a_box.type(curDir) == IndexType::NODE) {
            fcDir = curDir;
            break;
        }
    }

    // Are there any more nodal dirs?
    for (++curDir; curDir < SpaceDim; ++curDir) {
        if (a_box.type(curDir) == IndexType::NODE) {
            return -1;
        }
    }

    return fcDir;
}


};  // end namespace BCTools
