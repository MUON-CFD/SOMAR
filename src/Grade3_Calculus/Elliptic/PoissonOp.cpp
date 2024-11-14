#include "PoissonOp.H"

#include "AnisotropicRefinementTools.H"
#include "BCToolsF_F.H"
#include "CFInterp.H"
#include "CFInterpF_F.H"
#include "Comm.H"
#include "Debug.H"
#include "FABAlgebra.H"
#include "FiniteDiffF_F.H"
#include "Integral.H"
#include "LayoutTools.H"
#include "NeighborIterator.H"
#include "PoissonOpF_F.H"
#include "ProblemContext.H"
#include "ProjectorParameters.H"
#include "SOMAR_Constants.H"
#include "Subspace.H"
#include "LepticBoxTools.H"

// NOTE: The nanChecks in applyOp(...) take the most time by far!
#ifndef NDEBUG
#define nanCheck(x) checkForValidNAN(x)
#else
#define nanCheck(x)
#endif


namespace Elliptic
{


// -----------------------------------------------------------------------------
// Constructor for op at an AMR level, not in the MG depths.
PoissonOp::PoissonOp(const LevelGeometry&                       a_levGeo,
                     const DisjointBoxLayout&                   a_fineGrids,
                     const DisjointBoxLayout&                   a_crseGrids,
                     const int                                  a_numComps,
                     const std::shared_ptr<BCTools::BCFunction> a_bcFuncPtr,
                     const Real                                 a_alpha,
                     const Real                                 a_beta,
                     const LevelData<FluxBox>*                  a_JgupPtr)
: m_crseAMRGrids(a_crseGrids)
, m_geoSrc(a_levGeo.getGeoSource())
, m_whoMadeMe(1)
, m_activeDirs(IntVect::Unit)
, m_bcFuncPtr(a_bcFuncPtr)
, m_L(a_levGeo.getDomainLength())
{
    m_domain      = a_levGeo.getDomain();
    m_grids       = a_levGeo.getBoxes();
    m_dXi         = a_levGeo.getDXi();
    m_numComps    = a_numComps;
    m_relaxMethod = ProblemContext::getInstance()->proj.relaxMethod;

    // J, Jgup cache
    m_Jgup.define(m_grids, 1);
    m_J.define(m_grids, 1);
    if (a_JgupPtr == nullptr) {
        for (DataIterator dit(m_grids); dit.ok(); ++dit) {
            const FluxBox&   srcJgupFlub = a_levGeo.getFCJgup()[dit];
            const FArrayBox& srcJFAB     = a_levGeo.getCCJ()[dit];

            for (int d = 0; d < SpaceDim; ++d) {
                m_Jgup[dit][d].copy(srcJgupFlub[d]);
            }
            m_J[dit].copy(srcJFAB);
        }
    } else {
        for (DataIterator dit(m_grids); dit.ok(); ++dit) {
            const FluxBox&   srcJgupFlub = (*a_JgupPtr)[dit];
            const FArrayBox& srcJFAB     = a_levGeo.getCCJ()[dit];

            for (int d = 0; d < SpaceDim; ++d) {
                m_Jgup[dit][d].copy(srcJgupFlub[d]);
            }
            m_J[dit].copy(srcJFAB);
        }
    }
    nanCheck(m_Jgup);
    nanCheck(m_J);

    // AMR stuff...
    if (m_crseAMRGrids.isClosed()) {
        const Box&    domBox     = m_domain.domainBox();
        const Box&    crseDomBox = m_crseAMRGrids.physDomain().domainBox();
        const IntVect ref        = calculateRefinementRatio(crseDomBox, domBox);
        CH_assert(ref >= IntVect::Unit);
        CH_assert(ref.product() >= 1);
        m_amrCrseDXi = m_dXi * RealVect(ref);
    } else {
        m_amrCrseDXi = RealVect::Zero;
    }

    if (a_fineGrids.isClosed()) {
        const ProblemDomain& fineDomain = a_fineGrids.physDomain();
        const Box&           fineDomBox = fineDomain.domainBox();
        const Box&           domBox     = m_domain.domainBox();
        const IntVect        ref = calculateRefinementRatio(domBox, fineDomBox);
        CH_assert(ref >= IntVect::Unit);
        CH_assert(ref.product() > 1);

        // We shouldn't be refining in an inactive dir.
        CH_assert(m_activeDirs[0] == 1 || ref[0] == 1);
        CH_assert(m_activeDirs[1] == 1 || ref[1] == 1);
        CH_assert(m_activeDirs[SpaceDim - 1] == 1 || ref[SpaceDim - 1] == 1);

        m_fluxReg.define(a_fineGrids, m_grids, fineDomain, ref, a_numComps);
    }

    // BC stuff...
    SideArray activeSides = AllSides;
    for (int d = 0; d < SpaceDim; ++d) {
        if (!m_activeDirs[d]) activeSides[d] = { 0, 0 };
    }
    m_physBdryIter.define(m_grids, activeSides);
    m_cfiIter.define(m_grids, a_levGeo.getCFRegion(), activeSides);
    if (m_crseAMRGrids.isClosed()) {
        m_cfInterp.define(m_grids, m_dXi, m_crseAMRGrids);
    }

    m_exCopier.exchangeDefine(m_grids, m_activeDirs);
    m_exCopier.trimEdges(m_grids, m_activeDirs);

    // This calls cacheMatrixElements and sets m_hasNullSpace.
    this->setAlphaAndBeta(a_alpha, a_beta);
}


// -----------------------------------------------------------------------------
// LepticOperator constructor.
PoissonOp::PoissonOp(const LevelGeometry&     a_levGeo,
                     const DisjointBoxLayout& a_crseGrids,
                     const int                a_numComps)
: m_crseAMRGrids(a_crseGrids)
, m_geoSrc(a_levGeo.getGeoSource())
, m_whoMadeMe(4)
, m_activeDirs(IntVect::Unit)
, m_bcFuncPtr(new BCTools::HomogNeumBC)
, m_L(a_levGeo.getDomainLength())
{
    m_domain      = a_levGeo.getDomain();
    m_grids       = a_levGeo.getBoxes();
    m_dXi         = a_levGeo.getDXi();
    m_numComps    = a_numComps;
    m_relaxMethod = ProblemContext::getInstance()->proj.relaxMethod;

    // J, Jgup cache
    m_Jgup.define(m_grids, 1);
    m_J.define(m_grids, 1);
    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        const FluxBox&   srcJgupFlub = a_levGeo.getFCJgup()[dit];
        const FArrayBox& srcJFAB     = a_levGeo.getCCJ()[dit];

        for (int d = 0; d < SpaceDim; ++d) {
            if (m_activeDirs[d]) {
                m_Jgup[dit][d].copy(srcJgupFlub[d]);
            } else {
                m_Jgup[dit][d].setVal(0.0);
            }
        }
        m_J[dit].copy(srcJFAB);
    }
    nanCheck(m_Jgup);
    nanCheck(m_J);

    // AMR stuff...
    if (m_crseAMRGrids.isClosed()) {
        const Box&    domBox     = m_domain.domainBox();
        const Box&    crseDomBox = m_crseAMRGrids.physDomain().domainBox();
        const IntVect ref        = calculateRefinementRatio(crseDomBox, domBox);
        CH_assert(ref >= IntVect::Unit);
        CH_assert(ref.product() >= 1);
        m_amrCrseDXi = m_dXi * RealVect(ref);
    } else {
        m_amrCrseDXi = RealVect::Zero;
    }

    // m_fluxReg remains undefined.

    // BC stuff...
    SideArray activeSides = AllSides;
    for (int d = 0; d < SpaceDim; ++d) {
        if (!m_activeDirs[d]) activeSides[d] = { 0, 0 };
    }
    m_physBdryIter.define(m_grids, activeSides);
    m_cfiIter.define(m_grids, a_levGeo.getCFRegion(), activeSides);
    if (m_crseAMRGrids.isClosed()) {
        m_cfInterp.define(m_grids, m_dXi, m_crseAMRGrids);
    }

    m_exCopier.exchangeDefine(m_grids, m_activeDirs);
    m_exCopier.trimEdges(m_grids, m_activeDirs);

    // This calls cacheMatrixElements and sets m_hasNullSpace.
    this->setAlphaAndBeta(0.0, 1.0);
}


// -----------------------------------------------------------------------------
// Clones a_srcOp.
// Metric data will be aliased.
// -----------------------------------------------------------------------------
PoissonOp::PoissonOp(const PoissonOp& a_srcOp)
: m_crseAMRGrids(a_srcOp.m_crseAMRGrids)
, m_geoSrc(a_srcOp.m_geoSrc)
, m_whoMadeMe(2)
, m_activeDirs(a_srcOp.m_activeDirs)
, m_domain(a_srcOp.m_domain)
, m_grids(a_srcOp.m_grids)
, m_dXi(a_srcOp.m_dXi)
, m_numComps(a_srcOp.m_numComps)
, m_relaxMethod(a_srcOp.m_relaxMethod)
// J, Jgup copied in body
, m_amrCrseDXi(a_srcOp.m_amrCrseDXi)
// m_fluxReg remains undefined
, m_bcFuncPtr(a_srcOp.m_bcFuncPtr)
// m_physBdryIter, m_cfiIter, m_cfInterp copied in body
, m_exCopier(a_srcOp.m_exCopier)
, m_alpha(a_srcOp.m_alpha)
, m_beta(a_srcOp.m_beta)
// m_M, m_Dinv copied in body
, m_hasNullSpace(a_srcOp.m_hasNullSpace)
// Line relaxation stuff copied in body
, m_L(a_srcOp.m_L)
{
    aliasLevelData(m_Jgup,
                   const_cast<LevelData<FluxBox>*>(&a_srcOp.m_Jgup),
                   a_srcOp.m_Jgup.interval());
    aliasLevelData(m_J,
                   const_cast<LevelData<FArrayBox>*>(&a_srcOp.m_J),
                   a_srcOp.m_J.interval());

    SideArray activeSides = AllSides;
    for (int d = 0; d < SpaceDim; ++d) {
        if (!m_activeDirs[d]) activeSides[d] = { 0, 0 };
    }
    m_physBdryIter.define(m_grids, activeSides);
    m_cfiIter.define(m_grids, a_srcOp.m_cfiIter.getCFRegion(), activeSides);
    if (m_crseAMRGrids.isClosed()) {
        m_cfInterp.define(m_grids, m_dXi, m_crseAMRGrids);
    }

    // Alias the matrix elements
    for (int d = 0; d < SpaceDim; ++d) {
        m_M[d].define(a_srcOp.m_M[d].interval(), (FArrayBox&)(a_srcOp.m_M[d]));
    }

    aliasLevelData(m_Dinv,
                   const_cast<LevelData<FArrayBox>*>(&a_srcOp.m_Dinv),
                   a_srcOp.m_Dinv.interval());

    if (a_srcOp.m_vertTriDiagsLoBCs.isDefined()) {
        aliasLevelData(
            m_vertTriDiagsLoBCs,
            const_cast<LevelData<FArrayBox>*>(&a_srcOp.m_vertTriDiagsLoBCs),
            a_srcOp.m_vertTriDiagsLoBCs.interval());
    }

    if (a_srcOp.m_vertTriDiagsHiBCs.isDefined()) {
        aliasLevelData(
            m_vertTriDiagsHiBCs,
            const_cast<LevelData<FArrayBox>*>(&a_srcOp.m_vertTriDiagsHiBCs),
            a_srcOp.m_vertTriDiagsHiBCs.interval());
    }
}


// -----------------------------------------------------------------------------
PoissonOp::PoissonOp(const PoissonOp&         a_srcOp,
                     const DisjointBoxLayout& a_newGrids)
: m_crseAMRGrids(a_srcOp.m_crseAMRGrids)
, m_geoSrc(a_srcOp.m_geoSrc)
, m_whoMadeMe(6)
, m_activeDirs(a_srcOp.m_activeDirs)
, m_domain(a_srcOp.m_domain)
, m_grids(a_newGrids)
, m_dXi(a_srcOp.m_dXi)
, m_numComps(a_srcOp.m_numComps)
, m_relaxMethod(a_srcOp.m_relaxMethod)
// J, Jgup defined in body
, m_amrCrseDXi(a_srcOp.m_amrCrseDXi)
// m_fluxReg remains undefined
, m_bcFuncPtr(a_srcOp.m_bcFuncPtr)
// m_physBdryIter, m_cfiIter, m_cfInterp defined in body
// m_exCopier defined in body //TODO
// m_alpha and m_beta defined in body
// m_M, m_Dinv defined in body
// m_hasNullSpace defined in body
// Line relaxation stuff defined in body
, m_L(a_srcOp.m_L)
{
    // J, Jgup cache
    m_Jgup.define(m_grids, 1);
    m_J.define(m_grids, 1);
    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        for (int d = 0; d < SpaceDim; ++d) {
            m_geoSrc.fill_Jgup(m_Jgup[dit][d], 0, d, m_dXi);
        }
        m_geoSrc.fill_J(m_J[dit], 0, m_dXi);
    }
    nanCheck(m_Jgup);
    nanCheck(m_J);

    // BC stuff...
    SideArray activeSides = AllSides;
    for (int d = 0; d < SpaceDim; ++d) {
        if (!m_activeDirs[d]) activeSides[d] = { 0, 0 };
    }
    m_physBdryIter.define(m_grids, activeSides);

    CFRegion cfRegion(m_grids, m_domain);
    m_cfiIter.define(m_grids, cfRegion, activeSides);

    if (m_crseAMRGrids.isClosed()) {
        m_cfInterp.define(m_grids, m_dXi, m_crseAMRGrids);
    }

    m_exCopier.exchangeDefine(m_grids, m_activeDirs);
    m_exCopier.trimEdges(m_grids, m_activeDirs);

    // This calls cacheMatrixElements and sets m_hasNullSpace.
    this->setAlphaAndBeta(a_srcOp.m_alpha, a_srcOp.m_beta);
}



// -----------------------------------------------------------------------------
// Defines this as a coarsened version of a_srcOp.
// a_crseGrids must be compatible with a_srcOp::m_grids.
// This constructor is used in MGSolvers to create coarsened ops.
// BCs will be homogeneous.
// -----------------------------------------------------------------------------
PoissonOp::PoissonOp(const PoissonOp&         a_srcOp,
                     const DisjointBoxLayout& a_crseGrids,
                     const IntVect&           a_refRatio)
: m_crseAMRGrids()
, m_geoSrc(a_srcOp.m_geoSrc)
, m_whoMadeMe(3)
, m_activeDirs(a_srcOp.m_activeDirs)
// m_domain, m_grids, m_dXi coarsened in body
, m_numComps(a_srcOp.m_numComps)
, m_relaxMethod(a_srcOp.m_relaxMethod)
// m_Jgup, m_J coarsened in body
, m_amrCrseDXi(a_srcOp.m_amrCrseDXi)
// m_fluxReg remains undefined
, m_bcFuncPtr(a_srcOp.m_bcFuncPtr)
// remains undefined
// m_physBdryIter, m_cfiIter, m_exCopier coarsened in body
// m_cfInterp remains undefined
, m_alpha(a_srcOp.m_alpha)
, m_beta(a_srcOp.m_beta)
// m_M, m_Dinv defined in body
, m_hasNullSpace(a_srcOp.m_hasNullSpace)
// Line relaxation stuff coarsened in body
, m_L(a_srcOp.m_L)
{
    CH_assert(a_srcOp.m_grids.compatible(a_crseGrids));
    CH_assert(a_refRatio >= IntVect::Unit);
    CH_assert(a_refRatio != IntVect::Unit);

    m_domain = a_crseGrids.physDomain();
    m_grids  = a_crseGrids;
    m_dXi    = a_srcOp.m_dXi * RealVect(a_refRatio);

    m_Jgup.define(m_grids, 1);
    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        CFInterp::localCoarsen(m_Jgup[dit],
                               a_srcOp.m_Jgup[dit],
                               m_Jgup[dit].box(),
                               a_refRatio,
                               nullptr);
    }  // dit
    nanCheck(m_Jgup);

    m_J.define(m_grids, 1);
    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        CFInterp::localCoarsen(m_J[dit],
                               a_srcOp.m_J[dit],
                               m_J[dit].box(),
                               a_refRatio,
                               false,
                               nullptr);
    }  // dit
    nanCheck(m_J);

    SideArray activeSides = AllSides;
    for (int d = 0; d < SpaceDim; ++d) {
        if (!m_activeDirs[d]) activeSides[d] = { 0, 0 };
    }
    m_physBdryIter.define(m_grids, activeSides);
    {
        CFRegion cfRegion(a_srcOp.m_cfiIter.getCFRegion());
        ::coarsen(cfRegion, a_refRatio);
        m_cfiIter.define(m_grids, cfRegion, activeSides);
    }
    // m_cfInterp remains undefined

    m_exCopier = a_srcOp.m_exCopier;
    ::coarsen(m_exCopier, a_refRatio);

    this->cacheMatrixElements();

    m_HOProlongCornerCopier.define(
        m_grids, m_grids, m_domain, m_activeDirs, true);
}


// -----------------------------------------------------------------------------
// Define a horizontal-only op
PoissonOp::PoissonOp(const GeoSourceInterface&                  a_geoSrc,
                     const DisjointBoxLayout&                   a_hGrids,
                     const DisjointBoxLayout&                   a_hCrseGrids,
                     const int                                  a_numComps,
                     const RealVect&                            a_L,
                     const RealVect&                            a_dXi,
                     const std::shared_ptr<BCTools::BCFunction> a_bcFuncPtr,
                     const int                                  a_relaxMethod)
: m_crseAMRGrids(a_hCrseGrids) // TODO: Is this needed?
, m_geoSrc(a_geoSrc)
, m_whoMadeMe(5)
, m_activeDirs(s_hunit)
, m_bcFuncPtr(a_bcFuncPtr)
, m_L(a_L)
{
    m_grids       = a_hGrids;
    m_domain      = a_hGrids.physDomain();
    m_dXi         = a_dXi;
    m_dXi[SpaceDim - 1] = 1.0;
    m_numComps    = a_numComps;
    m_relaxMethod = a_relaxMethod;

    // J, Jgup cache
    m_Jgup.define(m_grids, 1);
    m_J.define(m_grids, 1);
    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        m_Jgup[dit][SpaceDim - 1].setVal(0.0);

        // These are vertical averages, so the dz/dZeta factors are 1.
        if constexpr (SpaceDim == 2) {
            // Jg^{00} = dXi/dx
            m_geoSrc.fill_dXidx(m_Jgup[dit][0], 0, 0, m_dXi);

            // J = dx/dXi
            m_geoSrc.fill_dxdXi(m_J[dit], 0, 0, m_dXi);
        } else {
            // Jg^{00} = dXi/dx * dy/dEta
            FArrayBox tmpFAB(m_Jgup[dit][0].box(), 1);
            m_geoSrc.fill_dXidx(m_Jgup[dit][0], 0, 0, m_dXi);
            m_geoSrc.fill_dxdXi(tmpFAB, 0, 1, m_dXi);
            m_Jgup[dit][0] *= tmpFAB;

            // Jg^{11} = dEta/dy * dx/dXi
            tmpFAB.define(m_Jgup[dit][1].box(), 1);
            m_geoSrc.fill_dXidx(m_Jgup[dit][1], 0, 1, m_dXi);
            m_geoSrc.fill_dxdXi(tmpFAB, 0, 0, m_dXi);
            m_Jgup[dit][1] *= tmpFAB;

            // J = dx/dXi * dy/dEta
            tmpFAB.define(m_J[dit].box(), 1);
            m_geoSrc.fill_dxdXi(m_J[dit], 0, 0, m_dXi);
            m_geoSrc.fill_dxdXi(tmpFAB, 0, 1, m_dXi);
            m_J[dit] *= tmpFAB;
        }
    }
    nanCheck(m_Jgup);
    nanCheck(m_J);

    // AMR stuff...
    if (m_crseAMRGrids.isClosed()) {
        const Box&    domBox     = m_domain.domainBox();
        const Box&    crseDomBox = m_crseAMRGrids.physDomain().domainBox();
        const IntVect ref        = calculateRefinementRatio(crseDomBox, domBox);
        CH_assert(ref >= IntVect::Unit);
        CH_assert(ref.product() >= 1);
        m_amrCrseDXi = m_dXi * RealVect(ref);
    } else {
        m_amrCrseDXi = RealVect::Zero;
    }

    // m_fluxReg remains undefined.

    // BC stuff...
    SideArray activeSides = AllSides;
    for (int d = 0; d < SpaceDim; ++d) {
        if (!m_activeDirs[d]) activeSides[d] = { 0, 0 };
    }

    m_physBdryIter.define(m_grids, activeSides);

    CFRegion cfRegion(m_grids, m_domain);
    m_cfiIter.define(m_grids, cfRegion, activeSides);

    if (m_crseAMRGrids.isClosed()) {
        m_cfInterp.define(m_grids, m_dXi, m_crseAMRGrids);
    }

    m_exCopier.exchangeDefine(m_grids, m_activeDirs);
    m_exCopier.trimEdges(m_grids, m_activeDirs);

    // This calls cacheMatrixElements and sets m_hasNullSpace.
    this->setAlphaAndBeta(0.0, 1.0);
}


// -----------------------------------------------------------------------------
// Redefines and recomputes m_invDiags = 1.0 / diagonal matrix elements.
// -----------------------------------------------------------------------------
void
PoissonOp::cacheMatrixElements()
{
    for (int d = 0; d < SpaceDim; ++d) {
        const IntVect e         = BASISV(d);
        const Box&    domainBox = m_domain.domainBox();
        const Box     ccRegion  = Subspace::flattenBox(domainBox, e);
        const Box     fcRegion  = surroundingNodes(ccRegion, d);

        FArrayBox fcDXidxFAB(fcRegion, 1);
        m_geoSrc.fill_dXidx(fcDXidxFAB, 0, d, m_dXi);

        FArrayBox ccDXidxFAB(ccRegion, 1);
        m_geoSrc.fill_dXidx(ccDXidxFAB, 0, d, m_dXi);

        m_M[d].define(ccRegion, 2);  // = 0 lower, 1 = upper diagonals
        if (m_activeDirs[d]) {
            FORT_POISSONOP_COMPUTEMATRIXELEMENTS(CHF_FRA(m_M[d]),
                                                 CHF_CONST_FRA1(fcDXidxFAB, 0),
                                                 CHF_CONST_FRA1(ccDXidxFAB, 0),
                                                 CHF_CONST_REAL(m_dXi[d]),
                                                 CHF_CONST_INT(d),
                                                 CHF_BOX(ccRegion));
        } else {
            m_M[d].setVal(0.0);
        }
    }

    m_Dinv.define(m_grids, 1);
    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        FArrayBox& DinvFAB = m_Dinv[dit];
        FORT_POISSONOP_COMPUTEDINV(CHF_FRA1(DinvFAB, 0),
                                   CHF_CONST_FRA(m_M[0]),
                                   CHF_CONST_FRA(m_M[1]),
                                   CHF_CONST_FRA(m_M[SpaceDim - 1]),
                                   CHF_CONST_FRA1(m_J[dit], 0),
                                   CHF_CONST_REAL(m_alpha),
                                   CHF_CONST_REAL(m_beta),
                                   CHF_BOX(DinvFAB.box()));
    }

    // Vertical line relaxation stuff.
    if (m_relaxMethod == ProjectorParameters::RelaxMethod::VERTLINE) {
        if (!m_activeDirs[SpaceDim - 1]) {
            MAYDAYERROR(
                "The vertical must be an active direction to use vertical line "
                "relaxation!");
        }

        // In fact, for now...
        CH_assert(m_activeDirs == IntVect::Unit);

        // Make sure grids are appropriate for line relaxation.
        NeighborIterator nit(m_grids);
        for (DataIterator dit(m_grids); dit.ok(); ++dit) {
            for (SideIterator sit; sit.ok(); ++sit) {
                const Box ghostBox = bdryBox(m_grids[dit], SpaceDim - 1, sit());
                for (nit.begin(dit()); nit.ok(); ++nit) {
                    if (!nit.box().intersectsNotEmpty(ghostBox)) continue;
                    MAYDAYERROR(
                        "Grids are not suitable for vertical line relaxation. "
                        "Try setting base.splitDirs = 1 1 0 in 3D or 1 0 in "
                        "2D.");
                }
            }
        }

        DisjointBoxLayout loBCGrids, hiBCGrids;
        adjCellLo(loBCGrids, m_grids, SpaceDim - 1, -1);
        adjCellHi(hiBCGrids, m_grids, SpaceDim - 1, -1);
        m_vertTriDiagsLoBCs.define(loBCGrids, m_numComps);
        m_vertTriDiagsHiBCs.define(hiBCGrids, m_numComps);

        for (DataIterator dit(m_grids); dit.ok(); ++dit) {
            FArrayBox&       loBCFAB = m_vertTriDiagsLoBCs[dit];
            FArrayBox&       hiBCFAB = m_vertTriDiagsHiBCs[dit];
            const FArrayBox& JFAB    = m_J[dit];
            const Box        valid   = m_grids[dit];

            FArrayBox xFAB(valid, SpaceDim);
            m_geoSrc.fill_physCoor(xFAB, m_dXi);

            FArrayBox      dummyFAB;
            constexpr Real dummyTime = 0.0;
            constexpr bool homogBCs  = true;

            FArrayBox phiFAB(valid, m_numComps);
            phiFAB.setVal(quietNAN);

            // Get important regions.
            const Box loBdry = bdryBox(valid, SpaceDim - 1, Side::Lo);
            FArrayBox loBCalphaFAB(loBdry, m_numComps);
            FArrayBox loBCbetaFAB(loBdry, m_numComps);
            (*m_bcFuncPtr)(loBCalphaFAB,
                           loBCbetaFAB,
                           dummyFAB,
                           phiFAB,
                           xFAB,
                           dit(),
                           SpaceDim - 1,
                           Side::Lo,
                           dummyTime,
                           homogBCs);
            const Real lodz = m_dXi[SpaceDim - 1];

            const Box hiBdry = bdryBox(valid, SpaceDim - 1, Side::Hi);
            FArrayBox hiBCalphaFAB(hiBdry, phiFAB.nComp());
            FArrayBox hiBCbetaFAB(hiBdry, phiFAB.nComp());
            (*m_bcFuncPtr)(hiBCalphaFAB,
                           hiBCbetaFAB,
                           dummyFAB,
                           phiFAB,
                           xFAB,
                           dit(),
                           SpaceDim - 1,
                           Side::Hi,
                           dummyTime,
                           homogBCs);
            const Real hidz = m_dXi[SpaceDim - 1];

            // We want to preserve the horizontal index to compute
            // the red-black ordering.
            IntVect validShift       = IntVect::Zero;
            validShift[SpaceDim - 1] = valid.smallEnd(SpaceDim - 1);

            if constexpr (SpaceDim == 2) {
                FORT_POISSONOP_DEFINEVERTLINERELAXBCS_2D(
                    CHF_CONST_FRA1_SHIFT(JFAB, 0, validShift),
                    CHF_CONST_FRA_SHIFT(m_M[SpaceDim - 1], validShift),
                    CHF_CONST_REAL(m_beta),
                    CHF_BOX_SHIFT(valid, validShift),
                    CHF_CONST_FRA_SHIFT(loBCalphaFAB, validShift),
                    CHF_CONST_FRA_SHIFT(loBCbetaFAB, validShift),
                    CHF_CONST_REAL(lodz),
                    CHF_CONST_FRA_SHIFT(hiBCalphaFAB, validShift),
                    CHF_CONST_FRA_SHIFT(hiBCbetaFAB, validShift),
                    CHF_CONST_REAL(hidz),
                    CHF_FRA_SHIFT(loBCFAB, validShift),
                    CHF_FRA_SHIFT(hiBCFAB, validShift));
            } else {
                FORT_POISSONOP_DEFINEVERTLINERELAXBCS_3D(
                    CHF_CONST_FRA1_SHIFT(JFAB, 0, validShift),
                    CHF_CONST_FRA_SHIFT(m_M[SpaceDim - 1], validShift),
                    CHF_CONST_REAL(m_beta),
                    CHF_BOX_SHIFT(valid, validShift),
                    CHF_CONST_FRA_SHIFT(loBCalphaFAB, validShift),
                    CHF_CONST_FRA_SHIFT(loBCbetaFAB, validShift),
                    CHF_CONST_REAL(lodz),
                    CHF_CONST_FRA_SHIFT(hiBCalphaFAB, validShift),
                    CHF_CONST_FRA_SHIFT(hiBCbetaFAB, validShift),
                    CHF_CONST_REAL(hidz),
                    CHF_FRA_SHIFT(loBCFAB, validShift),
                    CHF_FRA_SHIFT(hiBCFAB, validShift));
            }
        }  // dit
    }
}


// -----------------------------------------------------------------------------
bool
PoissonOp::checkForNullSpace() const
{
    LevelData<FArrayBox> ones(m_grids, 1, IntVect::Unit);
    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        ones[dit].setVal(1.0);
    }

    LevelData<FArrayBox> Lones(m_grids, 1);
    this->applyOp(Lones, ones, nullptr, 0.0, true, true);

    int hasNonZero;
    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        const FArrayBox& LonesFAB = Lones[dit];
        const Box&       valid    = m_grids[dit];

        FORT_POISSONOP_NONZEROSEARCH(CHF_INT(hasNonZero),
                                     CHF_CONST_FRA1(LonesFAB, 0),
                                     CHF_BOX(valid),
                                     CHF_CONST_REAL(smallReal));
        if (hasNonZero) break;
    }
    Comm::reduce(hasNonZero, MPI_SUM);

    return (hasNonZero == 0);
}


// -----------------------------------------------------------------------------
// Destructor for good measure.
// -----------------------------------------------------------------------------
PoissonOp::~PoissonOp()
{
}


// -----------------------------------------------------------------------------
void
PoissonOp::setAlphaAndBeta(const Real a_alpha, const Real a_beta)
{
    m_alpha = a_alpha;
    m_beta  = a_beta;

    this->cacheMatrixElements();

    if (RealCmp::neq(m_alpha, 0.0)) {
        m_hasNullSpace = false;
    } else {
        m_hasNullSpace = this->checkForNullSpace();
    }
}


// ======================== LinearOperator overrides ===========================

// -----------------------------------------------------------------------------
void
PoissonOp::applyBCsBegin(LevelData<FArrayBox>&       a_phi,
                         const LevelData<FArrayBox>* a_crsePhiPtr,
                         const Real                  a_time,
                         const bool                  a_homogPhysBCs,
                         const bool                  a_homogCFIBCs) const
{
    CH_assert(a_phi.ghostVect() >= m_activeDirs);
    CH_assert(a_phi.nComp() == m_numComps);
    CH_assert(a_phi.getBoxes().physDomain() == m_domain);
    CH_assert(a_phi.getBoxes().compatible(m_grids));
    nanCheck(a_phi);
    CH_assert(m_cfiIter.isDefined());

    // Begin exchange
    a_phi.exchangeBegin(m_exCopier);

    // CFI BCs
    if (!m_cfiIter.isEmpty()) {
        if (a_homogCFIBCs) {
            CFInterp::homogInterpAtCFI(a_phi, m_dXi, m_amrCrseDXi, m_cfiIter);
        } else {
            // if MG depth > 0, we should be calling homogInterpAtCFI!
            CH_assert(a_crsePhiPtr);
            CH_assert(a_crsePhiPtr->getBoxes() == m_crseAMRGrids);

            nanCheck(*a_crsePhiPtr);
            m_cfInterp.interpAtCFI(a_phi, *a_crsePhiPtr);
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
}


// -----------------------------------------------------------------------------
void
PoissonOp::applyOpNoBCs(LevelData<FArrayBox>&       a_lhs,
                        const LevelData<FArrayBox>& a_phi) const
{
    // I don't think we need these checks in all of the functions.
    // Putting them here should be enough to catch an error.
    CH_assert(a_lhs.getBoxes().physDomain() == m_domain);
    CH_assert(a_phi.getBoxes().physDomain() == m_domain);
    CH_assert(a_lhs.getBoxes().compatible(m_grids));
    CH_assert(a_phi.getBoxes().compatible(m_grids));
    CH_assert(a_lhs.nComp() == m_numComps);
    CH_assert(a_phi.nComp() == m_numComps);
    CH_assert(a_phi.ghostVect() >= m_activeDirs);
    nanCheck(a_phi);

    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        if (m_activeDirs == IntVect::Unit) {
            FORT_POISSONOP_APPLYOP(CHF_FRA(a_lhs[dit]),
                                   CHF_CONST_FRA(a_phi[dit]),
                                   CHF_CONST_FRA(m_M[0]),
                                   CHF_CONST_FRA(m_M[1]),
                                   CHF_CONST_FRA(m_M[SpaceDim - 1]),
                                   CHF_CONST_FRA1(m_Dinv[dit], 0),
                                   CHF_CONST_FRA1(m_J[dit], 0),
                                   CHF_CONST_REALVECT(m_dXi),
                                   CHF_CONST_REAL(m_beta),
                                   CHF_BOX(m_grids[dit]));
        } else {
            FORT_POISSONOP_APPLYOPDIRS(CHF_FRA(a_lhs[dit]),
                                       CHF_CONST_FRA(a_phi[dit]),
                                       CHF_CONST_FRA(m_M[0]),
                                       CHF_CONST_FRA(m_M[1]),
                                       CHF_CONST_FRA(m_M[SpaceDim - 1]),
                                       CHF_CONST_FRA1(m_Dinv[dit], 0),
                                       CHF_CONST_FRA1(m_J[dit], 0),
                                       CHF_CONST_REALVECT(m_dXi),
                                       CHF_CONST_REAL(m_beta),
                                       CHF_BOX(m_grids[dit]),
                                       CHF_CONST_INTVECT(m_activeDirs));
        }
    }  // dit

    nanCheck(a_lhs);
}


// -----------------------------------------------------------------------------
// If this operator has a null space, a_phi needs to be projected.
// In other words, if L[e^1] = 0, and phi = c_1*e^1 + c_2*e^2 + ...,
// then this function should remove c_1*e^1.
// (Typically, you would just remove the average.)
// -----------------------------------------------------------------------------
void
PoissonOp::removeKernel(LevelData<FArrayBox>& a_phi) const
{
    nanCheck(a_phi);

    if (!m_hasNullSpace) return;

    // Compute the average...
    Real vol = 0.0;
    Real sum = 0.0;
    {
        CH_assert(a_phi.nComp() == 1);
        CH_assert(a_phi.getBoxes().physDomain() == m_domain);
        CH_assert(a_phi.getBoxes().compatible(m_grids));

        Integral::sum(sum, vol, a_phi, &m_J, m_dXi, 0);
    }

    // Remove the average...
    const Real avgPhi = sum / vol;
    for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit) {
        a_phi[dit] -= avgPhi;
    }

    nanCheck(a_phi);
}


// -----------------------------------------------------------------------------
// Checks if L[phi] = rhs is a solvable problem.
// For this op, we check if Sum[rhs] = 0.
// -----------------------------------------------------------------------------
bool
PoissonOp::levelEquationIsConsistent(const StateType& a_phi,
                                     const StateType* /*a_crsePhiPtr*/,
                                     const StateType& a_rhs,
                                     const Real /*a_time*/,
                                     const bool /*a_homogPhysBCs*/,
                                     const bool /*a_homogCFIBCs*/) const
{
    CH_assert(a_rhs.nComp() == 1);
    CH_assert(a_rhs.getBoxes().physDomain() == m_domain);
    CH_assert(a_rhs.getBoxes().compatible(m_grids));

    if (RealCmp::eq(m_alpha, 0.0)) {
        // Compute the sums.
        Real vol    = 0.0;
        Real phiSum = 0.0;
        Real rhsSum = 0.0;

        // Notify the user...
        if (m_hasNullSpace) {
            Integral::sum(phiSum, vol, a_phi, &m_J, m_dXi, 0);
            pout() << "Sum[phi] = " << phiSum << endl;
        }

        // Do not scale rhs by J. Remember, this op's divergence already does
        // that.
        Integral::sum(rhsSum, vol, a_rhs, nullptr, m_dXi, 0);
        pout() << "Sum[rhs] = " << rhsSum << endl;
    }

    // ...but don't stop the solve.
    return true;
}


// -----------------------------------------------------------------------------
// Sets phi to a preconditioned value.
// Typically, this will be phi = rhs / invDiags, then relaxed, but you
// can ignore a_relaxIters if you have a better idea.
// -----------------------------------------------------------------------------
void
PoissonOp::preCond(LevelData<FArrayBox>&       a_phi,
                   const LevelData<FArrayBox>& a_rhs,
                   const Real                  a_time,
                   const int                   a_relaxIters) const
{
    nanCheck(a_rhs);

    this->assignLocal(a_phi, a_rhs);

    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        for (int comp = 0; comp < m_numComps; ++comp) {
            a_phi[dit].mult(m_Dinv[dit], m_grids[dit], 0, comp, 1);
        }
    }
    this->relax(a_phi, a_rhs, a_time, a_relaxIters);

    nanCheck(a_phi);
}


// -----------------------------------------------------------------------------
// Applies relaxation to the residual equation.
// If you don't plan to use this in preCond, or if you aren't creating
// an MG or AMR op, then this function can be a simple no-op.
// -----------------------------------------------------------------------------
void
PoissonOp::relax(LevelData<FArrayBox>&       a_cor,
                 const LevelData<FArrayBox>& a_res,
                 const Real                  a_time,
                 const int                   a_iters) const
{
    nanCheck(a_cor);
    nanCheck(a_res);

    switch (m_relaxMethod) {
        case ProjectorParameters::RelaxMethod::NONE:
            break;
        case ProjectorParameters::RelaxMethod::POINT:
            this->point_relax(a_cor, a_res, a_time, a_iters);
            break;
        case ProjectorParameters::RelaxMethod::JACOBI:
            this->jacobi_relax(a_cor, a_res, a_time, a_iters);
            break;
        case ProjectorParameters::RelaxMethod::JACOBIRB:
            this->jacobiRB_relax(a_cor, a_res, a_time, a_iters);
            break;
        case ProjectorParameters::RelaxMethod::GS:
            this->gs_relax(a_cor, a_res, a_time, a_iters);
            break;
        case ProjectorParameters::RelaxMethod::GSRB:
            this->gsrb_relax(a_cor, a_res, a_time, a_iters);
            break;
        case ProjectorParameters::RelaxMethod::VERTLINE:
            this->vertLineGSRB_relax(a_cor, a_res, a_time, a_iters);
            break;
        default:
            MAYDAYERROR("Unknown relaxation method.");
            break;
    }

    nanCheck(a_cor);
}


// =========================== MGOperator methods ==============================

// -----------------------------------------------------------------------------
// Factory method.
// Allocate + define a coarsened version of *this.
// Deletion is left to the caller.
//
// a_refRatio specifies how much *coarsening* will be needed to get from
// this to the new operator.
//
// The new operator will only be used with homogeneous BCs.
// -----------------------------------------------------------------------------
MGOperator<LevelData<FArrayBox>>*
PoissonOp::newMGOperator(const IntVect& a_refRatio) const
{
    PoissonOp* newOpPtr = nullptr;

    if (a_refRatio == IntVect::Unit) {
        newOpPtr = new PoissonOp(*this);
    } else {
        DisjointBoxLayout crseGrids;
        ::coarsen(crseGrids, m_grids, a_refRatio);

        newOpPtr = new PoissonOp(*this, crseGrids, a_refRatio);
    }

    return newOpPtr;
}


// -----------------------------------------------------------------------------
// Restrict to coarser MG depth: a_crseRes = I[h->2h](a_fineRes).
// This op is at the fine level.
// a_crseRes and a_fineRes should be defined over compatible grids.
// -----------------------------------------------------------------------------
void
PoissonOp::MGRestrict(StateType&       a_crseRes,
                      const StateType& a_fineRes,
                      const Real /*a_time*/,
                      const IntVect& a_refRatio,
                      const MGOpType& /*a_crseOp*/) const
{
    //CH_assert(m_activeDirs == IntVect::Unit);

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
// Prolong (interpolate) a_crse to finer MG depth and add to a_fine:
//   a_finePhi += I[2h->h](a_crseCor).
// This op is at the fine level.
// a_crse and a_fine should be defined over compatible grids.
//
// We do not overwrite a_finePhi! We ADD the correction.
//
// If a_interpOrder > 1, then all ghosts of a_crseCor must be extrapolated,
// including edge and vertex ghosts.
// -----------------------------------------------------------------------------
void
PoissonOp::MGProlong(StateType&      a_finePhi,
                     StateType&      a_crseCor,
                     const Real      a_time,
                     const IntVect&  a_refRatio,
                     const MGOpType& a_crseOp,
                     const int       a_interpOrder) const
{
    // CH_assert(m_activeDirs == IntVect::Unit);

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

        FORT_MGOPERATOR_PROLONG_CONSTANT(CHF_FRA_SHIFT(fineFAB, fiv),
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
        bool           noRoom = false;
        LayoutIterator lit    = crseGrids.layoutIterator();
        for (lit.reset(); lit.ok(); ++lit) {
            const IntVect& size = crseGrids[lit].size();
            if (!(size >= IntVect(D_DECL(4, 4, 4)))) {
                noRoom = true;
                break;
            }
        }
        if (noRoom) break;

        // Update with d^2/dx^2, d^2/dy^2, and d^2/dz^2 terms.
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox&       fineFAB   = a_finePhi[dit];
            const FArrayBox& crseFAB   = a_crseCor[dit];
            const Box&       crseValid = crseGrids[dit];

            FORT_MGOPERATOR_PROLONG_QUADUPGRADE1(CHF_FRA(fineFAB),
                                                 CHF_CONST_FRA(crseFAB),
                                                 CHF_BOX(crseValid),
                                                 CHF_BOX(refBox),
                                                 CHF_CONST_INTVECT(a_refRatio));
        }

        // Fill *ALL* ghosts. Some are already filled, the corners are not.
        for (dit.reset(); dit.ok(); ++dit) {
            BCTools::extrapCorners(a_crseCor[dit], crseGrids[dit], 2);
        }
        const auto& crsePoissonOp = static_cast<const PoissonOp&>(a_crseOp);
        CH_verify(crsePoissonOp.m_HOProlongCornerCopier.isDefined());
        a_crseCor.exchange(crsePoissonOp.m_HOProlongCornerCopier);

        // Update with mixed second derivative terms.
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox&       fineFAB   = a_finePhi[dit];
            const FArrayBox& crseFAB   = a_crseCor[dit];
            const Box&       crseValid = crseGrids[dit];

            FORT_MGOPERATOR_PROLONG_QUADUPGRADE2(CHF_FRA(fineFAB),
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

// ======================== AMRMGOperator overrides ============================

// -----------------------------------------------------------------------------
// Apply the AMR operator, including coarse-fine matching.
// The CFI BCs are always inhomogeneous in this function.
// -----------------------------------------------------------------------------
void
PoissonOp::AMROperator(LevelData<FArrayBox>& a_LofPhi,
                       LevelData<FArrayBox>& a_phiFine,  // You may set ghosts.
                       LevelData<FArrayBox>& a_phi,      // You may set ghosts.
                       const LevelData<FArrayBox>& a_phiCoarse,
                       const IntVect&              a_fineRefRatio,
                       const IntVect& /*a_crseRefRatio*/,
                       const Real         a_time,
                       const bool         a_homogPhysBCs,
                       const AMRMGOpType& a_finerOp) const
{
    CH_assert(m_activeDirs == IntVect::Unit);
    this->applyOp(a_LofPhi, a_phi, &a_phiCoarse, a_time, a_homogPhysBCs, false);
    this->reflux(a_LofPhi, a_phiFine, a_phi, a_fineRefRatio, a_finerOp);
}


// -----------------------------------------------------------------------------
// Apply the AMR operator, including coarse-fine matching.
// The CFI BCs are always inhomogeneous in this function.
// Assume no finer AMR level.
// -----------------------------------------------------------------------------
void
PoissonOp::AMROperatorNF(LevelData<FArrayBox>& a_LofPhi,
                         LevelData<FArrayBox>& a_phi,  // You may set ghosts.
                         const LevelData<FArrayBox>& a_phiCoarse,
                         const IntVect& /*a_crseRefRatio*/,
                         const Real a_time,
                         const bool a_homogPhysBCs) const
{
    CH_assert(m_activeDirs == IntVect::Unit);
    this->applyOp(a_LofPhi, a_phi, &a_phiCoarse, a_time, a_homogPhysBCs, false);
}


// -----------------------------------------------------------------------------
// Apply the AMR operator, including coarse-fine matching.
// Assume no coarser AMR level.
// -----------------------------------------------------------------------------
void
PoissonOp::AMROperatorNC(
    LevelData<FArrayBox>& a_LofPhi,
    LevelData<FArrayBox>& a_phiFine,  // You may set ghosts.
    LevelData<FArrayBox>& a_phi,      // You may set ghosts.
    const IntVect&        a_fineRefRatio,
    const Real            a_time,
    const bool            a_homogPhysBCs,
    const AMRMGOpType&    a_finerOp) const
{
    CH_assert(m_activeDirs == IntVect::Unit);
    this->applyOp(a_LofPhi, a_phi, nullptr, a_time, a_homogPhysBCs, false);
    this->reflux(a_LofPhi, a_phiFine, a_phi, a_fineRefRatio, a_finerOp);
}


// -----------------------------------------------------------------------------
// \brief
//  Compute the p-Norm of a_x only on where not covered by a_fineRes.
//  If a_p = 0, then we compute an inf-norm.
// \details
//  This function must perform all MPI communication.
//  We want the norm on this level to be comparable to the norms on
//  other levels, so we choose powScale = dx*dy*dz and compute the norm as
//   |a_res|_p = ( Sum |a_res|^p * powScale ) ^(1/p).
//
//  Sorry if the name is confusing. It means we want the norm of a single
//  level, but that level may be covered by another level in an AMR solve.
//  If a_fineResPtr is nullptr, then we just need plain ol' norm().
// -----------------------------------------------------------------------------
Real
PoissonOp::AMRNormLevel(const StateType& a_res,
                        const StateType* a_fineResPtr,
                        const IntVect&   a_refRatio,
                        const int        a_p) const
{
    CH_assert(m_activeDirs == IntVect::Unit);

    const Real powScale = this->getDXi().product();

    if (!a_fineResPtr) {
        return this->norm(a_res, a_p, powScale);
    }

    nanCheck(a_res);

    const int                numComps = a_res.nComp();
    const DisjointBoxLayout& grids    = a_res.getBoxes();
    DataIterator             dit      = grids.dataIterator();
    Real                     retVal   = 0.0;
    Real                     boxVal;

    DisjointBoxLayout coarsenedFineGrids;
    ::coarsen(coarsenedFineGrids, a_fineResPtr->getBoxes(), a_refRatio);
    LayoutIterator lit = coarsenedFineGrids.layoutIterator();

    for (dit.reset(); dit.ok(); ++dit) {
        const Box& validBox = grids[dit];

        // Create a copy of resFAB that is zero at covered cells.
        FArrayBox validFAB;
        {
            validFAB.define(validBox, numComps);
            validFAB.copy(a_res[dit]);

            Box overlap;
            for (lit.reset(); lit.ok(); ++lit) {
                overlap = coarsenedFineGrids[lit];
                overlap &= validBox;
                if (!overlap.isEmpty()) {
                    validFAB.setVal(0.0, overlap, 0, numComps);
                }
            }
        }

        // Compute norm of validFAB.
        boxVal = validFAB.norm(validBox, a_p, 0, numComps);
        if (a_p == 0) {
            retVal = max(retVal, boxVal);
        } else {
            retVal += pow(boxVal, a_p);
        }
    }

    if (a_p == 0) {
        Comm::reduce(retVal, MPI_MAX);
    } else {
        Comm::reduce(retVal, MPI_SUM);
        retVal = pow(retVal * powScale, 1.0 / Real(a_p));
    }

    return retVal;
}


// -----------------------------------------------------------------------------
// If L[a_phi] = a divergence of fluxes, then refluxing will be needed.
// This function provides the fluxes needed by the refluxing scheme.
// a_phi's ghosts must be filled prior to call.
// -----------------------------------------------------------------------------
void
PoissonOp::getFlux(FArrayBox&       a_fluxFAB,
                   const FArrayBox& a_phiFAB,
                   const Box&       a_fluxBox,
                   const DataIndex& a_di,
                   const int        a_fluxDir) const
{
    CH_assert(m_activeDirs == IntVect::Unit);

    // Gather references
    const FArrayBox& JgupFAB = m_Jgup[a_di][a_fluxDir];
    const Real       dXiDir  = m_dXi[a_fluxDir];

    // Sanity checks
    CH_assert(a_fluxFAB.nComp() == m_numComps);
    CH_assert(a_phiFAB.nComp() == m_numComps);
    CH_assert(surroundingNodes(m_grids[a_di], a_fluxDir).contains(a_fluxBox));

    // Compute derivs in each direction.
    for (int phiComp = 0; phiComp < m_numComps; ++phiComp) {
        FORT_FINITEDIFF_PARTIALD_CC2NC(CHF_FRA1(a_fluxFAB, phiComp),
                                       CHF_CONST_FRA1(a_phiFAB, phiComp),
                                       CHF_CONST_INT(a_fluxDir),
                                       CHF_CONST_REAL(dXiDir),
                                       CHF_BOX(a_fluxBox));

        a_fluxFAB.mult(JgupFAB, 0, phiComp, 1);
    }  // phiComp

    a_fluxFAB *= m_beta;

    checkForNAN(a_fluxFAB, a_fluxBox);
}


// -----------------------------------------------------------------------------
// Used by the AMR ops with a finer level.
// -----------------------------------------------------------------------------
void
PoissonOp::reflux(LevelData<FArrayBox>& a_res,
                  LevelData<FArrayBox>& a_finePhi,  // Ghost will be filled.
                  const LevelData<FArrayBox>& a_phi,
                  const IntVect&              a_fineRefRatio [[maybe_unused]],
                  const AMRMGOpType&          a_finerOp) const
{
    CH_assert(m_activeDirs == IntVect::Unit);

    // Gather references.
    const DisjointBoxLayout& fineGrids = a_finePhi.getBoxes();
    const Interval&          ivl       = a_phi.interval();

    // Sanity checks. Make sure the levGeo has the AMR, not the MG refRatio.
    CH_assert(coarsen(fineGrids.physDomain(), a_fineRefRatio) == m_domain);
    CH_assert(m_fluxReg.isDefined());
    CH_assert(a_phi.getBoxes() == m_grids);
    CH_assert(a_phi.nComp() == m_numComps);
    CH_assert(a_finePhi.nComp() == m_numComps);

    // 1. Initialize the flux register to zero.
    m_fluxReg.setToZero();

    // 2. Increment the coarse side.
    for (DataIterator ditc(m_grids); ditc.ok(); ++ditc) {
        const FArrayBox& crsePhi = a_phi[ditc];

        if (m_fluxReg.hasCF(ditc())) {
            for (int idir = 0; idir < SpaceDim; ++idir) {
                // Find the destination flux box.
                Box fluxBox = crsePhi.box();
                fluxBox.grow(1);
                fluxBox &= m_grids[ditc];
                fluxBox.surroundingNodes(idir);
                CH_assert(!fluxBox.isEmpty());

                // Calculate the fluxes
                FArrayBox crseFlux(fluxBox, m_numComps);
                this->getFlux(crseFlux, crsePhi, fluxBox, ditc(), idir);

                // Increment the register
                const Real scale = 1.0 / m_dXi[idir];
                m_fluxReg.incrementCoarse(
                    crseFlux, scale, ditc(), ivl, ivl, idir);
            }
        }
    }

    // 3. Increment the fine side.

    // This requires fine ghosts at CFI.
    const PoissonOp& finePoissonOp = static_cast<const PoissonOp&>(a_finerOp);

    CH_assert(finePoissonOp.m_cfInterp.isDefined());
    finePoissonOp.m_cfInterp.interpAtCFI(a_finePhi, a_phi);

    for (DataIterator ditf(fineGrids); ditf.ok(); ++ditf) {
        const FArrayBox& finePhi    = a_finePhi[ditf];
        const Box&       fineRegion = fineGrids[ditf];

        for (int idir = 0; idir < SpaceDim; ++idir) {
            SideIterator sit;
            for (sit.begin(); sit.ok(); sit.next()) {
                if (m_fluxReg.hasCF(ditf(), sit())) {
                    const Side::LoHiSide hiorlo = sit();
                    const Box fluxBox = bdryBox(fineRegion, idir, hiorlo, 1);

                    // Calculate the fluxes
                    FArrayBox fineFlux(fluxBox, m_numComps);
                    finePoissonOp.getFlux(
                        fineFlux, finePhi, fluxBox, ditf(), idir);

                    // Increment the register
                    const Real scale =
                        1.0 / m_dXi[idir];  // Yes, the coarse dXi!
                    m_fluxReg.incrementFine(
                        fineFlux, scale, ditf(), ivl, ivl, idir, hiorlo);
                }  // if has CFI
            }  // sit
        }  // idir
    }  // ditf

    // 4. Reflux.
    nanCheck(a_res);
    m_fluxReg.reflux(a_res);
    nanCheck(a_res);
}


// -----------------------------------------------------------------------------
// Used by the AMR ops with a finer level.
// -----------------------------------------------------------------------------
void
PoissonOp::reflux(LevelData<FArrayBox>&     a_div,
                  const LevelData<FluxBox>& a_flux,
                  const LevelData<FluxBox>& a_fineFlux) const
{
    CH_assert(m_activeDirs == IntVect::Unit);

    // Sanity checks.
    CH_assert(a_div.getBoxes().physDomain() == m_domain);
    CH_assert(a_div.getBoxes() == m_grids);
    CH_assert(m_fluxReg.isDefined());

    // Gather references.
    const DisjointBoxLayout& fineGrids = a_fineFlux.getBoxes();
    const Interval&          ivl       = a_div.interval();

    // 1. Initialize the flux register to zero.
    m_fluxReg.setToZero();

    // 2. Increment the coarse side.
    nanCheck(a_flux);
    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        if (m_fluxReg.hasCF(dit())) {
            for (int idir = 0; idir < SpaceDim; ++idir) {
                const FArrayBox& crseFluxFAB = a_flux[dit][idir];
                const Real       scale       = 1.0 / m_dXi[idir];

                m_fluxReg.incrementCoarse(
                    crseFluxFAB, scale, dit(), ivl, ivl, idir);
            }
        }
    }

    // 3. Increment the fine side.
    nanCheck(a_fineFlux);
    for (DataIterator dit(fineGrids); dit.ok(); ++dit) {
        for (int idir = 0; idir < SpaceDim; ++idir) {
            const FArrayBox& fineFluxFAB = a_fineFlux[dit][idir];
            const Real       scale = 1.0 / m_dXi[idir];  // Yes, the coarse dXi!

            for (SideIterator sit; sit.ok(); sit.next()) {
                if (m_fluxReg.hasCF(dit(), sit())) {
                    m_fluxReg.incrementFine(
                        fineFluxFAB, scale, dit(), ivl, ivl, idir, sit());
                }  // if has CFI
            }  // sit
        }  // idir
    }  // ditf

    // 4. Reflux.
    nanCheck(a_div);
    m_fluxReg.reflux(a_div);
    nanCheck(a_div);
}


// -----------------------------------------------------------------------------
// Set phi BCs and compute Grad[phi].
// Since we need to set CFI BCs, levelGrad = compGrad.
// Both are included for completeness.
// -----------------------------------------------------------------------------
void
PoissonOp::levelGradient(LevelData<FluxBox>&         a_gradPhi,
                         LevelData<FArrayBox>&       a_phi,
                         const LevelData<FArrayBox>* a_crsePhiPtr,
                         const Real                  a_time,
                         const bool                  a_homogPhysBCs,
                         const bool                  a_homogCFIBCs) const
{
    CH_assert(m_activeDirs == IntVect::Unit);

    // Sanity checks
    CH_assert(a_gradPhi.getBoxes().physDomain() == m_domain);
    CH_assert(a_phi.getBoxes().physDomain() == m_domain);
    CH_assert(a_gradPhi.getBoxes().compatible(m_grids));
    CH_assert(a_phi.getBoxes().compatible(m_grids));
    CH_assert(a_gradPhi.nComp() == m_numComps);
    CH_assert(a_phi.nComp() == m_numComps);
    CH_assert(a_phi.ghostVect() >= IntVect::Unit);

    // crse matching
    this->applyBCs(a_phi, a_crsePhiPtr, a_time, a_homogPhysBCs, a_homogCFIBCs);

    // grad
    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        const FArrayBox& phiFAB = a_phi[dit];

        for (int gradDir = 0; gradDir < SpaceDim; ++gradDir) {
            // Gather references
            FArrayBox&       gradFAB = a_gradPhi[dit][gradDir];
            const FArrayBox& JgupFAB = m_Jgup[dit][gradDir];
            const Real       dXiDir  = m_dXi[gradDir];

            // Compute destination box.
            Box fcValid = m_grids[dit];
            fcValid.surroundingNodes(gradDir);
            fcValid &= gradFAB.box();

            // Compute derivs in each direction.
            for (int phiComp = 0; phiComp < m_numComps; ++phiComp) {
                FORT_FINITEDIFF_PARTIALD_CC2NC(CHF_FRA1(gradFAB, phiComp),
                                               CHF_CONST_FRA1(phiFAB, phiComp),
                                               CHF_CONST_INT(gradDir),
                                               CHF_CONST_REAL(dXiDir),
                                               CHF_BOX(fcValid));

                gradFAB.mult(JgupFAB, 0, phiComp, 1);
            }  // phiComp
        }  // gradDir
    }  // dit

    if (RealCmp::neq(m_beta, 1.0)) {
        for (DataIterator dit(m_grids); dit.ok(); ++dit) {
            a_gradPhi[dit] *= m_beta;
        }
    }

    nanCheck(a_gradPhi);
}


// -----------------------------------------------------------------------------
// Set phi BCs and compute the composite Grad[phi].
// -----------------------------------------------------------------------------
void
PoissonOp::compGradient(LevelData<FluxBox>&         a_gradPhi,
                        LevelData<FArrayBox>&       a_phi,
                        const LevelData<FArrayBox>* a_crsePhiPtr,
                        const Real                  a_time,
                        const bool                  a_homogPhysBCs,
                        const bool                  a_homogCFIBCs) const
{
    CH_assert(m_activeDirs == IntVect::Unit);

    // crse matching + grad
    this->levelGradient(
        a_gradPhi, a_phi, a_crsePhiPtr, a_time, a_homogPhysBCs, a_homogCFIBCs);
}


// -----------------------------------------------------------------------------
// Divergence of fluxes without flux corrections from the finer level.
// -----------------------------------------------------------------------------
void
PoissonOp::levelDivergence(LevelData<FArrayBox>&     a_div,
                           const LevelData<FluxBox>& a_flux) const
{
    CH_assert(m_activeDirs == IntVect::Unit);

    // Sanity checks
    CH_assert(a_div.getBoxes().physDomain() == m_domain);
    CH_assert(a_flux.getBoxes().physDomain() == m_domain);
    CH_assert(a_div.getBoxes().compatible(m_grids));
    CH_assert(a_flux.getBoxes().compatible(m_grids));
    CH_assert(a_div.nComp() == a_flux.nComp());

    debugInitLevel(a_div);
    nanCheck(a_flux);

    const RealVect scale = RealVect::Unit;

    // div
    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        FArrayBox&     divFAB = a_div[dit];
        const FluxBox& fluxFB = a_flux[dit];
        const Box&     valid  = m_grids[dit];
        CH_assert(fluxFB.box().contains(valid));

        if (SpaceDim == 2) {
            FORT_FINITEDIFF_DIV2D(CHF_FRA(divFAB),
                                  CHF_CONST_FRA(fluxFB[0]),
                                  CHF_CONST_FRA(fluxFB[1]),
                                  CHF_BOX(valid),
                                  CHF_CONST_REALVECT(m_dXi),
                                  CHF_CONST_REALVECT(scale));
        } else {
            FORT_FINITEDIFF_DIV3D(CHF_FRA(divFAB),
                                  CHF_CONST_FRA(fluxFB[0]),
                                  CHF_CONST_FRA(fluxFB[1]),
                                  CHF_CONST_FRA(fluxFB[2]),
                                  CHF_BOX(valid),
                                  CHF_CONST_REALVECT(m_dXi),
                                  CHF_CONST_REALVECT(scale));
        }
    }  // end loop over grids (dit)

    nanCheck(a_div);
}


// -----------------------------------------------------------------------------
// Composite divergence of fluxes with refluxing.
// -----------------------------------------------------------------------------
void
PoissonOp::compDivergence(LevelData<FArrayBox>&     a_div,
                          const LevelData<FluxBox>& a_flux,
                          const LevelData<FluxBox>* a_fineFluxPtr) const
{
    CH_assert(m_activeDirs == IntVect::Unit);

    // div
    this->levelDivergence(a_div, a_flux);

    // fine matching
    if (a_fineFluxPtr) {
        this->reflux(a_div, a_flux, *a_fineFluxPtr);
    }
}


// ======================== LepticOperator overrides ===========================

// -----------------------------------------------------------------------------
std::shared_ptr<LevelOperator<LevelData<FArrayBox>>>
PoissonOp::createLevelOperator(const DisjointBoxLayout& a_grids) const
{
    return std::shared_ptr<LevelOperator<LevelData<FArrayBox>>>(
        new PoissonOp(*this, a_grids));
}


// -----------------------------------------------------------------------------
std::shared_ptr<MGOperator<LevelData<FArrayBox>>>
PoissonOp::createHorizontalMGOperator(const DisjointBoxLayout& a_hgrids) const
{
    DisjointBoxLayout hCrseGrids;
    if (m_crseAMRGrids.isClosed()) {
        const auto& horizBlockFactor = this->getBlockFactor() * s_hunit + s_vunit;
        const auto& crseDomain = m_crseAMRGrids.physDomain();
        const auto& crseDomBox = crseDomain.domainBox();

        const bool isPeriodic[SpaceDim] = {D_DECL(
            crseDomain.isPeriodic(0),
            crseDomain.isPeriodic(1),
            crseDomain.isPeriodic(2)
        )};

        // Calculate horizontal domain.
        const Box flatCrseDomBox = Subspace::flattenBox(crseDomBox, SpaceDim - 1);
        ProblemDomain hCrseDomain(flatCrseDomBox, isPeriodic);

        // Create the load-balanced horizontal grids.
        Vector<Box> hCrseBoxArray;
        LepticBoxTools::createHorizontalSolverGrids(hCrseBoxArray,
                                                    m_crseAMRGrids.boxArray(),
                                                    crseDomBox,
                                                    horizBlockFactor);
        CH_assert(hCrseBoxArray.size() > 0);
        hCrseGrids.defineAndLoadBalance(hCrseBoxArray, nullptr, hCrseDomain);
    }

    const int hRelaxMethod = ((m_relaxMethod == 6) ? 5 : m_relaxMethod);

    return std::shared_ptr<MGOperator<LevelData<FArrayBox>>>(
        new PoissonOp(m_geoSrc,
                      a_hgrids,
                      hCrseGrids,
                      m_numComps,
                      m_L,
                      m_dXi,
                      m_bcFuncPtr,
                      hRelaxMethod));
}


// =========================== Relaxation methods ==============================

// -----------------------------------------------------------------------------
// Point-relaxation scheme used in Martin's thesis.
// cor^{i+1} = cor^{i} + lambda * (L[cor^{i}] - rhs)
//           = phi^{i} - lambda * res
// lambda = omega * min(dx)^2 / 4
// -----------------------------------------------------------------------------
void
PoissonOp::point_relax(LevelData<FArrayBox>& /*a_cor*/,
                       const LevelData<FArrayBox>& /*a_res*/,
                       const Real /*a_time*/,
                       const int /*a_iters*/) const
{
    MAYDAYERROR("PoissonOp::point_relax is being removed.");
}


// -----------------------------------------------------------------------------
// Jacobi relaxation.
// -----------------------------------------------------------------------------
void
PoissonOp::jacobi_relax(LevelData<FArrayBox>&       a_phi,
                        const LevelData<FArrayBox>& a_rhs,
                        const Real                  a_time,
                        const int                   a_iters) const
{
    LevelData<FArrayBox> res(m_grids, m_numComps);

    for (int iter = 0; iter < a_iters; ++iter) {
        // Compute residual
        this->residual(res, a_phi, nullptr, a_rhs, a_time, true, true);

        // Relax
        for (DataIterator dit(m_grids); dit.ok(); ++dit) {
            FArrayBox& resFAB = res[dit];
            const Box& valid  = m_grids[dit];

            FORT_POISSONOP_JACOBI(CHF_FRA(a_phi[dit]),
                                  CHF_CONST_FRA(resFAB),
                                  CHF_CONST_FRA1(m_Dinv[dit], 0),
                                  CHF_BOX(valid));
        }
    }
}


// -----------------------------------------------------------------------------
// Jacobi relaxation.
// -----------------------------------------------------------------------------
void
PoissonOp::jacobiRB_relax(LevelData<FArrayBox>&       a_phi,
                          const LevelData<FArrayBox>& a_rhs,
                          const Real                  a_time,
                          const int                   a_iters) const
{
    LevelData<FArrayBox> res(m_grids, m_numComps);

    for (int iter = 0; iter < a_iters; ++iter) {
        this->residual(res, a_phi, nullptr, a_rhs, a_time, true, true);

        int whichPass = 0;
        for (DataIterator dit(m_grids); dit.ok(); ++dit) {
            FArrayBox& resFAB = res[dit];
            const Box& valid  = m_grids[dit];

            FORT_POISSONOP_JACOBIRB(CHF_FRA(a_phi[dit]),
                                    CHF_CONST_FRA(resFAB),
                                    CHF_CONST_FRA1(m_Dinv[dit], 0),
                                    CHF_CONST_INT(whichPass),
                                    CHF_BOX(valid));
        }

        a_phi.exchange(m_exCopier);

        whichPass = 1;
        for (DataIterator dit(m_grids); dit.ok(); ++dit) {
            FArrayBox& resFAB = res[dit];
            const Box& valid  = m_grids[dit];

            FORT_POISSONOP_JACOBIRB(CHF_FRA(a_phi[dit]),
                                    CHF_CONST_FRA(resFAB),
                                    CHF_CONST_FRA1(m_Dinv[dit], 0),
                                    CHF_CONST_INT(whichPass),
                                    CHF_BOX(valid));
        }
    }  // iter
}


// -----------------------------------------------------------------------------
// Gauss-Seidel relaxation.
// -----------------------------------------------------------------------------
void
PoissonOp::gs_relax(LevelData<FArrayBox>&       a_phi,
                    const LevelData<FArrayBox>& a_rhs,
                    const Real                  a_time,
                    const int                   a_iters) const
{
    if (m_activeDirs == IntVect::Unit) {
        for (int iter = 0; iter < a_iters; ++iter) {
            this->applyBCs(a_phi, nullptr, a_time, true, true);

            for (DataIterator dit(m_grids); dit.ok(); ++dit) {
                FORT_POISSONOP_GS(CHF_FRA(a_phi[dit]),
                                  CHF_CONST_FRA(a_rhs[dit]),
                                  CHF_CONST_FRA1(m_J[dit], 0),
                                  CHF_CONST_FRA(m_M[0]),
                                  CHF_CONST_FRA(m_M[1]),
                                  CHF_CONST_FRA(m_M[SpaceDim - 1]),
                                  CHF_CONST_FRA1(m_Dinv[dit], 0),
                                  CHF_CONST_REAL(m_beta),
                                  CHF_BOX(m_grids[dit]));
            }  // dit
        }  // iter

    } else if (m_activeDirs == s_hunit) {
        for (int iter = 0; iter < a_iters; ++iter) {
            this->applyBCs(a_phi, nullptr, a_time, true, true);

            for (DataIterator dit(m_grids); dit.ok(); ++dit) {
                FORT_POISSONOP_GS_HORIZ(CHF_FRA(a_phi[dit]),
                                        CHF_CONST_FRA(a_rhs[dit]),
                                        CHF_CONST_FRA1(m_J[dit], 0),
                                        CHF_CONST_FRA(m_M[0]),
                                        CHF_CONST_FRA(m_M[SpaceDim - 2]),
                                        CHF_CONST_FRA1(m_Dinv[dit], 0),
                                        CHF_CONST_REAL(m_beta),
                                        CHF_BOX(m_grids[dit]));
            }  // dit
        }  // iter

    } else {
        MAYDAYERROR(
            "gs_relax only works for activeDirs = all dirs or all dirs but "
            "vertical.");
    }
}


// -----------------------------------------------------------------------------
// Red-Black Gauss-Seidel relaxation.
// -----------------------------------------------------------------------------
void
PoissonOp::gsrb_relax(LevelData<FArrayBox>&       a_phi,
                      const LevelData<FArrayBox>& a_rhs,
                      const Real                  a_time,
                      const int                   a_iters) const
{
    if (m_activeDirs == IntVect::Unit) {
        for (int iter = 0; iter < a_iters; ++iter) {
            this->applyBCs(a_phi,
                           nullptr,
                           a_time,
                           true,
                           true);  // Bottleneck! (due to exchange)

            int whichPass = 0;
            for (DataIterator dit(m_grids); dit.ok(); ++dit) {
                FORT_POISSONOP_GSRB(CHF_FRA(a_phi[dit]),
                                    CHF_CONST_FRA(a_rhs[dit]),
                                    CHF_CONST_FRA1(m_J[dit], 0),
                                    CHF_CONST_FRA(m_M[0]),
                                    CHF_CONST_FRA(m_M[1]),
                                    CHF_CONST_FRA(m_M[SpaceDim - 1]),
                                    CHF_CONST_FRA1(m_Dinv[dit], 0),
                                    CHF_CONST_REAL(m_beta),
                                    CHF_BOX(m_grids[dit]),
                                    CHF_CONST_INT(whichPass));
            }  // dit

            a_phi.exchange(m_exCopier);  // Bottleneck!

            whichPass = 1;
            for (DataIterator dit(m_grids); dit.ok(); ++dit) {
                FORT_POISSONOP_GSRB(CHF_FRA(a_phi[dit]),
                                    CHF_CONST_FRA(a_rhs[dit]),
                                    CHF_CONST_FRA1(m_J[dit], 0),
                                    CHF_CONST_FRA(m_M[0]),
                                    CHF_CONST_FRA(m_M[1]),
                                    CHF_CONST_FRA(m_M[SpaceDim - 1]),
                                    CHF_CONST_FRA1(m_Dinv[dit], 0),
                                    CHF_CONST_REAL(m_beta),
                                    CHF_BOX(m_grids[dit]),
                                    CHF_CONST_INT(whichPass));
            }  // dit
        }  // iter

    } else if (m_activeDirs == s_hunit) {
        for (int iter = 0; iter < a_iters; ++iter) {
            this->applyBCs(a_phi,
                           nullptr,
                           a_time,
                           true,
                           true);  // Bottleneck! (due to exchange)

            int whichPass = 0;
            for (DataIterator dit(m_grids); dit.ok(); ++dit) {
                FORT_POISSONOP_GSRB_HORIZ(CHF_FRA(a_phi[dit]),
                                          CHF_CONST_FRA(a_rhs[dit]),
                                          CHF_CONST_FRA1(m_J[dit], 0),
                                          CHF_CONST_FRA(m_M[0]),
                                          CHF_CONST_FRA(m_M[SpaceDim - 2]),
                                          CHF_CONST_FRA1(m_Dinv[dit], 0),
                                          CHF_CONST_REAL(m_beta),
                                          CHF_BOX(m_grids[dit]),
                                          CHF_CONST_INT(whichPass));
            }  // dit

            a_phi.exchange(m_exCopier);  // Bottleneck!

            whichPass = 1;
            for (DataIterator dit(m_grids); dit.ok(); ++dit) {
                FORT_POISSONOP_GSRB_HORIZ(CHF_FRA(a_phi[dit]),
                                          CHF_CONST_FRA(a_rhs[dit]),
                                          CHF_CONST_FRA1(m_J[dit], 0),
                                          CHF_CONST_FRA(m_M[0]),
                                          CHF_CONST_FRA(m_M[SpaceDim - 2]),
                                          CHF_CONST_FRA1(m_Dinv[dit], 0),
                                          CHF_CONST_REAL(m_beta),
                                          CHF_BOX(m_grids[dit]),
                                          CHF_CONST_INT(whichPass));
            }  // dit
        }  // iter

    } else {
        MAYDAYERROR(
            "gsrb_relax only works for activeDirs = all dirs or all dirs but "
            "vertical.");
    }
}


// -----------------------------------------------------------------------------
// Red-Black Gauss-Seidel line relaxation.
// -----------------------------------------------------------------------------
#if 1 // 1 = async exchange, 0 = sync exchange.
void
PoissonOp::vertLineGSRB_relax(LevelData<FArrayBox>&       a_phi,
                              const LevelData<FArrayBox>& a_rhs,
                              const Real                  a_time,
                              const int                   a_iters) const
{
    CH_assert(m_activeDirs == IntVect::Unit);

    if (a_iters == 0) return;

    CH_assert(m_bcFuncPtr);
    FArrayBox dummyFAB;

    // Create workspace for Fortran
    int Nz = m_domain.size(SpaceDim - 1);
    Box workspace(IntVect::Unit, IntVect(D_DECL(Nz, 1, 1)));
    BaseFab<Real> D(workspace, 1);
    BaseFab<Real> B(workspace, 1);

    workspace.shift(0, 1);
    BaseFab<Real> DU(workspace, 1);  // Last element not used.

    workspace.shift(0, -1);
    BaseFab<Real> DL(workspace, 1);  // Last element not used.

    // We want to preserve the horizontal index to compute
    // the red-black ordering.
    const IntVect validShift = s_vunit * m_domain.domainBox().smallEnd(SpaceDim - 1);

    for (int iter = 0; iter < a_iters; ++iter) {
        for (int whichPass = 0; whichPass < 2; ++whichPass) {
            if (whichPass == 0) {
                this->applyBCs(a_phi,
                               nullptr,
                               a_time,
                               true,
                               true);  // Bottleneck! (due to exchange)
            } else {
                a_phi.exchange(m_exCopier);  // Bottleneck!
            }

            for (DataIterator dit(m_grids); dit.ok(); ++dit) {
                FArrayBox&       phiFAB  = a_phi[dit];
                const FArrayBox& rhsFAB  = a_rhs[dit];
                const FArrayBox& JFAB    = m_J[dit];
                const FArrayBox& DinvFAB = m_Dinv[dit];
                const FArrayBox& loBCFAB = m_vertTriDiagsLoBCs[dit];
                const FArrayBox& hiBCFAB = m_vertTriDiagsHiBCs[dit];
                const Box        valid   = m_grids[dit];

                if constexpr (SpaceDim == 2) {
                    FORT_POISSONOP_VERTLINEGSRB_2D(
                        CHF_FRA1_SHIFT(phiFAB, 0, validShift),
                        CHF_CONST_FRA1_SHIFT(rhsFAB, 0, validShift),
                        CHF_CONST_FRA1_SHIFT(JFAB, 0, validShift),
                        CHF_CONST_FRA(m_M[0]),
                        CHF_CONST_FRA_SHIFT(m_M[SpaceDim - 1], validShift),
                        CHF_CONST_FRA1_SHIFT(DinvFAB, 0, validShift),
                        CHF_CONST_REAL(m_beta),
                        CHF_BOX_SHIFT(valid, validShift),
                        CHF_CONST_INT(whichPass),
                        CHF_CONST_FRA1_SHIFT(loBCFAB, 0, validShift),
                        CHF_CONST_FRA1_SHIFT(hiBCFAB, 0, validShift),
                        CHF_FRA1(DU, 0),
                        CHF_FRA1(D, 0),
                        CHF_FRA1(DL, 0),
                        CHF_FRA1(B, 0));
                } else {
                    FORT_POISSONOP_VERTLINEGSRB_3D(
                        CHF_FRA1_SHIFT(phiFAB, 0, validShift),
                        CHF_CONST_FRA1_SHIFT(rhsFAB, 0, validShift),
                        CHF_CONST_FRA1_SHIFT(JFAB, 0, validShift),
                        CHF_CONST_FRA(m_M[0]),
                        CHF_CONST_FRA(m_M[1]),
                        CHF_CONST_FRA_SHIFT(m_M[SpaceDim - 1], validShift),
                        CHF_CONST_FRA1_SHIFT(DinvFAB, 0, validShift),
                        CHF_CONST_REAL(m_beta),
                        CHF_BOX_SHIFT(valid, validShift),
                        CHF_CONST_INT(whichPass),
                        CHF_CONST_FRA1_SHIFT(loBCFAB, 0, validShift),
                        CHF_CONST_FRA1_SHIFT(hiBCFAB, 0, validShift),
                        CHF_FRA1(DU, 0),
                        CHF_FRA1(D, 0),
                        CHF_FRA1(DL, 0),
                        CHF_FRA1(B, 0));
                }
            }  // dit
        }  // whichPass
    }  // iter
}

#else

void
PoissonOp::vertLineGSRB_relax(LevelData<FArrayBox>&       a_phi,
                              const LevelData<FArrayBox>& a_rhs,
                              const Real                  a_time,
                              const int                   a_iters) const
{
    CH_assert(m_activeDirs == IntVect::Unit);

    if (a_iters == 0) return;

    CH_assert(m_bcFuncPtr);
    FArrayBox dummyFAB;

    // Create workspace for Fortran
    int Nz = m_domain.size(SpaceDim - 1);
    Box workspace(IntVect::Unit, IntVect(D_DECL(Nz, 1, 1)));
    BaseFab<Real> D(workspace, 1);
    BaseFab<Real> B(workspace, 1);

    workspace.shift(0, 1);
    BaseFab<Real> DU(workspace, 1);  // Last element not used.

    workspace.shift(0, -1);
    BaseFab<Real> DL(workspace, 1);  // Last element not used.

    // We want to preserve the horizontal index to compute
    // the red-black ordering.
    IntVect validShift       = IntVect::Zero;
    validShift[SpaceDim - 1] = m_domain.domainBox().smallEnd(SpaceDim - 1);


    for (int iter = 0; iter < a_iters; ++iter) {
        for (int whichPass = 0; whichPass < 2; ++whichPass) {
            if (whichPass == 0) {
                this->applyBCsBegin(a_phi,
                                    nullptr,
                                    a_time,
                                    true,
                                    true);  // Bottleneck! (due to exchange)
            } else {
                a_phi.exchangeBegin(m_exCopier);  // Bottleneck!
            }

            for (DataIterator dit(m_grids); dit.ok(); ++dit) {
                FArrayBox&       phiFAB  = a_phi[dit];
                const FArrayBox& rhsFAB  = a_rhs[dit];
                const FArrayBox& JFAB    = m_J[dit];
                const FArrayBox& DinvFAB = m_Dinv[dit];
                const FArrayBox& loBCFAB = m_vertTriDiagsLoBCs[dit];
                const FArrayBox& hiBCFAB = m_vertTriDiagsHiBCs[dit];
                const Box        valid   = m_grids[dit].grow(-s_hunit);

                if (valid.isEmpty()) continue;

                if constexpr (SpaceDim == 2) {
                    FORT_POISSONOP_VERTLINEGSRB_2D(
                        CHF_FRA1_SHIFT(phiFAB, 0, validShift),
                        CHF_CONST_FRA1_SHIFT(rhsFAB, 0, validShift),
                        CHF_CONST_FRA1_SHIFT(JFAB, 0, validShift),
                        CHF_CONST_FRA(m_M[0]),
                        CHF_CONST_FRA_SHIFT(m_M[SpaceDim - 1], validShift),
                        CHF_CONST_FRA1_SHIFT(DinvFAB, 0, validShift),
                        CHF_CONST_REAL(m_beta),
                        CHF_BOX_SHIFT(valid, validShift),
                        CHF_CONST_INT(whichPass),
                        CHF_CONST_FRA1_SHIFT(loBCFAB, 0, validShift),
                        CHF_CONST_FRA1_SHIFT(hiBCFAB, 0, validShift),
                        CHF_FRA1(DU, 0),
                        CHF_FRA1(D, 0),
                        CHF_FRA1(DL, 0),
                        CHF_FRA1(B, 0));
                } else {
                    FORT_POISSONOP_VERTLINEGSRB_3D(
                        CHF_FRA1_SHIFT(phiFAB, 0, validShift),
                        CHF_CONST_FRA1_SHIFT(rhsFAB, 0, validShift),
                        CHF_CONST_FRA1_SHIFT(JFAB, 0, validShift),
                        CHF_CONST_FRA(m_M[0]),
                        CHF_CONST_FRA(m_M[1]),
                        CHF_CONST_FRA_SHIFT(m_M[SpaceDim - 1], validShift),
                        CHF_CONST_FRA1_SHIFT(DinvFAB, 0, validShift),
                        CHF_CONST_REAL(m_beta),
                        CHF_BOX_SHIFT(valid, validShift),
                        CHF_CONST_INT(whichPass),
                        CHF_CONST_FRA1_SHIFT(loBCFAB, 0, validShift),
                        CHF_CONST_FRA1_SHIFT(hiBCFAB, 0, validShift),
                        CHF_FRA1(DU, 0),
                        CHF_FRA1(D, 0),
                        CHF_FRA1(DL, 0),
                        CHF_FRA1(B, 0));
                }
            }  // dit

            if (whichPass == 0) {
                this->applyBCsEnd(a_phi);  // Bottleneck! (due to exchange)
            } else {
                a_phi.exchangeEnd();  // Bottleneck!
            }

            for (DataIterator dit(m_grids); dit.ok(); ++dit) {
                FArrayBox&       phiFAB  = a_phi[dit];
                const FArrayBox& rhsFAB  = a_rhs[dit];
                const FArrayBox& JFAB    = m_J[dit];
                const FArrayBox& DinvFAB = m_Dinv[dit];
                const FArrayBox& loBCFAB = m_vertTriDiagsLoBCs[dit];
                const FArrayBox& hiBCFAB = m_vertTriDiagsHiBCs[dit];

                Box fullValid = m_grids[dit];
                for (int dir = 0; dir < SpaceDim - 1; ++dir) {
                    for (SideIterator sit; sit.ok(); ++sit) {
                        const Box valid = adjCellBox(fullValid, dir, sit(), -1);
                        if (valid.isEmpty()) continue;

                        if constexpr (SpaceDim == 2) {
                            FORT_POISSONOP_VERTLINEGSRB_2D(
                                CHF_FRA1_SHIFT(phiFAB, 0, validShift),
                                CHF_CONST_FRA1_SHIFT(rhsFAB, 0, validShift),
                                CHF_CONST_FRA1_SHIFT(JFAB, 0, validShift),
                                CHF_CONST_FRA(m_M[0]),
                                CHF_CONST_FRA_SHIFT(m_M[SpaceDim - 1], validShift),
                                CHF_CONST_FRA1_SHIFT(DinvFAB, 0, validShift),
                                CHF_CONST_REAL(m_beta),
                                CHF_BOX_SHIFT(valid, validShift),
                                CHF_CONST_INT(whichPass),
                                CHF_CONST_FRA1_SHIFT(loBCFAB, 0, validShift),
                                CHF_CONST_FRA1_SHIFT(hiBCFAB, 0, validShift),
                                CHF_FRA1(DU, 0),
                                CHF_FRA1(D, 0),
                                CHF_FRA1(DL, 0),
                                CHF_FRA1(B, 0));
                        } else {
                            FORT_POISSONOP_VERTLINEGSRB_3D(
                                CHF_FRA1_SHIFT(phiFAB, 0, validShift),
                                CHF_CONST_FRA1_SHIFT(rhsFAB, 0, validShift),
                                CHF_CONST_FRA1_SHIFT(JFAB, 0, validShift),
                                CHF_CONST_FRA(m_M[0]),
                                CHF_CONST_FRA(m_M[1]),
                                CHF_CONST_FRA_SHIFT(m_M[SpaceDim - 1], validShift),
                                CHF_CONST_FRA1_SHIFT(DinvFAB, 0, validShift),
                                CHF_CONST_REAL(m_beta),
                                CHF_BOX_SHIFT(valid, validShift),
                                CHF_CONST_INT(whichPass),
                                CHF_CONST_FRA1_SHIFT(loBCFAB, 0, validShift),
                                CHF_CONST_FRA1_SHIFT(hiBCFAB, 0, validShift),
                                CHF_FRA1(DU, 0),
                                CHF_FRA1(D, 0),
                                CHF_FRA1(DL, 0),
                                CHF_FRA1(B, 0));
                        }
                    } // sit

                    fullValid.grow(dir, -1);
                } // dir
            }  // dit
        }  // whichPass
    }  // iter
}
#endif


};  // namespace Elliptic
