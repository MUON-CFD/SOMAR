/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2019
 *    Jefferson (Philadelphia University + Thomas Jefferson University) and
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
#include "LevelGeometry.H"
#include "LevelGeometryF_F.H"
#include "TensorComp.H"
#include "Format.H"
#include "LayoutTools.H"
#include "AnisotropicRefinementTools.H"
#include "CFInterp.H"
#include "Debug.H"
#include "Subspace.H"
#include <sstream>


// -----------------------------------------------------------------------------
// Full constructor (calls define)
// -----------------------------------------------------------------------------
LevelGeometry::LevelGeometry(const ProblemDomain&      a_domain,
                             const RealVect&           a_domainLength,
                             LevelGeometry*            a_crseLevGeoPtr,
                             const GeoSourceInterface* a_geoSourcePtr)
  : m_geoSourcePtr(nullptr), m_coarserPtr(nullptr), m_finerPtr(nullptr)
{
    this->define(a_domain, a_domainLength, a_crseLevGeoPtr, a_geoSourcePtr);
}


// -----------------------------------------------------------------------------
// Puts this Geometry object in a usable state. Called by the full constructor.
// a_crseLevGeoPtr links this with the coarser level, which sets the refRatios.
// The finer pointer will still be nullptr. It is expected that these defines will
// be called from coarsest to finest.
// If a_geoSourcePtr is nullptr, we will attempt to take a pointer from the coarser
// LevelGeometry.
// -----------------------------------------------------------------------------
void
LevelGeometry::define(const ProblemDomain&      a_domain,
                      const RealVect&           a_domainLength,
                      LevelGeometry*            a_crseLevGeoPtr,
                      const GeoSourceInterface* a_geoSourcePtr)
{
    // Sanity checks
    CH_assert(!a_domain.domainBox().isEmpty());
    D_TERM(
        CH_assert(a_domainLength[0] > 0.0);,
        CH_assert(a_domainLength[1] > 0.0);,
        CH_assert(a_domainLength[2] > 0.0);
    )

    // Just in case the user is redefining this LevelGeometry, we need to
    // remove the data stored over the old grids.
    this->clearMetricCache();

    // Level parameters
    m_domain = a_domain;
    m_domainLength = a_domainLength;
    m_dXi = m_domainLength / RealVect(m_domain.size());
    CH_assert(m_dXi >= RealVect::Zero);

    // Link to coarser LevelGeometry
    m_coarserPtr = nullptr;
    m_finerPtr = nullptr;
    m_crseRefRatio = IntVect::Unit;
    m_fineRefRatio = IntVect::Unit;

    if (a_crseLevGeoPtr != nullptr) {
        // Is the coarser level already linked to a finer level?
        if (a_crseLevGeoPtr->m_finerPtr != nullptr) {
            ERROR(
                "Coarser LevGeo is already linked to a finer LevGeo. Maybe you "
                "are regridding the finer level and mistakenly called "
                "LevelGeometry::define instead of "
                "LevelGeometry::createMetricCache.");
        }

        // Save the coarser object's pointer and compute the ref ratio.
        m_coarserPtr = a_crseLevGeoPtr;
        m_crseRefRatio = m_domain.size() / m_coarserPtr->m_domain.size();

        // Is this domain coarsenable to the coarser domain?
        CH_assert(m_crseRefRatio >= IntVect::Unit);
        CH_assert(m_crseRefRatio.product() >= 1);
        CH_assert(m_crseRefRatio * m_coarserPtr->m_domain.size() == m_domain.size());

        // Inform the coarser LevelGeometry object.
        m_coarserPtr->m_finerPtr = this;
        m_coarserPtr->m_fineRefRatio = m_crseRefRatio;
    }

    // If this level's m_geoSourcePtr is nullptr, take it from the coarser level.
    m_geoSourcePtr = a_geoSourcePtr;

    if ((m_geoSourcePtr == nullptr) && (m_coarserPtr != nullptr)) {
        m_geoSourcePtr = m_coarserPtr->m_geoSourcePtr;
    }
    
    if (m_geoSourcePtr == nullptr) {
        MayDay::Error("LevelGeometry must get a well-defined GeoSourceInterface from a_geoSourcePtr or a_crseLevGeoPtr.");
    }

    // We can cache the physical coordinates before we have grids.
    {
        // The cache should have at least refRatio ghosts for interpolations.
        const IntVect ghostVect(D_DECL(
            max(4, m_crseRefRatio[0]),
            max(4, m_crseRefRatio[1]),
            max(4, m_crseRefRatio[2])
        ));

        for (int dir = 0; dir < SpaceDim; ++dir) {
            Box xBox = m_domain.domainBox();
            xBox.grow(ghostVect);
            xBox = Subspace::flattenBox(xBox, BASISV(dir));

            m_xCellFAB[dir].define(xBox, 1);
            m_geoSourcePtr->fill_physCoor(m_xCellFAB[dir], 0, dir, m_dXi);

            xBox.surroundingNodes(dir);
            m_xNodeFAB[dir].define(xBox, 1);
            m_geoSourcePtr->fill_physCoor(m_xNodeFAB[dir], 0, dir, m_dXi);
        }
    }
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
LevelGeometry::~LevelGeometry ()
{
    // Clean up!
    this->clearMetricCache();

    // Inform the coarser LevelGeometry object.
    if (m_coarserPtr) {
        m_coarserPtr->m_finerPtr = nullptr;
        // Leave the coarser LevGeo's refRatio alone.
    }
}


// -----------------------------------------------------------------------------
// Creates a new version of this levGeo, but coarsened by a_crseRefRatio.
// It is up to you to delete the resulting object when you are through.
// If a_crseRefRatio == (1,1,1), then the cached metric will be aliased.
// The result will NOT be linked to coarser/finer levGeos.
// -----------------------------------------------------------------------------
LevelGeometry*
LevelGeometry::newCoarsenedLG(const IntVect& a_crseRefRatio) const
{
    CH_assert(a_crseRefRatio >= IntVect::Unit);

    LevelGeometry* crseLGPtr = nullptr;

    if (a_crseRefRatio == IntVect::Unit) {
        crseLGPtr = new LevelGeometry(
            m_domain, m_domainLength, nullptr, m_geoSourcePtr);

        if (this->hasGrids()) {
            crseLGPtr->m_hasGrids = true;
            crseLGPtr->m_grids    = this->m_grids;
            crseLGPtr->m_cfRegion = this->m_cfRegion;

            aliasLevelData(crseLGPtr->m_CCJCache,
                           const_cast<LevelData<FArrayBox>*>(&m_CCJCache),
                           m_CCJCache.interval());

            aliasLevelData(crseLGPtr->m_CCJinvCache,
                           const_cast<LevelData<FArrayBox>*>(&m_CCJinvCache),
                           m_CCJinvCache.interval());

            aliasLevelData(crseLGPtr->m_FCJgupCache,
                           const_cast<LevelData<FluxBox>*>(&m_FCJgupCache),
                           m_FCJgupCache.interval());
        }
    } else {
        ProblemDomain crseDomain;
        coarsen(crseDomain, m_domain, a_crseRefRatio);

        crseLGPtr = new LevelGeometry(
            crseDomain, m_domainLength, nullptr, m_geoSourcePtr);

        if (this->hasGrids()) {
            DisjointBoxLayout crseGrids;
            coarsen(crseGrids, m_grids, a_crseRefRatio);

            // Right now, I will not define ghosts. Creating them would be a
            // headache. If we find we need them, we will revisit the issue.
            const IntVect& ghostVect = IntVect::Zero;

            crseLGPtr->m_hasGrids = true;
            crseLGPtr->m_grids    = crseGrids;
            crseLGPtr->m_cfRegion.define(crseGrids, crseDomain);

            crseLGPtr->m_CCJCache.define(crseGrids, 1, ghostVect);
            crseLGPtr->m_CCJinvCache.define(crseGrids, 1, ghostVect);
            crseLGPtr->m_FCJgupCache.define(crseGrids, 1, ghostVect);

            CFInterp interpObj;
            interpObj.define(m_grids, m_dXi, crseGrids);
            interpObj.localCoarsen(crseLGPtr->m_CCJCache, m_CCJCache, false, nullptr);
            interpObj.localCoarsen(crseLGPtr->m_CCJinvCache, m_CCJinvCache, true, nullptr);
            interpObj.localCoarsen(crseLGPtr->m_FCJgupCache, m_FCJgupCache, nullptr);
        }
    }

    return crseLGPtr;
}


// -----------------------------------------------------------------------------
// Evaluates commonly used metric data over the specified grids.
// Before calling this function, this LevelGeometry has no knowledge of this
// level's grids.
// -----------------------------------------------------------------------------
void
LevelGeometry::createMetricCache (const DisjointBoxLayout& a_newGrids)
{
    // Do nothing if the grids haven't changed
    if(m_grids == a_newGrids) return;

    // Clear old metric cache
    this->clearMetricCache();

    // Redefine using the new grids
    m_grids = a_newGrids;
    CH_assert(m_grids.physDomain() == m_domain);

    m_cfRegion.define(m_grids, m_domain);

    // Create cache...
    // The cache should have at least refRatio ghosts for interpolations.
    const IntVect ghostVect(D_DECL(
        max(4, m_crseRefRatio[0]),
        max(4, m_crseRefRatio[1]),
        max(4, m_crseRefRatio[2])
    ));

    m_CCJCache.define(m_grids, 1, ghostVect);
    m_CCJinvCache.define(m_grids, 1, ghostVect);
    m_FCJgupCache.define(m_grids, 1, ghostVect);

    DataIterator dit = m_grids.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        m_geoSourcePtr->fill_J(m_CCJCache[dit], 0, m_dXi);

        m_geoSourcePtr->fill_Jinv(m_CCJinvCache[dit], 0, m_dXi);

        for (int adir = 0; adir < SpaceDim; ++adir) {
            m_geoSourcePtr->fill_Jgup(m_FCJgupCache[dit][adir], 0, adir, m_dXi);
        }
    }

    // the metric cache is ready for use.
    m_hasGrids = true;
}


// -----------------------------------------------------------------------------
// This copies an already established metric cache to a new set of grids.
// The new grids can have a new decomposition and set of processor
// assignments, but can not contain cells in new locations. For that, you
// will need to call createMetricCache() to recompute the metric quantities.
// -----------------------------------------------------------------------------
void
LevelGeometry::redistributeCache (const DisjointBoxLayout& a_newGrids)
{
    // We can only call this function after creating the cache.
    if (!m_hasGrids) {
        ERROR(
            "LevelGeometry::redistributeCache was called before cache was "
            "created");
    }

    // If the new grids are identical to the old grids, scram.
    if (a_newGrids == m_grids)
        return;

    // Copy the cache data to a temp holder (local operation).
    LevelData<FArrayBox> oldJ;
    oldJ.define(m_grids, m_CCJCache.nComp(), m_CCJCache.ghostVect());

    LevelData<FArrayBox> oldJinv;
    oldJinv.define(m_grids, m_CCJinvCache.nComp(), m_CCJinvCache.ghostVect());

    LevelData<FluxBox> oldJgup;
    oldJgup.define(m_grids, m_FCJgupCache.nComp(), m_FCJgupCache.ghostVect());

    for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit) {
        oldJ[dit]   .copy(m_CCJCache[dit]);
        oldJinv[dit].copy(m_CCJinvCache[dit]);
        oldJgup[dit].copy(m_FCJgupCache[dit]);
    }

    // Assign the new grids.
    m_grids = a_newGrids;
    m_cfRegion.define(m_grids, m_domain);

    // Copy the data to the new grids (non local operation).
    {
        Copier cp;
        LayoutTools::defineImpartialCopier(cp,
                                           oldJ.boxLayout(),
                                           m_grids,
                                           m_domain,
                                           oldJ.ghostVect(),
                                           oldJ.ghostVect());
        m_CCJCache.define(m_grids, oldJ.nComp(), oldJ.ghostVect());
        oldJ.copyTo(m_CCJCache, cp);
    }
    {
        Copier cp;
        LayoutTools::defineImpartialCopier(cp,
                                           oldJinv.boxLayout(),
                                           m_grids,
                                           m_domain,
                                           oldJinv.ghostVect(),
                                           oldJinv.ghostVect());
        m_CCJinvCache.define(m_grids, oldJinv.nComp(), oldJinv.ghostVect());
        oldJinv.copyTo(m_CCJinvCache, cp);
    }
    {
        std::array<StaggeredCopier, CH_SPACEDIM> cp;
        LayoutTools::defineImpartialCopier(cp,
                                           oldJgup.boxLayout(),
                                           m_grids,
                                           m_domain,
                                           oldJgup.ghostVect(),
                                           oldJgup.ghostVect());
        m_FCJgupCache.define(m_grids, oldJgup.nComp(), oldJgup.ghostVect());
        oldJgup.copyTo(m_FCJgupCache, cp);
    }
}


// // -----------------------------------------------------------------------------
// // Regrids the current geometry data with flat (vertically averaged) metrics.
// // -----------------------------------------------------------------------------
// void LevelGeometry::regridVertAvg (const DisjointBoxLayout& a_newGrids,
//                                    const Box&               a_fullDomainBox)
// {
//     CH_TIME("LevelGeometry::regridVertAvg");
//     CH_assert(this->isDefined());

//     // Get the new ProblemDomain and make sure it hasn't changed, not
//     // because it will cause problems, but because there is no reason
//     // for that to happpen and usually means something is ill-defined.
//     ProblemDomain domain(a_newGrids.physDomain());
//     CH_assert(m_grids.physDomain().isEmpty() || domain == m_grids.physDomain());

//     // Do nothing if the grids haven't changed
//     if(m_grids == a_newGrids) return;

//     // Clean up the map by removing all fields that still use oldGrids.
//     // It is assumed that these RefCountedPtrs are unique. It may be nice to
//     // check this by including a !isNonUnique() assert, but I'll save that
//     // for the day I begin dereferencing nullptr pointers.
//     s_CCJMap.erase(m_grids);
//     s_CCJinvMap.erase(m_grids);
//     s_FCgupMap.erase(m_grids);
//     s_FCJgupMap.erase(m_grids);
//     s_CCgdnMap.erase(m_grids);

//     // This flushes the entire cache. I could only erase elements that belong to the
//     // regridded levels, but I tried that and it requires a lot of code refactoring.
//     // Besides, this is efficient enough.
//     s_bdryNormMaps.clear();

//     // Redefine using the new grids
//     m_grids = a_newGrids;

//     // Check if the fields are in the field map.
//     // If not, create new ones.
//     m_CCJPtr      = t_CCJPtr(nullptr);
//     m_FCgupPtr    = t_FCgupPtr(nullptr);
//     m_CCgdnPtr    = t_CCgdnPtr(nullptr);
//     m_CCJinvPtr   = this->createVertAvgCCJinvPtr(m_grids, a_fullDomainBox);
//     m_FCJgupPtr   = this->createVertAvgFCJgupPtr(m_grids, a_fullDomainBox);
// }


// -----------------------------------------------------------------------------
// Clears the metric cache and removes all knowledge of this level's grids.
// -----------------------------------------------------------------------------
void
LevelGeometry::clearMetricCache ()
{
    m_hasGrids = false;
    m_grids = DisjointBoxLayout();
    m_cfRegion = CFRegion();
    m_CCJCache.clear();
    m_CCJinvCache.clear();
    m_FCJgupCache.clear();
}



// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
Vector<const LevelData<FArrayBox>*>
LevelGeometry::getAMRCCJ() const
{
    Vector<const LevelData<FArrayBox>*> vptrs(0);
    Vector<const LevelGeometry*> vLevGeoPtrs = getAMRLevGeos();
    for (unsigned int lev = 0; lev < vLevGeoPtrs.size(); ++lev) {
        vptrs.push_back(&vLevGeoPtrs[lev]->getCCJ());
    }
    return vptrs;
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
Vector<const LevelData<FArrayBox>*>
LevelGeometry::getAMRCCJinv() const
{
    Vector<const LevelData<FArrayBox>*> vptrs(0);
    Vector<const LevelGeometry*> vLevGeoPtrs = getAMRLevGeos();
    for (unsigned int lev = 0; lev < vLevGeoPtrs.size(); ++lev) {
        vptrs.push_back(&vLevGeoPtrs[lev]->getCCJinv());
    }
    return vptrs;
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
Vector<const LevelData<FluxBox>*>
LevelGeometry::getAMRFCJgup() const
{
    Vector<const LevelData<FluxBox>*> vptrs(0);
    Vector<const LevelGeometry*> vLevGeoPtrs = getAMRLevGeos();
    for (unsigned int lev = 0; lev < vLevGeoPtrs.size(); ++lev) {
        vptrs.push_back(&vLevGeoPtrs[lev]->getFCJgup());
    }
    return vptrs;
}


// -----------------------------------------------------------------------------
// Fills a FAB with displacements from Xi to physical locations.
// For use with VisIt's displace operator.
// -----------------------------------------------------------------------------
void
LevelGeometry::fill_displacement (FArrayBox& a_dest) const
{
    // Sanity checks
    CH_assert(a_dest.nComp() == SpaceDim);

    // Fill data holder
    this->fill_physCoor(a_dest);

    // Figure out the centering
    const Box destBox(a_dest.box());
    const IntVect boxType(destBox.type());

    // Calculate displacements
    FORT_LEVELGEOMETRY_FILLDISPLACEMENT (
        CHF_FRA(a_dest),
        CHF_CONST_FRA(a_dest),
        CHF_CONST_REALVECT(m_dXi),
        CHF_BOX(destBox),
        CHF_CONST_INTVECT(boxType));
}


// -----------------------------------------------------------------------------
// Fills a mapped box with physical locations
// a_dest must have SpaceDim comps.
// -----------------------------------------------------------------------------
void
LevelGeometry::fill_physCoor (FArrayBox& a_dest) const
{
    // Sanity checks
    CH_assert(a_dest.nComp() == SpaceDim);

    // Fill data holder
    m_geoSourcePtr->fill_physCoor(a_dest, m_dXi);
}


// -----------------------------------------------------------------------------
// Fills a mapped box with physical locations
// Single component version.
// -----------------------------------------------------------------------------
void
LevelGeometry::fill_physCoor(FArrayBox& a_dest,
                             const int  a_destComp,
                             const int  a_coorComp) const
{
    // Sanity checks
    CH_assert(0 <= a_destComp);
    CH_assert(a_destComp < a_dest.nComp());
    CH_assert(0 <= a_coorComp);
    CH_assert(a_coorComp < SpaceDim);

    // Fill data holder
    m_geoSourcePtr->fill_physCoor(a_dest, a_destComp, a_coorComp, m_dXi);
}


// -----------------------------------------------------------------------------
std::vector<Real>
LevelGeometry::getX(const int  a_coorDir,
                    const Box& a_bx) const
{
    size_t       idx = 0;
    IntVect      iv  = a_bx.smallEnd();
    const size_t N   = static_cast<size_t>(a_bx.size(a_coorDir));

    std::vector<Real> vecX(N);

    if (a_bx.type(a_coorDir) == IndexType::CELL) {
        while (idx < N) {
            vecX[idx] = this->getCellX(a_coorDir, iv);
            ++idx;
            ++iv[a_coorDir];
        };

    } else {
        while (idx < N) {
            vecX[idx] = this->getFaceX(a_coorDir, iv, a_coorDir);
            ++idx;
            ++iv[a_coorDir];
        };
    }

    return vecX;
}


// // -----------------------------------------------------------------------------
// // Ensures that metric data is the average of finer data.
// // -----------------------------------------------------------------------------
// void LevelGeometry::averageMetricsDown ()
// {
//     CH_TIME("LevelGeometry::averageMetricsDown");

//     Vector<LevelGeometry*> amrLevGeos = getAMRLevGeos();
//     for (int fineLev = amrLevGeos.size()-1; fineLev > 0; --fineLev) {
//         LevelGeometry* const crseLevGeoPtr = amrLevGeos[fineLev-1];

//         const LevelGeometry* const fineLevGeoPtr = amrLevGeos[fineLev];
//         const DisjointBoxLayout& fineGrids = fineLevGeoPtr->getBoxes();
//         const IntVect& refToCrse = fineLevGeoPtr->getCrseRefRatio();

//         MappedCoarseAverage crseAvgObj;
//         MappedCoarseAverageFace crseAvgFaceObj;
//         int ncomps;

//         { // CC J
//             t_CCJPtr crseMetricPtr = crseLevGeoPtr->m_CCJPtr;
//             const t_CCJPtr fineMetricPtr = fineLevGeoPtr->m_CCJPtr;
//             ncomps = fineMetricPtr->nComp();
//             crseAvgObj.define(fineGrids, ncomps, refToCrse);
//             crseAvgObj.averageToCoarse(*crseMetricPtr, *fineMetricPtr);
//         }

//         { // CC Jinv
//             t_CCJinvPtr crseMetricPtr = crseLevGeoPtr->m_CCJinvPtr;
//             const t_CCJinvPtr fineMetricPtr = fineLevGeoPtr->m_CCJinvPtr;
//             ncomps = fineMetricPtr->nComp();
//             crseAvgObj.define(fineGrids, ncomps, refToCrse);
//             crseAvgObj.averageToCoarseHarmonic(*crseMetricPtr, *fineMetricPtr);
//         }

//         { // CC gdn
//             t_CCgdnPtr crseMetricPtr = crseLevGeoPtr->m_CCgdnPtr;
//             const t_CCgdnPtr fineMetricPtr = fineLevGeoPtr->m_CCgdnPtr;
//             ncomps = fineMetricPtr->nComp();
//             crseAvgObj.define(fineGrids, ncomps, refToCrse);
//             crseAvgObj.averageToCoarse(*crseMetricPtr, *fineMetricPtr);
//         }

//         { // FC gup
//             t_FCgupPtr crseMetricPtr = crseLevGeoPtr->m_FCgupPtr;
//             const t_FCgupPtr fineMetricPtr = fineLevGeoPtr->m_FCgupPtr;
//             ncomps = fineMetricPtr->nComp();
//             crseAvgFaceObj.define(fineGrids, ncomps, refToCrse);
//             crseAvgFaceObj.averageToCoarseHarmonic(*crseMetricPtr, *fineMetricPtr);
//         }

//         { // FC Jgup
//             t_FCJgupPtr crseMetricPtr = crseLevGeoPtr->m_FCJgupPtr;
//             const t_FCJgupPtr fineMetricPtr = fineLevGeoPtr->m_FCJgupPtr;
//             ncomps = fineMetricPtr->nComp();
//             crseAvgFaceObj.define(fineGrids, ncomps, refToCrse);
//             crseAvgFaceObj.averageToCoarseHarmonic(*crseMetricPtr, *fineMetricPtr);
//         }
//     }
// }


// -----------------------------------------------------------------------------
// Returns the coarsest LevGeo in the hierarchy.
// -----------------------------------------------------------------------------
const LevelGeometry*
LevelGeometry::getCoarsestPtr () const
{
    const LevelGeometry* levGeoPtr = this;
    while (levGeoPtr->getCoarserPtr() != nullptr) {
        levGeoPtr = levGeoPtr->getCoarserPtr();
    }

    return levGeoPtr;
}


// -----------------------------------------------------------------------------
// Returns all levgGeos over the AMR hierarchy
// -----------------------------------------------------------------------------
Vector<const LevelGeometry*>
LevelGeometry::getAMRLevGeos () const
{
    // Find the coarsest levGeo
    const LevelGeometry* levGeoPtr = this->getCoarsestPtr();

    // Work up the levels, collecting pointers
    Vector<const LevelGeometry*> vLevGeo(0);
    while (levGeoPtr != nullptr) {
        vLevGeo.push_back(levGeoPtr);
        levGeoPtr = levGeoPtr->getFinerPtr();
    }

    return vLevGeo;
}


// -----------------------------------------------------------------------------
// Returns all fineRefRatios over the AMR hierarchy
// -----------------------------------------------------------------------------
Vector<IntVect>
LevelGeometry::getAMRRefRatios () const
{
    // Find the coarsest levGeo
    const LevelGeometry* levGeoPtr = this->getCoarsestPtr();

    // Work up the levels, collecting crseRefRatios
    Vector<IntVect> vRefRatios(0);
    while (levGeoPtr != nullptr) {
        vRefRatios.push_back(levGeoPtr->getFineRefRatio());
        levGeoPtr = levGeoPtr->getFinerPtr();
    }

    return vRefRatios;
}


// -----------------------------------------------------------------------------
// Returns the LevelGeometry that is at a_domBox's index space or nullptr
// if none found.
// -----------------------------------------------------------------------------
const LevelGeometry*
LevelGeometry::getLevGeoPtr (const Box& a_domBox) const
{
    Vector<const LevelGeometry*> vptrs = getAMRLevGeos();
    for (unsigned int lev = 0; lev < vptrs.size(); ++lev) {
        if (vptrs[lev]->getDomain().domainBox() == a_domBox)
            return vptrs[lev];
    }
    return nullptr;
}



// -----------------------------------------------------------------------------
// Returns all DisjointBoxLayouts over the AMR hierarchy
// -----------------------------------------------------------------------------
Vector<DisjointBoxLayout>
LevelGeometry::getAMRGrids () const
{
    // Find the coarsest levGeo
    const LevelGeometry* levGeoPtr = this->getCoarsestPtr();

    // Work up the levels, collecting crseRefRatios
    Vector<DisjointBoxLayout> vGrids(0);
    while (levGeoPtr != nullptr) {
        vGrids.push_back(levGeoPtr->getBoxes());
        levGeoPtr = levGeoPtr->getFinerPtr();
    }

    return vGrids;
}


// -----------------------------------------------------------------------------
// Returns the number of levels attached through getCoarserPtr/getFinerPtr.
// You can call this from any level.
// -----------------------------------------------------------------------------
int
LevelGeometry::getNumLevels() const
{
    // Find the coarsest levGeo
    const LevelGeometry* levGeoPtr = this->getCoarsestPtr();

    // Count the levels.
    int numLevels = 0;
    while (levGeoPtr != nullptr) {
        ++numLevels;
        levGeoPtr = levGeoPtr->getFinerPtr();
    }

    return numLevels;
}


// -----------------------------------------------------------------------------
// Sends a vector from the mapped basis to the Cartesian basis.
// uCart^{a} = [dx^{a} / dXi^{a}] * uMapped^{a}
// -----------------------------------------------------------------------------
void
LevelGeometry::sendToCartesianBasis (FluxBox& a_vectFlub) const
{
    CH_assert(a_vectFlub.nComp() == 1);

    for (int adir = 0; adir < SpaceDim; ++adir) {
        FArrayBox& vectFAB = a_vectFlub[adir];
        const Box& region = vectFAB.box();

        FArrayBox dxdXiFAB(region, 1);
        m_geoSourcePtr->fill_dxdXi(dxdXiFAB, 0, adir, m_dXi);

        vectFAB.mult(dxdXiFAB, region, 0, 0, 1);
    }
}


// -----------------------------------------------------------------------------
// Sends a vector from the Cartesian basis to the mapped basis.
// uMapped^{a} = [dXi^{a} / dx^{a}] * uCart^{a}
// -----------------------------------------------------------------------------
void
LevelGeometry::sendToMappedBasis (FluxBox& a_vectFlub) const
{
    CH_assert(a_vectFlub.nComp() == 1);

    for (int adir = 0; adir < SpaceDim; ++adir) {
        FArrayBox& vectFAB = a_vectFlub[adir];
        const Box& region = vectFAB.box();

        FArrayBox dXidxFAB(region, 1);
        m_geoSourcePtr->fill_dXidx(dXidxFAB, 0, adir, m_dXi);

        vectFAB.mult(dXidxFAB, region, 0, 0, 1);
    }
}


// -----------------------------------------------------------------------------
// Multiplies every element by J (Entire level version)
// -----------------------------------------------------------------------------
void
LevelGeometry::multByJ (LevelData<FArrayBox>& a_data) const
{
    CH_assert(a_data.getBoxes().physDomain().domainBox() ==
              this->getDomain().domainBox());

    // Just loop over this level's grids and call the single grid version.
    DataIterator dit = a_data.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        this->multByJ(a_data[dit], dit());
    }
}


// -----------------------------------------------------------------------------
// Multiplies every element by J (Single grid version)
// -----------------------------------------------------------------------------
void
LevelGeometry::multByJ (FArrayBox&       a_data,
                        const DataIndex& a_di,
                        int              a_comp) const
{
    CH_assert(m_hasGrids);
    CH_assert(m_grids.check(a_di));
    CH_assert(-1 <= a_comp && a_comp < a_data.nComp());

    // If a_comp is -1, then do all comps.
    int startcomp, numcomp;
    if (a_comp == -1) {
        startcomp = 0;
        numcomp = a_data.nComp();
    } else {
        startcomp = a_comp;
        numcomp = 1;
    }

    const Box& region = a_data.box();
    if (region.type() == IntVect::Zero && m_hasGrids) {
        // We can use the Jinv stored in the cache.
        const FArrayBox& JinvFAB = this->getCCJinv()[a_di];
        CH_assert(JinvFAB.box().contains(region));

        // Scale each component requested
        for (int n = startcomp; n < startcomp + numcomp; ++n) {
            a_data.divide(JinvFAB, 0, n, 1);
        }
    } else {
        // The cache does not support this centering, so we need
        // to fill a custom holder in this region with J.
        FArrayBox JFAB(region, 1);
        m_geoSourcePtr->fill_J(JFAB, 0, m_dXi);

        // Scale each component requested
        for (int n = startcomp; n < startcomp + numcomp; ++n) {
            a_data.mult(JFAB, 0, n, 1);
        }
    }
}


// -----------------------------------------------------------------------------
// Multiplies every element by J (Entire level version)
// -----------------------------------------------------------------------------
void
LevelGeometry::multByJ (LevelData<FluxBox>& a_data) const
{
    CH_assert(a_data.getBoxes().physDomain().domainBox() ==
              this->getDomain().domainBox());

    // Just loop over this level's grids and call the single grid version.
    DataIterator dit = a_data.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        this->multByJ(a_data[dit], dit());
    }
}


// -----------------------------------------------------------------------------
// Multiplies every element by J (Single grid version)
// -----------------------------------------------------------------------------
void
LevelGeometry::multByJ (FluxBox&         a_data,
                        const DataIndex& /*a_di*/,
                        int              a_comp) const
{
    CH_assert(-1 <= a_comp && a_comp < a_data.nComp());

    // If a_comp is -1, then do all comps.
    int startcomp, numcomp;
    if (a_comp == -1) {
        startcomp = 0;
        numcomp = a_data.nComp();
    } else {
        startcomp = a_comp;
        numcomp = 1;
    }

    // Loop over FABs in the FluxBox
    for (int FCdir = 0; FCdir < SpaceDim; ++FCdir) {
        FArrayBox& dataFAB = a_data[FCdir];

        // Fill a holder in this region with J.
        FArrayBox JFAB(dataFAB.box(), 1);
        m_geoSourcePtr->fill_J(JFAB, 0, m_dXi);

        // Scale each component requested
        for (int n = startcomp; n < startcomp + numcomp; ++n) {
            dataFAB.mult(JFAB, 0, n, 1);
        }
    }
}


// -----------------------------------------------------------------------------
// Divides every element by J (Entire level version)
// -----------------------------------------------------------------------------
void
LevelGeometry::divByJ (LevelData<FArrayBox>& a_data) const
{
    CH_assert(a_data.getBoxes().physDomain().domainBox() ==
              this->getDomain().domainBox());

    // Just loop over this level's grids and call the single grid version.
    DataIterator dit = a_data.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        this->divByJ(a_data[dit], dit());
    }
}


// -----------------------------------------------------------------------------
// Divides every element by J (Single grid version)
// -----------------------------------------------------------------------------
void
LevelGeometry::divByJ (FArrayBox&       a_data,
                       const DataIndex& a_di,
                       int              a_comp) const
{
    CH_assert(m_hasGrids);
    CH_assert(m_grids.check(a_di));
    CH_assert(-1 <= a_comp && a_comp < a_data.nComp());

    // If a_comp is -1, then do all comps.
    int startcomp, numcomp;
    if (a_comp == -1) {
        startcomp = 0;
        numcomp = a_data.nComp();
    } else {
        startcomp = a_comp;
        numcomp = 1;
    }

    const Box& region = a_data.box();
    if (region.type() == IntVect::Zero) {
        // We can use the Jinv stored in the cache.
        const FArrayBox& JinvFAB = this->getCCJinv()[a_di];
        CH_assert(JinvFAB.box().contains(region));

        // Scale each component requested
        for (int n = startcomp; n < startcomp + numcomp; ++n) {
            a_data.mult(JinvFAB, 0, n, 1);
        }
    } else {
        // The cache does not support this centering, so we need
        // to fill a custom holder in this region with Jinv.
        FArrayBox JinvFAB(region, 1);
        m_geoSourcePtr->fill_Jinv(JinvFAB, 0, m_dXi);

        // Scale each component requested
        for (int n = startcomp; n < startcomp + numcomp; ++n) {
            a_data.mult(JinvFAB, 0, n, 1);
        }
    }
}


// -----------------------------------------------------------------------------
// Divides every element by J (Entire level version)
// -----------------------------------------------------------------------------
void
LevelGeometry::divByJ (LevelData<FluxBox>& a_data) const
{
    CH_assert(a_data.getBoxes().physDomain().domainBox() ==
              this->getDomain().domainBox());

    // Just loop over this level's grids and call the single grid version.
    DataIterator dit = a_data.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        this->divByJ(a_data[dit], dit());
    }
}


// -----------------------------------------------------------------------------
// Divides every element by J (Single grid version)
// -----------------------------------------------------------------------------
void
LevelGeometry::divByJ (FluxBox&         a_data,
                       const DataIndex& /*a_di*/,
                       int              a_comp) const
{
    CH_assert(-1 <= a_comp && a_comp < a_data.nComp());

    // If a_comp is -1, then do all comps.
    int startcomp, numcomp;
    if (a_comp == -1) {
        startcomp = 0;
        numcomp = a_data.nComp();
    } else {
        startcomp = a_comp;
        numcomp = 1;
    }

    // Loop over FABs in the FluxBox
    for (int FCdir = 0; FCdir < SpaceDim; ++FCdir) {
        FArrayBox& dataFAB = a_data[FCdir];

        // Fill a holder in this region with Jinv.
        FArrayBox JinvFAB(dataFAB.box(), 1);
        m_geoSourcePtr->fill_Jinv(JinvFAB, 0, m_dXi);

        // Scale each component requested
        for (int n = startcomp; n < startcomp + numcomp; ++n) {
            dataFAB.mult(JinvFAB, 0, n, 1);
        }
    }
}


// -----------------------------------------------------------------------------
// Contracts a CC contra-vector with gdn, making it covariant.
// Single comp, single grid version.
// -----------------------------------------------------------------------------
void
LevelGeometry::makeCovariant(FArrayBox&       a_coVect,
                             const int        a_coVectComp,
                             const FArrayBox& a_contraVect,
                             const int        a_contraVectComp,
                             const Box&       a_region,
                             const DataIndex& /*a_di*/,
                             const int        a_mu) const
{
    CH_assert(0 <= a_coVectComp);
    CH_assert(a_coVectComp < a_coVect.nComp());

    CH_assert(0 <= a_contraVectComp);
    CH_assert(a_contraVectComp < a_contraVect.nComp());

    CH_assert(0 <= a_mu);
    CH_assert(a_mu < SpaceDim);

    CH_assert(a_coVect.dataPtr() != a_contraVect.dataPtr());

    CH_assert(a_coVect.box().type() == a_region.type());
    CH_assert(a_coVect.box().contains(a_region));

    CH_assert(a_contraVect.box().type() == a_region.type());
    CH_assert(a_contraVect.box().contains(a_region));

    m_geoSourcePtr->fill_gdn(a_coVect, a_coVectComp, a_mu, m_dXi);
    a_coVect.mult(a_contraVect, a_region, a_contraVectComp, a_coVectComp, 1);
}


// -----------------------------------------------------------------------------
// Contracts a CC contra-vector with gdn, making it covariant.
// (Single-grid version)
// -----------------------------------------------------------------------------
void
LevelGeometry::makeCovariant (FArrayBox&       a_coVect,
                              const FArrayBox& a_contraVect,
                              const Box&       a_region,
                              const DataIndex& a_di) const
{
    CH_assert(a_coVect.nComp() == SpaceDim);
    CH_assert(a_contraVect.nComp() == SpaceDim);

    for (int mu = 0; mu < SpaceDim; ++mu) {
        this->makeCovariant(a_coVect, mu, a_contraVect, mu, a_region, a_di, mu);
    }
}


// -----------------------------------------------------------------------------
// Contracts a CC covariant one-form with gup, making it a contravariant vector.
// Single comp, single grid version.
// -----------------------------------------------------------------------------
void
LevelGeometry::makeContravariant(FArrayBox&       a_contraVect,
                                 const int        a_contraVectComp,
                                 const FArrayBox& a_coVect,
                                 const int        a_coVectComp,
                                 const Box&       a_region,
                                 const DataIndex& /*a_di*/,
                                 const int        a_mu) const
{
    CH_assert(0 <= a_coVectComp);
    CH_assert(a_coVectComp < a_coVect.nComp());

    CH_assert(0 <= a_contraVectComp);
    CH_assert(a_contraVectComp < a_contraVect.nComp());

    CH_assert(0 <= a_mu);
    CH_assert(a_mu < SpaceDim);

    CH_assert(a_coVect.dataPtr() != a_contraVect.dataPtr());

    CH_assert(a_coVect.box().type() == a_region.type());
    CH_assert(a_coVect.box().contains(a_region));

    CH_assert(a_contraVect.box().type() == a_region.type());
    CH_assert(a_contraVect.box().contains(a_region));

    m_geoSourcePtr->fill_gup(a_contraVect, a_contraVectComp, a_mu, m_dXi);
    a_contraVect.mult(a_coVect, a_region, a_coVectComp, a_contraVectComp, 1);
}


// -----------------------------------------------------------------------------
// Contracts a CC covariant one-form with gup, making it a contravariant vector.
// (Single-grid version)
// -----------------------------------------------------------------------------
void
LevelGeometry::makeContravariant (FArrayBox&       a_contraVect,
                                  const FArrayBox& a_coVect,
                                  const Box&       a_region,
                                  const int        a_vectComp,
                                  const DataIndex& a_di) const
{
    CH_assert(a_coVect.nComp() == SpaceDim);
    CH_assert(a_contraVect.nComp() == SpaceDim);

    for (int mu = 0; mu < SpaceDim; ++mu) {
        this->makeContravariant(a_contraVect, mu, a_coVect, mu, a_region, a_di, mu);
    }
}


// -----------------------------------------------------------------------------
// Compute the magnitude of a vector
// (Single-grid version)
// -----------------------------------------------------------------------------
void
LevelGeometry::contractVectors (FArrayBox&       a_res,
                                const FArrayBox& a_vec1,
                                const FArrayBox& a_vec2,
                                const DataIndex& /*a_di*/,
                                const bool       a_isCartesian) const
{
    // Sanity checks
    CH_assert(a_res.nComp() == 1);
    CH_assert(a_vec1.nComp() == SpaceDim);
    CH_assert(a_vec2.nComp() == SpaceDim);

    const Box& region = a_res.box();
    CH_assert(a_vec1.box().type() == region.type());
    CH_assert(a_vec2.box().type() == region.type());

    CH_assert(a_vec1.box().contains(region));
    CH_assert(a_vec2.box().contains(region));

    if (a_isCartesian) {
        // Compute
        FORT_LEVELGEOMETRY_CONTRACTCARTESIANVECTORS (
            CHF_FRA1(a_res,0),
            CHF_CONST_FRA(a_vec1),
            CHF_CONST_FRA(a_vec2),
            CHF_BOX(region));

    } else {
        a_res.setVal(0.0);

        FArrayBox gdn(region, 1);
        for (int mu = 0; mu < SpaceDim; ++mu) {
            m_geoSourcePtr->fill_gdn(gdn, 0, mu, m_dXi);
            gdn.mult(a_vec1);
            gdn.mult(a_vec2);
            a_res.plus(gdn);
        }
    }
}


// -----------------------------------------------------------------------------
Real
LevelGeometry::getMinCellDx(const int a_coorDir) const
{
    const FArrayBox& xFAB = m_xNodeFAB[a_coorDir];
    const int        imax = xFAB.box().bigEnd(a_coorDir);
    const IntVect    e    = BASISV(a_coorDir);

    IntVect iv = xFAB.box().smallEnd();
    Real minDx = xFAB(iv + e) - xFAB(iv);

    ++iv[a_coorDir];
    while (iv[a_coorDir] < imax) {
        minDx = std::min(minDx, xFAB(iv + e) - xFAB(iv));
        ++iv[a_coorDir];
    }

    return minDx;
}


// -----------------------------------------------------------------------------
Real
LevelGeometry::getMaxCellDx(const int a_coorDir) const
{
    const FArrayBox& xFAB = m_xNodeFAB[a_coorDir];
    const int        imax = xFAB.box().bigEnd(a_coorDir);
    const IntVect    e    = BASISV(a_coorDir);

    IntVect iv = xFAB.box().smallEnd();
    Real maxDx = xFAB(iv + e) - xFAB(iv);

    ++iv[a_coorDir];
    while (iv[a_coorDir] < imax) {
        maxDx = std::max(maxDx, xFAB(iv + e) - xFAB(iv));
        ++iv[a_coorDir];
    }

    return maxDx;
}


// -----------------------------------------------------------------------------
int
LevelGeometry::findContainingCellIdxDir(const Real a_X, const int a_dir) const
{
    const Box& nodeBox = m_xNodeFAB[a_dir].box();
    IntVect ivmin = nodeBox.smallEnd();
    IntVect ivmax = nodeBox.bigEnd();

    if (m_xNodeFAB[a_dir](ivmin) <= a_X && a_X <= m_xNodeFAB[a_dir](ivmax)) {
        // X is within bounds. Perform binary search.
        IntVect iv(D_DECL(0, 0, 0));
        while (ivmax[a_dir] - ivmin[a_dir] > 1) {
            iv[a_dir] = (ivmin[a_dir] + ivmax[a_dir]) / 2;

            if (a_X < m_xNodeFAB[a_dir](iv)) {
                ivmax[a_dir] = iv[a_dir];
            } else {
                ivmin[a_dir] = iv[a_dir];
            }
        }
        return ivmin[a_dir];  // The left FC index is also the CC index.

    } else {
        // X is out of bounds. Compute estimate.
        return std::floor(a_X / m_dXi[a_dir]);
    }
}


// -----------------------------------------------------------------------------
IntVect
LevelGeometry::findContainingCellIdx(const RealVect& a_X) const
{
    IntVect retVal;

    for (int d = 0; d < SpaceDim; ++d) {
        const Box& nodeBox = m_xNodeFAB[d].box();
        IntVect    ivmin   = nodeBox.smallEnd();
        IntVect    ivmax   = nodeBox.bigEnd();

        if (m_xNodeFAB[d](ivmin) <= a_X[d] && a_X[d] <= m_xNodeFAB[d](ivmax)) {
            // X is within bounds. Perform binary search.
            IntVect iv(D_DECL(0, 0, 0));
            while (ivmax[d] - ivmin[d] > 1) {
                iv[d] = (ivmin[d] + ivmax[d]) / 2;

                if (a_X[d] < m_xNodeFAB[d](iv)) {
                    ivmax[d] = iv[d];
                } else {
                    ivmin[d] = iv[d];
                }
            }
            retVal[d] = ivmin[d];  // The left FC index is also the CC index.

            // Double check the iv.
            CH_assert(getFaceX(d, retVal, d) <= a_X[d]);
            CH_assert(a_X[d] <= getFaceX(d, retVal + BASISV(d), d));

        } else {
            // X is out of bounds. Compute estimate.
            retVal[d] = std::floor(a_X[d] / m_dXi[d]);
            // MAYDAYERROR("Out of bounds.");
        }
    }

    return retVal;
}


#ifndef NDEBUG
// -----------------------------------------------------------------------------
// Dumps debugging info to pout()
// -----------------------------------------------------------------------------
void LevelGeometry::dump () const
{
    pout() << "LevelGeometry dump:\n"
           << "\tm_finerPtr   = " << (void*)m_finerPtr << "\n"
           << "\tthis         = " << (void*)this << "\n"
           << "\tm_coarserPtr = " << (void*)m_coarserPtr << "\n"
           << "\tCoordinate map = " << getCoorMapName() << "\n"
           << "\tis uniform  = " << (isUniform()? "true": "false") << "\n"
           << "\tm_dXi            = " << m_dXi << "\n"
           << "\tm_fineRefRatio   = " << m_fineRefRatio << "\n"
           << "\tm_crseRefRatio   = " << m_crseRefRatio << "\n"
           << "\tm_domain = " << m_domain << "\n"
           << "\tgrids.physDomain() = " << m_grids.physDomain() << "\n"
           << "\tgrids = " << m_grids << std::flush;
}
#endif

