/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2024 Thomas Jefferson University and Arizona State University
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
 *  https://github.com/MUON-CFD/SOMAR.
 ******************************************************************************/
#include "AMRNSLevel.H"
#include "AMRNSLevelF_F.H"
#include "CFInterp.H"
#include "Debug.H"


// -----------------------------------------------------------------------------
// Returns a pointer to a specified level or nullptr if level DNE.
// -----------------------------------------------------------------------------
AMRNSLevel*
AMRNSLevel::getLevel(const int a_level)
{
    CH_assert(a_level >= 0);

    AMRNSLevel* levPtr = this;
    while (levPtr && a_level < levPtr->m_level) {
        levPtr = levPtr->crseNSPtr();
    }
    while (levPtr && a_level > levPtr->m_level) {
        levPtr = levPtr->fineNSPtr();
        if (!levPtr) return nullptr;
    }
    CH_verify(levPtr);
    return levPtr;
}


// -----------------------------------------------------------------------------
// Returns a pointer to a specified level or nullptr if level DNE.
// (const version)
// -----------------------------------------------------------------------------
const AMRNSLevel*
AMRNSLevel::getLevel(const ProblemDomain& a_domain) const
{
    return const_cast<AMRNSLevel*>(this)->getLevel(a_domain);
}


// -----------------------------------------------------------------------------
// Returns a pointer to a specified level or nullptr if level DNE.
// -----------------------------------------------------------------------------
AMRNSLevel*
AMRNSLevel::getLevel(const ProblemDomain& a_domain)
{
    AMRNSLevel* levPtr = coarsestNSPtr();
    while (levPtr) {
        const Box& levDomBox = levPtr->m_problem_domain.domainBox();
        if (levDomBox.size() == a_domain.size()) break;
        levPtr = levPtr->fineNSPtr();
    }
    return levPtr;
}


// -----------------------------------------------------------------------------
// Returns a pointer to a specified level or nullptr if level DNE.
// (const version)
// -----------------------------------------------------------------------------
const AMRNSLevel*
AMRNSLevel::getLevel(const int a_level) const
{
    return const_cast<AMRNSLevel*>(this)->getLevel(a_level);
}


// -----------------------------------------------------------------------------
// Casts the fine level ptr to an AMRNSLevel object.
// -----------------------------------------------------------------------------
AMRNSLevel*
AMRNSLevel::fineNSPtr() const
{
    AMRNSLevel* fine_ns_ptr;
    if (m_finer_level_ptr != nullptr) {
        CH_assert(m_finer_level_ptr != nullptr);

        fine_ns_ptr = dynamic_cast<AMRNSLevel*>(m_finer_level_ptr);
        if (fine_ns_ptr == nullptr) {
            MayDay::Error(
                "AMRNSLevel::fineNSPtr: fineptr not castable to AMRNSLevel*");
        }
    } else {
        fine_ns_ptr = nullptr;
    }

    // Do not return level if it is empty.
    if (fine_ns_ptr) {
        if (fine_ns_ptr->isEmpty()) {
            fine_ns_ptr = nullptr;
        }
    }

    return fine_ns_ptr;
}


// -----------------------------------------------------------------------------
// Casts the coarse level ptr to an AMRNSLevel object.
// -----------------------------------------------------------------------------
AMRNSLevel*
AMRNSLevel::crseNSPtr() const
{
    AMRNSLevel* crse_ns_ptr;
    if (m_coarser_level_ptr != nullptr) {
        CH_assert(m_coarser_level_ptr != nullptr);
        crse_ns_ptr = dynamic_cast<AMRNSLevel*>(m_coarser_level_ptr);
        if (crse_ns_ptr == nullptr) {
            MayDay::Error(
                "AMRNSLevel::crseNSPtr: crseptr not castable to AMRNSLevel*");
        }
    } else {
        crse_ns_ptr = nullptr;
    }

    return crse_ns_ptr;
}


// -----------------------------------------------------------------------------
// Returns a pointer to the finest level
// -----------------------------------------------------------------------------
AMRNSLevel*
AMRNSLevel::finestNSPtr()
{
    AMRNSLevel* levelNSPtr = this;
    while (!levelNSPtr->isFinestLevel()) {
        levelNSPtr = levelNSPtr->fineNSPtr();
    }

    CH_assert(levelNSPtr->isFinestLevel() == true);
    CH_assert(levelNSPtr->isEmpty() == false);

    return levelNSPtr;
}


// -----------------------------------------------------------------------------
// Returns a pointer to the finest level (const version)
// -----------------------------------------------------------------------------
const AMRNSLevel*
AMRNSLevel::finestNSPtr() const
{
    return ((AMRNSLevel*)this)->finestNSPtr();
}


// -----------------------------------------------------------------------------
// Returns a pointer to the coarsest level
// -----------------------------------------------------------------------------
AMRNSLevel*
AMRNSLevel::coarsestNSPtr()
{
    AMRNSLevel* levelNSPtr = this;
    while (levelNSPtr && levelNSPtr->m_level != 0) {
        levelNSPtr = levelNSPtr->crseNSPtr();
    }
    CH_verify(levelNSPtr);
    return levelNSPtr;
}


// -----------------------------------------------------------------------------
// Returns a pointer to the coarsest level (const version)
// -----------------------------------------------------------------------------
const AMRNSLevel*
AMRNSLevel::coarsestNSPtr() const
{
    return ((AMRNSLevel*)this)->coarsestNSPtr();
}


// -----------------------------------------------------------------------------
Real
AMRNSLevel::getMinAMRDx(const int a_coorDir) const
{
    static bool isComputed = false;
    static RealVect minDx(D_DECL(quietNAN, quietNAN, quietNAN));

    if (!isComputed) {
        const ProblemContext*     ctx       = ProblemContext::getInstance();
        const RealVect            L         = ctx->base.L;
        const GeoSourceInterface& geoSrc    = m_levGeoPtr->getGeoSource();
        const Vector<IntVect>&    refRatios = ctx->amr.refRatios;

        ProblemDomain fineDomain = ctx->base.domain;
        for (const auto& r : refRatios) {
            fineDomain.refine(r);
        }

        LevelGeometry fineLevGeo(fineDomain, L, nullptr, &geoSrc);

        for (int d = 0; d < SpaceDim; ++d) {
            minDx[d] = fineLevGeo.getMinCellDx(d);
        }
    }

    return minDx[a_coorDir];
}


// -----------------------------------------------------------------------------
// Is this this finest level?
// -----------------------------------------------------------------------------
bool
AMRNSLevel::isFinestLevel() const
{
    bool retval = true;

    if (!m_finer_level_ptr) return true;

    const AMRNSLevel* fineptr = this->fineNSPtr();
    if (fineptr) {
        retval = fineptr->isEmpty();
    }
    return retval;
}


// -----------------------------------------------------------------------------
// Is this an empty level?
// -----------------------------------------------------------------------------
bool
AMRNSLevel::isEmpty() const
{
    return (m_level_grids.size() == 0);
}


// -----------------------------------------------------------------------------
// Returns this level's DisjointBoxLayout.
// -----------------------------------------------------------------------------
const DisjointBoxLayout&
AMRNSLevel::getBoxes() const
{
    return m_levGeoPtr->getBoxes();
}


// -----------------------------------------------------------------------------
// Returns this level's ProblemDomain.
// -----------------------------------------------------------------------------
const ProblemDomain&
AMRNSLevel::getDomain() const
{
    return m_levGeoPtr->getDomain();
}


// -----------------------------------------------------------------------------
// Returns this level's ProblemDomain as a simple Box.
// -----------------------------------------------------------------------------
const Box&
AMRNSLevel::getDomainBox() const
{
    return m_levGeoPtr->getDomain().domainBox();
}


// -----------------------------------------------------------------------------
// Returns a pointer to the coarser grids. nullptr if this is level 0.
// -----------------------------------------------------------------------------
const DisjointBoxLayout*
AMRNSLevel::getCrseGridsPtr () const
{
    const AMRNSLevel* crsePtr = this->crseNSPtr();
    return (crsePtr ? &(crsePtr->getBoxes()) : nullptr);
}


// -----------------------------------------------------------------------------
// Returns a pointer to the finer grids. nullptr if this is the finest level.
// -----------------------------------------------------------------------------
const DisjointBoxLayout*
AMRNSLevel::getFineGridsPtr () const
{
    const AMRNSLevel* finePtr = this->fineNSPtr();
    return (finePtr ? &(finePtr->getBoxes()) : nullptr);
}


// -----------------------------------------------------------------------------
// Retrieves the State.
// -----------------------------------------------------------------------------
State&
AMRNSLevel::getState()
{
    return *m_statePtr;
}


// -----------------------------------------------------------------------------
// Retrieves the State.
// (const version)
// -----------------------------------------------------------------------------
const State&
AMRNSLevel::getState() const
{
    return *m_statePtr;
}


// -----------------------------------------------------------------------------
// Retrieves the LevelGeometry.
// -----------------------------------------------------------------------------
const LevelGeometry&
AMRNSLevel::getLevGeo() const
{
    return *m_levGeoPtr;
}


// -----------------------------------------------------------------------------
// Packages all LevGeos in the AMR hierarchy into a vector.
// -----------------------------------------------------------------------------
Vector<const LevelGeometry*>
AMRNSLevel::getAMRLevGeos () const
{
    Vector<const LevelGeometry*> amrLG(0);

    const AMRNSLevel* levPtr = this->coarsestNSPtr();
    while (levPtr) {
        amrLG.push_back(&levPtr->getLevGeo());
        levPtr = levPtr->fineNSPtr();
    }

    return amrLG;
}


// -----------------------------------------------------------------------------
// Allocates and defines m_qcc on the indicated levels.
//
// If a_lmax = -1, then we allocate up to the top level.
// All elements outside of the range [a_lbase, a_lmax] will be nullptr.
//
// The resulting vector's indices will correspond to the levels and the
// vector will have as many elements as there are extant levels.
//
// If the levels in the range [a_lbase, a_lmax] are not at a_time, an
// error will be thrown.
//
// It is the caller's responsibility to delete the allocated memory!
// -----------------------------------------------------------------------------
void
AMRNSLevel::allocateAndDefine(Vector<LevelData<FArrayBox>*>& a_amrData,
                              const int                      a_numComps,
                              const IntVect&                 a_ghostVect,
                              int                            a_lbase,
                              int                            a_lmax) const
{
    CH_assert(a_numComps > 0);

    const AMRNSLevel* levPtr;
    int lev;

    // Set default input parameters then make sure they are are well-defined.
    {
        // Find finest level.
        levPtr = this->finestNSPtr();
        lev = levPtr->m_level;

        // Set a_lmax if < 0.
        if (a_lmax < 0) a_lmax = lev;

        // Go the a_lmax and make sure level exists.
        levPtr = this->getLevel(a_lmax);
        lev = a_lmax;
        CH_assert(!levPtr->isEmpty());

        // a_lbase must be set correctly.
        CH_assert(a_lbase >= 0);
        CH_assert(a_lbase <= a_lmax);

        // All input parameters are well-defined for use.
        // levPtr and lev correspond to a_lmax.
    }

    // Now, go from a_lmax to 0 allocating, defining, etc...
    // We begin by setting all data to nullptr.
    this->deallocate(a_amrData);
    a_amrData = Vector<LevelData<FArrayBox>*>(lev + 1, nullptr);

    while (levPtr) {
        // We are sill at an extant level. Get lev.
        lev = levPtr->m_level;

        if (levPtr->m_level >= a_lbase) {
            // We are within [a_lbase, a_lmax].

            // Allocate and define.
            a_amrData[lev] = new LevelData<FArrayBox>();
            a_amrData[lev]->define(levPtr->getBoxes(), a_numComps, a_ghostVect);
        }

        // Move to next level down.
        levPtr = levPtr->crseNSPtr();
    }
}

// -----------------------------------------------------------------------------
// FC version of allocateAndDefine. See it's comments for details.
// -----------------------------------------------------------------------------
void
AMRNSLevel::allocateAndDefine(Vector<LevelData<FluxBox>*>& a_amrData,
                              const int                    a_numComps,
                              const IntVect&               a_ghostVect,
                              int                          a_lbase,
                              int                          a_lmax) const
{
    CH_assert(a_numComps > 0);

    const AMRNSLevel* levPtr;
    int lev;

    // Set default input parameters then make sure they are are well-defined.
    {
        // Find finest level.
        levPtr = this->finestNSPtr();
        lev = levPtr->m_level;

        // Set a_lmax if < 0.
        if (a_lmax < 0) a_lmax = lev;

        // Go the a_lmax and make sure level exists.
        levPtr = this->getLevel(a_lmax);
        lev = a_lmax;
        CH_assert(!levPtr->isEmpty());

        // a_lbase must be set correctly.
        CH_assert(a_lbase >= 0);
        CH_assert(a_lbase <= a_lmax);

        // All input parameters are well-defined for use.
        // levPtr and lev correspond to a_lmax.
    }

    // Now, go from a_lmax to 0 allocating, defining, etc...
    // We begin by setting all data to nullptr.
    this->deallocate(a_amrData);
    a_amrData = Vector<LevelData<FluxBox>*>(lev + 1, nullptr);

    while (levPtr) {
        // We are sill at an extant level. Get lev.
        lev = levPtr->m_level;

        if (levPtr->m_level >= a_lbase) {
            // We are within [a_lbase, a_lmax].

            // Allocate and define.
            a_amrData[lev] = new LevelData<FluxBox>();
            a_amrData[lev]->define(levPtr->getBoxes(), a_numComps, a_ghostVect);
        }

        // Move to next level down.
        levPtr = levPtr->crseNSPtr();
    }
}


// -----------------------------------------------------------------------------
// Allocates and aliases components of m_pPtr on the indicated levels.
//
// If a_ivl = Interval(), then we collect all comps.
// If a_lmax = -1, then we collect up to the top level.
// All elements outside of the range [a_lbase, a_lmax] will be nullptr.
//
// The resulting vector's indices will correspond to the levels and the
// vector will have as many elements as there are extant levels.
//
// If the levels in the range [a_lbase, a_lmax] are not at a_time, an
// error will be thrown.
//
// It is the caller's responsibility to delete the allocated memory!
// -----------------------------------------------------------------------------
void
AMRNSLevel::allocateAndAliasPressure(Vector<LevelData<FArrayBox>*>& a_amrData,
                                     const Real                     a_time,
                                     Interval                       a_ivl,
                                     int                            a_lbase,
                                     int                            a_lmax)
{
    AMRNSLevel* levPtr;
    int lev;

    // Set default input parameters then make sure they are are well-defined.
    {
        // Find finest level.
        levPtr = this->finestNSPtr();
        lev = levPtr->m_level;

        // Set a_lmax if < 0.
        if (a_lmax < 0) a_lmax = lev;

        // Go the a_lmax and make sure level exists.
        levPtr = this->getLevel(a_lmax);
        lev = a_lmax;
        CH_assert(!levPtr->isEmpty());

        // a_lbase must be set correctly.
        CH_assert(a_lbase >= 0);
        CH_assert(a_lbase <= a_lmax);

        // Set a_ivl to all comps if passed in empty.
        if (a_ivl == Interval()) {
            a_ivl = m_pPtr->interval();
        }

        // At this point, a_ivl should be well-defined.
        CH_assert(a_ivl.begin() >= m_pPtr->interval().begin());
        CH_assert(a_ivl.end() <= m_pPtr->interval().end());

        // All input parameters are well-defined for use.
        // levPtr and lev correspond to a_lmax.
    }

    // Now, go from a_lmax to 0 allocating, collecting aliases, etc...
    // We begin by setting all data to nullptr.
    this->deallocate(a_amrData);
    a_amrData = Vector<LevelData<FArrayBox>*>(lev + 1, nullptr);

    while (levPtr) {
        // We are sill at an extant level. Get lev.
        lev = levPtr->m_level;

        if (levPtr->m_level >= a_lbase) {
            // We are within [a_lbase, a_lmax].

            // Make sure levPtr's data is at a_time.
            if (RealCmp::neq(a_time, levPtr->m_time)) {
                MAYDAYERROR("Level " << lev << " at time " << levPtr->m_time
                            << ". Expected " << a_time << ".");
            }

            // Allocate and alias data.
            a_amrData[lev] = new LevelData<FArrayBox>();
            aliasLevelData(*a_amrData[lev], levPtr->m_pPtr, a_ivl);
        }

        // Move to next level down.
        levPtr = levPtr->crseNSPtr();
    }
}


// -----------------------------------------------------------------------------
// m_qPtr version of allocateAndAlias. See it's comments for details.
// -----------------------------------------------------------------------------
void
AMRNSLevel::allocateAndAliasScalars(Vector<LevelData<FArrayBox>*>& a_amrData,
                                    const Real                     a_time,
                                    Interval                       a_ivl,
                                    int                            a_lbase,
                                    int                            a_lmax)
{
    AMRNSLevel* levPtr;
    int lev;

    // Set default input parameters then make sure they are are well-defined.
    {
        // Find finest level.
        levPtr = this->finestNSPtr();
        lev = levPtr->m_level;

        // Set a_lmax if < 0.
        if (a_lmax < 0) a_lmax = lev;

        // Go the a_lmax and make sure level exists.
        levPtr = this->getLevel(a_lmax);
        lev = a_lmax;
        CH_assert(!levPtr->isEmpty());

        // a_lbase must be set correctly.
        CH_assert(a_lbase >= 0);
        CH_assert(a_lbase <= a_lmax);

        // Set a_ivl to all comps if passed in empty.
        if (a_ivl == Interval()) {
            a_ivl = m_qPtr->interval();
        }

        // At this point, a_ivl should be well-defined.
        CH_assert(a_ivl.begin() >= m_qPtr->interval().begin());
        CH_assert(a_ivl.end() <= m_qPtr->interval().end());

        // All input parameters are well-defined for use.
        // levPtr and lev correspond to a_lmax.
    }

    // Now, go from a_lmax to 0 allocating, collecting aliases, etc...
    // We begin by setting all data to nullptr.
    this->deallocate(a_amrData);
    a_amrData = Vector<LevelData<FArrayBox>*>(lev + 1, nullptr);

    while (levPtr) {
        // We are sill at an extant level. Get lev.
        lev = levPtr->m_level;

        if (levPtr->m_level >= a_lbase) {
            // We are within [a_lbase, a_lmax].

            // Make sure levPtr's data is at a_time.
            if (RealCmp::neq(a_time, levPtr->m_time)) {
                MAYDAYERROR("Level " << lev << " at time " << levPtr->m_time
                            << ". Expected " << a_time << ".");
            }

            // Allocate and alias data.
            a_amrData[lev] = new LevelData<FArrayBox>();
            aliasLevelData(*a_amrData[lev], levPtr->m_qPtr, a_ivl);
        }

        // Move to next level down.
        levPtr = levPtr->crseNSPtr();
    }
}


// -----------------------------------------------------------------------------
// m_velPtr of allocateAndAlias. See it's comments for details.
// -----------------------------------------------------------------------------
void
AMRNSLevel::allocateAndAliasVel(Vector<LevelData<FluxBox>*>& a_amrData,
                                const Real                   a_time,
                                Interval                     a_ivl,
                                int                          a_lbase,
                                int                          a_lmax)
{
    AMRNSLevel* levPtr;
    int lev;

    // Set default input parameters then make sure they are are well-defined.
    {
        // Find finest level.
        levPtr = this->finestNSPtr();
        lev = levPtr->m_level;

        // Set a_lmax if < 0.
        if (a_lmax < 0) a_lmax = lev;

        // Go the a_lmax and make sure level exists.
        levPtr = this->getLevel(a_lmax);
        lev = a_lmax;
        CH_assert(!levPtr->isEmpty());

        // a_lbase must be set correctly.
        CH_assert(a_lbase >= 0);
        CH_assert(a_lbase <= a_lmax);

        // Set a_ivl to all comps if passed in empty.
        if (a_ivl == Interval()) {
            a_ivl = m_velPtr->interval();
        }

        // At this point, a_ivl should be well-defined.
        CH_assert(a_ivl.begin() >= m_velPtr->interval().begin());
        CH_assert(a_ivl.end() <= m_velPtr->interval().end());

        // All input parameters are well-defined for use.
        // levPtr and lev correspond to a_lmax.
    }

    // Now, go from a_lmax to 0 allocating, collecting aliases, etc...
    // We begin by setting all data to nullptr.
    this->deallocate(a_amrData);
    a_amrData = Vector<LevelData<FluxBox>*>(lev + 1, nullptr);

    while (levPtr) {
        // We are sill at an extant level. Get lev.
        lev = levPtr->m_level;

        if (levPtr->m_level >= a_lbase) {
            // We are within [a_lbase, a_lmax].

            // Make sure levPtr's data is at a_time.
            if (RealCmp::neq(a_time, levPtr->m_time)) {
                MAYDAYERROR("Level " << lev << " at time " << levPtr->m_time
                            << ". Expected " << a_time << ".");
            }

            // Allocate and alias data.
            a_amrData[lev] = new LevelData<FluxBox>();
            aliasLevelData(*a_amrData[lev], levPtr->m_velPtr, a_ivl);
        }

        // Move to next level down.
        levPtr = levPtr->crseNSPtr();
    }
}


// -----------------------------------------------------------------------------
// Deallocates Vectors created by allocateAndAlias. CC version.
// This would be better as a template that works on any Vector<T*>.
// -----------------------------------------------------------------------------
void
AMRNSLevel::deallocate(Vector<LevelData<FArrayBox>*>& a_amrData) const
{
    for (unsigned int idx = 0; idx < a_amrData.size(); ++idx) {
        delete a_amrData[idx];
        a_amrData[idx] = nullptr;
    }
}


// -----------------------------------------------------------------------------
// Deallocates Vectors created by allocateAndAlias. FC version.
// This would be better as a template that works on any Vector<T*>.
// -----------------------------------------------------------------------------
void
AMRNSLevel::deallocate(Vector<LevelData<FluxBox>*>& a_amrData) const
{
    for (unsigned int idx = 0; idx < a_amrData.size(); ++idx) {
        delete a_amrData[idx];
        a_amrData[idx] = nullptr;
    }
}


// -----------------------------------------------------------------------------
// Average data down. CC version.
// This function starts at the highest, non-null entry and ends at a_lmin
// or the coarsest non-null level.  It is okay to let a_lmin be < 0.
// The indices of a_amrData MUST correspond to the level numbers.
// -----------------------------------------------------------------------------
void
AMRNSLevel::averageDown(Vector<LevelData<FArrayBox>*>& a_amrData,
                        int                            a_lmin,
                        const bool                     a_useJWeighting) const
{
    if (a_lmin < 0) a_lmin = 0;
    const int numLevels = (int)a_amrData.size();

    for (int lev = numLevels - 1; lev > a_lmin; --lev) {
        // lev is the finer level. We average down to lev-1.

        // Is this level null?
        if (a_amrData[lev] == nullptr) continue;

        // Is the coarser level null? If so, we are done.
        if (a_amrData[lev - 1] == nullptr) break;

        // Gather data for the averaging utility.
        const AMRNSLevel*           fineLevPtr    = getLevel(lev);
        const CFInterp&             interpObj     = *fineLevPtr->m_cfInterpPtr;
        const LevelData<FArrayBox>& fineData      = *a_amrData[lev];
        LevelData<FArrayBox>&       crseData      = *a_amrData[lev - 1];
        const bool                  doHarmonicAvg = false;

        const LevelData<FArrayBox>* JPtr = nullptr;
        if (a_useJWeighting) JPtr = &(fineLevPtr->m_levGeoPtr->getCCJ());

        // Sanity checks.
        CH_assert(fineData.getBoxes() == fineLevPtr->getBoxes());
        CH_assert(crseData.getBoxes() == *fineLevPtr->getCrseGridsPtr());

        // Do it!
        interpObj.coarsen(crseData, fineData, doHarmonicAvg, JPtr);
    }
}


// -----------------------------------------------------------------------------
// Average data down. FC version.
// See the CC version's comments for details.
// -----------------------------------------------------------------------------
void
AMRNSLevel::averageDown(Vector<LevelData<FluxBox>*>& a_amrData,
                        int                          a_lmin,
                        const bool                   a_useJWeighting) const
{
    if (a_lmin < 0) a_lmin = 0;

    for (int lev = a_amrData.size() - 1; lev > a_lmin; --lev) {
        // lev is the finer level. We average down to lev-1.

        // Is this level null?
        if (a_amrData[lev] == nullptr) continue;

        // Is the coarser level null? If so, we are done.
        if (a_amrData[lev - 1] == nullptr) break;

        // Gather data for the averaging utility.
        const AMRNSLevel*         fineLevPtr = getLevel(lev);
        const CFInterp&           interpObj  = *fineLevPtr->m_cfInterpPtr;
        const LevelData<FluxBox>& fineData   = *a_amrData[lev];
        LevelData<FluxBox>&       crseData   = *a_amrData[lev - 1];

        const LevelData<FluxBox>* fineJgupPtr = nullptr;
        if (a_useJWeighting) {
            fineJgupPtr = &(fineLevPtr->m_levGeoPtr->getFCJgup());
        }

        // Do it!
        // pout() << "Averaging from level " << lev << " to " << lev - 1 << endl;
        interpObj.coarsen(crseData, fineData, fineJgupPtr);
    }
}


// -----------------------------------------------------------------------------
void
AMRNSLevel::averageDownToThis(const int a_lmax)
{
    // const int lbase = (m_level > 0)? (m_level - 1): 0;
    const int lmin = m_level;
    const int lmax = (a_lmax < 0) ? finestNSPtr()->m_level : a_lmax;

    // Is there any work to do?
    if (m_level == lmax) return;
    CH_assert(lmax > m_level);

    // Velocity
    {
        Vector<LevelData<FluxBox>*> amrVel;
        this->allocateAndAliasVel(
            amrVel, m_time, m_statePtr->velInterval, lmin, lmax);
        this->averageDown(amrVel, lmin, true);
        this->deallocate(amrVel);
    }

    // User-defined scalars
    if (this->numScalars() > 0) {
        Vector<LevelData<FArrayBox>*> amrS;
        this->allocateAndAliasScalars(
            amrS, m_time, m_statePtr->scalarsInterval, lmin, lmax);
        this->averageDown(amrS, lmin, true);
        this->deallocate(amrS);
    }

    // Temperature
    {
        Vector<LevelData<FArrayBox>*> amrT;
        this->allocateAndAliasScalars(
            amrT, m_time, m_statePtr->TInterval, lmin, lmax);
        this->averageDown(amrT, lmin, true);
        this->deallocate(amrT);
    }

    // Salinity
    {
        Vector<LevelData<FArrayBox>*> amrS;
        this->allocateAndAliasScalars(
            amrS, m_time, m_statePtr->SInterval, lmin, lmax);
        this->averageDown(amrS, lmin, true);
        this->deallocate(amrS);
    }
}


// -----------------------------------------------------------------------------
// Set BCs on *m_statePtr on all levels from a_lmax to *this.
// See averageDownToThis for details on a_lmax.
// -----------------------------------------------------------------------------
void
AMRNSLevel::setBCsDownToThis(const int a_lmax)
{
    const int lmax = (a_lmax < 0) ? finestNSPtr()->m_level : a_lmax;
    CH_assert(lmax >= m_level);

    AMRNSLevel* levPtr = this;
    while (levPtr) {
        levPtr->setBC(*levPtr->m_statePtr, m_time);
        if (levPtr->m_level == lmax) break;
        levPtr = levPtr->fineNSPtr();
    };
}
