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
 *  https://github.com/somarhub.
 ******************************************************************************/
#include "State.H"
#include "BoxLayoutData.H"
#include "Debug.H"
#include "LoHiSide.H"


// -----------------------------------------------------------------------------
// Basic constructor.
// Sets intervals, does not attach to Q.
// -----------------------------------------------------------------------------
State::State(const int a_numScalars)
: m_velptr(nullptr)
, m_pptr(nullptr)
, m_qptr(nullptr)
//
, numScalars(a_numScalars)
, numQComps(numScalars + 3)
//
, velComp(0)
, velInterval(0, 0)
//
, pComp(0)
, pInterval(0, 0)
//
, TComp(numScalars)
, SComp(numScalars + 1)
, eddyNuComp(numScalars + 2)
//
, scalarsInterval(0, numScalars - 1)
, TInterval(TComp, TComp)
, SInterval(SComp, SComp)
, eddyNuInterval(eddyNuComp, eddyNuComp)
//
, TSInterval(TComp, SComp)
{
    // Do nothing else.
}


// -----------------------------------------------------------------------------
// Copy constructor.
// Sets intervals, a_src and *this will be attached to the same Q.
// -----------------------------------------------------------------------------
State::State(const State& a_src)
: m_velptr(a_src.m_velptr)
, m_pptr(a_src.m_pptr)
, m_qptr(a_src.m_qptr)
//
, numScalars(a_src.numScalars)
, numQComps(numScalars + 3)
//
, velComp(0)
, velInterval(0, 0)
//
, pComp(0)
, pInterval(0, 0)
//
, TComp(numScalars)
, SComp(numScalars + 1)
, eddyNuComp(numScalars + 2)
//
, scalarsInterval(0, numScalars - 1)
, TInterval(TComp, TComp)
, SInterval(SComp, SComp)
, eddyNuInterval(eddyNuComp, eddyNuComp)
//
, TSInterval(TComp, SComp)
{
    this->attachToQ(*a_src.m_velptr, *a_src.m_pptr, *a_src.m_qptr);
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
State::~State()
{
    this->detachFromQ();
}


// -----------------------------------------------------------------------------
// Assignment.
// Intervals already set. Both States will be attached to the same Q.
// -----------------------------------------------------------------------------
const State&
State::operator=(const State& a_src)
{
    this->attachToQ(*a_src.m_velptr, *a_src.m_pptr, *a_src.m_qptr);
    return *this;
}


// -----------------------------------------------------------------------------
// Intervals already set. Calls define on Q, ensuring the correct number
// of comps, then creates aliases to various intervals.
// -----------------------------------------------------------------------------
void
State::defineAndAttachToQ(LevelData<FluxBox>&      a_vel,
                          LevelData<FArrayBox>&    a_p,
                          LevelData<FArrayBox>&    a_q,
                          const DisjointBoxLayout& a_grids,
                          const IntVect&           a_velGhostVect,
                          const IntVect&           a_pGhostVect,
                          const IntVect&           a_qGhostVect)
{
    // I don't want to redefine (reallocate) Q while attached.
    // That sounds like a segfault to me.
    this->detachFromQ();

    a_vel.define(a_grids, 1, a_velGhostVect);
    a_p.define(a_grids, 1, a_pGhostVect);
    a_q.define(a_grids, numQComps, a_qGhostVect);

    this->attachToQ(a_vel, a_p, a_q);
}


// -----------------------------------------------------------------------------
// Intervals already set. This redefines the aliasing LevelDatas,
// grids, and dit.
// -----------------------------------------------------------------------------
void
State::attachToQ(LevelData<FluxBox>&   a_vel,
                 LevelData<FArrayBox>& a_p,
                 LevelData<FArrayBox>& a_q)
{
    this->detachFromQ();

    // This is a mandatory check. We don't just want compatibility.
    CH_verify(a_p.getBoxes() == a_vel.getBoxes());
    CH_verify(a_q.getBoxes() == a_vel.getBoxes());

    // Pointers to user's holders
    m_velptr = &a_vel;
    m_pptr = &a_p;
    m_qptr = &a_q;

    // Layout
    grids = a_vel.getBoxes();
    dit = grids.dataIterator();

    this->velDefineExCopier(velExCopier, a_vel);
    this->pDefineExCopier(pExCopier, a_p);
    this->qDefineExCopier(qExCopier, a_q);

    this->velDefineExCornerCopier1(velExCornerCopier1, a_vel);
    this->velDefineExCornerCopier2(velExCornerCopier2, a_vel);
    this->pDefineExCornerCopier(pExCornerCopier, a_p);
    this->qDefineExCornerCopier(qExCornerCopier, a_q);

    aliasLevelData(vel, m_velptr, m_velptr->interval());
    aliasLevelData(p, m_pptr, m_pptr->interval());
    aliasLevelData(q, m_qptr, m_qptr->interval());
    if (numScalars > 0) {
        aliasLevelData(scalars, m_qptr, scalarsInterval);
    }
    aliasLevelData(T, m_qptr, TInterval);
    aliasLevelData(S, m_qptr, SInterval);
    aliasLevelData(eddyNu, m_qptr, eddyNuInterval);
}

// -----------------------------------------------------------------------------
// Intervals will remain set. Aliasing LevelDatas, grids, and dit undefined.
// -----------------------------------------------------------------------------
void
State::detachFromQ()
{
    grids = DisjointBoxLayout();
    dit = DataIterator();

    D_TERM(
    velExCopier[0].clear();,
    velExCopier[1].clear();,
    velExCopier[2].clear();)
    pExCopier.clear();
    qExCopier.clear();

    D_TERM(
    velExCornerCopier2[0].clear();,
    velExCornerCopier2[1].clear();,
    velExCornerCopier2[2].clear();)
    D_TERM(
    velExCornerCopier1[0].clear();,
    velExCornerCopier1[1].clear();,
    velExCornerCopier1[2].clear();)
    pExCornerCopier.clear();
    qExCornerCopier.clear();

    vel.clear();
    p.clear();
    q.clear();
    scalars.clear();
    T.clear();
    S.clear();
    eddyNu.clear();

    // This MUST be last!
    m_velptr = nullptr;
    m_pptr = nullptr;
    m_qptr = nullptr;
}


// -----------------------------------------------------------------------------
// Checks if this State object is attached to a_vel.
// -----------------------------------------------------------------------------
bool
State::isAttachedTo(const LevelData<FluxBox>& a_vel) const
{
    // Careful! If !m_velptr, then we are not attached to anything, yet
    // nullptr == nullptr would return true. Hence the additional && m_velptr.
    return (&a_vel == m_velptr) && m_velptr;
}


// -----------------------------------------------------------------------------
// Checks if this State object is attached to a_q.
// You can send in a p also.
// -----------------------------------------------------------------------------
bool
State::isAttachedTo(const LevelData<FArrayBox>& a_q) const
{
    // Careful! If !m_qptr, then we are not attached to anything, yet
    // nullptr == nullptr would return true. Hence the additional && m_qptr.
    return ((&a_q == m_qptr) || (&a_q == m_pptr)) && m_qptr;
}


// -----------------------------------------------------------------------------
// Makes an alias to a single component of q.
// -----------------------------------------------------------------------------
void
State::alias(LevelData<FArrayBox>& a_q, const int a_comp)
{
    CH_assert(!a_q.isDefined());
    aliasLevelData(a_q, m_qptr, Interval(a_comp, a_comp));
}


// -----------------------------------------------------------------------------
// Makes an alias to an interval of q.
// -----------------------------------------------------------------------------
void
State::alias(LevelData<FArrayBox>& a_q, const Interval& a_ivl)
{
    CH_assert(!a_q.isDefined());
    aliasLevelData(a_q, m_qptr, a_ivl);
}


// -----------------------------------------------------------------------------
// Define an exchange copier. This ensures we always create the copier in
// the same way. exchangeDefines do different things than regular defines!
// -----------------------------------------------------------------------------
void
State::velDefineExCopier(std::array<StaggeredCopier, CH_SPACEDIM>& a_exCopier,
                         const LevelData<FluxBox>&                 a_vel) const
{
    for (int d = 0; d < SpaceDim; ++d) {
        constexpr bool doValidCorners = true;
        a_exCopier[d].defineValidExchange(
            a_vel.getBoxes(), a_vel.ghostVect(), d, doValidCorners);
    }
}


// -----------------------------------------------------------------------------
// Define an exchange copier. This ensures we always create the copier in
// the same way. exchangeDefines do different things than regular defines!
// -----------------------------------------------------------------------------
void
State::pDefineExCopier(Copier&                     a_exCopier,
                       const LevelData<FArrayBox>& a_p) const
{
    a_exCopier.exchangeDefine(a_p.getBoxes(), a_p.ghostVect());
    a_exCopier.trimEdges(a_p.getBoxes(), a_p.ghostVect());
}


// -----------------------------------------------------------------------------
// Define an exchange copier. This ensures we always create the copier in
// the same way. exchangeDefines do different things than regular defines!
// -----------------------------------------------------------------------------
void
State::qDefineExCopier(Copier&                     a_exCopier,
                       const LevelData<FArrayBox>& a_q) const
{
    a_exCopier.exchangeDefine(a_q.getBoxes(), a_q.ghostVect());
    a_exCopier.trimEdges(a_q.getBoxes(), a_q.ghostVect());
}


// -----------------------------------------------------------------------------
// CornerCopier version.
// -----------------------------------------------------------------------------
void
State::velDefineExCornerCopier1(
    std::array<StaggeredCopier, CH_SPACEDIM>& a_exCopier,
    const LevelData<FluxBox>&                 a_vel) const
{
    for (int d = 0; d < SpaceDim; ++d) {
        a_exCopier[d].defineInvalidCornerExchange1(
            a_vel.getBoxes(), a_vel.ghostVect(), d);
    }
}


// -----------------------------------------------------------------------------
// CornerCopier version.
// -----------------------------------------------------------------------------
void
State::velDefineExCornerCopier2(
    std::array<StaggeredCopier, CH_SPACEDIM>& a_exCopier,
    const LevelData<FluxBox>&                 a_vel) const
{
    for (int d = 0; d < SpaceDim; ++d) {
        a_exCopier[d].defineInvalidCornerExchange2(
            a_vel.getBoxes(), a_vel.ghostVect(), d);
    }
}


// -----------------------------------------------------------------------------
// CornerCopier version.
// -----------------------------------------------------------------------------
void
State::pDefineExCornerCopier(CornerCopier&               a_exCopier,
                             const LevelData<FArrayBox>& a_p) const
{
    // The underscores prevent shadowing member variables.
    const DisjointBoxLayout& _grids     = a_p.getBoxes();
    const ProblemDomain&     _domain    = _grids.physDomain();
    const IntVect&           _ghostVect = a_p.ghostVect();
    a_exCopier.define(_grids, _grids, _domain, _ghostVect, true);
}


// -----------------------------------------------------------------------------
// CornerCopier version.
// -----------------------------------------------------------------------------
void
State::qDefineExCornerCopier(CornerCopier&               a_exCopier,
                             const LevelData<FArrayBox>& a_q) const
{
    // The underscores prevent shadowing member variables.
    const DisjointBoxLayout& _grids     = a_q.getBoxes();
    const ProblemDomain&     _domain    = _grids.physDomain();
    const IntVect&           _ghostVect = a_q.ghostVect();
    a_exCopier.define(_grids, _grids, _domain, _ghostVect, true);
}
