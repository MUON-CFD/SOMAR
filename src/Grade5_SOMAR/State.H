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
#ifndef ___State_H__INCLUDED___
#define ___State_H__INCLUDED___

#include "LevelData.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "CornerCopier.H"
#include "StaggeredCopier.H"


// A wrapper for AMRNSLevel's state variables.
// The FC velocity and CC pressure are special, so they are excluded
// from the q LevelData. Everything else is within q -- user-defined
// scalars, temperature, salinity, and eddyNu. I will use a
// capital Q when referring to the entire set of variables Q = [vel, p, q].
class State
{
public:
    // Basic constructor.
    // Sets intervals, does not attach to Q.
    State(const int a_numScalars);

    // Copy constructor.
    // Sets intervals, both States will be attached to the same Q.
    State(const State& a_src);

    // Destructor
    ~State();

    // Assignment.
    // Intervals already set. Both States will be attached to the same Q.
    const State&
    operator=(const State& a_src);

    // Intervals already set. Calls define on a_q, ensuring the correct number
    // of comps, then creates aliases to various intervals.
    void
    defineAndAttachToQ(LevelData<FluxBox>&      a_vel,
                       LevelData<FArrayBox>&    a_p,
                       LevelData<FArrayBox>&    a_q,
                       const DisjointBoxLayout& a_grids,
                       const IntVect&           a_velGhostVect,
                       const IntVect&           a_pGhostVect,
                       const IntVect&           a_qGhostVect);

    // Intervals already set. This redefines the aliasing LevelDatas,
    // grids, and dit.
    void
    attachToQ (LevelData<FluxBox>&   a_vel,
               LevelData<FArrayBox>& a_p,
               LevelData<FArrayBox>& a_q);

    // Intervals will remain set. Aliasing LevelDatas, grids, and dit undefined.
    void
    detachFromQ ();

    // Checks if this State object is attached to a_vel.
    bool
    isAttachedTo (const LevelData<FluxBox>& a_vel) const;

    // Checks if this State object is attached to a_q.
    // You can send in a p also.
    bool
    isAttachedTo (const LevelData<FArrayBox>& a_q) const;

    // Makes an alias to a single component of q.
    void
    alias(LevelData<FArrayBox>& a_q, const int a_comp);

    // Makes an alias to an interval of q.
    void
    alias(LevelData<FArrayBox>& a_q, const Interval& a_ivl);

    // Define an exchange copier. This ensures we always create the copier in
    // the same way. exchangeDefines do different things than regular defines!
    static void
    velDefineExCopier(std::array<StaggeredCopier, CH_SPACEDIM>& a_exCopier,
                      const LevelData<FluxBox>&                 a_vel);

    // Define an exchange copier. This ensures we always create the copier in
    // the same way. exchangeDefines do different things than regular defines!
    static void
    pDefineExCopier(Copier& a_exCopier, const LevelData<FArrayBox>& a_p);

    // Define an exchange copier. This ensures we always create the copier in
    // the same way. exchangeDefines do different things than regular defines!
    static void
    qDefineExCopier(Copier& a_exCopier, const LevelData<FArrayBox>& a_q);

    // CornerCopier version.
    static void
    velDefineExCornerCopier1(
        std::array<StaggeredCopier, CH_SPACEDIM>& a_exCopier,
        const LevelData<FluxBox>&                 a_vel);

    // CornerCopier version.
    static void
    velDefineExCornerCopier2(
        std::array<StaggeredCopier, CH_SPACEDIM>& a_exCopier,
        const LevelData<FluxBox>&                 a_vel);

    // CornerCopier version.
    static void
    pDefineExCornerCopier(CornerCopier&               a_exCopier,
                          const LevelData<FArrayBox>& a_p);

    // CornerCopier version.
    static void
    qDefineExCornerCopier(CornerCopier&               a_exCopier,
                          const LevelData<FArrayBox>& a_q);

    // The state data...
    // We do not own this data, AMRNSLevel will handle memory management.
    LevelData<FluxBox>*   m_velptr;
    LevelData<FArrayBox>* m_pptr;
    LevelData<FArrayBox>* m_qptr;

    DisjointBoxLayout                        grids;
    DataIterator                             dit;
    std::array<StaggeredCopier, CH_SPACEDIM> velExCopier;
    Copier                                   pExCopier;
    Copier                                   qExCopier;
    std::array<StaggeredCopier, CH_SPACEDIM> velExCornerCopier1;
    std::array<StaggeredCopier, CH_SPACEDIM> velExCornerCopier2;
    CornerCopier                             pExCornerCopier;
    CornerCopier                             qExCornerCopier;

    // q interval assignments and aliases...
    const int numScalars;
    const int numQComps; // This includes u, (v in 3D), w, p, and q.

    // vel comps and intervals
    const int velComp;
    const Interval velInterval;

    // p comps and intervals
    const int pComp;
    const Interval pInterval;

    // q comps and intervals
    const int TComp;
    const int SComp;

    const Interval scalarsInterval;
    const Interval TInterval;
    const Interval SInterval;

    const Interval TSInterval;   // Composite

    LevelData<FluxBox>   vel;
    LevelData<FArrayBox> p;
    LevelData<FArrayBox> q;
    LevelData<FArrayBox> scalars;
    LevelData<FArrayBox> T;
    LevelData<FArrayBox> S;
};


#endif //!___State_H__INCLUDED___
