/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2019
 *    Jefferson University and
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
#ifndef ___BdryIter_H__INCLUDED___
#define ___BdryIter_H__INCLUDED___

#include "DisjointBoxLayout.H"
#include "SideArray.H"
#include "CFRegion.H"
#include "Lockable.H"


struct EmptyStruct {};


// -----------------------------------------------------------------------------
// This data will be available at each physical boundary while iterating.
// See PhysBdryIter.
// -----------------------------------------------------------------------------
struct PhysBdryElement final
{
    DataIndex      di;
    Box            ccValidBox;
    Box            fcBdryBox;
    Box            ccAdjGhostBox;
    int            dir;
    Side::LoHiSide side;
    int            iside;
    int            isign;
    Real           rsign;
};


// -----------------------------------------------------------------------------
// An iterator over the physical boundaries.
// This is not a std::iterator, so don't use it in the std algorithms.
// None of the members are virtual, so I would avoid inheritance.
//
// Typical usage:
//   PhysBdryIter it(grids);
//   for (it.reset(); it.ok(); ++it) {
//       const DataIndex&     di       = it->di;
//       const Box&           ccValid  = it->valid;
//       const int            bdryDir  = it->dir;
//       const Side::LoHiSide side     = it->side;
//       FArrayBox&           stateFAB = state[di];
//
//       // Do something to stateFAB at this boundary...
//   }
// -----------------------------------------------------------------------------
class PhysBdryIter final: public Lockable
{
public:
    PhysBdryIter();

    PhysBdryIter(const PhysBdryIter& a_src)
    : Lockable()
    , m_isDefined(a_src.m_isDefined)
    , m_grids(a_src.m_grids)
    , m_elems(a_src.m_elems)
    , m_numElems(a_src.m_numElems)
    , m_idx(a_src.m_idx)
    , m_isEmpty(a_src.m_isEmpty)
    {}

    PhysBdryIter(const DisjointBoxLayout& a_grids,
                 const SideArray&         a_activeSides = AllSides);

    void
    define(const DisjointBoxLayout& a_grids,
           const SideArray&         a_activeSides = AllSides);

    void
    clear();

    inline bool
    isDefined() const
    {
        return m_isDefined;
    }

    virtual inline bool
    isCompatible(const DisjointBoxLayout& a_grids) const
    {
        return m_isDefined && a_grids.compatible(m_grids);
    }

    // Careful! Rest are not virtual for speed.

    inline bool
    isEmpty() const
    {
        CH_assert(m_isDefined);
        return m_isEmpty;
    }

    inline void
    reset()
    {
        CH_assert(m_isDefined);
        CH_assert(this->isUnlocked());
        m_idx = 0;
    }

    inline void
    reset(const Key& a_key)
    {
        CH_assert(m_isDefined);
        CH_verify(this->accessGranted(a_key));
        m_idx = 0;
    }

    inline bool
    ok() const
    {
        CH_assert(m_isDefined);
        return m_idx < m_numElems;
    }

    inline void
    operator++()
    {
        this->next();
    }

    inline void
    next()
    {
        CH_assert(m_isDefined);
        CH_assert(this->isUnlocked());
        ++m_idx;
    }

    inline void
    next(const Key& a_key)
    {
        CH_assert(m_isDefined);
        CH_verify(this->accessGranted(a_key));
        ++m_idx;
    }

    inline const PhysBdryElement*
    operator->() const
    {
        CH_assert(this->ok());
        return &m_elems[m_idx];
    }

private:
    bool                         m_isDefined;
    DisjointBoxLayout            m_grids;
    std::vector<PhysBdryElement> m_elems;
    unsigned int                 m_numElems;
    unsigned int                 m_idx;
    bool                         m_isEmpty;
};


// -----------------------------------------------------------------------------
// This data will be available at the CFI while iterating.
// See CFIIter.
// -----------------------------------------------------------------------------
struct CFIElement final
{
    CFIElement(const CFIVS& a_cfivs) : cfivs(a_cfivs){};

    DataIndex      di;
    Box            ccValidBox;
    Box            fcBdryBox;
    Box            ccAdjGhostBox; // Overlaps, but not equal to CFI ghosts!!!
    int            dir;
    Side::LoHiSide side;
    int            iside;
    int            isign;
    Real           rsign;
    const CFIVS&   cfivs;
};


// -----------------------------------------------------------------------------
// An iterator over the CFI.
// This is not a std::iterator, so don't use it in the std algorithms.
// None of the members are virtual, so I would avoid inheritance.
// *this exists at the fine-grid level, not the coarse!
//
// Typical usage:
//   CFIIter it(grids, cfRegion);
//   for (it.reset(); it.ok(); ++it) {
//       const DataIndex&     di       = it->di;
//       const Box&           ccValid  = it->ccValidBox;
//       const Box&           cfi      = it->fcBdryBox;
//       const int            bdryDir  = it->dir;
//       const Side::LoHiSide side     = it->side;
//       FArrayBox&           stateFAB = state[di];
//
//       // Do something to stateFAB at this boundary...
//   }
// -----------------------------------------------------------------------------
class CFIIter final: public Lockable
{
public:
    CFIIter();

    CFIIter(const CFIIter& a_src)
    : Lockable()
    , m_isDefined(a_src.m_isDefined)
    , m_grids(a_src.m_grids)
    , m_cfRegion(a_src.m_cfRegion)
    , m_elems(a_src.m_elems)
    , m_numElems(a_src.m_numElems)
    , m_idx(a_src.m_idx)
    , m_isEmpty(a_src.m_isEmpty)
    {}

    CFIIter(const DisjointBoxLayout& a_grids,
            const CFRegion&          a_cfRegion,
            const SideArray&         a_activeSides = AllSides);

    virtual void
    define(const DisjointBoxLayout& a_grids,
           const CFRegion&          a_cfRegion,
           const SideArray&         a_activeSides = AllSides);

    virtual void
    clear();

    virtual inline bool
    isDefined() const
    {
        return m_isDefined;
    }

    virtual inline bool
    isCompatible(const DisjointBoxLayout& a_grids) const
    {
        return m_isDefined && a_grids.compatible(m_grids);
    }

    virtual inline const CFRegion&
    getCFRegion() const
    {
        return m_cfRegion;
    }

    // Careful! Rest are not virtual.

    inline void
    reset()
    {
        CH_assert(m_isDefined);
        CH_assert(this->isUnlocked());
        m_idx = 0;
    }

    inline bool
    isEmpty() const
    {
        CH_assert(m_isDefined);
        return m_isEmpty;
    }

    inline void
    reset(const Key& a_key)
    {
        CH_assert(m_isDefined);
        CH_verify(this->accessGranted(a_key));
        m_idx = 0;
    }

    inline bool
    ok() const
    {
        CH_assert(m_isDefined);
        return m_idx < m_numElems;
    }

    inline void
    operator++()
    {
        this->next();
    }

    inline void
    next()
    {
        CH_assert(m_isDefined);
        CH_assert(this->isUnlocked());
        ++m_idx;
    }

    inline void
    next(const Key& a_key)
    {
        CH_assert(m_isDefined);
        CH_verify(this->accessGranted(a_key));
        ++m_idx;
    }

    inline const CFIElement*
    operator->() const
    {
        CH_assert(this->ok());
        return &m_elems[m_idx];
    }

private:
    bool                    m_isDefined;
    DisjointBoxLayout       m_grids;
    CFRegion                m_cfRegion;
    std::vector<CFIElement> m_elems;
    unsigned int            m_numElems;
    unsigned int            m_idx;
    bool                    m_isEmpty;
};


#endif //!___BdryIter_H__INCLUDED___
