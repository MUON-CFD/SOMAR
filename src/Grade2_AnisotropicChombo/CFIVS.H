#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CFIVS_H_
#define _CFIVS_H_

#include "REAL.H"
#include "DataIndex.H"
#include "IntVect.H"
#include "Box.H"
#include "IntVectSet.H"
#include "NamespaceHeader.H"
class DisjointBoxLayout;
class ProblemDomain;
class LoHiSide;



/// Internal class to find parts of a box outside the valid region of a level
/**
 *  Typically used to find the ghost cells between this and a coarser level.
 *  These cells are represented by an IntVectSet and, if possible, packed into
 *  a box.  This class should only be used through a CFRegion which breaks down
 *  the parts of a box into regions that have a high probability of being
 *  packed.
 */

class CFIVS
{
    friend void coarsen (CFIVS&         a_cfivs,
                         const IntVect& a_refRatio);

public:
    /// Explicit define based on an IntVectSet.
    /// Using this means the cells to store in the IVS have been identified
    /// elsewhere!
    void
    define(const IntVectSet& a_IVS);

    /// @brief General define for any box using NeighborIterators.
    /// @param a_dataIndex
    /// @param a_grids
    /// @param a_ghostBox The staggering type will be deduced from this.
    ///
    /// @details
    ///  You don't need to worry about valid cells in a_ghostBox. They will be
    ///  detected and removed for you.
    ///
    /// @warning
    ///  To use this, neighbors must be defined in the DBL.
    ///  Otherwise continue to use the legacy define functions.
    void
    define(const DataIndex&         a_dataIndex,
           const DisjointBoxLayout& a_grids,
           const ProblemDomain&     a_domain,
           const Box&               a_ghostBox);

    /// Coarsen the stored cells
    void
    coarsen(const IntVect& a_ref);

    /// Returns packed box for when isPacked() is true.
    inline const Box&
    packedBox() const;


    /// Returns minbox for when isPacked() is false
    /// (but can also be used if true)
    inline const Box&
    minBox() const;

    /// Get the set of IntVects
    inline const IntVectSet&
    getIVS() const;

    /// Returns "true" if this CF IntVectSet can be represented as a simple Box.
    inline bool
    isPacked() const;

    /// Returns "true" if the defined IVS is empty
    inline bool
    isEmpty() const;

    /// Returns "true" if the object has been defined
    inline bool
    isDefined() const;

protected:
    /// Decrement counts (during redefine)
    void
    decrementCounts();

    /// Pack the IVS
    void
    packIVS();

public:
    static long long s_packCount, s_sparseCount;

protected:
    /// IntVects that are outside the level.
    IntVectSet m_IVS = IntVectSet();

    /// If m_packed is true, the Box with which the coarse-fine IntVectSet can
    /// be represented. If m_packed is false, same as m_IVS.minBox().
    Box m_packedBox = Box();

    /// Can the coarse-fine IntVectSet be represented as a Box?
    /// Set to false if the IVS is empty.
    bool m_packed = false;

    /// Is the IVS empty - cached here for fast lookup
    bool m_empty = true;

    /// Has this object been defined?
    bool m_defined = false;
};


// -----------------------------------------------------------------------------
inline const Box&
CFIVS::packedBox() const
{
    CH_assert(this->isPacked());
    return m_packedBox;
}


// -----------------------------------------------------------------------------
inline const Box&
CFIVS::minBox() const
{
    return m_packedBox;
}


// -----------------------------------------------------------------------------
inline const IntVectSet&
CFIVS::getIVS() const
{
    return m_IVS;
}


// -----------------------------------------------------------------------------
inline bool
CFIVS::isPacked() const
{
    CH_assert(m_defined);
    return m_packed;
}


// -----------------------------------------------------------------------------
inline bool
CFIVS::isEmpty() const
{
    CH_assert(m_defined);
    return m_empty;
}


// -----------------------------------------------------------------------------
inline bool
CFIVS::isDefined() const
{
    return m_defined;
}


#include "NamespaceFooter.H"
#endif
