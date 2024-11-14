#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CFREGION_H_
#define _CFREGION_H_

#include <array>
#include "IntVect.H"
#include "LayoutData.H"
#include "CFIVS.H"
#include "NamespaceHeader.H"
class Box;
class ProblemDomain;
class DisjointBoxLayout;


/// \brief Class that represents the edge region around a DisjointBoxLayout
/// \details The coarsen functions have been removed. They only work with CC IVs.
class CFRegion final
{
    friend void
    coarsen(CFRegion& a_cfRegion, const IntVect& a_refRatio);

public:
    // Rule of five stuff. Nontrivial because LayoutData<T> is not copyable.
    CFRegion()                            = default;
    ~CFRegion()                           = default;
    CFRegion(CFRegion&& a_src) noexcept   = default;
    CFRegion& operator=(CFRegion&& a_src) = default;
    CFRegion(const CFRegion& a_src) {*this = a_src;};
    CFRegion& operator=(const CFRegion& a_src);

    /// @brief Defining constructor.  Constructs a valid object.
    //    Equivalent to default construction followed by define().
    /// @param a_grids
    /// @param a_domain
    /// @param a_boxType
    CFRegion(const DisjointBoxLayout& a_grids,
             const ProblemDomain&     a_domain,
             const IntVect&           a_boxType = IntVect::Zero)
    {
        this->define(a_grids, a_domain, a_boxType);
    }

    /// @brief Defines this object.  Existing information is overriden.
    /// @param a_grids
    /// @param a_domain
    /// @param a_boxType
    void
    define(const DisjointBoxLayout& a_grids,
           const ProblemDomain&     a_domain,
           const IntVect&           a_boxType = IntVect::Zero);

    /// @brief Is there a CFI?
    bool
    isEmpty() const noexcept
    {
        return m_isEmpty;
    }

    const CFIVS&
    getCFIVS(const DataIndex&      a_dit,
             const int             a_bdryDir,
             const Side::LoHiSide& a_side) const
    {
        CH_assert(m_defined);
        if (a_side == Side::Lo) {
            return m_loCFIVS[a_bdryDir][a_dit];
        } else {
            return m_hiCFIVS[a_bdryDir][a_dit];
        }
    }

    const CFIVS&
    loCFIVS(const DataIndex& a_dit, const int a_bdryDir) const
    {
        CH_assert(m_defined);
        return (m_loCFIVS[a_bdryDir])[a_dit];
    }

    const CFIVS&
    hiCFIVS(const DataIndex& a_dit, const int a_bdryDir) const
    {
        CH_assert(m_defined);
        return (m_hiCFIVS[a_bdryDir])[a_dit];
    }

protected:
    static std::vector<Box>
    buildPeriodicVector(const ProblemDomain&     a_fineDomain,
                        const DisjointBoxLayout& a_fineBoxes,
                        const IntVect&           a_boxType);

    bool              m_defined = false;
    IntVect           m_boxType = IntVect::Zero;
    bool              m_isEmpty = true;

    std::array<LayoutData<CFIVS>, SpaceDim> m_loCFIVS;
    std::array<LayoutData<CFIVS>, SpaceDim> m_hiCFIVS;
};


#include "NamespaceFooter.H"
#endif
