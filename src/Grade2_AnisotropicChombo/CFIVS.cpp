#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "LayoutIterator.H"
#include "DataIterator.H"

#include "CFIVS.H"
#include "NeighborIterator.H"
#include "CH_Timer.H"
#include "LoHiSide.H"
#include "ProblemDomain.H"
#include "DisjointBoxLayout.H"
#include "AnisotropicRefinementTools.H"
#include "NamespaceHeader.H"

long long CFIVS::s_packCount   = 0;
long long CFIVS::s_sparseCount = 0;


// -----------------------------------------------------------------------------
void
CFIVS::define(const IntVectSet & a_IVS)
{
    this->decrementCounts();
    if (a_IVS.isEmpty())  // Avoid compacting
    {
        m_IVS.define();
        m_empty  = true;
        m_packed = false;
    } else {
        m_IVS.define(a_IVS);
        this->packIVS();
    }
    m_defined = true;
}


// -----------------------------------------------------------------------------
void
CFIVS::define(const DataIndex&         a_dataIndex,
              const DisjointBoxLayout& a_grids,
              const ProblemDomain&     a_domain,
              const Box&               a_ghostBox)
{
    CH_assert(a_grids.isSorted());  //**Not sufficient.  Need neighbours defined!
    this->decrementCounts();

    const Box domBox = [&]() {
        Box bx = a_domain.domainBox();
        for (int dir = 0; dir < SpaceDim; ++dir) {
            if (!a_domain.isPeriodic(dir)) continue;

            const int numGhosts = a_ghostBox.size(dir) - a_ghostBox.type()[dir];
            bx.grow(dir, numGhosts);
        }
        bx.convert(a_ghostBox.type());
        return bx;
    }();

    m_IVS.define(a_ghostBox & domBox);
    if (!m_IVS.isEmpty()) {
        m_IVS -= a_grids[a_dataIndex].convert(a_ghostBox.type());
        NeighborIterator nit(a_grids);
        for (nit.begin(a_dataIndex); nit.ok(); ++nit) {
            m_IVS -= nit.box().convert(a_ghostBox.type());
        }
    }

    this->packIVS();
    m_defined = true;
}


// -----------------------------------------------------------------------------
void
CFIVS::decrementCounts()
{
    if (m_defined && !m_empty) {
        if (m_packed) {
            --s_packCount;
        } else {
            --s_sparseCount;
        }
    }
    CH_assert(s_packCount >= 0);
    CH_assert(s_sparseCount >= 0);
}


// -----------------------------------------------------------------------------
void
CFIVS::packIVS()
{
    m_IVS.compact();
    if (m_IVS.isEmpty()) {
        m_empty  = true;
        m_packed = false;
    } else {
        m_empty     = false;
        m_packedBox = m_IVS.minBox();
        if (m_packedBox.numPts() == m_IVS.numPts()) {
            ++s_packCount;
            m_packed = true;
        } else {
            ++s_sparseCount;
            m_packed = false;
        }
    }
}


#include "NamespaceFooter.H"
