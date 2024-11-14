#include "BdryIter.H"
#include "DataIterator.H"
#include "Comm.H"


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
PhysBdryIter::PhysBdryIter()
: Lockable()
, m_isDefined(false)
{}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
PhysBdryIter::PhysBdryIter(const DisjointBoxLayout& a_grids,
                           const SideArray&         a_activeSides)
: m_isDefined(false)
{
    this->define(a_grids, a_activeSides);
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void
PhysBdryIter::define(const DisjointBoxLayout& a_grids,
                     const SideArray&         a_activeSides)

{
    m_isDefined = true;
    m_grids     = a_grids;
    m_idx       = 0;

    // Gather data structures.
    const ProblemDomain& domain = m_grids.physDomain();
    const Box&           domBox = domain.domainBox();
    DataIterator         dit    = m_grids.dataIterator();

    // Find the domain boundaries.
    Box domBdry[CH_SPACEDIM][2];
    for (int dir = 0; dir < SpaceDim; ++dir) {
        domBdry[dir][0] = bdryBox(domBox, dir, Side::Lo);
        domBdry[dir][1] = bdryBox(domBox, dir, Side::Hi);
    }
    Box overlap;

    // Capacity estimate
    {
        Real numGrids         = Real(m_grids.size());
        Real power            = 1.0 / Real(SpaceDim);
        int  capacityEstimate = int(pow(numGrids + 1.0, power));
        m_elems.reserve(capacityEstimate);
    }

    // Find boundaries
    for (dit.reset(); dit.ok(); ++dit) {
        const Box& ccValid = m_grids[dit];

        for (int dir = 0; dir < SpaceDim; ++dir) {
            if (domain.isPeriodic(dir)) continue;

            for (SideIterator sit; sit.ok(); ++sit) {
                if (a_activeSides[dir][sit()] == 0) continue;

                // Are we at a boundary?
                overlap = bdryBox(ccValid, dir, sit());
                overlap &= domBdry[dir][sit()];
                if (overlap.isEmpty()) continue;

                PhysBdryElement newElem;
                newElem.di            = dit();
                newElem.ccValidBox    = ccValid;
                newElem.fcBdryBox     = bdryBox(ccValid, dir, sit());
                newElem.ccAdjGhostBox = adjCellBox(ccValid, dir, sit(), 1);
                newElem.dir           = dir;
                newElem.side          = sit();
                newElem.iside         = int(newElem.side);
                newElem.isign         = sign(newElem.side);
                newElem.rsign         = Real(newElem.isign);

                m_elems.push_back(newElem);
            } // end loop over boundary sides (sit)
        } // end loop over boundary directions (dir)
    } // end loop over grids (dit)

    m_numElems = m_elems.size();

    int globalSize = m_numElems;
    Comm::reduce(globalSize, MPI_SUM); // Bottleneck!

    m_isEmpty = (globalSize == 0);
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void
PhysBdryIter::clear()
{
    m_isDefined = false;
    m_elems.clear();
    m_numElems = 0;
    m_idx      = (unsigned int)(-1);
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
CFIIter::CFIIter()
: Lockable()
, m_isDefined(false)
{}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
CFIIter::CFIIter(const DisjointBoxLayout& a_grids,
                 const CFRegion&          a_cfRegion,
                 const SideArray&         a_activeSides)
: m_isDefined(false)
{
    this->define(a_grids, a_cfRegion, a_activeSides);
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void
CFIIter::define(const DisjointBoxLayout& a_grids,
                const CFRegion&          a_cfRegion,
                const SideArray&         a_activeSides)
{
    m_isDefined = true;
    m_grids     = a_grids;
    m_cfRegion  = a_cfRegion;
    m_idx       = 0;

    // Capacity estimate. Unless you have a better idea...
    m_elems.reserve(m_grids.size());

    // Find boundaries
    DataIterator dit = m_grids.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        const Box& ccValid = m_grids[dit];

        for (int dir = 0; dir < SpaceDim; ++dir) {
            for (SideIterator sit; sit.ok(); ++sit) {
                if (a_activeSides[dir][sit()] == 0) continue;

                const CFIVS& cfivs =
                    ((sit() == Side::Lo) ? m_cfRegion.loCFIVS(dit(), dir)
                                         : m_cfRegion.hiCFIVS(dit(), dir));

                if (cfivs.isEmpty()) continue;

                CFIElement newElem(cfivs);
                newElem.di            = dit();
                newElem.ccValidBox    = ccValid;
                newElem.fcBdryBox     = bdryBox(ccValid, dir, sit());
                newElem.ccAdjGhostBox = adjCellBox(ccValid, dir, sit(), 1);
                newElem.dir           = dir;
                newElem.side          = sit();
                newElem.iside         = int(newElem.side);
                newElem.isign         = sign(newElem.side);
                newElem.rsign         = Real(newElem.isign);

                m_elems.push_back(newElem);
            } // end loop over boundary sides (sit)
        } // end loop over boundary directions (dir)
    } // end loop over grids (dit)

    m_numElems = m_elems.size();

    int globalSize = m_numElems;
    Comm::reduce(globalSize, MPI_SUM);

    m_isEmpty = (globalSize == 0);
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void
CFIIter::clear()
{
    CH_assert(this->isUnlocked());
    m_isDefined = false;
    m_cfRegion  = CFRegion();
    m_elems.clear();
    m_numElems = 0;
    m_idx      = (unsigned int)(-1);
}
