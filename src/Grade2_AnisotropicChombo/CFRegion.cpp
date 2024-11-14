#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CFRegion.H"
#include "Debug.H"
#include "Comm.H"

#include "NamespaceHeader.H"



// -----------------------------------------------------------------------------
void
CFRegion::define(const DisjointBoxLayout& a_grids,
                 const ProblemDomain&     a_domain,
                 const IntVect&           a_boxType)
{
    int localSize = 0;
    for (int bdryDir = 0; bdryDir < SpaceDim; ++bdryDir) {
        m_loCFIVS[bdryDir].define(a_grids);
        m_hiCFIVS[bdryDir].define(a_grids);

        for (DataIterator dit(a_grids); dit.ok(); ++dit) {
            Box ghostBox = adjCellBox(a_grids[dit], bdryDir, Side::Lo, 1)
                          .convert(a_boxType);
            m_loCFIVS[bdryDir][dit].define(dit(), a_grids, a_domain, ghostBox);
            localSize += m_loCFIVS[bdryDir][dit].getIVS().numPts();


            ghostBox = adjCellBox(a_grids[dit], bdryDir, Side::Hi, 1)
                      .convert(a_boxType);
            m_hiCFIVS[bdryDir][dit].define(dit(), a_grids, a_domain, ghostBox);
            localSize += m_hiCFIVS[bdryDir][dit].getIVS().numPts();
        }
    }

    int globalSize = localSize;
    Comm::reduce(globalSize, MPI_SUM);

    m_isEmpty = (globalSize == 0);
    m_boxType = a_boxType;
    m_defined = true;
}

// -----------------------------------------------------------------------------
CFRegion&
CFRegion::operator=(const CFRegion& a_rhs)
{
    if (a_rhs.m_defined) {
        const BoxLayout& grids = (a_rhs.m_loCFIVS[0]).boxLayout();
        for (int bdryDir = 0; bdryDir < SpaceDim; ++bdryDir) {
            m_loCFIVS[bdryDir].define(grids);
            m_hiCFIVS[bdryDir].define(grids);

            for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit) {
                m_loCFIVS[bdryDir][dit] = a_rhs.m_loCFIVS[bdryDir][dit];
                m_hiCFIVS[bdryDir][dit] = a_rhs.m_hiCFIVS[bdryDir][dit];
            }
        }
        m_boxType = a_rhs.m_boxType;
        m_isEmpty = a_rhs.m_isEmpty;
        m_defined = true;
    } else {
        m_boxType = IntVect::Zero;
        m_isEmpty = true;
        m_defined = false;
    }
    return *this;
}


// -----------------------------------------------------------------------------
Vector<Box>
CFRegion::buildPeriodicVector(const ProblemDomain&     a_fineDomain,
                              const DisjointBoxLayout& a_fineBoxes,
                              const IntVect&           a_boxType)
{
    // Local utility:
    // bx is at a periodic boundary if test.atPeriodicBoundary(bx) is true.
    class PeriodicBoundaryTest
    {
    public:
        PeriodicBoundaryTest(const ProblemDomain& aa_domain,
                             const IntVect&       aa_boxType)
        : m_isPeriodic(aa_domain.isPeriodic())
        , m_periodicTestBox(aa_domain.domainBox())
        {
            m_periodicTestBox.convert(aa_boxType);

            if (!m_isPeriodic) return;
            for (int idir = 0; idir < SpaceDim; ++idir) {
                if (!aa_domain.isPeriodic(idir)) continue;
                m_periodicTestBox.grow(idir, -1);
            }
        }

        inline bool
        atPeriodicBoundary(const Box& a_box) const
        {
            return m_isPeriodic && !m_periodicTestBox.contains(a_box);
        }

    private:
        bool m_isPeriodic;
        Box  m_periodicTestBox;
    };
    const PeriodicBoundaryTest test(a_fineDomain, a_boxType);


    // We are ready. Build periodicTestBox.
    std::vector<Box> periodicVector;
    periodicVector.reserve(a_fineBoxes.size());

    for (LayoutIterator lit = a_fineBoxes.layoutIterator(); lit.ok(); ++lit) {
        // Start with the valid FC box.
        const Box& box = a_fineBoxes[lit()].convert(a_boxType);
        periodicVector.push_back(box);

        // Then add periodic images, if needed.
        if (test.atPeriodicBoundary(box)) {
            ShiftIterator shiftIt = a_fineDomain.shiftIterator();
            IntVect shiftMult(a_fineDomain.domainBox().size());
            Box shiftedBox(box);

            for (shiftIt.begin(); shiftIt.ok(); ++shiftIt) {
                IntVect shiftVect = shiftMult * shiftIt();
                shiftedBox.shift(shiftVect);
                periodicVector.push_back(shiftedBox);
                shiftedBox.shift(-shiftVect);
            }
        }
    }

    std::sort(periodicVector.begin(), periodicVector.end());
    return periodicVector;
}


#include "NamespaceFooter.H"
