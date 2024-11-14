#include "LevelData.H"
#ifdef USE_LEVELDATAFLUXBOXMOD

#if defined(CH_USE_THRUST) || defined(__NVCC__)
#   include "CudaFluxBox.hpp" // keep this here to prevent issues with macros
#endif

// #include <cstdlib>
// #include <algorithm>
// using std::sort;

// #include "parstream.H"
// #include "CH_Timer.H"
// #include <float.h>

#include "NeighborIterator.H"
#include "Comm.H"

#include "NamespaceHeader.H"

// =============================================================================
// Construction / destruction

// -----------------------------------------------------------------------------
void
LevelData<FluxBox>::clear()
{
    D_TERM(
    m_alias[0].clear();,
    m_alias[1].clear();,
    m_alias[2].clear();)

    D_TERM(
    m_exCornerCopier2[0].clear();,
    m_exCornerCopier2[1].clear();,
    m_exCornerCopier2[2].clear();)

    D_TERM(
    m_exCornerCopier1[0].clear();,
    m_exCornerCopier1[1].clear();,
    m_exCornerCopier1[2].clear();)

    D_TERM(
    m_exCopier[0].clear();,
    m_exCopier[1].clear();,
    m_exCopier[2].clear();)

    m_ghost             = IntVect(D_DECL(-1, -1, -1));
    m_disjointBoxLayout = DisjointBoxLayout();
    BoxLayoutData<FluxBox>::clear();
}


// -----------------------------------------------------------------------------
LevelData<FluxBox>::LevelData()
{
    this->clear();
}


// -----------------------------------------------------------------------------
LevelData<FluxBox>::~LevelData()
{
    this->clear();
}


// -----------------------------------------------------------------------------
LevelData<FluxBox>::LevelData(const DisjointBoxLayout&    dp,
                              int                         comps,
                              const IntVect&              ghost,
                              const DataFactory<FluxBox>& a_factory)
: m_disjointBoxLayout(dp), m_ghost(ghost)
{
#ifdef CH_MPI
    this->numSends    = 0;
    this->numReceives = 0;
#endif
    this->m_boxLayout = dp;
    this->m_comps     = comps;
    this->m_isdefined = true;

    if (!dp.isClosed()) {
        MayDay::Error(
            "non-disjoint DisjointBoxLayout: "
            "LevelData<FluxBox>::LevelData(const DisjointBoxLayout& dp, int "
            "comps)");
    }

    Interval interval(0, comps - 1);
    this->allocateGhostVector(a_factory, ghost);
    this->setVector(*this, interval, interval);  // Does nothing.

    for (int d = 0; d < SpaceDim; ++d) {
        constexpr bool doValidCorners = true;
        m_exCopier[d].defineValidExchange(
            m_disjointBoxLayout, m_ghost, d, doValidCorners);

        m_exCornerCopier1[d].defineInvalidCornerExchange1(
            m_disjointBoxLayout, m_ghost, d);

        m_exCornerCopier2[d].defineInvalidCornerExchange2(
            m_disjointBoxLayout, m_ghost, d);

        FABAliasFlBxDataFactory factory(this, this->interval(), d);
        m_alias[d].define(
            m_disjointBoxLayout, this->nComp(), m_ghost, factory);
    }
}


// -----------------------------------------------------------------------------
void
LevelData<FluxBox>::define(const DisjointBoxLayout&    dp,
                           int                         comps,
                           const IntVect&              ghost,
                           const DataFactory<FluxBox>& a_factory)
{
    CH_TIME("LevelData<FluxBox>::define(dbl,comps,ghost,factory)");
    // clear exchange copier if it's already been defined
    if (this->m_isdefined) {
        D_TERM(
        m_alias[0].clear();,
        m_alias[1].clear();,
        m_alias[2].clear();)

        D_TERM(
        m_exCornerCopier2[0].clear();,
        m_exCornerCopier2[1].clear();,
        m_exCornerCopier2[2].clear();)

        D_TERM(
        m_exCornerCopier1[0].clear();,
        m_exCornerCopier1[1].clear();,
        m_exCornerCopier1[2].clear();)

        D_TERM(
        m_exCopier[0].clear();,
        m_exCopier[1].clear();,
        m_exCopier[2].clear();)
    }

    this->m_isdefined = true;
    if (!dp.isClosed()) {
        MayDay::Error(
            "non-disjoint DisjointBoxLayout: LevelData<FluxBox>::define(const "
            "DisjointBoxLayout& dp,....)");
    }
    if (comps <= 0) {
        MayDay::Error(
            "LevelData::LevelData(const BoxLayout& dp, int comps)  comps<=0");
    }
    this->m_comps     = comps;
    this->m_boxLayout = dp;

    m_disjointBoxLayout = dp;
    m_ghost             = ghost;

    // Interval interval(0, comps-1);
    this->allocateGhostVector(a_factory, ghost);
    //  this->setVector(*this, interval, interval);

    for (int d = 0; d < SpaceDim; ++d) {
        constexpr bool doValidCorners = true;
        m_exCopier[d].defineValidExchange(
            m_disjointBoxLayout, m_ghost, d, doValidCorners);

        m_exCornerCopier1[d].defineInvalidCornerExchange1(
            m_disjointBoxLayout, m_ghost, d);

        m_exCornerCopier2[d].defineInvalidCornerExchange2(
            m_disjointBoxLayout, m_ghost, d);

        FABAliasFlBxDataFactory factory(this, this->interval(), d);
        m_alias[d].define(
            m_disjointBoxLayout, this->nComp(), m_ghost, factory);
    }
}


// -----------------------------------------------------------------------------
void
LevelData<FluxBox>::define(const LevelData<FluxBox>&   da,
                           const DataFactory<FluxBox>& a_factory)
{
    CH_TIME("LevelData<FluxBox>::define(LevelData<FluxBox>,factory)");
    // clear exchange copier if it's already been defined
    if (this->m_isdefined) {
        D_TERM(
        m_alias[0].clear();,
        m_alias[1].clear();,
        m_alias[2].clear();)

        D_TERM(
        m_exCornerCopier2[0].clear();,
        m_exCornerCopier2[1].clear();,
        m_exCornerCopier2[2].clear();)

        D_TERM(
        m_exCornerCopier1[0].clear();,
        m_exCornerCopier1[1].clear();,
        m_exCornerCopier1[2].clear();)

        D_TERM(
        m_exCopier[0].clear();,
        m_exCopier[1].clear();,
        m_exCopier[2].clear();)
    }
    this->m_isdefined = true;
    if (this == &da) return;
    m_disjointBoxLayout = da.m_disjointBoxLayout;
    this->m_boxLayout   = da.m_disjointBoxLayout;
    this->m_comps       = da.m_comps;
    m_ghost             = da.m_ghost;

    Interval srcAnddest(0, this->m_comps - 1);

    this->allocateGhostVector(a_factory, m_ghost);
    da.copyTo(*this);

    for (int d = 0; d < SpaceDim; ++d) {
        constexpr bool doValidCorners = true;
        m_exCopier[d].defineValidExchange(
            m_disjointBoxLayout, m_ghost, d, doValidCorners);

        m_exCornerCopier1[d].defineInvalidCornerExchange1(
            m_disjointBoxLayout, m_ghost, d);

        m_exCornerCopier2[d].defineInvalidCornerExchange2(
            m_disjointBoxLayout, m_ghost, d);

        FABAliasFlBxDataFactory factory(this, this->interval(), d);
        m_alias[d].define(
            m_disjointBoxLayout, this->nComp(), m_ghost, factory);
    }
}


// -----------------------------------------------------------------------------
void
LevelData<FluxBox>::define(const LevelData<FluxBox>&   da,
                           const Interval&             comps,
                           const DataFactory<FluxBox>& a_factory)
{
    CH_TIME("LevelData<FluxBox>::define(LevelData<FluxBox>,comps,factory)");
    // clear exchange copier if it's already been defined
    if (this->m_isdefined) {
        D_TERM(
        m_alias[0].clear();,
        m_alias[1].clear();,
        m_alias[2].clear();)

        D_TERM(
        m_exCornerCopier2[0].clear();,
        m_exCornerCopier2[1].clear();,
        m_exCornerCopier2[2].clear();)

        D_TERM(
        m_exCornerCopier1[0].clear();,
        m_exCornerCopier1[1].clear();,
        m_exCornerCopier1[2].clear();)

        D_TERM(
        m_exCopier[0].clear();,
        m_exCopier[1].clear();,
        m_exCopier[2].clear();)
    }
    this->m_isdefined = true;
    if (this == &da) {
        MayDay::Error(
            " LevelData<FluxBox>::define(const LevelData<FluxBox>& da, const "
            "Interval& comps) called with 'this'");
    }
    CH_assert(comps.size() > 0);
    // this line doesn't make any sense!
    // CH_assert(comps.end()<=this->m_comps);
    CH_assert(comps.begin() >= 0);

    m_disjointBoxLayout = da.m_disjointBoxLayout;
    this->m_boxLayout   = da.m_disjointBoxLayout;

    this->m_comps = comps.size();

    m_ghost = da.m_ghost;

    Interval dest(0, this->m_comps - 1);
    this->allocateGhostVector(a_factory, m_ghost);
    this->setVector(da, comps, dest);


    for (int d = 0; d < SpaceDim; ++d) {
        constexpr bool doValidCorners = true;
        m_exCopier[d].defineValidExchange(
            m_disjointBoxLayout, m_ghost, d, doValidCorners);

        m_exCornerCopier1[d].defineInvalidCornerExchange1(
            m_disjointBoxLayout, m_ghost, d);

        m_exCornerCopier2[d].defineInvalidCornerExchange2(
            m_disjointBoxLayout, m_ghost, d);

        FABAliasFlBxDataFactory factory(this, this->interval(), d);
        m_alias[d].define(
            m_disjointBoxLayout, this->nComp(), m_ghost, factory);
    }
}


// =============================================================================
// LevelData --> LevelData copyTo functions

// -----------------------------------------------------------------------------
void
LevelData<FluxBox>::copyTo(LevelData<FluxBox>& a_dest) const
{
    CH_assert(this->nComp() == a_dest.nComp());
    this->copyTo(this->interval(), a_dest, a_dest.interval());
}


// -----------------------------------------------------------------------------
void
LevelData<FluxBox>::copyTo(const Interval&     a_srcComps,
                           LevelData<FluxBox>& a_dest,
                           const Interval&     a_destComps) const
{
    if (this == &a_dest) {
        MayDay::Error(
            "src == dest in copyTo function. Perhaps you want exchange?");
    }

    if (this->boxLayout() == a_dest.boxLayout() &&
        a_dest.ghostVect() == IntVect::Zero) {
        // parallel direct copy here, no communication issues
        for (DataIterator dit(this->dataIterator()); dit.ok(); ++dit) {
            a_dest[dit].copy(this->box(dit),
                             a_destComps,
                             this->box(dit),
                             this->operator[](dit),
                             a_srcComps);
        }
    } else {
        CopierT copier;
        for (int d = 0; d < SpaceDim; ++d) {
            const BoxLayout&     srcLayout  = this->boxLayout();
            const BoxLayout&     destLayout = a_dest.boxLayout();
            const ProblemDomain& domain     = m_disjointBoxLayout.physDomain();
            const IntVect&       destGhost  = a_dest.ghostVect();

            copier[d].define(srcLayout, destLayout, domain, destGhost, d);
        }

        this->copyTo(a_srcComps, a_dest, a_destComps, copier);
    }
}


// -----------------------------------------------------------------------------
void
LevelData<FluxBox>::copyTo(LevelData<FluxBox>& a_dest,
                           const CopierT&      a_copier,
                           const LDOpT&        a_op) const
{
    CH_assert(this->nComp() == a_dest.nComp());
    this->copyTo(this->interval(), a_dest, a_dest.interval(), a_copier, a_op);
}


// -----------------------------------------------------------------------------
void
LevelData<FluxBox>::copyTo(const Interval&     a_srcComps,
                           LevelData<FluxBox>& a_dest,
                           const Interval&     a_destComps,
                           const CopierT&      a_copier,
                           const LDOpT&        a_op) const
{
    if ((BoxLayoutData<FluxBox>*)this == &a_dest) return;

    if (this->boxLayout() == a_dest.boxLayout() &&
        a_dest.ghostVect() == IntVect::Zero) {
        // parallel direct copy here, no communication issues
        for (DataIterator dit(this->dataIterator()); dit.ok(); ++dit) {
            a_dest[dit].copy(this->box(dit),
                             a_destComps,
                             this->box(dit),
                             this->operator[](dit),
                             a_srcComps);
        }
    } else {
        std::array<BoxLayoutData<FArrayBox>, SpaceDim> destAlias;

        for (int d = 0; d < SpaceDim; ++d) {
            CH_assert(m_alias[d].isDefined());
            CH_assert(m_alias[d].getBoxes() == this->getBoxes());
            CH_assert(m_alias[d].ghostVect() == this->ghostVect());

            FABAliasFlBxDataFactory factory(&a_dest, a_dest.interval(), d);
            destAlias[d].define(a_dest.boxLayout(), a_dest.nComp(), factory);

            // m_alias[d].copyTo(
            //     a_srcComps, destAlias, a_destComps, a_copier[d], a_op[d]);
            if (a_op[d]) {
                m_alias[d].copyToBegin(
                    a_srcComps, destAlias[d], a_destComps, a_copier[d], *a_op[d]);
            } else {
                m_alias[d].copyToBegin(
                    a_srcComps, destAlias[d], a_destComps, a_copier[d]);
            }
        }

        for (int d = 0; d < SpaceDim; ++d) {
            if (a_op[d]) {
                m_alias[d].copyToEnd(
                    destAlias[d], a_destComps, a_copier[d], *a_op[d]);
            } else {
                m_alias[d].copyToEnd(
                    destAlias[d], a_destComps, a_copier[d]);
            }
        }
    }
}


// =============================================================================
// LevelData --> BoxLayoutData copyTo functions

// -----------------------------------------------------------------------------
void
LevelData<FluxBox>::copyTo(BoxLayoutData<FluxBox>& a_dest) const
{
    CH_assert(this->nComp() == a_dest.nComp());
    this->copyTo(this->interval(), a_dest, a_dest.interval());
}


// -----------------------------------------------------------------------------
void
LevelData<FluxBox>::copyTo(const Interval&         a_srcComps,
                           BoxLayoutData<FluxBox>& a_dest,
                           const Interval&         a_destComps) const
{
    if ((BoxLayoutData<FluxBox>*)this == &a_dest) return;

    if (this->boxLayout() == a_dest.boxLayout()) {
        // parallel direct copy here, no communication issues
        for (DataIterator dit(this->dataIterator()); dit.ok(); ++dit) {
            a_dest[dit].copy(this->box(dit),
                             a_destComps,
                             this->box(dit),
                             this->operator[](dit),
                             a_srcComps);
        }
    } else {
        CopierT copier;
        for (int d = 0; d < SpaceDim; ++d) {
            const BoxLayout&     srcLayout  = this->boxLayout();
            const BoxLayout&     destLayout = a_dest.boxLayout();
            const ProblemDomain& domain     = m_disjointBoxLayout.physDomain();
            const IntVect&       destGhost  = IntVect::Zero;

            copier[d].define(srcLayout, destLayout, domain, destGhost, d);
        }

        this->copyTo(a_srcComps, a_dest, a_destComps, copier);
    }
}


// -----------------------------------------------------------------------------
void
LevelData<FluxBox>::copyTo(BoxLayoutData<FluxBox>& a_dest,
                           const CopierT&          a_copier,
                           const LDOpT&            a_op) const
{
    CH_assert(this->nComp() == a_dest.nComp());
    this->copyTo(this->interval(), a_dest, a_dest.interval(), a_copier, a_op);
}


// -----------------------------------------------------------------------------
void
LevelData<FluxBox>::copyTo(const Interval&         a_srcComps,
                           BoxLayoutData<FluxBox>& a_dest,
                           const Interval&         a_destComps,
                           const CopierT&          a_copier,
                           const LDOpT&            a_op) const
{
    if ((BoxLayoutData<FluxBox>*)this == &a_dest) return;

    if (this->boxLayout() == a_dest.boxLayout()) {
        // parallel direct copy here, no communication issues
        for (DataIterator dit(this->dataIterator()); dit.ok(); ++dit) {
            a_dest[dit].copy(this->box(dit),
                             a_destComps,
                             this->box(dit),
                             this->operator[](dit),
                             a_srcComps);
        }
    } else {
        std::array<BoxLayoutData<FArrayBox>, SpaceDim> destAlias;

        for (int d = 0; d < SpaceDim; ++d) {
            CH_assert(m_alias[d].isDefined());
            CH_assert(m_alias[d].getBoxes() == this->getBoxes());
            CH_assert(m_alias[d].ghostVect() == this->ghostVect());

            FABAliasFlBxDataFactory factory(&a_dest, a_dest.interval(), d);
            destAlias[d].define(a_dest.boxLayout(), a_dest.nComp(), factory);

            if (a_op[d]) {
                m_alias[d].copyToBegin(
                    a_srcComps, destAlias[d], a_destComps, a_copier[d], *a_op[d]);
            } else {
                m_alias[d].copyToBegin(
                    a_srcComps, destAlias[d], a_destComps, a_copier[d]);
            }
        }

        for (int d = 0; d < SpaceDim; ++d) {
            if (a_op[d]) {
                m_alias[d].copyToEnd(
                    destAlias[d], a_destComps, a_copier[d], *a_op[d]);
            } else {
                m_alias[d].copyToEnd(
                    destAlias[d], a_destComps, a_copier[d]);
            }
        }
    }
}


// =============================================================================
// Exchange functions

// -----------------------------------------------------------------------------
void
LevelData<FluxBox>::exchange(void)
{
    this->exchange(this->interval());
}


// -----------------------------------------------------------------------------
void
LevelData<FluxBox>::exchange(const Interval& comps)
{
    if (comps == this->interval()) {
        for (int d = 0; d < SpaceDim; ++d) {
            CH_assert(m_alias[d].isDefined());
            CH_assert(m_alias[d].getBoxes() == this->getBoxes());
            CH_assert(m_alias[d].ghostVect() == this->ghostVect());

            m_alias[d].exchangeBegin(m_exCopier[d]);
        }
        for (int d = 0; d < SpaceDim; ++d) {
            m_alias[d].exchangeEnd();
        }

        for (int d = 0; d < SpaceDim; ++d) {
            m_alias[d].exchangeBegin(m_exCornerCopier1[d]);
        }
        for (int d = 0; d < SpaceDim; ++d) {
            m_alias[d].exchangeEnd();
        }

        for (int d = 0; d < SpaceDim; ++d) {
            m_alias[d].exchangeBegin(m_exCornerCopier2[d]);
        }
        for (int d = 0; d < SpaceDim; ++d) {
            m_alias[d].exchangeEnd();
        }

    } else {
        for (int d = 0; d < SpaceDim; ++d) {
            CH_assert(m_alias[d].isDefined());
            CH_assert(m_alias[d].getBoxes() == this->getBoxes());
            CH_assert(m_alias[d].ghostVect() == this->ghostVect());

            m_alias[d].exchange(comps, m_exCopier[d]);
            m_alias[d].exchange(comps, m_exCornerCopier1[d]);
            m_alias[d].exchange(comps, m_exCornerCopier2[d]);
        }
    }
}


// -----------------------------------------------------------------------------
void
LevelData<FluxBox>::exchange(const Interval& /*a_comps*/,
                             const CopierT& a_copier)
{
    for (int d = 0; d < SpaceDim; ++d) {
        CH_assert(m_alias[d].isDefined());
        CH_assert(m_alias[d].getBoxes() == this->getBoxes());
        CH_assert(m_alias[d].ghostVect() == this->ghostVect());
        CH_assert(a_copier[d].getGhostVect() <= this->ghostVect());

        m_alias[d].exchange(this->interval(), a_copier[d]);
    }
}


// -----------------------------------------------------------------------------
void
LevelData<FluxBox>::exchange(const CopierT& a_copier)
{
    this->exchange(this->interval(), a_copier);
}


// -----------------------------------------------------------------------------
void
LevelData<FluxBox>::exchangeNoOverlap(const CopierT& a_copier)
{
    for (int d = 0; d < SpaceDim; ++d) {
        CH_assert(m_alias[d].isDefined());
        CH_assert(m_alias[d].getBoxes() == this->getBoxes());
        CH_assert(m_alias[d].ghostVect() == this->ghostVect());
        CH_assert(a_copier[d].getGhostVect() <= this->ghostVect());

        m_alias[d].exchangeNoOverlap(a_copier[d]);
    }
}


// -----------------------------------------------------------------------------
void
LevelData<FluxBox>::exchangeBegin(const CopierT& a_copier)
{
    for (int d = 0; d < SpaceDim; ++d) {
        CH_assert(m_alias[d].isDefined());
        CH_assert(m_alias[d].getBoxes() == this->getBoxes());
        CH_assert(m_alias[d].ghostVect() == this->ghostVect());
        CH_assert(a_copier[d].getGhostVect() <= this->ghostVect());

        m_alias[d].exchangeBegin(a_copier[d]);
    }
}


// -----------------------------------------------------------------------------
void
LevelData<FluxBox>::exchangeEnd()
{
    for (int d = 0; d < SpaceDim; ++d) {
        CH_assert(m_alias[d].isDefined());
        CH_assert(m_alias[d].getBoxes() == this->getBoxes());
        CH_assert(m_alias[d].ghostVect() == this->ghostVect());

        m_alias[d].exchangeEnd();
    }
}


// =============================================================================
// Miscellaneous

// -----------------------------------------------------------------------------
void
LevelData<FluxBox>::apply(void (*a_Function)(const Box& box,
                                             int        comps,
                                             FluxBox&   t))
{
    for (DataIterator it(this->dataIterator()); it.ok(); ++it) {
        a_Function(m_disjointBoxLayout.get(it()),
                   this->m_comps,
                   *(this->m_vector[it().datInd()]));
    }
}


// -----------------------------------------------------------------------------
void
LevelData<FluxBox>::apply(const ApplyFunctor& f)
{
    for (DataIterator it(this->dataIterator()); it.ok(); ++it) {
        //     unsigned int index = this->m_boxLayout.lindex(it());
        f(m_disjointBoxLayout.get(it()),
          this->m_comps,
          *(this->m_vector[it().datInd()]));
    }
}


// -----------------------------------------------------------------------------
void
LevelData<FluxBox>::degenerate(LevelData<FluxBox>& a_to,
                               const SliceSpec&    a_sliceSpec) const
{
    DisjointBoxLayout toDBL;
    m_disjointBoxLayout.degenerate(toDBL, a_sliceSpec);
    IntVect toGhost;
    for (int i = 0; i < CH_SPACEDIM; ++i) {
        if (i != a_sliceSpec.direction) {
            toGhost[i] = m_ghost[i];
        } else {
            toGhost[i] = 0;
        }
    }
    a_to.define(toDBL, this->nComp(), toGhost);
    copyTo(a_to);
}


// =============================================================================
// Overidden virtual functions

// -----------------------------------------------------------------------------
void
LevelData<FluxBox>::define(const BoxLayout&            /*dp*/,
                           int                         /*comps*/,
                           const DataFactory<FluxBox>& /*a_factory*/)
{
    MayDay::Error("LevelData<FluxBox>::define called with BoxLayout input");
}


// -----------------------------------------------------------------------------
void
LevelData<FluxBox>::define(const BoxLayout& /*dp*/)
{
    MayDay::Error("LevelData<FluxBox>::define called with BoxLayout input");
}


// -----------------------------------------------------------------------------
void
LevelData<FluxBox>::define(const BoxLayoutData<FluxBox>& /*da*/,
                           const DataFactory<FluxBox>&   /*a_factory*/)
{
    MayDay::Error("LevelData<FluxBox>::define called with BoxLayout input");
}


// -----------------------------------------------------------------------------
void
LevelData<FluxBox>::define(const BoxLayoutData<FluxBox>& /*da*/,
                           const Interval&               /*comps*/,
                           const DataFactory<FluxBox>&   /*a_factory*/)
{
    MayDay::Error("LevelData<FluxBox>::define called with BoxLayout input");
}
//-----------------------------------------------------------------------
#ifdef CH_USE_PYTHON
  PyObject* LevelData<FluxBox>::pack()  {

    Box Domain = this->disjointBoxLayout().physDomain().domainBox();
    PyObject *pDomain = Domain.pack();

    DataIterator dit = this->dataIterator();
    PyObject *pLDFABs = PyTuple_New(dit.size());
    PyObject *pLDBoxes = PyTuple_New(this->getBoxes().size());
    PyObject *pIDs = PyTuple_New(this->getBoxes().size());

    int i = 0;
    for (dit.reset(); dit.ok(); ++dit, ++i) {
      FluxBox &A = this->operator[](dit);
      PyObject *pFAB = A.pack();
      PyTuple_SetItem(pLDFABs, i, pFAB);
    }
    LayoutIterator lit = this->getBoxes().layoutIterator();
    i = 0;
    for (lit.reset(); lit.ok(); ++lit, ++i) {
      Box B = this->getBoxes()[lit];
      PyObject *pBox = B.pack();
      int PID=this->getBoxes().procID(lit());
      PyObject *pId= PyTuple_New(2);
      PyTuple_SetItem(pId, 0, PyLong_FromLong((int)PID));
      PyTuple_SetItem(pId, 1, PyUnicode_FromString(std::string("int").c_str()));


      PyTuple_SetItem(pLDBoxes, i, pBox);
      PyTuple_SetItem(pIDs, i, pId);
    }

    PyObject *pArg = PyTuple_New(4 + 1);
    PyTuple_SetItem(pArg, 0, pLDFABs);
    PyTuple_SetItem(pArg, 1, pLDBoxes);
    PyTuple_SetItem(pArg, 2, pIDs);
    PyTuple_SetItem(pArg, 3, pDomain);

    std::string Label = "LevelDataFluxBox";
    PyObject *pLabel = PyUnicode_FromString(Label.c_str());
    PyTuple_SetItem(pArg, 4, pLabel);
    return pArg;
  }


  PyObject* LevelData<FluxBox>::pack() const  {

    Box Domain = this->disjointBoxLayout().physDomain().domainBox();
    PyObject *pDomain = Domain.pack();

    DataIterator dit = this->dataIterator();
    PyObject *pLDFABs = PyTuple_New(dit.size());
    PyObject *pLDBoxes = PyTuple_New(this->getBoxes().size());
    PyObject *pIDs = PyTuple_New(this->getBoxes().size());

    int i = 0;
    for (dit.reset(); dit.ok(); ++dit, ++i) {
      const FluxBox &A = this->operator[](dit);
      PyObject *pFAB = A.pack();
      PyTuple_SetItem(pLDFABs, i, pFAB);
    }
    LayoutIterator lit = this->getBoxes().layoutIterator();
    i = 0;
    for (lit.reset(); lit.ok(); ++lit, ++i) {
      Box B = this->getBoxes()[lit];
      PyObject *pBox = B.pack();
      int PID=this->getBoxes().procID(lit());
      PyObject *pId= PyTuple_New(2);
      PyTuple_SetItem(pId, 0, PyLong_FromLong((int)PID));
      PyTuple_SetItem(pId, 1, PyUnicode_FromString(std::string("int").c_str()));


      PyTuple_SetItem(pLDBoxes, i, pBox);
      PyTuple_SetItem(pIDs, i, pId);
    }

    PyObject *pArg = PyTuple_New(4 + 1);
    PyTuple_SetItem(pArg, 0, pLDFABs);
    PyTuple_SetItem(pArg, 1, pLDBoxes);
    PyTuple_SetItem(pArg, 2, pIDs);
    PyTuple_SetItem(pArg, 3, pDomain);
    std::string Label = "LevelDataFluxBox";
    PyObject *pLabel = PyUnicode_FromString(Label.c_str());
    PyTuple_SetItem(pArg, 4, pLabel);
    return pArg;
  }
#endif
// =============================================================================
// Aliasing functions
// =============================================================================

// -----------------------------------------------------------------------------
// This is modeled after aliasLevelData, which takes a pointer for the
// second argument.
// -----------------------------------------------------------------------------
#if defined(CH_USE_THRUST) || defined(__NVCC__)
template <typename Memory>
void
aliasFluxBox(CudaFluxBox<Memory>& a_new, CudaFluxBox<Memory>* a_oldPtr, const Interval& a_ivl)
#else
void aliasFluxBox(FluxBox& a_new, FluxBox* a_oldPtr, const Interval& a_ivl)
#endif
{
    CH_assert(a_new.m_nvar == -1); // Check if undefined
    CH_assert(a_oldPtr != nullptr);
    CH_assert(a_ivl.begin() >= 0);
    CH_assert(a_ivl.end() < a_oldPtr->nComp());
    CH_assert(a_ivl.size() > 0);

    a_new.m_bx   = a_oldPtr->box();
    a_new.m_nvar = a_ivl.size();

    if (a_new.m_fluxes.size() == 0) {
        a_new.m_fluxes.resize(SpaceDim, nullptr);
    }

    for (int dir = 0; dir < SpaceDim; dir++) {
        if (a_new.m_fluxes[dir] != nullptr) {
            delete a_new.m_fluxes[dir];
        }
        a_new.m_fluxes[dir] = new FArrayBox(a_ivl, (*a_oldPtr)[dir]);
    }
}


// -----------------------------------------------------------------------------
// The LevelData<FluxBox> version of aliasLevelData.
// -----------------------------------------------------------------------------
template <>
void
aliasLevelData<FluxBox>(LevelData<FluxBox>& a_alias,
                        LevelData<FluxBox>* a_original,
                        const Interval&     a_interval)
{
    CH_assert(a_original != nullptr);
    CH_assert(a_original->interval().contains(a_interval.begin()));
    CH_assert(a_original->interval().contains(a_interval.end()));

    AliasDataFactory<FluxBox> factory(a_original, a_interval);
    a_alias.define(a_original->disjointBoxLayout(),
                   a_interval.size(),
                   a_original->ghostVect(),
                   factory);
}


// -----------------------------------------------------------------------------
// This is modeled after aliasLevelData, which takes a pointer for the
// second argument. *a_oldPtr must have SpaceDim comps.
// -----------------------------------------------------------------------------
#if defined(CH_USE_THRUST) || defined(__NVCC__)
template <typename Memory>
void
aliasNormalFluxes(CudaFluxBox<Memory>& a_new, CudaFluxBox<Memory>* a_oldPtr)
#else
void
aliasNormalFluxes(FluxBox& a_new, FluxBox* a_oldPtr)
#endif
{
    CH_assert(a_new.m_nvar == -1); // Check if undefined
    CH_assert(a_oldPtr != nullptr);
    CH_assert(a_oldPtr->nComp() == SpaceDim);

    a_new.m_bx   = a_oldPtr->box();
    a_new.m_nvar = 1;

    if (a_new.m_fluxes.size() == 0) {
        a_new.m_fluxes.resize(SpaceDim, nullptr);
    }

    for (int dir = 0; dir < SpaceDim; dir++) {
        if (a_new.m_fluxes[dir] != nullptr) {
            delete a_new.m_fluxes[dir];
        }
        a_new.m_fluxes[dir] =
            new FArrayBox(Interval(dir, dir), (*a_oldPtr)[dir]);
    }
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
AliasNormalDataFactory::
AliasNormalDataFactory(BoxLayoutData<FluxBox>* a_origPtr)
: m_origPtr(a_origPtr)
{
    // TODO: Sanity checks.
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
AliasNormalDataFactory::~AliasNormalDataFactory()
{
}


// -----------------------------------------------------------------------------
// factory function.  creates a new 'T' object using an aliased dataPtr for T
// creates a new 'T' object and returns a pointer to it.  Responsiblitly
// for calling operator 'delete' on this pointer is passed to the user.
// -----------------------------------------------------------------------------
FluxBox*
AliasNormalDataFactory::create(const Box&       /*a_box*/,
                               int              a_nComps [[maybe_unused]],
                               const DataIndex& a_di) const
{
    CH_assert(a_nComps == 1);
    FluxBox* newPtr = new FluxBox;
    aliasNormalFluxes(*newPtr, &((*m_origPtr)[a_di]));
    return newPtr;
}


// -----------------------------------------------------------------------------
// Similar to aliasLevelData. *a_original must have SpaceDim comps.
// -----------------------------------------------------------------------------
#if defined(CH_USE_THRUST) || defined(__NVCC__)
template <>
#endif
void
aliasNormalFluxes(LevelData<FluxBox>& a_alias,
                  LevelData<FluxBox>* a_original)
{
    CH_assert(a_original != nullptr);
    CH_assert(a_original->interval() == Interval(0, SpaceDim - 1));

    AliasNormalDataFactory factory(a_original);
    a_alias.define(a_original->disjointBoxLayout(),
                   1,
                   a_original->ghostVect(),
                   factory);
}


#include "NamespaceFooter.H"
#endif  // USE_LEVELDATAFLUXBOXMOD
