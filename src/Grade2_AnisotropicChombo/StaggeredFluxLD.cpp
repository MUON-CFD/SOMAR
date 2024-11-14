#include "StaggeredFluxLD.H"
#include "MiscUtils.H"
#include "Debug.H"


// -----------------------------------------------------------------------------
// Default constructor. Leaves object undefined.
// -----------------------------------------------------------------------------
StaggeredFluxLD::StaggeredFluxLD()
: m_isDefined(false)
{
}


// -----------------------------------------------------------------------------
// Full constructor. Leaves object in a usable state.
// -----------------------------------------------------------------------------
StaggeredFluxLD::StaggeredFluxLD(const DisjointBoxLayout& a_grids,
                                 const bool               a_noDiags)
: m_isDefined(false)
{
    this->define(a_grids, a_noDiags);
}


// -----------------------------------------------------------------------------
// Full virtual constructor. Leaves object in a usable state.
// -----------------------------------------------------------------------------
void
StaggeredFluxLD::define(const DisjointBoxLayout& a_grids,
                        const bool               a_noDiags)
{
    for (int i = 0; i < SpaceDim; ++i) {
        // Diagonal element is CC.
        if (!a_noDiags) {
            (*this)[i][i].define(a_grids, 1, BASISV(i));
            m_exCopier[i].exchangeDefine(a_grids, BASISV(i));
        }

        // Off diagonal element is nodal in both i and j dirs.
        for (int j = i + 1; j < SpaceDim; ++j) {
            NCDataFactory<FArrayBox> ncFact(BASISV(i) + BASISV(j));
            (*this)[i][j].define(a_grids, 1, IntVect::Zero, ncFact);
            (*this)[j][i].define(a_grids, 1, IntVect::Zero, ncFact);
        }
    }

    m_grids     = a_grids;
    m_noDiags   = a_noDiags;
    m_isDefined = true;
}


// -----------------------------------------------------------------------------
// Free memory. Leaves object undefined.
// -----------------------------------------------------------------------------
void
StaggeredFluxLD::clear()
{
    m_isDefined = false;
    m_noDiags   = false;
    m_grids     = DisjointBoxLayout();

    for (int i = 0; i < SpaceDim; ++i) {
        m_exCopier[i].clear();

        for (int j = 0; j < SpaceDim; ++j) {
            (*this)[i][j].clear();
        }
    }
}


// -----------------------------------------------------------------------------
// MPI data exchange between neighboring boxes. This is blocking.
// -----------------------------------------------------------------------------
void
StaggeredFluxLD::exchange()
{
    if (m_noDiags) return;

    for (int i = 0; i < SpaceDim; ++i) {
        (*this)[i][i].exchangeBegin(m_exCopier[i]);
    }
    for (int i = 0; i < SpaceDim; ++i) {
        (*this)[i][i].exchangeEnd();
    }
}


// -----------------------------------------------------------------------------
// Sets all data to a single value.
// -----------------------------------------------------------------------------
void
StaggeredFluxLD::setVal(const Real a_val)
{
    DataIterator dit = m_grids.dataIterator();

    for (dit.reset(); dit.ok(); ++dit) {
        for (int i = 0; i < SpaceDim; ++i) {
            if (!m_noDiags) {
                (*this)[i][i][dit].setVal(a_val);
            }
            for (int j = i + 1; j < SpaceDim; ++j) {
                (*this)[i][j][dit].setVal(a_val);
                (*this)[j][i][dit].setVal(a_val);
            }
        } // i
    } // dit
}
