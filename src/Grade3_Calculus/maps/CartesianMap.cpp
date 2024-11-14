#include "CartesianMap.H"
#include "CartesianMapF_F.H"


// -----------------------------------------------------------------------------
CartesianMap::~CartesianMap()
{
    ;
}


// -----------------------------------------------------------------------------
const char*
CartesianMap::getCoorMapName() const
{
    return "Cartesian";
}


// -----------------------------------------------------------------------------
bool
CartesianMap::isUniform() const
{
    return true;
}


// -----------------------------------------------------------------------------
void
CartesianMap::interp(Vector<Real>&       a_x,
                     const Vector<Real>& a_xi,
                     const int           /*a_mu*/) const
{
    a_x = a_xi;
}


// -----------------------------------------------------------------------------
void
CartesianMap::fill_physCoor(Vector<Real>& a_x,
                            const int     a_mu,
                            const Real    a_dXi,
                            const Box&    a_box) const
{
    const size_t Nx = a_box.size(a_mu);

    // Resize a_x, if needed.
    if (a_x.size() != Nx) {
        a_x.resize(Nx);
    }

    // Compute xi at each box location.
    {
        // Compute initial xi.
        Real xi = a_box.smallEnd(a_mu) * a_dXi;
        if (a_box.type(a_mu) == 0) {
            xi += 0.5 * a_dXi;
        }

        // Step along, computing xi at each location.
        size_t idx = 0;
        while (idx < Nx) {
            a_x[idx] = xi;
            ++idx;
            xi += a_dXi;
        }
    }
}


// -----------------------------------------------------------------------------
void
CartesianMap::fill_physCoor(FArrayBox&      a_dest,
                            const int       a_destComp,
                            const int       a_mu,
                            const RealVect& a_dXi) const
{
    // Sanity checks
    CH_assert(0 <= a_destComp);
    CH_assert(a_destComp < a_dest.nComp());
    CH_assert(0 <= a_mu);
    CH_assert(a_mu < SpaceDim);

    const Real& dXimu       = a_dXi[a_mu];
    const Box&  destBox     = a_dest.box();
    const int   destBoxType = destBox.type()[a_mu];

    FORT_CARTESIAN_FILL_PHYSCOOR(
        CHF_FRA1(a_dest, a_destComp),
        CHF_CONST_INT(a_mu),
        CHF_CONST_REAL(dXimu),
        CHF_BOX(destBox),
        CHF_CONST_INT(destBoxType));
}


// -----------------------------------------------------------------------------
void
CartesianMap::fill_physCoor(FArrayBox&      a_dest,
                            const RealVect& a_dXi,
                            const RealVect  a_scale) const
{
    const RealVect scaledDXi   = a_scale * a_dXi;
    const Box&     destBox     = a_dest.box();
    const IntVect& destBoxType = destBox.type();

    CH_assert(a_dest.nComp() == SpaceDim);

    FORT_CARTESIAN_FILL_PHYSCOOR_ALL_COMPS(
        CHF_FRA(a_dest),
        CHF_CONST_REALVECT(scaledDXi),
        CHF_BOX(destBox),
        CHF_CONST_INTVECT(destBoxType));
}


// -----------------------------------------------------------------------------
void
CartesianMap::fill_dxdXi(FArrayBox&      a_dest,
                         const int       a_destComp,
                         const int       a_mu   [[maybe_unused]],
                         const RealVect& a_dXi  [[maybe_unused]],
                         const Real      a_scale) const
{
    // Sanity checks
    CH_assert(a_dest.interval().contains(a_destComp));
    CH_assert((0 <= a_mu) && (a_mu < SpaceDim));

    a_dest.setVal(a_scale, a_destComp);
}


// -----------------------------------------------------------------------------
void
CartesianMap::fill_J(FArrayBox&      a_dest,
                     const int       a_destComp,
                     const RealVect& /*a_dXi*/,
                     const Real      a_scale) const
{
    // Sanity check
    CH_assert(a_dest.interval().contains(a_destComp));

    a_dest.setVal(a_scale, a_destComp);
}


// -----------------------------------------------------------------------------
void
CartesianMap::fill_dXidx(FArrayBox&      a_dest,
                         const int       a_destComp,
                         const int       a_mu  [[maybe_unused]],
                         const RealVect& a_dXi [[maybe_unused]],
                         const Real      a_scale) const
{
    // Sanity checks
    CH_assert(a_dest.interval().contains(a_destComp));
    CH_assert((0 <= a_mu) && (a_mu < SpaceDim));

    a_dest.setVal(a_scale, a_destComp);
}


// -----------------------------------------------------------------------------
void
CartesianMap::fill_Jinv(FArrayBox&      a_dest,
                        const int       a_destComp,
                        const RealVect& a_dXi [[maybe_unused]],
                        const Real      a_scale) const
{
    // Sanity check
    CH_assert(a_dest.interval().contains(a_destComp));

    a_dest.setVal(a_scale, a_destComp);
}


// -----------------------------------------------------------------------------
void
CartesianMap::fill_gdn(FArrayBox&      a_dest,
                       const int       a_destComp,
                       const int       a_mu  [[maybe_unused]],
                       const RealVect& a_dXi [[maybe_unused]],
                       const Real      a_scale) const
{
    // Sanity checks
    CH_assert(a_dest.interval().contains(a_destComp));
    CH_assert((0 <= a_mu) && (a_mu < SpaceDim));

    a_dest.setVal(a_scale, a_destComp);
}


// -----------------------------------------------------------------------------
void
CartesianMap::fill_gup(FArrayBox&      a_dest,
                       const int       a_destComp,
                       const int       a_mu  [[maybe_unused]],
                       const RealVect& a_dXi [[maybe_unused]],
                       const Real      a_scale) const
{
    // Sanity checks
    CH_assert(a_dest.interval().contains(a_destComp));
    CH_assert((0 <= a_mu) && (a_mu < SpaceDim));

    a_dest.setVal(a_scale, a_destComp);
}


// -----------------------------------------------------------------------------
void
CartesianMap::fill_Jgup(FArrayBox&      a_dest,
                        const int       a_destComp,
                        const int       a_mu  [[maybe_unused]],
                        const RealVect& a_dXi [[maybe_unused]],
                        const Real      a_scale) const
{
    // Sanity checks
    CH_assert(a_dest.interval().contains(a_destComp));
    CH_assert((0 <= a_mu) && (a_mu < SpaceDim));

    a_dest.setVal(a_scale, a_destComp);
}
