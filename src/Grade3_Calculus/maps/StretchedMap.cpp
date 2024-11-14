#include "StretchedMap.H"
#include "SOMAR_Constants.H"
#include "Debug.H"


// -----------------------------------------------------------------------------
StretchedMap::StretchedMap(const RealVect& a_xmin,
                           const RealVect& a_xmax,
                           const RealVect& a_ampl)
: m_xmin(a_xmin)
, m_xmax(a_xmax)
, m_ampl(a_ampl)
, m_k(2.0 * Pi / (a_xmax - a_xmin))
{
    CH_verify(a_xmin < a_xmax);
}


// -----------------------------------------------------------------------------
StretchedMap::~StretchedMap()
{
}


// -----------------------------------------------------------------------------
void
StretchedMap::interp(Vector<Real>&       a_x,
                     const Vector<Real>& a_xi,
                     const int           a_mu) const
{
    const Real A  = m_ampl[a_mu];
    const Real k  = m_k[a_mu];
    const Real x0 = m_xmin[a_mu];

    for (size_t idx = 0; idx < a_x.size(); ++idx) {
        a_x[idx] = a_xi[idx] + A * sin(k * (a_xi[idx] - x0));
    }
}
