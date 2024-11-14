/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2014 University of North Carolina at Chapel Hill
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
#include "CubicSpline.H"
#include "CubicSplineF_F.H"
#include "Debug.H"


// -----------------------------------------------------------------------------
// Default constructor
// -----------------------------------------------------------------------------
CubicSpline::CubicSpline ()
{;}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
CubicSpline::~CubicSpline ()
{
    this->clear();
}


// -----------------------------------------------------------------------------
// Free memory and make object unusable.
// -----------------------------------------------------------------------------
void
CubicSpline::clear ()
{
    m_x.clear();
    m_f.clear();
    m_d2f.clear();
}


// -----------------------------------------------------------------------------
// Computes the *natural* spline coefficients from the nodal data, (a_x, a_f)
// if a_lofbc and a_hifbc are ommitted. Otherwise, we will use those
// parameters as Neumann BC values.
// -----------------------------------------------------------------------------
void
CubicSpline::solve(const std::vector<Real>& a_f,
                   const std::vector<Real>& a_x,
                   const Real               a_lofbc,
                   const Real               a_hifbc)
{
    // Start clean.
    this->clear();

    // Make sure input vectors are the same size
    CH_assert(a_x.size() == a_f.size());

    // Make sure the a_x vector is sorted and ascending.
    CH_verify(std::is_sorted(a_x.begin(), a_x.end()));

    // Create our own copies of the input data.
    m_x = a_x;
    m_f = a_f;

    // Create space for the second derivatives
    m_d2f.resize(a_f.size());

    // Create workspace for the solver
    std::vector<Real> workspace(a_f.size());

    // Solve for the second derivatives
    FORT_CUBICSPLINE_SOLVESECONDDERIV (
        CHF_VR(m_d2f),
        CHF_CONST_VR(m_x),
        CHF_CONST_VR(m_f),
        CHF_CONST_REAL(a_lofbc),
        CHF_CONST_REAL(a_hifbc),
        CHF_VR(workspace));
}


// -----------------------------------------------------------------------------
// Fills a_f with interpolated values at positions indicated by a_x.
// -----------------------------------------------------------------------------
void
CubicSpline::interp(std::vector<Real>&       a_f,
                    const std::vector<Real>& a_x,
                    const bool               a_linearInterp) const
{
    // Make sure the solve has been called.
    CH_assert(m_d2f.size() > 0);

    // Make sure the output vector has room for the results.
    const int Nax = a_x.size();
    a_f.resize(Nax);

    // Loop over a_x and compute f(x).
    Real x, xlo, xhi, dx;
    Real A, B, C, D;
    int klo, khi, k;
    const int Nmx = m_x.size();

    for (int i = 0; i < Nax; ++i) {
        // Search for x(i) in m_x.
        x = a_x[i];
        klo = 0;
        khi = Nmx-1;

        while (khi-klo > 1) {
            k = (khi + klo) >> 1;
            if (m_x[k] > x) {
                khi = k;
            } else {
                klo = k;
            }
        }

        // Compute cell extents and width.
        xlo = m_x[klo];
        xhi = m_x[khi];
        dx = xhi - xlo;
        CH_assert(dx != 0.0);

        // Predict value using linear interpolant.
        A = (xhi - x) / dx;
        B = (x - xlo) / dx;
        a_f[i] = A*m_f[klo] + B*m_f[khi];

        // Apply cubic correction.
        if (!a_linearInterp) {
            C = A*(A*A-1.0);
            D = B*(B*B-1.0);
            a_f[i] += (C*m_d2f[klo] + D*m_d2f[khi]) * dx*dx/6.0;
        }
    }
}


// -----------------------------------------------------------------------------
// Returns the interpolated value at the position indicated by a_x.
// -----------------------------------------------------------------------------
Real
CubicSpline::interp (const Real& a_x) const
{
    // Make sure the solve has been called.
    CH_assert(m_d2f.size() > 0);

    // Just call the vector version.
    std::vector<Real> xVec(1, a_x);
    std::vector<Real> fVec(1);
    this->interp(fVec, xVec);

    // Return the one and only result.
    return fVec[0];
}


// -----------------------------------------------------------------------------
// Fills a_df with the interpolated first derivatives at positions
// indicated by a_x.
// -----------------------------------------------------------------------------
void
CubicSpline::interpFirstDeriv(std::vector<Real>&       a_df,
                              const std::vector<Real>& a_x) const
{
    // Make sure the solve has been called.
    CH_assert(m_d2f.size() > 0);

    // Make sure the output vector has room for the results.
    const int Nax = a_x.size();
    a_df.resize(Nax);

    // Loop over a_x and compute dx/dx.
    Real x, xlo, xhi, dx, df;
    Real A, B;
    int klo, khi, k;
    const int Nmx = m_x.size();

    for (int i = 0; i < Nax; ++i) {
        // Search for x(i) in m_x.
        x = a_x[i];
        klo = 0;
        khi = Nmx-1;

        while (khi-klo > 1) {
            k = (khi + klo) >> 1;
            if (m_x[k] > x) {
                khi = k;
            } else {
                klo = k;
            }
        }

        // Compute cell extents and width.
        xlo = m_x[klo];
        xhi = m_x[khi];
        dx = xhi - xlo;
        CH_assert(dx != 0.0);

        // Compute finite difference approximation to the derivative.
        df = m_f[khi] - m_f[klo];
        a_df[i] = df/dx;

        // Apply higher-order corrections
        A = (xhi - x) / dx;
        B = (x - xlo) / dx;
        a_df[i] += (- (3.0*A*A-1.0) * m_d2f[klo]
                    + (3.0*B*B-1.0) * m_d2f[khi]) * (dx / 6.0);
    }
}


// -----------------------------------------------------------------------------
Real
CubicSpline::interpFirstDeriv(const Real& a_x) const
{
    // Make sure the solve has been called.
    CH_assert(m_d2f.size() > 0);

    // Just call the vector version.
    std::vector<Real> xVec(1, a_x);
    std::vector<Real> fVec(1);
    this->interpFirstDeriv(fVec, xVec);

    // Return the one and only result.
    return fVec[0];
}


// -----------------------------------------------------------------------------
// Fills a_d2f with the interpolated second derivatives at positions
// indicated by a_x.
// -----------------------------------------------------------------------------
void
CubicSpline::interpSecondDeriv(std::vector<Real>&       a_d2f,
                               const std::vector<Real>& a_x) const
{
    // Make sure the solve has been called.
    CH_assert(m_d2f.size() > 0);

    // Make sure the output vector has room for the results.
    const int Nax = a_x.size();
    a_d2f.resize(Nax);

    // Loop over a_x and compute second derivative.
    Real x, xlo, xhi, dx;
    Real A, B;
    int klo, khi, k;
    const int Nmx = m_x.size();

    for (int i = 0; i < Nax; ++i) {
        // Search for x(i) in m_x.
        x = a_x[i];
        klo = 0;
        khi = Nmx-1;

        while (khi-klo > 1) {
            k = (khi + klo) >> 1;
            if (m_x[k] > x) {
                khi = k;
            } else {
                klo = k;
            }
        }

        // Compute cell extents and width.
        xlo = m_x[klo];
        xhi = m_x[khi];
        dx = xhi - xlo;
        CH_assert(dx != 0.0);

        // Second derivative is just a linear combination of values at nodes.
        A = (xhi - x) / dx;
        B = (x - xlo) / dx;
        a_d2f[i] = A * m_d2f[klo] + B * m_d2f[khi];
    }
}


// -----------------------------------------------------------------------------
// Fills a_d2f with the interpolated second derivatives at the position
// indicated by a_x
// -----------------------------------------------------------------------------
Real
CubicSpline::interpSecondDeriv (const Real& a_x) const
{
    // Make sure the solve has been called.
    CH_assert(m_d2f.size() > 0);

    // Just call the vector version.
    std::vector<Real> xVec(1, a_x);
    std::vector<Real> fVec(1);
    this->interpSecondDeriv(fVec, xVec);

    // Return the one and only result.
    return fVec[0];
}


// -----------------------------------------------------------------------------
// Uses a pre-computed set of data.
// -----------------------------------------------------------------------------
void
CubicSpline::useSolution(const std::vector<Real>& a_x,
                         const std::vector<Real>& a_f,
                         const std::vector<Real>& a_d2f)
{
    CH_assert(a_f.size() == a_x.size());
    CH_assert(a_d2f.size() == a_x.size());

    m_x   = a_x;
    m_f   = a_f;
    m_d2f = a_d2f;
}
