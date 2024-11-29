/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2024 Thomas Jefferson University and Arizona State University
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
#ifndef CubicSpline_H__INCLUDED__
#define CubicSpline_H__INCLUDED__

#include <vector>
#include "REAL.H"


// -----------------------------------------------------------------------------
// Utility to compute and store spline data.
// -----------------------------------------------------------------------------
class CubicSpline
{
public:
    // Default constructor
    CubicSpline ();

    // Destructor
    ~CubicSpline ();

    // Free memory and make object unusable.
    void clear ();

    // Computes the *natural* spline coefficients from the nodal data, (a_x, a_f)
    // if a_lofbc and a_hifbc are ommitted. Otherwise, we will use those
    // parameters as Neumann BC values.
    void
    solve(const std::vector<Real>& a_f,
          const std::vector<Real>& a_x,
          const Real               a_lofbc = 1.0e50,
          const Real               a_hifbc = 1.0e50);

    // Fills a_f with interpolated values at positions indicated by a_x.
    void
    interp(std::vector<Real>&       a_f,
           const std::vector<Real>& a_x,
           const bool               a_linearInterp = false) const;

    // Returns the interpolated value at the position indicated by a_x.
    Real
    interp(const Real& a_x) const;

    // Fills a_df with the interpolated first derivatives at positions
    // indicated by a_x.
    // These inputs must have the same centering.
    void
    interpFirstDeriv(std::vector<Real>&       a_df,
                     const std::vector<Real>& a_x) const;

    // Fills a_df with the interpolated first derivatives at the position
    // indicated by a_x
    Real
    interpFirstDeriv(const Real& a_x) const;

    // Fills a_d2f with the interpolated second derivatives at positions
    // indicated by a_x.
    void
    interpSecondDeriv(std::vector<Real>&       a_d2f,
                      const std::vector<Real>& a_x) const;

    // Fills a_d2f with the interpolated second derivatives at the position
    // indicated by a_x
    Real
    interpSecondDeriv(const Real& a_x) const;

    // Uses a pre-computed set of data.
    void
    useSolution(const std::vector<Real>& a_x,
                const std::vector<Real>& a_f,
                const std::vector<Real>& a_d2f);

protected:
    // The nodal coordinates
    std::vector<Real> m_x;

    // The nodal function values
    std::vector<Real> m_f;

    // The nodal second derivates needed to create the splines.
    std::vector<Real> m_d2f;
};


#endif //!CubicSpline_H__INCLUDED__
