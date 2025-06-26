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
 *  https://github.com/MUON-CFD/SOMAR.
 ******************************************************************************/
#ifndef __CONSTANTS_H__INCLUDED__
#define __CONSTANTS_H__INCLUDED__

#include <limits>
#include <cmath>
#include <algorithm>

#include "Box.H"
#include "REAL.H"


/// NaN
#define _USE_NAN_
#ifndef _USE_NAN_
    static constexpr Real quietNAN = 4.321e100;
    static constexpr Real loudNAN  = 4.321e100;
#else // !_USE_NAN_
    static constexpr Real quietNAN = std::numeric_limits<Real>::quiet_NaN();
    static constexpr Real loudNAN  = std::numeric_limits<Real>::signaling_NaN();
#endif // !_USE_NAN_


/// Flags
static constexpr Real BOGUS_TIME = 1.0e300; // NaN causes issues in VisIt.
static constexpr Real timeEps = 1.0e-8;     // See also: smallReal


/// Mathematical constants
static constexpr Real E     = 2.718281828459045e0;
static constexpr Real SQRT2 = 1.414213562373095e0;
static constexpr Real SQRT3 = 1.732050807568877e0;
static constexpr Real Pi    = 3.14159265358979323846264338327950288e0;
static constexpr int TWOPOWD = D_TERM(2, *2, *2);

/// Limiting values of the Real type.
static constexpr Real minReal = std::numeric_limits<Real>::min();
static constexpr Real maxReal = std::numeric_limits<Real>::max();

/// Numbers smaller than this will be considered zero.
static constexpr Real smallReal = 1.0e4 * std::numeric_limits<Real>::epsilon();

/// Some simple functions for comparing Reals.
namespace RealCmp
{
    /// Returns true if a number is considered zero.
    inline bool isZero (Real a_val)
    {
        return (std::abs(a_val) < smallReal);
    }

    /// Returns true if a number is considered NAN.
    inline bool isNAN (Real a_val)
    {
        return ((a_val != a_val) || (a_val <= -1.0e10)  || (1.0e10 <= a_val));
    }


    /**
     * @brief      Returns if two numbers are considered equal.
     * @detail
     *  Checks if |a1 - a2| <= smallReal * max(|a1|, |a2|).
     *  I do this because the simpler check |a1-a2| <= smallReal can produce
     *  the wrong result if a1 or a2 are large. For example, suppose
     *  a1 = 123456.789012342,
     *  a2 = 123456.789012349, and
     *  smallReal = 1.0e-12.
     *  |a1-a2| = 7.0e-9, which is larger than smallReal. So |a1-a2|<=smallReal
     *  would return false and it would appear a1 != a2. But this is not true
     *  since the two numbers are equal in the first 12 significant digits.
     */
    inline bool eq (Real a_1, Real a_2)
    {
        const Real diff = std::abs(a_1 - a_2);
        const Real sm = smallReal * std::max(std::abs(a_1), std::abs(a_2));
        return diff <= sm;
    }

    /// Returns if two numbers are considered unequal.
    inline bool neq (Real a_1, Real a_2)
    {
        return !eq(a_1, a_2);
    }
}

// The canonical empty box.
const Box EMPTYBOX;


#endif //!__CONSTANTS_H__INCLUDED__

