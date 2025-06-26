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
#ifndef ___SideArray_H__INCLUDED___
#define ___SideArray_H__INCLUDED___
#include <array>

// -----------------------------------------------------------------------------
// A SideArray contains one int for each of the domain's boundaries.
// In 3D, this is 6 elements -- a lower and upper x-boundary, a lower
// and upper y-boundary, and a lower and upper z-boundary.
// This can be used to tell a BC-setting function which ghosts to fill.
// For example, if we create
//    SideArray sides = {{ {1,1}, {0,0}, {1,0} }};
// to the BCTools::diri function, we would be setting Dirichlet BCs on
// the lower and upper x-boundaries and the lower z-boundary.
// WARNING: Note the extra pair of brackets!!! This issue was fixed in
//  C++14, but not in C++11.
// -----------------------------------------------------------------------------
typedef std::array<std::array<int, 2>, CH_SPACEDIM> SideArray;


// I'm not thrilled that these are in the global namespace. I'd like to have
// SideArray inherit from std::array<std::array<int, 2>, CH_SPACEDIM>
// and these are members. Perhaps another day.
#if CH_SPACEDIM == 2
    static constexpr SideArray  NoSides = {{ {0,0}, {0,0} }};
    static constexpr SideArray AllSides = {{ {1,1}, {1,1} }};
    static constexpr SideArray   XSides = {{ {1,1}, {0,0} }};
    static constexpr SideArray  XYSides = {{ {1,1}, {0,0} }};
    static constexpr SideArray   ZSides = {{ {0,0}, {1,1} }};
    static constexpr SideArray LoZSide  = {{ {0,0}, {1,0} }};
    static constexpr SideArray HiZSide  = {{ {0,0}, {0,1} }};
    static constexpr std::array<SideArray, 2>
                                  SideV = {{ XSides, ZSides }};
#else
    static constexpr SideArray  NoSides = {{ {0,0}, {0,0}, {0,0} }};
    static constexpr SideArray AllSides = {{ {1,1}, {1,1}, {1,1} }};
    static constexpr SideArray   XSides = {{ {1,1}, {0,0}, {0,0} }};
    static constexpr SideArray   YSides = {{ {0,0}, {1,1}, {0,0} }};
    static constexpr SideArray  XYSides = {{ {1,1}, {1,1}, {0,0} }};
    static constexpr SideArray   ZSides = {{ {0,0}, {0,0}, {1,1} }};
    static constexpr SideArray LoZSide  = {{ {0,0}, {0,0}, {1,0} }};
    static constexpr SideArray HiZSide  = {{ {0,0}, {0,0}, {0,1} }};
    static constexpr std::array<SideArray, 3>
                                  SideV = {{ XSides, YSides, ZSides }};
#endif


// -----------------------------------------------------------------------------
/// @brief Checks if a SideArray is properly defined.
/// @param a
/// @return true if all elements are 0 or 1.
constexpr bool
check(const SideArray& a) {
    return (a[0][0] == 0 || a[0][0] == 1)
        && (a[0][1] == 0 || a[0][1] == 1)
        && (a[1][0] == 0 || a[1][0] == 1)
        && (a[1][1] == 0 || a[1][1] == 1)
#if CH_SPACEDIM > 2
        && (a[2][0] == 0 || a[2][0] == 1)
        && (a[2][1] == 0 || a[2][1] == 1)
#endif
    ;
}


// -----------------------------------------------------------------------------
/// @brief Subtraction op. Useful when you need to specify all sides except one.
/// @param a
/// @param b
constexpr SideArray
operator-(const SideArray& a, const SideArray& b) {
    CH_assert(check(a));
    CH_assert(check(b));
    return {{
         {a[0][0] - b[0][0], a[0][1] - b[0][1]}
        ,{a[1][0] - b[1][0], a[1][1] - b[1][1]}
#if CH_SPACEDIM > 2
        ,{a[2][0] - b[2][0], a[2][1] - b[2][1]}
#endif
    }};
}

#endif //!___SideArray_H__INCLUDED___
