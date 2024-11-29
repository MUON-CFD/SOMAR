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
#ifndef ___EllipticStateOps_H__INCLUDED___
#define ___EllipticStateOps_H__INCLUDED___

#include "REAL.H"

namespace Elliptic {


class NullType {};


// You need to specialize this template for your data holder.
// If your data holder does not use Copiers, leave it NullType and ignore the
// copiers that are passed into some of your operator's functions.
template<class T>
struct StateTraits {
    typedef NullType GridsType;
    typedef NullType CopierType;
};


// This must be specialized to your data holder.
// See LDFABOps for details.
template <class T, class TTraits = StateTraits<T>>
class StateOps {};



}; // namespace Elliptic
#endif //___EllipticStateOps_H__INCLUDED___
