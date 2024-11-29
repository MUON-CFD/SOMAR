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
 *  https://github.com/MUON-CFD/somar.
 ******************************************************************************/

// this structure is used to test if a class T defines the a linear operator
// over a Hilber space suitable for use in an iterative solver.
// isLinearOperatorOverHilbert<T,H>::value returns true if T has the right
// functions (right signature and all);
#pragma once
// Note the use of pointers below, rather than the object itself. This is so
// that isLinearOperatorOverHilbert<T,H>::value will work even if T and/or H are
// abstract classes and/or have private default constructors.
#include "isLinearOperatorOverHilbert.hpp"
namespace Elliptic {

namespace traits {
template <typename LinearOperator, typename Hilbert, int M = 1000>
struct isLinearOperatorOverHilbertMGCompatible {
private:
  static Hilbert *H;
  static const Hilbert *Hc;
  static const Real r;
  static const bool b;
  static const int i;
  static IntVect Iv;
  static LinearOperator *L;
  typedef char no[M];

  // check for restrictResidual (void)
  template <typename C>
  static auto test_for_restrictResidual(C *a)
      -> decltype(std::declval<C>().restrictResidual(*H, *Hc, r, Iv));
  template <typename C> static no &test_for_restrictResidual(...);
  // check for prolongIncrement (void)
  template <typename C>
  static auto test_for_prolongIncrement(C *a)
      -> decltype(std::declval<C>().prolongIncrement(*H, *H, r, Iv, i));
  template <typename C> static no &test_for_prolongIncrement(...);
  // check for createCoarsened (void)
  template <typename C>
  static auto test_for_createCoarsened(C *a)
      -> decltype(std::declval<C>().createCoarsened(*H, *Hc, *H, *Hc));
  template <typename C> static no &test_for_createCoarsened(...);

  LinearOperator *t;

public:
  static const bool value =
      std::is_same<decltype(test_for_restrictResidual<LinearOperator>(t)),
                   void>::value &&
      std::is_same<decltype(test_for_prolongIncrement<LinearOperator>(t)),
                   void>::value &&
      std::is_same<decltype(test_for_createCoarsened<LinearOperator>(t)),
                   void>::value &&
      Elliptic::traits::isLinearOperatorOverHilbert<LinearOperator,
                                                    Hilbert>::value;
};
} // namespace traits
} // namespace Elliptic
