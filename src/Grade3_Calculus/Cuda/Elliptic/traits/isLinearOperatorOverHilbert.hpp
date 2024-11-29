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
// this structure is used to test if a class T defines the a linear operator
// over a Hilber space suitable for use in an iterative solver.
// isLinearOperatorOverHilbert<T,H>::value returns true if T has the right
// functions (right signature and all);
#pragma once
#include "FArrayBox.H"
#include "LevelData.H"
// Note the use of pointers below, rather than the object itself. This is so
// that isLinearOperatorOverHilbert<T,H>::value will work even if T and/or H are
// abstract classes.
namespace Elliptic {
namespace traits {

template <typename LinearOperator, typename Hilbert, int M = 1000>
struct isLinearOperatorOverHilbert {
private:
  static const LevelData<FArrayBox> *LDc;
  static LevelData<FArrayBox> *LD;
  static Hilbert *H;
  static const Hilbert *Hc;
  static const Real r;
  static const bool b;
  static const int i;

  typedef char no[M];
  // check that solve exists (void).
  template <typename C>
  static auto test_for_applyOp(C *a)
      -> decltype(std::declval<C>().applyOp(*H, *H, r, b));
  template <typename C> static no &test_for_applyOp(...);

  // check for preCond (void).
  template <typename C>
  static auto test_for_preCond(C *a)
      -> decltype(std::declval<C>().preCond(*H, *Hc, r, i));
  template <typename C> static no &test_for_preCond(...);

  // check for residual() (void)
  template <typename C>
  static auto test_for_residual(C *a)
      -> decltype(std::declval<C>().residual(*H, *H, *Hc, r, b));
  template <typename C> static no &test_for_residual(...);

  // check for relax() (void)
  template <typename C>
  static auto test_for_relax(C *a)
      -> decltype(std::declval<C>().relax(*H, *Hc, r, i, r));
  template <typename C> static no &test_for_relax(...);

  // check for LDToContainer (void)
  template <typename C>
  static auto test_for_LDToContainer(C *a)
      -> decltype(std::declval<C>().LDToContainer(*LDc, *H));
  template <typename C> static no &test_for_LDToContainer(...);

  // check for ContainerToLD (void)
  template <typename C>
  static auto test_for_ContainerToLD(C *a)
      -> decltype(std::declval<C>().ContainerToLD(*Hc, *LD));
  template <typename C> static no &test_for_ContainerToLD(...);
  LinearOperator *t;

public:
  static const bool value =
      std::is_same<decltype(test_for_applyOp<LinearOperator>(t)),
                   void>::value &&
      std::is_same<decltype(test_for_preCond<LinearOperator>(t)),
                   void>::value &&
      std::is_same<decltype(test_for_residual<LinearOperator>(t)),
                   void>::value &&
      std::is_same<decltype(test_for_relax<LinearOperator>(t)), void>::value &&
      std::is_same<decltype(test_for_LDToContainer<LinearOperator>(t)),
                   void>::value &&
      std::is_same<decltype(test_for_ContainerToLD<LinearOperator>(t)),
                   void>::value;
};
} // namespace traits
} // namespace Elliptic
