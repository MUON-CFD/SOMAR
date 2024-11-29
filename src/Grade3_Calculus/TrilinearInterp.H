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
#ifndef TrilinearInterp_H__INCLUDED__
#define TrilinearInterp_H__INCLUDED__

#include "FArrayBox.H"


// -----------------------------------------------------------------------------
// Simple trilinear interpolation utility.
// -----------------------------------------------------------------------------
void
TrilinearInterp3D(FArrayBox&          a_fInterp,
                  const FArrayBox&    a_xInterp,
                  const FArrayBox&    a_yInterp,
                  const FArrayBox&    a_zInterp,
                  const Box&          a_interpBox,
                  const int           a_xdir,
                  const int           a_ydir,
                  const int           a_zdir,
                  const Vector<Real>& a_x,
                  const Vector<Real>& a_y,
                  const Vector<Real>& a_z,
                  const FArrayBox&    a_f);


#if CH_SPACEDIM == 3
#ifndef NDEBUG
// -----------------------------------------------------------------------------
// A test of TrilinearInterp3D. By refining inBox by 2, we should get an
// error drop of 4. This returns the inf-norm of the error.
// -----------------------------------------------------------------------------
Real
TrilinearInterp3D_Test(std::function<Real(Real, Real, Real)> a_func,
                       const IntVect& a_ref = IntVect::Unit);
#endif // !NDEBUG
#endif // CH_SPACEDIM == 3


#endif //!TrilinearInterp_H__INCLUDED__
