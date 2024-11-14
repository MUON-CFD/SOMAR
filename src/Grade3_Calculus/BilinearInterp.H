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
#ifndef BilinearInterp_H__INCLUDED__
#define BilinearInterp_H__INCLUDED__

#include <vector>
#include "FArrayBox.H"


///
void
LinearInterp(FArrayBox&                                     a_destfFAB,
             const std::array<std::vector<Real>, SpaceDim>& a_destx,
             const Box&                                     a_destBox,
             const FArrayBox&                               a_srcfFAB,
             const std::array<std::vector<Real>, SpaceDim>& a_srcx,
             const Box&                                     a_srcBox,
             const int                                      a_srcExtrapRadius = 1);









/// \brief Simple bilinear interpolation utility on a FAB.
///
/// \param[out] a_destf
/// \param[ in] a_destx      Must be *exactly* the same length as a_destBox.
/// \param[ in] a_desty      Must be *exactly* the same length as a_destBox.
/// \param[ in] a_destBox
/// \param[ in] a_srcf       The data to be interpolated.
/// \param[ in] a_srcx       Must be *exactly* the same length as a_srcBox.
/// \param[ in] a_srcy       Must be *exactly* the same length as a_srcBox.
/// \param[ in] a_srcBox     This is where we are allowed to use a_f's data.
void
BilinearInterp2D(FArrayBox&               a_destf,
                 const std::vector<Real>& a_destx,
                 const std::vector<Real>& a_desty,
                 const Box&               a_destBox,
                 const FArrayBox&         a_srcf,
                 std::vector<Real>        a_srcx,
                 std::vector<Real>        a_srcy,
                 Box                      a_srcBox);


/// \brief Simple bilinear interpolation utility on a FAB.
///
/// \param[out] a_fInterp
/// \param[ in] a_xInterp
/// \param[ in] a_yInterp
/// \param[ in] a_interpBox
/// \param[ in] a_xdir
/// \param[ in] a_ydir
/// \param[ in] a_x          The x-coords of a_f's cell-centers.
/// \param[ in] a_y          The y-coords of a_f's cell-centers.
/// \param[ in] a_f          The data to be interpolated.
///
/// \warning a_f must have exactly the same dimensions as a_x (X) a_y,
///          where (X) is an outer product.
void
BilinearInterp2D(FArrayBox&               a_fInterp,
                 const FArrayBox&         a_xInterp,
                 const FArrayBox&         a_yInterp,
                 const Box&               a_interpBox,
                 const int                a_xdir,
                 const int                a_ydir,
                 const FArrayBox&         a_f,
                 const std::vector<Real>& a_x,
                 const std::vector<Real>& a_y);


#endif //!BilinearInterp_H__INCLUDED__
