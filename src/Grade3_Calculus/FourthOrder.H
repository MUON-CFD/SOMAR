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
 *  https://github.com/UNC-CFD/LES-SOMAR.
 ******************************************************************************/
#ifndef ___FourthOrder_H__INCLUDED___
#define ___FourthOrder_H__INCLUDED___

#include "FArrayBox.H"
#include "FluxBox.H"
#include "LevelData.H"


class FourthOrder
{
public:
    /// \name Finite differerence <--> finite volume
    /// \{

    /// @brief Converts nodal to averaged data in all centered directions.
    /// @param a_dataFAB
    /// @param a_conversionBox
    static void
    nodalToAvg(FArrayBox& a_dataFAB, const Box& a_conversionBox);

    /// @brief Converts averaged to nodal data in all centered directions.
    /// @param a_dataFAB
    /// @param a_conversionBox
    static void
    avgToNodal(FArrayBox& a_dataFAB, const Box& a_conversionBox);

    /// @brief Converts nodal to averaged data in conversionDir.
    /// @param a_dataFAB
    /// @param a_conversionBox
    /// @param a_conversionDir
    /// @param a_scale Leave default. Not for public use.
    static void
    nodalToAvg(FArrayBox& a_dataFAB,
               const Box& a_conversionBox,
               const int  a_conversionDir,
               const Real a_scale = 1.0);

    /// @brief Converts averaged to nodal data in conversionDir.
    /// @param a_dataFAB
    /// @param a_conversionBox
    /// @param a_conversionDir
    static inline void
    avgToNodal(FArrayBox& a_dataFAB,
               const Box& a_conversionBox,
               const int  a_conversionDir)
    {
        FourthOrder::nodalToAvg(
            a_dataFAB, a_conversionBox, a_conversionDir, -1.0);
    }

    /// @brief Converts face averages to face-centered nodal values.
    /// @param a_data
    /// @details Does not exchange or fill ghosts. That is up to the caller.
    ///  This only converts valid faces.
    static void
    avgToNodal(LevelData<FluxBox>& a_data);

    /// @brief Converts face-centered nodal values to face averages.
    /// @param a_data
    /// @details Does not exchange or fill ghosts. That is up to the caller.
    ///  This only converts valid faces.
    static void
    nodalToAvg(LevelData<FluxBox>& a_data);
    /// \}

    /// \name Recentering
    /// \{

    /// @brief Conservative interpolation between nodal and zone-averaged data.
    /// @param a_destFAB  Must have opposite type as a_srcFAB only in a_interpDir.
    /// @param a_destBox  Must have same type as a_destFAB.
    /// @param a_srcFAB   Must have opposite type as a_destFAB only in a_interpDir.
    /// @param a_interpDir
    /// @details This function only interpolates in one direction.
    ///
    /// If srcFAB is zonal in interpDir, then destFAB and destBox must
    /// be nodal in interpDir, have identical centerings in all other dirs,
    /// and srcFAB must at least be defined on enclosedCells(destBox, interpDir).
    ///
    /// If srcFAB is nodal in interpDir, then destFAB and destBox must
    /// be zonal in interpDir, have identical centerings in all other dirs,
    /// and srcFAB must at least be defined on surroundingNodes(destBox, interpDir).
    static void
    conservativeInterp(FArrayBox&       a_destFAB,
                       const int        a_destFABComp,
                       const Box&       a_destBox,
                       const FArrayBox& a_srcFAB,
                       const int        a_srcFABComp,
                       const int        a_interpDir);

    /// @brief Interpolation between nodal and zone-centered data.
    /// @param a_destFAB  Must have opposite type as a_srcFAB only in a_interpDir.
    /// @param a_destBox  Must have same type as a_destFAB.
    /// @param a_srcFAB   Must have opposite type as a_destFAB only in a_interpDir.
    /// @param a_interpDir
    /// @details This function only interpolates in one direction.
    ///
    /// If srcFAB is zonal in interpDir, then destFAB and destBox must be nodal
    /// in interpDir, have identical centerings in all other dirs, and srcFAB
    /// must at least be defined on enclosedCells(destBox, interpDir).grow(interpDir, 1).
    /// NOTE: This is a different requirement than conservativeInterp to avoid extrapolation!
    ///
    /// If srcFAB is nodal in interpDir, then destFAB and destBox must
    /// be zonal in interpDir, have identical centerings in all other dirs,
    /// and srcFAB must at least be defined on surroundingNodes(destBox, interpDir).
    static void
    interp(FArrayBox&       a_destFAB,
           const int        a_destFABComp,
           const Box&       a_destBox,
           const FArrayBox& a_srcFAB,
           const int        a_srcFABComp,
           const int        a_interpDir);
    /// \}


    // /// @brief Given the cell/face/edge averages a_f and a_g, we compute the
    // ///        cell/face/edge average a_f*a_g to fourth order.
    // ///        This stays within a_validPhibox by using biased stencils.
    // /// @param a_fg
    // /// @param a_fgComp
    // /// @param a_f
    // /// @param a_fComp
    // /// @param a_g
    // /// @param a_gComp
    // /// @param a_validBox
    // static void
    // multiplyFourthOrderAverages(FArrayBox&       a_fg,
    //                             const int        a_fgComp,
    //                             const FArrayBox& a_f,
    //                             const int        a_fComp,
    //                             const FArrayBox& a_g,
    //                             const int        a_gComp,
    //                             const Box&       a_validBox);

protected:
    static bool
    transverseSameType(const Box& a_box1, const Box& a_box2, const int a_dir);

    static inline bool
    oppositeType(const Box& a_box1, const Box& a_box2, const int a_dir)
    {
        return a_box1.type(a_dir) != a_box2.type(a_dir);
    }

    static void
    getDestBoxes(Box&       a_lolo,
                 Box&       a_lo,
                 Box&       a_mid,
                 Box&       a_hi,
                 Box&       a_hihi,
                 const Box& a_srcBox,
                 const Box& a_destBox,
                 const int  a_interpDir);
};


#endif  //!___FourthOrder_H__INCLUDED___
