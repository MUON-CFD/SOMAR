/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2017
 *    Jefferson (Philadelphia University + Thomas Jefferson University) and
 *    University of North Carolina at Chapel Hill
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
#ifndef ___FABAlgebra_H__INCLUDED___
#define ___FABAlgebra_H__INCLUDED___

#include "FArrayBox.H"


class FABAlgebra
{
public:
    // Calculates dest = dest + p0 * p1
    static void
    AddProd2(FArrayBox&       a_dest,
             const int        a_destComp,
             const FArrayBox& a_p0,
             const int        a_p0Comp,
             const FArrayBox& a_p1,
             const int        a_p1Comp,
             const Box&       a_destBox);

    // Calculates dest = dest - p0 * p1
    static void
    SubProd2(FArrayBox&       a_dest,
             const int        a_destComp,
             const FArrayBox& a_p0,
             const int        a_p0Comp,
             const FArrayBox& a_p1,
             const int        a_p1Comp,
             const Box&       a_destBox);

    // Just like FArrayBox::axby, except this works on the intersection of
    // all the FABs.
    static void
    axby(FArrayBox&       a_dest,
         const FArrayBox& a_x,
         const FArrayBox& a_y,
         const Real       a_a,
         const Real       a_b);

    // Calculates dest = dest^pow
    static void
    pow(FArrayBox& a_dest,
        const Real a_pow);

    // Calculates dest = dest^pow
    static void
    pow(FArrayBox& a_dest,
        const int  a_destComp,
        const Box& a_destBox,
        const Real a_pow);

    // Computes EC = EC * (alpha + beta*Av[CC])
    static void
    ECmultCC(FArrayBox&       a_ecFAB,
             const int        a_ecFABComp,
             const Box&       a_ecBox,
             const FArrayBox& a_ccFAB,
             const int        a_ccFABComp,
             const Real       a_alpha      = 0.0,
             const Real       a_beta       = 1.0);

    // Computes FC = FC * (alpha + beta*Av[CC])
    static void
    FCmultCC(FArrayBox&       a_fcFAB,
             const int        a_fcFABComp,
             const Box&       a_fcBox,
             const FArrayBox& a_ccFAB,
             const int        a_ccFABComp,
             const Real       a_alpha      = 0.0,
             const Real       a_beta       = 1.0);

    // Computes CC1 = CC1 * (alpha + beta*CC2)
    static void
    CCmultCC(FArrayBox&       a_cc1FAB,
             const int        a_cc1FABComp,
             const Box&       a_region,
             const FArrayBox& a_cc2FAB,
             const int        a_cc2FABComp,
             const Real       a_alpha      = 0.0,
             const Real       a_beta       = 1.0);

    // Computes FC = FC + scale * Av[CC].
    static void
    FCaddCC(FArrayBox&       a_fcFAB,
            const int        a_fcFABComp,
            const Box&       a_fcBox,
            const FArrayBox& a_ccFAB,
            const int        a_ccFABComp,
            const Real       a_scale = 1.0);

    // Computes CC = CC + scale * Av[FC].
    static void
    CCaddFC(FArrayBox&       a_ccFAB,
            const int        a_ccFABComp,
            const Box&       a_ccBox,
            const FArrayBox& a_fcFAB,
            const int        a_fcFABComp,
            const Real       a_scale = 1.0);
};


#endif //!___FABAlgebra_H__INCLUDED___
