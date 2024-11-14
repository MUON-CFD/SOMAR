/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2018
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
#include "PARKCoeffs.H"


// ----------------------------------------------------------------------------
// ForwardEuler_Coeffs
// ----------------------------------------------------------------------------
const Real ForwardEuler_Coeffs::c[] ;

const Real ForwardEuler_Coeffs::aE[numStages][numStages] ;

const Real ForwardEuler_Coeffs::aI[numStages][numStages] ;

const Real ForwardEuler_Coeffs::bE[numStages];
const Real ForwardEuler_Coeffs::bI[numStages];

const Real ForwardEuler_Coeffs::bEhat[numStages];
const Real ForwardEuler_Coeffs::bIhat[numStages] ;

const Real ForwardEuler_Coeffs::bstar[numDenseCoeffs][numStages] ;


// ----------------------------------------------------------------------------
// Midpoint Method
// ----------------------------------------------------------------------------
const Real Midpoint_Coeffs::c[] ;

const Real Midpoint_Coeffs::aE[numStages][numStages] ;

const Real Midpoint_Coeffs::aI[numStages][numStages] ;

const Real Midpoint_Coeffs::bE[numStages];
const Real Midpoint_Coeffs::bI[numStages] ;

const Real Midpoint_Coeffs::bEhat[numStages] ;
const Real Midpoint_Coeffs::bIhat[numStages] ;

const Real Midpoint_Coeffs::bstar[numDenseCoeffs][numStages] ;


// ----------------------------------------------------------------------------
// IMEXRKCB2
// ----------------------------------------------------------------------------
const Real IMEXRKCB2_Coeffs::c[] ;

const Real IMEXRKCB2_Coeffs::aE[numStages][numStages] ;

const Real IMEXRKCB2_Coeffs::aI[numStages][numStages] ;

const Real IMEXRKCB2_Coeffs::bE[numStages] ;
const Real IMEXRKCB2_Coeffs::bI[numStages];

const Real IMEXRKCB2_Coeffs::bEhat[numStages] ;
const Real IMEXRKCB2_Coeffs::bIhat[numStages] ;

const Real IMEXRKCB2_Coeffs::bstar[numDenseCoeffs][numStages] ;


// ----------------------------------------------------------------------------
// IMEXRKCB3b (Might work well with gravity theta method)
// ----------------------------------------------------------------------------
const Real IMEXRKCB3b_Coeffs::c[] ;

const Real IMEXRKCB3b_Coeffs::aE[numStages][numStages];

const Real IMEXRKCB3b_Coeffs::aI[numStages][numStages];

const Real IMEXRKCB3b_Coeffs::bE[numStages] ;
const Real IMEXRKCB3b_Coeffs::bI[numStages];

const Real IMEXRKCB3b_Coeffs::bstar[numDenseCoeffs][numStages] ;


// ----------------------------------------------------------------------------
// IMEXRKCB3c
// ----------------------------------------------------------------------------
const Real IMEXRKCB3c_Coeffs::c[];

const Real IMEXRKCB3c_Coeffs::aE[numStages][numStages];

const Real IMEXRKCB3c_Coeffs::aI[numStages][numStages] ;

const Real IMEXRKCB3c_Coeffs::bE[numStages];
const Real IMEXRKCB3c_Coeffs::bI[numStages];

const Real IMEXRKCB3c_Coeffs::bEhat[numStages];
const Real IMEXRKCB3c_Coeffs::bIhat[numStages];

const Real IMEXRKCB3c_Coeffs::bstar[numDenseCoeffs][numStages];


// ----------------------------------------------------------------------------
// IMEXRKCB3e
// ----------------------------------------------------------------------------
const Real IMEXRKCB3e_Coeffs::c[];

const Real IMEXRKCB3e_Coeffs::aE[numStages][numStages];

const Real IMEXRKCB3e_Coeffs::aI[numStages][numStages];

const Real IMEXRKCB3e_Coeffs::bE[numStages];
const Real IMEXRKCB3e_Coeffs::bI[numStages];

const Real IMEXRKCB3e_Coeffs::bstar[numDenseCoeffs][numStages];


// ----------------------------------------------------------------------------
// IMEXRKCB3f
// ----------------------------------------------------------------------------
const Real IMEXRKCB3f_Coeffs::c[] ;

const Real IMEXRKCB3f_Coeffs::aE[numStages][numStages] ;

const Real IMEXRKCB3f_Coeffs::aI[numStages][numStages];

const Real IMEXRKCB3f_Coeffs::bE[numStages];
const Real IMEXRKCB3f_Coeffs::bI[numStages];

const Real IMEXRKCB3f_Coeffs::bEhat[numStages] ;
const Real IMEXRKCB3f_Coeffs::bIhat[numStages] ;

const Real IMEXRKCB3f_Coeffs::bstar[numDenseCoeffs][numStages] ;



// ----------------------------------------------------------------------------
// CN_RKW3
// ----------------------------------------------------------------------------
const Real CN_RKW3_Coeffs::c[] ;

const Real CN_RKW3_Coeffs::aE[numStages][numStages];

const Real CN_RKW3_Coeffs::aI[numStages][numStages];

const Real CN_RKW3_Coeffs::bE[numStages];

const Real CN_RKW3_Coeffs::bI[numStages] ;


const Real CN_RKW3_Coeffs::bEhat[numStages];

const Real CN_RKW3_Coeffs::bIhat[numStages];


const Real CN_RKW3_Coeffs::bstar[numDenseCoeffs][numStages];

const Real CN_RKW3_Coeffs::bIstar[numDenseCoeffs][numStages] ;


// ----------------------------------------------------------------------------
// Classical RK4
// ----------------------------------------------------------------------------
const Real RK4_Coeffs::c[] ;

const Real RK4_Coeffs::aE[numStages][numStages];

const Real RK4_Coeffs::aI[numStages][numStages] ;

const Real RK4_Coeffs::bE[numStages];
const Real RK4_Coeffs::bI[numStages];

const Real RK4_Coeffs::bEhat[numStages];

const Real RK4_Coeffs::bIhat[numStages];

const Real RK4_Coeffs::bstar[numDenseCoeffs][numStages] ;


// ----------------------------------------------------------------------------
// ARK3_2_4L_2_SA
// ----------------------------------------------------------------------------
const Real ARK3_2_4L_2_SA_Coeffs::c[];

const Real ARK3_2_4L_2_SA_Coeffs::aE[numStages][numStages];

const Real ARK3_2_4L_2_SA_Coeffs::aI[numStages][numStages];


const Real ARK3_2_4L_2_SA_Coeffs::bE[numStages];
const Real ARK3_2_4L_2_SA_Coeffs::bI[numStages];

const Real ARK3_2_4L_2_SA_Coeffs::bEhat[numStages] ;
const Real ARK3_2_4L_2_SA_Coeffs::bIhat[numStages];

const Real ARK3_2_4L_2_SA_Coeffs::bstar[numDenseCoeffs][numStages] ;


// -----------------------------------------------------------------------------
// ARK4_3_6L_2_SA
// ----------------------------------------------------------------------------
const Real ARK4_3_6L_2_SA_Coeffs::c[] ;

const Real ARK4_3_6L_2_SA_Coeffs::aE[numStages][numStages] ;

const Real ARK4_3_6L_2_SA_Coeffs::aI[numStages][numStages];

const Real ARK4_3_6L_2_SA_Coeffs::bE[numStages];
const Real ARK4_3_6L_2_SA_Coeffs::bI[numStages] ;

const Real ARK4_3_6L_2_SA_Coeffs::bEhat[numStages] ;
const Real ARK4_3_6L_2_SA_Coeffs::bIhat[numStages] ;

const Real ARK4_3_6L_2_SA_Coeffs::bstar[numDenseCoeffs][numStages] ;

