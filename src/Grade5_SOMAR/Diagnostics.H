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
#ifndef ___Diagnostics_H__INCLUDED___
#define ___Diagnostics_H__INCLUDED___
#include "LevelData.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "LevelGeometry.H"
#include "StaggeredFluxLD.H"

namespace Diagnostics{

void 
N2(LevelData<FArrayBox> & a_N2, 
            const LevelData<FArrayBox> & a_b,
            const LevelGeometry* const a_levGeoPtr);
void 
S2(LevelData<FArrayBox> & a_S2, const LevelData<FluxBox> & a_vel, 
            const LevelGeometry* const a_levGeoPtr); 

void 
Ri(      LevelData<FArrayBox>& a_Ri, 
            const LevelData<FluxBox> & a_vel, 
            const LevelData<FArrayBox> & a_b,
            const LevelGeometry* const a_levGeoPtr);
RealVect
    calculatePressureForce(const LevelData<FArrayBox>& a_press,
                           const LevelData<FArrayBox>& a_buoyancy,
                           const RealVect&             a_accel);

RealVect
    calculateStressForce(const StaggeredFluxLD& a_S);

    

};

#endif //___Diagnostics_H__INCLUDED___

