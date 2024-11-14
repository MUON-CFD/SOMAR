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
#ifndef ___UserMain_H__INCLUDED___
#define ___UserMain_H__INCLUDED___

#include "GeoSourceInterface.H"
#include "AMRNSLevel.H"
#include "ProblemContext.H"


// As far as the general user is concerned, this takes the place of the main()
// function. There is a true main() function, but it contains ugly boilerplate
// code that no one should care about.
class UserMain
{
public:
    // As the name suggests, this will be called just one time during the
    // initialization process. MPI will already be setup, so don't worry about
    // those kinds of details. Instead, just focus on tweaking parameters (if
    // needed) and choosing a coordinate system.
    //
    // When creating the geometry, be sure to use "new" keyword to allocate the
    // object on the heap just as I have done in the example code. You don't
    // need to worry about freeing memory. SOMAR will take care of that.
    static GeoSourceInterface*
    oneTimeSetup ();

    // This will be called an unknown number of times, but after oneTimeSetup().
    // Basically, SOMAR doesn't know what user-defined physics class you plan to
    // use or how to set it up. Getting that ready is your job here.
    //
    // When creating your physics class, be sure you use "new" keyword to
    // allocate your object on the heap just as I have done in the example code.
    // Don't need to worry about freeing memory. SOMAR will take care of that.
    static AMRNSLevel*
    createPhysics ();
};


#endif //!___UserMain_H__INCLUDED___
