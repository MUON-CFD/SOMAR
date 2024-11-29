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
 *  https://github.com/MUON-CFD/SOMAR.
 ******************************************************************************/
#ifdef CH_USE_CUDA_PROJECTOR
#include <cusp/multiply.h>
#endif

#include "UserMain.H"
#include "CartesianMap.H"
#include "StretchedMap.H"
#include "ParmParse.H"
#include "LockExchangePhysics.H"


// -----------------------------------------------------------------------------
// As the name suggests, this will be called just one time during the
// initialization process. MPI will already be setup, so don't worry about
// those kinds of details. Instead, just focus on tweaking parameters (if
// needed) and choosing a coordinate system.
//
// When creating the geometry, be sure to use "new" keyword to allocate the
// object on the heap just as I have done in the example code. You don't
// need to worry about freeing memory. SOMAR will take care of that.
// -----------------------------------------------------------------------------
GeoSourceInterface*
UserMain::oneTimeSetup ()
{
    // Gather input file settings related to the geometry.
    ProblemContext* __nowarn_unused ctx = ProblemContext::getNonConstInstance();
    // const RealVect& L    = ctx->base.L;
    // const IntVect&  nx   = ctx->base.nx;
    // const IntVect&  imin = ctx->base.nxOffset;

    // // Compute the domain extents.
    // RealVect dx   = L / RealVect(nx);
    // RealVect xmin = RealVect(imin) * dx;
    // RealVect xmax = xmin + L;

    // // The streching amplitudes.
    // RealVect ampl;
    // ampl[0]            = 2.5;   // Resolve mid-domain, at the initial interface.
    // ampl[1]            = 0.0;
    // ampl[SpaceDim - 1] = -0.25; // Resolve at the top and bottom boundaries.

    // Create the geometry class.
    CartesianMap* geoPtr = new CartesianMap();
    // StretchedMap* geoPtr = new StretchedMap(xmin, xmax, ampl);

    // We are done. Return the geometry class pointer.
    return geoPtr;
}


// -----------------------------------------------------------------------------
// This will be called an unknown number of times, but after oneTimeSetup().
// Basically, SOMAR doesn't know what user-defined physics class you plan to
// use or how to set it up. Getting that ready is your job here.
//
// When creating your physics class, be sure you use "new" keyword to
// allocate your object on the heap just as I have done in the example code.
// Don't need to worry about freeing memory. SOMAR will take care of that.
// -----------------------------------------------------------------------------
AMRNSLevel*
UserMain::createPhysics ()
{
    // Create and return the appropriate user-defined Physics class.
    // Again, use the "new" keyword to allocate and don't worry about cleanup.
    LockExchangePhysics* physPtr = new LockExchangePhysics();
    return physPtr;
}

