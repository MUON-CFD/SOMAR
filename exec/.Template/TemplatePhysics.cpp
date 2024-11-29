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
#include "T_USER_PHYSICS_NAME.H"
#include "SetValLevel.H"
#include "Debug.H"


//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
T_USER_PHYSICS_NAME::T_USER_PHYSICS_NAME()
: AMRNSLevel::AMRNSLevel()
{
}


//------------------------------------------------------------------------------
// Virtual destructor for good measure.
//------------------------------------------------------------------------------
T_USER_PHYSICS_NAME::~T_USER_PHYSICS_NAME()
{
}


// -----------------------------------------------------------------------------
// Returns the names of the scalars. This will be used to put scalar names
// in the plot files.
// -----------------------------------------------------------------------------
std::string
T_USER_PHYSICS_NAME::getScalarName(const int a_comp) const
{
    switch (a_comp) {
        // If you have scalars, add cases here.
        default:
            MayDay::Error(
                "T_USER_PHYSICS_NAME::getScalarName: "
                "a_comp out of range.");
    }

    return "Undefined";
}


//------------------------------------------------------------------------------
// Set the ICs on all fields at once. This includes the velocity,
// total temperature, total salinity, eddy viscosity and diffusivity, all
// user-defined scalars, and pressure.
//
// By default, everything is set to zero except T and S which are set to
// thier background stratification values (T' = S' = 0).
//
// Also, you do not need to specify the pressure. If rhs.computeInitPressure
// = 1, SOMAR will automatically figure out the pressure. However, if you
// have a good guess, SOMAR can use that as an initial guess in the Poisson
// solver.
//
// By the time this function is called, *m_levGeoPtr will be ready for use.
//------------------------------------------------------------------------------
void
T_USER_PHYSICS_NAME::setICs(State& a_state)
{
}
