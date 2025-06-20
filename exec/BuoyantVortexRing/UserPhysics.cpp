/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2022
 *    Jefferson University and University of North Carolina at Chapel Hill
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
#include "UserPhysics.H"
#include "SetValLevel.H"
#include "Debug.H"
#include "ParmParse.H"
#include "PyGlue.H"

Real UserPhysics::R = 1.0;
Real UserPhysics::a = 0.1;
Real UserPhysics::alpha = 0.0;
Real UserPhysics::Gamma = 1.0;
Real UserPhysics::B = 0.0; 
bool UserPhysics::isDefined = false;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
UserPhysics::UserPhysics()
: AMRNSLevel::AMRNSLevel()
{
    if (!UserPhysics::isDefined){
        UserPhysics::isDefined = true;
        ParmParse pp("bvr");
        pp.get("R", UserPhysics::R);
        pp.get("a", UserPhysics::a);
        pp.get("alpha", UserPhysics::alpha);
        pp.get("Gamma", UserPhysics::Gamma);
        pp.get("B", UserPhysics::B);
    }
}


//------------------------------------------------------------------------------
// Virtual destructor for good measure.
//------------------------------------------------------------------------------
UserPhysics::~UserPhysics()
{
}


// -----------------------------------------------------------------------------
// Returns the names of the scalars. This will be used to put scalar names
// in the plot files.
// -----------------------------------------------------------------------------
std::string
UserPhysics::getScalarName(const int a_comp) const
{
    switch (a_comp) {
        // If you have scalars, add cases here.
        default:
            MayDay::Error(
                "UserPhysics::getScalarName: "
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
UserPhysics::setICs(State& a_state)
{
    setValLevel(a_state.q, 0.0);
    setValLevel(a_state.p, 0.0);
    setValLevel(a_state.q, 0.0);
    const DisjointBoxLayout &grids = this->getBoxes();
    DataIterator dit = grids.dataIterator();
    for (dit.reset(); dit.ok();++dit)
    {
        FArrayBox &u = a_state.vel[dit][0];
        FArrayBox &v = a_state.vel[dit][1];
        FArrayBox &w = a_state.vel[dit][2];
        FArrayBox &b = a_state.T[dit];
        FArrayBox Ucoords(u.box(),SpaceDim);
        FArrayBox Vcoords(v.box(),SpaceDim);
        FArrayBox Wcoords(w.box(),SpaceDim);
        FArrayBox Bcoords(b.box(),SpaceDim);
        m_levGeoPtr->fill_physCoor(Ucoords);
        m_levGeoPtr->fill_physCoor(Vcoords);
        m_levGeoPtr->fill_physCoor(Wcoords);
        m_levGeoPtr->fill_physCoor(Bcoords);
        Py::PythonFunction("InitVortex_vect", "initialize", 
                            u,
                            v,
                            w,
                            b,
                            Ucoords,
                            Vcoords,
                            Wcoords,
                            Bcoords, 
                            UserPhysics::R,
                            UserPhysics::a,
                            UserPhysics::Gamma,
                            UserPhysics::alpha);
    }
    LayoutTools::averageOverlappingValidFaces(a_state.vel);
}
