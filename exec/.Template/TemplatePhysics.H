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
#ifndef ___T_USER_PHYSICS_NAME_H__INCLUDED___
#define ___T_USER_PHYSICS_NAME_H__INCLUDED___

#include "AMRNSLevel.H"


class T_USER_PHYSICS_NAME: public AMRNSLevel
{
public:
    // Constructor
    T_USER_PHYSICS_NAME ();

    // Virtual destructor for good measure.
    virtual ~T_USER_PHYSICS_NAME ();


    // Custom scalars ----------------------------------------------------------
    // Custom scalars will automatically be advected and, if requested,
    // diffused. If your needs are more specialized than that, then it
    // is up to you to define a forcing function (see below).

    // This tells the Navier-Stokes class how many additional scalars
    // to create. Do not include T, S, or eddyNu. This is inline so that it is
    // just as fast as reading from a const int.
    virtual int
    numScalars() const
    {
        return 0;
    }

    // Returns the names of the scalars. This will be used to put scalar names
    // in the plot files.
    virtual std::string
    getScalarName(const int a_comp) const;


    // ICs ---------------------------------------------------------------------
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
    virtual void
    setICs(State& a_state);

protected:

};


#endif //!___T_USER_PHYSICS_NAME_H__INCLUDED___
