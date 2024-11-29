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
#ifndef ___ProblemContext_H__INCLUDED___
#define ___ProblemContext_H__INCLUDED___

#include <memory>
#include "Registry.H"

#include "BaseParameters.H"
#include "TimeParameters.H"
#include "OutputParameters.H"
#include "AMRParameters.H"
#include "RHSParameters.H"
#include "ProjectorParameters.H"

#include "Debug.H"


// This is an easy way to set all parameters for the run (once) and be able
// to access them from anywhere within the code. This is also a Registry, so
// if you want to store an unforseen parameter on-the-fly, you can.
class ProblemContext: public Registry
{
private:
    // This is a singleton class -- private constructor, static pointer!
    ProblemContext();

    static std::unique_ptr<ProblemContext> s_ctx;

public:
    // Destructor
    ~ProblemContext();

    // This returns the single ProblemContext object. (const version)
    static const ProblemContext*
    getInstance();

    // This returns the single ProblemContext object.
    static ProblemContext*
    getNonConstInstance();

    // This deletes the static pointer in hopes valgrind will take notice.
    static void
    freeMemory();

    // The user-accessible members...
    BaseParameters      base;   // Level 0 geometry and discretization.
    TimeParameters      time;   // Timestepping stuff.
    OutputParameters    output; // also needed by AnisotropicAMR. Plot stuff.
    AMRParameters       amr;    // everything needed by AnisotropicAMR.
    RHSParameters       rhs;    // Model parameters.
    ProjectorParameters proj;   // for the projector.
    // ...
};


#endif  //!___ProblemContext_H__INCLUDED___
