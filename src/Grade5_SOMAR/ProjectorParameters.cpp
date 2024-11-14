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
#include "ProjectorParameters.H"
#include "Debug.H"
#include "Format.H"
#include "ParmParse.H"


//-----------------------------------------------------------------------
// Static variable definitions for ProjectorParameters.
//-----------------------------------------------------------------------
std::unique_ptr<ProjectorParameters> ProjectorParameters::s_defPtr;
bool ProjectorParameters::s_constructorLock = false;

//-----------------------------------------------------------------------
// The default constructor sets default / reads parameters.
//-----------------------------------------------------------------------
ProjectorParameters::ProjectorParameters()
{
    if (s_constructorLock) return;

    if (s_defPtr == NULL) {
        s_constructorLock = true;

        s_defPtr.reset(new ProjectorParameters);
        createDefaults();

        s_constructorLock = false;
    }

    // Copy default values.
    *this = *s_defPtr;
}


//-----------------------------------------------------------------------
// For the bean counters like me.
//-----------------------------------------------------------------------
void
ProjectorParameters::freeMemory()
{
    s_defPtr.reset();
}


//-----------------------------------------------------------------------
// It's nice to be able to see these parameters in pout.*.
//-----------------------------------------------------------------------
void
ProjectorParameters::dump() const
{
    pout() << "ProjectorParameters:\n" << Format::indent() << std::flush;

    pout() << "doLevelProj = " << (doLevelProj ? "true" : "false") << "\n";
    pout() << "doInitProj = " << (doInitProj ? "true" : "false") << "\n";
    pout() << "doSyncProj = " << (doSyncProj ? "true" : "false") << "\n";

    pout() << "absTol = " << absTol << "\n";
    pout() << "relTol = " << relTol << "\n";
    pout() << "numSmoothDown = " << numSmoothDown << "\n";
    pout() << "numSmoothUp = " << numSmoothUp << "\n";
    pout() << "numSmoothBottom = " << numSmoothBottom << "\n";
    pout() << "numSmoothPrecond = " << numSmoothPrecond << "\n";
    pout() << "prolongOrder = " << prolongOrder << "\n";
    pout() << "prolongOrderFMG = " << prolongOrderFMG << "\n";
    pout() << "numSmoothUpFMG = " << numSmoothUpFMG << "\n";
    pout() << "maxDepth = " << maxDepth << "\n";
    pout() << "numCycles = " << numCycles << "\n";
    pout() << "maxIters = " << maxIters << "\n";
    pout() << "hang = " << hang << "\n";
    pout() << "normType = " << normType << "\n";
    pout() << "verbosity = " << verbosity << "\n";
    pout() << "relaxMethod = ";
    switch (relaxMethod) {
        case RelaxMethod::NONE:     pout() << "NONE\n";     break;
        case RelaxMethod::POINT:    pout() << "POINT\n";    break;
        case RelaxMethod::JACOBI:   pout() << "JACOBI\n";   break;
        case RelaxMethod::JACOBIRB: pout() << "JACOBIRB\n"; break;
        case RelaxMethod::GS:       pout() << "GS\n";       break;
        case RelaxMethod::GSRB:     pout() << "GSRB\n";     break;
        case RelaxMethod::VERTLINE: pout() << "VERTLINE\n"; break;
        default:                    pout() << "UNKNOWN\n";  break;
    }

    pout() << "bottom_absTol = " << bottom_absTol << "\n";
    pout() << "bottom_relTol = " << bottom_relTol << "\n";
    pout() << "bottom_small = " << bottom_small << "\n";
    pout() << "bottom_hang = " << bottom_hang << "\n";
    pout() << "bottom_maxIters = " << bottom_maxIters << "\n";
    pout() << "bottom_maxRestarts = " << bottom_maxRestarts << "\n";
    pout() << "bottom_normType = " << bottom_normType << "\n";
    pout() << "bottom_verbosity = " << bottom_verbosity << "\n";
    pout() << "bottom_numSmoothPrecond = " << bottom_numSmoothPrecond << "\n";

    pout() << Format::unindent << std::endl;
}


//-----------------------------------------------------------------------
// Fills the *s_defPtr object.
//-----------------------------------------------------------------------
void
ProjectorParameters::createDefaults()
{
    ParmParse    pp("proj");
    Vector<int>  vint(SpaceDim);
    Vector<Real> vreal(SpaceDim);


    s_defPtr->doLevelProj = true;
    pp.query("doLevelProj", s_defPtr->doLevelProj);

    s_defPtr->doInitProj = true;
    pp.query("doInitProj", s_defPtr->doInitProj);

    s_defPtr->doSyncProj = true;
    pp.query("doSyncProj", s_defPtr->doSyncProj);


    // Multigrid settings...
    s_defPtr->absTol = 1.0e-12;
    pp.query("absTol", s_defPtr->absTol);

    s_defPtr->relTol = 1.0e-6;
    pp.query("relTol", s_defPtr->relTol);

    s_defPtr->numSmoothDown = 16;
    pp.query("numSmoothDown", s_defPtr->numSmoothDown);

    s_defPtr->numSmoothUp = 16;
    pp.query("numSmoothUp", s_defPtr->numSmoothUp);

    s_defPtr->numSmoothBottom = 2;
    pp.query("numSmoothBottom", s_defPtr->numSmoothBottom);

    s_defPtr->numSmoothPrecond = 2;
    pp.query("numSmoothPrecond", s_defPtr->numSmoothPrecond);

    s_defPtr->prolongOrder = 1;
    pp.query("prolongOrder", s_defPtr->prolongOrder);

    s_defPtr->prolongOrderFMG = 3;
    pp.query("prolongOrderFMG", s_defPtr->prolongOrderFMG);

    s_defPtr->numSmoothUpFMG = 2;
    pp.query("numSmoothUpFMG", s_defPtr->numSmoothUpFMG);

    s_defPtr->maxDepth = -1;
    pp.query("maxDepth", s_defPtr->maxDepth);

    s_defPtr->numCycles = -1;
    pp.query("numCycles", s_defPtr->numCycles);

    s_defPtr->maxIters = 10;
    pp.query("maxIters", s_defPtr->maxIters);

    s_defPtr->hang = 1.0e-7;
    pp.query("hang", s_defPtr->hang);

    s_defPtr->normType = 2;
    pp.query("normType", s_defPtr->normType);

    s_defPtr->verbosity = 4;
    pp.query("verbosity", s_defPtr->verbosity);

    s_defPtr->relaxMethod = RelaxMethod::GSRB;
    pp.query("relaxMethod", s_defPtr->relaxMethod);
    CH_verify(0 <= s_defPtr->relaxMethod);
    CH_verify(s_defPtr->relaxMethod < RelaxMethod::_NUM_RELAXMETHODS);

    // Bottom solver settings...
    s_defPtr->bottom_absTol = 1.0e-6;
    pp.query("bottom_absTol", s_defPtr->bottom_absTol);

    s_defPtr->bottom_relTol = 1.0e-4;
    pp.query("bottom_relTol", s_defPtr->bottom_relTol);

    s_defPtr->bottom_small = 1.0e-30;
    pp.query("bottom_small", s_defPtr->bottom_small);

    s_defPtr->bottom_hang = 1.0e-7;
    pp.query("bottom_hang", s_defPtr->bottom_hang);

    s_defPtr->bottom_maxIters = 80;
    pp.query("bottom_maxIters", s_defPtr->bottom_maxIters);

    s_defPtr->bottom_maxRestarts = 5;
    pp.query("bottom_maxRestarts", s_defPtr->bottom_maxRestarts);

    s_defPtr->bottom_normType = 2;
    pp.query("bottom_normType", s_defPtr->bottom_normType);

    s_defPtr->bottom_verbosity = 0;
    pp.query("bottom_verbosity", s_defPtr->bottom_verbosity);

    s_defPtr->bottom_numSmoothPrecond = 2;
    pp.query("bottom_numSmoothPrecond", s_defPtr->bottom_numSmoothPrecond);

    // Send defaults to pout.
    s_defPtr->dump();
}
