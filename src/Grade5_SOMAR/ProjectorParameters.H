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
#ifndef ___ProjectorParameters_H__INCLUDED___
#define ___ProjectorParameters_H__INCLUDED___

#include <memory>
#include "RealVect.H"


class ProjectorParameters
{
public:
    // Default constructor
    ProjectorParameters();

    bool doLevelProj;
    bool doInitProj;
    bool doSyncProj;

    Real absTol;            // Solver tolerance
    Real relTol;            // Solver relative tolerance
    int  numSmoothDown;     // MG pre smoothing iters
    int  numSmoothUp;       // MG post smoothing iters
    int  numSmoothBottom;   // MG bottom smoothing iters
    int  numSmoothPrecond;  // MG preconditioner smoothing iters
    int  prolongOrder;      // Interpolator used in V-cycles. (0 = const, 1 = linear, ...)
    int  prolongOrderFMG;   // Interpolator used by FMG _after_ V-cycles.
    int  numSmoothUpFMG;    // FMG post-interpolation smoothing iters.
    int  maxDepth;          // Max MG depth (-1 for as deep as possible)
    int  numCycles;         // 1=V-cycle, 2=W-cycle, etc...
    int  maxIters;          // Max number of V-cycles.
    Real hang;              // Stop when residual drops less than this fraction.
    int  normType;          // 0 = inf-norm, etc...
    int  verbosity;         //

    struct RelaxMethod {
        enum {
            NONE     = 0,
            POINT    = 1,
            JACOBI   = 2,
            JACOBIRB = 3,
            GS       = 4,
            GSRB     = 5,
            VERTLINE = 6,
            _NUM_RELAXMETHODS
        };
    };
    int relaxMethod;

    Real bottom_absTol;       // Solver tolerance
    Real bottom_relTol;       // Solver relative tolerance
    Real bottom_small;        //
    Real bottom_hang;         //
    int  bottom_maxIters;     // Max BiCGStab iterations
    int  bottom_maxRestarts;  // Max BiCGStab restarts
    int  bottom_normType;     //
    int  bottom_verbosity;    //
    int  bottom_numSmoothPrecond; // Smoothing iters in preconditioner.

    // You shouldn't need to call this. AnisotropicAMR will do it for you.
    static void
    freeMemory();

    // It's nice to be able to see these parameters in pout.*.
    void
    dump() const;

private:
    static void
    createDefaults();

    static std::unique_ptr<ProjectorParameters> s_defPtr;
    static bool                                 s_constructorLock;
};


#endif //!___ProjectorParameters_H__INCLUDED___
