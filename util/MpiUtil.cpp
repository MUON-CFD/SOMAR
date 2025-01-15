/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2014 University of North Carolina at Chapel Hill
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

// OutPut optimal number of processes.

// Standard headers
#ifdef CH_USE_CUDA_PROJECTOR
#include <cusp/multiply.h> // leave this at first line. There are macros in CHOMBO that interfere with it. 
#endif
#include <sys/utsname.h>
#include <iostream>
#ifdef CH_USE_PYTHON
#include <Python.h> 
#endif


#ifdef CH_MPI
#include "mpi.h"
#endif

// Chombo headers
#include "CH_Attach.H"
#include "ParmParse.H"
#include "memusage.H"
#include "parstream.H"
#include "AnisotropicMeshRefine.H"

// Do we need this?
#ifdef CH_USE_MEMORY_TRACKING
#include "Pool.H"
#endif

#ifndef CH_NTIMER
#include "OldTimer.H"
#endif

// Project headers
#ifdef CH_USE_PYTHON
#include "PyGlue.H"
#endif

#include "Debug.H"
#include "Format.H"
#include "IO.H"
#include "ProblemContext.H"
#include "AMRNSLevelFactory.H"
#include "AnisotropicAMR.H"



// // #define TRAP_FPE  //(should be off by default)
// #ifdef TRAP_FPE
// // Previous versions of glibc require the following code:
// #include "parstream.H"
// extern "C" {
// #include <fpu_control.h>
// }

// /* IM: Invalid operation mask
//  * DM: Denormalized operand mask
//  * ZM: Zero-divide mask
//  * OM: Overflow mask
//  * UM: Underflow mask
//  * PM: Precision (inexact result) mask */
// static void __attribute__((constructor)) trapfpe(void)
// {
//     pout() << " Turning on floating-point traps! " << std::endl;
//     // fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_ZM | _FPU_MASK_OM |
//     // _FPU_MASK_UM); fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_ZM);
//     // fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_OM |
//     // _FPU_MASK_UM);
//     fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_UM);
//     // fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_ZM |
//     // _FPU_MASK_OM); fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_IM |
//     // _FPU_MASK_ZM | _FPU_MASK_OM | _FPU_MASK_DM | _FPU_MASK_UM);
//     fpu_control_t
//     // cw = _FPU_DEFAULT;
//     _FPU_SETCW(cw);
//     /* On x86, this expands to: */
//     /* unsigned int cw = 0x037f & ~(0x01 | 0x04 | 0x08); */
//     /* __asm__ ("fldcw %0" : : "m" (*&cw));              */
// }
// #endif


// -----------------------------------------------------------------------------
// Setup and shutdown routines. This is the main()
// function, but the real work is done in mainLoop().
// -----------------------------------------------------------------------------
#ifdef CH_USE_PYTHON
// this starts the Python interpreter. It shuts down when Python goes out 
// of scope. Being a global variable, that means at the very end
Py Python; 
#endif
int
main(int argc, char* argv[])
{
#ifdef CH_MPI
    MPI_Init(&argc, &argv);
#endif

    // Reset the stdout color
    std::cout << Format::none << std::flush;

    // Open the input file
    char*     in_file = argv[1];
    ParmParse pp(argc - 2, argv + 2, NULL, in_file);
    pout() << "Input file: " << in_file << endl;

#ifdef CH_MPI
#ifdef CH_AIX
    H5dont_atexit();
#endif
#endif

    // #ifdef TRAP_FPE
    //     trapfpe();
    // #endif

#ifndef CH_NTIMER
    // Start yer timers.
    OldTimer Everything;
    Everything.start();
#endif

#ifdef CH_MPI
    pout() << "Using MPI." << endl;

    int chomboRank     = procID();
    int chomboNumRanks = numProc();
    int procid         = getpid();
    pout() << "chomboRank = " << chomboRank
           << "\nchomboNumRanks = " << chomboNumRanks << "\nPID = " << procid
           << endl;
#else
    pout() << "Not using MPI." << endl;
#endif

    // Make sure the user knows if we are in debug mode.
    if (debugMode) {
        pout() << "*** DEBUG MODE ***" << endl;
        IO::tout(0) << Format::hired << "*** DEBUG MODE ***" << Format::none
                    << endl;
    } else {
        pout() << "*** RELEASE MODE ***" << endl;
        IO::tout(0) << Format::higreen << "*** RELEASE MODE ***" << Format::none
                    << endl;
    }

    pout() << "SpaceDim = " << SpaceDim << endl;

#ifdef CH_USE_MEMORY_TRACKING
    pout() << "Using memory tracking." << endl;
#endif

    // Setting this to true will only register the debugger in DEBUG mode
    // if you are using <= 4 procs!
    bool useDebugger = true;
    if (useDebugger && debugMode && chomboNumRanks <= 4) {
        IO::tout(0) << "Registering debugger..." << flush;
        registerDebugger();
        IO::tout(0) << "done." << endl;
    }


    // Setup AMR and run the simulation
    pout() << "\n\n" << flush;
    {

        
        // This sets up default parameters and performs the one and only
        // read from the input file.
        const ProblemContext* ctx = ProblemContext::getInstance();

        // NavierStokes factory setup. This object is in charge of creating
        // level solvers when needed.
        AnisotropicAMRLevelFactory* nsFactPtr = new AMRNSLevelFactory();

        // The AMR driver...
       

        AnisotropicAMR amr(nsFactPtr, ctx->amr);

        

        // Free memory
        delete nsFactPtr;
        nsFactPtr = NULL;
    }
    pout() << endl;
    barrier();
    
    // Free statically allocated memory.
    ProblemContext::freeMemory();
    AnisotropicMeshRefine::deleteBuffer();
    
    
#ifndef CH_NTIMER
    // Stop timing.
    Everything.stop();

    Real        end_memory = get_memory_usage_from_OS();
    std::string end_time   = Format::textTime(Everything.wc_time());

    pout() << "Everything done.\n"
           << "mem usage: " << end_memory << "MB\n"
           << "elapsed time: " << end_time << endl;
    IO::tout(0) << "Elapsed time = " << end_time << "." << endl;
#endif

#ifdef CH_USE_MEMORY_TRACKING
    dumpmemoryatexit();
#endif

#ifndef CH_NTIMER
    CH_TIMER_REPORT();
#endif

#ifdef CH_MPI
    MPI_Finalize();
#endif
    
      

    return 0;
}
