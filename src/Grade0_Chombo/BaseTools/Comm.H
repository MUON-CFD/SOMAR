/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2017
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
#ifndef ___Comm_H__INCLUDED___
#define ___Comm_H__INCLUDED___

#include <vector>
#include <ostream>
#include "REAL.H"
#include "Box.H"

#ifdef CH_MPI
#   include "mpi.h"
#   include "SPMD.H"
#else
    typedef int MPI_Comm;
    typedef int MPI_Op;
    typedef int MPI_Request;
    typedef int MPI_Status;

    static constexpr int MPI_INT     = 0;
    static constexpr int MPI_DOUBLE  = 0;
    static constexpr int MPI_FLOAT   = 0;
    static constexpr int MPI_BYTE    = 0;

    static constexpr MPI_Op    MPI_MAX     = 0;
    static constexpr MPI_Op    MPI_MIN     = 0;
    static constexpr MPI_Op    MPI_SUM     = 0;
    static constexpr MPI_Op    MPI_PROD    = 0;
    static constexpr MPI_Op    MPI_LAND    = 0;
    static constexpr MPI_Op    MPI_BAND    = 0;
    static constexpr MPI_Op    MPI_LOR     = 0;
    static constexpr MPI_Op    MPI_BOR     = 0;
    static constexpr MPI_Op    MPI_LXOR    = 0;
    static constexpr MPI_Op    MPI_BXOR    = 0;
    static constexpr MPI_Op    MPI_MINLOC  = 0;
    static constexpr MPI_Op    MPI_MAXLOC  = 0;
    static constexpr MPI_Op    MPI_REPLACE = 0;
    static constexpr MPI_Op    MPI_NO_OP   = 0;
#endif

namespace Comm {


#ifdef CH_MPI
    static constexpr MPI_Comm& primaryComm = Chombo_MPI::comm;
    static constexpr int       primaryRank = 0;
#else
    static constexpr MPI_Comm  primaryComm = 0;
    static constexpr int       primaryRank = 0;
#endif


// Forward declarations ========================================================

/// What is the main proc?
inline int
rootRank();

/// Is this the main proc?
inline bool
iAmRoot();

/// Reduces a set of integers over many processors.
void
reduce(int& a_val, MPI_Op a_mpiOp);

/// Reduces a set of real numbers over many processors.
void
reduce(Real& a_val, MPI_Op a_mpiOp);

/// All procs -> destProc
template <class T>
inline void
gather(std::vector<T>& a_out,
       const T&        a_in,
       const int       a_destRank = primaryRank,
       const MPI_Comm  a_comm     = primaryComm);

/// All procs -> rootProc -> all procs
template <class T>
inline void
allgather(std::vector<T>& a_out,
          const T&        a_in,
          const int       a_rootProc = primaryRank,
          const MPI_Comm  a_comm     = primaryComm);

/// All procs -> rootProc -> all procs
/// in-place version. Vector will be grown to fit all elements.
template <class T>
void
allgather2(std::vector<T>& a_inAndOut,
           const int       a_rootProc = primaryRank,
           const MPI_Comm  a_comm     = primaryComm);

/// srcProc -> all procs
template <class T>
inline void
broadcast(T&             a_inAndOut,
          const int      a_srcProc = primaryRank,
          const MPI_Comm a_comm    = primaryComm);

/// wait for all
inline void
barrier(MPI_Comm a_comm = primaryComm);


// Inline definitions ==========================================================

// -----------------------------------------------------------------------------
inline int
rootRank()
{
    return Comm::primaryRank;
}


// -----------------------------------------------------------------------------
inline bool
iAmRoot()
{
    return procID() == Comm::rootRank();
}


// -----------------------------------------------------------------------------
template <class T>
inline void
gather(std::vector<T>& a_out,
       const T&        a_in,
       const int       a_destRank,
       const MPI_Comm  a_comm)
{
    // ::gather(a_out, a_in, a_destProc);

#ifdef CH_MPI
    if (a_comm == MPI_COMM_NULL) return;

    const int commSize = [&a_comm]() {
        int s = -1;
        MPI_Comm_size(a_comm, &s);
        return s;
    }();

    const int myRank = [&a_comm]() {
        int r = -1;
        MPI_Comm_rank(a_comm, &r);
        return r;
    }();

    CH_assert(a_destRank >= 0);
    CH_assert(a_destRank < commSize);
    // now THIS size lives on THIS processor
    int isize = linearSize(a_in);

    // make stuff for linearout
    void* localBuf = mallocMT(isize);
    if (localBuf == NULL) MayDay::Error("out of memory in gather 1");

    // put linearized T into its proper buffer
    linearOut(localBuf, a_in);

    int sendCount = 1;
    int recdCount = 1;

    // need to gather isizes onto processor a_destRank
    int*  vectSize = NULL;
    int*  vectDisp = NULL;
    void* sendBuf  = static_cast<void*>(&isize);
    // allocate received buffer
    if (myRank == a_destRank) {
        vectSize = new int[commSize];
        vectDisp = new int[commSize];
    }

    const int result1 = MPI_Gather(sendBuf,
                                   sendCount,
                                   MPI_INT,
                                   vectSize,
                                   recdCount,
                                   MPI_INT,
                                   a_destRank,
                                   a_comm);

    if (result1 != MPI_SUCCESS)
        MayDay::Error("Comm::gather<T> failed in MPI_Gather 1");

    // make memory for gather, linearin
    void* recdBuf = NULL;
    if (myRank == a_destRank) {
        int itotsize = 0;
        for (int iproc = 0; iproc < commSize; iproc++) {
            vectDisp[iproc] = itotsize;
            itotsize += vectSize[iproc];
        }
        recdBuf = mallocMT(itotsize);
        if (recdBuf == NULL)
            MayDay::Error("Comm::gather<T> is out of memory in gather 2");
    }

    // gather data
    const int result2 = MPI_Gatherv(localBuf,
                                    isize,
                                    MPI_BYTE,
                                    recdBuf,
                                    vectSize,
                                    vectDisp,
                                    MPI_BYTE,
                                    a_destRank,
                                    a_comm);
    if (result2 != MPI_SUCCESS)
        MayDay::Error("Comm::gather<T> failed in MPI_Gather 2");

    if (myRank == a_destRank) {
        // calculate offset into array for current processor
        int ioffset = 0;
        a_out.resize(commSize);
        // need to cast to char* to do pointer arithmetic
        char* arithPtr = (char*)recdBuf;
        for (int iproc = 0; iproc < commSize; iproc++) {
            ioffset           = vectDisp[iproc];
            char* thisProcBuf = arithPtr + ioffset;
            linearIn(a_out[iproc], thisProcBuf);
        }

        // delete memory for dest-specific arrays
        delete[] vectSize;
        delete[] vectDisp;
        freeMT(recdBuf);
    }

    // delete memory for local buffer
    freeMT(localBuf);

#else
    a_out.resize(1);
    a_out[0] = a_in;
#endif
}


// -----------------------------------------------------------------------------
template <class T>
inline void
allgather(std::vector<T>& a_out,
          const T&        a_in,
          const int       a_rootProc,
          const MPI_Comm  a_comm)
{
    // TODO: Use MPI_Allgather
    Comm::gather(a_out, a_in, a_rootProc, a_comm);
    Comm::broadcast(a_out, a_rootProc, a_comm);
}


// -----------------------------------------------------------------------------
template <class T>
void
allgather2(std::vector<T>& a_inAndOut,
           const int       a_rootProc,
           const MPI_Comm  a_comm)
{
    std::vector<std::vector<T>> allElems;
    Comm::gather(allElems, a_inAndOut, a_rootProc, a_comm);
    Comm::broadcast(allElems, a_rootProc, a_comm);

    a_inAndOut.clear();
    for (auto& v : allElems) {
        a_inAndOut.insert(a_inAndOut.end(), v.begin(), v.end());
    }
}


// -----------------------------------------------------------------------------
template <class T>
inline void
broadcast(T&             a_inAndOut,
          const int      a_srcRank,
          const MPI_Comm a_comm)
{
    // ::broadcast(a_inAndOut, a_srcProc);

#ifdef CH_MPI
    if (a_comm == MPI_COMM_NULL) return;

#ifndef NDEBUG
    const int commSize = [&a_comm]() {
        int s = -1;
        MPI_Comm_size(a_comm, &s);
        return s;
    }();
#endif

    const int myRank = [&a_comm]() {
        int r = -1;
        MPI_Comm_rank(a_comm, &r);
        return r;
    }();

    CH_assert(a_srcRank >= 0);
    CH_assert(a_srcRank < commSize);
    int isize = -1;
    if (myRank == a_srcRank) {
        isize = linearSize(a_inAndOut);
    }

    MPI_Bcast(&isize, 1, MPI_INT, a_srcRank, a_comm);
    CH_assert(isize >= 0);

    void* broadBuf = mallocMT(isize);

    if (broadBuf == NULL) {
        MayDay::Error("out of memory in Comm::broadcast");
    }

    // take inAndOut from src and put it into broadBuf
    if (myRank == a_srcRank) {
        linearOut(broadBuf, a_inAndOut);
    }

    // broadcast broadBuf to all procs
    MPI_Bcast(broadBuf, isize, MPI_BYTE, a_srcRank, a_comm);

    if (myRank == a_srcRank) {
        CH_MaxMPISendSize = Max<long long>(CH_MaxMPISendSize, isize);
    } else {
        CH_MaxMPIRecvSize = Max<long long>(CH_MaxMPIRecvSize, isize);
    }
    // take broadBuf and put back into inAndOut if not src
    if (myRank != a_srcRank) {
        linearIn(a_inAndOut, broadBuf);
    }

    // delete memory for buffer
    freeMT(broadBuf);
#endif
}


// -----------------------------------------------------------------------------
inline void
barrier(MPI_Comm a_comm)
{
#ifdef CH_MPI
    MPI_Barrier(a_comm);
#else
    (void)a_comm;
#endif
}


}; // end namespace Comm

#endif //!___Comm_H__INCLUDED___
