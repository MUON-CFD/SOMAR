#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//
// These are the inevitable functions that people
// can't avoid using when writing a SPMD program.
// It is a minimal set from which more elaborate
// functionality can be generated.  As always, a
// user is free to utilize the entire MPI programming
// on their own platform.  The functions are
// assured to work on all machines supported.
//

#ifndef _SPMD_H_
#define _SPMD_H_

#include <array>
#include <string>

#include "REAL.H"
#include "Vector.H"
#include "MayDay.H"
#include "Misc.H"
#include <type_traits>
#ifdef CH_MPI
#include <mpi.h>
#endif

#include "BaseNamespaceHeader.H"

namespace std{
  template<class...>
  using void_t=void;
}

// =============================================================================
// Begin linearization functions

// -----------------------------------------------------------------------------
template <typename T, typename = void>
struct HasMemberLinearSize : std::false_type
{};

template <typename T>
struct HasMemberLinearSize<
    T,
    std::void_t<decltype(std::declval<T>().linearSize())>> : std::true_type
{};


template <typename T>
int
linearSize(const T& a_input);


template <typename T>
int
linearSize(const std::vector<T>& a_input);


// -----------------------------------------------------------------------------
template <typename T, typename = void>
struct HasMemberLinearOut : std::false_type
{};

template <typename T>
struct HasMemberLinearOut<
    T,
    std::void_t<decltype(std::declval<T>().linearOut(nullptr))>> : std::true_type
{};


template <typename T>
void
linearOut(void* const a_outBuf, const T& a_input);

template <typename T>
void
linearOut(void* const a_outBuf, const std::vector<T>& a_input);

// -----------------------------------------------------------------------------
template <typename T, typename = void>
struct HasMemberLinearIn : std::false_type
{};

template <typename T>
struct HasMemberLinearIn<
    T,
    std::void_t<decltype(std::declval<T>().linearIn(nullptr))>> : std::true_type
{};


template <typename T>
void
linearIn(T& a_output, const void* const a_inBuf);

template <typename T>
void
linearIn(std::vector<T>& a_output, const void* const a_inBuf);

// End linearization functions.
// =============================================================================


#ifdef CH_MPI
struct Chombo_MPI
{
  static MPI_Comm comm;
};

#else
// this can be changed for debugging parallel code in serial
extern int num_procs ;
#endif

extern long long CH_MAX_MPI_MESSAGE_SIZE;

extern long long CH_MaxMPISendSize;
extern long long CH_MaxMPIRecvSize;


int reportMPIStats();

/// local process ID
/**
    Returns the ID of the locally running process in the range 0 <=
    procID() < numProc().  This has no relation to the operating system
    pid.  There is always a procID() == 0.  */
int procID();

inline int CHprocID()
{
  return procID();
}

/// number of parallel processes
/**
    Returns the number of parallel processes running.
    Always returns at least 1. */
unsigned int numProc();

///all ranks wait here to sync-up
/**
   All MPI ranks wait here to sync-up.  Calls MPI_Barrier(comm)  */
void barrier(void);


#ifdef CH_MPI
///gather from multiple processes
/**
  Gathers a_input of every process onto Vector<T> a_outVec on process
  number a_dest.  It is an error if a_dest < 0 or a_dest >= numProc().
  a_outVec is numProc() long with the value of a_input on every
  process in its elements.\\

T must have have the following functions:
  \\

  int linearSize<T>(const T& inputT); \\
  void linearIn<T>(T& a_outputT, const void* const inBuf); \\
  void linearOut<T>(void* a_outBuf, const T& inputT); \\

*/
template <class T>
void gather(Vector<T>& a_outVec, const T& a_input, int a_dest);

// this has to be here so that linearSize<IntVectSet>, linearIn<IntVectSet>,
//  and linearOut<IntVectSet> are defined.
// template < >
// void gather(Vector<IntVectSet>& a_outVec,
//             const IntVectSet& a_input,
//             int a_dest);
#endif

///broadcast to every process
/**
  Broadcasts a_inAndOut to every process from process number a_src.
  It is an error if a_src < 0 or a_src >= numProc().  \\

  T must have have the following functions: \\

  int linearSize<T>(const T& inputT); \\
  void linearIn<T>(T& a_outputT, const void* const inBuf); \\
  void linearOut<T>(void* a_outBuf, const T& inputT); \\ */
template <class T>
void broadcast(T& a_inAndOut,  int a_src);


class SerialTask
{
public:
    enum task
    {
      compute=0
    };
};

int GetPID(int rank);
int GetRank(int pid);

// return id of unique processor for special serial tasks
int
uniqueProc(const SerialTask::task& a_task);

#include "BaseNamespaceFooter.H"

#include "SPMDI.H"

#endif
