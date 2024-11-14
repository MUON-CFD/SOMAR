#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
#include <cstring>
// #extern "C" {      // The #extern "C" might have been here for a reason...
#include <unistd.h>
// }

#include "SPMD.H"
#include "parstream.H"
#include "BaseNamespaceHeader.H"

using std::cout;
using std::endl;

// try a 30 Mbyte max message size and see if that helps.

long long CH_MAX_MPI_MESSAGE_SIZE = 30*1024*1024;
long long CH_MaxMPISendSize = 0;
long long CH_MaxMPIRecvSize  = 0;

int reportMPIStats()
{
  pout()<<"Chombo message size limit:"<< CH_MAX_MPI_MESSAGE_SIZE<<"\n"
        <<"Max send message size:"<<CH_MaxMPISendSize<<"\n"
        <<"Max recv message size:"<<CH_MaxMPIRecvSize<<std::endl;
  return 0;
}

Vector<int> pids;

#ifndef CH_MPI

int procID()
{
  return 0;
}

// reset this to fool serial code into thinking its parallel
int num_procs = 1 ;

int GetRank(int /*pid*/)
{
  return 0;
}

int GetPID(int rank)
{
  CH_assert(rank == 0);
  return getpid();
}

unsigned int numProc()
{
  return num_procs;
}

#else // CH_MPI version

int GetPID(int rank)
{
  if (rank<0) return -1;
  if ((unsigned int)rank>=pids.size()) return -2;
  return pids[rank];
}

int GetRank(int pid)
{
  for (unsigned int i=0; i<pids.size(); ++i)
    {
      if (pids[i]== pid) return (int)i;
    }
  return -1;
}

int procID()
{
  static int ret = -1;
  if (ret==-1)
  {
//     int proc = getpid();
//     gather(pids, proc, 0);
//     broadcast(pids, 0);
    MPI_Comm_rank(Chombo_MPI::comm, &ret);
  }
  return ret;
}

unsigned int numProc()
{
  static int ret = -1;
  if (ret == -1)
  {
    MPI_Comm_size(Chombo_MPI::comm, &ret);
  }
  return ret;
}

// hopefully static copy of opaque handles
MPI_Comm Chombo_MPI::comm = MPI_COMM_WORLD;

#endif // CH_MPI




// template < >
// void gather<IntVectSet>(Vector<IntVectSet>& a_outVec,
//             const IntVectSet& a_input,
//             int a_dest)
// {
//   Vector<Box> boxes = a_input.boxes();
//   Vector<Vector<Box> > all_boxes;

//   gather (all_boxes, boxes, a_dest);

//   const int num_vecs = all_boxes.size();
//   a_outVec.resize(num_vecs);

//   for (int i = 0; i < num_vecs; ++i)
//   {
//     IntVectSet& ivs = a_outVec[i];
//     const Vector<Box>& vb = all_boxes[i];
//     for (int ibox = 0; ibox < vb.size(); ++ibox)
//     {
//       ivs |= vb[ibox];
//     }
//   }
// }

// return id of unique processor for special serial tasks
int
uniqueProc(const SerialTask::task& a_task)
{
#ifdef NDEBUG
    switch (a_task)
    {
    case SerialTask::compute:
    default:
        return (0);
        //break;  // unreachable break can generate compiler warning
    }
#else
    // in mpi, the debugger attaches to process 0
    (void)a_task;
    return (0);
#endif
}

#include "BaseNamespaceFooter.H"
