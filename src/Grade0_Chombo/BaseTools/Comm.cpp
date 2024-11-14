#include "Comm.H"

#ifdef CH_MPI
#   include "SPMD.H"
#endif

namespace Comm {


// -----------------------------------------------------------------------------
// Reduces a set of integers over many processors.
// -----------------------------------------------------------------------------
void
reduce (int&   a_val,
        MPI_Op a_mpiOp)
{
#ifdef CH_MPI
    int localVal = a_val;
    int result = MPI_Allreduce(&localVal, &a_val, 1, MPI_INT,
                               a_mpiOp, Chombo_MPI::comm);

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in Comm:reduce");
    }
#else
    (void)a_val;
    (void)a_mpiOp;
#endif
}


// -----------------------------------------------------------------------------
// Reduces a set of real numbers over many processors.
// -----------------------------------------------------------------------------
void
reduce (Real&  a_val,
        MPI_Op a_mpiOp)
{
#ifdef CH_MPI
    Real localVal = a_val;
    int result = MPI_Allreduce(&localVal, &a_val, 1, MPI_CH_REAL,
                               a_mpiOp, Chombo_MPI::comm);

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in Comm:reduce");
    }
#else
    (void)a_val;
    (void)a_mpiOp;
#endif
}


}; // end namespace Comm
