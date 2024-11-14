#include "BasicIO.H"
#include "SPMD.H"
#include <iomanip>


BasicIO::basic_onullstream<char> BasicIO::s_null_ostream;


// -----------------------------------------------------------------------------
// Print to terminal from one or all procs
// -----------------------------------------------------------------------------
std::ostream&
BasicIO::tout (int a_procid)
{
#ifdef CH_MPI
    static unsigned int width = ceil(log10(numProc()));

    // Print...
    if(a_procid == -1) {
        // ...from all procs with a nice rank label.
        unsigned int saveWidth = std::cout.width();
        return std::cout << "Proc "
                         << std::setw(width)
                         << procID()
                         << std::setw(saveWidth) << ": ";

    } else if(a_procid <= -2) {
        // ... from all procs with no label.
        return std::cout;

    } else if(procID() != a_procid) {
        // ... from a single proc, and this ain't it.
        return s_null_ostream;
    }
#else
    (void)a_procid;
#endif

    // ...from a single proc, and this is it.
    return std::cout;
}
