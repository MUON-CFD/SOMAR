#include "BaseParameters.H"
#include "Format.H"
#include "ParmParse.H"
#include "AnisotropicRefinementTools.H"
#include "BasicIO.H"
#include "Debug.H"


//-----------------------------------------------------------------------
// Static variable definitions for BaseParameters.
//-----------------------------------------------------------------------
std::unique_ptr<BaseParameters> BaseParameters::s_defPtr;
bool                            BaseParameters::s_constructorLock = false;

//-----------------------------------------------------------------------
// The default constructor sets default / reads parameters.
//-----------------------------------------------------------------------
BaseParameters::BaseParameters()
{
    if (s_constructorLock) return;

    if (s_defPtr == NULL) {
        s_constructorLock = true;

        s_defPtr.reset(new BaseParameters);
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
BaseParameters::freeMemory()
{
    s_defPtr.reset();
}


//-----------------------------------------------------------------------
// It's nice to be able to see these parameters in pout.*.
//-----------------------------------------------------------------------
void
BaseParameters::dump() const
{
    pout() << "BaseParameters:\n" << Format::indent() << std::flush;

    pout() << "L = " << L << '\n';
    pout() << "nx = " << nx << '\n';
    pout() << "dx = " << L / RealVect(nx) << '\n';
    pout() << "nxOffset = " << nxOffset << '\n';
    {  // isPeriodic block
        pout() << "isPeriodic = (" << (isPeriodic[0] ? "true" : "false");
        for (int dir = 1; dir < SpaceDim; ++dir)
            pout() << ", " << (isPeriodic[dir] ? "true" : "false");
        pout() << ")\n";
    }
    pout() << "base domain = " << domain << '\n';

    pout() << "splitDirs = " << splitDirs << '\n';
    pout() << "maxBaseGridSize = " << maxBaseGridSize << '\n';
    pout() << "blockFactor = " << blockFactor << '\n';

    pout() << Format::unindent << std::endl;
}


//-----------------------------------------------------------------------
// Fills the *s_defPtr object.
//-----------------------------------------------------------------------
void
BaseParameters::createDefaults()
{
    ParmParse pp("base");
    Vector<int> vint(SpaceDim);
    Vector<Real> vreal(SpaceDim);

    pp.getarr("L", vreal, 0, SpaceDim);
    s_defPtr->L = RealVect(vreal);

    pp.getarr("nx", vint, 0, SpaceDim);
    s_defPtr->nx = IntVect(vint);

    vint = Vector<int>(SpaceDim, 0);
    pp.queryarr("nxOffset", vint, 0, SpaceDim);
    s_defPtr->nxOffset = IntVect(vint);

    vint = Vector<int>(SpaceDim, 0);
    pp.queryarr("isPeriodic", vint, 0, SpaceDim);
    for (int dir = 0; dir < SpaceDim; ++dir)
        s_defPtr->isPeriodic[dir] = (vint[dir] != 0);

    vint = Vector<int>(SpaceDim, 1);
    pp.queryarr("splitDirs", vint, 0, SpaceDim);
    for (int dir = 0; dir < SpaceDim; ++dir)
        s_defPtr->splitDirs[dir] = (vint[dir] == 0)? 0: 1;

    s_defPtr->domain.define(s_defPtr->nxOffset,
                            s_defPtr->nxOffset + s_defPtr->nx - IntVect::Unit,
                            s_defPtr->isPeriodic);

    // Set automated values...
    {
        // maxBaseGridSize
        chooseMaxBaseGridSize(s_defPtr->maxBaseGridSize,
                              s_defPtr->domain,
                              s_defPtr->L);

        // blockFactor
        s_defPtr->blockFactor = 32;
        for(int dir = 0; dir < SpaceDim; ++dir) {
            s_defPtr->blockFactor =
                min(s_defPtr->blockFactor, s_defPtr->maxBaseGridSize[dir]);
        }
    }

    // ...then, give the user a chance to override.
    {
        if (pp.contains("maxBaseGridSize")) {
            vint = Vector<int>(SpaceDim, 0);
            pp.queryarr("maxBaseGridSize", vint, 0, SpaceDim);
            for (int dir = 0; dir < SpaceDim; ++dir) {
                s_defPtr->maxBaseGridSize[dir] = vint[dir];
                CH_verify(s_defPtr->maxBaseGridSize[dir] >= 0);
            }
        }

        pp.query("blockFactor", s_defPtr->blockFactor);
        CH_assert(s_defPtr->blockFactor >= 0);
    }

    pout() << endl;

    // Send defaults to pout.
    s_defPtr->dump();
}


//------------------------------------------------------------------------------
// Attempts to find the best values.
//------------------------------------------------------------------------------
void
BaseParameters::chooseMaxBaseGridSize(IntVect&             a_maxBaseGridSize,
                                      const ProblemDomain& a_domain,
                                      const RealVect&      a_L)
{
    int nproc = numProc();

    // If we are only using 1 proc, then we consider it a special case.
    if (nproc == 1) {
       a_maxBaseGridSize = a_domain.size();
        return;
    }

    // Coarsen the domain until the dXi is isotropic.
    ProblemDomain isoDomain;
    IntVect refToIsotropy;
    coarsenToIsotropy(isoDomain, refToIsotropy, a_domain, a_L);

    int idealNumProc;
    powOf2Strategy(a_maxBaseGridSize, idealNumProc, isoDomain, nproc);
    a_maxBaseGridSize *= refToIsotropy;
    // IntVect splits = a_domain.size() / a_maxBaseGridSize;

    // pout() << "Setting maxBaseGridSize = " << a_maxBaseGridSize << endl;
    // pout() << "This splits the domain into " << splits << " parts." << endl;
    // pout() << "The ideal number of MPI ranks to use is " << idealNumProc
    //             << "." << endl;

    // BasicIO::tout(0) << "Setting maxBaseGridSize = " << a_maxBaseGridSize << endl;
    // BasicIO::tout(0) << "This splits the domain into " << splits << " parts." << endl;
    // BasicIO::tout(0) << "The ideal number of MPI ranks to use is " << idealNumProc
    //             << "." << endl;

    // // We are done with maxBaseGridSize.
    // a_maxGridSize = 16 * a_refToIsotropy;


    // D_TERM(
    // a_maxGridSize[0] = min(32, a_maxBaseGridSize;
    // )

    // a_blockFactor = 32;
    // for(int dir = 0; dir < SpaceDim; ++dir) {
    //     a_blockFactor = min(a_blockFactor, a_maxBaseGridSize[dir]);
    // }
    // pout() << "Setting blockFactor = " << a_blockFactor << endl;
}


//------------------------------------------------------------------------------
// Attempts to find the best values.
//------------------------------------------------------------------------------
void
BaseParameters::powOf2Strategy(IntVect&             a_maxGridSize,
                               int&                 a_idealNumProc,
                               const ProblemDomain& a_domain,
                               const int            a_nproc)
{
    // Find common power of 2.
    int b = 2;
    int p = 12345;
    for (int dir = 0; dir < SpaceDim; ++dir) {
        int n = a_domain.size(dir);
        int thisp = 0;
        while (n % b == 0) {
            ++thisp;
            n /= b;
        }
        p = min(p, thisp);
    }

    int l = pow((double)b, p);
    IntVect splits = a_domain.size() / l;

    while (splits.product() < a_nproc) {
        splits *= 2;
    }

    const int vol2 = D_TERM(2,*2,*2);
    while (splits.product() > a_nproc * vol2 && splits > IntVect::Unit) {
        splits /= 2;
    }

    a_maxGridSize = a_domain.size() / splits;
    a_idealNumProc = splits.product();

    if (splits < IntVect::Unit) {
        MayDay::Error("Domain deconposition could not be found.");
    }

    if (splits.product() % a_nproc != 0) {
        ostringstream msg;
        msg << "Domain does not split evenly over MPI ranks. Try using "
            << splits.product() << " MPI ranks.";
        MayDay::Warning(msg.str().c_str());
    }
}
