#include "AMRParameters.H"
#include "Format.H"
#include "ParmParse.H"
#include "AnisotropicRefinementTools.H"
#include "BasicIO.H"
#include "Debug.H"


//-----------------------------------------------------------------------
// Static variable definitions for AMRParameters.
//-----------------------------------------------------------------------
std::unique_ptr<AMRParameters> AMRParameters::s_defPtr;
bool                           AMRParameters::s_constructorLock = false;


//-----------------------------------------------------------------------
// The default constructor sets default / reads parameters.
//-----------------------------------------------------------------------
AMRParameters::AMRParameters(const BaseParameters& a_baseParams)
{
    if (s_constructorLock) return;

    if (s_defPtr == NULL) {
        s_constructorLock = true;

        s_defPtr.reset(new AMRParameters(a_baseParams));
        createDefaults(a_baseParams);

        s_constructorLock = false;
    }

    // Copy default values.
    *this = *s_defPtr;
}


//-----------------------------------------------------------------------
// For the bean counters like me.
//-----------------------------------------------------------------------
void
AMRParameters::freeMemory()
{
    s_defPtr.reset();
}


//-----------------------------------------------------------------------
// It's nice to be able to see these parameters in pout.*.
//-----------------------------------------------------------------------
void
AMRParameters::dump() const
{
    pout() << "AMRParameters:\n" << Format::indent() << std::flush;


    pout() << "maxLevel = " << maxLevel << "\n";
    pout() << "numLevels = " << numLevels << "\n";
    pout() << "refRatios = " << refRatios << "\n";
    pout() << "maxGridSize = " << maxGridSize << "\n";
    pout() << "bufferSize = " << bufferSize << "\n";
    pout() << "fillRatio = " << fillRatio << "\n";

    pout() << "regridIntervals = " << regridIntervals << "\n";
    pout() << "useSubcycling = " << (useSubcycling ? "true" : "false") << "\n";

    pout() << "tagIB = " << (tagIB ? "true" : "false") << "\n";
    pout() << "velTagTol = " << velTagTol << "\n";
    pout() << "bTagTol = " << bTagTol << "\n";
    pout() << "TTagTol = " << TTagTol << "\n";
    pout() << "STagTol = " << STagTol << "\n";
    pout() << "bpertTagTol = " << bpertTagTol << "\n";
    pout() << "TpertTagTol = " << TpertTagTol << "\n";
    pout() << "SpertTagTol = " << SpertTagTol << "\n";
    pout() << "scalarsTagTol = " << scalarsTagTol << "\n";
    pout() << "growTags = " << growTags << "\n";

    pout() << Format::unindent << std::endl;
}


//-----------------------------------------------------------------------
// Fills the *s_defPtr object.
//-----------------------------------------------------------------------
void
AMRParameters::createDefaults(const BaseParameters& a_baseParams)
{
    ParmParse pp("amr");
    std::vector<int> vint(SpaceDim);

    pp.get("maxLevel", s_defPtr->maxLevel);
    CH_verify(s_defPtr->maxLevel >= 0);
    s_defPtr->numLevels = s_defPtr->maxLevel + 1;

    vint = std::vector<int>(SpaceDim, 4);
    pp.queryarr("refRatio", vint, 0, SpaceDim);
    s_defPtr->refRatios =
        std::vector<IntVect>(s_defPtr->numLevels, IntVect(vint));
    //
    for (int lev = 0; lev < s_defPtr->maxLevel; ++lev) {
        std::ostringstream str;
        str << "refRatio_lev" << lev;
        bool levSet = pp.queryarr(str.str().c_str(), vint, 0, SpaceDim);

        if (levSet) {
            s_defPtr->refRatios[lev] = IntVect(D_DECL(vint[0], vint[1], vint[2]));
        }
    }
    s_defPtr->refRatios[s_defPtr->maxLevel] = IntVect::Unit; // This is needed by Chombo.
    for (unsigned int lev = 0; lev < s_defPtr->refRatios.size(); ++lev)
        for (int dir = 0; dir < SpaceDim; ++dir)
            CH_verify(s_defPtr->refRatios[lev][dir] >= 1);


    // // Set automated values...
    // {
    //     const IntVect        maxBaseGridSize = a_baseParams.maxBaseGridSize;
    //     const ProblemDomain& baseDomain      = a_baseParams.domain;

    //     // maxGridSize
    //     if (s_defPtr->maxLevel > 0) {
            // chooseMaxGridSize(s_defPtr->maxGridSize,
            //                   s_defPtr->refRatios,
            //                   a_baseParams);

    //     } else {
    //         s_defPtr->maxGridSize = maxBaseGridSize;
    //     }
    // }
    s_defPtr->maxGridSize = a_baseParams.maxBaseGridSize; // For now.

    // ...then, give the user a chance to override.
    if (pp.contains("maxGridSize")) {
        vint = std::vector<int>(SpaceDim, 0);
        pp.queryarr("maxGridSize", vint, 0, SpaceDim);
        for (int dir = 0; dir < SpaceDim; ++dir) {
            s_defPtr->maxGridSize[dir] = vint[dir];
            CH_verify(s_defPtr->maxGridSize[dir] >= 0);
        }
    }

    s_defPtr->bufferSize = 1;
    pp.query("bufferSize", s_defPtr->bufferSize);
    CH_verify(s_defPtr->bufferSize > 0);

    s_defPtr->fillRatio = 0.80;
    pp.query("fillRatio", s_defPtr->fillRatio);

    const int numReadLevels = Max(s_defPtr->maxLevel, 1);
    std::vector<int> vintLevels(numReadLevels, 10);
    s_defPtr->regridIntervals = vintLevels;
    if (s_defPtr->maxLevel > 0) {
        pp.queryarr(
            "regridIntervals", s_defPtr->regridIntervals, 0, numReadLevels);
    }

    s_defPtr->useSubcycling = true;
    pp.query("useSubcycling", s_defPtr->useSubcycling);

    // Tag tol
    s_defPtr->tagIB         = false;
    s_defPtr->velTagTol     = -1.0;
    s_defPtr->bTagTol       = -1.0;
    s_defPtr->TTagTol       = -1.0;
    s_defPtr->STagTol       = -1.0;
    s_defPtr->bpertTagTol   = -1.0;
    s_defPtr->TpertTagTol   = -1.0;
    s_defPtr->SpertTagTol   = -1.0;
    s_defPtr->scalarsTagTol = std::vector<Real>(0);
    s_defPtr->growTags      = 0;

    if (s_defPtr->maxLevel > 0) {
        bool tolSpecified = false;
        tolSpecified |= pp.query("tagIB", s_defPtr->tagIB);
        tolSpecified |= pp.query("velTagTol", s_defPtr->velTagTol);
        tolSpecified |= pp.query("bTagTol", s_defPtr->bTagTol);
        tolSpecified |= pp.query("TTagTol", s_defPtr->TTagTol);
        tolSpecified |= pp.query("STagTol", s_defPtr->STagTol);
        tolSpecified |= pp.query("bpertTagTol", s_defPtr->bpertTagTol);
        tolSpecified |= pp.query("TpertTagTol", s_defPtr->TpertTagTol);
        tolSpecified |= pp.query("SpertTagTol", s_defPtr->SpertTagTol);

        int numVals = pp.countval("scalarsTagTol");
        if (numVals > 0) {
            s_defPtr->scalarsTagTol = std::vector<Real>(numVals, -1.0);
            pp.getarr("scalarsTagTol", s_defPtr->scalarsTagTol, 0, numVals);
            tolSpecified = true;
        }

        if (!tolSpecified) {
            MAYDAYWARNING(
                "No dynamic tagging criteria specified. Static refinement "
                "assumed.");
        }

        pp.query("growTags", s_defPtr->growTags);
    }

    pout() << endl;

    // Send defaults to pout.
    s_defPtr->dump();
}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void
AMRParameters::chooseMaxGridSize(IntVect&                    a_maxGridSize,
                                 const std::vector<IntVect>& a_refRatio,
                                 const BaseParameters&       a_baseParams)
{
    const int localNumLevels = a_refRatio.size();
    std::vector<IntVect> vMaxGridSizes(localNumLevels, a_baseParams.maxBaseGridSize);

    IntVect ref = IntVect::Unit;
    for(int lev = 1; lev < localNumLevels; ++lev) {
        ref *= a_refRatio[lev - 1];
        chooseMaxGridSize(vMaxGridSizes[lev], ref, a_baseParams);
    }

    pout() << "Possible maxGridSizes = " << vMaxGridSizes << endl;
    BasicIO::tout(0) << "Possible maxGridSizes = " << vMaxGridSizes << endl;

    // Use the top-level grid size since we need to help it be more efficient
    // than any other level.
    a_maxGridSize = vMaxGridSizes.back();
}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void
AMRParameters::chooseMaxGridSize(IntVect&              a_maxGridSize,
                                 const IntVect&        a_refRatio,
                                 const BaseParameters& a_baseParams)
{
    ProblemDomain fineDomain;
    refine(fineDomain, a_baseParams.domain, a_refRatio);

    ProblemDomain isoFineDomain;
    IntVect refToIsotropy;
    coarsenToIsotropy(
        isoFineDomain, refToIsotropy, a_baseParams.domain, a_baseParams.L);

    a_maxGridSize = refToIsotropy;

    // Cut everything until the smallest side is not > 8.
    {
        D_TERM(
        int minSize = a_maxGridSize[0];,
        minSize = min(minSize, a_maxGridSize[1]);,
        minSize = min(minSize, a_maxGridSize[2]);)

       while (minSize > 8) {
            // int l = std::accumulate(&a_maxGridSize[0], &a_maxGridSize[0] + SpaceDim, 1, lcm);
            // pout() << "lcm of " << a_maxGridSize << " is " << l << endl;

            // Find least common factor.
            int f = 2;
            while (f < minSize) {
                if (D_TERM(   a_maxGridSize[0] % f == 0,
                           && a_maxGridSize[1] % f == 0,
                           && a_maxGridSize[2] % f == 0   )) {
                    a_maxGridSize /= f;
                    minSize /= f;
                    break;
                }
                ++f; // Not very efficient, but not a problem here.
            }

            // Were we able to cut the box?
            if (f == minSize) {
                break; // Nope. Give up this tactic.
            }
        }
    }

    // Double everything until the smallest side is not < 8.
    {
        D_TERM(
        int minSize = a_maxGridSize[0];,
        minSize = min(minSize, a_maxGridSize[1]);,
        minSize = min(minSize, a_maxGridSize[2]);)

        int p = 1;
        while (minSize < 8) {
            p *= 2;
            minSize *= 2;
        }
        a_maxGridSize *= p;
    }

    // Double everything until the largest side is not < 32.
    {
        D_TERM(
        int maxSize = a_maxGridSize[0];,
        maxSize = max(maxSize, a_maxGridSize[1]);,
        maxSize = max(maxSize, a_maxGridSize[2]);)

        int p = 1;
        while (maxSize < 32) {
            p *= 2;
            maxSize *= 2;
        }
        a_maxGridSize *= p;
    }
}
