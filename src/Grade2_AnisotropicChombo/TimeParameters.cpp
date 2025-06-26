#include "TimeParameters.H"
#include "Format.H"
#include "ParmParse.H"
#include "AnisotropicRefinementTools.H"
#include "BasicIO.H"
#include "Debug.H"


//-----------------------------------------------------------------------
// Static variable definitions for TimeParameters.
//-----------------------------------------------------------------------
std::unique_ptr<TimeParameters> TimeParameters::s_defPtr;
bool                            TimeParameters::s_constructorLock = false;


//-----------------------------------------------------------------------
// The default constructor sets default / reads parameters.
//-----------------------------------------------------------------------
TimeParameters::TimeParameters()
{
    if (s_constructorLock) return;

    if (s_defPtr == NULL) {
        s_constructorLock = true;

        s_defPtr.reset(new TimeParameters());
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
TimeParameters::freeMemory()
{
    s_defPtr.reset();
}


//-----------------------------------------------------------------------
// It's nice to be able to see these parameters in pout.*.
//-----------------------------------------------------------------------
void
TimeParameters::dump() const
{
    pout() << "TimeParameters:\n" << Format::indent() << std::flush;

    pout() << "stopTime = " << stopTime << '\n';
    pout() << "maxSteps = " << maxSteps << '\n';
    pout() << "fixedDt = " << fixedDt << '\n';
    pout() << "maxDt = " << maxDt << '\n';
    pout() << "initDtMult = " << initDtMult << '\n';
    pout() << "dtMult = " << dtMult << '\n';
    pout() << "maxDtGrow = " << maxDtGrow << '\n';

    pout() << "useElementaryController = " << (useElementaryController ? "true" : "false") << '\n';
    pout() << "usePIController = " << (usePIController ? "true" : "false") << '\n';
    pout() << "usePIDController = " << (usePIDController ? "true" : "false") << '\n';
    if (useElementaryController || usePIController || usePIDController) {
        pout() << "absTol = " << absTol << '\n';
        pout() << "relTol = " << relTol << '\n';
    }

    if (isRestart) {
        pout() << "restartFile = " << restartFile << '\n';
    } else {
        pout() << "Not restarting from checkpoint.\n";
    }

    pout() << Format::unindent << std::endl;
}


//-----------------------------------------------------------------------
// Fills the *s_defPtr object.
//-----------------------------------------------------------------------
void
TimeParameters::createDefaults()
{
    ParmParse pp("time");
    Vector<int>  vint(SpaceDim);
    Vector<Real> vreal(SpaceDim);

    pp.get("stopTime", s_defPtr->stopTime);
    CH_verify(s_defPtr->stopTime > 0.0);

    pp.get("maxSteps", s_defPtr->maxSteps);
    CH_verify(s_defPtr->maxSteps >= 0);

    s_defPtr->fixedDt = -1.0;
    pp.query("fixedDt", s_defPtr->fixedDt);

    s_defPtr->maxDt = 1.0e8;
    pp.query("maxDt", s_defPtr->maxDt);

    s_defPtr->initDtMult = 0.1;
    pp.query("initDtMult", s_defPtr->initDtMult);

    s_defPtr->dtMult = 0.80;
    pp.query("dtMult", s_defPtr->dtMult);
    if (s_defPtr->dtMult < 0.0) {
        MAYDAYERROR("dtMult must be a positive number. Default = 0.80.");
    }

    s_defPtr->maxDtGrow = 1.5;
    pp.query("maxDtGrow", s_defPtr->maxDtGrow);
    if (s_defPtr->maxDtGrow <= 0.0) s_defPtr->maxDtGrow = 1.0e8;


    s_defPtr->useElementaryController = false;
    pp.query("useElementaryController", s_defPtr->useElementaryController);

    s_defPtr->usePIController = false;
    pp.query("usePIController", s_defPtr->usePIController);

    s_defPtr->usePIDController = false;
    pp.query("usePIDController", s_defPtr->usePIDController);

    // Read controller error tolerances...
    if (s_defPtr->useElementaryController || s_defPtr->usePIController ||
        s_defPtr->usePIDController) {
        if (!pp.contains("absTol") && !pp.contains("relTol")) {
            MAYDAYERROR(
                "If using an error controller, you must specify either "
                "time.absTol or time.relTol.");
        }

        s_defPtr->absTol = 0.0;
        pp.query("absTol", s_defPtr->absTol);

        s_defPtr->relTol = 0.0;
        pp.query("relTol", s_defPtr->relTol);
    }


    // Read checkpoint stuff...
    s_defPtr->isRestart = pp.contains("restartFile");
    if (s_defPtr->isRestart) {
        pp.get("restartFile", s_defPtr->restartFile);
    } else {
        s_defPtr->restartFile = std::string("");
    }

    pout() << endl;

    // Send defaults to pout.
    s_defPtr->dump();
}
