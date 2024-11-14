#include "OutputParameters.H"
#include "Format.H"
#include "ParmParse.H"
#include "Debug.H"


//-----------------------------------------------------------------------
// Static variable definitions for OutputParameters.
//-----------------------------------------------------------------------
std::unique_ptr<OutputParameters> OutputParameters::s_defPtr;
bool                              OutputParameters::s_constructorLock = false;


//-----------------------------------------------------------------------
// The default constructor sets default / reads parameters.
//-----------------------------------------------------------------------
OutputParameters::OutputParameters()
{
    if (s_constructorLock) return;

    if (s_defPtr == NULL) {
        s_constructorLock = true;

        s_defPtr.reset(new OutputParameters);
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
OutputParameters::freeMemory()
{
    s_defPtr.reset();
}


//-----------------------------------------------------------------------
// It's nice to be able to see these parameters in pout.*.
//-----------------------------------------------------------------------
void
OutputParameters::dump() const
{
    pout() << "OutputParameters:\n" << Format::indent() << std::flush;

    // pout() << "verbosity = " << verbosity << "\n";
    // pout() << "doFlowchart = " << (doFlowchart ? "true" : "false") << "\n";

    pout() << "plotInterval = " << plotInterval << "\n";
    pout() << "plotPeriod = " << plotPeriod << "\n";
    pout() << "plotPrefix = " << plotPrefix << "\n";
    pout() << "checkpointInterval = " << checkpointInterval << "\n";
    pout() << "checkpointPrefix = " << checkpointPrefix << "\n";

    pout() << Format::unindent << std::endl;
}


//-----------------------------------------------------------------------
// Fills the *s_defPtr object.
//-----------------------------------------------------------------------
void
OutputParameters::createDefaults()
{
    ParmParse pp("output");
    Vector<int> vint(SpaceDim);

    s_defPtr->verbosity = 1;
    pp.query("verbosity", s_defPtr->verbosity);
    CH_verify(s_defPtr->verbosity >= 0);

    s_defPtr->doFlowchart = false;
    pp.query("doFlowchart", s_defPtr->doFlowchart);

    { // plot block
        bool plotScheduled = false;
        s_defPtr->plotInterval = -1;
        s_defPtr->plotPeriod = -1.0;
        s_defPtr->plotPrefix = std::string("plot_");
        if (pp.query("plotInterval", s_defPtr->plotInterval)) {
            plotScheduled = true;
        }
        if (pp.query("plotPeriod", s_defPtr->plotPeriod)) {
            plotScheduled = true;
        }
        if (plotScheduled) {
            pp.query("plotPrefix", s_defPtr->plotPrefix);
        } else {
            MAYDAYERROR("No plots scheduled. You must set either "
                        "output.plotInterval or output.plotPeriod");
        }
    }

    { // checkpoint block
        bool checkpointScheduled = false;
        s_defPtr->checkpointInterval = -1;
        s_defPtr->checkpointPrefix = std::string("chkpt_");

        if (pp.query("checkpointInterval", s_defPtr->checkpointInterval)) {
            checkpointScheduled = true;
        }
        if (checkpointScheduled) {
            pp.query("checkpointPrefix", s_defPtr->checkpointPrefix);
        } else {
            MAYDAYWARNING("No checkpoints scheduled");
        }
    }

    pout() << endl;

    // Send defaults to pout.
    s_defPtr->dump();
}
