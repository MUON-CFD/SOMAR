#include <iostream>
using std::cerr;
using std::cin;
using std::cout;
using std::endl;

#include "Box.H"
#include "LayoutIterator.H"
#include "Vector.H"
#include "parstream.H"

#include "AnisotropicAMRLevel.H"


int AnisotropicAMRLevel::s_verbosity = 0;

//-----------------------------------------------------------------------
bool
AnisotropicAMRLevel::isDefined() const
{
    return m_isDefined;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
AnisotropicAMRLevel::~AnisotropicAMRLevel()
{
    // Make sure we delete levels from the top down.
    if (m_finer_level_ptr) {
        MayDay::Error("You must delete levels from top to bottom!");
    }

    // Make sure we don't try to access this level ever again!
    if (m_coarser_level_ptr) {
        m_coarser_level_ptr->finerLevelPtr(nullptr);
        m_coarser_level_ptr = nullptr;
    }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
AnisotropicAMRLevel::AnisotropicAMRLevel()
{
    m_coarser_level_ptr     = nullptr;
    m_finer_level_ptr       = nullptr;
    m_isDefined             = false;
    m_level                 = 0;
    m_fineRefRatio          = IntVect::Unit;
    m_time                  = 0;
    m_dt                    = 0;
    m_initial_dt_multiplier = 0.1;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
AnisotropicAMRLevel::define(AnisotropicAMRLevel* a_coarser_level_ptr,
                            const Box&           a_problem_domain,
                            int                  a_level,
                            const IntVect&       a_ref_ratio)
{
    ProblemDomain physDomain(a_problem_domain);
    define(a_coarser_level_ptr, physDomain, a_level, a_ref_ratio);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
AnisotropicAMRLevel::define(AnisotropicAMRLevel* a_coarser_level_ptr,
                            const ProblemDomain& a_problem_domain,
                            int                  a_level,
                            const IntVect&       a_ref_ratio)
{
    if (s_verbosity >= 3) {
        pout() << "AnisotropicAMRLevel::define" << endl;
    }

    m_coarser_level_ptr = a_coarser_level_ptr;
    m_problem_domain    = a_problem_domain;
    m_level             = a_level;
    m_fineRefRatio      = a_ref_ratio;
    m_finer_level_ptr   = nullptr;
    m_isDefined         = true;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
AnisotropicAMRLevel::finerLevelPtr(AnisotropicAMRLevel* a_finer_level_ptr)
{
    m_finer_level_ptr = a_finer_level_ptr;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
AnisotropicAMRLevel::dt(Real a_dt)
{
    m_dt = a_dt;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real
AnisotropicAMRLevel::dt() const
{
    return (m_dt);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
const ProblemDomain&
AnisotropicAMRLevel::problemDomain() const
{
    return (m_problem_domain);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Vector<Box>
AnisotropicAMRLevel::boxes() const
{
    return (m_level_grids);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
bool
AnisotropicAMRLevel::hasCoarserLevel() const
{
    return ((m_coarser_level_ptr != nullptr) &&
            (m_coarser_level_ptr->m_isDefined) &&
            (m_coarser_level_ptr->m_level_grids.size() > 0));
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
bool
AnisotropicAMRLevel::hasFinerLevel() const
{
    return ((m_finer_level_ptr != nullptr) &&
            (m_finer_level_ptr->m_isDefined) &&
            (m_finer_level_ptr->m_level_grids.size() > 0));
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
int
AnisotropicAMRLevel::level() const
{
    return m_level;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
const IntVect&
AnisotropicAMRLevel::getFineRefRatio() const
{
    return (m_fineRefRatio);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
const IntVect&
AnisotropicAMRLevel::getCrseRefRatio() const
{
    return this->hasCoarserLevel() ? m_coarser_level_ptr->getFineRefRatio()
                                   : IntVect::Unit;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
AnisotropicAMRLevel::time(Real a_time)
{
    m_time = a_time;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real
AnisotropicAMRLevel::time() const
{
    return m_time;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
AnisotropicAMRLevel::initialDtMultiplier(Real a_initial_dt_multiplier)
{
    m_initial_dt_multiplier = a_initial_dt_multiplier;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real
AnisotropicAMRLevel::initialDtMultiplier() const
{
    return m_initial_dt_multiplier;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// static
void
AnisotropicAMRLevel::verbosity(int a_verbosity)
{
    s_verbosity = a_verbosity;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// static
int
AnisotropicAMRLevel::verbosity()
{
    return (s_verbosity);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
AnisotropicAMRLevel::preRegrid(int                         /*a_base_level*/,
                               const Vector<Vector<Box> >& /*a_new_grids*/)
{
    if (s_verbosity >= 3) {
        pout() << "AnisotropicAMRLevel::preRegrid" << endl;
    }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
AnisotropicAMRLevel::postRegrid(int /*a_base_level*/)
{
    if (s_verbosity >= 3) {
        pout() << "AnisotropicAMRLevel::postRegrid" << endl;
    }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
AnisotropicAMRLevel::postInitialGrid(const bool /*a_restart*/)
{
    if (s_verbosity >= 3) {
        pout() << "AnisotropicAMRLevel::postInitialGrid" << endl;
    }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Vector<AnisotropicAMRLevel*>
AnisotropicAMRLevel::getAMRLevelHierarchy()
{
    Vector<AnisotropicAMRLevel*> retval;
    // First go to level 0
    AnisotropicAMRLevel* levelPtr = this;
    while (levelPtr->hasCoarserLevel()) {
        levelPtr = levelPtr->m_coarser_level_ptr;
    }

    // Now can accumulate the pointers by chasing finer level
    retval.push_back(levelPtr);
    while (levelPtr->hasFinerLevel()) {
        levelPtr = levelPtr->m_finer_level_ptr;
        retval.push_back(levelPtr);
    }

    return retval;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
AnisotropicAMRLevel::writeCustomPlotFile(const std::string& /*a_prefix*/,
                                         int                /*a_step*/) const
{
    // By default, this does nothing.
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
AnisotropicAMRLevel::conclude(int /*a_step*/) const
{
    // By default, this does nothing.
}
//-----------------------------------------------------------------------
