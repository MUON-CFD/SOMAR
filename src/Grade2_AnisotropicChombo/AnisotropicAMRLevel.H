/*******************************************************************************
 *    FILE: AnisotropicAMRLevel.H
 *    DESC: Same as AnisotropicAMRLevel, but allows for anisotropic refinement.
 *    DATE: Wed 14 Aug 2013 06:58:45 PM EDT
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __AnisotropicAMRLevel_H__INCLUDED__
#define __AnisotropicAMRLevel_H__INCLUDED__

#include "Box.H"
#include "CH_HDF5.H"
#include "DataIterator.H"
#include "IntVectSet.H"
#include "ProblemDomain.H"
#include "REAL.H"
#include "Vector.H"
#include "HeaderData.H"
// class HDF5Handle;
// class IntVectSet;
// template <class TTYPE> class Vector;

/// Abstract base class for time-dependent data at a level of refinement.
/**
   This is an abstract base class for data at the same level of
   refinement within a hierarchy of levels.  The concrete class
   derived from AnisotropicAMRLevel is called a  physics class.  The domain
   of a level is a disjoint union of rectangles in a logically
   rectangular index space.  Data is defined within this domain.
   There is also a problem domain, which may be larger, within which
   data can, in theory, be interpolated from some coarser level.

   AnisotropicAMRLevel is intimately entangled with the AMR class.
   The AMR contains a collection of AnisotropicAMRLevel objects which
   represent a hierarchy of levels of refinement.  The class  AMR
   is responsible for calling the correct member functions of
   AnisotropicAMRLevel.  The user is responsible for implementing the physics
   class, and might never call any AnisotropicAMRLevel functions described
   here.
*/

class AnisotropicAMRLevel
{
public:
    ///
    /**
       Default constructor.
    */
    AnisotropicAMRLevel();

    ///
    /**
       Destructor.
    */
    virtual ~AnisotropicAMRLevel();

    ///
    /**
       Defines this AnisotropicAMRLevel.

       <ul>
         <li>a_coarser_level_ptr (not modified): pointer to next coarser level
             object.
         <li>a_problem_domain (not modified): problem domain of this level.
         <li>a_level (not modified): index of this level.  The base level is
             zero.
         <li>a_ref_ratio (not modified): the refinement ratio between this level
             and the next finer level.
       </ul>
    */
    virtual void
    define(AnisotropicAMRLevel* a_coarser_level_ptr,
           const Box&           a_problem_domain,
           int                  a_level,
           const IntVect&       a_ref_ratio);

    ///
    /**
       Defines this AnisotropicAMRLevel.

       <ul>
         <li>a_coarser_level_ptr (not modified): pointer to next coarser level
             object.
         <li>a_problem_domain (not modified): problem domain of this level.
         <li>a_level (not modified): index of this level.  The base level is
             zero.
         <li>a_ref_ratio (not modified): the refinement ratio between this level
             and the next finer level.
       </ul>
    */
    virtual void
    define(AnisotropicAMRLevel* a_coarser_level_ptr,
           const ProblemDomain& a_problem_domain,
           int                  a_level,
           const IntVect&       a_ref_ratio);

    ///
    /**
       Advances this level by one time step.  Returns an estimate of the
       new time step.

       This is a pure virtual function and MUST be defined in the derived
       class.

    */
    virtual Real
    advance() = 0;

    ///
    /**
       Things to do after advancing this level by one time step.

       This is a pure virtual function and  MUST be defined in the derived
       class.

       a_step is the step number that is being completed.
       Unlike in AnisotropicAMR.cpp, we start counting from 1.

    */
    virtual void
    postLevelTimeStep(int a_step) = 0;

    ///
    /**
       Creates tagged cells for dynamic mesh refinement.

       This is a pure virtual function and  MUST be defined in the derived
       class.

    */
    virtual void
    tagCells(IntVectSet& a_tags) = 0;

    ///
    /**
       Creates tagged cells for mesh refinement at initialization.

       This is a pure virtual function and MUST be defined in the derived
       class.

    */
    virtual void
    tagCellsInit(IntVectSet& a_tags) = 0;

    ///
    /**
        Performs any pre-regridding operations which are necessary.

        This is not a pure virtual function to preserve compatibility
        with earlier versions of AnisotropicAMRLevel.  The
       AnisotropicAMRLevel::preRegrid() instantiation is a no-op.
    */
    virtual void
    preRegrid(int a_base_level, const Vector<Vector<Box> >& a_new_grids);

    ///
    /**
       Redefines this level to have the specified domain a_new_grids.

       This is a pure virtual function and MUST be defined in the derived
       class.

    */
    virtual void
    regrid(const Vector<Box>& a_new_grids) = 0;

    ///
    /**
        Performs any post-regridding operations which are necessary.

        This is not a pure virtual function to preserve compatibility
        with earlier versions of AnisotropicAMRLevel.  The
       AnisotropicAMRLevel::postRegrid() instantiation is a no-op.
    */
    virtual void
    postRegrid(int a_base_level);

    ///
    /**
       Initializes this level to have the specified domain a_new_grids.

       This is a pure virtual function and MUST be defined in the derived
       class.

    */
    virtual void
    initialGrid(const Vector<Box>& a_new_grids) = 0;

    ///
    /**
       Performs operations required after the grid has been defined but
       before data initialization.  This will also be called after
       readCheckpointLevel during a restart procedure with argument
       a_restart set to true.

       Levels are accessed from finest to coarsest.  The
       AnisotropicAMRLevel::postInitialGrid() instantiation is a no-op.

    */
    virtual void
    postInitialGrid(const bool a_restart);

    ///
    /**
       Initializes data.

       This is a pure virtual function and  MUST be defined in the derived
       class.

    */
    virtual void
    initialData() = 0;

    ///
    /**
       Things to do after initialization.

       This is a pure virtual function and  MUST be defined in the derived
       class.

    */
    virtual void
    postInitialize() = 0;

    //! Override this method to have an AnisotropicAMRLevel subclass perform
    //! some operation upon the conclusion of a simulation. This is called when
    //! AMR::conclude() is called. The final step is passed to the method.
    //! \param a_step The last step in the simulation.
    virtual void
    conclude(int a_step) const;

    /**
       \name I/O functions
    */
    /**@{*/


    ///
    /**
       Writes checkpoint header.

       This is a pure virtual function and  MUST be defined in the derived
       class.

    */
    virtual void
    writeCheckpointHeader(HeaderData & a_header, const std::string & a_filename) const = 0;

    ///
    /**
       Write checkpoint data for this level.

       This is a pure virtual function and  MUST be defined in the derived
       class.

    */
    virtual void
    writeCheckpointLevel(const std::string & a_filename, int level) const = 0;


    virtual void
    openFile(const std::string & a_filename, const bool checkpoint) const = 0;

    virtual void
    closeFile(const std::string & a_filename) const = 0;

    virtual void
    writePlotHeader(HeaderData & a_handle, const std::string & a_filename) const = 0;

    ///
    /**
       Write plot file for this level.

       This is a pure virtual function and MUST be defined in the derived
       class.
    */

    virtual void
    writePlotLevel(const std::string & a_filename, int level) const = 0;

///
    /**
       Reads checkpoint header.

       This is a pure virtual function and  MUST be defined in the derived
       class.

    */
    virtual void
    readCheckpointHeader(const std::string & a_fileName) = 0;

    ///
    /**
       Reads checkpoint data for this level.

       This is a pure virtual function and MUST be defined in the derived
       class.

    */
    virtual void
    readCheckpointLevel(const std::string & a_fileName) = 0;

    ///
    /**
       Writes plot header.

       This is a pure virtual function and MUST be defined in the derived
       class.
*/
    //! This allows one to write a plot file in a non-HDF5 format. It is called
    //! only at refinement level 0, so AMR data will have to be handled by the
    //! implementer. \param a_prefix A prefix for the custom plot file name.
    //! \param a_step The current time step.
    virtual void
    writeCustomPlotFile(const std::string& a_prefix, int a_step) const;

    /**@}*/

    /**
       \name Parameter-setting functions
    */
    /**@{*/

    ///
    /**
       Sets the pointer-to-finer-level member to a_finer_level_ptr.
    */
    virtual void
    finerLevelPtr(AnisotropicAMRLevel* a_finer_level_ptr);

    ///
    /**
       Sets the time step to a_dt.
    */
    virtual void
    dt(Real a_dt);

    ///
    /**
       Sets the time to a_time.

    */
    virtual void
    time(Real a_time);

    ///
    /**
       Sets the initial dt multiplier to a_initial_dt_multiplier.
    */
    virtual void
    initialDtMultiplier(Real a_initial_dt_multiplier);

    /**@}*/

    /**
       \name Access functions
    */
    /**@{*/

    ///
    /**
       Returns the current value of the time step.

    */
    virtual Real
    dt() const;

    ///
    /**
       Returns the current value of the time on this level.
    */
    virtual Real
    time() const;

    ///
    /**
       Returns the initial dt multiplier.
    */
    virtual Real
    initialDtMultiplier() const;

    ///
    /**
       Returns the problem domain of this level.

    */
    virtual const ProblemDomain&
    problemDomain() const;

    ///
    /**
       Returns the domain of this level.

    */
    virtual Vector<Box>
    boxes() const;

    ///
    /**
       Returns true if any AnisotropicAMRLevel::define function has been called,
       false otherwise.
    */
    bool
    isDefined() const;

    ///
    /**
       Returns true if a coarser level exists, is defined, and has a grid.
    */
    bool
    hasCoarserLevel() const;

    ///
    /**
       Returns true if a finer level exists, is defined, and has a grid.
    */
    bool
    hasFinerLevel() const;

    ///
    /**
       Returns the index of this level
    */
    virtual int
    level() const;

    ///
    /**
       Gets the refinement ratio between this level and the next finer level.
       Returns IntVect::Unit if this is the finest level.
    */
    virtual const IntVect&
    getFineRefRatio() const;

    ///
    /**
       Gets the refinement ratio between this level and the next coarser level.
       Returns IntVect::Unit if this is the coarsest level.
    */
    virtual const IntVect&
    getCrseRefRatio() const;

    ///
    /**
       Returns maximum stable time step for this level.

       This is a pure virtual function and  MUST be defined in the derived
       class.

    */
    virtual Real
    computeDt() = 0;

    ///
    /**
        Returns maximum stable time step for this level with initial data.

        This is a pure virtual function and  MUST be defined in the derived
        class.

    */
    virtual Real
    computeInitialDt() = 0;

    //! Retrieve an array of all of the AnisotropicAMRLevel objects in the
    //! entire hierarchy.
    Vector<AnisotropicAMRLevel*>
    getAMRLevelHierarchy();

    /**@}*/

    ///
    /**
       Returns current verbosity level.  Minimum verbosity is 0, for
       which nothing is printed.

    */
    static int
    verbosity();

    ///
    /**
       Sets verbosity level to a_verbosity.  Minimum verbosity is 0, for
       which nothing is printed.
    */
    static void
    verbosity(int a_verbosity);

protected:
    // verbosity level
    static int s_verbosity;

    // the problem domain
    ProblemDomain m_problem_domain;

    //
    Vector<Box> m_level_grids;

    // the level
    int m_level;

    // refinement ratios
    IntVect m_fineRefRatio;

    // initial time step multipier
    Real m_initial_dt_multiplier;

    // time step
    Real m_dt;

    // time
    Real m_time;

    // pointer to next coarser level
    AnisotropicAMRLevel* m_coarser_level_ptr;

    // pointer to next finer level
    AnisotropicAMRLevel* m_finer_level_ptr;

    bool m_isDefined;
};


#endif  //!__AnisotropicAMRLevel_H__INCLUDED__
