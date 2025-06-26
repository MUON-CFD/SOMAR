/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2024 Thomas Jefferson University and Arizona State University
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
 *  USA
 *
 *  For up-to-date contact information, please visit the repository homepage,
 *  https://github.com/MUON-CFD/SOMAR.
 ******************************************************************************/
#ifndef ___LayoutTools_H__INCLUDED___
#define ___LayoutTools_H__INCLUDED___

#include "LevelData.H"
#include "BaseFab.H"
#include "FluxBox.H"
#include "Copier.H"

namespace LayoutTools {


// Takes a DBL and returns a new one whose boxes are merged in a_dir.
void
mergeLayout(DisjointBoxLayout&       a_newLayout,
            const DisjointBoxLayout& a_origLayout,
            const int                a_dir);

// Splits a box into a load balanced set of boxes that are not split in a_dir.
void
lineLayout(DisjointBoxLayout&   a_newLayout,
           const ProblemDomain& a_domain,
           const int            a_dir);

// Define a dbl with just one box.
// This function must be called on ALL procs, but a_box only needs to be
// defined on a_srcProc.
void
defineOneProcGrids(DisjointBoxLayout&   a_grids,
                   const ProblemDomain& a_domain,
                   Box                  a_box,
                   const int a_srcProc = uniqueProc(SerialTask::compute));


/// This will define a copier that does not care about valid vs invalid data -
/// it will just copy everything.
/// This is a CC version.
void
defineImpartialCopier(Copier&              a_copier,
                      const BoxLayout&     a_srcLayout,
                      const BoxLayout&     a_destLayout,
                      const ProblemDomain& a_domain,
                      const IntVect&       a_srcGhostVect,
                      const IntVect&       a_destGhostVect,
                      const IntVect&       a_shift = IntVect::Zero);

/// This will define a copier that does not care about valid vs invalid data -
/// it will just copy everything.
/// This is a staggered version.
void
defineImpartialCopier(StaggeredCopier&     a_copier,
                      const int            a_fcDir,
                      const BoxLayout&     a_srcLayout,
                      const BoxLayout&     a_destLayout,
                      const ProblemDomain& a_domain,
                      const IntVect&       a_srcGhostVect,
                      const IntVect&       a_destGhostVect,
                      const IntVect&       a_shift = IntVect::Zero);


/// This will define a copier that does not care about valid vs invalid data -
/// it will just copy everything.
/// This is a staggered version.
void
defineImpartialCopier(std::array<StaggeredCopier, CH_SPACEDIM>& a_copier,
                      const BoxLayout&                          a_srcLayout,
                      const BoxLayout&                          a_destLayout,
                      const ProblemDomain&                      a_domain,
                      const IntVect&                            a_srcGhostVect,
                      const IntVect&                            a_destGhostVect,
                      const IntVect& a_shift = IntVect::Zero);

///
std::vector<Real>
domainDecompCoordVec(const std::vector<Real>& a_globalCoords,
                     const int                a_dir,
                     Box                      a_valid,
                     const Box&               a_domBox);


/// Grows each Box in a BoxLayout while preserving proc assignments, etc.
class GrowTransform: public BaseTransform
{
public:
    GrowTransform(const IntVect& a_grow);
    virtual Box operator()(const Box& a_inputBox);

protected:
    IntVect m_grow;
};


/// Changes the range in a single direction of all boxes in a BoxLayout
/// while preserving proc assignments, etc.
class RangeTransform: public BaseTransform
{
public:
    RangeTransform (const int a_dir, const int a_loEnd, const int a_numCells = 1);
    virtual ~RangeTransform ();
    virtual Box operator() (const Box& a_inputBox);

protected:
    int m_dir;
    int m_loEnd;
    int m_numCells;
};


/// Shifts all boxes in a BoxLayout while preserving proc assignments, etc.
class ShiftTransform: public BaseTransform
{
public:
    ShiftTransform (const IntVect& a_shift);
    ShiftTransform (const int a_dir, const int a_shift);
    virtual ~ShiftTransform ();
    virtual Box operator() (const Box& a_inputBox);

protected:
    IntVect m_shift;
};


/// Shifts all boxes in a BoxLayout while preserving proc assignments, etc.
class ShiftHalfTransform: public BaseTransform
{
public:
    ShiftHalfTransform (const IntVect& a_shift);
    ShiftHalfTransform (const int a_dir, const int a_shift);
    virtual ~ShiftHalfTransform ();
    virtual Box operator() (const Box& a_inputBox);

protected:
    IntVect m_shift;
};


/// Slices all boxes in a BoxLayout while preserving proc assignments, etc.
class SliceTransform: public BaseTransform
{
public:
    SliceTransform (const int a_dir, const int a_pos);
    virtual ~SliceTransform ();
    virtual Box operator() (const Box& a_inputBox);

protected:
    int m_dir, m_pos;
};


/// Takes each Box in the BoxLayout and converts it to it's counterpart in a_newBoxes.
class SwapTransform : public BaseTransform
{
public:
    /// @brief The two input vectors must be 1-to-1. a_oldBoxes[i] will be
    /// transformed into a_newBoxes[i].
    /// @param a_oldBoxes
    /// @param a_newBoxes
    SwapTransform(const std::vector<Box>& a_oldBoxes,
                  const std::vector<Box>& a_newBoxes);
    virtual ~SwapTransform() = default;
    virtual Box operator()(const Box& a_inputBox) override;

protected:
    std::vector<Box> m_oldBoxes;
    std::vector<Box> m_newBoxes;
};



/// This shifts shiftFAB so that its lower vertical end coincides with the lower
/// vertical end of target*. The centering will also be changed. You MUST call
/// restore() before this object goes out of scope or an error will be thrown.
template<typename T>
class VertShifter
{
    typedef BaseFab<T> DataType;

public:
    VertShifter (DataType&       a_shiftFAB,
                 const DataType& a_targetFAB);
    VertShifter (DataType&  a_shiftFAB,
                 const Box& a_targetBox);
    ~VertShifter ();
    void restore ();

private:
    // This is private so that we can only use this object once.
    void shift (DataType&  a_shiftFAB,
                const Box& a_targetBox);

    int  m_totalShift;
    int  m_typeIsChanged;
    DataType* m_dataPtr;
};


/// \brief Checks overlapping faces to ensure data is non-NaN and identical.
///
/// \param[in] a_data The staggered data to be checked.
/// \param[in] a_file Pass in __FILE__ so the error msg pinpoints where the error occurred.
/// \param[in] a_line Pass in __LINE__ so the error msg pinpoints where the error occurred.
void
checkValidFaceOverlap(const LevelData<FluxBox>& a_data,
                      const char*               a_file,
                      const size_t              a_line);

/// \brief Checks overlapping faces to ensure data is non-NaN and identical.
///
/// \param[in] a_data  A face-centered LevelData.
/// \param[in] a_fcDir The face-centered direction.
/// \param[in] a_file  Pass in __FILE__ so the error msg pinpoints where the error occurred.
/// \param[in] a_line  Pass in __LINE__ so the error msg pinpoints where the error occurred.
void
checkValidFaceOverlap(const LevelData<FArrayBox>& a_data,
                      const int                   a_fcDir,
                      const char*                 a_file,
                      const size_t                a_line);

/// @brief This version writes discrepancies to HDF5.
/// @param a_filename
/// @param a_data
void
checkValidFaceOverlap(const LevelData<FluxBox>& a_data,
                      LevelData<FluxBox>&       a_destDiscrepancies);

/// Small discrepancies (the size of rounding errors) at overlying faces on
/// different grids can lead to large discrepancies after many timesteps.
/// If left unchecked, these errors magnify and take down some simulations.
/// Averaging the overlying values seems to help.
void
averageOverlappingValidFaces(LevelData<FluxBox>& a_data);


// Include templated definitions
#define Me949b7e6032829e355cb3ebad4d68ff6
#   include "LayoutToolsI.H"
#undef Me949b7e6032829e355cb3ebad4d68ff6

}; // end LayoutTools namespace
#endif //!___LayoutTools_H__INCLUDED___
