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
 *  https://github.com/MUON-CFD/somar.
 ******************************************************************************/
#ifndef ___IO_H__INCLUDED___
#define ___IO_H__INCLUDED___

#include <fstream>
#include "BasicIO.H"
#include "SOMAR_Constants.H"
#include "LevelGeometry.H"


class IO: public BasicIO
{
public:
    // Reads data from an ASCII file and puts it in the user's holder.
    // Data must be stored with space and/or newline separators.
    // There must be at least 1 entry to be read.
    //
    // It is up to you to check if the filestream is open. If the stream is
    // closed, this function will quietly return.
    static void
    readASCII(int&           a_data,
              std::ifstream& a_ifstream);

    // Reads data from an ASCII file and puts it in the user's holder.
    // Data must be stored with space and/or newline separators.
    // There must be at least 1 entry to be read.
    //
    // It is up to you to check if the filestream is open. If the stream is
    // closed, this function will quietly return.
    static void
    readASCII(size_t&        a_data,
              std::ifstream& a_ifstream);

    // Reads data from an ASCII file and puts it in the user's holder.
    // Data must be stored with space and/or newline separators.
    // There must be at least 1 entry to be read.
    //
    // It is up to you to check if the filestream is open. If the stream is
    // closed, this function will quietly return.
    static void
    readASCII(Real&          a_data,
              std::ifstream& a_ifstream);

    // Reads data from an ASCII file and puts it in the user's holder.
    // Data must be stored with space and/or newline separators.
    // There must be at least a_data.size() entries to be read.
    //
    // It is up to you to check if the filestream is open and vector is properly
    // sized. If the stream is closed or the vector has zero length, this
    // function will quietly return.
    //
    // Data from a Matlab/Octave vector, v, can be read via this function if
    // written to file like this...
    //
    // fid = fopen('filename.txt', 'W');
    // fprintf(fid, '%f ', v(:));
    // fclose(fid);
    static void
    readASCII(std::vector<Real>& a_data,
              std::ifstream&     a_ifstream);

    // Reads data from an ASCII file and puts it in the user's holder.
    // Data must be stored with space and/or newline separators.
    // There must be at least a_dataBox.size() entries to be read.
    //
    // It is up to you to check if the filestream is open. If it isn't, this
    // function will quietly return.
    //
    // Data from a Matlab/Octave matrix, A, can be read via this function if
    // written to file like this...
    //
    // fid = fopen('filename.txt', 'W');
    // for k = 1:NZ
    //   for j = 1:NY   # if 3D
    //     fprintf(fid, '%f ', A'(:,j,k));
    //     fprintf(fid, '\n');
    //   end            # if 3D
    // end
    // fclose(fid);
    static void
    readASCII(FArrayBox&     a_dataFAB,
              const int      a_dataComp,
              const Box&     a_dataBox,
              std::ifstream& a_ifstream);

    // Reads data from an ASCII file and puts it in the user's holder.
    // Data must be stored with space and/or newline separators.
    // There must be at least a_data.getBoxes().physDomain().size() entries to
    // be read.
    //
    // It is up to you to check if the filestream is open. If it isn't, this
    // function will quietly return.
    //
    // Data from a Matlab/Octave matrix, A, can be read via this function if
    // written to file like this...
    //
    // fid = fopen('filename.txt', 'W');
    // for k = 1:NZ
    //   for j = 1:NY   # if 3D
    //     fprintf(fid, '%f ', A'(:,j,k));
    //     fprintf(fid, '\n');
    //   end            # if 3D
    // end
    // fclose(fid);
    static void
    readASCII(LevelData<FArrayBox>& a_data,
              const int             a_dataComp,
              std::ifstream&        a_ifstream);

    // Writes hierarchy of levels in HDF5 format.  Only available if the
    // preprocessor macro HDF5 is defined at compilation.
    //
    // Arguments:
    //   a_filename  : file to output to.
    //   a_vectGrids : grids at each level.
    //   a_vectData  : data at each level.
    //   a_compNames : names of variables.
    //   a_domBox    : domain at coarsest level.
    //   a_dXi       : grid spacing in each direction at coarsest level.
    //   a_dt        : time step at coarsest level.
    //   a_time      : the current time.
    //   a_refRatios : refinement ratio in each direction at all levels
    //                 (ith entry is refinement ratio in each direction
    //                 between levels i and i + 1).
    //   a_numLevels : number of levels to output.
    //
    // This is blocking.
    static void
    writeHDF5 (std::string                           a_filename,
               const Vector<DisjointBoxLayout>&      a_vectGrids,
               const Vector<LevelData<FArrayBox>* >& a_vectData,
               const Vector<std::string>&            a_compNames,
               const Box&                            a_domBox,
               const RealVect&                       a_dXi,
               const Real                            a_dt,
               const Real                            a_time,
               const Vector<IntVect>&                a_refRatios,
               const int                             a_numLevels);

    // Writes hierarchy of levels in HDF5 format.  Only available if the
    // preprocessor macro HDF5 is defined at compilation.
    //
    // This is the less general version that you will probably want to use
    // most of the time for multi-level data.
    //
    // Arguments:
    //   a_filename  : file to output to.
    //   a_vectData  : data at each level.
    //   a_levGeo    : A LevelGeometry object at any level in the hierarchy.
    //                 It must be able to produce all other LevelGeometry
    //                 objects via getCoarserPtr() and getFinerPtr().
    //   a_dt        : time step at coarsest level.
    //   a_time      : the current time.
    //   a_compNames : names of variables.
    //
    // This is blocking.
    static void
    writeHDF5 (std::string                           a_filename,
               const Vector<LevelData<FArrayBox>* >& a_vectData,
               const LevelGeometry&                  a_levGeo,
               const Real                            a_dt        = BOGUS_TIME,
               const Real                            a_time      = BOGUS_TIME,
               const Vector<std::string>&            a_compNames = Vector<std::string>());

    // Writes a single level of data in HDF5 format.  Only available if the
    // preprocessor macro HDF5 is defined at compilation.
    //
    // This is the less general version that you will probably want to use
    // most of the time for single-level data.
    //
    // Arguments:
    //   a_filename  : file to output to.
    //   a_data      : a single level of data.
    //   a_levGeo    : this level's LevelGeometry object.
    //   a_dt        : time step.
    //   a_time      : the current time.
    //   a_compNames : names of variables.
    //
    // This is blocking.
    static void
    writeHDF5 (std::string                 a_filename,
               const LevelData<FArrayBox>& a_data,
               const LevelGeometry&        a_levGeo,
               const Real                  a_dt        = BOGUS_TIME,
               const Real                  a_time      = BOGUS_TIME,
               const Vector<std::string>&  a_compNames = Vector<std::string>());


    // Writes hierarchy of levels in HDF5 format.  Only available if the
    // preprocessor macro HDF5 is defined at compilation.
    //
    // This converts a LevelData<FluxBox> to cell-centered and calls the
    // LevelData<FArrayBox> version.
    //
    // This is the less general version that you will probably want to use
    // most of the time for multi-level, face-centered data.
    //
    // Arguments:
    //   a_filename  : file to output to.
    //   a_vectData  : data at each level.
    //   a_levGeo    : A LevelGeometry object at any level in the hierarchy.
    //                 It must be able to produce all other LevelGeometry
    //                 objects via getCoarserPtr() and getFinerPtr().
    //   a_dt        : time step at coarsest level.
    //   a_time      : the current time.
    //   a_compNames : names of variables.
    //
    // This is blocking.
    static void
    writeHDF5 (std::string                         a_filename,
               const Vector<LevelData<FluxBox>* >& a_vectData,
               const LevelGeometry&                a_levGeo,
               const Real                          a_dt        = BOGUS_TIME,
               const Real                          a_time      = BOGUS_TIME,
               const Vector<std::string>&          a_compNames = Vector<std::string>());

    // Writes a single level of data in HDF5 format.  Only available if the
    // preprocessor macro HDF5 is defined at compilation.
    //
    // This converts a LevelData<FluxBox> to cell-centered and calls the
    // LevelData<FArrayBox> version.
    //
    // This is the less general version that you will probably want to use
    // most of the time for single-level, face-centered data.
    //
    // Arguments:
    //   a_filename  : file to output to.
    //   a_data      : a single level of data.
    //   a_levGeo    : this level's LevelGeometry object.
    //   a_dt        : time step.
    //   a_time      : the current time.
    //   a_compNames : names of variables.
    //
    // This is blocking.
    static void
    writeHDF5 (std::string                a_filename,
               const LevelData<FluxBox>&  a_data,
               const LevelGeometry&       a_levGeo,
               const Real                 a_dt        = BOGUS_TIME,
               const Real                 a_time      = BOGUS_TIME,
               const Vector<std::string>& a_compNames = Vector<std::string>());

    // // This version does not require a levGeo.
    // static void
    // writeHDF5 (std::string                a_filename,
    //            const LevelData<FluxBox>&  a_data,
    //            const Real                 a_dt        = BOGUS_TIME,
    //            const Real                 a_time      = BOGUS_TIME,
    //            const Vector<std::string>& a_compNames = Vector<std::string>());
};



#endif //!___IO_H__INCLUDED___
