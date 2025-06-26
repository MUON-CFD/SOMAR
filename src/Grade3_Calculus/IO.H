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
#ifndef ___IO_H__INCLUDED___
#define ___IO_H__INCLUDED___

#include <fstream>
#include "AMRIO.H"
#include "BasicIO.H"
#include "Debug.H"
#include "HeaderData.H"
#include "LevelGeometry.H"
#include "SOMAR_Constants.H"

#ifdef CH_USE_PYTHON
#   include "PyGlue.H"
#endif


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

    /// @brief Reads data from an ASCII file and puts it in the user's holder(s).
    /// @param a_ifstream
    /// @param a_v0
    /// @param a_vs...
    /// @details
    /// As an example, if the data within a file is:
    ///     a0 b0 c0
    ///     a1 b1 c1
    ///     a2 b2 c2
    ///     ...
    /// or like this
    ///     a0 b0 c0 a1 b1 c1 a2 b2 c2 ...
    /// then calling IO::readASCII(ifstream, a, b, c); would result in
    ///   a = {a0, a1, a2, ...}
    ///   b = {b0, b1, b2, ...}
    ///   c = {c0, c1, c2, ...}
    ///
    /// Data is read to the list of vectors {a, b, c, ...} in successive sweeps.
    /// If a sweep reaches the end of the file before all vectors in the list
    /// can be updated, then the entire sweep is discarded. This ensures all
    /// vectors are exactly the same size upon return.
    ///
    /// Data must be stored with space and/or newline separators, and this
    /// function makes no distinction between the two.
    /// There must be at least 1 entry to be read.
    ///
    /// It is up to you to check if the filestream is open. If the stream is
    /// closed, this function will quietly return.
    template <typename T0, typename... vTs>
    static void
    readASCII(std::ifstream&      a_ifstream,
              std::vector<T0>&    a_v0,
              vTs&...             a_vs);

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
    template <typename T>
    static void
    writeHDF5(std::string                      a_filename,
              const Vector<DisjointBoxLayout>& a_vectGrids,
              const Vector<T*>&                a_vectData,
              const Vector<std::string>&       a_compNames,
              const Box&                       a_domBox,
              const RealVect&                  a_dXi,
              const Real                       a_dt,
              const Real                       a_time,
              const Vector<IntVect>&           a_refRatios,
              const int                        a_numLevels);

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

// private:
//     template <typename _T>
//     class HasGhostVect
//     {
//         template <typename U>
//         static auto
//         check(int) -> decltype(std::declval<U>().ghostVect(), std::true_type{});
//         template <typename U>
//         static std::false_type
//         check(...);

//     public:
//         static constexpr bool value = decltype(check<_T>(0))::value;
//     };
};


// -----------------------------------------------------------------------------
template <typename T0, typename... vTs>
void
IO::readASCII(std::ifstream&      a_ifstream,
              std::vector<T0>&    a_v0,
              vTs&...             a_vs)
{
    // If you pass in N vectors, this function will read 1 element into each,
    // for a total of N reads.
    // It returns if there was an EOF error.
    auto _readOneLine = [](auto&&             self,
                           std::ifstream&     ifstream,
                           std::vector<auto>& v0,
                           auto&...           vs) -> bool {
        Real data;
        if (!(ifstream >> data)) {
            if (ifstream.eof()) return true;
            MAYDAYERROR("Unknown error while reading file.");
        }
        v0.push_back(static_cast<T0>(data));

        if constexpr (sizeof...(vs) > 0) return self(self, ifstream, vs...);
        return false;
    };

    // Crops all vectors to size N.
    // This assumes all vectors are at least size N.
    auto _trim = [](auto&&             self,
                    const size_t       N,
                    std::vector<auto>& v0,
                    auto&...           vs) {
        if (v0.size() > N) v0.resize(N);
        if constexpr (sizeof...(vs) > 0) self(self, N, vs...);
    };

    // This is the actual function...
    _trim(_trim, 0, a_v0, a_vs...);
    if (!a_ifstream.is_open()) return;
    while (!_readOneLine(_readOneLine, a_ifstream, a_v0, a_vs...)) {}
    const size_t minSize = std::min({a_v0.size(), a_vs.size()...});
    _trim(_trim, minSize, a_v0, a_vs...);
}


// -----------------------------------------------------------------------------
template <typename T> void
IO::writeHDF5 (std::string                           a_filename,
               const Vector<DisjointBoxLayout>&      /*a_vectGrids*/,
               const Vector<T* >& a_vectData,
               const Vector<std::string>&            a_compNames,
               const Box&                            a_domBox,
               const RealVect&                       a_dXi,
               const Real                            a_dt,
               const Real                            a_time,
               const Vector<IntVect>&                a_refRatios,
               const int                             a_numLevels)
{
    CH_assert(a_numLevels > 0);
    CH_assert(static_cast<int>(a_vectData.size())  >= a_numLevels);
    CH_assert(static_cast<int>(a_refRatios.size()) >= a_numLevels - 1);
    // if file exists, delete it
#ifdef CH_USE_PYTHON
    Py::PythonFunction("IO", "OpenFileForWrite", a_filename);
#endif
    // write some Chombo boilerplate (do we need it?)
    {
        HeaderData temp;
        temp.m_int["SpaceDim"] = SpaceDim;
        temp.m_real["testReal"] = (Real)0.0;
        temp.writeToFile(a_filename, std::string("Chombo_global"));
    }

    // Write file's header.
    HeaderData header;
    int nComp = a_compNames.size();

    string filedescriptor("VanillaAMRFileType");
    header.m_string["filetype"]       = filedescriptor;
    header.m_int   ["num_levels"]     = a_numLevels;
    header.m_int   ["num_components"] = nComp;
    header.m_int   ["max_level"]      = a_numLevels - 1;
    header.m_real  ["time"]           = a_time;

    for (int ivar = 0; ivar < nComp; ivar++) {
        char labelChSt[100];
        sprintf(labelChSt, "component_%d", ivar);
        string label(labelChSt);
        header.m_string[label] = a_compNames[ivar];
    }
    header.writeToFile(a_filename, std::string("/"));

    // Write data from level 0, up.
    Box domainLevel = a_domBox;
    Real dtLevel = a_dt;
    RealVect dXiLevel = a_dXi;
    for (int ilev = 0; ilev < a_numLevels; ilev++) {
        // Get ref ratio between this and finer levels.
        IntVect refLevel = IntVect::Unit;
        if (ilev != a_numLevels - 1) {
            refLevel = a_refRatios[ilev];
        }

        // Compute this levels domain and dXi by refining.
        if (ilev != 0) {
            domainLevel.refine(a_refRatios[ilev - 1]);
            dtLevel  /= a_refRatios[ilev - 1][0]; // HACK - just use 0 dir ref ratio
            dXiLevel /= a_refRatios[ilev - 1];
        }

        // Set this level's group name.
        // char levelName[20];
        // sprintf(levelName, "%d", ilev);
        // const std::string label = std::string("level_") + levelName;
        const std::string label = "level_" + std::to_string(ilev);

        // Write this level's header.
        HeaderData meta;
        meta.m_realvect["vec_dx"]      = dXiLevel;
        meta.m_real    ["dt"]          = dtLevel;
        meta.m_real    ["time"]        = a_time;
        meta.m_box     ["prob_domain"] = domainLevel;
        meta.m_intvect ["ref_ratio"]   = refLevel;

        meta.writeToFile(a_filename, label);


        // Write this level's BoxLayout.
#ifdef CH_USE_PYTHON
        CH_assert(a_vectData[ilev] != NULL);
        const T& dataLevel = *a_vectData[ilev];
        std::string name="data";
        std::string levelNameStr(label);
        if constexpr (std::is_same<T, LevelData<FArrayBox>>::value || std::is_same<T, LevelData<FluxBox>>::value) {
            Py::PythonFunction("IO",
                                  "WriteLevelDataFAB",
                                  a_filename,
                                  levelNameStr,
                                  name,
                                  dataLevel,
                                  dataLevel.ghostVect());
        }
        else
        {
        Py::PythonFunction("IO",
                              "WriteLevelDataFAB",
                              a_filename,
                              levelNameStr,
                              name,
                              dataLevel,
                              IntVect::Zero);
        }
#endif
    }
#ifdef CH_USE_PYTHON
    Py::PythonFunction("IO", "CloseFile", a_filename);
#endif
}


#endif //!___IO_H__INCLUDED___
