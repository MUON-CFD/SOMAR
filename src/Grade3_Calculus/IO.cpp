/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2017
 *    Jefferson (Philadelphia University + Thomas Jefferson University) and
 *    University of North Carolina at Chapel Hill
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
 *  https://github.com/somarhub.
 ******************************************************************************/
#include <sstream>
#include <fstream>

#include "IO.H"
#include "Format.H"
#include "AMRIO.H"
#include "HeaderData.H"
#include "Convert.H"
#include "LayoutTools.H"
#include "Debug.H"

#ifdef CH_USE_PYTHON
#include "PyGlue.H"
#endif


// -----------------------------------------------------------------------------
void
IO::readASCII(int& a_data, std::ifstream& a_ifstream)
{
    Real tmp;
    IO::readASCII(tmp, a_ifstream);
    a_data = static_cast<int>(tmp);
}


// -----------------------------------------------------------------------------
void
IO::readASCII(size_t& a_data, std::ifstream& a_ifstream)
{
    Real tmp;
    IO::readASCII(tmp, a_ifstream);
    CH_verify(tmp >= 0.0);
    a_data = static_cast<size_t>(tmp);
}


// -----------------------------------------------------------------------------
void
IO::readASCII(Real& a_data, std::ifstream& a_ifstream)
{
    if (!a_ifstream.is_open()) return;

    if (!(a_ifstream >> a_data)) {
        if (a_ifstream.eof()) {
            MAYDAYERROR("File does not contain the requested Real.");
        } else {
            MAYDAYERROR("Unknown error while reading file.");
        }
    }
}


// -----------------------------------------------------------------------------
void
IO::readASCII(std::vector<Real>& a_data, std::ifstream& a_ifstream)
{
    if (!a_ifstream.is_open() || a_data.empty()) return;

    for (size_t i = 0; i < a_data.size(); ++i) {
        if (!(a_ifstream >> a_data[i])) {
            if (a_ifstream.eof()) {
                MAYDAYERROR(
                    "File does not contain enough data to fill vector.");
            } else {
                MAYDAYERROR("Unknown error while reading file.");
            }
        }
    }
}

// -----------------------------------------------------------------------------
void
IO::readASCII(FArrayBox&     a_dataFAB,
              const int      a_dataComp,
              const Box&     a_dataBox,
              std::ifstream& a_ifstream)
{
    if (!a_ifstream.is_open()) return;

    CH_assert(0 <= a_dataComp && a_dataComp < a_dataFAB.nComp());
    CH_assert(a_dataFAB.box().contains(a_dataBox));

    const IntVect& sm = a_dataBox.smallEnd();
    const IntVect& bg = a_dataBox.bigEnd();

    IntVect iv = sm;
#if CH_SPACEDIM > 2
    for (iv[2] = sm[2]; iv[2] <= bg[2]; ++iv[2]) {
#endif
        for (iv[1] = sm[1]; iv[1] <= bg[1]; ++iv[1]) {
            for (iv[0] = sm[0]; iv[0] <= bg[0]; ++iv[0]) {
                if (!(a_ifstream >> a_dataFAB(iv, a_dataComp))) {
                    if (a_ifstream.eof()) {
                        MAYDAYERROR(
                            "File does not contain enough data to fill FAB.");
                    } else {
                        MAYDAYERROR("Unknown error while reading file.");
                    }
                }
            }
        }
#if CH_SPACEDIM > 2
    }
#endif
}

// -----------------------------------------------------------------------------
void
IO::readASCII(LevelData<FArrayBox>& a_data,
              const int             a_dataComp,
              std::ifstream&        a_ifstream)
{
    if (!a_ifstream.is_open()) return;

    const DisjointBoxLayout& grids  = a_data.getBoxes();
    const ProblemDomain&     domain = grids.physDomain();
    const Box&               domBox = domain.domainBox();

    // Create space on a single proc.
    DisjointBoxLayout oneProcGrids;
    LayoutTools::defineOneProcGrids(oneProcGrids, domain, domBox);
    LevelData<FArrayBox> oneProcData(oneProcGrids, 1);

    // Read data on that single proc.
    for (DataIterator dit(oneProcGrids); dit.ok(); ++dit) {
        IO::readASCII(oneProcData[dit], a_dataComp, domBox, a_ifstream);
    }

    // Copy to user's holder / Distribute data.
    oneProcData.copyTo(Interval(0,0), a_data, Interval(a_dataComp, a_dataComp));
}


// -----------------------------------------------------------------------------
// Writes hierarchy of levels in HDF5 format.  Only available if the
// preprocessor macro HDF5 is defined at compilation.
//
// Arguments:
//   a_filename  : name of output file.
//   a_vectGrids : grids at each level.
//   a_vectData  : data at each level.
//   a_compNames : names of variables.
//   a_domBox    : domain at coarsest level.
//   a_dXi       : grid spacing in each direction at coarsest level.
//   a_dt        : time step at coarsest level.
//   a_time      : the current time.
//   a_vectRatio : refinement ratio in each direction at all levels
//                 (ith entry is refinement ratio in each direction
//                 between levels i and i + 1).
//   a_numLevels : number of levels to output.
//
// This is blocking.
// -----------------------------------------------------------------------------
void
IO::writeHDF5 (std::string                           a_filename,
               const Vector<DisjointBoxLayout>&      /*a_vectGrids*/,
               const Vector<LevelData<FArrayBox>* >& a_vectData,
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
        const LevelData<FArrayBox>& dataLevel = *a_vectData[ilev];
        std::string name="data";
        std::string levelNameStr(label);

        Py::PythonFunction("IO",
                              "WriteLevelDataFAB",
                              a_filename,
                              levelNameStr,
                              name,
                              dataLevel,
                              a_vectData[0]->ghostVect());
#endif
    }
#ifdef CH_USE_PYTHON
    Py::PythonFunction("IO", "CloseFile", a_filename);
#endif
}


// -----------------------------------------------------------------------------
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
// -----------------------------------------------------------------------------
void
IO::writeHDF5 (std::string                           a_filename,
               const Vector<LevelData<FArrayBox>* >& a_vectData,
               const LevelGeometry&                  a_levGeo,
               const Real                            /*a_dt*/,
               const Real                            a_time,
               const Vector<std::string>&            a_compNames)
{
    // Gather geometric data
    const int                          finestLevel = a_vectData.size() - 1;
    const Vector<const LevelGeometry*> vLevGeo     = a_levGeo.getAMRLevGeos();
    const Vector<DisjointBoxLayout>    vGrids      = a_levGeo.getAMRGrids();
    const Vector<IntVect>              vRefRatio   = a_levGeo.getAMRRefRatios();
    const ProblemDomain&               lev0Domain  = vLevGeo[0]->getDomain();
    const RealVect                     lev0dXi     = vLevGeo[0]->getDXi();

    // These used to be input parameters, but I think they just made things
    // complicated. I'm leaving them hard-coded to write the entire hierarchy.
    // Maybe someday I'll have a need to make these function inputs again.
    const int a_lmin = 0;
    const int a_lmax = finestLevel;
    CH_assert(0 <= a_lmin);
    CH_assert(a_lmin <= a_lmax);
    CH_assert(a_lmax <= finestLevel);

    // How many comps will we need to copy? How many ghosts are there?
    int numComps = -1;
    IntVect ghostVect = IntVect::Unit; // At least one for proper display of mapping.
    for (int lev = 0; lev <= finestLevel; ++lev) {
        if (a_vectData[lev] == NULL) continue;

        int thisNumComps = a_vectData[lev]->nComp();
        if (numComps == -1) numComps = thisNumComps;
        CH_assert(thisNumComps == numComps);

        IntVect thisGhostVect = a_vectData[lev]->ghostVect();
        D_TERM(ghostVect[0] = Max(ghostVect[0], thisGhostVect[0]);,
               ghostVect[1] = Max(ghostVect[1], thisGhostVect[1]);,
               ghostVect[2] = Max(ghostVect[2], thisGhostVect[2]);)
    }
    CH_assert(numComps > 0);

    // Package output into one CC holder
    Vector<LevelData<FArrayBox>*> vOutput(finestLevel+1, NULL);
    for (int lev = 0; lev <= finestLevel; ++lev) {
        const DisjointBoxLayout& grids = vGrids[lev];
        DataIterator dit = grids.dataIterator();

        // No matter what, initialize output and fill the displacement field.
        vOutput[lev] = new LevelData<FArrayBox>(grids, numComps+SpaceDim, ghostVect);
        Interval dispInt(numComps, numComps+SpaceDim-1);
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox& outputFAB = (*vOutput[lev])[dit];
            outputFAB.setVal(quietNAN);

            FArrayBox dispFAB(dispInt, (*vOutput[lev])[dit]);
            vLevGeo[lev]->fill_displacement(dispFAB);
        }

        // Do we want this level's data?
        if (a_vectData[lev] == NULL) continue;
        if (lev < a_lmin) continue;
        if (a_lmax < lev) continue;

        // Copy the data
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox&       outputFAB = (*vOutput[lev])[dit];
            const FArrayBox& dataFAB   = (*a_vectData[lev])[dit];
            const Interval   ivl(0, numComps - 1);

            Box destBox = dataFAB.box();
            destBox.enclosedCells();
            destBox &= outputFAB.box();

            Convert::Simple(outputFAB, ivl, destBox, dataFAB, ivl);
        }
    }

    // Write to file
    {
        Vector<string> newCompNames(numComps + SpaceDim);
        int comp = 0;
        if (a_compNames.size() == 0) {
            if (numComps == SpaceDim) {
                // It's probably a vector.
                D_TERM(newCompNames[comp++] = "x_component";,
                       newCompNames[comp++] = "y_component";,
                       newCompNames[comp++] = "z_component";)
            } else {
                // It's probably not a vector.
                while (comp < numComps) {
                    char compName[40];
                    sprintf(compName, "comp_%d", comp);
                    newCompNames[comp++] = std::string(compName);
                }
            }
        } else {
            // The caller provided names.
            while (comp < numComps) {
                newCompNames[comp] = a_compNames[comp];
                ++comp;
            }
        }
        D_TERM(
        newCompNames[comp++] = "x_Displacement";,
        newCompNames[comp++] = "y_Displacement";,
        newCompNames[comp++] = "z_Displacement";)
        CH_assert(comp == numComps+SpaceDim);

        Real dt = 1.0;
        Real dummyTime = ((a_time < 0.0)? 0.0: a_time);

        IO::writeHDF5(a_filename,
                      vGrids,
                      vOutput,
                      newCompNames,
                      lev0Domain.domainBox(),
                      lev0dXi,
                      dt,
                      dummyTime,
                      vRefRatio,
                      vOutput.size());

        barrier();

    }

    // Free memory
    for (int lev = 0; lev <= finestLevel; ++lev) {
        delete vOutput[lev];
    }

}


// -----------------------------------------------------------------------------
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
// -----------------------------------------------------------------------------
void
IO::writeHDF5 (std::string                 a_filename,
               const LevelData<FArrayBox>& a_data,
               const LevelGeometry&        a_levGeo,
               const Real                  a_dt,
               const Real                  a_time,
               const Vector<std::string>&  a_compNames)
{
    // We need to find where we are in the AMR hierarchy.
    const Vector<const LevelGeometry*> vLevGeo = a_levGeo.getAMRLevGeos();
    const int numLevels = a_levGeo.getNumLevels();

    int lev;
    for (lev = 0; lev < numLevels; ++lev) {
        const Box dataBox = a_data.getBoxes().physDomain().domainBox();
        const Box geoBox = vLevGeo[lev]->getDomain().domainBox();
        if (geoBox == dataBox) {
            if (a_data.getBoxes().compatible(vLevGeo[lev]->getBoxes())) break;
        }
    }
    if (lev == numLevels) {
        MayDay::Error("Could not find levelGeometry at data's index space.");
    }

    // I need to do a const_cast. I promise IO::writeHDF5 will not alter a_data
    // in any way. I am doing this because C++ has no good way of handling
    // a vector of pointers to const objects. Specifically, this doesn't work:
    //
    // void foo (Vector<const T*> a_data) {;}
    // int main () {
    //     Vector<T*> bar;
    //     foo(bar);  // ERROR!
    //
    //     // The reason this fails is because T* cannot be cast to const T*
    //     // when it is the template parameter. To understand why this type of
    //     // cast is disallowed in C++, take a look at
    //     // https://stackoverflow.com/questions/2102244/vector-and-const.
    //     // Too bad for us.
    // }
    LevelData<FArrayBox>* castDataPtr = const_cast<LevelData<FArrayBox>*>(&a_data);

    // Now that that's over with, we can package the level's data into a
    // data vector with only one non-NULL level.
    Vector<LevelData<FArrayBox>*> vData(numLevels, NULL);
    vData[lev] = castDataPtr;

    // Then, call the multi-level version of writeHDF5().
    IO::writeHDF5(a_filename, vData, a_levGeo, a_dt, a_time, a_compNames);
}


// -----------------------------------------------------------------------------
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
// -----------------------------------------------------------------------------
void
IO::writeHDF5 (std::string                         a_filename,
               const Vector<LevelData<FluxBox>* >& a_vectData,
               const LevelGeometry&                a_levGeo,
               const Real                          a_dt,
               const Real                          a_time,
               const Vector<std::string>&          a_compNames)
{
    // Gather geometric data
    const int                          finestLevel = a_vectData.size() - 1;
    const Vector<const LevelGeometry*> vLevGeo     = a_levGeo.getAMRLevGeos();
    const Vector<DisjointBoxLayout>    vGrids      = a_levGeo.getAMRGrids();
    const Vector<IntVect>              vRefRatio   = a_levGeo.getAMRRefRatios();
    const ProblemDomain& __attribute__((unused))   lev0Domain  = vLevGeo[0]->getDomain();
    const RealVect __attribute__((unused)) lev0dXi     = vLevGeo[0]->getDXi();

    // These used to be input parameters, but I think they just made things
    // complicated. I'm leaving them hard-coded to write the entire hierarchy.
    // Maybe someday I'll have a need to make these function inputs again.
    const int a_lmin = 0;
    const int a_lmax = finestLevel;
    CH_assert(0 <= a_lmin);
    CH_assert(a_lmin <= a_lmax);
    CH_assert(a_lmax <= finestLevel);

    // How many comps will we need to copy? How many ghosts are there?
    int numComps = -1;
    IntVect ghostVect = IntVect::Unit; // At least one for proper display of mapping.
    for (int lev = 0; lev <= finestLevel; ++lev) {
        if (a_vectData[lev] == NULL) continue;

        int thisNumComps = a_vectData[lev]->nComp();
        if (numComps == -1) numComps = thisNumComps;
        CH_assert(thisNumComps == numComps);

        IntVect thisGhostVect = a_vectData[lev]->ghostVect();
        D_TERM(ghostVect[0] = Max(ghostVect[0], thisGhostVect[0]);,
               ghostVect[1] = Max(ghostVect[1], thisGhostVect[1]);,
               ghostVect[2] = Max(ghostVect[2], thisGhostVect[2]);)
    }
    CH_assert(numComps > 0);

    // Did the user provide comp names? If not, generate default names.
    Vector<std::string> origCompNames(numComps);
    if (a_compNames.size() == 0) {
        int comp = 0;
        if (numComps == SpaceDim) {
            // It's probably a vector.
            D_TERM(origCompNames[comp++] = "x_component";,
                   origCompNames[comp++] = "y_component";,
                   origCompNames[comp++] = "z_component";)
        } else {
            // It's probably not a vector.
            while (comp < numComps) {
                char compName[40];
                sprintf(compName, "comp_%d", comp);
                origCompNames[comp++] = std::string(compName);
            }
        }
    } else {
        // The caller provided names.
        for (int comp = 0; comp < numComps; ++comp) {
            origCompNames[comp] = a_compNames[comp];
        }
    }

    // We will need 2*SpaceDim FArrayBox comps for each FluxBox comp.
    // left x, right x, left y, right y, ...
    const int newNumComps = numComps * 2 * SpaceDim;
    Vector<string> newCompNames(newNumComps);
    for (int newComp = 0; newComp < newNumComps; newComp += 2*SpaceDim) {
        int origComp = newComp / (2*SpaceDim);
        D_TERM(
        {
            ostringstream newCompName;
            newCompName << origCompNames[origComp] << " - left x comp " << origComp;
            newCompNames[newComp] = newCompName.str();
        }
        {
            ostringstream newCompName;
            newCompName << origCompNames[origComp] << " - right x comp " << origComp;
            newCompNames[newComp+1] = newCompName.str();
        },
        {
            ostringstream newCompName;
            newCompName << origCompNames[origComp] << " - left y comp " << origComp;
            newCompNames[newComp+2] = newCompName.str();
        }
        {
            ostringstream newCompName;
            newCompName << origCompNames[origComp] << " - right y comp " << origComp;
            newCompNames[newComp+3] = newCompName.str();
        },
        {
            ostringstream newCompName;
            newCompName << origCompNames[origComp] << " - left z comp " << origComp;
            newCompNames[newComp+4] = newCompName.str();
        }
        {
            ostringstream newCompName;
            newCompName << origCompNames[origComp] << " - right z comp " << origComp;
            newCompNames[newComp+5] = newCompName.str();
        })
    }

    Vector<LevelData<FArrayBox>*> ccData(finestLevel+1, nullptr);
    for (int lev = 0; lev <= finestLevel; ++lev) {
        const DisjointBoxLayout& grids = vGrids[lev];
        DataIterator dit = grids.dataIterator();

        ccData[lev] = new LevelData<FArrayBox>(grids, newNumComps, ghostVect);
        for (dit.reset(); dit.ok(); ++dit) {
            (*ccData[lev])[dit].setVal(quietNAN);
        }

        if (a_vectData[lev] == nullptr) continue;
        if (lev < a_lmin) continue;
        if (a_lmax < lev) continue;

        for (int newComp = 0; newComp < newNumComps; newComp += 2*SpaceDim) {
            int oldComp = newComp / (2*SpaceDim);

            for (int FCdir = 0; FCdir < SpaceDim; ++FCdir) {
                for (dit.reset(); dit.ok(); ++dit) {
                    FArrayBox& oldFAB = (FArrayBox&)((*a_vectData[lev])[dit][FCdir]);
                    FArrayBox& newFAB = (*ccData[lev])[dit];

                    newFAB.setVal(quietNAN, newComp + 2*FCdir);
                    newFAB.setVal(quietNAN, newComp + 2*FCdir + 1);

                    oldFAB.shiftHalf(FCdir,1);
                    newFAB.copy(oldFAB, oldComp, newComp + 2*FCdir, 1);

                    oldFAB.shiftHalf(FCdir,-2);
                    newFAB.copy(oldFAB, oldComp, newComp + 2*FCdir + 1, 1);

                    oldFAB.shiftHalf(FCdir,1);
                }
            }
        }
    }

    // Write to file.
    IO::writeHDF5(a_filename, ccData, a_levGeo, a_dt, a_time, newCompNames);

    // Free memory
    for (int lev = 0; lev <= finestLevel; ++lev) {
        delete ccData[lev];
    }
}


// -----------------------------------------------------------------------------
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
// -----------------------------------------------------------------------------
void
IO::writeHDF5 (std::string                a_filename,
               const LevelData<FluxBox>&  a_data,
               const LevelGeometry&       a_levGeo,
               const Real                 a_dt,
               const Real                 a_time,
               const Vector<std::string>& a_compNames)
{
    // We need to find where we are in the AMR hierarchy.
    const Vector<const LevelGeometry*> vLevGeo = a_levGeo.getAMRLevGeos();
    const int numLevels = a_levGeo.getNumLevels();

    int lev;
    for (lev = 0; lev < numLevels; ++lev) {
        const Box dataBox = a_data.getBoxes().physDomain().domainBox();
        const Box geoBox = vLevGeo[lev]->getDomain().domainBox();
        if (geoBox == dataBox) {
            if (a_data.getBoxes().compatible(vLevGeo[lev]->getBoxes())) break;
        }
    }
    if (lev == numLevels) {
        MayDay::Error("Could not find levelGeometry at data's index space.");
    }

    // I need to do a const_cast. I promise IO::writeHDF5 will not alter a_data
    // in any way. I am doing this because C++ has no good way of handling
    // a vector of pointers to const objects. Specifically, this doesn't work:
    //
    // void foo (Vector<const T*> a_data) {;}
    // int main () {
    //     Vector<T*> bar;
    //     foo(bar);  // ERROR!
    //
    //     // The reason this fails is because T* cannot be cast to const T*
    //     // when it is the template parameter. To understand why this type of
    //     // cast is disallowed in C++, take a look at
    //     // https://stackoverflow.com/questions/2102244/vector-and-const.
    //     // Too bad for us.
    // }
    LevelData<FluxBox>* castDataPtr = const_cast<LevelData<FluxBox>*>(&a_data);

    // Now that that's over with, we can package the level's data into a
    // data vector with only one non-NULL level.
    Vector<LevelData<FluxBox>*> vData(numLevels, nullptr);
    vData[lev] = castDataPtr;

    // Then, call the multi-level version of writeHDF5().
    IO::writeHDF5(a_filename, vData, a_levGeo, a_dt, a_time, a_compNames);
}


// // -----------------------------------------------------------------------------
// void
// IO::writeHDF5 (std::string                a_filename,
//                const LevelData<FluxBox>&  a_data,
//                const Real                 a_dt,
//                const Real                 a_time,
//                const Vector<std::string>& a_compNames)
// {
//     const DisjointBoxLayout& grids = a_data.getBoxes();
//     const ProblemDomain& domain = grids.physDomain();
//     const RealVect L(domain.size());

//     CartesianMap cartMap;
//     LevelGeometry levGeo(domain, L, nullptr, &cartMap);
//     levGeo.createMetricCache(grids);

//     IO::writeHDF5(a_filename, a_data, levGeo, a_dt, a_time, a_compNames);
// }
