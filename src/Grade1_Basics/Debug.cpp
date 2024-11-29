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
#include <iomanip>
#include "Debug.H"
#include "Format.H"
#include "FluxBox.H"
#include "AMRIO.H"
#include "CFRegion.H"
#include "SOMAR_Constants.H"
#include "SetValLevel.H"
#include "NeighborIterator.H"

#include <execinfo.h>
#include <cxxabi.h>
// #include <stdio.h>
// #include <stdlib.h>
// #include <unistd.h>


#ifdef CH_MULTIDIM
    using Chombo::pout;
#endif



// -----------------------------------------------------------------------------
// Private!
// Write a LevelData<FluxBox> to HDF5
// -----------------------------------------------------------------------------
void _debug_writeLevelHDF5 (const LevelData<FluxBox>& a_data,
                            const char*               a_filename,
                            Real                      a_time,
                            bool                      a_oneGhost)
{
#ifdef USE_CH_HDF5
    const DisjointBoxLayout& grids = a_data.getBoxes();
    const Box& domainBox = grids.physDomain().domainBox();
    const int nOldComps = a_data.nComp();
    const int nNewComps = nOldComps * 2 * SpaceDim;
    DataIterator dit = a_data.dataIterator();
    const IntVect& ghostVect = a_oneGhost? IntVect::Unit: a_data.ghostVect();

    Vector<LevelData<FArrayBox>*> vData(1);
    vData[0] = new LevelData<FArrayBox>(grids, nNewComps, ghostVect);

    Vector<DisjointBoxLayout> vGrids(1);
    vGrids[0] = grids;

    Vector<string> vNames(nNewComps);
    for (int newComp = 0; newComp < nNewComps; newComp += 2*SpaceDim) {
        int oldComp = newComp / (2*SpaceDim);
        D_TERM(
        {
            ostringstream newCompName;
            newCompName << "left x comp " << oldComp;
            vNames[newComp] = newCompName.str();
        }
        {
            ostringstream newCompName;
            newCompName << "right x comp " << oldComp;
            vNames[newComp+1] = newCompName.str();
        },
        {
            ostringstream newCompName;
            newCompName << "left y comp " << oldComp;
            vNames[newComp+2] = newCompName.str();
        }
        {
            ostringstream newCompName;
            newCompName << "right y comp " << oldComp;
            vNames[newComp+3] = newCompName.str();
        },
        {
            ostringstream newCompName;
            newCompName << "left z comp " << oldComp;
            vNames[newComp+4] = newCompName.str();
        }
        {
            ostringstream newCompName;
            newCompName << "right z comp " << oldComp;
            vNames[newComp+5] = newCompName.str();
        })

        for (int FCdir = 0; FCdir < SpaceDim; ++FCdir) {
            for (dit.reset(); dit.ok(); ++dit) {
                FArrayBox& oldFAB = (FArrayBox&)(a_data[dit][FCdir]);
                FArrayBox& newFAB = (*vData[0])[dit];

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

    Vector<int> refRatios(1,1);
    int numLevels = 1;
    Real dx = 1.0;
    Real dt = 1.0;
#ifdef CH_USE_HDF5
    WriteAMRHierarchyHDF5(string(a_filename),
                          vGrids, vData, vNames,
                          domainBox, dx, dt, a_time,
                          refRatios, numLevels);
#endif
    delete vData[0];
#else
    (void)a_data;
    (void)a_filename;
    (void)a_time;
    (void)a_oneGhost;
#endif
}


// -----------------------------------------------------------------------------
// Private!
// Write a LevelData<FArrayBox> to HDF5
// -----------------------------------------------------------------------------
void _debug_writeLevelHDF5 (const LevelData<FArrayBox>& a_data,
                            const char*                 a_filename,
                            Real                        a_time,
                            bool                        a_oneGhost)
{
#ifdef USE_CH_HDF5
    {
        // Figure out the centering of the data
        IntVect dataType = IntVect::Zero;
        DataIterator dit = a_data.dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            const Box thisBox = a_data[dit].box();
            if (thisBox.isEmpty()) continue;
            dataType = thisBox.type();
            break;
        }

#if CH_MPI
        // It's possible this rank had no grids. In that case, we must get the
        // dataType from the other ranks.
        for (int dir = 0; dir < SpaceDim; ++dir) {
            int localDataTypeDir = dataType[dir];
            int dataTypeDir;
            int ierr = MPI_Allreduce(&localDataTypeDir, &dataTypeDir, 1, MPI_INT, MPI_SUM, Chombo_MPI::comm);
            if (ierr != MPI_SUCCESS) {
                std::ostringstream errmsg;
                errmsg << "MPI_Allreduce failed. Error " << ierr << std::endl;
                MayDay::Error(errmsg.str().c_str());
            }
            dataType[dir] = ((dataTypeDir == 0)? 0: 1);
        }
#endif

        if (dataType == IntVect::Zero) {
            // Do nothing special
        } else if (dataType.sum() == 1) {
            // Call the FluxBox version
            LevelData<FluxBox> fluxData(a_data.getBoxes(), a_data.nComp(), a_data.ghostVect());
            D_TERM(int FCdir = 0;,
                   if (dataType[1] == 1) FCdir = 1;,
                   if (dataType[2] == 1) FCdir = 2;)

            for (dit.reset(); dit.ok(); ++dit) {
                fluxData[dit].setVal(0.0);
                fluxData[dit][FCdir].copy(a_data[dit]);
            }

            _debug_writeLevelHDF5(fluxData, a_filename, a_time, a_oneGhost);
            return;

        } else {
            // Throw an error
            pout() << "dataType = " << dataType << endl;
            MayDay::Error("_debug_writeLevelHDF5 only works with CC or FC data and needs comm code for dataType");
        }
    }

    Vector<LevelData<FArrayBox>*> vData(1);
    if (a_oneGhost) {
        const DisjointBoxLayout& grids = a_data.getBoxes();
        vData[0] = new LevelData<FArrayBox>(grids, a_data.nComp(), IntVect::Unit);

        DataIterator dit = grids.dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            (*vData[0])[dit].copy(a_data[dit]);
        }
        // Copier oneGhostCopier(grids, grids, grids.physDomain(), IntVect::Unit, false);
        // a_data.copyTo(a_data.interval(), *vData[0], vData[0]->interval(), oneGhostCopier);
    } else {
        vData[0] = (LevelData<FArrayBox>*)(&a_data);
    }

    Vector<DisjointBoxLayout> vGrids(1);
    vGrids[0] = a_data.getBoxes();
    const Box& domainBox = vGrids[0].physDomain().domainBox();

    Vector<string> vNames(vData[0]->nComp());
    for (int comp = 0; comp < vData[0]->nComp(); ++comp) {
        ostringstream compName;
        compName << "comp " << comp;
        vNames[comp] = compName.str();
    }

    Vector<int> refRatios(1,1);
    int numLevels = 1;
    Real dx = 1.0;
    Real dt = 1.0;
#ifdef CH_USE_HDF5
    WriteAMRHierarchyHDF5(string(a_filename),
                          vGrids, vData, vNames,
                          domainBox, dx, dt, a_time,
                          refRatios, numLevels);
#endif
    if (a_oneGhost) {
        delete vData[0];
    }

#else
    (void)a_data;
    (void)a_filename;
    (void)a_time;
    (void)a_oneGhost;
#endif
}


// -----------------------------------------------------------------------------
// Puts a checkpoint into pout.*
// -----------------------------------------------------------------------------
void _WriteCheckpoint(std::string   /*a_func*/,
                      std::string   a_file,
                      unsigned int  a_line,
                      std::string   a_mark) {
    using namespace std;
    // string::size_type   start, end, center;

    // center = a_func.find_first_of("::");
    // start = a_func.rfind(" ", center-1) + 1;
    // end   = a_func.find_first_of("(", center+1) - 1;
    // a_func = a_func.substr(start, end - start + 1) + "()";

    // start = a_file.find_last_of("/") + 1;
    // if(start != string::npos && start < a_file.length())
    //     a_file = a_file.substr(start, a_file.length() - start);

    // if(a_mark.length() > 0) pout() << setfill('.');
    // pout() << "CHECK:\t" << setiosflags(ios::left) << setw(80) << a_func << "\t" << setw(30) << a_file;
    // if(a_mark.length() > 0) pout() << setfill(' ');
    // pout() << "[" << setw(4) << a_line << "] " << a_mark << endl;

    pout() << "CHECK: " << a_file << " [" << a_line << "] " << a_mark
           << endl;
}


// -----------------------------------------------------------------------------
// Used to make TODO() output readable.
// -----------------------------------------------------------------------------
static std::string stripFuncName(const std::string& a_func) {
    std::string::size_type  start, end, center;
    center = a_func.find_last_of("::");
    start = a_func.rfind(" ", center-1) + 1;
    end   = a_func.find_first_of("(", center+1) - 1;
    return a_func.substr(start, end - start + 1) + "()";
}


// -----------------------------------------------------------------------------
// Used to make TODO() output readable.
// -----------------------------------------------------------------------------
static std::string stripFileName(const std::string& a_file) {
    std::string::size_type start = a_file.find_last_of("/") + 1;
    if(start != std::string::npos && start < a_file.length())
        return a_file.substr(start, a_file.length() - start);
    return a_file;
}


// -----------------------------------------------------------------------------
// Write a note to the terminal.
// -----------------------------------------------------------------------------
void
_alwaysNote(const char*       a_filename,
            const char*       a_funcname,
            const int         a_linenumber,
            const std::string a_msg,
            const std::string a_premsg)
{
    if (procID() == 0) {
        ostringstream ss;
        if (a_premsg.length() == 0) {
            ss << Format::brown << "NOTE: ";
        } else {
            ss << a_premsg;
        }

        ss << Format::none
           << stripFuncName(a_funcname) << " in "
           << stripFileName(a_filename) << ":" << a_linenumber;

        if (a_msg.length() > 0) {
            ss << Format::hiwhite << ": " << a_msg << Format::none;
        }
        std::cout << ss.str() << std::endl;
    }
}


// -----------------------------------------------------------------------------
// Send a "needs testing" message to stdout
// -----------------------------------------------------------------------------
#ifndef NDEBUG
void _test(const char* a_filename, const char* a_funcname, const int a_linenumber) {
        if (procID() == 0) {
            std::cout << Format::brown << "TODO:"
                      << Format::none << " Test "
                      << stripFuncName(a_funcname)
                      << " in " << stripFileName(a_filename)
                      << " [" << a_linenumber << "]" << std::endl;
        }
}
#endif // End debug code


// -----------------------------------------------------------------------------
// Throws an error if a NAN is found in the testBox.
// -----------------------------------------------------------------------------
void
_checkForNAN(const FArrayBox&  a_data,
             const Box&        a_testBox,
             const std::string a_file,
             const int         a_line,
             naninfo_type*     a_info)
{
    CH_assert(a_data.box().sameType(a_testBox));

    if (a_info != NULL) {
        a_info->problemFound = false;
    }
    const int proc = procID();

    Box region = a_testBox;
    region &= a_data.box();

    for (int c = 0; c < a_data.nComp(); ++c) {
        BoxIterator bit(region);
        for (bit.reset(); bit.ok(); ++bit) {
            const IntVect& i = bit();
            Real v = a_data(i,c);

            if (v != v || v <= -1.0e10  || 1.0e10 <= v) {
                ostringstream str;
                    str << "\ncheckForNAN at " << a_file << ":" << a_line
                        << " found an error.\n"
                        << "FArrayBox(" << i << ", " << c
                        << ") = " << v << " on proc " << proc;
                pout() << str.str().c_str() << std::endl;

                if (a_info == NULL) {
                    MayDay::Error(str.str().c_str());
                } else {
                    a_info->problemFound = true;
                    a_info->pos = i;
                    a_info->comp = c;
                    a_info->box = region;
                    a_info->val = v;
                    a_info->msg = str.str();
                    return;
                }
            }
        }
    }
}


// -----------------------------------------------------------------------------
// Throws an error if a NAN is found in the testBox.
// -----------------------------------------------------------------------------
void
_checkForNAN(const FluxBox&    a_data,
             const Box&        a_testBox,
             const std::string a_file,
             const int         a_line,
             naninfo_type*     a_info)
{
    CH_assert(a_testBox.type() == IntVect::Zero);

    if (a_info != NULL) {
        a_info->problemFound = false;
    }
    const int proc = procID();

    Box CCregion = a_testBox;
    CCregion &= a_data.box();

    for (int d = 0; d < SpaceDim; ++d) {
        Box region = CCregion;
        region.surroundingNodes(d);
        region &= a_data[d].box();

        for (int c = 0; c < a_data.nComp(); ++c) {
            BoxIterator bit(region);
            for (bit.reset(); bit.ok(); ++bit) {
                const IntVect& i = bit();
                Real v = a_data[d](i,c);

                if (v != v || v <= -1.0e10  || 1.0e10 <= v) {
                    ostringstream str;
                    str << "\ncheckForNAN at " << a_file << ":" << a_line
                        << " found an error.\n"
                        << "FluxBox[" << d << "](" << i << ", " << c
                        << ") = " << v << " on proc " << proc;
                    pout() << str.str().c_str() << std::endl;

                    if (a_info == NULL) {
                        MayDay::Error(str.str().c_str());
                    } else {
                        a_info->problemFound = true;
                        a_info->pos = i;
                        a_info->comp = c;
                        a_info->box = region;
                        a_info->val = v;
                        a_info->msg = str.str();
                        return;
                    }
                }
            }
        }
    }
}


// -----------------------------------------------------------------------------
// Throws an error if a NAN is found in the valid regions.
// -----------------------------------------------------------------------------
void
_checkForValidNAN(const LevelData<FArrayBox>& a_data,
                  const std::string           a_file,
                  const int                   a_line)
{
    const DisjointBoxLayout& grids = a_data.getBoxes();
    DataIterator dit = grids.dataIterator();

    naninfo_type info;
    info.problemFound = false;

    for (dit.reset(); dit.ok(); ++dit) {
        const FArrayBox& d = a_data[dit];
        Box valid = grids[dit];
        valid &= d.box();

        _checkForNAN(d, valid, a_file, a_line, &info);
        if (info.problemFound == true) break;
    }

    Vector<int> vpf(numProc(), 0);
    int pf = info.problemFound? 1: 0;
    const int srcProc = uniqueProc(SerialTask::compute);
    gather(vpf, pf, srcProc);
    if (procID() == srcProc) {
        for (unsigned int idx = 0; idx < vpf.size(); ++idx) {
            if (vpf[idx] == 1) {
                pf = 1;
                break;
            }
        }
    }
    broadcast(pf, srcProc);

    if (pf == 1) {
        _debug_writeLevelHDF5(a_data, "NANData.hdf5", 0.0, false);
        ostringstream msg;
        msg << Format::hired << "BUG: " << Format::none
            << "NAN found at " << a_file.c_str() << ":" << a_line
            << ". Data written to NANData.hdf5."
            << info.msg.c_str();
        MayDay::Error(msg.str().c_str());
    }
}


// -----------------------------------------------------------------------------
// Throws an error if a NAN is found in the valid regions.
// -----------------------------------------------------------------------------
void
_checkForValidNAN(const LevelData<FluxBox>& a_data,
                  const std::string         a_file,
                  const int                 a_line)
{
    const DisjointBoxLayout& grids = a_data.getBoxes();
    DataIterator dit = grids.dataIterator();

    naninfo_type info;
    info.problemFound = false;

    for (dit.reset(); dit.ok(); ++dit) {
        const FluxBox& d = a_data[dit];
        Box valid = grids[dit];
        valid &= d.box();

        _checkForNAN(d, valid, a_file, a_line, &info);
        if (info.problemFound == true) break;
    }

    Vector<int> vpf(numProc(), 0);
    int pf = info.problemFound? 1: 0;
    const int srcProc = uniqueProc(SerialTask::compute);
    gather(vpf, pf, srcProc);
    if (procID() == srcProc) {
        for (unsigned int idx = 0; idx < vpf.size(); ++idx) {
            if (vpf[idx] == 1) {
                pf = 1;
                break;
            }
        }
    }
    broadcast(pf, srcProc);

    if (pf == 1) {
        _debug_writeLevelHDF5(a_data, "NANData.hdf5", 0.0, false);
        ostringstream msg;
        msg << Format::hired << "BUG: " << Format::none
            << "NAN found at " << a_file.c_str() << ":" << a_line
            << ". Data written to NANData.hdf5."
            << info.msg.c_str();
        MayDay::Error(msg.str().c_str());
    }
}


// -----------------------------------------------------------------------------
// Initializes data holders to NAN
// -----------------------------------------------------------------------------
#if defined(USE_INIT_FUNCTIONS)
    // -------------------------------------------------------------------------
    // Set a_levels[a_min:a_max] to NAN.
    // -------------------------------------------------------------------------
    void
    debugInitLevels(Vector<LevelData<FArrayBox>*>& a_levels, int a_min, int a_max)
    {
        setValLevels(a_levels, a_min, a_max, quietNAN);
    }


    // -------------------------------------------------------------------------
    // Set a_level to NAN.
    // -------------------------------------------------------------------------
    void debugInitLevel (LevelData<FArrayBox>& a_level)
    {
        setValLevel(a_level, quietNAN);
    }


    // -------------------------------------------------------------------------
    // Set a_level to NAN.
    // -------------------------------------------------------------------------
    void
    debugInitLevel(BoxLayoutData<FArrayBox>& a_level)
    {
        setValLevel(a_level, quietNAN);
    }


    // -------------------------------------------------------------------------
    // Set a_level to NAN.
    // -------------------------------------------------------------------------
    void
    debugInitLevel(LevelData<FluxBox>& a_level)
    {
        setValLevel(a_level, quietNAN);
    }


    // -------------------------------------------------------------------------
    // Set a_level to NAN.
    // -------------------------------------------------------------------------
    void
    debugInitLevel(LevelData<NodeFArrayBox>& a_level)
    {
        setValLevel(a_level, quietNAN);
    }


    // -------------------------------------------------------------------------
    // Set FArrayBox to NAN.
    // -------------------------------------------------------------------------
    void
    debugInit (FArrayBox& a_fab)
    {
        a_fab.setVal(quietNAN);
    }


    // -------------------------------------------------------------------------
    // Set FluxBox to NAN.
    // -------------------------------------------------------------------------
    void
    debugInit(FluxBox& a_flub)
    {
        a_flub.setVal(quietNAN);
    }


    // -------------------------------------------------------------------------
    // Set ghosts of a_level to NAN.
    // -------------------------------------------------------------------------
    void
    debugInitLevelGhosts(LevelData<FArrayBox>& a_level)
    {
        const DisjointBoxLayout& grids = a_level.getBoxes();
        DataIterator dit = a_level.dataIterator();

        for (dit.reset(); dit.ok(); ++dit) {
            CH_assert(a_level[dit].box().type() == IntVect::Zero);
            debugInitGhosts(a_level[dit], grids[dit]);
        }
    }


    // -------------------------------------------------------------------------
    // Set ghosts of a_level to NAN.
    // -------------------------------------------------------------------------
    void
    debugInitLevelGhosts(LevelData<FluxBox>& a_level)
    {
        const DisjointBoxLayout& grids = a_level.getBoxes();
        DataIterator dit = a_level.dataIterator();

        for (dit.reset(); dit.ok(); ++dit) {
            debugInitGhosts(a_level[dit], grids[dit]);
        }
    }


    // -------------------------------------------------------------------------
    // Set cells of a_fab that lie outside of a_ccValid to NAN.
    // -------------------------------------------------------------------------
    void debugInitGhosts(FArrayBox& a_fab, Box a_valid)
    {
        CH_assert(a_fab.box().type() == a_valid.type());

        a_valid &= a_fab.box();

        if (a_valid.isEmpty()) {
            a_fab.setVal(quietNAN);
        } else {
            FArrayBox tmpFAB(a_valid, a_fab.nComp());
            tmpFAB.copy(a_fab);
            a_fab.setVal(quietNAN);
            a_fab.copy(tmpFAB);
        }
    }


    // -------------------------------------------------------------------------
    // Set faces of a_fab that lie outside of a_valid to NAN.
    // -------------------------------------------------------------------------
    void
    debugInitGhosts(FluxBox& a_flub, const Box& a_ccValid)
    {
        for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
            const Box fcValid = surroundingNodes(a_ccValid, fcDir);
            debugInitGhosts(a_flub[fcDir], fcValid);
        }
    }


    // -------------------------------------------------------------------------
    // Set physical and CFI boundary faces of a_level to NAN.
    // -------------------------------------------------------------------------
    void
    debugInitLevelBoundaryFaces(LevelData<FluxBox>& a_level)
    {
        const DisjointBoxLayout& grids  = a_level.getBoxes();
        const ProblemDomain&     domain = grids.physDomain();
        DataIterator             dit    = a_level.dataIterator();
        CFRegion                 cfRegion(grids, domain);

        for (dit.reset(); dit.ok(); ++dit) {
            for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
                FArrayBox&    levelFAB = a_level[dit][fcDir];
                const IntVect e        = BASISV(fcDir);

                for (SideIterator sit; sit.ok(); ++sit) {
                    // Are we at a physical boundary?
                    if (!domain.isPeriodic(fcDir)) {
                        Box physBdry =
                            bdryBox(domain.domainBox(), fcDir, sit());
                        physBdry &= levelFAB.box();
                        if (!physBdry.isEmpty()) {
                            levelFAB.setVal(
                                quietNAN, physBdry, 0, levelFAB.nComp());
                        }
                    }

                    // Are we at the CFI?
                    const CFIVS& cfivs =
                        (sit() == Side::Lo ? cfRegion.loCFIVS(dit(), fcDir)
                                           : cfRegion.hiCFIVS(dit(), fcDir));
                    if (!cfivs.isEmpty()) {
                        IntVectSet ivs(cfivs.getIVS());
                        if (sit() == Side::Lo) {
                            ivs.shift(e);
                        }
                        IVSIterator ivsit(ivs);

                        for (int comp = 0; comp < levelFAB.nComp(); ++comp) {
                            for (ivsit.reset(); ivsit.ok(); ++ivsit) {
                                const IntVect& iv  = ivsit();
                                levelFAB(iv, comp) = quietNAN;
                            }
                        }
                    }
                }
            }
        }
    }
#endif
