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
#include "AMRNSLevel.H"
#include "Convert.H"
#include "Debug.H"
#include "FABAlgebra.H"
#include "HeaderData.H"
#include "Subspace.H"
#include "Masks.H"
#include "ProblemContext.H"
#ifdef CH_USE_PYTHON
#include "PyGlue.H"
#endif


// #define WRITE_METRIC_TO_HDF5
// #define WRITE_GRADP_TO_HDF5
// #define WRITE_Sij_TO_HDF5
// #define WRITE_VORTICITY_TO_HDF5

// -----------------------------------------------------------------------------
// write checkpoint header
// -----------------------------------------------------------------------------
void
AMRNSLevel::writeCheckpointHeader(HeaderData&        header,
                                  const std::string& a_filename) const
{
    BEGIN_FLOWCHART();

    const ProblemContext* ctx = ProblemContext::getInstance();
    char comp_str[80];

    // Create the header. This will only store metadata
    // about the number of fields and thier names.

    // write some Chombo boilerplate (do we need it?)
    {
        HeaderData temp;
        temp.m_int["SpaceDim"]  = SpaceDim;
        temp.m_real["testReal"] = (Real)0.0;
        temp.writeToFile(a_filename, std::string("Chombo_global"));
    }

    // AMR parameters
    {
        header.m_real["time.fixedDt"]     = ctx->time.fixedDt;
        header.m_real["time.dtMult"]      = ctx->time.dtMult;
        header.m_real["time.maxDt"]       = ctx->time.maxDt;
        header.m_int["amr.useSubcycling"] = ctx->amr.useSubcycling ? 1 : 0;
        header.m_int["amr.maxLevel"]      = ctx->amr.maxLevel;
        header.m_int["base.blockFactor"]  = ctx->base.blockFactor;
        header.m_intvect["base.nx"]       = ctx->base.nx;
        header.m_intvect["base.nxOffset"] = ctx->base.nxOffset;
        header.m_int["base.isPeriodic_0"] = (ctx->base.isPeriodic[0] ? 1 : 0);
        header.m_int["base.isPeriodic_1"] = (ctx->base.isPeriodic[1] ? 1 : 0);
        if (SpaceDim > 2) {
            header.m_int["base.isPeriodic_2"] =
                (ctx->base.isPeriodic[SpaceDim - 1] ? 1 : 0);
        }
    }

    // RHS parameters
    {
        header.m_realvect["base.L"]  = ctx->base.L;
        header.m_real["rhs.nu"]     = ctx->rhs.nu;
        header.m_real["rhs.TKappa"] = ctx->rhs.TKappa;
        header.m_real["rhs.SKappa"] = ctx->rhs.SKappa;

        header.m_int["rhs.scalarsKappaSize"] = this->numScalars();
        for (int idx = 0; idx < this->numScalars(); ++idx) {
            sprintf(comp_str, "rhs.scalarsKappa_%d", idx);
            header.m_real[comp_str] = ctx->rhs.getScalarsKappa(idx);
        }

        header.m_realvect["rhs.coriolisF"] = ctx->rhs.coriolisF;
    }

    // LES stuff
    {
        const Vector<int>& eddyViscMethod = ctx->rhs.eddyViscMethod;
        header.m_int["rhs.eddyViscMethodSize"] = eddyViscMethod.size();
        for (unsigned int idx = 0; idx < eddyViscMethod.size(); ++idx) {
            sprintf(comp_str, "rhs.eddyViscMethod_%d", idx);
            header.m_int[comp_str] = eddyViscMethod[idx];
        }

        header.m_real["rhs.eddyPrandtlT"] = ctx->rhs.eddyPrandtlT;
        header.m_real["rhs.eddyPrandtlS"] = ctx->rhs.eddyPrandtlS;

        header.m_int["rhs.eddyPrandtlScalarsSize"] = this->numScalars();
        for (int idx = 0; idx < this->numScalars(); ++idx) {
            sprintf(comp_str, "rhs.eddyPrandtlScalars_%d", idx);
            header.m_real[comp_str] = ctx->rhs.getEddyPrandtlScalars(idx);
        }
    }

    // Scalar names
    {
        const int nScal = (int)this->numScalars();
        header.m_int["numScalars"] = nScal;
        for (int comp = 0; comp < nScal; ++comp) {
            sprintf(comp_str, "scalarComponent_%d", comp);
            header.m_string[comp_str] = this->getScalarName(comp);
        }
    }

    // State components.
    header.m_int["numQComps"] = m_statePtr->numQComps;

    // Write the metadata to HDF5 and pout.*
    header.writeToFile(a_filename, std::string("/"));
    if (s_verbosity >= 6) {
        pout() << header << endl;
    }
}


// -----------------------------------------------------------------------------
// write checkpoint data for this level
// -----------------------------------------------------------------------------
void
AMRNSLevel::writeCheckpointLevel(const std::string& a_fileName, int /*level*/) const
{
    BEGIN_FLOWCHART();

    // Set group for this level.
    char level_str[20];
    sprintf(level_str, "%d", m_level);
    const std::string label = std::string("level_") + level_str;

    // Create the header for this group.
    // This will store the level's metadata.
    HeaderData header;

    // Collect all metadata that will be needed at restart.
    // header.m_int     ["step_number"] = s_step_number;
    header.m_intvect["ref_ratio"]     = this->getFineRefRatio();
    header.m_realvect["vec_dx"]       = m_levGeoPtr->getDXi();
    header.m_real["dt"]               = m_dt;
    header.m_real["time"]             = m_time;
    header.m_int["isEmpty"]           = (this->isEmpty() ? 1 : 0);
    header.m_box["domainBox"]         = this->getDomainBox();
    header.m_int["finestExtantLevel"] = this->finestNSPtr()->m_level;

    // Write the metadata to file and pout.*
    header.writeToFile(a_fileName, label);
    if (s_verbosity >= 6) {
        pout() << header << endl;
    }

    // If this level has valid data, we need to write it to HDF5.
    if (!this->isEmpty()) {
#ifdef CH_USE_PYTHON
        std::string name = "velData";
        Py::PythonFunction("IO",
                              "WriteCheckPointLevelDataFluxBox",
                              a_fileName,
                              label,
                              name,
                              *m_velPtr,
                              IntVect::Unit);
        name = "pData";
        Py::PythonFunction("IO",
                              "WriteCheckPointLevelDataFAB",
                              a_fileName,
                              label,
                              name,
                              *m_pPtr,
                              IntVect::Unit);
        name = "qData";
        Py::PythonFunction("IO",
                              "WriteCheckPointLevelDataFAB",
                              a_fileName,
                              label,
                              name,
                              *m_qPtr,
                              IntVect::Unit);
#endif
    }
}

// -----------------------------------------------------------------------------
// read checkpoint header
// -----------------------------------------------------------------------------
void
AMRNSLevel::readCheckpointHeader(const std::string& a_fileName)
{
    BEGIN_FLOWCHART();

    const ProblemContext* ctx = ProblemContext::getInstance();
    char comp_str[80];

    // Get the checkpoint's metadata
    HeaderData header(a_fileName, s_verbosity);

    // time.fixedDt
    if (header.m_real.find("time.fixedDt") == header.m_real.end()) {
        MayDay::Warning("Checkfile does not have time.fixedDt");
    } else {
        Real fixedDt = header.m_real["time.fixedDt"];
        if (abs(fixedDt - ctx->time.fixedDt) > smallReal) {
            ostringstream msg;
            msg << "time.fixedDt changed from " << fixedDt << " to "
                << ctx->time.fixedDt;
            MayDay::Warning(msg.str().c_str());
        }
    }

    // time.dtMult
    if (header.m_real.find("time.dtMult") == header.m_real.end()) {
        MayDay::Warning("Checkfile does not have time.dtMult");
    } else {
        Real dtMult = header.m_real["time.dtMult"];
        if (abs(dtMult - ctx->time.dtMult) > smallReal) {
            ostringstream msg;
            msg << "time.dtMult changed from " << dtMult << " to "
                << ctx->time.dtMult;
            MayDay::Warning(msg.str().c_str());
        }
    }

    // time.maxDt
    if (header.m_real.find("time.maxDt") == header.m_real.end()) {
        MayDay::Warning("Checkfile does not have time.maxDt");
    } else {
        Real maxDt = header.m_real["time.maxDt"];
        if (abs(maxDt - ctx->time.maxDt) > smallReal) {
            ostringstream msg;
            msg << "time.maxDt changed from " << maxDt << " to "
                << ctx->time.maxDt;
            MayDay::Warning(msg.str().c_str());
        }
    }

    // amr.useSubcycling
    if (header.m_int.find("amr.useSubcycling") == header.m_int.end()) {
        MayDay::Warning("Checkfile does not have amr.useSubcycling");
    } else {
        int useSubcycling = header.m_int["amr.useSubcycling"];
        if (useSubcycling != ctx->amr.useSubcycling) {
            ostringstream msg;
            msg << "amr.useSubcycling changed from " << useSubcycling << " to "
                << ctx->amr.useSubcycling;
            MayDay::Warning(msg.str().c_str());
        }
    }

    // amr.maxLevel
    if (header.m_int.find("amr.maxLevel") == header.m_int.end()) {
        MayDay::Warning("Checkfile does not have amr.maxLevel");
    } else {
        int maxLevel = header.m_int["amr.maxLevel"];
        if (maxLevel != ctx->amr.maxLevel) {
            ostringstream msg;
            msg << "amr.maxLevel changed from " << maxLevel << " to "
                << ctx->amr.maxLevel;
            MayDay::Warning(msg.str().c_str());
        }
    }

    // base.blockFactor
    if (header.m_int.find("base.blockFactor") == header.m_int.end()) {
        MayDay::Warning("Checkfile does not have base.blockFactor");
    } else {
        int blockFactor = header.m_int["base.blockFactor"];
        if (blockFactor != ctx->base.blockFactor) {
            ostringstream msg;
            msg << "base.blockFactor changed from " << blockFactor << " to "
                << ctx->base.blockFactor
                << ". If you get an obscure runtime error"
                << " this may be why.";
            MayDay::Warning(msg.str().c_str());
        }
    }

    // base.nx
    if (header.m_intvect.find("base.nx") == header.m_intvect.end()) {
        MayDay::Error("Checkfile does not have base.nx");
    } else {
        IntVect base_nx = header.m_intvect["base.nx"];
        if (base_nx != ctx->base.nx) {
            ostringstream msg;
            msg << "base.nx changed from " << base_nx << " to " << ctx->base.nx;
            MayDay::Error(msg.str().c_str());
        }
    }

    // base.nxOffset
    if (header.m_intvect.find("base.nxOffset") == header.m_intvect.end()) {
        MayDay::Error("Checkfile does not have base.nxOffset");
    } else {
        IntVect base_nxOffset = header.m_intvect["base.nxOffset"];
        if (base_nxOffset != ctx->base.nxOffset) {
            ostringstream msg;
            msg << "base.nxOffset changed from " << base_nxOffset << " to "
                << ctx->base.nxOffset;
            MayDay::Error(msg.str().c_str());
        }
    }

    // base.isPeriodic_0
    if (header.m_int.find("base.isPeriodic_0") == header.m_int.end()) {
        MayDay::Warning("Checkfile does not have base.isPeriodic_0");
    } else {
        int isPeriodic_0 = header.m_int["base.isPeriodic_0"];
        if (isPeriodic_0 != ctx->base.isPeriodic[0]) {
            ostringstream msg;
            msg << "base.isPeriodic[0] changed from " << isPeriodic_0 << " to "
                << ctx->base.isPeriodic[0];
            MayDay::Warning(msg.str().c_str());
        }
    }

    // base.isPeriodic_1
    if (header.m_int.find("base.isPeriodic_1") == header.m_int.end()) {
        MayDay::Warning("Checkfile does not have base.isPeriodic_1");
    } else {
        int isPeriodic_1 = header.m_int["base.isPeriodic_1"];
        if (isPeriodic_1 != ctx->base.isPeriodic[1]) {
            ostringstream msg;
            msg << "base.isPeriodic[1] changed from " << isPeriodic_1 << " to "
                << ctx->base.isPeriodic[1];
            MayDay::Warning(msg.str().c_str());
        }
    }

    // base.isPeriodic_2
    if (SpaceDim > 2) {
        if (header.m_int.find("base.isPeriodic_2") == header.m_int.end()) {
            MayDay::Warning("Checkfile does not have isPeriodic_2");
        } else {
            int isPeriodic_2 = header.m_int["isPeriodic_2"];
            if (isPeriodic_2 != ctx->base.isPeriodic[SpaceDim - 1]) {
                ostringstream msg;
                msg << "base.isPeriodic[2] changed from " << isPeriodic_2
                    << " to " << ctx->base.isPeriodic[SpaceDim - 1];
                MayDay::Warning(msg.str().c_str());
            }
        }
    }

    // base.L
    if (header.m_realvect.find("base.L") == header.m_realvect.end()) {
        MayDay::Error("Checkfile does not have base.L");
    } else {
        RealVect L = header.m_realvect["base.L"];
        if (L - ctx->base.L >=
            RealVect(D_DECL(smallReal, smallReal, smallReal))) {
            ostringstream msg;
            msg << "base.L changed from " << L << " to " << ctx->base.L;
            MayDay::Error(msg.str().c_str());
        }
    }

    // rhs.nu
    if (header.m_real.find("rhs.nu") == header.m_real.end()) {
        MayDay::Warning("Checkfile does not have rhs.nu");
    } else {
        Real nu = header.m_real["rhs.nu"];
        if (abs(nu - ctx->rhs.nu) > smallReal) {
            ostringstream msg;
            msg << "rhs.nu changed from " << nu << " to " << ctx->rhs.nu;
            MayDay::Warning(msg.str().c_str());
        }
    }

    // rhs.TKappa
    if (header.m_real.find("rhs.TKappa") == header.m_real.end()) {
        MayDay::Warning("Checkfile does not have rhs.TKappa");
    } else {
        Real TKappa = header.m_real["rhs.TKappa"];
        if (abs(TKappa - ctx->rhs.TKappa) > smallReal) {
            ostringstream msg;
            msg << "rhs.TKappa changed from " << TKappa << " to "
                << ctx->rhs.TKappa;
            MayDay::Warning(msg.str().c_str());
        }
    }

    // rhs.SKappa
    if (header.m_real.find("rhs.SKappa") == header.m_real.end()) {
        MayDay::Warning("Checkfile does not have rhs.SKappa");
    } else {
        Real SKappa = header.m_real["rhs.SKappa"];
        if (abs(SKappa - ctx->rhs.SKappa) > smallReal) {
            ostringstream msg;
            msg << "rhs.SKappa changed from " << SKappa << " to "
                << ctx->rhs.SKappa;
            MayDay::Warning(msg.str().c_str());
        }
    }

    // rhs.scalarsKappaSize
    int scalarsKappaSize = 0;
    if (header.m_int.find("rhs.scalarsKappaSize") == header.m_int.end()) {
        MayDay::Warning("Checkfile does not have rhs.scalarsKappaSize");
    } else {
        scalarsKappaSize = header.m_int["rhs.scalarsKappaSize"];
    }
    if (scalarsKappaSize != this->numScalars()) {
        MAYDAYERROR("There are " << this->numScalars()
                                 << " scalars, but rhs.scalarsKappaSize is "
                                 << scalarsKappaSize << ".");
    }

    // rhs.scalarsKappa
    for (int idx = 0; idx < scalarsKappaSize; ++idx) {
        sprintf(comp_str, "rhs.scalarsKappa_%d", idx);

        if (header.m_real.find(comp_str) == header.m_real.end()) {
            ostringstream msg;
            msg << "Checkfile does not have " << comp_str;
            MayDay::Warning(msg.str().c_str());
        } else {
            Real scalarsKappaComp = header.m_real[comp_str];
            if (RealCmp::neq(scalarsKappaComp, ctx->rhs.getScalarsKappa(idx))) {
                ostringstream msg;
                msg << "rhs.scalarsKappa[" << idx << "] changed from "
                    << scalarsKappaComp << " to " << ctx->rhs.getScalarsKappa(idx);
                MayDay::Warning(msg.str().c_str());
            }
        }
    }

    // rhs.coriolisF
    if (header.m_realvect.find("rhs.coriolisF") == header.m_realvect.end()) {
        MayDay::Warning("Checkfile does not have rhs.coriolisF");
    } else {
        const RealVect coriolisF = header.m_realvect["rhs.coriolisF"];

        if (D_TERM( RealCmp::neq(coriolisF[0], ctx->rhs.coriolisF[0]),
                 || RealCmp::neq(coriolisF[1], ctx->rhs.coriolisF[1]),
                 || RealCmp::neq(coriolisF[2], ctx->rhs.coriolisF[2]) )) {
            ostringstream msg;
            msg << "rhs.coriolisF changed from " << coriolisF << " to "
                << ctx->rhs.coriolisF;
            MayDay::Warning(msg.str().c_str());
        }
    }

    // rhs.eddyViscMethodSize
    unsigned int eddyViscMethodSize = 0;
    if (header.m_int.find("rhs.eddyViscMethodSize") == header.m_int.end()) {
        MayDay::Warning("Checkfile does not have rhs.eddyViscMethodSize");
    } else {
        eddyViscMethodSize = header.m_int["rhs.eddyViscMethodSize"];
        if (eddyViscMethodSize != ctx->rhs.eddyViscMethod.size()) {
            ostringstream msg;
            msg << "rhs.eddyViscMethodSize changed from " << eddyViscMethodSize
                << " to " << ctx->rhs.eddyViscMethod.size();
            MayDay::Warning(msg.str().c_str());
        }
    }

    // rhs.eddyViscMethod
    CH_assert(ctx->rhs.eddyViscMethod.size() == eddyViscMethodSize);
    for (unsigned int idx = 0; idx < eddyViscMethodSize; ++idx) {
        sprintf(comp_str, "rhs.eddyViscMethod_%d", idx);

        if (header.m_int.find(comp_str) == header.m_int.end()) {
            ostringstream msg;
            msg << "Checkfile does not have " << comp_str;
            MayDay::Warning(msg.str().c_str());
        } else {
            Real eddyViscMethodComp = header.m_int[comp_str];
            if (eddyViscMethodComp != ctx->rhs.eddyViscMethod[idx]) {
                ostringstream msg;
                msg << "rhs.eddyViscMethod[" << idx << "] changed from "
                    << eddyViscMethodComp << " to "
                    << ctx->rhs.eddyViscMethod[idx];
                MayDay::Warning(msg.str().c_str());
            }
        }
    }

    TODONOTE("Read eddyPrandtl* from checkpoint?");

    // numScalars
    if (header.m_int.find("numScalars") == header.m_int.end()) {
        MayDay::Error("Checkfile does not have numScalars");
    } else {
        int nScal = header.m_int["numScalars"];
        if (nScal != (int)this->numScalars()) {
            ostringstream msg;
            msg << "numScalars changed from " << nScal << " to "
                << this->numScalars() << ".";
            MayDay::Error(msg.str().c_str());
        }
    }

    // scalarComponent_*
    for (int comp = 0; comp < this->numScalars(); ++comp) {
        sprintf(comp_str, "scalarComponent_%d", comp);

        if (header.m_string.find(comp_str) == header.m_string.end()) {
            ostringstream msg;
            msg << "Checkfile does not have " << comp_str;
            MayDay::Error(msg.str().c_str());
        }
        std::string scalarName = header.m_string[comp_str];
        if (scalarName.compare(this->getScalarName(comp)) != 0) {
            ostringstream msg;
            msg << "getScalarName(" << comp << ") changed from " << scalarName
                << " to " << this->getScalarName(comp) << ".";
            MayDay::Warning(msg.str().c_str());
        }
    }

    // numQComps
    if (header.m_int.find("numQComps") == header.m_int.end()) {
        MayDay::Error("Checkfile does not have numQComps");
    } else {
        int numQComps = header.m_int["numQComps"];
        if (numQComps != m_statePtr->numQComps) {
            ostringstream msg;
            msg << "numQComps changed from " << numQComps << " to "
                << m_statePtr->numQComps << ".";
            MayDay::Warning(msg.str().c_str());
        }
    }
}

// -----------------------------------------------------------------------------
// read checkpoint data for this level
// -----------------------------------------------------------------------------
void
AMRNSLevel::readCheckpointLevel(const std::string& a_filename)
{
    BEGIN_FLOWCHART();

#ifdef CH_USE_PYTHON
    const ProblemContext* ctx = ProblemContext::getInstance();

    // Open this level's group.

    // Read this group's metadata
    HeaderData header(a_filename, m_level, s_verbosity);
    if (s_verbosity >= 1) {
        pout() << "hdf5 header data:" << endl;
        pout() << header << endl;
    }
    // ref_ratio
    if (header.m_intvect.find("ref_ratio") == header.m_intvect.end()) {
        MayDay::Error("Checkfile does not have ref_ratio");
    } else {
        IntVect ref_ratio = header.m_intvect["ref_ratio"];
        if (ref_ratio != this->getFineRefRatio()) {
            // If ref_ratio (old value) is (1,1,1), then we must be adding a
            // level. This is okay.
            if (ref_ratio == IntVect::Unit) {
                // We must be adding a level. This is okay.

            } else if (this->getFineRefRatio() == IntVect::Unit) {
                // We must be removing a level. For now, we do not allow this.
                ostringstream msg;
                msg << "ref_ratio changed from " << ref_ratio << " to "
                    << this->getFineRefRatio() << ".\n"
                    << "If you are trying to remove a level, you can just "
                       "avoid "
                    << "tagging cells." << endl;
                MayDay::Error(msg.str().c_str());

            } else {
                // We must be changing the refinemnent ratios. This is NOT okay.
                ostringstream msg;
                msg << "ref_ratio changed from " << ref_ratio << " to "
                    << this->getFineRefRatio() << ".";
                MayDay::Error(msg.str().c_str());
            }
        }
    }

    // vec_dx
    if (header.m_realvect.find("vec_dx") == header.m_realvect.end()) {
        MayDay::Error("Checkfile does not have vec_dx");
    }
    const RealVect vec_dx = header.m_realvect["vec_dx"];  // Check later.

    // dt
    if (header.m_real.find("dt") == header.m_real.end()) {
        MayDay::Error("Checkfile does not have dt");
    }
    m_dt = header.m_real["dt"];

    // time
    if (header.m_real.find("time") == header.m_real.end()) {
        MayDay::Error("Checkfile does not have time");
    }
    m_time = header.m_real["time"];

    // isEmpty
    if (header.m_int.find("isEmpty") == header.m_int.end()) {
        MayDay::Error("Checkfile does not have isEmpty");
    }
    const bool __attribute__((unused)) isEmpty = (header.m_int["isEmpty"] != 0);

    // domainBox
    if (header.m_box.find("domainBox") == header.m_box.end()) {
        MayDay::Error("Checkfile does not have domainBox");
    }
    const Box domainBox = header.m_box["domainBox"];
    m_problem_domain    = ProblemDomain(domainBox, ctx->base.isPeriodic);

    // finestExtantLevel
    if (header.m_int.find("finestExtantLevel") == header.m_int.end()) {
        MayDay::Error("Checkfile does not have finestExtantLevel");
    }
    // const int finestExtantLevel = header.m_int["finestExtantLevel"];

    // Read and create level grids

    const int status =
        Py::PythonReturnFunction<int>("IO", "OpenFileForRead", a_filename);
    if (status < 0) MayDay::Error("We can't seem to read data from a_filename");

    const int grid_size =
      Py::PythonReturnFunction<int>("IO", "SizeOfBoxLayout", a_filename, m_level);
    if (grid_size == 0) {
        MayDay::Warning("Checkfile does not contain a disjointBoxLayout");
    }
    Vector<Box> boxArrayFromFile(grid_size);
    for (unsigned int i = 0; i < boxArrayFromFile.size(); ++i) {
        int j = (int)i;
        boxArrayFromFile[i] =
	  Py::PythonReturnFunction<Box>("IO", "GetBox", a_filename, m_level, j);
      }
      Vector<int> procMap(grid_size);
      for (unsigned int i = 0; i < boxArrayFromFile.size(); ++i)
      {
          int j = (int)i;
          procMap[i] =
	    Py::PythonReturnFunction<int>("IO", "getProcID", a_filename, m_level, j);
      }

          // If this level has valid data, read it now.
          if (boxArrayFromFile.size() > 0) {
              // Allocate and define everything on this level.
              // These must be called from bottom to top level.

              this->activateLevel(boxArrayFromFile, procMap);


              // Final checks...
              for (int dir = 0; dir < SpaceDim; ++dir) {
                  if (RealCmp::neq(vec_dx[dir], m_levGeoPtr->getDXi(dir))) {
                      ostringstream msg;
                      msg << "vec_dx[" << dir << "] changed from "
                          << vec_dx[dir] << " to " << m_levGeoPtr->getDXi(dir)
                          << ".";
                      MayDay::Error(msg.str().c_str());
                  }
              }
              // Read the data
              std::string vel = std::string("velData");
              std::string p   = std::string("pData");
              std::string q   = std::string("qData");
              Py::PythonFunction("IO", "ReadFB", a_filename, m_level, *m_velPtr, vel);
              Py::PythonFunction("IO", "ReadFAB", a_filename,  m_level, *m_pPtr, p);
              Py::PythonFunction("IO", "ReadFAB", a_filename, m_level, *m_qPtr, q);

              Py::PythonFunction("IO", "CloseFile", a_filename);
              m_velPtr->exchange();
              m_qPtr->exchange();
              m_pPtr->exchange();

              // Reset BCs.
              this->setBC(*m_statePtr, m_time);

          } else {
              // Level is empty.

	    Py::PythonFunction("IO", "CloseFile", a_filename);
              this->deactivateLevel();
          }
#endif
}


// -----------------------------------------------------------------------------
// Open plotfile. Only called on level 0.
// -----------------------------------------------------------------------------
void
AMRNSLevel::openFile(const std::string& a_filename, const bool checkpoint) const
{
#ifdef CH_USE_PYTHON
// if file exists, delete it
    Py::PythonFunction("IO", "DeleteIfFileExists", a_filename);
    Py::PythonFunction("IO", "OpenFileForWrite", a_filename, checkpoint);
#endif
}


// -----------------------------------------------------------------------------
// Open plotfile. Only called on level 0.
// -----------------------------------------------------------------------------
void
AMRNSLevel::closeFile(const std::string& a_filename) const
{
#ifdef CH_USE_PYTHON
    Py::PythonFunction("IO", "CloseFile", a_filename);
#endif
}


// -----------------------------------------------------------------------------
// Write plotfile header. Only called on level 0.
// -----------------------------------------------------------------------------
void
AMRNSLevel::writePlotHeader(HeaderData&        a_header,
                            const std::string& a_filename) const
{
    BEGIN_FLOWCHART();

#ifdef CH_USE_PYTHON

    // write some Chombo boilerplate to it (do we need it?)
    {
        HeaderData temp;
        temp.m_int["SpaceDim"]  = SpaceDim;
        temp.m_real["testReal"] = (Real)0.0;
        temp.writeToFile(a_filename, std::string("Chombo_global"));
    }
    auto& header = a_header;
    char comp_str[30];

    // This MUST match numComps in writePlotLevel.
    const int numComps = this->numScalars()
                       + SpaceDim  // u,(v),w
                       + 1         // divVel
                       + 1         // p
                       + 3         // T,S,b totals
                       + 3         // T,S,b perturbations
                       + 1         // eddy viscosity
#ifdef WRITE_GRADP_TO_HDF5
                       + SpaceDim  // gradP
#endif //WRITE_GRADP_TO_HDF5
#ifdef WRITE_METRIC_TO_HDF5
                       + SpaceDim  // x^i(Xi)
                       + SpaceDim  // dx^i/dXi^i
                       + SpaceDim  // dXi^i/dx^i
                       + SpaceDim  // Jg^{ii}
                       + SpaceDim  // g^{ii}
                       + SpaceDim  // g_{ii}
                       + 1         // J
                       + 1         // 1/J
#endif //WRITE_METRIC_TO_HDF5
#ifdef WRITE_Sij_TO_HDF5
                       + SpaceDim * (SpaceDim + 1) / 2 // Sij
#endif //WRITE_Sij_TO_HDF5
#ifdef WRITE_VORTICITY_TO_HDF5
                       + ((SpaceDim == 2) ? 1 : 3)
#endif //WRITE_VORTICITY_TO_HDF5
                       + SpaceDim; // Displacement
    header.m_int["num_components"] = numComps;

    // User-defined scalars
    int comp = 0;
    for (int sc = 0; sc < this->numScalars(); ++sc) {
        sprintf(comp_str, "component_%d", comp);
        header.m_string[comp_str] = this->getScalarName(sc);
        comp++;
    }

    // Velocity
    D_TERM(sprintf(comp_str, "component_%d", comp);
           header.m_string[comp_str] = "x_vel";
           comp++;
           , sprintf(comp_str, "component_%d", comp);
           header.m_string[comp_str] = "y_vel";
           comp++;
           , sprintf(comp_str, "component_%d", comp);
           header.m_string[comp_str] = "z_vel";
           comp++;)

    // divVel
    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "div_vel";
    comp++;

    // Pressure
    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "pressure";
    comp++;

    // Total temperature
    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "T_total";
    comp++;

    // Total salinity
    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "S_total";
    comp++;

    // Total buoyancy
    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "b_total";
    comp++;

    // Temperature perturbation
    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "T_pert";
    comp++;

    // Salinity perturbation
    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "S_pert";
    comp++;

    // Buoyancy perturbation
    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "b_pert";
    comp++;

    // Eddy viscosity
    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "eddyNu";
    comp++;

#ifdef WRITE_GRADP_TO_HDF5
    // gradP
    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "x_gradP";
    comp++;

    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "y_gradP";
    comp++;

    if (SpaceDim > 2) {
        sprintf(comp_str, "component_%d", comp);
        header.m_string[comp_str] = "z_gradP";
        comp++;
    }
#endif //WRITE_GRADP_TO_HDF5

#ifdef WRITE_METRIC_TO_HDF5
    // physCoor
    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "x_physCoor";
    comp++;

    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "y_physCoor";
    comp++;

    if (SpaceDim > 2) {
        sprintf(comp_str, "component_%d", comp);
        header.m_string[comp_str] = "z_physCoor";
        comp++;
    }

    // dxdXi
    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "x_dxdXi";
    comp++;

    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "y_dxdXi";
    comp++;

    if (SpaceDim > 2) {
        sprintf(comp_str, "component_%d", comp);
        header.m_string[comp_str] = "z_dxdXi";
        comp++;
    }

    // dXidx
    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "x_dXidx";
    comp++;

    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "y_dXidx";
    comp++;

    if (SpaceDim > 2) {
        sprintf(comp_str, "component_%d", comp);
        header.m_string[comp_str] = "z_dXidx";
        comp++;
    }

    // Jgup
    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "x_Jgup";
    comp++;

    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "y_Jgup";
    comp++;

    if (SpaceDim > 2) {
        sprintf(comp_str, "component_%d", comp);
        header.m_string[comp_str] = "z_Jgup";
        comp++;
    }

    // gup
    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "x_gup";
    comp++;

    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "y_gup";
    comp++;

    if (SpaceDim > 2) {
        sprintf(comp_str, "component_%d", comp);
        header.m_string[comp_str] = "z_gup";
        comp++;
    }

    // gdn
    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "x_gdn";
    comp++;

    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "y_gdn";
    comp++;

    if (SpaceDim > 2) {
        sprintf(comp_str, "component_%d", comp);
        header.m_string[comp_str] = "z_gdn";
        comp++;
    }

    // J
    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "J";
    comp++;

    // Jinv
    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "Jinv";
    comp++;
#endif //WRITE_METRIC_TO_HDF5

#ifdef WRITE_Sij_TO_HDF5
    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "S_xx";
    comp++;

    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "S_yy";
    comp++;

#if CH_SPACEDIM > 2
    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "S_zz";
    comp++;
#endif

    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "S_xy";
    comp++;

#if CH_SPACEDIM > 2
    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "S_zx";
    comp++;

    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "S_yz";
    comp++;
#endif
#endif //WRITE_Sij_TO_HDF5

#ifdef WRITE_VORTICITY_TO_HDF5
#if CH_SPACEDIM > 2
    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "vorticity_x";
    comp++;

    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "vorticity_y";
    comp++;
#endif

    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "vorticity_z";
    comp++;
#endif //WRITE_VORTICITY_TO_HDF5

    // Displacement
    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "x_Displacement";
    comp++;

    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "y_Displacement";
    comp++;

    if constexpr (SpaceDim > 2) {
        sprintf(comp_str, "component_%d", comp);
        header.m_string[comp_str] = "z_Displacement";
        comp++;
    }

    CH_assert(comp == numComps);

    header.writeToFile(a_filename, std::string("/"));

    if (s_verbosity >= 5) {
        pout() << header << endl;
    }
#endif //CH_USE_PYTHON
}


// -----------------------------------------------------------------------------
// Write plotfile data for this level. Called from bottom->top.
// -----------------------------------------------------------------------------
void
AMRNSLevel::writePlotLevel(const std::string& a_filename, int /*level*/) const
{
    BEGIN_FLOWCHART();

    char level_str[20];
    sprintf(level_str, "%d", m_level);
    const std::string label = std::string("level_") + level_str;

    HeaderData header;

    header.m_intvect["ref_ratio"] = this->getFineRefRatio();
    header.m_realvect["vec_dx"]   = m_levGeoPtr->getDXi();
    header.m_real["dt"]           = m_dt;
    header.m_real["time"]         = m_time;
    header.m_box["prob_domain"]   = m_problem_domain.domainBox();
    header.writeToFile(a_filename, label);

    // if (s_verbosity >= 5) {
    //   pout() << header << endl;
    // }

    // This MUST match numComps in writePlotHeader.
    const int numComps = this->numScalars()
                       + SpaceDim  // u,(v),w
                       + 1         // divVel
                       + 1         // p
                       + 3         // T,S,b totals
                       + 3         // T,S,b perturbations
                       + 1         // eddy viscosity
#ifdef WRITE_GRADP_TO_HDF5
                       + SpaceDim  // gradP
#endif //WRITE_GRADP_TO_HDF5
#ifdef WRITE_METRIC_TO_HDF5
                       + SpaceDim  // x^i(Xi)
                       + SpaceDim  // dx^i/dXi^i
                       + SpaceDim  // dXi^i/dx^i
                       + SpaceDim  // Jg^{ii}
                       + SpaceDim  // g^{ii}
                       + SpaceDim  // g_{ii}
                       + 1         // J
                       + 1         // 1/J
#endif //WRITE_METRIC_TO_HDF5
#ifdef WRITE_Sij_TO_HDF5
                       + SpaceDim * (SpaceDim + 1) / 2
#endif //WRITE_Sij_TO_HDF5
#ifdef WRITE_VORTICITY_TO_HDF5
                       + ((SpaceDim == 2) ? 1 : 3)
#endif //WRITE_VORTICITY_TO_HDF5
                       + SpaceDim; // Displacement

    const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();
    DataIterator dit = grids.dataIterator();

    // One ghost helps VisIt's interpolation
    LevelData<FArrayBox> plotData(grids, numComps, IntVect::Unit);

    // Two ghosts help us convert the FC vel to CC with fourth-order accuracy.
    const LevelData<FluxBox>& cartVel = m_statePtr->vel;
    LevelData<FluxBox> tmpVel(grids, 1, 2 * IntVect::Unit);

    // Set BCs
    this->setBC(*m_statePtr, m_time);


    // User-defined scalars
    int comp = 0;
    if (this->numScalars() > 0) {
        LevelData<FArrayBox> dest;
        const Interval destIvl(comp, comp + this->numScalars() - 1);
        aliasLevelData(dest, &plotData, destIvl);

        for (dit.reset(); dit.ok(); ++dit) {
            dest[dit].copy(m_statePtr->scalars[dit]);
        }

        comp += this->numScalars();
    }

    // Velocity
    {
        LevelData<FArrayBox> dest;
        aliasLevelData(dest, &plotData, Interval(comp, comp + SpaceDim - 1));

        // cartVel.copyTo(tmpVel);
        // m_levGeoPtr->multByJ(tmpVel);
        // Convert::FacesToCells(dest, tmpVel);
        // m_levGeoPtr->divByJ(dest);

        for (dit.reset(); dit.ok(); ++dit) {
            for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
                FArrayBox& tmpFAB = tmpVel[dit][fcDir];
                const FArrayBox& srcFAB = cartVel[dit][fcDir];
                const Box srcBox = grids[dit].grow(1).surroundingNodes(fcDir);

                tmpFAB.copy(srcFAB, srcBox);
                m_levGeoPtr->multByJ(tmpFAB, dit());
            }
        }
        BCTools::extrapAllGhosts(tmpVel, 2, IntVect::Unit);
        tmpVel.exchange();
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox&     destFAB = dest[dit];
            const Box      destBox = grids[dit];
            constexpr Real alpha   = 7.0 / 12.0;
            constexpr Real beta    = -1.0 / 12.0;

            for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
                FArrayBox&     tmpFAB = tmpVel[dit][fcDir];
                const IntVect& e      = BASISV(fcDir);

                for (BoxIterator bit(destBox); bit.ok(); ++bit) {
                    const IntVect& cc = bit();
                    destFAB(cc, fcDir) = alpha * (tmpFAB(cc +   e) + tmpFAB(cc    ))
                                       +  beta * (tmpFAB(cc + 2*e) + tmpFAB(cc - e));
                }
            }

            m_levGeoPtr->divByJ(destFAB, dit());
        }

        comp += SpaceDim;
    }

    // divVel
    {
        LevelData<FArrayBox> dest;
        aliasLevelData(dest, &plotData, Interval(comp, comp));


        // Gather fine vel and convert to advecting vel.
        LevelData<FluxBox>* fineVelPtr = nullptr;
        AMRNSLevel* fineLevPtr = this->fineNSPtr();
        if (fineLevPtr) {
            fineVelPtr = &(fineLevPtr->m_statePtr->vel);
            fineLevPtr->sendToAdvectingVelocity(*fineVelPtr, *fineVelPtr);
        }
        this->sendToAdvectingVelocity(tmpVel, cartVel);

        // Compute divergence.
        if (m_projOpPtr) {
            m_projOpPtr->compDivergence(dest, tmpVel, fineVelPtr);
            m_levGeoPtr->divByJ(dest);
        } else {
            const bool scaleByJinv = true;
            m_finiteDiffPtr->levelDivergenceMAC(dest, tmpVel, scaleByJinv);
        }

        // Restore Cartesian basis.
        if (fineLevPtr) {
            fineLevPtr->sendToCartesianVelocity(*fineVelPtr, *fineVelPtr);
            fineVelPtr = nullptr;
        }
        tmpVel.clear(); // Done with this.

        // Don't let invalid values make VisIt think max(div) is large.
        if (!this->isFinestLevel()) {
            Masks::zeroInvalid(dest, this->getFineGridsPtr());
        }

        comp += 1;
    }

    // pressure
    {
        LevelData<FArrayBox> dest;
        aliasLevelData(dest, &plotData, Interval(comp, comp));

        for (dit.reset(); dit.ok(); ++dit) {
            dest[dit].copy(m_statePtr->p[dit]);
        }

        comp += 1;
    }

    // Total temperature
    {
        LevelData<FArrayBox> dest;
        aliasLevelData(dest, &plotData, Interval(comp, comp));

        for (dit.reset(); dit.ok(); ++dit) {
            dest[dit].copy(m_statePtr->T[dit]);
        }

        comp += 1;
    }

    // Total salinity
    {
        LevelData<FArrayBox> dest;
        aliasLevelData(dest, &plotData, Interval(comp, comp));

        for (dit.reset(); dit.ok(); ++dit) {
            dest[dit].copy(m_statePtr->S[dit]);
        }

        comp += 1;
    }

    // Total buoyancy
    LevelData<FArrayBox> b; // We're gonna need this later.
    {
        aliasLevelData(b, &plotData, Interval(comp, comp));

        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox zFAB(b[dit].box(), 1);
            m_levGeoPtr->fill_physCoor(zFAB, 0, SpaceDim - 1);
            this->equationOfState(
                b[dit], m_statePtr->T[dit], m_statePtr->S[dit], zFAB);
        }

        comp += 1;
    }

    // Temperature perturbation
    {
        LevelData<FArrayBox> dest;
        aliasLevelData(dest, &plotData, Interval(comp, comp));

        for (dit.reset(); dit.ok(); ++dit) {
            dest[dit].copy(m_statePtr->T[dit]);
            Subspace::addHorizontalExtrusion(
                dest[dit], 0, *m_TbarPtr, 0, 1, -1.0);
        }

        comp += 1;
    }

    // Salinity perturbation
    {
        LevelData<FArrayBox> dest;
        aliasLevelData(dest, &plotData, Interval(comp, comp));

        for (dit.reset(); dit.ok(); ++dit) {
            dest[dit].copy(m_statePtr->S[dit]);
            Subspace::addHorizontalExtrusion(
                dest[dit], 0, *m_SbarPtr, 0, 1, -1.0);
        }

        comp += 1;
    }

    // Buoyancy perturbation
    {
        LevelData<FArrayBox> dest;
        aliasLevelData(dest, &plotData, Interval(comp, comp));

        for (dit.reset(); dit.ok(); ++dit) {
            dest[dit].copy(b[dit]);
            Subspace::addHorizontalExtrusion(
                dest[dit], 0, *m_bbarPtr, 0, 1, -1.0);
        }

        comp += 1;
    }

    // eddyNu
    {
        // In incremental PARK, eddyNu is already filled, but in standard form,
        // eddyNu is all zeros. This happens because RK in standard form always
        // starts with q^{n} and adds forces to it. But there are no forces that
        // add to the eddy viscosity, it's just something we calculate at each
        // stage. So, when it comes time to assemble q^{n+1} = q^{n} + dt * forces,
        // the eddy viscosity never gets copied over from the RK stages.
        //
        // Don't worry, a new eddyNu is computed each time it is used dynamically,
        // so your simulation is fine in either incremental or standard mode.
        // We just need to remember to fill eddyNu for post-processing.
        this->computeEddyNu(m_statePtr->eddyNu, cartVel, m_time);

        LevelData<FArrayBox> dest;
        aliasLevelData(dest, &plotData, Interval(comp, comp));

        for (dit.reset(); dit.ok(); ++dit) {
            dest[dit].copy(m_statePtr->eddyNu[dit]);
        }

        comp += 1;
    }

#ifdef WRITE_GRADP_TO_HDF5
    // gradP
    {
        LevelData<FArrayBox> dest;
        aliasLevelData(dest, &plotData, Interval(comp, comp + SpaceDim - 1));

        // Since the BCs are been set, this is the same as a composite gradient.
        LevelData<FluxBox> gradP(grids, 1);
        m_finiteDiffPtr->levelGradientMAC(gradP, m_statePtr->p);
        this->sendToCartesianVelocity(gradP, gradP);
        m_levGeoPtr->multByJ(gradP);
        Convert::FacesToCells(dest, gradP);
        m_levGeoPtr->divByJ(dest);

        // If you want to know what the level-projector sees, comment this out.
        // if (!this->isFinestLevel()) {
        //     Masks::zeroInvalid(dest, m_levGeoPtr->getFineGridsPtr());
        // }

        comp += SpaceDim;
    }
#endif //WRITE_GRADP_TO_HDF5

#ifdef WRITE_METRIC_TO_HDF5
    // physCoor
    {
        LevelData<FArrayBox> dest;
        aliasLevelData(dest, &plotData, Interval(comp, comp + SpaceDim - 1));

        for (dit.reset(); dit.ok(); ++dit) {
            m_levGeoPtr->fill_physCoor(dest[dit]);
        }

        comp += SpaceDim;
    }

    // dxdXi
    {
        LevelData<FArrayBox> dest;
        aliasLevelData(dest, &plotData, Interval(comp, comp + SpaceDim - 1));

        for (dit.reset(); dit.ok(); ++dit) {
            for (int d = 0; d < SpaceDim; ++d) {
                m_levGeoPtr->getGeoSource().fill_dxdXi(
                    dest[dit], d, d, m_levGeoPtr->getDXi());
            }
        }

        comp += SpaceDim;
    }

    // dXidx
    {
        LevelData<FArrayBox> dest;
        aliasLevelData(dest, &plotData, Interval(comp, comp + SpaceDim - 1));

        for (dit.reset(); dit.ok(); ++dit) {
            for (int d = 0; d < SpaceDim; ++d) {
                m_levGeoPtr->getGeoSource().fill_dXidx(
                    dest[dit], d, d, m_levGeoPtr->getDXi());
            }
        }

        comp += SpaceDim;
    }

    // Jgup
    {
        LevelData<FArrayBox> dest;
        aliasLevelData(dest, &plotData, Interval(comp, comp + SpaceDim - 1));
        Convert::FacesToCells(dest, m_levGeoPtr->getFCJgup());
        // for (dit.reset(); dit.ok(); ++dit) {
        //     for (int d = 0; d < SpaceDim; ++d) {
        //         m_levGeoPtr->getGeoSource().fill_Jgup(
        //             dest[dit], d, d, m_levGeoPtr->getDXi());
        //     }
        // }

        comp += SpaceDim;
    }

    // gup
    {
        LevelData<FArrayBox> dest;
        aliasLevelData(dest, &plotData, Interval(comp, comp + SpaceDim - 1));

        for (dit.reset(); dit.ok(); ++dit) {
            for (int d = 0; d < SpaceDim; ++d) {
                m_levGeoPtr->getGeoSource().fill_gup(
                    dest[dit], d, d, m_levGeoPtr->getDXi());
            }
        }

        comp += SpaceDim;
    }

    // gdn
    {
        LevelData<FArrayBox> dest;
        aliasLevelData(dest, &plotData, Interval(comp, comp + SpaceDim - 1));

        for (dit.reset(); dit.ok(); ++dit) {
            for (int d = 0; d < SpaceDim; ++d) {
                m_levGeoPtr->getGeoSource().fill_gdn(
                    dest[dit], d, d, m_levGeoPtr->getDXi());
            }
        }

        comp += SpaceDim;
    }

    // J
    {
        LevelData<FArrayBox> dest;
        aliasLevelData(dest, &plotData, Interval(comp, comp));

        for (dit.reset(); dit.ok(); ++dit) {
            dest[dit].copy(m_levGeoPtr->getCCJ()[dit]);
        }

        comp += 1;
    }

    // Jinv
    {
        LevelData<FArrayBox> dest;
        aliasLevelData(dest, &plotData, Interval(comp, comp));

        for (dit.reset(); dit.ok(); ++dit) {
            dest[dit].copy(m_levGeoPtr->getCCJinv()[dit]);
        }

        comp += 1;
    }
#endif //WRITE_METRIC_TO_HDF5


#ifdef WRITE_Sij_TO_HDF5
    {
        const RealVect nuDummy(D_DECL(0.5, 0.5, 0.5));

        LevelData<FArrayBox> eddyNuDummy(grids, 1, IntVect::Unit);
        setValLevel(eddyNuDummy, 0.0);

        LevelData<FluxBox> kvel(grids, 1);
        StaggeredFluxLD Sij(grids);
        setValLevel(kvel, 0.0);
        Sij.setVal(0.0);
        this->computeMomentumDiffusion(kvel, Sij, cartVel, nuDummy, eddyNuDummy);

        for (size_t i = 0; i < SpaceDim; ++i) {
            LevelData<FArrayBox> dest;
            aliasLevelData(dest, &plotData, Interval(comp, comp));

            for (dit.reset(); dit.ok(); ++dit) {
                dest[dit].copy(Sij[i][i][dit]);
            }
            ++comp;
        }

        for (size_t i = 0; i < SpaceDim; ++i) {
            for (size_t j = i + 1; j < SpaceDim; ++j) {
                LevelData<FArrayBox> dest;
                aliasLevelData(dest, &plotData, Interval(comp, comp));

                for (dit.reset(); dit.ok(); ++dit) {
                    Convert::Simple(dest[dit], grids[dit], Sij[i][j][dit]);
                }
                ++comp;
            }
        }
    }
#endif //WRITE_Sij_TO_HDF5

#ifdef WRITE_VORTICITY_TO_HDF5
    if (SpaceDim == 2) {
        LevelData<FArrayBox> dest;
        aliasLevelData(dest, &plotData, Interval(comp, comp));
        ++comp;

        const IntVect   ex  = BASISV(0);
        const IntVect   ey  = BASISV(1);
        const RealVect& dXi = m_levGeoPtr->getDXi();

        for (dit.reset(); dit.ok(); ++dit) {
            const Box ccBox   = grids[dit];
            const Box ecBox01 = surroundingNodes(ccBox, 0).surroundingNodes(1);
            const Box fcBox0  = enclosedCells(ecBox01, 1).grow(ey);
            const Box fcBox1  = enclosedCells(ecBox01, 0).grow(ex);

            const FArrayBox& contraUFAB = cartVel[dit][0];
            const FArrayBox& contraVFAB = cartVel[dit][1];
            CH_assert(contraUFAB.box().contains(fcBox0));
            CH_assert(contraVFAB.box().contains(fcBox1));

            FArrayBox coUFAB(fcBox0, 1);
            m_levGeoPtr->getGeoSource().fill_gdn(coUFAB, 0, 0, dXi);
            coUFAB.mult(contraUFAB, fcBox0, 0, 0, 1);

            FArrayBox coVFAB(fcBox1, 1);
            m_levGeoPtr->getGeoSource().fill_gdn(coVFAB, 0, 1, dXi);
            coVFAB.mult(contraVFAB, fcBox1, 0, 0, 1);

            FArrayBox dvelFAB(ecBox01, 1);
            FiniteDiff::partialD(dvelFAB, 0, ecBox01, coUFAB, 0, 1, dXi[1]);
            FiniteDiff::partialD(dvelFAB, 0, ecBox01, coVFAB, 0, 0, dXi[0], true);

            Convert::Simple(dest[dit], ccBox, dvelFAB);
            m_levGeoPtr->divByJ(dest[dit], dit());
        }

    } else {
        UNDEFINED_FUNCTION();
    }
#endif //WRITE_VORTICITY_TO_HDF5

    // Displacement
    {
        LevelData<FArrayBox> dest;
        aliasLevelData(dest, &plotData, Interval(comp, comp + SpaceDim - 1));

        for (dit.reset(); dit.ok(); ++dit) {
            m_levGeoPtr->fill_displacement(dest[dit]);
        }

        comp += SpaceDim;
    }

    CH_assert(comp == numComps);
    std::string name = "data";

    BCTools::extrapAllGhosts(plotData, 2);
    plotData.exchange();

#ifdef CH_USE_PYTHON
    // Bottleneck! (But we don't need to call this often.)
    Py::PythonFunction("IO",
                       "WriteLevelDataFAB",
                       a_filename,
                       label,
                       name,
                       plotData,
                       IntVect::Unit);
#endif
}
