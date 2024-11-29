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
#include "RHSParameters.H"
#include "AMRParameters.H"
#include "Debug.H"
#include "Format.H"
#include "ParmParse.H"
#include "SOMAR_Constants.H"


//------------------------------------------------------------------------------
// Static variable definitions for RHSParameters.
//------------------------------------------------------------------------------
std::unique_ptr<RHSParameters> RHSParameters::s_defPtr;
bool                           RHSParameters::s_constructorLock = false;


//------------------------------------------------------------------------------
// The default constructor sets default / reads parameters.
//------------------------------------------------------------------------------
RHSParameters::RHSParameters(const BaseParameters& a_baseParams,
                             const int             a_numLevels)
{
    if (s_constructorLock) return;

    if (!s_defPtr) {
        s_constructorLock = true;

        s_defPtr.reset(new RHSParameters(a_baseParams, a_numLevels));
        createDefaults(a_baseParams, a_numLevels);

        s_constructorLock = false;
    }

    // Copy default values.
    *this = *s_defPtr;
}


//------------------------------------------------------------------------------
// For the bean counters like me.
//------------------------------------------------------------------------------
void
RHSParameters::freeMemory()
{
    s_defPtr.reset();
}


//------------------------------------------------------------------------------
// It's nice to be able to see these parameters in pout.*.
//------------------------------------------------------------------------------
void
RHSParameters::dump() const
{
    pout() << "RHSParameters:\n" << Format::indent() << std::flush;

    pout() << "nu = " << nu << "\n";
    pout() << "TKappa = " << TKappa << "\n";
    pout() << "SKappa = " << SKappa << "\n";
    if (scalarsKappa.size() > 0) {
        pout() << "scalarsKappa = " << scalarsKappa << "\n";
    } else {
        pout() << "scalarsKappa = ()\n";
    }
    pout() << "coriolisF = " << coriolisF << "\n";


    pout() << "doImplicitDiffusion = " << doImplicitDiffusion << '\n';

    pout() << "velBCTypeLo = (";
    for (int dir = 0; dir < SpaceDim; ++dir) {
        switch (velBCTypeLo[dir]) {
            case VelBCType::UNDEFINED: pout() << "UNDEFINED"; break;
            case VelBCType::PERIODIC : pout() << " PERIODIC"; break;
            case VelBCType::NO_SLIP  : pout() << "  NO_SLIP"; break;
            case VelBCType::FREE_SLIP: pout() << "FREE_SLIP"; break;
            case VelBCType::CUSTOM   : pout() << "   CUSTOM"; break;
            case VelBCType::TIDAL    : pout() << "    TIDAL"; break;
            case VelBCType::OUTFLOW  : pout() << "  OUTFLOW"; break;
            default:
                MAYDAYERROR("velBCTypeLo[" << dir << "] = "
                            << velBCTypeLo[dir] << " not recognized.");
        }
        if (dir < SpaceDim - 1) pout() << ", ";
    }
    pout() << ")" << endl;

    pout() << "velBCTypeHi = (";
    for (int dir = 0; dir < SpaceDim; ++dir) {
        switch (velBCTypeHi[dir]) {
            case VelBCType::UNDEFINED: pout() << "UNDEFINED"; break;
            case VelBCType::PERIODIC : pout() << " PERIODIC"; break;
            case VelBCType::NO_SLIP  : pout() << "  NO_SLIP"; break;
            case VelBCType::FREE_SLIP: pout() << "FREE_SLIP"; break;
            case VelBCType::CUSTOM   : pout() << "   CUSTOM"; break;
            case VelBCType::TIDAL    : pout() << "    TIDAL"; break;
            case VelBCType::OUTFLOW  : pout() << "  OUTFLOW"; break;
            default:
                MAYDAYERROR("velBCTypeHi[" << dir << "] = "
                            << velBCTypeHi[dir] << " not recognized.");
        }
        if (dir < SpaceDim - 1) pout() << ", ";
    }
    pout() << ")" << endl;

    pout() << "tidalU0 = " << tidalU0
           << "\ntidalOmega = " << tidalOmega
           << "\ntidalInitPhase = " << tidalInitPhase
           << endl;


    pout() << "tempBCTypeLo = (";
    for (int dir = 0; dir < SpaceDim; ++dir) {
        switch (tempBCTypeLo[dir]) {
            case TempBCType::UNDEFINED        : pout() << "        UNDEFINED"; break;
            case TempBCType::PERIODIC         : pout() << "         PERIODIC"; break;
            case TempBCType::ZERO_NEUM        : pout() << "        ZERO_NEUM"; break;
            case TempBCType::ZERO_NEUM_ON_PERT: pout() << "ZERO_NEUM_ON_PERT"; break;
            case TempBCType::CUSTOM           : pout() << "           CUSTOM"; break;
            case TempBCType::TIDAL            : pout() << "            TIDAL"; break;
            case TempBCType::OUTFLOW          : pout() << "          OUTFLOW"; break;
            default:
                MAYDAYERROR("tempBCTypeLo[" << dir << "] = "
                            << tempBCTypeLo[dir] << " not recognized.");
        }
        if (dir < SpaceDim - 1) pout() << ", ";
    }
    pout() << ")" << endl;

    pout() << "tempBCTypeHi = (";
    for (int dir = 0; dir < SpaceDim; ++dir) {
        switch (tempBCTypeHi[dir]) {
            case TempBCType::UNDEFINED        : pout() << "        UNDEFINED"; break;
            case TempBCType::PERIODIC         : pout() << "         PERIODIC"; break;
            case TempBCType::ZERO_NEUM        : pout() << "        ZERO_NEUM"; break;
            case TempBCType::ZERO_NEUM_ON_PERT: pout() << "ZERO_NEUM_ON_PERT"; break;
            case TempBCType::CUSTOM           : pout() << "           CUSTOM"; break;
            case TempBCType::TIDAL            : pout() << "            TIDAL"; break;
            case TempBCType::OUTFLOW          : pout() << "          OUTFLOW"; break;
            default:
                MAYDAYERROR("tempBCTypeHi[" << dir << "] = "
                            << tempBCTypeHi[dir] << " not recognized.");
        }
        if (dir < SpaceDim - 1) pout() << ", ";
    }
    pout() << ")" << endl;


    pout() << "salinityBCTypeLo = (";
    for (int dir = 0; dir < SpaceDim; ++dir) {
        switch (salinityBCTypeLo[dir]) {
            case SalinityBCType::UNDEFINED        : pout() << "        UNDEFINED"; break;
            case SalinityBCType::PERIODIC         : pout() << "         PERIODIC"; break;
            case SalinityBCType::ZERO_NEUM        : pout() << "        ZERO_NEUM"; break;
            case SalinityBCType::ZERO_NEUM_ON_PERT: pout() << "ZERO_NEUM_ON_PERT"; break;
            case SalinityBCType::CUSTOM           : pout() << "           CUSTOM"; break;
            case SalinityBCType::TIDAL            : pout() << "            TIDAL"; break;
            case SalinityBCType::OUTFLOW          : pout() << "          OUTFLOW"; break;
            default:
                MAYDAYERROR("salinityBCTypeLo[" << dir << "] = "
                            << salinityBCTypeLo[dir] << " not recognized.");
        }
        if (dir < SpaceDim - 1) pout() << ", ";
    }
    pout() << ")" << endl;

    pout() << "salinityBCTypeHi = (";
    for (int dir = 0; dir < SpaceDim; ++dir) {
        switch (salinityBCTypeHi[dir]) {
            case SalinityBCType::UNDEFINED        : pout() << "        UNDEFINED"; break;
            case SalinityBCType::PERIODIC         : pout() << "         PERIODIC"; break;
            case SalinityBCType::ZERO_NEUM        : pout() << "        ZERO_NEUM"; break;
            case SalinityBCType::ZERO_NEUM_ON_PERT: pout() << "ZERO_NEUM_ON_PERT"; break;
            case SalinityBCType::CUSTOM           : pout() << "           CUSTOM"; break;
            case SalinityBCType::TIDAL            : pout() << "            TIDAL"; break;
            case SalinityBCType::OUTFLOW          : pout() << "          OUTFLOW"; break;
            default:
                MAYDAYERROR("salinityBCTypeHi[" << dir << "] = "
                            << salinityBCTypeHi[dir] << " not recognized.");
        }
        if (dir < SpaceDim - 1) pout() << ", ";
    }
    pout() << ")" << endl;


    for (int dir = 0; dir < SpaceDim; ++dir) {
        for (SideIterator sit; sit.ok(); ++sit) {
            const int iside = int(sit());
            pout() << "spongeWidth[" << dir << "]["
                   << (iside == 0 ? "Lo" : "Hi")
                   << "] = " << spongeWidth[dir][iside] << "\n";
        }
    }
    pout() << "spongeTimeCoeff = " << spongeTimeCoeff << "\n";

    pout() << "eddyViscMethod = (";
    for (unsigned int lev = 0; lev < s_defPtr->eddyViscMethod.size(); ++lev) {
        if (lev > 0) pout() << ", ";
        switch (s_defPtr->eddyViscMethod[lev]) {
            case EddyViscMethods::AVG_DOWN    : pout() << "AVG_DOWN"    ; break;
            case EddyViscMethods::SET_TO_ZERO : pout() << "SET_TO_ZERO" ; break;
            case EddyViscMethods::STATIC_SMAG : pout() << "STATIC_SMAG" ; break;
            case EddyViscMethods::DYNAMIC_SMAG: pout() << "DYNAMIC_SMAG"; break;
            case EddyViscMethods::DUCROS      : pout() << "DUCROS"      ; break;
        }
    }
    pout() << ")\n";
    pout() << "eddyPrandtlT = " << eddyPrandtlT << "\n";
    pout() << "eddyPrandtlS = " << eddyPrandtlS << "\n";
    pout() << "eddyPrandtlScalars = " << eddyPrandtlScalars << "\n";
    pout() << "eddyScale = " << eddyScale << "\n";

    if (velReconstruction == Reconstruction::SecondOrder) {
        pout() << "velReconstruction = SecondOrder\n";
    } else if (velReconstruction == Reconstruction::FourthOrder) {
        pout() << "velReconstruction = FourthOrder\n";
    } else {
        MAYDAYERROR("velReconstruction setting not recognized.");
    }

    if (scalReconstruction == Reconstruction::SecondOrder) {
        pout() << "scalReconstruction = SecondOrder\n";
    } else if (scalReconstruction == Reconstruction::FourthOrder) {
        pout() << "scalReconstruction = FourthOrder\n";
    } else {
        MAYDAYERROR("scalReconstruction setting not recognized.");
    }

    pout() << "momAdvSkewness = " << momAdvSkewness << '\n';

    pout() << "doMomentumAdvection = "
           << (doMomentumAdvection ? "true" : "false") << "\n";

    pout() << "doTemperatureAdvection = "
           << (doTemperatureAdvection ? "true" : "false") << "\n";

    pout() << "doSalinityAdvection = "
           << (doSalinityAdvection ? "true" : "false") << "\n";

    pout() << "doScalarAdvection = " << (doScalarAdvection ? "true" : "false")
           << "\n";

    pout() << "doViscousForcing = " << (doViscousForcing ? "true" : "false")
           << "\n";

    pout() << "doTemperatureDiffusion = "
           << (doTemperatureDiffusion ? "true" : "false") << "\n";

    pout() << "doSalinityDiffusion = "
           << (doSalinityDiffusion ? "true" : "false") << "\n";

    pout() << "doScalarDiffusion = " << (doScalarDiffusion ? "true" : "false")
           << "\n";

    pout() << "doGravityForcing = " << (doGravityForcing ? "true" : "false")
           << "\n";

    pout() << "doCoriolisForcing = " << (doCoriolisForcing ? "true" : "false")
           << "\n";

    pout() << "doTidalForcing = " << (doTidalForcing ? "true" : "false")
           << "\n";

    pout() << "doSpongeForcing = " << (doSpongeForcing ? "true" : "false")
           << "\n";

    pout() << "doMomAdvRefluxing = " << (doMomAdvRefluxing ? "true" : "false")
           << "\n";

    pout() << "doTemperatureAdvRefluxing = "
           << (doTemperatureAdvRefluxing ? "true" : "false") << "\n";

    pout() << "doSalnityAdvRefluxing = "
           << (doSalinityAdvRefluxing ? "true" : "false") << "\n";

    pout() << "doScalarAdvRefluxing = "
           << (doScalarAdvRefluxing ? "true" : "false") << "\n";

    pout() << "doViscousRefluxing = " << (doViscousRefluxing ? "true" : "false")
           << "\n";

    pout() << "doTemperatureDiffusiveRefluxing = "
           << (doTemperatureDiffusiveRefluxing ? "true" : "false") << "\n";

    pout() << "doSalinityDiffusiveRefluxing = "
           << (doSalinityDiffusiveRefluxing ? "true" : "false") << "\n";

    pout() << "doScalarDiffusiveRefluxing = "
           << (doScalarDiffusiveRefluxing ? "true" : "false") << "\n";


    pout() << "computeInitPressure = "
           << (computeInitPressure ? "true" : "false") << "\n";

    pout() << Format::unindent << std::endl;
}


//------------------------------------------------------------------------------
Real
RHSParameters::getScalarsKappa(const unsigned int a_comp) const
{
    Real retVal;
    if (scalarsKappa.size() == 1) {
        retVal = scalarsKappa[0];
    } else if (scalarsKappa.size() > a_comp) {
        retVal = scalarsKappa[a_comp];
    } else {
        MAYDAYERROR("scalarsKappa must have 1 or numScalars elements."
                    << "\n\tscalarsKappa.size() = " << scalarsKappa.size()
                    << "\n\ta_comp = " << a_comp);
        retVal = quietNAN; // Just to quiet a compiler warning.
    }
    return retVal;
}


//------------------------------------------------------------------------------
Real
RHSParameters::getEddyPrandtlScalars(const unsigned int a_comp) const
{
    Real retVal;
    if (eddyPrandtlScalars.size() == 1) {
        retVal = eddyPrandtlScalars[0];
    } else if (eddyPrandtlScalars.size() > a_comp) {
        retVal = eddyPrandtlScalars[a_comp];
    } else {
        MAYDAYERROR("eddyPrandtlScalars must have 1 or numScalars elements."
                    << "\n\teddyPrandtlScalars.size() = " << eddyPrandtlScalars.size()
                    << "\n\ta_comp = " << a_comp);
        retVal = quietNAN; // Just to quiet a compiler warning.
    }
    return retVal;
}


//------------------------------------------------------------------------------
// Fills the *s_defPtr object.
//------------------------------------------------------------------------------
void
RHSParameters::createDefaults(const BaseParameters& a_baseParams,
                              const int             a_numLevels)
{
    ParmParse pp("rhs");
    Vector<int> vint(SpaceDim);
    Vector<Real> vreal(SpaceDim);

    s_defPtr->nu = 0.0;
    pp.query("nu", s_defPtr->nu);

    s_defPtr->TKappa = 0.0;
    pp.query("TKappa", s_defPtr->TKappa);

    s_defPtr->SKappa = 0.0;
    pp.query("SKappa", s_defPtr->SKappa);

    int numVals = pp.countval("scalarsKappa");
    if (numVals > 0) {
        s_defPtr->scalarsKappa = Vector<Real>(numVals, 0.0);
        pp.queryarr("scalarsKappa", s_defPtr->scalarsKappa, 0, numVals);
    } else {
        s_defPtr->scalarsKappa.resize(0);
    }

    s_defPtr->doImplicitDiffusion = false;
    pp.query("doImplicitDiffusion", s_defPtr->doImplicitDiffusion);

    if (pp.queryarr("coriolisF", vreal, 0, SpaceDim)) {
        s_defPtr->coriolisF = RealVect(vreal);
    } else {
        s_defPtr->coriolisF = RealVect::Zero;
    }

    bool tidalBCPresent = false;
    vint = Vector<int>(SpaceDim, VelBCType::NO_SLIP);
    pp.queryarr("velBCTypeLo", vint, 0, SpaceDim);
    for (int dir = 0; dir < SpaceDim; ++dir) {
        if (a_baseParams.isPeriodic[dir]) {
            vint[dir] = VelBCType::PERIODIC;
        }

        if (vint[dir] == VelBCType::TIDAL) {
            tidalBCPresent = true;
        }

        if (vint[dir] < 0 || VelBCType::_NUM_BCTYPES <= vint[dir]) {
            MAYDAYERROR("velBCTypeLo[" << vint[dir] << "] = " << vint[dir]
                                       << " not recognized.");
        }
        s_defPtr->velBCTypeLo[dir] = vint[dir];
        s_defPtr->velBCType[dir][0] = vint[dir];
    }

    vint = Vector<int>(SpaceDim, VelBCType::NO_SLIP);
    pp.queryarr("velBCTypeHi", vint, 0, SpaceDim);
    for (int dir = 0; dir < SpaceDim; ++dir) {
        if (a_baseParams.isPeriodic[dir]) {
            vint[dir] = VelBCType::PERIODIC;
        }

        if (vint[dir] == VelBCType::TIDAL) {
            tidalBCPresent = true;
        }

        if (vint[dir] < 0 || VelBCType::_NUM_BCTYPES <= vint[dir]) {
            MAYDAYERROR("velBCTypeHi[" << dir << "] = " << vint[dir]
                                       << " not recognized.");
        }
        s_defPtr->velBCTypeHi[dir] = vint[dir];
        s_defPtr->velBCType[dir][1] = vint[dir];
    }


    vint = Vector<int>(SpaceDim, TempBCType::ZERO_NEUM);
    pp.queryarr("tempBCTypeLo", vint, 0, SpaceDim);
    for (int dir = 0; dir < SpaceDim; ++dir) {
        if (a_baseParams.isPeriodic[dir]) {
            vint[dir] = TempBCType::PERIODIC;
        }

        if (vint[dir] == TempBCType::TIDAL) {
            CH_verify(tidalBCPresent);
        }

        if (vint[dir] < 0 || TempBCType::_NUM_BCTYPES <= vint[dir]) {
            MAYDAYERROR("tempBCTypeLo[" << vint[dir] << "] = " << vint[dir]
                                        << " not recognized.");
        }
        s_defPtr->tempBCTypeLo[dir] = vint[dir];
        s_defPtr->tempBCType[dir][0] = vint[dir];
    }

    vint = Vector<int>(SpaceDim, TempBCType::ZERO_NEUM);
    pp.queryarr("tempBCTypeHi", vint, 0, SpaceDim);
    for (int dir = 0; dir < SpaceDim; ++dir) {
        if (a_baseParams.isPeriodic[dir]) {
            vint[dir] = TempBCType::PERIODIC;
        }

        if (vint[dir] == TempBCType::TIDAL) {
            CH_verify(tidalBCPresent);
        }

        if (vint[dir] < 0 || TempBCType::_NUM_BCTYPES <= vint[dir]) {
            MAYDAYERROR("tempBCTypeHi[" << dir << "] = " << vint[dir]
                                        << " not recognized.");
        }
        s_defPtr->tempBCTypeHi[dir] = vint[dir];
        s_defPtr->tempBCType[dir][1] = vint[dir];
    }


    vint = Vector<int>(SpaceDim, SalinityBCType::ZERO_NEUM);
    pp.queryarr("salinityBCTypeLo", vint, 0, SpaceDim);
    for (int dir = 0; dir < SpaceDim; ++dir) {
        if (a_baseParams.isPeriodic[dir]) {
            vint[dir] = SalinityBCType::PERIODIC;
        }

        if (vint[dir] == SalinityBCType::TIDAL) {
            CH_verify(tidalBCPresent);
        }

        if (vint[dir] < 0 || SalinityBCType::_NUM_BCTYPES <= vint[dir]) {
            MAYDAYERROR("salinityBCTypeLo[" << vint[dir] << "] = " << vint[dir]
                                            << " not recognized.");
        }
        s_defPtr->salinityBCTypeLo[dir] = vint[dir];
        s_defPtr->salinityBCType[dir][0] = vint[dir];
    }

    vint = Vector<int>(SpaceDim, SalinityBCType::ZERO_NEUM);
    pp.queryarr("salinityBCTypeHi", vint, 0, SpaceDim);
    for (int dir = 0; dir < SpaceDim; ++dir) {
        if (a_baseParams.isPeriodic[dir]) {
            vint[dir] = SalinityBCType::PERIODIC;
        }

        if (vint[dir] == SalinityBCType::TIDAL) {
            CH_verify(tidalBCPresent);
        }

        if (vint[dir] < 0 || SalinityBCType::_NUM_BCTYPES <= vint[dir]) {
            MAYDAYERROR("salinityBCTypeHi[" << dir << "] = " << vint[dir]
                                            << " not recognized.");
        }
        s_defPtr->salinityBCTypeHi[dir] = vint[dir];
        s_defPtr->salinityBCType[dir][1] = vint[dir];
    }

    CH_verify(a_numLevels > 0);
    s_defPtr->eddyViscMethod =
        Vector<int>(a_numLevels, EddyViscMethods::SET_TO_ZERO);
    if (pp.contains("eddyViscMethod")) {
        pp.queryarr("eddyViscMethod", s_defPtr->eddyViscMethod, 0, a_numLevels);

        bool usesEddyVisc = false;
        for (int lev = 0; lev < a_numLevels; ++lev) {
            CH_verify(s_defPtr->eddyViscMethod[lev] >= -1);
            CH_verify(s_defPtr->eddyViscMethod[lev] <
                      EddyViscMethods::_NUM_EDDYVISCMETHODS);

            usesEddyVisc |= (s_defPtr->eddyViscMethod[lev] != 0);
        }

        if (usesEddyVisc) {
            pp.get("eddyPrandtlT", s_defPtr->eddyPrandtlT);
            pp.get("eddyPrandtlS", s_defPtr->eddyPrandtlS);

            numVals = pp.countval("eddyPrandtlScalars");
            if (numVals > 0) {
                s_defPtr->eddyPrandtlScalars = Vector<Real>(numVals, quietNAN);
                pp.getarr("eddyPrandtlScalars", s_defPtr->eddyPrandtlScalars, 0, numVals);
            } else {
                s_defPtr->eddyPrandtlScalars.resize(0);
            }

            vreal = Vector<Real>(SpaceDim, 1.0);
            pp.queryarr("eddyScale", vreal, 0, SpaceDim);
            s_defPtr->eddyScale = RealVect(vreal);
        } else {
            s_defPtr->eddyPrandtlT = 1.0;
            s_defPtr->eddyPrandtlS = 1.0;
            s_defPtr->eddyPrandtlScalars = Vector<Real>(1, 1.0);
            s_defPtr->eddyScale =
                RealVect(D_DECL(quietNAN, quietNAN, quietNAN));
        }
    } else {
        s_defPtr->eddyPrandtlT = 1.0;
        s_defPtr->eddyPrandtlS = 1.0;
        s_defPtr->eddyPrandtlScalars = Vector<Real>(1, 1.0);
        s_defPtr->eddyScale =
            RealVect(D_DECL(quietNAN, quietNAN, quietNAN));
    }

    s_defPtr->velReconstruction = Reconstruction::FourthOrder;
    {
        int order = 4;
        pp.query("velReconstruction", order);
        if (order == 2) {
            s_defPtr->velReconstruction = Reconstruction::SecondOrder;
        } else if (order == 4) {
            s_defPtr->velReconstruction = Reconstruction::FourthOrder;
        } else {
            MAYDAYERROR("rhs.velReconstruction = "
                        << order << " not a valid option. Try 2 or 4.");
        }
    }

    s_defPtr->scalReconstruction = Reconstruction::FourthOrder;
    {
        int order = 4;
        pp.query("scalReconstruction", order);
        if (order == 2) {
            s_defPtr->scalReconstruction = Reconstruction::SecondOrder;
        } else if (order == 4) {
            s_defPtr->scalReconstruction = Reconstruction::FourthOrder;
        } else {
            MAYDAYERROR("rhs.scalReconstruction = "
                        << order << " not a valid option. Try 2 or 4.");
        }
    }

    s_defPtr->momAdvSkewness = 0.0; // Start with conservative form.
    pp.query("momAdvSkewness", s_defPtr->momAdvSkewness);
    if (s_defPtr->momAdvSkewness < 0.0 || 1.0 < s_defPtr->momAdvSkewness) {
        MAYDAYWARNING("rhs.momAdvSkewness = "
                      << s_defPtr->momAdvSkewness
                      << ". This is typically between 0 (conservative form) "
                         "and 1 (advective form).");
    }

    s_defPtr->doMomentumAdvection = true;
    pp.query("doMomentumAdvection", s_defPtr->doMomentumAdvection);

    s_defPtr->doTemperatureAdvection = true;
    pp.query("doTemperatureAdvection", s_defPtr->doTemperatureAdvection);

    s_defPtr->doSalinityAdvection = true;
    pp.query("doSalinityAdvection", s_defPtr->doSalinityAdvection);

    s_defPtr->doScalarAdvection = true;
    pp.query("doScalarAdvection", s_defPtr->doScalarAdvection);

    s_defPtr->doViscousForcing = true;
    pp.query("doViscousForcing", s_defPtr->doViscousForcing);

    s_defPtr->doTemperatureDiffusion = true;
    pp.query("doTemperatureDiffusion", s_defPtr->doTemperatureDiffusion);

    s_defPtr->doSalinityDiffusion = true;
    pp.query("doSalinityDiffusion", s_defPtr->doSalinityDiffusion);

    s_defPtr->doScalarDiffusion = true;
    pp.query("doScalarDiffusion", s_defPtr->doScalarDiffusion);

    s_defPtr->doGravityForcing = true;
    pp.query("doGravityForcing", s_defPtr->doGravityForcing);

    s_defPtr->doCoriolisForcing = pp.contains("coriolisF");
    pp.query("doCoriolisForcing", s_defPtr->doCoriolisForcing);
    if (s_defPtr->doCoriolisForcing) {
        if (s_defPtr->coriolisF.vectorLength() <= smallReal) {
            MAYDAYWARNING(
                "You set rhs.doCoriolisForcing = true but set rhs.coriolisF = "
                "all zeros. Are you sure you want this?");
        }
    }

    s_defPtr->doTidalForcing = pp.contains("tidalU0");
    s_defPtr->doTidalForcing |= pp.contains("tidalOmega");
    s_defPtr->doTidalForcing |= pp.contains("tidalInitPhase");
    pp.query("doTidalForcing", s_defPtr->doTidalForcing);

    s_defPtr->doSpongeForcing = pp.contains("spongeWidthLo");
    s_defPtr->doSpongeForcing |= pp.contains("spongeWidthHi");
    s_defPtr->doSpongeForcing |= pp.contains("spongeTimeCoeff");
    pp.query("doSpongeForcing", s_defPtr->doSpongeForcing);


    s_defPtr->doMomAdvRefluxing = false;
    pp.query("doMomAdvRefluxing", s_defPtr->doMomAdvRefluxing);

    s_defPtr->doTemperatureAdvRefluxing = false;
    pp.query("doTemperatureAdvRefluxing", s_defPtr->doTemperatureAdvRefluxing);

    s_defPtr->doSalinityAdvRefluxing = false;
    pp.query("doSalinityAdvRefluxing", s_defPtr->doSalinityAdvRefluxing);

    s_defPtr->doScalarAdvRefluxing = false;
    pp.query("doScalarAdvRefluxing", s_defPtr->doScalarAdvRefluxing);

    s_defPtr->doViscousRefluxing = false;
    pp.query("doViscousRefluxing", s_defPtr->doViscousRefluxing);

    s_defPtr->doTemperatureDiffusiveRefluxing = false;
    pp.query("doTemperatureDiffusiveRefluxing",
             s_defPtr->doTemperatureDiffusiveRefluxing);

    s_defPtr->doSalinityDiffusiveRefluxing = false;
    pp.query("doSalinityDiffusiveRefluxing",
             s_defPtr->doSalinityDiffusiveRefluxing);

    s_defPtr->doScalarDiffusiveRefluxing = false;
    pp.query("doScalarDiffusiveRefluxing",
             s_defPtr->doScalarDiffusiveRefluxing);


    s_defPtr->computeInitPressure = true;
    pp.query("computeInitPressure", s_defPtr->computeInitPressure);


    // Tidal parameters
    if (s_defPtr->doTidalForcing || tidalBCPresent) {
        vreal = Vector<Real>(SpaceDim, 0.0);
        pp.getarr("tidalU0", vreal, 0, SpaceDim);
        s_defPtr->tidalU0 = RealVect(vreal);

        vreal = Vector<Real>(SpaceDim, 0.0);
        pp.getarr("tidalOmega", vreal, 0, SpaceDim);
        s_defPtr->tidalOmega = RealVect(vreal);

        vreal = Vector<Real>(SpaceDim, 0.0);
        pp.getarr("tidalInitPhase", vreal, 0, SpaceDim);
        s_defPtr->tidalInitPhase = RealVect(vreal);

    } else {
        vreal = Vector<Real>(SpaceDim, 0.0);
        s_defPtr->tidalU0 = RealVect(vreal);

        vreal = Vector<Real>(SpaceDim, 0.0);
        s_defPtr->tidalOmega = RealVect(vreal);

        vreal = Vector<Real>(SpaceDim, 0.0);
        s_defPtr->tidalInitPhase = RealVect(vreal);
    }


    // Sponge parameters
    if (s_defPtr->doSpongeForcing) {
        vreal = Vector<Real>(SpaceDim, 0.0);
        pp.getarr("spongeWidthLo", vreal, 0, SpaceDim);
        for (int d = 0; d < SpaceDim; ++d) {
            s_defPtr->spongeWidth[d][0] = vreal[d];
        }

        vreal = Vector<Real>(SpaceDim, 0.0);
        pp.getarr("spongeWidthHi", vreal, 0, SpaceDim);
        for (int d = 0; d < SpaceDim; ++d) {
            s_defPtr->spongeWidth[d][1] = vreal[d];
        }

        s_defPtr->spongeTimeCoeff = 15.0;
        pp.get("spongeTimeCoeff", s_defPtr->spongeTimeCoeff);

    } else {
        for (int d = 0; d < SpaceDim; ++d) {
            s_defPtr->spongeWidth[d][0] = 0.0;
            s_defPtr->spongeWidth[d][1] = 0.0;
        }
        s_defPtr->spongeTimeCoeff = 1.0e300;
    }

    // Send defaults to pout.
    s_defPtr->dump();
}
