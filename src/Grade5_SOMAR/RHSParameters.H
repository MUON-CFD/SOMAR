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
 *  https://github.com/somarhub.
 ******************************************************************************/
#ifndef ___RHSParameters_H__INCLUDED___
#define ___RHSParameters_H__INCLUDED___

#include <memory>
#include "RealVect.H"
#include "BaseParameters.H"


class RHSParameters
{
public:
    // Default constructor
    RHSParameters(const BaseParameters& a_baseParams,
                  const int             a_numLevels);

    Real         nu;
    Real         TKappa;
    Real         SKappa;
    // scalarsKappa is private. See below.

    bool doImplicitDiffusion;

    RealVect     coriolisF;

    // Not the same as EllipticOpBC::BCType! Perhaps they should be
    // consolidated somehow. *BCType has 4 values in 2D and
    // 6 values in 3D, one for each dir/side pair.
    struct VelBCType {
        enum {
            UNDEFINED = -1, // Can't imagine why we would need this.
            PERIODIC  =  0, // Automatically set when isPeriodic is true.
            NO_SLIP   =  1,
            FREE_SLIP =  2,
            CUSTOM    =  3, // Must be overridden.
            TIDAL     =  4,
            OUTFLOW   =  5,
            _NUM_BCTYPES
        };
    };
    int velBCTypeLo[CH_SPACEDIM];
    int velBCTypeHi[CH_SPACEDIM];
    int velBCType[CH_SPACEDIM][2]; // For convenience

    RealVect tidalU0;
    RealVect tidalOmega;
    RealVect tidalInitPhase; // 0 for sin...0.5*Pi for cos.

    struct TempBCType {
        enum {
            UNDEFINED = -1,  // Can't imagine why we would need this.
            PERIODIC  = 0,   // Automatically set when isPeriodic is true.
            ZERO_NEUM = 1,
            ZERO_NEUM_ON_PERT = 2,
            CUSTOM            = 3,  // Must be overridden.
            TIDAL             = 4,
            OUTFLOW           = 5,
            _NUM_BCTYPES
        };
    };
    int tempBCTypeLo[CH_SPACEDIM];
    int tempBCTypeHi[CH_SPACEDIM];
    int tempBCType[CH_SPACEDIM][2]; // For convenience

    struct SalinityBCType {
        enum {
            UNDEFINED = -1,  // Can't imagine why we would need this.
            PERIODIC  = 0,   // Automatically set when isPeriodic is true.
            ZERO_NEUM = 1,
            ZERO_NEUM_ON_PERT = 2,
            CUSTOM            = 3,  // Must be overridden.
            TIDAL             = 4,
            OUTFLOW           = 5,
            _NUM_BCTYPES
        };
    };
    int salinityBCTypeLo[CH_SPACEDIM];
    int salinityBCTypeHi[CH_SPACEDIM];
    int salinityBCType[CH_SPACEDIM][2]; // For convenience

    // Sponge stuff. Set in input file via spongeWidthLo and spongeWidthHi.
    // bool useSponge;
    Real spongeWidth[CH_SPACEDIM][2];
    Real spongeTimeCoeff;

    struct EddyViscMethods {
        enum {
            AVG_DOWN = -1,
            SET_TO_ZERO = 0,
            STATIC_SMAG,
            DYNAMIC_SMAG,
            DUCROS,
            _NUM_EDDYVISCMETHODS
        };
    };
    Vector<int> eddyViscMethod; // One elem for each level.
    Real eddyPrandtlT;
    Real eddyPrandtlS;
    // eddyPrandtlScalars is private. See below.
    RealVect eddyScale;         //(1,1,0) = 4 pt method, (1,1,1) = 6 pt method

    // Advection term settings.
    enum class Reconstruction {
        SecondOrder,
        FourthOrder,
    };
    Reconstruction velReconstruction;
    Reconstruction scalReconstruction;
    Real momAdvSkewness; // 0 = conservative, 1 = advective

    bool doMomentumAdvection;
    bool doTemperatureAdvection;
    bool doSalinityAdvection;
    bool doScalarAdvection;
    bool doViscousForcing;
    bool doTemperatureDiffusion;
    bool doSalinityDiffusion;
    bool doScalarDiffusion;
    bool doGravityForcing;
    bool doCoriolisForcing;
    bool doTidalForcing;
    bool doSpongeForcing;

    bool doMomAdvRefluxing;
    bool doTemperatureAdvRefluxing;
    bool doSalinityAdvRefluxing;
    bool doScalarAdvRefluxing;
    bool doViscousRefluxing;
    bool doTemperatureDiffusiveRefluxing;
    bool doSalinityDiffusiveRefluxing;
    bool doScalarDiffusiveRefluxing;

    bool computeInitPressure;

    // You shouldn't need to call this. AnisotropicAMR will do it for you.
    static void
    freeMemory();

    // It's nice to be able to see these parameters in pout.*.
    void
    dump() const;

    // This retrieves the Vectors and makes sure they are the correct size.
    Real
    getScalarsKappa(const unsigned int a_comp) const;

    Real
    getEddyPrandtlScalars(const unsigned int a_comp) const;

private:
    static void
    createDefaults(const BaseParameters& a_baseParams,
                   const int             a_numLevels);

    static std::unique_ptr<RHSParameters> s_defPtr;
    static bool                           s_constructorLock;

    Vector<Real> scalarsKappa;
    Vector<Real> eddyPrandtlScalars;
};


#endif //!___RHSParameters_H__INCLUDED___
