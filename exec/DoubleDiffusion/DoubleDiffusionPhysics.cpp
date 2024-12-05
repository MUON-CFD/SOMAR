/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2018
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
#include "DoubleDiffusionPhysics.H"
#include "SetValLevel.H"
#include "Debug.H"
#include "ParmParse.H"
#include <random>
#include "AMRNSLevelF_F.H"
#include "MiscUtils.H"


//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
DoubleDiffusionPhysics::DoubleDiffusionPhysics()
: AMRNSLevel::AMRNSLevel()
{
}


//------------------------------------------------------------------------------
// Virtual destructor for good measure.
//------------------------------------------------------------------------------
DoubleDiffusionPhysics::~DoubleDiffusionPhysics()
{
}


// -----------------------------------------------------------------------------
// Returns the names of the scalars. This will be used to put scalar names
// in the plot files.
// -----------------------------------------------------------------------------
std::string
DoubleDiffusionPhysics::getScalarName(const int a_comp) const
{
    switch (a_comp) {
        // If you have scalars, add cases here.
        default:
            MayDay::Error(
                "DoubleDiffusionPhysics::getScalarName: "
                "a_comp out of range.");
    }

    return "Undefined";
}


//------------------------------------------------------------------------------
// Set the ICs on all fields at once. This includes the velocity,
// total temperature, total salinity, eddy viscosity and diffusivity, all
// user-defined scalars, and pressure.
//
// By default, everything is set to zero except T and S which are set to
// thier background stratification values (T' = S' = 0).
//
// Also, you do not need to specify the pressure. If rhs.computeInitPressure
// = 1, SOMAR will automatically figure out the pressure. However, if you
// have a good guess, SOMAR can use that as an initial guess in the Poisson
// solver.
//
// By the time this function is called, *m_levGeoPtr will be ready for use.
//------------------------------------------------------------------------------
void
DoubleDiffusionPhysics::setICs(State& a_state)
{
    const ProblemDomain&     domain = m_levGeoPtr->getDomain();
    const DisjointBoxLayout& grids  = m_levGeoPtr->getBoxes();
    DataIterator             dit    = grids.dataIterator();

    ParmParse pp("strat");

    Real Tmin, Tmax;
    pp.get("Tmin", Tmin);
    pp.get("Tmax", Tmax);

    Real Smin, Smax;
    pp.get("Smin", Smin);
    pp.get("Smax", Smax);

    if constexpr (SpaceDim == 2) {
        Real z0, deltaz, pertA;
        Real pertK1, pertK2, pertK3;
        pp.get("z0", z0);
        pp.get("deltaz", deltaz);
        pp.get("pertA", pertA);
        pp.get("pertK", pertK2);
        pertK1 = pertK2 - 12.0;
        pertK3 = pertK2 + 12.0;
        pertK1 *= 2.0 * Pi / m_levGeoPtr->getDomainLength(0);
        pertK2 *= 2.0 * Pi / m_levGeoPtr->getDomainLength(0);
        pertK3 *= 2.0 * Pi / m_levGeoPtr->getDomainLength(0);

        const Real avgT    = 0.5 * (Tmax + Tmin);
        const Real deltaT  = 0.5 * (Tmax - Tmin);

        const Real avgS    = 0.5 * (Smax + Smin);
        const Real deltaS  = 0.5 * (Smax - Smin);

        const Real zScale  = 1.0 / deltaz;

        std::vector<Real> randos(domain.size(0) + 2);
        generateRandomNumbers(randos);
        const int ivx0 = domain.domainBox().smallEnd(0) - 1;


        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox& TFAB = a_state.T[dit];
            FArrayBox& SFAB = a_state.S[dit];
            const Box& region = TFAB.box();

            FArrayBox xFAB(region, 2);
            m_levGeoPtr->fill_physCoor(xFAB, 0, 0);
            m_levGeoPtr->fill_physCoor(xFAB, 1, SpaceDim - 1);

            for (BoxIterator bit(region); bit.ok(); ++bit) {
                const IntVect& iv = bit();
                // const Real x = xFAB(iv, 0);
                // const Real z = xFAB(iv, 1) - pertA * sin(pertK2 * x);

                // const Real z = xFAB(iv, 1) - pertA * sin(pertK1 * x) * 0.5
                //                            - pertA * sin(pertK2 * x)
                //                            - pertA * sin(pertK3 * x) * 0.5;

                Real r = randos[iv[0] - ivx0];
                r = 2.0 * pertA * (r - 0.5);
                const Real z = xFAB(iv, 1) + r;

                const Real tanhArg = tanh(zScale * (z - z0));
                TFAB(iv) = avgT + deltaT * tanhArg;
                SFAB(iv) = avgS + deltaS * tanhArg;
            }
        }


        // Parameters defined in Singh and Srinivasan, Phys. Fluids 26, 062104 (2014).
        const ProblemContext* ctx = ProblemContext::getInstance();
        const Real nu     = ctx->rhs.nu;
        const Real kappaT = ctx->rhs.TKappa;
        const Real kappaS = ctx->rhs.SKappa;
        const Real H      = ctx->base.L[SpaceDim - 1];
        const Real dT     = Tmax - Tmin;
        const Real dS     = Smax - Smin;

        // Initital density stability ratio (alpha = betaT)
        const Real Rp0 = s_alphag * dT / (s_betag * dS);
        // Thermal Rayleigh number
        const Real RaT = s_alphag * dT * H * H * H / (nu * kappaT);
        // Salinity Rayleigh number
        const Real RaS = RaT / Rp0;
        // Timescale ratio
        const Real tau = kappaT / kappaS;

        POUT(Rp0);
        POUT(RaT);
        POUT(RaS);
        POUT(tau);
        pout() << "Stable IC if 1.0 << " << Rp0 << " << " << std::pow(tau, -3.0/2.0) << std::endl;

    } else {
        Real z0, deltaz;

        constexpr size_t numPerts = 3;
        std::vector<Real> pertAx(numPerts, 0.0);
        std::vector<Real> pertKx(numPerts, 0.0);
        std::vector<Real> pertAy(numPerts, 0.0);
        std::vector<Real> pertKy(numPerts, 0.0);

        pp.get("z0", z0);
        pp.get("deltaz", deltaz);
        pp.getarr("pertAx", pertAx, 0, numPerts);
        pp.getarr("pertKx", pertKx, 0, numPerts);
        pp.getarr("pertAy", pertAy, 0, numPerts);
        pp.getarr("pertKy", pertKy, 0, numPerts);

        for (size_t idx = 0; idx < numPerts; ++idx) {
            pertKx[idx] *= 2.0 * Pi / m_levGeoPtr->getDomainLength(0);
            pertKy[idx] *= 2.0 * Pi / m_levGeoPtr->getDomainLength(1);
        }

        std::mt19937       gen(procID());  // seed the generator
        std::uniform_real_distribution<> distr(-pertAx[1], pertAx[1]);

        const Real avgT    = 0.5 * (Tmax + Tmin);
        const Real deltaT  = 0.5 * (Tmax - Tmin);

        const Real avgS    = 0.5 * (Smax + Smin);
        const Real deltaS  = 0.5 * (Smax - Smin);

        const Real zScale  = 1.0 / deltaz;

        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox& TFAB = a_state.T[dit];
            FArrayBox& SFAB = a_state.S[dit];
            const Box& region = TFAB.box();

            FArrayBox xFAB(region, SpaceDim);
            m_levGeoPtr->fill_physCoor(xFAB);

            for (BoxIterator bit(region); bit.ok(); ++bit) {
                const IntVect& iv = bit();
                // const Real x = xFAB(iv, 0);
                // const Real y = xFAB(iv, 1);
                // z = xFAB(iv, 2) - pertAx[0] * sin(pertKx[0] * x)
                //                 - pertAx[1] * sin(pertKx[1] * x)
                //                 - pertAx[2] * sin(pertKx[2] * x)
                //                 - pertAy[0] * sin(pertKy[0] * y)
                //                 - pertAy[1] * sin(pertKy[1] * y)
                //                 - pertAy[2] * sin(pertKy[2] * y);

                const Real z = xFAB(iv, 2) + distr(gen);

                const Real tanhArg = tanh(zScale * (z - z0));
                TFAB(iv) = avgT + deltaT * tanhArg;
                SFAB(iv) = avgS + deltaS * tanhArg;
            }
        }

    }
}


// -----------------------------------------------------------------------------
// The equation of state.
// Fills a data holder with the buoyancy given temperature and salinity.
//
// The FABs will be at this level's index space, but may be defined over a
// patch that is not in this level's DisjointBoxLayout. Therefore, you
// can compare the FAB's boxes to this level's ProblemDomain to figure out
// where the ghosts are, but you should not loop over the grids to find
// an associated DataIndex. For any reasonable EoS, this should not be
// a problem.
//
// T and S must be the total temperature and salinity, not just the deviation
// from the background. Likewise, this function returns the total buoyancy.
// -----------------------------------------------------------------------------
void
DoubleDiffusionPhysics::equationOfState(FArrayBox&       a_bFAB,
                                        const FArrayBox& a_TFAB,
                                        const FArrayBox& a_SFAB,
                                        const FArrayBox& a_zFAB) const
{
    const Box& region = a_bFAB.box();

    FORT_DEFAULTLINEAREOS(
        CHF_FRA1(a_bFAB, 0),
        CHF_CONST_FRA1(a_TFAB, 0),
        CHF_CONST_FRA1(a_SFAB, 0),
        CHF_BOX(region),
        CHF_CONST_REAL(s_alphag),
        CHF_CONST_REAL(s_betag));
}


// // -----------------------------------------------------------------------------
// // Sets the T and S background stratification at specified heights.
// // Provide as many elements as you need to resolve the stratification.
// // This function will be called just once on level 0.
// // a_zmin and a_zmax are provided so that you know where the extrapolation
// // will take place. Keep in mind this range may include some space for
// // ghost values.
// // -----------------------------------------------------------------------------
// void
// DoubleDiffusionPhysics::setStratification(Vector<Real>& a_vTbar,
//                                           Vector<Real>& a_vSbar,
//                                           Vector<Real>& a_vz,
//                                           const Real    a_zmin,
//                                           const Real    a_zmax) const
// {
//     const ProblemContext* ctx = ProblemContext::getInstance();
//     const Real L = a_zmax - a_zmin;

//     // If no stratification, set to 0 and exit.
//     if (ctx->rhs.stratType == RHSParameters::StratType::NONE) {
//         // No stratification.
//         // Setting min and max values to 0 should will interp to 0 everywhere.
//         a_vz.resize(2);
//         a_vTbar.resize(2);
//         a_vSbar.resize(2);

//         a_vz[0] = a_zmin;
//         a_vTbar[0] = 0.0;
//         a_vSbar[0] = 0.0;

//         a_vz[1] = a_zmax;
//         a_vTbar[1] = 0.0;
//         a_vSbar[1] = 0.0;

//     } else if (ctx->rhs.stratType == RHSParameters::StratType::LINEAR) {
//         // Linear stratification
//         a_vz.resize(2);
//         a_vTbar.resize(2);
//         a_vSbar.resize(2);

//         a_vz[0] = a_zmin;
//         a_vTbar[0] = ctx->rhs.T0;
//         a_vSbar[0] = ctx->rhs.S0;

//         a_vz[1] = a_zmax;
//         a_vTbar[1] = ctx->rhs.T0 + L * ctx->rhs.dTdz;
//         a_vSbar[1] = ctx->rhs.S0 + L * ctx->rhs.dSdz;

//     } else if (ctx->rhs.stratType == RHSParameters::StratType::QUADRATIC) {
//         MAYDAYERROR("StratType::QUADRATIC not complete");

//     } else if (ctx->rhs.stratType == RHSParameters::StratType::TANH) {
//         MAYDAYERROR("StratType::TANH not complete");

//     } else if (ctx->rhs.stratType == RHSParameters::StratType::USER) {
//         MAYDAYERROR(
//             "rhs.statType = -1 means user-defined stratification. In this "
//             "case, you must override AMRNSLevel::setStratification in your "
//             "custom physics class.");

//     } else {
//         MAYDAYERROR(
//             "rhs.stratType = "
//             << ctx->rhs.stratType << " not recognized. "
//             << "Choose a pre-defined type or override "
//                "AMRNSLevel::setStratification in your custom physics class.");
//     }
// }
