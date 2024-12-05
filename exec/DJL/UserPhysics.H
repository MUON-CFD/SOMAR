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
#ifndef ___UserPhysics_H__INCLUDED___
#define ___UserPhysics_H__INCLUDED___

#include "AMRNSLevel.H"


class UserPhysics: public AMRNSLevel
{
public:
    // Constructor
    UserPhysics ();

    // Virtual destructor for good measure.
    virtual ~UserPhysics ();


    // Custom scalars ----------------------------------------------------------
    // Custom scalars will automatically be advected and, if requested,
    // diffused. If your needs are more specialized than that, then it
    // is up to you to define a forcing function (see below).

    // This tells the Navier-Stokes class how many additional scalars
    // to create. Do not include T, S, or eddyNu. This is inline so that it is
    // just as fast as reading from a const int.
    virtual int
    numScalars() const
    {
        return 0;
    }

    // Returns the names of the scalars. This will be used to put scalar names
    // in the plot files.
    virtual std::string
    getScalarName(const int a_comp) const;


    // ICs ---------------------------------------------------------------------
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
    virtual void
    setICs(State& a_state);


    // Other -------------------------------------------------------------------

    /// Our custom EoS. This simply sets b = T.
    virtual void
    equationOfState(FArrayBox&       a_bFAB,
                    const FArrayBox& a_TFAB,
                    const FArrayBox& a_SFAB,
                    const FArrayBox& a_zFAB) const override
    {
        a_bFAB.copy(a_TFAB);
    }

    /// Tag cells during runtime.
    virtual void
    tagCells(IntVectSet &a_tags);


    virtual void
    computeVelSoln(LevelData<FluxBox>& a_vel,
                   const Real          a_time) const;

    virtual void
    computeBSoln(LevelData<FArrayBox>& a_b,
                 const Real            a_time) const;

    virtual void
    postStep(const LevelData<FluxBox>&   a_vel  __nowarn_unused,
             const LevelData<FArrayBox>& a_p    __nowarn_unused,
             const LevelData<FArrayBox>& a_q    __nowarn_unused,
             const Real                  a_time __nowarn_unused) const
    {
        // Do nothing for now.
    }

protected:
    static void
    setStaticMembers();

    static std::vector<Real>
    getLayer_eta(const Real a_z);

    static std::vector<Real>
    getLayer_eta_z(const Real a_z);

    static inline Real
    extrusionEnvelope(const Real a_yp)
    {
        return 0.5 * (tanh(s_envSlope * (a_yp + 0.5 * s_envWidth)) -
                      tanh(s_envSlope * (a_yp - 0.5 * s_envWidth)));
    }

#if CH_SPACEDIM == 2
    static constexpr Real s_rotAngle    = 90.0 * Pi / 180.0;
    static constexpr Real s_xOffset     = 1.0 * 181.0;
    static constexpr Real s_yOffset     = 0.0 * 96.0;
    static constexpr Real s_envSlope    = 1.0 / 6.4;
    static constexpr Real s_envWidth    = 19200.0;
    static constexpr Real s_bProjTagTol = 0.01;
    static constexpr Real s_pTagTol     = 0.05;
#else
    static constexpr Real s_rotAngle    = 45.0 * Pi / 180.0;
    static constexpr Real s_xOffset     = 128.0;
    static constexpr Real s_yOffset     = 128.0;
    static constexpr Real s_envSlope    = 1.0 / 6.4;
    static constexpr Real s_envWidth    = 192.0;
    static constexpr Real s_bProjTagTol = 0.1;//0.05;
    static constexpr Real s_pTagTol     = 0.05;
#endif

    struct SrcData {
        SrcData(): L(-1.0) {}
        bool isDefined() const {return L > 0.0;}

        Real              L, H, A, c;
        size_t            nx, nz;
        std::vector<Real> x, z;
        FArrayBox         etaFAB;
    };
    static SrcData s_src;
};


#endif //!___UserPhysics_H__INCLUDED___
