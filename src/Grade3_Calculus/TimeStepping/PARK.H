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
#ifndef ___PARK_H__INCLUDED___
#define ___PARK_H__INCLUDED___

#include "PARKRHS.H"
#include "PARKCoeffs.H" // the user will need this.
#include "RealVect.H"
#include "TimeParameters.H"


/// This is a stand-alone utility function. It is not a member of PARK.
/// It returns the order of an RK method and which test failed.
/// This is very helpful when creating new RK methods...it's very easy to
/// type just one number wrong in the Butcher table!
template <int numStages>
constexpr int
RKOrderChecks(const Real a[numStages][numStages],
              const Real b[numStages],
              const Real c[numStages],
              const bool a_quietMode = true);

template <int numStages>
constexpr int
RKStageOrderChecks(const Real a[numStages][numStages],
                   const Real bE[numStages],
                   const Real bI[numStages],
                   const Real c[numStages],
                   const bool a_quietMode = true);


/**
 * \class   PARK
 * \brief   A template for all Projected Additive Runge-Kutta methods.
 * \date    2019/12/06
 * \author  Edward.Santilli@jefferson.edu
 * \details
 *  A template for all Projected Additive Runge-Kutta methods.
 *  The structs that can be passed in as RKC are defined in ARKCoeffs.H.
 */
template <class RKC>
class PARK
{
public:
    /// Constructor
    PARK(const DisjointBoxLayout& a_grids,
         const int                a_velNumComps,
         const IntVect&           a_velGhostVect,
         const int                a_pNumComps,
         const IntVect&           a_pGhostVect,
         const int                a_qNumComps,
         const IntVect&           a_qGhostVect);

    /// Destructor
    ~PARK();

    /// Accessors
    Real oldTime() const {return m_oldTime;}
    Real newTime() const {return m_newTime;}
    Real dt() const {return m_newTime - m_oldTime;}
    bool hasEmbeddedScheme() const {return RKC::hasEmbeddedScheme;}

    /// The main timestepper.
    void
    advance(LevelData<FluxBox>&   a_vel,
            LevelData<FArrayBox>& a_p,
            LevelData<FArrayBox>& a_q,
            const Real            a_oldTime,
            const Real            a_dt,
            PARKRHS*              a_rhsPtr);

    /// A cheap Forward Euler timestepper.
    /// This is useful for generating rough estimates.
    void
    FEadvance(LevelData<FluxBox>&   a_vel,
              LevelData<FArrayBox>& a_p,
              LevelData<FArrayBox>& a_q,
              const Real            a_oldTime,
              const Real            a_dt,
              PARKRHS*              a_rhsPtr);

    inline constexpr Real
    ERKStabilityRe() const
    {
        return RKC::ERKStabilityRe;
    }

    inline constexpr Real
    ERKStabilityIm() const
    {
        return RKC::ERKStabilityIm;
    }

    /// @brief Computed the next dt based on error estimates.
    /// @param a_tol         Velocity error tolerance.
    /// @param a_useImplicit If false, we will consider the implicit part of the
    ///                      RK scheme to be unused.
    /// @return The new, limited dt.
    Real
    controllerDt(const Real a_tol,
                 const bool a_useImplicit,
                 const bool a_useElementary,
                 const bool a_usePI,
                 const bool a_usePID) const;

    /// Returns true if we can use the time interp functions.
    inline bool
    readyToInterp() const {return m_interpDataReady;}

    /// Interpolate vel in time.
    Real
    velTimeInterp(LevelData<FluxBox>& a_vel,
                  const Real          a_time,
                  int                 a_srcComp = 0,
                  int                 a_destComp = 0,
                  int                 a_numComp  = -1) const;

    /// Interpolate p in time.
    Real
    pTimeInterp(LevelData<FArrayBox>& a_p,
                const Real            a_time,
                int                   a_srcComp = 0,
                int                   a_destComp = 0,
                int                   a_numComp  = -1) const;

    /// Interpolate q in time.
    Real
    qTimeInterp(LevelData<FArrayBox>& a_q,
                const Real            a_time,
                int                   a_srcComp = 0,
                int                   a_destComp = 0,
                int                   a_numComp  = -1) const;

protected:
    /// Computes theta = (a_time - oldTime) / dt.
    /// If a_time is close to oldTime or newTime, this will snap to 0 or 1
    /// values to prevent extrapolation when not expected.
    Real
    computeTheta(const Real a_time) const;

    /// Interpolate data in time using a linear approx.
    /// This will use piecewise linear interpolation and return how far along the
    /// initial and final data the interpolation took place (0 = old time,
    /// 1 = new time). If this returns a number outside of [0,1], then we
    /// extrapolated and you need to worry about stability.
    Real
    velLinearInterp(LevelData<FluxBox>& a_vel,
                    const Real          a_time,
                    int                 a_srcComp = 0,
                    int                 a_destComp = 0,
                    int                 a_numComp  = -1) const;

    /// Interpolate data in time using a linear approx.
    /// This will use piecewise linear interpolation and return how far along the
    /// initial and final data the interpolation took place (0 = old time,
    /// 1 = new time). If this returns a number outside of [0,1], then we
    /// extrapolated and you need to worry about stability.
    Real
    pLinearInterp(LevelData<FArrayBox>& a_p,
                  const Real            a_time,
                  int                   a_srcComp = 0,
                  int                   a_destComp = 0,
                  int                   a_numComp  = -1) const;

    /// Interpolate data in time using a linear approx.
    /// This will use piecewise linear interpolation and return how far along the
    /// initial and final data the interpolation took place (0 = old time,
    /// 1 = new time). If this returns a number outside of [0,1], then we
    /// extrapolated and you need to worry about stability.
    Real
    qLinearInterp(LevelData<FArrayBox>& a_q,
                  const Real            a_time,
                  int                   a_srcComp = 0,
                  int                   a_destComp = 0,
                  int                   a_numComp  = -1) const;

    /// Interpolate data in time using the dense output formula.
    /// This will use the dense output formula and return how far along the
    /// initial and final data the interpolation took place (0 = old time,
    /// 1 = new time). If this returns a number outside of [0,1], then we
    /// extrapolated and you need to worry about stability.
    Real
    velDenseOutput(LevelData<FluxBox>& a_vel,
                   const Real          a_time,
                   int                 a_srcComp  = 0,
                   int                 a_destComp = 0,
                   int                 a_numComp  = -1) const;

    /// Interpolate data in time using the dense output formula.
    /// This will use the dense output formula and return how far along the
    /// initial and final data the interpolation took place (0 = old time,
    /// 1 = new time). If this returns a number outside of [0,1], then we
    /// extrapolated and you need to worry about stability.
    Real
    qDenseOutput(LevelData<FArrayBox>& a_q,
                 const Real            a_time,
                 int                   a_srcComp  = 0,
                 int                   a_destComp = 0,
                 int                   a_numComp  = -1) const;


    /// Interpolate data in time using a Lagrange polynomial.
    /// This will use a Lagrange polynomial to interpolate between initial and
    /// final data, returning how far along the timestep the interpolated state
    /// is (0 = old time, 1 = new time). If this returns a number outside of
    /// [0,1], then we extrapolated and you need to worry about stability.
    Real
    velLagrangeInterp(LevelData<FluxBox>& a_vel,
                      const Real          a_time,
                      int                 a_srcComp  = 0,
                      int                 a_destComp = 0,
                      int                 a_numComp  = -1) const;

    /// Interpolate data in time using a Lagrange polynomial.
    /// This will use a Lagrange polynomial to interpolate between initial and
    /// final data, returning how far along the timestep the interpolated state
    /// is (0 = old time, 1 = new time). If this returns a number outside of
    /// [0,1], then we extrapolated and you need to worry about stability.
    Real
    qLagrangeInterp(LevelData<FArrayBox>& a_q,
                    const Real            a_time,
                    int                   a_srcComp  = 0,
                    int                   a_destComp = 0,
                    int                   a_numComp  = -1) const;

    /// Interpolate data in time using a Lagrange polynomial.
    /// This will use a Lagrange polynomial to interpolate between initial and
    /// final data, returning how far along the timestep the interpolated state
    /// is (0 = old time, 1 = new time). If this returns a number outside of
    /// [0,1], then we extrapolated and you need to worry about stability.
    Real
    pLagrangeInterp(LevelData<FArrayBox>& a_p,
                    const Real            a_time,
                    int                   a_srcComp  = 0,
                    int                   a_destComp = 0,
                    int                   a_numComp  = -1) const;

    /// Resets the data. The interp functions will be able to copy vel, p, and q
    /// if called with a_time = a_oldTime.
    void
    setOldQ(const LevelData<FluxBox>&   a_velOld,
            const LevelData<FArrayBox>& a_pOld,
            const LevelData<FArrayBox>& a_qOld,
            const Real                  a_oldTime);

    /// \name Static utilities
    /// \{

    /// Sets x = y on all comps.
    /// x and y must have the same number of comps.
    static void
    copy(LevelData<FluxBox>&       a_x,
         const LevelData<FluxBox>& a_y);

    /// Sets x = y on all comps.
    /// x and y must have the same number of comps.
    static void
    copy(LevelData<FArrayBox>&       a_x,
         const LevelData<FArrayBox>& a_y);

    /// Sets x += b*y on all comps.
    /// x and y must have the same number of comps.
    static void
    plus(LevelData<FluxBox>&       a_x,
         const Real                a_b,
         const LevelData<FluxBox>& a_y);

    /// Sets x += b*y on all comps.
    /// x and y must have the same number of comps.
    static void
    plus(LevelData<FArrayBox>&       a_x,
         const Real                  a_b,
         const LevelData<FArrayBox>& a_y);

    /// Sets x = a*x + b*y on all comps.
    /// x and y must have the same number of comps.
    static void
    axby(const Real                a_a,
         LevelData<FluxBox>&       a_x,
         const Real                a_b,
         const LevelData<FluxBox>& a_y);

    /// Sets x = a*x + b*y on all comps.
    /// x and y must have the same number of comps.
    static void
    axby(const Real                  a_a,
         LevelData<FArrayBox>&       a_x,
         const Real                  a_b,
         const LevelData<FArrayBox>& a_y);
    /// \}

    DisjointBoxLayout m_grids;
    const int         m_velNumComps;
    const IntVect     m_velGhostVect;
    const int         m_pNumComps;
    const IntVect     m_pGhostVect;
    const int         m_qNumComps;
    const IntVect     m_qGhostVect;

    int  m_prevIdxForProj[RKC::numStages + 1];
    Real m_zetaE[RKC::numStages + 1][RKC::numStages];
    Real m_zetaI[RKC::numStages + 1][RKC::numStages];
    bool m_needsAssembly;
    bool m_needsAssemblyProjection;

// #define PARK__USE_INCREMENTAL_RK
    static const bool s_extrapStagePressure    = true;
    static const bool s_extrapAssemblyPressure = true;
    static const bool s_removeDivFromK         = false;
    static const bool s_useDenseOutput         = false;

    // The stage data.
    LevelData<FluxBox>   m_vel[RKC::numStages];
    LevelData<FluxBox>   m_kvelE[RKC::numStages];
    LevelData<FluxBox>   m_kvelI[RKC::numStages];
    LevelData<FArrayBox> m_p[RKC::numStages];
    LevelData<FArrayBox> m_q[RKC::numStages];
    LevelData<FArrayBox> m_kqE[RKC::numStages];
    LevelData<FArrayBox> m_kqI[RKC::numStages];

    Real m_oldTime;
    Real m_newTime;
    Real m_dt;

    bool m_Q0Ready;
    bool m_interpDataReady;


    /// \name Error controller functions and data
    /// \{
    struct ErrorHistoryElem
    {
        ErrorHistoryElem()
        : hasValue(false), time(quietNAN), dt(quietNAN), localError(quietNAN) {}

        ErrorHistoryElem(const Real a_time, const Real a_dt, const Real a_localError)
        : hasValue(true), time(a_time), dt(a_dt), localError(a_localError) {}

        bool hasValue;
        Real time;
        Real dt;
        Real localError;
    };
    mutable std::array<ErrorHistoryElem, 3> m_errorHistory; // array idx 0 is newest.
    mutable std::array<Real, 3>             m_controllerHistory; // For limiting.

    /// @brief Pushes a new element to a stack, throwing away the oldest element.
    /// @tparam E          The type of the new element.
    /// @tparam N          The number of elements on the fixed-size stack.
    /// @param a_newElem   The new element to put in a_hsitory[0].
    /// @param a_history   The stack.
    /// @param a_isNewTime If this is false, a_history[0] will be overwritten and
    ///                    nothing else will change -- there will be no "push."
    template<typename E, size_t N>
    static void
    pushToHistory(const E&          a_newElem,
                  std::array<E, N>& a_history,
                  const bool        a_isNewTime);

    /// @brief The controller limiter of Doderlind and Wang,
    /// doi:10.1016/j.cam.2005.03.008, eq. 10.
    /// @param a_proposedController        The controller to be limited.
    /// @param a_lastController            The controller that was used in the last timestep.
    /// @param a_maxGrowth                 Roughly the maximum amount of controller growth/shrinkage.
    /// @return The new, limited controller.
    static Real
    softLimiter(Real       a_proposedController,
                const Real a_lastController,
                const Real a_maxGrowth = 1.0);

    /// @brief A simple controller limiter restricts the growth/shrinkage to be
    /// within [a_maxShrinkage, a_maxGrowth].
    /// @param a_proposedController The controller to be limited.
    /// @param a_lastController     The controller that was used in the last timestep.
    /// @param a_maxShrinkage
    /// @param a_maxGrowth
    /// @return The new, limited controller.
    static Real
    hardLimiter(Real       a_proposedController,
                const Real a_lastController,
                const Real a_maxShrinkage = 0.5,
                const Real a_maxGrowth    = 2.0);

    /// @brief
    /// @param a_proposedController  The controller to be limited.
    /// @param a_oldController       The controller that was used in the last timestep.
    /// @param a_reallyOldController The controller that was used in two timesteps ago.
    /// @param a_dampingFactor       <= 0 produces no damping. 1 = very strong damping.
    /// @return The new, limited controller.
    static Real
    oscillationLimiter(Real       a_proposedController,
                       const Real a_oldController,
                       const Real a_reallyOldController,
                       const Real a_dampingFactor = 0.5);
    /// \}

    static const int s_orderE = RKOrderChecks<RKC::numStages>(RKC::aE, RKC::bE, RKC::c);
    static const int s_orderI = RKOrderChecks<RKC::numStages>(RKC::aI, RKC::bI, RKC::c);
    static const int s_orderEmbeddedE = RKOrderChecks<RKC::numStages>(RKC::aE, RKC::bEhat, RKC::c);
    static const int s_orderEmbeddedI = RKOrderChecks<RKC::numStages>(RKC::aI, RKC::bIhat, RKC::c);
    static const int s_stageOrderE = RKStageOrderChecks<RKC::numStages>(RKC::aE, RKC::bE, RKC::bI, RKC::c);
    static const int s_stageOrderI = RKStageOrderChecks<RKC::numStages>(RKC::aI, RKC::bE, RKC::bI, RKC::c);
};


// Include templated definitions
#define M0091158a93725eb57192ee021f0139bf
#   include "PARKI.H"
#undef M0091158a93725eb57192ee021f0139bf

#endif  //!___PARK_H__INCLUDED___
