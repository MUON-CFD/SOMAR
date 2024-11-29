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
 *  https://github.com/UNC-CFD/LES-SOMAR.
 ******************************************************************************/
#ifndef ___SDC_H__INCLUDED___
#define ___SDC_H__INCLUDED___

#include "PARKRHS.H"

// Make generic:
//  Low-order base integrator. (Can use PARK.)
//  Quadrature scheme (nodes and integrator)
//  Interpolation polynomial  ...   use cmath's legendre and consider long double version


template <size_t NumNodes>
class LobattoQuadrature
{
public:
    /// Will be on the interval [a, b].
    inline static Real getNode(const size_t n, const Real a, const Real b) {
        CH_assert(n < NumNodes);
        return 0.5 * (b - a) * s_nodes[n] + 0.5 * (b + a);
    }

    inline static Real getWeight(const size_t n, const Real a, const Real b) {
        CH_assert(n < NumNodes);
        return 0.5 * (b - a) * s_weights[n];
        // if (a_n == 0 || a_n == NumNodes - 1) {
        //     return 2.0 / (a_n * (a_n + 1.0));
        // } else {
        //     const Real l = std::legendre(a_n, getNode(a_n));
        //     return 2.0 / (a_n * (a_n + 1.0) * l * l);
        // }
    }

private:
    static const std::array<Real, NumNodes> s_nodes;
    static const std::array<Real, NumNodes> s_weights;
};






class SDC
{
public:
    static constexpr size_t NumNodes = 4;
    using Quadrature = LobattoQuadrature<NumNodes>;

    /// Constructor
    SDC(const DisjointBoxLayout& a_grids,
        const int                a_velNumComps,
        const IntVect&           a_velGhostVect,
        const int                a_pNumComps,
        const IntVect&           a_pGhostVect,
        const int                a_qNumComps,
        const IntVect&           a_qGhostVect);

    /// Destructor
    ~SDC() {}

    /// Accessors
    Real oldTime() const {return m_oldTime;}
    Real newTime() const {return m_newTime;}
    Real dt() const {return m_newTime - m_oldTime;}

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
        return 1.0;
    }

    inline constexpr Real
    ERKStabilityIm() const
    {
        return 1.0;
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
                 const bool a_usePID) const
    {return 0.0;}


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
    /// @brief Advanced the solution. No error checks. No setup.
    /// @param a_vel
    /// @param a_p
    /// @param a_q
    /// @param a_kvelE
    /// @param a_kqE
    /// @param a_kvelI
    /// @param a_kqI
    /// @param a_vel0
    /// @param a_p0
    /// @param a_q0
    /// @param a_time0
    /// @param a_dt
    /// @param a_refluxDt
    /// @param a_rhsPtr
    static void
    basicAdvance(LevelData<FluxBox>&   a_vel,
                 LevelData<FArrayBox>& a_p,
                 LevelData<FArrayBox>& a_q,
                 LevelData<FluxBox>&   a_kvelE,
                 LevelData<FArrayBox>& a_kqE,
                 LevelData<FluxBox>&   a_kvelI,
                 LevelData<FArrayBox>& a_kqI,
                 LevelData<FluxBox>&   a_vel0,
                 LevelData<FArrayBox>& a_p0,
                 LevelData<FArrayBox>& a_q0,
                 const Real            a_time0,
                 const Real            a_dt,
                 const Real            a_refluxDt,
                 PARKRHS*              a_rhsPtr);

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

    /// Resets the data. The interp functions will be able to copy vel, p, and q
    /// if called with a_time = a_oldTime.
    void
    setOldQ(const LevelData<FluxBox>&   a_velOld,
            const LevelData<FArrayBox>& a_pOld,
            const LevelData<FArrayBox>& a_qOld,
            const Real                  a_oldTime);

    // /// @brief Sets m_dt, m_newTime, m_quadNodes, m_quadWeights.
    // ///  Must call this after calling setOldQ.
    // /// @param a_dt
    // void
    // setQuadratureNodes(const Real a_dt);


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

    std::array<LevelData<FluxBox>, NumNodes>   m_vel;
    std::array<LevelData<FArrayBox>, NumNodes> m_p;
    std::array<LevelData<FArrayBox>, NumNodes> m_q;

    std::array<LevelData<FluxBox>, NumNodes>   m_kvelE;
    std::array<LevelData<FArrayBox>, NumNodes> m_kqE;
    std::array<LevelData<FluxBox>, NumNodes>   m_kvelI;
    std::array<LevelData<FArrayBox>, NumNodes> m_kqI;

    std::array<LevelData<FluxBox>, NumNodes>   m_kvelOld;
    std::array<LevelData<FArrayBox>, NumNodes> m_kqOld;

    Real m_oldTime;
    Real m_newTime;
    Real m_dt;

    std::array<std::array<Real, NumNodes>, NumNodes - 1> m_correctionWeights;

    bool m_Q0Ready;
    bool m_finalForcesReady;
    bool m_interpDataReady;
};


#endif //!___SDC_H__INCLUDED___
