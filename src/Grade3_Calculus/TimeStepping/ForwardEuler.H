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
#ifndef ___ForwardEuler_H__INCLUDED___
#define ___ForwardEuler_H__INCLUDED___

#include "PARKRHS.H"


class ForwardEuler
{
public:
    // Constructor
    ForwardEuler(const DisjointBoxLayout& a_grids,
                 const int                a_velNumComps,
                 const IntVect&           a_velGhostVect,
                 const int                a_pNumComps,
                 const IntVect&           a_pGhostVect,
                 const int                a_qNumComps,
                 const IntVect&           a_qGhostVect);

    // Destructor
    ~ForwardEuler();

    // The main timestepper.
    void
    advance(LevelData<FluxBox>&   a_vel,
            LevelData<FArrayBox>& a_p,
            LevelData<FArrayBox>& a_q,
            const Real            a_oldTime,
            const Real            a_dt,
            PARKRHS*              a_rhsPtr);

    //
    inline bool
    readyToInterp()
    {
        return m_interpDataReady;
    }

    // Interpolate data in time using a linear approx.
    // This will use piecewise linear interpolation and return how far along the
    // initial and final data the interpolation took place (0 = old time,
    // 1 = new time). If this returns a number outside of [0,1], then we
    // extrapolated and you need to worry about stability.
    Real
    velTimeInterp(LevelData<FluxBox>& a_vel,
                  const Real          a_time,
                  int                 a_srcComp  = 0,
                  int                 a_destComp = 0,
                  int                 a_numComp  = -1) const;

    // p version.
    Real
    pTimeInterp(LevelData<FArrayBox>& a_p,
                const Real            a_time,
                int                   a_srcComp  = 0,
                int                   a_destComp = 0,
                int                   a_numComp  = -1) const;

    // q version.
    Real
    qTimeInterp(LevelData<FArrayBox>& a_q,
                const Real            a_time,
                int                   a_srcComp  = 0,
                int                   a_destComp = 0,
                int                   a_numComp  = -1) const;

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

    void
    FEadvance(LevelData<FluxBox>&   a_vel,
              LevelData<FArrayBox>& a_p,
              LevelData<FArrayBox>& a_q,
              const Real            a_oldTime,
              const Real            a_dt,
              PARKRHS*              a_rhsPtr)
    {
        this->advance(a_vel, a_p, a_q, a_oldTime, a_dt, a_rhsPtr);
    }

    Real
    controllerDt(const Real /* a_tol */,
                 const bool /* a_useImplicit */,
                 const bool /* a_useElementary */,
                 const bool /* a_usePI */,
                 const bool /* a_usePID */) const
    {
        return 1.0e100;
    }

protected:
    // Resets the data. The interp functions will be able to copy vel, p, and q
    // if called with a_time = a_oldTime.
    void
    setOldQ(const LevelData<FluxBox>&   a_velOld,
            const LevelData<FArrayBox>& a_pOld,
            const LevelData<FArrayBox>& a_qOld,
            const Real                  a_oldTime);

    // Sets x = y on all comps.
    // x and y must have the same number of comps.
    static void
    copy(LevelData<FluxBox>& a_x, const LevelData<FluxBox>& a_y);

    // Sets x = y on all comps.
    // x and y must have the same number of comps.
    static void
    copy(LevelData<FArrayBox>& a_x, const LevelData<FArrayBox>& a_y);

    // Sets x += b*y on all comps.
    // x and y must have the same number of comps.
    static void
    plus(LevelData<FluxBox>&       a_x,
         const Real                a_b,
         const LevelData<FluxBox>& a_y);

    // Sets x += b*y on all comps.
    // x and y must have the same number of comps.
    static void
    plus(LevelData<FArrayBox>&       a_x,
         const Real                  a_b,
         const LevelData<FArrayBox>& a_y);

    // Sets x = a*x + b*y on all comps.
    // x and y must have the same number of comps.
    static void
    axby(const Real                a_a,
         LevelData<FluxBox>&       a_x,
         const Real                a_b,
         const LevelData<FluxBox>& a_y);

    // Sets x = a*x + b*y on all comps.
    // x and y must have the same number of comps.
    static void
    axby(const Real                  a_a,
         LevelData<FArrayBox>&       a_x,
         const Real                  a_b,
         const LevelData<FArrayBox>& a_y);

    DisjointBoxLayout m_grids;
    const int         m_velNumComps;
    const IntVect     m_velGhostVect;
    const int         m_pNumComps;
    const IntVect     m_pGhostVect;
    const int         m_qNumComps;
    const IntVect     m_qGhostVect;

    // The stage data.
    LevelData<FluxBox> m_velOld, m_velNew;
    LevelData<FluxBox> m_kvelE, m_kvelI;

    LevelData<FArrayBox> m_pOld, m_pNew;

    LevelData<FArrayBox> m_qOld, m_qNew;
    LevelData<FArrayBox> m_kqE, m_kqI;

    Real m_dt;
    Real m_oldTime;
    Real m_newTime;

    bool m_Q0Ready;
    bool m_interpDataReady;
};


#endif  //!___ForwardEuler_H__INCLUDED___
