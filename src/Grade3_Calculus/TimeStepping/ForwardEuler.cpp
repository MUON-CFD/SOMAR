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
#include "SOMAR_Constants.H"
#include "Debug.H"
#include "ForwardEuler.H"
#include "SetValLevel.H"


// #ifndef NDEBUG
// // Debug mode
// #define nanCheck(x) checkForValidNAN(x)
// #else
// Release mode
#define nanCheck(x)
// #endif

#define timeEps (1.0e-6)


// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
ForwardEuler::ForwardEuler(const DisjointBoxLayout& a_grids,
                           const int                a_velNumComps,
                           const IntVect&           a_velGhostVect,
                           const int                a_pNumComps,
                           const IntVect&           a_pGhostVect,
                           const int                a_qNumComps,
                           const IntVect&           a_qGhostVect)
: m_grids(a_grids)
, m_velNumComps(a_velNumComps)
, m_velGhostVect(a_velGhostVect)
, m_pNumComps(a_pNumComps)
, m_pGhostVect(a_pGhostVect)
, m_qNumComps(a_qNumComps)
, m_qGhostVect(a_qGhostVect)
, m_dt(quietNAN)
, m_oldTime(quietNAN)
, m_newTime(quietNAN)
, m_Q0Ready(false)
, m_interpDataReady(false)
{
    m_velOld.define(m_grids, m_velNumComps, m_velGhostVect);
    m_velNew.define(m_grids, m_velNumComps, m_velGhostVect);
    m_kvelE.define(m_grids, m_velNumComps);
    m_kvelI.define(m_grids, m_velNumComps);

    debugInitLevel(m_velOld);
    debugInitLevel(m_velNew);
    debugInitLevel(m_kvelE);
    debugInitLevel(m_kvelI);

    m_pOld.define(m_grids, m_pNumComps, m_pGhostVect);
    m_pNew.define(m_grids, m_pNumComps, m_pGhostVect);

    debugInitLevel(m_pOld);
    debugInitLevel(m_pNew);

    m_qOld.define(m_grids, m_qNumComps, m_qGhostVect);
    m_qNew.define(m_grids, m_qNumComps, m_qGhostVect);
    m_kqE.define(m_grids, m_qNumComps);
    m_kqI.define(m_grids, m_qNumComps);

    debugInitLevel(m_qOld);
    debugInitLevel(m_qNew);
    debugInitLevel(m_kqE);
    debugInitLevel(m_kqI);
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
ForwardEuler::~ForwardEuler()
{
}


// -----------------------------------------------------------------------------
// The main timestepper.
// -----------------------------------------------------------------------------
void
ForwardEuler::advance(LevelData<FluxBox>&   a_vel,
                      LevelData<FArrayBox>& a_p,
                      LevelData<FArrayBox>& a_q,
                      const Real            a_oldTime,
                      const Real            a_dt,
                      PARKRHS*              a_rhsPtr)
{
    // Sanity checks
    CH_assert(a_vel.getBoxes().compatible(m_grids));
    CH_assert(a_p.getBoxes().compatible(m_grids));
    CH_assert(a_q.getBoxes().compatible(m_grids));

    CH_assert(a_vel.nComp() == m_velNumComps);
    CH_assert(a_p.nComp() == m_pNumComps);
    CH_assert(a_q.nComp() == m_qNumComps);

    nanCheck(a_vel);
    nanCheck(a_p);
    nanCheck(a_q);

    // Initialization. This sets Q[0].
    this->setOldQ(a_vel, a_p, a_q, a_oldTime);
    m_dt = a_dt;


    // --- Begin stage ---
    pout() << Format::fixed << "FE timestep..."
           << "\n"
           << Format::indent() << flush;

    // Evaluate forces.
    a_rhsPtr->setExplicitRHS(m_kvelE, m_kqE, a_vel, a_p, a_q, a_oldTime, a_dt);
    a_rhsPtr->setImplicitRHS(m_kvelI, m_kqI, a_vel, a_p, a_q, a_oldTime, 0.5 * a_dt);

    nanCheck(m_kvelE);
    nanCheck(m_kqE);
    nanCheck(m_kvelI);
    nanCheck(m_kqI);

    // Put all (scaled) stage forces kqE.
    // kQE <-- kQE + 0.5*kQI
    ForwardEuler::plus(m_kvelE, 0.5, m_kvelI);
    ForwardEuler::plus(m_kqE, 0.5, m_kqI);

    // Update q (Ri)
    // Q += dt*kQE
    ForwardEuler::plus(a_vel, a_dt, m_kvelE);
    ForwardEuler::plus(a_q, a_dt, m_kqE);

    // Remove lagged pressure gradient
    a_rhsPtr->projectPredict(a_vel, a_p, a_oldTime, a_dt);

    // Solve implicit eqs for this stage's state.
    a_rhsPtr->solveImplicit(a_vel, a_q, 0.5 * a_dt, a_oldTime, 0.5 * a_dt);

    // Correct the projected state.
    a_rhsPtr->projectCorrect(a_vel, a_p, a_oldTime, a_dt);

    nanCheck(a_vel);
    nanCheck(a_p);
    nanCheck(a_q);
    pout() << Format::unindent << flush;
    // --- End stage ---


    // Save state data
    m_oldTime = a_oldTime;
    m_newTime = a_oldTime + a_dt;

    ForwardEuler::copy(m_velNew, a_vel);
    ForwardEuler::copy(m_pNew, a_p);
    ForwardEuler::copy(m_qNew, a_q);

    m_Q0Ready         = true;
    m_interpDataReady = true;

    pout() << Format::unindent << flush;
}


// -----------------------------------------------------------------------------
// Interpolate data in time using a linear approx.
// This will use piecewise linear interpolation and return how far along the
// initial and final data the interpolation took place (0 = old time,
// 1 = new time). If this returns a number outside of [0,1], then we
// extrapolated and you need to worry about stability.
// -----------------------------------------------------------------------------
Real
ForwardEuler::velTimeInterp(LevelData<FluxBox>& a_vel,
                            const Real          a_time,
                            int                 a_srcComp,
                            int                 a_destComp,
                            int                 a_numComp) const
{
    CH_assert(a_numComp <= m_velNumComps);
    CH_assert(a_vel.getBoxes().compatible(m_grids));
    CH_assert(a_vel.getBoxes().physDomain().size() ==
              m_grids.physDomain().size());

    DataIterator dit = m_grids.dataIterator();

    if (a_numComp == -1) a_numComp = m_velNumComps - a_srcComp;

    // If we are at m_times[0], just copy the old data.
    CH_assert(m_Q0Ready);
    if (abs(a_time - m_oldTime) < timeEps) {
        for (dit.reset(); dit.ok(); ++dit) {
            a_vel[dit].copy(m_velOld[dit], a_srcComp, a_destComp, a_numComp);
        }
        return 0.0;
    }

    // Linear interpolation...
    CH_assert(m_interpDataReady);

    // Compute theta and if nearby, snap into 0 or 1 values to prevent
    // extrapolation when not expected.
    Real theta = (a_time - m_oldTime) / m_dt;

    if (abs(theta) < timeEps) {
        // A bit redundant...
        for (dit.reset(); dit.ok(); ++dit) {
            a_vel[dit].copy(m_velOld[dit], a_srcComp, a_destComp, a_numComp);
        }
        return 0.0;
    } else if (abs(theta - 1.0) < timeEps) {
        for (dit.reset(); dit.ok(); ++dit) {
            a_vel[dit].copy(m_velNew[dit], a_srcComp, a_destComp, a_numComp);
        }
        return 1.0;
    }

    // Perform interpolation.
    for (dit.reset(); dit.ok(); ++dit) {
        for (int dir = 0; dir < SpaceDim; ++dir) {
            a_vel[dit][dir].setVal(0.0);

            a_vel[dit][dir].plus(m_velOld[dit][dir],
                                 1.0 - theta,
                                 a_srcComp,
                                 a_destComp,
                                 a_numComp);

            a_vel[dit][dir].plus(
                m_velNew[dit][dir], theta, a_srcComp, a_destComp, a_numComp);
        }
    }

    CH_assert(0.0 <= theta && theta <= 1.0);
    return theta;
}


// -----------------------------------------------------------------------------
// p version.
// -----------------------------------------------------------------------------
Real
ForwardEuler::pTimeInterp(LevelData<FArrayBox>& a_p,
                          const Real            a_time,
                          int                   a_srcComp,
                          int                   a_destComp,
                          int                   a_numComp) const
{
    CH_assert(a_numComp <= m_pNumComps);
    CH_assert(a_p.getBoxes().compatible(m_grids));
    CH_assert(a_p.getBoxes().physDomain().size() ==
              m_grids.physDomain().size());

    DataIterator dit = m_grids.dataIterator();

    if (a_numComp == -1) a_numComp = m_pNumComps - a_srcComp;

    // If we are at m_times[0], just copy the old data.
    CH_assert(m_Q0Ready);
    if (abs(a_time - m_oldTime) < timeEps) {
        for (dit.reset(); dit.ok(); ++dit) {
            a_p[dit].copy(m_pOld[dit], a_srcComp, a_destComp, a_numComp);
        }
        return 0.0;
    }

    // Linear interpolation...
    CH_assert(m_interpDataReady);

    // Compute theta and if nearby, snap into 0 or 1 values to prevent
    // extrapolation when not expected.
    Real theta = (a_time - m_oldTime) / m_dt;

    if (abs(theta) < timeEps) {
        // A bit redundant...
        for (dit.reset(); dit.ok(); ++dit) {
            a_p[dit].copy(m_pOld[dit], a_srcComp, a_destComp, a_numComp);
        }
        return 0.0;
    } else if (abs(theta - 1.0) < timeEps) {
        for (dit.reset(); dit.ok(); ++dit) {
            a_p[dit].copy(m_pNew[dit], a_srcComp, a_destComp, a_numComp);
        }
        return 1.0;
    }

    // Perform interpolation.
    for (dit.reset(); dit.ok(); ++dit) {
        a_p[dit].setVal(0.0);
        a_p[dit].plus(
            m_pOld[dit], 1.0 - theta, a_srcComp, a_destComp, a_numComp);
        a_p[dit].plus(m_pNew[dit], theta, a_srcComp, a_destComp, a_numComp);
    }

    CH_assert(0.0 <= theta && theta <= 1.0);
    return theta;
}


// -----------------------------------------------------------------------------
// q version.
// -----------------------------------------------------------------------------
Real
ForwardEuler::qTimeInterp(LevelData<FArrayBox>& a_q,
                          const Real            a_time,
                          int                   a_srcComp,
                          int                   a_destComp,
                          int                   a_numComp) const
{
    CH_assert(a_numComp <= m_qNumComps);
    CH_assert(a_q.getBoxes().compatible(m_grids));
    CH_assert(a_q.getBoxes().physDomain().size() ==
              m_grids.physDomain().size());

    DataIterator dit = m_grids.dataIterator();

    if (a_numComp == -1) a_numComp = m_qNumComps - a_srcComp;

    // If we are at m_times[0], just copy the old data.
    CH_assert(m_Q0Ready);
    if (abs(a_time - m_oldTime) < timeEps) {
        for (dit.reset(); dit.ok(); ++dit) {
            a_q[dit].copy(m_qOld[dit], a_srcComp, a_destComp, a_numComp);
        }
        return 0.0;
    }

    // Linear interpolation...
    CH_assert(m_interpDataReady);

    // Compute theta and if nearby, snap into 0 or 1 values to prevent
    // extrapolation when not expected.
    Real theta = (a_time - m_oldTime) / m_dt;

    if (abs(theta) < timeEps) {
        // A bit redundant...
        for (dit.reset(); dit.ok(); ++dit) {
            a_q[dit].copy(m_qOld[dit], a_srcComp, a_destComp, a_numComp);
        }
        return 0.0;
    } else if (abs(theta - 1.0) < timeEps) {
        for (dit.reset(); dit.ok(); ++dit) {
            a_q[dit].copy(m_qNew[dit], a_srcComp, a_destComp, a_numComp);
        }
        return 1.0;
    }

    // Perform interpolation.
    for (dit.reset(); dit.ok(); ++dit) {
        a_q[dit].setVal(0.0);
        a_q[dit].plus(
            m_qOld[dit], 1.0 - theta, a_srcComp, a_destComp, a_numComp);
        a_q[dit].plus(m_qNew[dit], theta, a_srcComp, a_destComp, a_numComp);
    }

    CH_assert(0.0 <= theta && theta <= 1.0);
    return theta;
}


// -----------------------------------------------------------------------------
// Resets the data. The interp functions will be able to copy vel, p, and q
// if called with a_time = a_oldTime.
// -----------------------------------------------------------------------------
void
ForwardEuler::setOldQ(const LevelData<FluxBox>&   a_velOld,
                      const LevelData<FArrayBox>& a_pOld,
                      const LevelData<FArrayBox>& a_qOld,
                      const Real                  a_oldTime)
{
    CH_assert(a_velOld.nComp() == m_velNumComps);
    CH_assert(a_velOld.getBoxes().compatible(m_grids));
    CH_assert(a_velOld.getBoxes().physDomain().size() ==
              m_grids.physDomain().size());

    CH_assert(a_pOld.nComp() == m_pNumComps);
    CH_assert(a_pOld.getBoxes().compatible(m_grids));
    CH_assert(a_pOld.getBoxes().physDomain().size() ==
              m_grids.physDomain().size());

    CH_assert(a_qOld.nComp() == m_qNumComps);
    CH_assert(a_qOld.getBoxes().compatible(m_grids));
    CH_assert(a_qOld.getBoxes().physDomain().size() ==
              m_grids.physDomain().size());

    debugInitLevel(m_velOld);
    debugInitLevel(m_velNew);
    debugInitLevel(m_kvelE);
    debugInitLevel(m_kvelI);

    debugInitLevel(m_pOld);
    debugInitLevel(m_pNew);

    debugInitLevel(m_qOld);
    debugInitLevel(m_qNew);
    debugInitLevel(m_kqE);
    debugInitLevel(m_kqI);

    m_oldTime = quietNAN;
    m_newTime = quietNAN;
    m_dt      = quietNAN;

    DataIterator dit = m_grids.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        m_velOld[dit].copy(a_velOld[dit]);
        m_pOld[dit].copy(a_pOld[dit]);
        m_qOld[dit].copy(a_qOld[dit]);
    }
    m_oldTime = a_oldTime;

    m_Q0Ready         = true;
    m_interpDataReady = false;
}


// -----------------------------------------------------------------------------
// Sets x = y on all comps.
// x and y must have the same number of comps.
// -----------------------------------------------------------------------------
void
ForwardEuler::copy(LevelData<FluxBox>& a_x, const LevelData<FluxBox>& a_y)
{
    CH_assert(a_x.getBoxes() == a_y.getBoxes());
    CH_assert(a_x.nComp() == a_y.nComp());

    DataIterator dit = a_x.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        for (int dir = 0; dir < SpaceDim; ++dir) {
            a_x[dit][dir].copy(a_y[dit][dir]);
        }
    }
}


// -----------------------------------------------------------------------------
// Sets x = y on all comps.
// x and y must have the same number of comps.
// -----------------------------------------------------------------------------
void
ForwardEuler::copy(LevelData<FArrayBox>& a_x, const LevelData<FArrayBox>& a_y)
{
    CH_assert(a_x.getBoxes() == a_y.getBoxes());
    CH_assert(a_x.nComp() == a_y.nComp());

    DataIterator dit = a_x.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        a_x[dit].copy(a_y[dit]);
    }
}


// -----------------------------------------------------------------------------
// Sets x += b*y on all comps.
// x and y must have the same number of comps.
// -----------------------------------------------------------------------------
void
ForwardEuler::plus(LevelData<FluxBox>&       a_x,
                   const Real                a_b,
                   const LevelData<FluxBox>& a_y)
{
    CH_assert(a_x.getBoxes() == a_y.getBoxes());
    CH_assert(a_x.nComp() == a_y.nComp());

    DataIterator dit = a_x.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        for (int dir = 0; dir < SpaceDim; ++dir) {
            a_x[dit][dir].plus(a_y[dit][dir], a_b);
        }
    }
}


// -----------------------------------------------------------------------------
// Sets x += b*y on all comps.
// x and y must have the same number of comps.
// -----------------------------------------------------------------------------
void
ForwardEuler::plus(LevelData<FArrayBox>&       a_x,
                   const Real                  a_b,
                   const LevelData<FArrayBox>& a_y)
{
    CH_assert(a_x.getBoxes() == a_y.getBoxes());
    CH_assert(a_x.nComp() == a_y.nComp());

    DataIterator dit = a_x.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        a_x[dit].plus(a_y[dit], a_b);
    }
}


// -----------------------------------------------------------------------------
// Sets x = a*x + b*y on all comps.
// x and y must have the same number of comps.
// -----------------------------------------------------------------------------
void
ForwardEuler::axby(const Real                a_a,
                   LevelData<FluxBox>&       a_x,
                   const Real                a_b,
                   const LevelData<FluxBox>& a_y)
{
    CH_assert(a_x.getBoxes() == a_y.getBoxes());
    CH_assert(a_x.nComp() == a_y.nComp());

    DataIterator dit = a_x.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        for (int dir = 0; dir < SpaceDim; ++dir) {
            a_x[dit][dir].mult(a_a);
            a_x[dit][dir].plus(a_y[dit][dir], a_b);
        }
    }
}


// -----------------------------------------------------------------------------
// Sets x = a*x + b*y on all comps.
// x and y must have the same number of comps.
// -----------------------------------------------------------------------------
void
ForwardEuler::axby(const Real                  a_a,
                   LevelData<FArrayBox>&       a_x,
                   const Real                  a_b,
                   const LevelData<FArrayBox>& a_y)
{
    CH_assert(a_x.getBoxes() == a_y.getBoxes());
    CH_assert(a_x.nComp() == a_y.nComp());

    DataIterator dit = a_x.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        a_x[dit].mult(a_a);
        a_x[dit].plus(a_y[dit], a_b);
    }
}
