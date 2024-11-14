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
#include "RKW3CN.H"
#include "SOMAR_Constants.H"
#include "Debug.H"
#include "SetValLevel.H"


#ifndef NDEBUG
// Debug mode
#define nanCheck(x) checkForValidNAN(x)
#else
// Release mode
#define nanCheck(x)
#endif



// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
RKW3CN::RKW3CN(const DisjointBoxLayout& a_grids,
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
RKW3CN::~RKW3CN()
{
}


#if 0
// -----------------------------------------------------------------------------
// The main timestepper.
// -----------------------------------------------------------------------------
void
RKW3CN::advance(LevelData<FluxBox>&   a_vel,
                LevelData<FArrayBox>& a_p,
                LevelData<FArrayBox>& a_q,
                const Real            a_oldTime,
                const Real            a_dt,
                PARKRHS*              a_rhsPtr)
{
    // Sanity checks
    CH_assert(a_vel.getBoxes().compatible(m_grids));
    CH_assert(a_p  .getBoxes().compatible(m_grids));
    CH_assert(a_q  .getBoxes().compatible(m_grids));

    CH_assert(a_vel.nComp() == m_velNumComps);
    CH_assert(a_p  .nComp() == m_pNumComps);
    CH_assert(a_q  .nComp() == m_qNumComps);

    checkForValidNAN(a_vel);
    checkForValidNAN(a_p  );
    checkForValidNAN(a_q  );

    DataIterator dit = m_grids.dataIterator();

    // Initialization. This sets Q[0].
    this->setOldQ(a_vel, a_p, a_q, a_oldTime);
    m_dt = a_dt;

    // Collect coefficients
    const Real hb[3] = {
        (8./15.) * a_dt,
        (2./15.) * a_dt,
        (1./ 3.) * a_dt
    };
    static constexpr Real betab[3] = {1.,  25./8.,  9./4.};
    static constexpr Real zetab[3] = {0., -17./8., -5./4.};
    static constexpr Real c[3]     = {0.,  8./15.,  2./3.};
    static constexpr bool projEveryStage = true;


    // ------ Main loop ------
    LevelData<FluxBox>   velm2(m_grids, a_vel.nComp(), a_vel.ghostVect());
    LevelData<FArrayBox> pm2  (m_grids, a_p  .nComp(), a_p  .ghostVect());
    LevelData<FArrayBox> qm2  (m_grids, a_q  .nComp(), a_q  .ghostVect());

    // Loop over the stages to compute the forces m_kqE, m_kqI, and m_Gp.
    for (int i = 0; i < 3; ++i) {
        pout() << Format::fixed;
        pout() << "RK Stage " << i << "\n";
        pout() << Format::indent() << flush;

        // Set time values for this stage.
        const Real stageTime = a_oldTime + a_dt * c[i];
        const Real gammaDt   = 0.5 * hb[i];

        // Evaluate this stage's forces.
        a_rhsPtr->setExplicitRHS(m_kvelE, m_kqE, a_vel, a_p, a_q, stageTime, hb[i] * betab[i]);
        a_rhsPtr->setImplicitRHS(m_kvelI, m_kqI, a_vel, a_p, a_q, stageTime, hb[i] - gammaDt);

        checkForValidNAN(m_kvelE);
        checkForValidNAN(m_kqE  );
        checkForValidNAN(m_kvelI);
        checkForValidNAN(m_kqI  );

        // Put all (scaled) stage forces kqE.
        // kQE <-- betab*kQE + 0.5*kQI
        RKW3CN::axby(betab[i], m_kvelE, 0.5, m_kvelI);
        RKW3CN::axby(betab[i], m_kqE  , 0.5, m_kqI  );

        if (i > 0) {
            // Use m_kqI to hold old-stage kqE.
            const Real oldStageTime = a_oldTime + a_dt * c[i - 1];
            a_rhsPtr->setExplicitRHS(m_kvelI, m_kqI, velm2, pm2, qm2, oldStageTime, hb[i] * zetab[i]);

            // Add contribution to kQE.
            // kQE += zetab * kQI
            RKW3CN::plus(m_kvelE, zetab[i], m_kvelI);
            RKW3CN::plus(m_kqE  , zetab[i], m_kqI  );
        }

        // Keep a copy of old state.
        if (i < 2) {
            RKW3CN::copy(velm2, a_vel);
            RKW3CN::copy(pm2  , a_p  );
            RKW3CN::copy(qm2  , a_q  );
        }

        // Update q (Ri)
        // Q += hb*kQE
        RKW3CN::plus(a_vel, hb[i], m_kvelE);
        RKW3CN::plus(a_q  , hb[i], m_kqE  );

        // Implicit solves...
        if (projEveryStage) {
            // proj predict -> viscous solve -> proj correct
            a_rhsPtr->projectPredict(a_vel, a_p, stageTime, hb[i]             );
            a_rhsPtr->solveImplicit (a_vel, a_q,   gammaDt, stageTime, gammaDt);
            a_rhsPtr->projectCorrect(a_vel, a_p, stageTime, hb[i]             );
        } else {
            // proj predict -> viscous solve -> unproject and wait until end
            // a_rhsPtr->projectPredict(a_vel, a_p, stageTime, hb[i]    );
            a_rhsPtr->solveImplicit (a_vel, a_q,   gammaDt, stageTime, gammaDt);
            // a_rhsPtr->projectPredict(a_vel, a_p, stageTime, -hb[i]   );
        }

        checkForValidNAN(a_vel);
        checkForValidNAN(a_p);
        checkForValidNAN(a_q);
        pout() << Format::unindent << flush;

    }  // end loop over stages


    // Save state data
    m_oldTime = a_oldTime;
    m_newTime = a_oldTime + a_dt;

    // Project if needed.
    if (!projEveryStage) {
        a_rhsPtr->projectPredict(a_vel, a_p, m_newTime, a_dt);
        a_rhsPtr->projectCorrect(a_vel, a_p, m_newTime, a_dt);
    }


    RKW3CN::copy(m_velNew, a_vel);
    RKW3CN::copy(m_pNew  , a_p  );
    RKW3CN::copy(m_qNew  , a_q  );

    // Allow user to perform postStep operations (e.g., preparing plots)
    a_rhsPtr->postStep(m_velNew, m_pNew, m_qNew, m_newTime);

    // Also, save a copy of the (properly scaled) RHS.
    {
        a_rhsPtr->setImplicitRHS(
            m_kvelI, m_kqI, m_velNew, m_pNew, m_qNew, m_newTime, 0.0);
        RKW3CN::plus(m_kvelE, 0.5, m_kvelI);
        RKW3CN::plus(  m_kqE, 0.5,   m_kqI);
    }

    m_Q0Ready         = true;
    m_interpDataReady = true;

    pout() << Format::unindent << flush;
}

#else

void
RKW3CN::advance(LevelData<FluxBox>&   a_vel,
                LevelData<FArrayBox>& a_p,
                LevelData<FArrayBox>& a_q,
                const Real            a_oldTime,
                const Real            a_dt,
                PARKRHS*              a_rhsPtr)
{
    // Sanity checks
    CH_assert(a_vel.getBoxes().compatible(m_grids));
    CH_assert(a_p  .getBoxes().compatible(m_grids));
    CH_assert(a_q  .getBoxes().compatible(m_grids));

    CH_assert(a_vel.nComp() == m_velNumComps);
    CH_assert(a_p  .nComp() == m_pNumComps);
    CH_assert(a_q  .nComp() == m_qNumComps);

    checkForValidNAN(a_vel);
    checkForValidNAN(a_p  );
    checkForValidNAN(a_q  );

    DataIterator dit = m_grids.dataIterator();

    // Initialization. This sets Q[0].
    this->setOldQ(a_vel, a_p, a_q, a_oldTime);
    m_dt = a_dt;

    // Collect coefficients
    const Real hb[3] = {
        (8./15.) * a_dt,
        (2./15.) * a_dt,
        (1./ 3.) * a_dt
    };
    static constexpr Real betab[3] = {1.,  25./8.,  9./4.};
    static constexpr Real zetab[3] = {0., -17./8., -5./4.};
    static constexpr Real c[3]     = {0.,  8./15.,  2./3.};
    static constexpr bool projEveryStage = true;

    static constexpr Real refluxDt = 0.0;


    // ------ Main loop ------
    LevelData<FluxBox>   kvelE2(m_grids, m_kvelE.nComp(), m_kvelE.ghostVect());
    LevelData<FArrayBox> kqE2  (m_grids, m_kqE  .nComp(), m_kqE  .ghostVect());

    // Loop over the stages to compute the forces m_kqE, m_kqI, and m_Gp.
    for (int i = 0; i < 3; ++i) {
        pout() << Format::fixed;
        pout() << "RK Stage " << i << "\n";
        pout() << Format::indent() << flush;

        // Set time values for this stage.
        const Real stageTime = a_oldTime + a_dt * c[i];

        // Evaluate this stage's forces.
        a_rhsPtr->setExplicitRHS(m_kvelE, m_kqE, a_vel, a_p, a_q, stageTime, refluxDt);
        a_rhsPtr->setImplicitRHS(m_kvelI, m_kqI, a_vel, a_p, a_q, stageTime, refluxDt);

        checkForValidNAN(m_kvelE);
        checkForValidNAN(m_kqE  );
        checkForValidNAN(m_kvelI);
        checkForValidNAN(m_kqI  );

        // Add current k*E
        RKW3CN::plus(a_vel, hb[i] * betab[i], m_kvelE);
        RKW3CN::plus(a_q  , hb[i] * betab[i], m_kqE  );

        // Add previous k*E
        if (i > 0) {
            RKW3CN::plus(a_vel, hb[i] * zetab[i], kvelE2);
            RKW3CN::plus(a_q  , hb[i] * zetab[i], kqE2  );
        }

        // Add current k*I
        RKW3CN::plus(a_vel, hb[i] * 0.5, m_kvelI);
        RKW3CN::plus(a_q  , hb[i] * 0.5, m_kqI  );

        // Save k*E for later stages, if needed.
        if (i < 2) {
            RKW3CN::copy(kvelE2, m_kvelE);
            RKW3CN::copy(kqE2  , m_kqE  );
        }

        // Implicit solves...
        if (projEveryStage) {
            // proj predict -> viscous solve -> proj correct
            a_rhsPtr->projectPredict(a_vel, a_p, stageTime, hb[i]);
            a_rhsPtr->solveImplicit (a_vel, a_q, 0.5 * hb[i], stageTime, refluxDt);
            a_rhsPtr->projectCorrect(a_vel, a_p, stageTime, hb[i]);
        } else {
            // proj predict -> viscous solve -> unproject and wait until end
            // a_rhsPtr->projectPredict(a_vel, a_p, stageTime, hb[i]    );
            a_rhsPtr->solveImplicit (a_vel, a_q, 0.5 * hb[i], stageTime, refluxDt);
            // a_rhsPtr->projectPredict(a_vel, a_p, stageTime, -hb[i]   );
        }

        checkForValidNAN(a_vel);
        checkForValidNAN(a_p);
        checkForValidNAN(a_q);
        pout() << Format::unindent << flush;

    }  // end loop over stages


    // Save state data
    m_oldTime = a_oldTime;
    m_newTime = a_oldTime + a_dt;

    // Project if needed.
    if (!projEveryStage) {
        a_rhsPtr->projectPredict(a_vel, a_p, m_newTime, a_dt);
        a_rhsPtr->projectCorrect(a_vel, a_p, m_newTime, a_dt);
    }


    RKW3CN::copy(m_velNew, a_vel);
    RKW3CN::copy(m_pNew  , a_p  );
    RKW3CN::copy(m_qNew  , a_q  );

    // debugInitLevel(m_kvelE);
    // debugInitLevel(m_kvelI);
    // debugInitLevel(m_kqE);
    // debugInitLevel(m_kqI);
    setValLevel(m_kvelE, quietNAN);
    setValLevel(m_kvelI, quietNAN);
    setValLevel(m_kqE, quietNAN);
    setValLevel(m_kqE, quietNAN);

    // Allow user to perform postStep operations (e.g., preparing plots)
    a_rhsPtr->postStep(m_velNew, m_pNew, m_qNew, m_newTime);

    m_Q0Ready         = true;
    m_interpDataReady = true;

    pout() << Format::unindent << flush;
}
#endif

// -----------------------------------------------------------------------------
// A cheap Forward Euler timestepper.
// This is useful for generating rough estimates.
// -----------------------------------------------------------------------------
void
RKW3CN::FEadvance(LevelData<FluxBox>&   a_vel,
                  LevelData<FArrayBox>& a_p,
                  LevelData<FArrayBox>& a_q,
                  const Real            a_oldTime,
                  const Real            a_dt,
                  PARKRHS*              a_rhsPtr)
{
    // Sanity checks
    CH_assert(a_vel.getBoxes().compatible(m_grids));
    CH_assert(a_p  .getBoxes().compatible(m_grids));
    CH_assert(a_q  .getBoxes().compatible(m_grids));

    CH_assert(a_vel.nComp() == m_velNumComps);
    CH_assert(a_p  .nComp() == m_pNumComps);
    CH_assert(a_q  .nComp() == m_qNumComps);

    checkForValidNAN(a_vel);
    checkForValidNAN(a_p  );
    checkForValidNAN(a_q  );

    this->setOldQ(a_vel, a_p, a_q, a_oldTime);
    m_dt = a_dt;


    // --- Begin stage ---
    pout() << Format::fixed << "FE timestep..."
           << "\n"
           << Format::indent() << flush;

    // Evaluate forces.
    a_rhsPtr->setExplicitRHS(m_kvelE, m_kqE, a_vel, a_p, a_q, a_oldTime, a_dt);
    a_rhsPtr->setImplicitRHS(m_kvelI, m_kqI, a_vel, a_p, a_q, a_oldTime, a_dt);

    checkForValidNAN(m_kvelE);
    checkForValidNAN(m_kqE);
    checkForValidNAN(m_kvelI);
    checkForValidNAN(m_kqI);

    // Update q (Ri)
    // Q += dt*(kQE + kQI)
    this->plus(a_vel, a_dt, m_kvelE);
    this->plus(a_vel, a_dt, m_kvelI);
    this->plus(a_q, a_dt, m_kqE);
    this->plus(a_q, a_dt, m_kqI);

    // Remove lagged pressure gradient
    a_rhsPtr->projectPredict(a_vel, a_p, a_oldTime, a_dt);

    // Solve implicit eqs for this stage's state.
    a_rhsPtr->solveImplicit(a_vel, a_q, 0.5 * a_dt, a_oldTime, 0.5 * a_dt);

    // Correct the projected state.
    a_rhsPtr->projectCorrect(a_vel, a_p, a_oldTime, a_dt);

    checkForValidNAN(a_vel);
    checkForValidNAN(a_p);
    checkForValidNAN(a_q);
    pout() << Format::unindent << flush;
    // --- End stage ---


    // Save state data
    m_oldTime = a_oldTime;
    m_newTime = a_oldTime + a_dt;

    // The initial state is assumed to be accurate...
    m_Q0Ready = true;
    // ...but the final state is not. We do not want to use it for BC
    // interpolations.
    m_interpDataReady = false;
    // ...also, we will not call postStep. I doubt you want to prepare
    // plots using this final state.
}


// -----------------------------------------------------------------------------
// Interpolate data in time using a linear approx.
// This will use piecewise linear interpolation and return how far along the
// initial and final data the interpolation took place (0 = old time,
// 1 = new time). If this returns a number outside of [0,1], then we
// extrapolated and you need to worry about stability.
// -----------------------------------------------------------------------------
Real
RKW3CN::velTimeInterp(LevelData<FluxBox>& a_vel,
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

            a_vel[dit][dir].plus(m_velNew[dit][dir],
                                 theta,
                                 a_srcComp,
                                 a_destComp,
                                 a_numComp);
        }
    }

    CH_assert(0.0 <= theta && theta <= 1.0);
    return theta;
}


// -----------------------------------------------------------------------------
// p version.
// -----------------------------------------------------------------------------
Real
RKW3CN::pTimeInterp(LevelData<FArrayBox>& a_p,
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
        a_p[dit].plus(m_pOld[dit], 1.0-theta, a_srcComp, a_destComp, a_numComp);
        a_p[dit].plus(m_pNew[dit], theta, a_srcComp, a_destComp, a_numComp);
    }

    CH_assert(0.0 <= theta && theta <= 1.0);
    return theta;
}


// -----------------------------------------------------------------------------
// q version.
// -----------------------------------------------------------------------------
Real
RKW3CN::qTimeInterp(LevelData<FArrayBox>& a_q,
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
        a_q[dit].plus(m_qOld[dit], 1.0-theta, a_srcComp, a_destComp, a_numComp);
        a_q[dit].plus(m_qNew[dit], theta, a_srcComp, a_destComp, a_numComp);
    }

    CH_assert(0.0 <= theta && theta <= 1.0);
    return theta;
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void
RKW3CN::velTimeDeriv(LevelData<FluxBox>& a_deriv,
                     const Real          /*a_time*/,
                     int                 a_srcComp,
                     int                 a_destComp,
                     int                 a_numComp) const
{
    CH_assert(a_numComp <= m_velNumComps);
    CH_assert(a_deriv.getBoxes().compatible(m_grids));
    CH_assert(a_deriv.getBoxes().physDomain().size() ==
              m_grids.physDomain().size());

    CH_assert(m_Q0Ready);
    CH_assert(m_interpDataReady);

    if (a_numComp == -1) a_numComp = m_velNumComps - a_srcComp;

    // Perform interpolation.
    DataIterator dit = m_grids.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        for (int dir = 0; dir < SpaceDim; ++dir) {
            a_deriv[dit][dir].setVal(0.0);

            a_deriv[dit][dir].plus(m_velNew[dit][dir],
                                   1.0 / m_dt,
                                   a_srcComp,
                                   a_destComp,
                                   a_numComp);

            a_deriv[dit][dir].plus(m_velOld[dit][dir],
                                   -1.0 / m_dt,
                                   a_srcComp,
                                   a_destComp,
                                   a_numComp);
        }
    }
}


// -----------------------------------------------------------------------------
// Resets the data. The interp functions will be able to copy vel, p, and q
// if called with a_time = a_oldTime.
// -----------------------------------------------------------------------------
void
RKW3CN::setOldQ(const LevelData<FluxBox>&   a_velOld,
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
    m_dt = quietNAN;

    DataIterator dit = m_grids.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        m_velOld[dit].copy(a_velOld[dit]);
        m_pOld[dit].copy(a_pOld[dit]);
        m_qOld[dit].copy(a_qOld[dit]);
    }
    m_oldTime = a_oldTime;

    m_Q0Ready = true;
    m_interpDataReady = false;
}


// -----------------------------------------------------------------------------
// Sets x = y on all comps.
// x and y must have the same number of comps.
// -----------------------------------------------------------------------------
void
RKW3CN::copy(LevelData<FluxBox>&       a_x,
             const LevelData<FluxBox>& a_y)
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
RKW3CN::copy(LevelData<FArrayBox>&       a_x,
             const LevelData<FArrayBox>& a_y)
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
RKW3CN::plus(LevelData<FluxBox>&       a_x,
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
RKW3CN::plus(LevelData<FArrayBox>&       a_x,
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
RKW3CN::axby(const Real                a_a,
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
RKW3CN::axby(const Real                  a_a,
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
