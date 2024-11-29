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
#ifndef M0091158a93725eb57192ee021f0139bf
#error PARKI.H is not for public use! Try including PARK.H instead.
#endif

#include "Analysis.H"
#include "Debug.H"
#include "SOMAR_Constants.H"
#include "SetValLevel.H"

#ifndef NDEBUG
// Debug mode
#define nanCheck(x) checkForValidNAN(x)
#else
// Release mode
#define nanCheck(x)
#endif


template <class RKC> const int PARK<RKC>::s_orderE;
template <class RKC> const int PARK<RKC>::s_orderI;
template <class RKC> const int PARK<RKC>::s_orderEmbeddedE;
template <class RKC> const int PARK<RKC>::s_orderEmbeddedI;
template <class RKC> const int PARK<RKC>::s_stageOrderE;
template <class RKC> const int PARK<RKC>::s_stageOrderI;


// -----------------------------------------------------------------------------
template <class RKC>
PARK<RKC>::PARK(const DisjointBoxLayout& a_grids,
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
, m_needsAssembly(true)
, m_needsAssemblyProjection(true)
, m_oldTime(quietNAN)
, m_newTime(quietNAN)
, m_dt(quietNAN)
, m_Q0Ready(false)
, m_interpDataReady(false)
, m_errorHistory()
, m_controllerHistory{-1.0, -1.0, -1.0}
{
    for (int i = 0; i < RKC::numStages; ++i) {
        m_vel[i].define(a_grids, a_velNumComps, m_velGhostVect);
        m_kvelE[i].define(a_grids, a_velNumComps);
        m_kvelI[i].define(a_grids, a_velNumComps);

        m_p[i].define(a_grids, a_pNumComps, a_pGhostVect);

        m_q[i].define(a_grids, a_qNumComps, a_qGhostVect);
        m_kqE[i].define(a_grids, a_qNumComps);
        m_kqI[i].define(a_grids, a_qNumComps);

        debugInitLevel(m_vel[i]);
        debugInitLevel(m_kvelE[i]);
        debugInitLevel(m_kvelI[i]);

        debugInitLevel(m_p[i]);

        debugInitLevel(m_q[i]);
        debugInitLevel(m_kqE[i]);
        debugInitLevel(m_kqI[i]);
    }

    // PARK methods have 3 strict requirements.
    CH_verify(abs(RKC::c[0]) < smallReal);
    CH_verify(abs(RKC::aI[0][0]) < smallReal);
    CH_verify(RKC::numStages == 1 || abs(RKC::c[1]) > smallReal);

    pout() << "Explicit order = " << s_orderE << '\n';
    pout() << "Implicit order = " << s_orderI << '\n';
    if (this->hasEmbeddedScheme()) {
        pout() << "Explicit embedded order = " << s_orderEmbeddedE << '\n';
        pout() << "Implicit embedded order = " << s_orderEmbeddedI << '\n';
    } else {
        pout() << "No embedded scheme.\n";
    }
    pout() << "Explicit stage order = " << s_stageOrderE << '\n';
    pout() << "Implicit stage order = " << s_stageOrderI << '\n';
    pout() << std::endl;

    // m_prevIdxForProj
    {
        int lastIdx = 0;
        m_prevIdxForProj[0] = -1; // Should not be used.

        for (int r = 1; r < RKC::numStages; ++r) {
            m_prevIdxForProj[r] = lastIdx;

            if (r < RKC::numStages - 1) {
                if (RealCmp::neq(RKC::c[r + 1], RKC::c[r])) lastIdx = r;
            } else if (r == RKC::numStages - 1) {
                if (RealCmp::neq(1.0, RKC::c[r])) lastIdx = r;
            }
        }
        m_prevIdxForProj[RKC::numStages] = lastIdx;
    }

    // m_zetaE, m_zetaI
    {
        constexpr int r = 0;
        for (int j = 0; j < RKC::numStages; ++j) {
            m_zetaE[r][j] = RKC::aE[r][j];
            m_zetaI[r][j] = RKC::aI[r][j];
        }
    }
    for (int r = 1; r < RKC::numStages; ++r) {
        const int rp = m_prevIdxForProj[r];
        for (int j = 0; j <= rp; ++j) {
            m_zetaE[r][j] = RKC::aE[r][j] - RKC::aE[rp][j];
            m_zetaI[r][j] = RKC::aI[r][j] - RKC::aI[rp][j];
        }
        for (int j = rp + 1; j < RKC::numStages; ++j) {
            m_zetaE[r][j] = RKC::aE[r][j];
            m_zetaI[r][j] = RKC::aI[r][j];
        }
    }
    {
        const int s = RKC::numStages;
        const int rp = m_prevIdxForProj[s];
        for (int j = 0; j <= rp; ++j) {
            m_zetaE[s][j] = RKC::bE[j] - RKC::aE[rp][j];
            m_zetaI[s][j] = RKC::bI[j] - RKC::aI[rp][j];
        }
        for (int j = rp + 1; j < s; ++j) {
            m_zetaE[s][j] = RKC::bE[j];
            m_zetaI[s][j] = RKC::bI[j];
        }
    }

    // m_needsAssembly, m_needsAssemblyProjection
    {
        // If the last row of RKC::a equals RKC::b, then we don't need
        // to do anything!
        m_needsAssembly = false;
        for (int r = 0; r < RKC::numStages; ++r) {
            if (RealCmp::neq(RKC::aE[RKC::numStages - 1][r], RKC::bE[r]))
                m_needsAssembly = true;
            if (RealCmp::neq(RKC::aI[RKC::numStages - 1][r], RKC::bI[r]))
                m_needsAssembly = true;
        }
        // CH_verify(!m_needsAssembly || RealCmp::eq(RKC::c[RKC::numStages - 1], 1.0));

        m_needsAssemblyProjection = false;
        if (m_needsAssembly) {
            const int s = RKC::numStages;
            for (int j = 0; j < s; ++j) {
                m_needsAssemblyProjection |= !RealCmp::isZero(m_zetaE[s][j]);
                m_needsAssemblyProjection |= !RealCmp::isZero(m_zetaI[s][j]);
            }
        }
    }
    TODONOTE("m_needsAssemblyProjection is hardcoded to true.");
    m_needsAssembly = true;
    m_needsAssemblyProjection = true;

    CH_verify(!s_removeDivFromK || m_needsAssemblyProjection);

    // Debugging code...
    if (0) {
        const Real a_oldTime = 0.0;
        const Real a_dt = 1.0;

        pout() << "Standard RK coeffs:" << Format::indent() << '\n';
        for (int i = 1; i < RKC::numStages; ++i) {
            const Real stageTime    = a_oldTime + a_dt * RKC::c[i];
            const Real gammaDt      = a_dt * RKC::aI[i][i];
            const Real projDt       = a_dt * RKC::c[i];

            pout() << Format::fixed;
            pout() << "RK Stage " << i << ":\n";
            pout() << Format::indent() << flush;

            POUT(stageTime);
            POUT(gammaDt);
            POUT(projDt);

            pout() << "aE[i][j] = [ ";
            Real sum = 0.0;
            for (int j = 0; j < i; ++j) {
                pout() << RKC::aE[i][j] << " ";
                sum += RKC::aE[i][j];
            }
            pout() << "]\n";
            pout() << "Sum_j[aE[i][j]] = " << sum << '\n';

            pout() << "aI[i][j] = [ ";
            sum = 0.0;
            for (int j = 0; j <= i; ++j) {
                pout() << RKC::aI[i][j] << " ";
                sum += RKC::aI[i][j];
            }
            pout() << "]\n";
            pout() << "Sum_j[aI[i][j]] = " << sum << '\n';

            pout() << Format::unindent << flush;
        }
        pout() << Format::unindent << std::endl;

        pout() << "Incremental RK coeffs:" << Format::indent() << '\n';
        int prevIdxForProj = 0;
        for (int i = 1; i < RKC::numStages; ++i) {
            const Real stageTime    = a_oldTime + a_dt * RKC::c[i];
            const Real gammaDt      = a_dt * RKC::aI[i][i];
            const Real projDt       = a_dt * (RKC::c[i] - RKC::c[prevIdxForProj]);

            pout() << Format::fixed;
            pout() << "RK Stage " << i << ":\n";
            pout() << Format::indent() << flush;

            POUT(prevIdxForProj);
            POUT(stageTime);
            POUT(gammaDt);
            POUT(projDt);

            pout() << "zetaE[i][j] = [ ";
            for (int j = 0; j <= prevIdxForProj; ++j) {
                Real zetaE = RKC::aE[i][j] - RKC::aE[i - 1][j];
                pout() << zetaE << ' ';
            }

            for (int j = prevIdxForProj + 1; j < i; ++j) {
                pout() << "aE ";
            }
            pout() << "]\n";

            pout() << "zetaI[i][j] = [ ";
            for (int j = 0; j <= prevIdxForProj; ++j) {
                Real zetaI = RKC::aI[i][j] - RKC::aI[i - 1][j];
                pout() << zetaI << ' ';
            }

            for (int j = prevIdxForProj + 1; j <= i; ++j) {
                pout() << "aI ";
            }
            pout() << "]\n";


            // If c[i+1] != c[i] or if c[last i] != 1, then the next stage will be
            // an increment to this stage.
            if (i < RKC::numStages - 1) {
                if (abs(RKC::c[i + 1] - RKC::c[i]) > smallReal) prevIdxForProj = i;
            } else if (i == RKC::numStages - 1) {
                if (abs(1.0 - RKC::c[i]) > smallReal) prevIdxForProj = i;
            }

            pout() << Format::unindent << flush;
        }
        do {
            // If the last row of RKC::a equals RKC::b, then we don't need
            // to do anything!
            bool lastARowEqualsB = true;
            for (int i = 0; i < RKC::numStages; ++i) {
                if (abs(RKC::aE[RKC::numStages - 1][i] - RKC::bE[i]) > smallReal)
                    lastARowEqualsB = false;
                if (abs(RKC::aI[RKC::numStages - 1][i] - RKC::bI[i]) > smallReal)
                    lastARowEqualsB = false;
            }
            if (lastARowEqualsB) {
                CH_assert(abs(RKC::c[RKC::numStages - 1] - 1.0) < smallReal);
                break;
            }

            const int i = RKC::numStages;

            bool needsAssemblyProjection = false;
            for (int j = 0; j <= prevIdxForProj; ++j) {
                Real zetaE = RKC::bE[j] - RKC::aE[prevIdxForProj][j];
                if (abs(zetaE) > smallReal) {
                    pout() << "zetaE[" << j << "] = " << zetaE << '\n';
                    needsAssemblyProjection = true;
                }

                Real zetaI = RKC::bI[j] - RKC::aI[prevIdxForProj][j];
                if (abs(zetaI) > smallReal) {
                    pout() << "zetaI[" << j << "] = " << zetaI << '\n';
                    needsAssemblyProjection = true;
                }
            }
            for (int j = prevIdxForProj + 1; j < i; ++j) {
                if (abs(RKC::bE[j]) > smallReal) {
                    pout() << "bE[" << j << "] = " << RKC::bE[j] << '\n';
                    needsAssemblyProjection = true;
                }

                if (abs(RKC::bI[j]) > smallReal) {
                    pout() << "bI[" << j << "] = " << RKC::bI[j] << '\n';
                    needsAssemblyProjection = true;
                }
            }
            if (needsAssemblyProjection) {
                pout() << Format::fixed;
                pout() << "Assembly:\n";
                pout() << Format::indent() << flush;


                pout() << Format::unindent << flush;
            }
        } while (0);



        pout() << "New incremental RK coeffs:" << Format::indent() << '\n';
        for (int i = 1; i < RKC::numStages; ++i) {
            const Real stageTime    = a_oldTime + a_dt * RKC::c[i];
            const Real gammaDt      = a_dt * m_zetaI[i][i];
            const Real projDt       = a_dt * (RKC::c[i] - RKC::c[m_prevIdxForProj[i]]);

            pout() << Format::fixed;
            pout() << "RK Stage " << i << ":\n";
            pout() << Format::indent() << flush;

            POUT(m_prevIdxForProj[i]);
            POUT(stageTime);
            POUT(gammaDt);
            POUT(projDt);

            pout() << "zetaE[i][j] = [ ";
            for (int j = 0; j <= m_prevIdxForProj[i]; ++j) {
                pout() << m_zetaE[i][j] << ' ';
            }

            for (int j = m_prevIdxForProj[i] + 1; j < i; ++j) {
                pout() << "aE ";
            }
            pout() << "]\n";

            pout() << "zetaI[i][j] = [ ";
            for (int j = 0; j <= m_prevIdxForProj[i]; ++j) {
                pout() << m_zetaI[i][j] << ' ';
            }

            for (int j = m_prevIdxForProj[i] + 1; j <= i; ++j) {
                pout() << "aI ";
            }
            pout() << "]\n";

            pout() << Format::unindent << flush;
        }
        {
            const int s = RKC::numStages;
            pout() << "zetaE = (" << m_zetaE[s][0];
            for (int j = 1; j < s; ++j) {
                pout() << ' ' << m_zetaE[s][j];
            }
            pout() << ")\n";

            pout() << "zetaI = (" << m_zetaI[s][0];
            for (int j = 1; j < s; ++j) {
                pout() << ' ' << m_zetaI[s][j];
            }
            pout() << ")\n";
        }

        pout() << Format::unindent << std::endl;
    }
}


// -----------------------------------------------------------------------------
template <class RKC>
PARK<RKC>::~PARK()
{
}


#ifdef PARK__USE_INCREMENTAL_RK
// -----------------------------------------------------------------------------
template <class RKC>
void
PARK<RKC>::advance(LevelData<FluxBox>&   a_vel,
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

    DataIterator dit = m_grids.dataIterator();

    // Initialization. This sets Q[0].
    this->setOldQ(a_vel, a_p, a_q, a_oldTime);
    m_dt = a_dt;

    // ------ Main loop ------

    // Stage 1 (r = 0)
    {
        const int  r            = 0;
        const Real stageTime    = a_oldTime;
        const Real refluxScaleE = a_dt * m_zetaE[RKC::numStages][r];
        const Real refluxScaleI = a_dt * m_zetaI[RKC::numStages][r];

        pout() << Format::fixed;
        pout() << "RK Stage " << r << ":\n";
        pout() << Format::indent() << flush;

        // Stage 1 data already saved.
        // Evaluate this stage's forces.
        a_rhsPtr->setExplicitRHS(m_kvelE[r], m_kqE[r], a_vel, a_p, a_q, stageTime, refluxScaleE);
        a_rhsPtr->setImplicitRHS(m_kvelI[r], m_kqI[r], a_vel, a_p, a_q, stageTime, refluxScaleI);

        pout() << Format::unindent << flush;
    }

    // Stages 2 to RKC::numStages  (r = 1 to RKC::numStages - 1)
    for (int r = 1; r < RKC::numStages; ++r) {
        const int  rp            = m_prevIdxForProj[r];
        const Real stageTime     = a_oldTime + a_dt * RKC::c[r];
        const Real gammaDt       = a_dt * m_zetaI[r][r];
        const Real projDt        = a_dt * (RKC::c[r] - RKC::c[rp]);
        const Real refluxScaleE  = a_dt * m_zetaE[RKC::numStages][r];
        const Real refluxScaleI  = a_dt * m_zetaI[RKC::numStages][r];

        pout() << Format::fixed;
        pout() << "RK Stage " << r << ":\n";
        pout() << Format::indent() << flush;

        // Construct this stage's state.
        if (rp != r - 1) {
            this->copy(a_vel, m_vel[rp]);
            this->copy(  a_q, m_q[rp]);
            this->copy(  a_p, m_p[rp]);
        }

        for (int j = 0; j < r; ++j) {
            if (!RealCmp::isZero(m_zetaE[r][j])) {
                this->plus(a_vel, a_dt * m_zetaE[r][j], m_kvelE[j]);
                this->plus(  a_q, a_dt * m_zetaE[r][j],   m_kqE[j]);
            }
            if (!RealCmp::isZero(m_zetaI[r][j])) {
                this->plus(a_vel, a_dt * m_zetaI[r][j], m_kvelI[j]);
                this->plus(  a_q, a_dt * m_zetaI[r][j],   m_kqI[j]);
            }
        }

        a_rhsPtr->projectPredict(a_vel, a_p, stageTime, projDt);
        a_rhsPtr->solveImplicit(a_vel, a_q, gammaDt, stageTime, 0.0);  // refluxing scale???
        a_rhsPtr->projectCorrect(a_vel, a_p, stageTime, projDt, gammaDt);

        // Extrapolate p_r to t_r.
        if (s_extrapStagePressure) {
            const Real theta = (RKC::c[r] - RKC::c[rp]) / (RKC::c[r] + RKC::c[rp]);
            this->axby(1.0 + theta, a_p, -theta, m_p[0]);

            // this->axby(2.0, a_p, -1.0, m_p[rp]);
        }

        // Evaluate this stage's forces.
        a_rhsPtr->setExplicitRHS(m_kvelE[r], m_kqE[r], a_vel, a_p, a_q, stageTime, refluxScaleE);
        a_rhsPtr->setImplicitRHS(m_kvelI[r], m_kqI[r], a_vel, a_p, a_q, stageTime, refluxScaleI);

        // Save stage data (except last stage, which will hold Q^{n+1})
        this->copy(m_vel[r], a_vel); // Technically, not needed at last stage.
        this->copy(m_q[r], a_q);     // Technically, not needed at last stage.
        this->copy(m_p[r], a_p);

        pout() << Format::unindent << flush;
    }

    // Assembly of new state.
    if (m_needsAssembly) {
        const int  s         = RKC::numStages;
        const int  rp        = m_prevIdxForProj[s];
        const Real stageTime = a_oldTime + a_dt;
        const Real projDt    = a_dt * (1.0 - RKC::c[rp]);

        // Assemble starting from the last state with c != 1.
        if (rp != s - 1) {
            this->copy(a_vel, m_vel[rp]);
            this->copy(  a_q,   m_q[rp]);
            this->copy(  a_p,   m_p[rp]);
        }
        for (int j = 0; j < s; ++j) {
            if (abs(m_zetaE[s][j]) > smallReal) {
                this->plus(a_vel, a_dt * m_zetaE[s][j], m_kvelE[j]);
                this->plus(  a_q, a_dt * m_zetaE[s][j],   m_kqE[j]);
            }
            if (abs(m_zetaI[s][j]) > smallReal) {
                this->plus(a_vel, a_dt * m_zetaI[s][j], m_kvelI[j]);
                this->plus(  a_q, a_dt * m_zetaI[s][j],   m_kqI[j]);
            }
        }

        // // Assemble starting from initial state.
        // if (rp != s - 1) {
        //     this->copy(a_vel, m_vel[0]);
        //     this->copy(  a_q,   m_q[0]);
        //     this->copy(  a_p,   m_p[0]);
        // }
        // for (int j = 0; j < s; ++j) {
        //     if (abs(RKC::bE[j]) > smallReal) {
        //         this->plus(a_vel, a_dt * RKC::bE[j], m_kvelE[j]);
        //         this->plus(  a_q, a_dt * RKC::bE[j],   m_kqE[j]);
        //     }
        //     if (abs(RKC::bI[j]) > smallReal) {
        //         this->plus(a_vel, a_dt * RKC::bI[j], m_kvelI[j]);
        //         this->plus(  a_q, a_dt * RKC::bI[j],   m_kqI[j]);
        //     }
        // }

        if (m_needsAssemblyProjection) {
            pout() << Format::fixed;
            pout() << "Assembly:\n";
            pout() << Format::indent() << flush;

            a_rhsPtr->projectPredict(a_vel, a_p, stageTime, projDt);
            a_rhsPtr->projectCorrect(a_vel, a_p, stageTime, projDt);

            // Extrapolate p^{n+1} to t^{n+1}.
            if (s_extrapAssemblyPressure) {
                const Real theta = (1.0 - RKC::c[rp]) / (1.0 + RKC::c[rp]);
                this->axby(1.0 + theta, a_p, -theta, m_p[0]);

                // this->axby(2.0, a_p, -1.0, m_p[rp]);

                // this->axby(2.0, a_p, -1.0, m_p[0]);
            }

            pout() << Format::unindent << flush;
        }
    }

    // Remove divergence from kE.
    if (s_removeDivFromK) {
        for (int r = 1; r < RKC::numStages; ++r) {
            this->copy(m_kvelE[r-1],  m_vel[r]);
            this->axby(1.0 / (a_dt * m_zetaE[r][r-1]), m_kvelE[r-1], -1.0 / (a_dt * m_zetaE[r][r-1]), m_vel[0]);
            for (int i = 0; i <= r - 2; ++i) {
                this->plus(m_kvelE[r-1], -m_zetaE[r][i] / m_zetaE[r][r-1], m_kvelE[i]);
            }
        }
        {
            const int s = RKC::numStages - 1;
            this->copy(m_kvelE[s], a_vel);
            this->axby(1.0 / (a_dt * m_zetaE[s+1][s]), m_kvelE[s], -1.0 / (a_dt * m_zetaE[s+1][s]), m_vel[0]);
            for (int i = 0; i <= s - 1; ++i) {
                this->plus(m_kvelE[s], -m_zetaE[s+1][i] / m_zetaE[s+1][s], m_kvelE[i]);
            }
        }

        for (int r = 0; r < RKC::numStages; ++r) {
            setValLevel(m_kvelI[r], 0.0);
        }
    }

    // Save new state.
    this->copy(m_vel[RKC::numStages - 1], a_vel);
    this->copy(m_q[RKC::numStages - 1], a_q);
    this->copy(m_p[RKC::numStages - 1], a_p);

    m_newTime         = m_oldTime + m_dt;
    m_Q0Ready         = true;
    m_interpDataReady = true;

    // Allow user to perform postStep operations (e.g., preparing plots)
    a_rhsPtr->postStep(a_vel, a_p, a_q, m_newTime);


    // Debugging code...
    if (this->hasEmbeddedScheme()) {
        // Assemble error approximation.
        LevelData<FluxBox>   dvel(m_grids, m_velNumComps, m_velGhostVect);
        LevelData<FArrayBox>   dq(m_grids,   m_qNumComps,   m_qGhostVect);
        setValLevel(dvel, 0.0);
        setValLevel(dq, 0.0);

        for (int j = 0; j < RKC::numStages; ++j) {
            this->plus(dvel, m_dt * (RKC::bE[j] - RKC::bEhat[j]), m_kvelE[j]);
            this->plus(  dq, m_dt * (RKC::bE[j] - RKC::bEhat[j]),   m_kqE[j]);

            this->plus(dvel, m_dt * (RKC::bI[j] - RKC::bIhat[j]), m_kvelI[j]);
            this->plus(  dq, m_dt * (RKC::bI[j] - RKC::bIhat[j]),   m_kqI[j]);
        }

        const RealVect velErrNorm0 = Analysis::pNorm(dvel, 2);
        const Real       qErrNorm0 = Analysis::pNorm(  dq, 2);

        pout() << Format::scientific;
        POUT(velErrNorm0);
        POUT(qErrNorm0);

        // static std::vector<Real> verr0(0);
        // verr0.push_back(std::max(velErrNorm0[0], velErrNorm0[1]));
        // POUT(verr0);

        LevelData<FArrayBox> p(m_grids, 1, IntVect::Unit);
        setValLevel(p, 0.0);
        a_rhsPtr->projectPredict(dvel, p, m_newTime, m_dt);
        a_rhsPtr->projectCorrect(dvel, p, m_newTime, m_dt);

        const RealVect velErrNorm1 = Analysis::pNorm(dvel, 2);
        POUT(velErrNorm1);

        // static std::vector<Real> verr1(0);
        // verr1.push_back(std::max(velErrNorm1[0], velErrNorm1[1]));
        // POUT(verr1);
    }

    if (0 && this->hasEmbeddedScheme()) {
        // E_projected
        LevelData<FluxBox> vel(m_grids, 1, m_velGhostVect);
        this->copy(vel, a_vel);
        LevelData<FArrayBox> E_projected(m_grids, 1);
        for (DataIterator dit(m_grids); dit.ok(); ++dit) {
            for (BoxIterator bit(m_grids[dit]); bit.ok(); ++bit) {
                const IntVect cc = bit();
                const RealVect ccvel(D_DECL(
                    0.5 * (vel[dit][0](cc) + vel[dit][0](cc + BASISV(0))),
                    0.5 * (vel[dit][1](cc) + vel[dit][1](cc + BASISV(1))),
                    0.5 * (vel[dit][2](cc) + vel[dit][2](cc + BASISV(2)))
                ));
                E_projected[dit](cc) = 0.5 * ccvel.vectorLengthSq();
            }
        }
        const Real E_projected_scale = Analysis::pNorm(E_projected, 2);
        POUT(E_projected_scale);


        // E
        this->copy(vel, m_vel[0]);
        for (int j = 0; j < RKC::numStages; ++j) {
            this->plus(vel, m_dt * RKC::bE[j], m_kvelE[j]);
            this->plus(vel, m_dt * RKC::bI[j], m_kvelI[j]);
        }
        LevelData<FArrayBox> E(m_grids, 1);
        for (DataIterator dit(m_grids); dit.ok(); ++dit) {
            for (BoxIterator bit(m_grids[dit]); bit.ok(); ++bit) {
                const IntVect cc = bit();
                const RealVect ccvel(D_DECL(
                    0.5 * (vel[dit][0](cc) + vel[dit][0](cc + BASISV(0))),
                    0.5 * (vel[dit][1](cc) + vel[dit][1](cc + BASISV(1))),
                    0.5 * (vel[dit][2](cc) + vel[dit][2](cc + BASISV(2)))
                ));
                E[dit](cc) = 0.5 * ccvel.vectorLengthSq();
            }
        }
        const Real E_scale = Analysis::pNorm(E, 2);
        POUT(E_scale);


        // Ehat
        this->copy(vel, m_vel[0]);
        for (int j = 0; j < RKC::numStages; ++j) {
            this->plus(vel, m_dt * RKC::bEhat[j], m_kvelE[j]);
            this->plus(vel, m_dt * RKC::bIhat[j], m_kvelI[j]);
        }
        LevelData<FArrayBox> Ehat(m_grids, 1);
        for (DataIterator dit(m_grids); dit.ok(); ++dit) {
            for (BoxIterator bit(m_grids[dit]); bit.ok(); ++bit) {
                const IntVect cc = bit();
                const RealVect ccvel(D_DECL(
                    0.5 * (vel[dit][0](cc) + vel[dit][0](cc + BASISV(0))),
                    0.5 * (vel[dit][1](cc) + vel[dit][1](cc + BASISV(1))),
                    0.5 * (vel[dit][2](cc) + vel[dit][2](cc + BASISV(2)))
                ));
                Ehat[dit](cc) = 0.5 * ccvel.vectorLengthSq();
            }
        }
        const Real Ehat_scale = Analysis::pNorm(Ehat, 2);
        POUT(Ehat_scale);


        // Ehat_semiprojected
        LevelData<FArrayBox> p(m_grids, 1, m_pGhostVect);
        this->copy(p, a_p);
        a_rhsPtr->projectPredict(vel, p, m_newTime, m_dt);
        LevelData<FArrayBox> Ehat_semiprojected(m_grids, 1);
        for (DataIterator dit(m_grids); dit.ok(); ++dit) {
            for (BoxIterator bit(m_grids[dit]); bit.ok(); ++bit) {
                const IntVect cc = bit();
                const RealVect ccvel(D_DECL(
                    0.5 * (vel[dit][0](cc) + vel[dit][0](cc + BASISV(0))),
                    0.5 * (vel[dit][1](cc) + vel[dit][1](cc + BASISV(1))),
                    0.5 * (vel[dit][2](cc) + vel[dit][2](cc + BASISV(2)))
                ));
                Ehat_semiprojected[dit](cc) = 0.5 * ccvel.vectorLengthSq();
            }
        }
        const Real Ehat_semiprojected_scale = Analysis::pNorm(Ehat_semiprojected, 2);
        POUT(Ehat_semiprojected_scale);


        // Ehat_projected
        a_rhsPtr->projectCorrect(vel, p, m_newTime, m_dt);
        LevelData<FArrayBox> Ehat_projected(m_grids, 1);
        for (DataIterator dit(m_grids); dit.ok(); ++dit) {
            for (BoxIterator bit(m_grids[dit]); bit.ok(); ++bit) {
                const IntVect cc = bit();
                const RealVect ccvel(D_DECL(
                    0.5 * (vel[dit][0](cc) + vel[dit][0](cc + BASISV(0))),
                    0.5 * (vel[dit][1](cc) + vel[dit][1](cc + BASISV(1))),
                    0.5 * (vel[dit][2](cc) + vel[dit][2](cc + BASISV(2)))
                ));
                Ehat_projected[dit](cc) = 0.5 * ccvel.vectorLengthSq();
            }
        }
        const Real Ehat_projected_scale = Analysis::pNorm(Ehat_projected, 2);
        POUT(Ehat_projected_scale);


        // E - Ehat
        LevelData<FArrayBox> dE(m_grids, 1);
        this->copy(dE, E);
        this->plus(dE, -1.0, Ehat);
        pout() << "|E - E_hat| = " << Analysis::pNorm(dE, 2) << '\n';

        // E_projected - Ehat_semiprojected
        this->copy(dE, E_projected);
        this->plus(dE, -1.0, Ehat_semiprojected);
        pout() << "|E_projected - Ehat_semiprojected| = " << Analysis::pNorm(dE, 2) << '\n';

        // E_projected - Ehat_projected
        this->copy(dE, E_projected);
        this->plus(dE, -1.0, Ehat_projected);
        pout() << "|E_projected - Ehat_projected| = " << Analysis::pNorm(dE, 2) << '\n';

        // E_projected - E
        this->copy(dE, E_projected);
        this->plus(dE, -1.0, E);
        pout() << "|E_projected - E| = " << Analysis::pNorm(dE, 2) << '\n';

        // Ehat_projected - Ehat
        this->copy(dE, Ehat_projected);
        this->plus(dE, -1.0, Ehat);
        pout() << "|Ehat_projected - Ehat| = " << Analysis::pNorm(dE, 2) << '\n';
    }
}

#else // Use Standard RK...

// -----------------------------------------------------------------------------
template <class RKC>
void
PARK<RKC>::advance(LevelData<FluxBox>&   a_vel,
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

    DataIterator dit = m_grids.dataIterator();

    // Initialization. This sets Q[0].
    this->setOldQ(a_vel, a_p, a_q, a_oldTime);
    m_dt = a_dt;

    // ------ Main loop ------

    // Stage 1 (i = 0)
    {
        constexpr int  i        = 0;
        const Real stageTime    = a_oldTime;
        const Real refluxScaleE = a_dt * RKC::bE[i];
        const Real refluxScaleI = a_dt * RKC::bI[i];

        pout() << Format::fixed;
        pout() << "RK Stage " << i << ":\n";
        pout() << Format::indent() << flush;

        // Stage 1 data already saved.
        // Evaluate this stage's forces.
        a_rhsPtr->setExplicitRHS(m_kvelE[i], m_kqE[i], m_vel[i], m_p[i], m_q[i], stageTime, refluxScaleE);
        a_rhsPtr->setImplicitRHS(m_kvelI[i], m_kqI[i], m_vel[i], m_p[i], m_q[i], stageTime, refluxScaleI);

        pout() << Format::unindent << flush;
    }

    // Stages 2 to RKC::numStages  (i = 1 to RKC::numStages - 1)
    for (int i = 1; i < RKC::numStages; ++i) {
        const Real stageTime    = a_oldTime + a_dt * RKC::c[i];
        const Real gammaDt      = a_dt * RKC::aI[i][i];
        const Real projDt       = a_dt * RKC::c[i];
        const Real refluxScaleE = a_dt * RKC::bE[i];
        const Real refluxScaleI = a_dt * RKC::bI[i];

        pout() << Format::fixed;
        pout() << "RK Stage " << i << ":\n";
        pout() << Format::indent() << flush;

        // Construct this stage's state.
        this->copy(m_vel[i], m_vel[0]);
        this->copy(m_q[i], m_q[0]);
        this->copy(m_p[i], m_p[0]);

        for (int j = 0; j < i; ++j) {
            if (abs(RKC::aE[i][j]) > smallReal) {
                this->plus(m_vel[i], a_dt * RKC::aE[i][j], m_kvelE[j]);
                this->plus(  m_q[i], a_dt * RKC::aE[i][j],   m_kqE[j]);
            }

            if (abs(RKC::aI[i][j]) > smallReal) {
                this->plus(m_vel[i], a_dt * RKC::aI[i][j], m_kvelI[j]);
                this->plus(  m_q[i], a_dt * RKC::aI[i][j],   m_kqI[j]);
            }
        }

        a_rhsPtr->projectPredict(m_vel[i], m_p[i], stageTime, projDt);
        a_rhsPtr->solveImplicit(m_vel[i], m_q[i], gammaDt, stageTime, 0.0);  // refluxing scale???
        a_rhsPtr->projectCorrect(m_vel[i], m_p[i], stageTime, projDt);


        // Extrapolate p to stage time.
        if (s_extrapStagePressure) {
            this->axby(2.0, a_p, -1.0, m_p[0]);
        }

        // Evaluate this stage's forces.
        a_rhsPtr->setExplicitRHS(m_kvelE[i], m_kqE[i], m_vel[i], m_p[i], m_q[i], stageTime, refluxScaleE);
        a_rhsPtr->setImplicitRHS(m_kvelI[i], m_kqI[i], m_vel[i], m_p[i], m_q[i], stageTime, refluxScaleI);

        pout() << Format::unindent << flush;
    }

    // Assembly of new state.
    do {
        // // If the last row of RKC::a equals RKC::b, then we don't need
        // // to do anything!
        // bool lastARowEqualsB = true;
        // for (int i = 0; i < RKC::numStages; ++i) {
        //     if (abs(RKC::aE[RKC::numStages - 1][i] - RKC::bE[i]) > smallReal)
        //         lastARowEqualsB = false;
        //     if (abs(RKC::aI[RKC::numStages - 1][i] - RKC::bI[i]) > smallReal)
        //         lastARowEqualsB = false;
        // }
        // if (lastARowEqualsB) {
        //     CH_assert(abs(RKC::c[RKC::numStages - 1] - 1.0) < smallReal);
        //     break;
        // }

        const Real stageTime = a_oldTime + a_dt;
        const Real projDt    = a_dt;

        // a_vel, etc still contain the initial state...
        for (int j = 0; j < RKC::numStages; ++j) {
            if (abs(RKC::bE[j]) > smallReal) {
                this->plus(a_vel, a_dt * RKC::bE[j], m_kvelE[j]);
                this->plus(  a_q, a_dt * RKC::bE[j],   m_kqE[j]);
            }

            if (abs(RKC::bI[j]) > smallReal) {
                this->plus(a_vel, a_dt * RKC::bI[j], m_kvelI[j]);
                this->plus(  a_q, a_dt * RKC::bI[j],   m_kqI[j]);
            }
        }

        pout() << Format::fixed;
        pout() << "Assembly:\n";
        pout() << Format::indent() << flush;

        a_rhsPtr->projectPredict(a_vel, a_p, stageTime, projDt);
        a_rhsPtr->projectCorrect(a_vel, a_p, stageTime, projDt);

        // Extrapolate p^{n+1} to t^{n+1}.
        if (s_extrapAssemblyPressure) {
            this->axby(2.0, a_p, -1.0, m_p[0]);
        }

        pout() << Format::unindent << flush;
    } while (0);


    // Remove divergence from kE.
    if (s_removeDivFromK) {
        for (int r = 1; r < RKC::numStages; ++r) {
            this->copy(m_kvelE[r-1],  m_vel[r]);
            this->axby(1.0 / (a_dt * RKC::aE[r][r-1]), m_kvelE[r-1], -1.0 / (a_dt * RKC::aE[r][r-1]), m_vel[0]);
            for (int i = 0; i <= r - 2; ++i) {
                this->plus(m_kvelE[r-1], -RKC::aE[r][i] / RKC::aE[r][r-1], m_kvelE[i]);
            }
            setValLevel(m_kvelI[r-1], 0.0);
        }
        {
            const int s = RKC::numStages - 1;
            this->copy(m_kvelE[s], a_vel);
            this->axby(1.0 / (a_dt * RKC::bE[s]), m_kvelE[s], -1.0 / (a_dt * RKC::bE[s]), m_vel[0]);
            for (int i = 0; i <= s - 1; ++i) {
                this->plus(m_kvelE[s], -RKC::bE[i] / RKC::bE[s], m_kvelE[i]);
            }
            setValLevel(m_kvelI[s], 0.0);
        }
    }

    // Save new state.
    this->copy(m_vel[RKC::numStages - 1], a_vel);
    this->copy(m_q[RKC::numStages - 1], a_q);
    this->copy(m_p[RKC::numStages - 1], a_p);


    m_newTime         = m_oldTime + m_dt;
    m_Q0Ready         = true;
    m_interpDataReady = true;

    // Allow user to perform postStep operations (e.g., preparing plots)
    a_rhsPtr->postStep(a_vel, a_p, a_q, m_newTime);

    if (0 & this->hasEmbeddedScheme()) {
        // Assemble error approximation.
        LevelData<FluxBox>   dvel(m_grids, m_velNumComps, m_velGhostVect);
        LevelData<FArrayBox>   dq(m_grids,   m_qNumComps,   m_qGhostVect);
        setValLevel(dvel, 0.0);
        setValLevel(dq, 0.0);

        for (int j = 0; j < RKC::numStages; ++j) {
            this->plus(dvel, m_dt * (RKC::bE[j] - RKC::bEhat[j]), m_kvelE[j]);
            this->plus(  dq, m_dt * (RKC::bE[j] - RKC::bEhat[j]),   m_kqE[j]);

            this->plus(dvel, m_dt * (RKC::bI[j] - RKC::bIhat[j]), m_kvelI[j]);
            this->plus(  dq, m_dt * (RKC::bI[j] - RKC::bIhat[j]),   m_kqI[j]);
        }

        const RealVect velErrNorm0 = Analysis::pNorm(dvel, 2);
        const Real       qErrNorm0 = Analysis::pNorm(  dq, 2);

        pout() << Format::scientific;
        POUT(velErrNorm0);
        POUT(qErrNorm0);

        // POUT(velErrNorm0 * 0.0159194573264);

        // static std::vector<Real> verr0(0);
        // verr0.push_back(std::max(velErrNorm0[0], velErrNorm0[1]));
        // POUT(verr0);


        // LevelData<FArrayBox> p(m_grids, 1, IntVect::Unit);
        // setValLevel(p, 0.0);
        // a_rhsPtr->projectPredict(dvel, p, m_newTime, m_dt);
        // a_rhsPtr->projectCorrect(dvel, p, m_newTime, m_dt);

        // const RealVect velErrNorm1 = Analysis::pNorm(dvel, 2);
        // POUT(velErrNorm1);

        // static std::vector<Real> verr1(0);
        // verr1.push_back(std::max(velErrNorm1[0], velErrNorm1[1]));
        // POUT(verr1);
    }
}
#endif


// -----------------------------------------------------------------------------
template <class RKC>
void
PARK<RKC>::FEadvance(LevelData<FluxBox>&   a_vel,
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

    this->setOldQ(a_vel, a_p, a_q, a_oldTime);
    m_dt = a_dt;


    // --- Begin stage ---
    pout() << Format::fixed << "FE timestep..."
           << "\n"
           << Format::indent() << flush;

    // Evaluate forces.
    a_rhsPtr->setExplicitRHS(
        m_kvelE[0], m_kqE[0], a_vel, a_p, a_q, a_oldTime, a_dt);
    a_rhsPtr->setImplicitRHS(
        m_kvelI[0], m_kqI[0], a_vel, a_p, a_q, a_oldTime, a_dt);

    nanCheck(m_kvelE[0]);
    nanCheck(m_kqE[0]);
    nanCheck(m_kvelI[0]);
    nanCheck(m_kqI[0]);

    // Update q (Ri)
    // Q += dt*(kQE + kQI)
    this->plus(a_vel, a_dt, m_kvelE[0]);
    this->plus(a_vel, a_dt, m_kvelI[0]);
    this->plus(a_q, a_dt, m_kqE[0]);
    this->plus(a_q, a_dt, m_kqI[0]);

    // Project.
    a_rhsPtr->projectPredict(a_vel, a_p, a_oldTime, a_dt);
    a_rhsPtr->projectCorrect(a_vel, a_p, a_oldTime, a_dt);

    nanCheck(a_vel);
    nanCheck(a_p);
    nanCheck(a_q);
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
template <class RKC>
Real
PARK<RKC>::controllerDt(const Real a_tol,
                        const bool a_useImplicit,
                        const bool a_useElementary,
                        const bool a_usePI,
                        const bool a_usePID) const
{
    if (!this->hasEmbeddedScheme()) return maxReal;
    if (!m_interpDataReady) return maxReal;

    // Sometimes, AMR will re-run a timestep to set up ICs or produce pressure
    // estimates. We don't want to push repeated values to the history.
    // WARNING: Use m_oldTime to tag history elements, not m_newTime!
    const Real timeTag = m_oldTime;
    const bool isNewTime = !m_errorHistory[0].hasValue || RealCmp::neq(m_errorHistory[0].time, timeTag);

    // |u - uhat| ~ c*dt^order.
    static const int order =
          (a_useImplicit
             ? std::min({s_orderE, s_orderEmbeddedE, s_orderI, s_orderEmbeddedI})
             : std::min({s_orderE, s_orderEmbeddedE}));

    const Real localError = [&]() {
        LevelData<FluxBox> dvel(m_grids, m_velNumComps, m_velGhostVect);

        setValLevel(dvel, 0.0);

        for (int j = 0; j < RKC::numStages; ++j) {
            this->plus(dvel, m_dt * (RKC::bE[j] - RKC::bEhat[j]), m_kvelE[j]);
            if (a_useImplicit) {
                this->plus(dvel, m_dt * (RKC::bI[j] - RKC::bIhat[j]), m_kvelI[j]);
            }
        }

        const RealVect norms = Analysis::pNorm(dvel, 2);
        return std::max({D_DECL(norms[0], norms[1], norms[2])});
    }();

    pout() << "|vel local temporal error| ~ " << localError << '\n';

    pushToHistory(
        ErrorHistoryElem(timeTag, m_dt, localError), m_errorHistory, isNewTime);

    Real retVal = maxReal;
    if (a_useElementary) {
        const Real eps      = a_tol / localError;
        const Real elemCtrl = std::pow(eps,  1.0 / Real(order + 1.0));
        retVal = std::min(retVal, elemCtrl);
    }

    if (a_usePI) {
        // Soderlind's PI11 controller
        const Real eps      = a_tol / localError;
        const Real r_ratio  = m_errorHistory[1].localError / localError;
        const Real dt_ratio = m_dt / m_errorHistory[1].dt;

        const Real piCtrl = std::pow(eps,     1.0 / Real(order + 1))
                          * std::pow(r_ratio, 1.0 / Real(order + 1))
                          * dt_ratio;

        retVal = std::min(retVal, piCtrl);
    }

    if (a_usePID) {
        Real pidCtrl;
        if (m_errorHistory[2].hasValue) {
            // From Kennedy & Carpenter
            const Real eps0 = a_tol / localError;
            const Real eps1 = a_tol / m_errorHistory[1].localError;
            const Real eps2 = a_tol / m_errorHistory[2].localError;
            pidCtrl = std::pow(eps0,  0.49 / Real(order + 1))
                    * std::pow(eps1, -0.34 / Real(order + 1))
                    * std::pow(eps2,  0.10 / Real(order + 1));

        } else if (m_errorHistory[1].hasValue) {
            // PI42 controller from Ranocha
            const Real eps0 = a_tol / localError;
            const Real eps1 = a_tol / m_errorHistory[1].localError;
            pidCtrl = std::pow(eps0,  0.60 / Real(order + 1))
                    * std::pow(eps1, -0.20 / Real(order + 1));

        } else {
            // Elementary controller
            const Real eps = a_tol / localError;
            pidCtrl = std::pow(eps,  1.0 / Real(order + 1.0));
        }
        retVal = std::min(retVal, pidCtrl);
    }

    // // Limiting
    // retVal = hardLimiter(retVal, m_controllerHistory[0]);
    // // retVal = softLimiter(retVal, m_controllerHistory[0]);
    // // retVal = oscillationLimiter(retVal, m_controllerHistory[0], m_controllerHistory[1]);
    // pushToHistory(retVal, m_controllerHistory, isNewTime);

    return retVal * m_dt;
}


// -----------------------------------------------------------------------------
template <class RKC>
Real
PARK<RKC>::velTimeInterp(LevelData<FluxBox>& a_vel,
                         const Real          a_time,
                         int                 a_srcComp,
                         int                 a_destComp,
                         int                 a_numComp) const
{
    // if (s_useDenseOutput && s_removeDivFromK) {
    if (s_useDenseOutput) {
        return this->velDenseOutput(a_vel, a_time, a_srcComp, a_destComp, a_numComp);
    } else {
        // return this->velLagrangeInterp(a_vel, a_time, a_srcComp, a_destComp, a_numComp);
        return this->velLinearInterp(a_vel, a_time, a_srcComp, a_destComp, a_numComp);
    }
}


// -----------------------------------------------------------------------------
template <class RKC>
Real
PARK<RKC>::pTimeInterp(LevelData<FArrayBox>& a_p,
                       const Real            a_time,
                       int                   a_srcComp,
                       int                   a_destComp,
                       int                   a_numComp) const
{
    // return this->pLagrangeInterp(a_p, a_time, a_srcComp, a_destComp, a_numComp);
    return this->pLinearInterp(a_p, a_time, a_srcComp, a_destComp, a_numComp);
}


// -----------------------------------------------------------------------------
template <class RKC>
Real
PARK<RKC>::qTimeInterp(LevelData<FArrayBox>& a_q,
                       const Real            a_time,
                       int                   a_srcComp,
                       int                   a_destComp,
                       int                   a_numComp) const
{
    if (s_useDenseOutput) {
        return this->qDenseOutput(a_q, a_time, a_srcComp, a_destComp, a_numComp);
    } else {
        // return this->qLagrangeInterp(a_q, a_time, a_srcComp, a_destComp, a_numComp);
        return this->qLinearInterp(a_q, a_time, a_srcComp, a_destComp, a_numComp);
    }
}


// -----------------------------------------------------------------------------
template <class RKC>
Real
PARK<RKC>::computeTheta(const Real a_time) const
{
    const Real theta = (a_time - m_oldTime) / m_dt;
    if (abs(theta - 0.0) < timeEps) return 0.0;
    if (abs(theta - 1.0) < timeEps) return 1.0;
    return theta;
}


// //
// -----------------------------------------------------------------------------
// // Interpolate kq = dq/dt in time using the dense output formula.
// // This will use the dense output formula and return how far along the
// // initial and final data the interpolation took place (0 = old time,
// // 1 = new time). If this returns a number outside of [0,1], then we
// // extrapolated and you need to worry about stability.
// //
// -----------------------------------------------------------------------------
// template <class RKC>
// Real
// PARK<RKC>::denseOutputQFirstDeriv(LevelData<FArrayBox>& a_kqtheta,
//                                   const Real            a_time,
//                                   int                   a_srcComp,
//                                   int                   a_destComp,
//                                   int                   a_numComp) const
// {
//     // Sanity checks
//     CH_assert(a_kqtheta.nComp() <= m_kqE[0].nComp());
//     CH_assert(a_kqtheta.getBoxes().compatible(m_grids));
//     CH_assert(a_kqtheta.getBoxes().physDomain().size() ==
//               m_grids.physDomain().size());


//     if (a_numComp == -1) a_numComp = m_kqE[0].nComp() - a_srcComp;

//     // About to do dense output interpolation. We need the kq's.
//     CH_assert(m_interpDataReady);

//     // Compute theta and if nearby, snap into 0 or 1 values to prevent
//     // extrapolation when not expected.
//     Real theta = (a_time - m_times[0]) / m_dt;

//     if (abs(theta) < timeEps) {
//         theta = 0.0;
//     } else if (abs(theta - 1.0) < timeEps) {
//         theta = 1.0;
//     }

//     // Compute powers of theta.
//     Real thetaPow[RKC::numDenseCoeffs];
//     for (int j = 0; j < RKC::numDenseCoeffs; ++j) {
//         thetaPow[j] = Real(j + 1) * pow(theta, j);
//     }

//     // Compute d/dtheta[bstar_i(theta)].
//     Real bstar[RKC::numStages];
//     for (int i = 0; i < RKC::numStages; ++i) {
//         bstar[i] = 0.0;
//         for (int j = 0; j < RKC::numDenseCoeffs; ++j) {
//             bstar[i] += RKC::bstar[j][i] * thetaPow[j];
//         }
//     }

//     // Perform interpolation.
//     setValLevel(a_kqtheta, 0.0);
//     for (int i = 0; i < RKC::numStages; ++i) {
//         plus(a_kqtheta, m_kqE[i], bstar[i], a_srcComp, a_destComp,
//         a_numComp); plus(a_kqtheta, m_kqI[i], bstar[i], a_srcComp,
//         a_destComp, a_numComp);
//     }

//     CH_assert(0.0 <= theta && theta <= 1.0);
//     return theta;
// }


// -----------------------------------------------------------------------------
template <class RKC>
Real
PARK<RKC>::velLinearInterp(LevelData<FluxBox>& a_vel,
                           const Real          a_time,
                           int                 a_srcComp,
                           int                 a_destComp,
                           int                 a_numComp) const
{
    CH_assert(a_numComp <= m_velNumComps);
    CH_assert(a_vel.getBoxes().compatible(m_grids));
    CH_assert(a_vel.getBoxes().physDomain().size() ==
              m_grids.physDomain().size());

    DataIterator              dit    = m_grids.dataIterator();
    const LevelData<FluxBox>& oldVel = m_vel[0];
    const LevelData<FluxBox>& newVel = m_vel[RKC::numStages - 1];

    if (a_numComp == -1) a_numComp = m_velNumComps - a_srcComp;

    // If we are at m_oldTime, just copy the old data.
    CH_assert(m_Q0Ready);
    if (abs(a_time - m_oldTime) < timeEps) {
        for (dit.reset(); dit.ok(); ++dit) {
            for (int dir = 0; dir < SpaceDim; ++dir) {
                a_vel[dit][dir].copy(
                    oldVel[dit][dir], a_srcComp, a_destComp, a_numComp);
            }
        }
        return 0.0;
    }

    // Linear interpolation...
    CH_assert(m_interpDataReady);

    // Compute theta and if nearby, snap into 0 or 1 values to prevent
    // extrapolation when not expected.
    const Real theta = this->computeTheta(a_time);

    if (abs(theta) < timeEps) {
        // A bit redundant...
        for (dit.reset(); dit.ok(); ++dit) {
            for (int dir = 0; dir < SpaceDim; ++dir) {
                a_vel[dit][dir].copy(
                    oldVel[dit][dir], a_srcComp, a_destComp, a_numComp);
            }
        }
        return 0.0;
    } else if (abs(theta - 1.0) < timeEps) {
        for (dit.reset(); dit.ok(); ++dit) {
            for (int dir = 0; dir < SpaceDim; ++dir) {
                a_vel[dit][dir].copy(
                    newVel[dit][dir], a_srcComp, a_destComp, a_numComp);
            }
        }
        return 1.0;
    }

    // Perform interpolation.
    for (dit.reset(); dit.ok(); ++dit) {
        for (int dir = 0; dir < SpaceDim; ++dir) {
            a_vel[dit][dir].setVal(0.0);

            a_vel[dit][dir].plus(oldVel[dit][dir],
                                 1.0 - theta,
                                 a_srcComp,
                                 a_destComp,
                                 a_numComp);

            a_vel[dit][dir].plus(
                newVel[dit][dir], theta, a_srcComp, a_destComp, a_numComp);
        }
    }

    CH_assert(0.0 <= theta && theta <= 1.0);
    return theta;
}


// -----------------------------------------------------------------------------
template <class RKC>
Real
PARK<RKC>::pLinearInterp(LevelData<FArrayBox>& a_p,
                         const Real            a_time,
                         int                   a_srcComp,
                         int                   a_destComp,
                         int                   a_numComp) const
{
    CH_assert(a_numComp <= m_pNumComps);
    CH_assert(a_p.getBoxes().compatible(m_grids));
    CH_assert(a_p.getBoxes().physDomain().size() ==
              m_grids.physDomain().size());

    DataIterator                dit  = m_grids.dataIterator();
    const LevelData<FArrayBox>& oldP = m_p[0];
    const LevelData<FArrayBox>& newP = m_p[RKC::numStages - 1];

    if (a_numComp == -1) a_numComp = m_pNumComps - a_srcComp;

    // If we are at m_oldTime, just copy the old data.
    CH_assert(m_Q0Ready);
    if (abs(a_time - m_oldTime) < timeEps) {
        for (dit.reset(); dit.ok(); ++dit) {
            a_p[dit].copy(oldP[dit], a_srcComp, a_destComp, a_numComp);
        }
        return 0.0;
    }

    // Linear interpolation...
    CH_assert(m_interpDataReady);

    // Compute theta and if nearby, snap into 0 or 1 values to prevent
    // extrapolation when not expected.
    const Real theta = this->computeTheta(a_time);

    if (abs(theta) < timeEps) {
        // A bit redundant...
        for (dit.reset(); dit.ok(); ++dit) {
            a_p[dit].copy(oldP[dit], a_srcComp, a_destComp, a_numComp);
        }
        return 0.0;
    } else if (abs(theta - 1.0) < timeEps) {
        for (dit.reset(); dit.ok(); ++dit) {
            a_p[dit].copy(newP[dit], a_srcComp, a_destComp, a_numComp);
        }
        return 1.0;
    }

    // Perform interpolation.
    for (dit.reset(); dit.ok(); ++dit) {
        a_p[dit].setVal(0.0);
        a_p[dit].plus(oldP[dit], 1.0 - theta, a_srcComp, a_destComp, a_numComp);
        a_p[dit].plus(newP[dit], theta, a_srcComp, a_destComp, a_numComp);
    }

    CH_assert(0.0 <= theta && theta <= 1.0);
    return theta;
}


// -----------------------------------------------------------------------------
template <class RKC>
Real
PARK<RKC>::qLinearInterp(LevelData<FArrayBox>& a_q,
                         const Real            a_time,
                         int                   a_srcComp,
                         int                   a_destComp,
                         int                   a_numComp) const
{
    CH_assert(a_numComp <= m_qNumComps);
    CH_assert(a_q.getBoxes().compatible(m_grids));
    CH_assert(a_q.getBoxes().physDomain().size() ==
              m_grids.physDomain().size());

    DataIterator                dit  = m_grids.dataIterator();
    const LevelData<FArrayBox>& oldQ = m_q[0];
    const LevelData<FArrayBox>& newQ = m_q[RKC::numStages - 1];

    if (a_numComp == -1) a_numComp = m_qNumComps - a_srcComp;

    // If we are at m_oldTime, just copy the old data.
    CH_assert(m_Q0Ready);
    if (abs(a_time - m_oldTime) < timeEps) {
        for (dit.reset(); dit.ok(); ++dit) {
            a_q[dit].copy(oldQ[dit], a_srcComp, a_destComp, a_numComp);
        }
        return 0.0;
    }

    // Linear interpolation...
    CH_assert(m_interpDataReady);

    // Compute theta and if nearby, snap into 0 or 1 values to prevent
    // extrapolation when not expected.
    const Real theta = this->computeTheta(a_time);

    if (abs(theta) < timeEps) {
        // A bit redundant...
        for (dit.reset(); dit.ok(); ++dit) {
            a_q[dit].copy(oldQ[dit], a_srcComp, a_destComp, a_numComp);
        }
        return 0.0;
    } else if (abs(theta - 1.0) < timeEps) {
        for (dit.reset(); dit.ok(); ++dit) {
            a_q[dit].copy(newQ[dit], a_srcComp, a_destComp, a_numComp);
        }
        return 1.0;
    }

    // Perform interpolation.
    for (dit.reset(); dit.ok(); ++dit) {
        a_q[dit].setVal(0.0);
        a_q[dit].plus(oldQ[dit], 1.0 - theta, a_srcComp, a_destComp, a_numComp);
        a_q[dit].plus(newQ[dit], theta, a_srcComp, a_destComp, a_numComp);
    }

    CH_assert(0.0 <= theta && theta <= 1.0);
    return theta;
}


// -----------------------------------------------------------------------------
template <class RKC>
Real
PARK<RKC>::velDenseOutput(LevelData<FluxBox>& a_vel,
                          const Real          a_time,
                          int                 a_srcComp,
                          int                 a_destComp,
                          int                 a_numComp) const
{
    CH_assert(RKC::hasDenseOutput);
    CH_assert(a_numComp <= m_velNumComps);
    CH_assert(a_vel.getBoxes().compatible(m_grids));
    CH_assert(a_vel.getBoxes().physDomain().size() ==
              m_grids.physDomain().size());

    DataIterator dit = m_grids.dataIterator();
    if (a_numComp == -1) a_numComp = m_velNumComps - a_srcComp;

    // No matter what, we need to start here.
    // And if we are at m_oldTime, we are done.
    CH_assert(m_Q0Ready);
    for (dit.reset(); dit.ok(); ++dit) {
        a_vel[dit].copy(m_vel[0][dit], a_srcComp, a_destComp, a_numComp);
    }
    if (abs(a_time - m_oldTime) < timeEps) return 0.0;

    // Interpolation...
    CH_assert(m_interpDataReady);

    // Compute theta and if nearby, snap into 0 or 1 values to prevent
    // extrapolation when not expected.
    const Real theta = this->computeTheta(a_time);

    // Compute powers of theta.
    Real thetaPow[RKC::numDenseCoeffs];
    for (int j = 0; j < RKC::numDenseCoeffs; ++j) {
        thetaPow[j] = pow(theta, j + 1);
    }

    // Compute dt * bstar_i(theta).
    Real dtbstar[RKC::numStages];
    for (int i = 0; i < RKC::numStages; ++i) {
        dtbstar[i] = 0.0;
        for (int j = 0; j < RKC::numDenseCoeffs; ++j) {
            dtbstar[i] += RKC::bstar[j][i] * thetaPow[j];
        }
        dtbstar[i] *= m_dt;
    }

    // Perform interpolation.
    for (int i = 0; i < RKC::numStages; ++i) {
        for (dit.reset(); dit.ok(); ++dit) {
            for (int dir = 0; dir < SpaceDim; ++dir) {
                a_vel[dit][dir].plus(m_kvelE[i][dit][dir],
                                     dtbstar[i],
                                     a_srcComp,
                                     a_destComp,
                                     a_numComp);
                a_vel[dit][dir].plus(m_kvelI[i][dit][dir],
                                     dtbstar[i],
                                     a_srcComp,
                                     a_destComp,
                                     a_numComp);
            }
        }
    }

    CH_assert(0.0 <= theta && theta <= 1.0);  // For now.
    return theta;
}


// -----------------------------------------------------------------------------
template <class RKC>
Real
PARK<RKC>::qDenseOutput(LevelData<FArrayBox>& a_q,
                        const Real            a_time,
                        int                   a_srcComp,
                        int                   a_destComp,
                        int                   a_numComp) const
{
    CH_assert(RKC::hasDenseOutput);
    CH_assert(a_numComp <= m_qNumComps);
    CH_assert(a_q.getBoxes().compatible(m_grids));
    CH_assert(a_q.getBoxes().physDomain().size() ==
              m_grids.physDomain().size());

    DataIterator dit = m_grids.dataIterator();
    if (a_numComp == -1) a_numComp = m_qNumComps - a_srcComp;

    // No matter what, we need to start here.
    // And if we are at m_oldTime, we are done.
    CH_assert(m_Q0Ready);
    for (dit.reset(); dit.ok(); ++dit) {
        a_q[dit].copy(m_q[0][dit], a_srcComp, a_destComp, a_numComp);
    }
    if (abs(a_time - m_oldTime) < timeEps) return 0.0;

    // Interpolation...
    CH_assert(m_interpDataReady);

    // Compute theta and if nearby, snap into 0 or 1 values to prevent
    // extrapolation when not expected.
    const Real theta = this->computeTheta(a_time);

    // Compute powers of theta.
    Real thetaPow[RKC::numDenseCoeffs];
    for (int j = 0; j < RKC::numDenseCoeffs; ++j) {
        thetaPow[j] = pow(theta, j + 1);
    }

    // Compute dt * bstar_i(theta).
    Real dtbstar[RKC::numStages];
    for (int i = 0; i < RKC::numStages; ++i) {
        dtbstar[i] = 0.0;
        for (int j = 0; j < RKC::numDenseCoeffs; ++j) {
            dtbstar[i] += RKC::bstar[j][i] * thetaPow[j];
        }
        dtbstar[i] *= m_dt;
    }

    // Perform interpolation.
    for (int i = 0; i < RKC::numStages; ++i) {
        for (dit.reset(); dit.ok(); ++dit) {
            a_q[dit].plus(
                m_kqE[i][dit], dtbstar[i], a_srcComp, a_destComp, a_numComp);
            a_q[dit].plus(
                m_kqI[i][dit], dtbstar[i], a_srcComp, a_destComp, a_numComp);
        }
    }

    CH_assert(0.0 <= theta && theta <= 1.0);  // For now.
    return theta;
}


// -----------------------------------------------------------------------------
template <class RKC>
Real
PARK<RKC>::velLagrangeInterp(LevelData<FluxBox>& a_vel,
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

    // If we are at m_oldTime, just copy the old data.
    CH_assert(m_Q0Ready);
    if (abs(a_time - m_oldTime) < timeEps) {
        for (dit.reset(); dit.ok(); ++dit) {
            a_vel[dit].copy(m_vel[0][dit], a_srcComp, a_destComp, a_numComp);
        }
        return 0.0;
    }

    CH_assert(m_interpDataReady);

    // Compute theta and if nearby, snap into 0 or 1 values to prevent
    // extrapolation when not expected.
    const Real theta = this->computeTheta(a_time);

    if (abs(theta) < timeEps) {
        // A bit redundant...
        for (dit.reset(); dit.ok(); ++dit) {
            a_vel[dit].copy(m_vel[0][dit], a_srcComp, a_destComp, a_numComp);
        }
        return 0.0;
    } else if (abs(theta - 1.0) < timeEps) {
        for (dit.reset(); dit.ok(); ++dit) {
            a_vel[dit].copy(m_vel[RKC::numStages - 1][dit], a_srcComp, a_destComp, a_numComp);
        }
        return 1.0;
    }

    // Perform interpolation.
    std::array<Real, RKC::numStages> lag;
    for (size_t j = 0; j < RKC::numStages; ++j) {
        const Real tj = m_oldTime + RKC::c[j] * m_dt;

        lag[j] = 1.0;
        for (size_t r = 0; r < RKC::numStages; ++r) {
            if (r == j) continue;

            const Real tr = m_oldTime + RKC::c[r] * m_dt;
            lag[j] *= (a_time - tr) / (tj - tr);
        }
    }
    for (dit.reset(); dit.ok(); ++dit) {
        a_vel[dit].setVal(0.0);
        for (size_t j = 0; j < RKC::numStages; ++j) {
            for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                a_vel[dit][velComp].plus(m_vel[j][dit][velComp], lag[j], a_srcComp, a_destComp, a_numComp);
            }
        }
    }

    CH_assert(0.0 <= theta && theta <= 1.0);
    return theta;
}


// -----------------------------------------------------------------------------
template <class RKC>
Real
PARK<RKC>::qLagrangeInterp(LevelData<FArrayBox>& a_q,
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

    // If we are at m_oldTime, just copy the old data.
    CH_assert(m_Q0Ready);
    if (abs(a_time - m_oldTime) < timeEps) {
        for (dit.reset(); dit.ok(); ++dit) {
            a_q[dit].copy(m_q[0][dit], a_srcComp, a_destComp, a_numComp);
        }
        return 0.0;
    }

    CH_assert(m_interpDataReady);

    // Compute theta and if nearby, snap into 0 or 1 values to prevent
    // extrapolation when not expected.
    const Real theta = this->computeTheta(a_time);

    if (abs(theta) < timeEps) {
        // A bit redundant...
        for (dit.reset(); dit.ok(); ++dit) {
            a_q[dit].copy(m_q[0][dit], a_srcComp, a_destComp, a_numComp);
        }
        return 0.0;
    } else if (abs(theta - 1.0) < timeEps) {
        for (dit.reset(); dit.ok(); ++dit) {
            a_q[dit].copy(m_q[RKC::numStages - 1][dit], a_srcComp, a_destComp, a_numComp);
        }
        return 1.0;
    }

    // Perform interpolation.
    std::array<Real, RKC::numStages> lag;
    for (size_t j = 0; j < RKC::numStages; ++j) {
        const Real tj = m_oldTime + RKC::c[j] * m_dt;

        lag[j] = 1.0;
        for (size_t r = 0; r < RKC::numStages; ++r) {
            if (r == j) continue;

            const Real tr = m_oldTime + RKC::c[r] * m_dt;
            lag[j] *= (a_time - tr) / (tj - tr);
        }
    }
    for (dit.reset(); dit.ok(); ++dit) {
        a_q[dit].setVal(0.0);
        for (size_t j = 0; j < RKC::numStages; ++j) {
            a_q[dit].plus(m_q[j][dit], lag[j], a_srcComp, a_destComp, a_numComp);
        }
    }

    CH_assert(0.0 <= theta && theta <= 1.0);
    return theta;
}


// -----------------------------------------------------------------------------
template <class RKC>
Real
PARK<RKC>::pLagrangeInterp(LevelData<FArrayBox>& a_p,
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

    // If we are at m_oldTime, just copy the old data.
    CH_assert(m_Q0Ready);
    if (abs(a_time - m_oldTime) < timeEps) {
        for (dit.reset(); dit.ok(); ++dit) {
            a_p[dit].copy(m_p[0][dit], a_srcComp, a_destComp, a_numComp);
        }
        return 0.0;
    }

    CH_assert(m_interpDataReady);

    // Compute theta and if nearby, snap into 0 or 1 values to prevent
    // extrapolation when not expected.
    const Real theta = this->computeTheta(a_time);

    if (abs(theta) < timeEps) {
        // A bit redundant...
        for (dit.reset(); dit.ok(); ++dit) {
            a_p[dit].copy(m_p[0][dit], a_srcComp, a_destComp, a_numComp);
        }
        return 0.0;
    } else if (abs(theta - 1.0) < timeEps) {
        for (dit.reset(); dit.ok(); ++dit) {
            a_p[dit].copy(m_p[RKC::numStages - 1][dit], a_srcComp, a_destComp, a_numComp);
        }
        return 1.0;
    }

    // Perform interpolation.
    std::array<Real, RKC::numStages> lag;
    for (size_t j = 0; j < RKC::numStages; ++j) {
        const Real tj = m_oldTime + RKC::c[j] * m_dt;

        lag[j] = 1.0;
        for (size_t r = 0; r < RKC::numStages; ++r) {
            if (r == j) continue;

            const Real tr = m_oldTime + RKC::c[r] * m_dt;
            lag[j] *= (a_time - tr) / (tj - tr);
        }
    }
    for (dit.reset(); dit.ok(); ++dit) {
        a_p[dit].setVal(0.0);
        for (size_t j = 0; j < RKC::numStages; ++j) {
            a_p[dit].plus(m_p[j][dit], lag[j], a_srcComp, a_destComp, a_numComp);
        }
    }

    CH_assert(0.0 <= theta && theta <= 1.0);
    return theta;
}


// -----------------------------------------------------------------------------
template <class RKC>
void
PARK<RKC>::setOldQ(const LevelData<FluxBox>&   a_velOld,
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

#ifndef NDEBUG
    for (int i = 0; i < RKC::numStages; ++i) {
        debugInitLevel(m_vel[i]);
        debugInitLevel(m_kvelE[i]);
        debugInitLevel(m_kvelI[i]);

        debugInitLevel(m_p[i]);

        debugInitLevel(m_q[i]);
        debugInitLevel(m_kqE[i]);
        debugInitLevel(m_kqI[i]);
    }

    m_oldTime = quietNAN;
    m_newTime = quietNAN;
    m_dt      = quietNAN;
#endif

    DataIterator dit = m_grids.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        for (int dir = 0; dir < SpaceDim; ++dir) {
            m_vel[0][dit][dir].copy(a_velOld[dit][dir]);
        }
        m_p[0][dit].copy(a_pOld[dit]);
        m_q[0][dit].copy(a_qOld[dit]);
    }
    m_oldTime = a_oldTime;

    m_Q0Ready         = true;
    m_interpDataReady = false;
}


// -----------------------------------------------------------------------------
template <class RKC>
void
PARK<RKC>::copy(LevelData<FluxBox>& a_x, const LevelData<FluxBox>& a_y)
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
template <class RKC>
void
PARK<RKC>::copy(LevelData<FArrayBox>& a_x, const LevelData<FArrayBox>& a_y)
{
    CH_assert(a_x.getBoxes() == a_y.getBoxes());
    CH_assert(a_x.nComp() == a_y.nComp());

    DataIterator dit = a_x.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        a_x[dit].copy(a_y[dit]);
    }
}


// -----------------------------------------------------------------------------
template <class RKC>
void
PARK<RKC>::plus(LevelData<FluxBox>&       a_x,
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
template <class RKC>
void
PARK<RKC>::plus(LevelData<FArrayBox>&       a_x,
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
template <class RKC>
void
PARK<RKC>::axby(const Real                a_a,
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
template <class RKC>
void
PARK<RKC>::axby(const Real                  a_a,
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


// -----------------------------------------------------------------------------
template <class RKC>
template <typename E, size_t N>
void
PARK<RKC>::pushToHistory(const E&          a_newElem,
                         std::array<E, N>& a_history,
                         const bool        a_isNewTime)
{
    if (a_isNewTime) {
        for (size_t i = N - 1; i > 0; --i) {
            a_history[i] = a_history[i - 1];
        }
    }
    a_history[0] = a_newElem;
}


// -----------------------------------------------------------------------------
template <class RKC>
Real
PARK<RKC>::softLimiter(Real       a_proposedController,
                       const Real a_lastController,
                       const Real a_maxGrowth)
{
    if (a_lastController > 0.0) {
        const Real arg = (a_proposedController - 1.0) / a_maxGrowth;
        a_proposedController *= 1.0 + a_maxGrowth * std::atan(arg);
    }
    return a_proposedController;
}


// -----------------------------------------------------------------------------
template <class RKC>
Real
PARK<RKC>::hardLimiter(Real       a_proposedController,
                       const Real a_lastController,
                       const Real a_minShrinkage,
                       const Real a_maxGrowth)
{
    if (a_lastController > 0.0) {
        a_proposedController =
            std::max(a_proposedController, a_minShrinkage * a_lastController);
        a_proposedController =
            std::min(a_proposedController, a_maxGrowth * a_lastController);
    }
    return a_proposedController;
}


// -----------------------------------------------------------------------------
template <class RKC>
Real
PARK<RKC>::oscillationLimiter(Real       a_proposedController,
                              const Real a_oldController,
                              const Real a_reallyOldController,
                              const Real a_dampingFactor)
{
    if (a_dampingFactor <= 0.0) {
        return a_proposedController;
    }

    if (a_oldController > 0.0 && a_reallyOldController > 0.0) {
        const Real oldSlope = a_oldController - a_reallyOldController;
        const Real newSlope = a_proposedController - a_oldController;

        if (oldSlope * newSlope < 0.0) {
            const Real avg = 0.5 * (a_oldController + a_reallyOldController);

            a_proposedController =
                    (1.0 - a_dampingFactor) * a_proposedController
                  + (      a_dampingFactor) * avg;
        }
    }

    return a_proposedController;
}


// -----------------------------------------------------------------------------
template <int numStages>
constexpr int
RKOrderChecks(const Real a[numStages][numStages],
              const Real b[numStages],
              const Real c[numStages],
              const bool a_quietMode)
{
    const Real small = timeEps;
    int        order = 0;

    do {
        // Consistency checks
        for (int i = 0; i < numStages; ++i) {
            Real sum = 0.0;
            for (int j = 0; j < numStages; ++j) {
                sum += a[i][j];
            }

            Real error = abs(sum - c[i]);
            if (error > small) {
                if (!a_quietMode) {
                    pout() << "RK method is inconsistent:"
                           << "\nSum[a_" << i << "][j] = " << sum << " != c["
                           << i << "] = " << c[i] << endl;
                }
                break;
            }
        }
        order = 0;

        // 1st order checks
        {
            Real sum = 0.0;
            for (int i = 0; i < numStages; ++i) {
                sum += b[i];
            }

            Real error = abs(sum - 1.0);
            if (error > small) {
                if (!a_quietMode) {
                    pout() << "RK order check failed:"
                           << "\nSum[b_i] = " << sum << " != 1.0" << endl;
                }
                break;
            }
        }
        order = 1;

        // 2nd order checks
        {
            Real sum = 0.0;
            for (int i = 0; i < numStages; ++i) {
                sum += b[i] * c[i];
            }

            Real error = abs(sum - 0.5);
            if (error > small) {
                if (!a_quietMode) {
                    pout() << "RK order check failed:"
                           << "\nSum[b_i * c_i] = " << sum << " != 0.5" << endl;
                }
                break;
            }
        }
        order = 2;

        // 3rd order checks
        {
            Real sum = 0.0;
            for (int i = 0; i < numStages; ++i) {
                sum += b[i] * c[i] * c[i];
            }

            Real error = abs(sum - 1.0 / 3.0);
            if (error > small) {
                if (!a_quietMode) {
                    pout() << "RK order check failed:"
                           << "\nSum[b_i * c_i^2] = " << sum << " != 1.0/3.0"
                           << endl;
                }
                break;
            }
        }
        {
            Real sum = 0.0;
            for (int i = 0; i < numStages; ++i) {
                for (int j = 0; j < numStages; ++j) {
                    sum += b[i] * a[i][j] * c[j];
                }
            }

            Real error = abs(sum - 1.0 / 6.0);
            if (error > small) {
                if (!a_quietMode) {
                    pout() << "RK order check failed:"
                           << "\nSum[b_i * a_ij * c_j] = " << sum
                           << " != 1.0/6.0" << endl;
                }
                break;
            }
        }
        order = 3;

        // 4th order checks
        {
            Real sum = 0.0;
            for (int i = 0; i < numStages; ++i) {
                sum += b[i] * c[i] * c[i] * c[i];
            }

            Real error = abs(sum - 0.25);
            if (error > small) {
                if (!a_quietMode) {
                    pout() << "RK order check failed:"
                           << "\nSum[b_i * c_i^3] = " << sum << " != 0.25"
                           << endl;
                }
                break;
            }
        }
        {
            Real sum = 0.0;
            for (int i = 0; i < numStages; ++i) {
                for (int j = 0; j < numStages; ++j) {
                    sum += b[i] * c[i] * a[i][j] * c[j];
                }
            }

            Real error = abs(sum - 1.0 / 8.0);
            if (error > small) {
                if (!a_quietMode) {
                    pout() << "RK order check failed:"
                           << "\nSum[b_i * c_i * a_ij * c_j] = " << sum
                           << " != 1.0/8.0" << endl;
                }
                break;
            }
        }
        {
            Real sum = 0.0;
            for (int i = 0; i < numStages; ++i) {
                for (int j = 0; j < numStages; ++j) {
                    sum += b[i] * a[i][j] * c[j] * c[j];
                }
            }

            Real error = abs(sum - 1.0 / 12.0);
            if (error > small) {
                if (!a_quietMode) {
                    pout() << "RK order check failed:"
                           << "\nSum[b_i * a_ij * c_j^2] = " << sum
                           << " != 1.0/12.0" << endl;
                }
                break;
            }
        }
        {
            Real sum = 0.0;
            for (int i = 0; i < numStages; ++i) {
                for (int j = 0; j < numStages; ++j) {
                    for (int k = 0; k < numStages; ++k) {
                        sum += b[i] * a[i][j] * a[j][k] * c[k];
                    }
                }
            }

            Real error = abs(sum - 1.0 / 24.0);
            if (error > small) {
                if (!a_quietMode) {
                    pout() << "RK order check failed:"
                           << "\nSum[b_i * a_ij * a_jk * c_k] = " << sum
                           << " != 1.0/24.0" << endl;
                }
                break;
            }
        }
        order = 4;
    } while (0);

    if (!a_quietMode) {
        if (order < 4) {
            pout() << "RK order is " << order << endl;
        } else {
            pout() << "RK order is at least " << order << endl;
        }
    }

    return order;
}


// -----------------------------------------------------------------------------
template <int numStages>
constexpr int
RKStageOrderChecks(const Real a[numStages][numStages],
                   const Real bE[numStages],
                   const Real bI[numStages],
                   const Real c[numStages],
                   const bool a_quietMode)
{
    const Real small = timeEps;
    Real       error = 0.0;
    int        order = 0;

    do {
        // 1st order checks
        for (int i = 0; i < numStages; ++i) {
            Real sum = 0.0;
            for (int j = 0; j < numStages; ++j) {
                sum += a[i][j];
            }
            error = abs(sum - c[i]);
            if (error > small) {
                if (!a_quietMode) {
                    pout() << "RK stage order 1 test failed:"
                           << "\nSum[a_" << i << "][j] = " << sum << " != c["
                           << i << "] = " << c[i] << endl;
                }
                break;
            }
        }
        if (error > small) break;
        order = 1;

        // 2nd order checks
        for (int i = 0; i < numStages; ++i) {
            error = std::abs(bE[i] - bI[i]);
            if (error > small) {
                if (!a_quietMode) {
                    pout() << "RK stage order 2 check failed:"
                           << "\nbE[" << i << "] = " << bE[i] << " != bI[" << i
                           << "] = " << bI[i] << endl;
                }
                break;
            }
        }
        if (error > small) break;

        for (int i = 0; i < numStages; ++i) {
            Real sum = 0.0;
            for (int j = 0; j < numStages; ++j) {
                sum += a[i][j] * c[j];
            }
            error = std::abs(sum - 0.5 * c[i] * c[i]);
            if (error > small) {
                if (!a_quietMode) {
                    pout() << "RK stage order 2 check failed:"
                           << "\nSum[a_{" << i << " j} * c_j] = " << sum
                           << " != 0.5 * c[" << i
                           << "]^2 = " << 0.5 * c[i] * c[i] << endl;
                }
                break;
            }
        }
        if (error > small) break;
        order = 2;
    } while (0);

    if (!a_quietMode) {
        if (order < 2) {
            pout() << "RK stage order is " << order << endl;
        } else {
            pout() << "RK stage order is at least " << order << endl;
        }
    }
    return order;
}


// Cleanup
#undef nanCheck
