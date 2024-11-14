#include "SDC.H"
#include "SOMAR_Constants.H"
#include "Debug.H"
#include "SetValLevel.H"
#include "Analysis.H"


#ifndef NDEBUG
    // Debug mode
#   define nanCheck(x) checkForValidNAN(x)
#else
    // Release mode
#   define nanCheck(x)
#endif


template<>
const std::array<Real, 3> LobattoQuadrature<3>::s_nodes =
    { -1.0
    ,  0.0
    ,  1.0};
template<>
const std::array<Real, 3> LobattoQuadrature<3>::s_weights =
    {  1.0 / 3.0
    ,  4.0 / 3.0
    ,  1.0 / 3.0};

template<>
const std::array<Real, 4> LobattoQuadrature<4>::s_nodes =
    { -1.0
    , -std::sqrt(5.0) / 5.0
    ,  std::sqrt(5.0) / 5.0
    ,  1.0};
template<>
const std::array<Real, 4> LobattoQuadrature<4>::s_weights =
    {  1.0 / 6.0
    ,  5.0 / 6.0
    ,  5.0 / 6.0
    ,  1.0 / 6.0};

template<>
const std::array<Real, 5> LobattoQuadrature<5>::s_nodes =
    { -1.0
    , -std::sqrt(21.0) / 7.0
    ,  0.0
    ,  std::sqrt(21.0) / 7.0
    ,  1.0};
template<>
const std::array<Real, 5> LobattoQuadrature<5>::s_weights =
    {  1.0 / 10.0
    ,  49.0 / 90.0
    ,  32.0 / 45.0
    ,  49.0 / 90.0
    ,  1.0 / 10.0};

template<>
const std::array<Real, 6> LobattoQuadrature<6>::s_nodes =
    { -1.0
    , -std::sqrt((7.0 + 2.0*std::sqrt(7.0)) / 21.0)
    , -std::sqrt((7.0 - 2.0*std::sqrt(7.0)) / 21.0)
    ,  std::sqrt((7.0 - 2.0*std::sqrt(7.0)) / 21.0)
    ,  std::sqrt((7.0 + 2.0*std::sqrt(7.0)) / 21.0)
    ,  1.0};
template<>
const std::array<Real, 6> LobattoQuadrature<6>::s_weights =
    {  1.0 / 15.0
    ,  (14.0 - std::sqrt(7.0)) / 30.0
    ,  (14.0 + std::sqrt(7.0)) / 30.0
    ,  (14.0 + std::sqrt(7.0)) / 30.0
    ,  (14.0 - std::sqrt(7.0)) / 30.0
    ,  1.0 / 15.0};

template<>
const std::array<Real, 7> LobattoQuadrature<7>::s_nodes =
    { -1.0
    , -std::sqrt((5.0 + 2.0 * std::sqrt(5.0 / 3.0)) / 11.0)
    , -std::sqrt((5.0 - 2.0 * std::sqrt(5.0 / 3.0)) / 11.0)
    ,  0.0
    ,  std::sqrt((5.0 - 2.0 * std::sqrt(5.0 / 3.0)) / 11.0)
    ,  std::sqrt((5.0 + 2.0 * std::sqrt(5.0 / 3.0)) / 11.0)
    ,  1.0};
template<>
const std::array<Real, 7> LobattoQuadrature<7>::s_weights =
    {  1.0 / 21.0
    ,  (124.0 - 7.0 * std::sqrt(15.0)) / 350.0
    ,  (124.0 + 7.0 * std::sqrt(15.0)) / 350.0
    ,  256.0 / 525.0
    ,  (124.0 + 7.0 * std::sqrt(15.0)) / 350.0
    ,  (124.0 - 7.0 * std::sqrt(15.0)) / 350.0
    ,  1.0 / 21.0};



// -----------------------------------------------------------------------------
SDC::SDC(const DisjointBoxLayout& a_grids,
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
// State variables
, m_oldTime(quietNAN)
, m_newTime(quietNAN)
, m_dt(quietNAN)
, m_Q0Ready(false)
, m_finalForcesReady(false)
, m_interpDataReady(false)
{
    for (size_t m = 0; m < NumNodes; ++m) {
        m_vel[m].define(a_grids, a_velNumComps, a_velGhostVect);
        m_p[m].define(a_grids, a_pNumComps, a_pGhostVect);
        m_q[m].define(a_grids, a_qNumComps, a_qGhostVect);
        debugInitLevel(m_vel[m]);
        debugInitLevel(m_p[m]);
        debugInitLevel(m_q[m]);

        m_kvelE[m].define(a_grids, a_velNumComps, a_velGhostVect);
        m_kqE[m].define(a_grids, a_qNumComps, a_qGhostVect);
        debugInitLevel(m_kvelE[m]);
        debugInitLevel(m_kqE[m]);

        m_kvelI[m].define(a_grids, a_velNumComps, a_velGhostVect);
        m_kqI[m].define(a_grids, a_qNumComps, a_qGhostVect);
        debugInitLevel(m_kvelI[m]);
        debugInitLevel(m_kqI[m]);

        m_kvelOld[m].define(a_grids, a_velNumComps, a_velGhostVect);
        m_kqOld[m].define(a_grids, a_qNumComps, a_qGhostVect);
        debugInitLevel(m_kvelOld[m]);
        debugInitLevel(m_kqOld[m]);
    }

    for (size_t m = 1; m < NumNodes; ++m) {
        const Real xim  = Quadrature::getNode(m, -1.0, 1.0);
        const Real xim1 = Quadrature::getNode(m-1, -1.0, 1.0);

        for (size_t l = 0; l < NumNodes; ++l) {
            const Real xil = Quadrature::getNode(l, -1.0, 1.0);

            m_correctionWeights[m-1][l] = 0.0;
            for (size_t a = 0; a < NumNodes; ++a) {
                const Real xia = Quadrature::getNode(a, -1.0, 1.0);
                const Real wa  = Quadrature::getWeight(a, -1.0, 1.0);

                Real prod = 1.0;
                for (size_t j = 0; j < NumNodes; ++j) {
                    if (j == l) continue;
                    const Real xij = Quadrature::getNode(j, -1.0, 1.0);
                    const Real xiaStar = 0.5 * ((xim - xim1) * xia + xim + xim1);
                    prod *= (xiaStar - xij) / (xil - xij);
                }

                m_correctionWeights[m-1][l] += wa * prod;
            }
            m_correctionWeights[m-1][l] *= 0.25 * (xim - xim1);
        } // l
    } // m


    POUT(m_correctionWeights);
}


// -----------------------------------------------------------------------------
void
SDC::advance(LevelData<FluxBox>&   a_vel,
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

    nanCheck(a_vel);
    nanCheck(a_p  );
    nanCheck(a_q  );

    const std::array<Real, NumNodes> stageTime = [&a_oldTime, &a_dt]() {
        std::array<Real, NumNodes> ret;
        for (size_t m = 0; m < NumNodes; ++m) {
            ret[m] = Quadrature::getNode(m, a_oldTime, a_oldTime + a_dt);
        }
        return ret;
    }();

    const std::array<Real, NumNodes - 1> stageDt = [&stageTime]() {
        std::array<Real, NumNodes - 1> ret;
        for (size_t m = 0; m < NumNodes - 1; ++m) {
            ret[m] = stageTime[m + 1] - stageTime[m];
        }
        return ret;
    }();

    // if (m_finalForcesReady) {
    //     SDC::copy(m_kvelE.back(), m_kvelE.front());
    //     SDC::copy(m_kvelI.back(), m_kvelI.front());
    // } else {
    //     // TODO
    // }

    // Initialization. This sets Q[0].
    this->setOldQ(a_vel, a_p, a_q, a_oldTime);
    m_dt = a_dt;

    // ------------------------------------------------------------
    // Predictor. Computes Q[m] and kQ[0][m] for all stages m >= 1.
    size_t k = 0;

    // Compute initial implicit forces.
    {
        constexpr Real refluxDt = 0.0;
        constexpr size_t m = 0;
        a_rhsPtr->setImplicitRHS(m_kvelI[m], m_kqI[m], m_vel[m], m_p[m], m_q[m], stageTime[m], refluxDt);
    }
    // Compute stages >= 1.
    for (size_t m = 1; m < NumNodes; ++m) {
        pout() << Format::fixed << "FE timestep (predictor k = " << k
               << ", subinterval m = " << m << ")\n"
               << Format::indent() << flush;

        constexpr Real refluxDt = 0.0;
        SDC::basicAdvance(m_vel[m], m_p[m], m_q[m],
                          m_kvelE[m-1], m_kqE[m-1],
                          m_kvelI[m], m_kqI[m],
                          m_vel[m-1], m_p[m-1], m_q[m-1],
                          stageTime[m-1], stageDt[m-1], refluxDt,
                          a_rhsPtr);

        pout() << Format::unindent << flush;
    }
    // Compute final explicit forces.
    {
        const Real refluxDt = 0.0;
        constexpr size_t m = NumNodes - 1;
        a_rhsPtr->setExplicitRHS(m_kvelE[m], m_kqE[m], m_vel[m], m_p[m], m_q[m], stageTime[m], refluxDt);
    }

    LevelData<FluxBox> velErr(m_grids, 1);
    SDC::copy(velErr, m_vel.back());

    // -------------------------------------------------------
    // Deferred corrections
    constexpr size_t numCorrections = NumNodes - 1;
    for (k = 1; k <= numCorrections; ++k) {
        // Save k-1 forces.
        for (size_t m = 0; m < NumNodes; ++m) {
            SDC::copy(m_kvelOld[m], m_kvelE[m]);
            SDC::plus(m_kvelOld[m], 1.0, m_kvelI[m]);

            SDC::copy(m_kqOld[m], m_kqE[m]);
            SDC::plus(m_kqOld[m], 1.0, m_kqI[m]);
        }

        // Correct stages >= 1.
        for (size_t m = 1; m < NumNodes; ++m) {
            pout() << Format::fixed << "FE timestep (corrector k = " << k
                << ", subinterval m = " << m << ")\n"
                << Format::indent() << flush;

            SDC::copy(m_vel[m], m_vel[m-1]);
            SDC::plus(m_vel[m], -stageDt[m-1], m_kvelE[m-1]);
            SDC::plus(m_vel[m], -stageDt[m-1], m_kvelI[m]);

            SDC::copy(m_p[m], m_p[m-1]);

            SDC::copy(m_q[m], m_q[m-1]);
            SDC::plus(m_q[m], -stageDt[m-1], m_kqE[m-1]);
            SDC::plus(m_q[m], -stageDt[m-1], m_kqI[m]);

            for (size_t l = 0; l < NumNodes; ++l) {
                SDC::plus(m_vel[m], a_dt * m_correctionWeights[m-1][l], m_kvelOld[l]);
                SDC::plus(  m_q[m], a_dt * m_correctionWeights[m-1][l],   m_kqOld[l]);
            }

            Real sum = 0;
            for (size_t l = 0; l < NumNodes; ++l) {
                sum += m_correctionWeights[m-1][l];
            }
            CH_verify(RealCmp::eq(sum*a_dt, stageDt[m-1]));

            const Real refluxDt = 0.0;
            SDC::basicAdvance(m_vel[m], m_p[m], m_q[m],
                              m_kvelE[m-1], m_kqE[m-1],
                              m_kvelI[m], m_kqI[m],
                              m_vel[m], m_p[m], m_q[m],
                              stageTime[m-1], stageDt[m-1], refluxDt,
                              a_rhsPtr);

            pout() << Format::unindent << flush;
        } // m
        // Compute final explicit forces.
        {
            const Real refluxDt = 0.0;
            constexpr size_t m = NumNodes - 1;
            a_rhsPtr->setExplicitRHS(m_kvelE[m], m_kqE[m], m_vel[m], m_p[m], m_q[m], stageTime[m], refluxDt);
        }

        // Error analysis
        SDC::plus(velErr, -1.0, m_vel.back());
        const RealVect velErrNorm = Analysis::pNorm(velErr, 2);
        pout() << "|velErr|_2 ~ " << Format::pushFlags << Format::scientific
               << velErrNorm << Format::popFlags << '\n';
        SDC::copy(velErr, m_vel.back());
    } // k

    // Update user's state.
    this->copy(a_vel, m_vel.back());
    this->copy(  a_p,   m_p.back());
    this->copy(  a_q,   m_q.back());

    m_newTime          = m_oldTime + m_dt;
    m_Q0Ready          = true;
    m_finalForcesReady = true;
    m_interpDataReady  = true;

    // Allow user to perform postStep operations (e.g., preparing plots)
    a_rhsPtr->postStep(a_vel, a_p, a_q, m_newTime);
}


// -----------------------------------------------------------------------------
void
SDC::basicAdvance(LevelData<FluxBox>&   a_vel,
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
                  PARKRHS*              a_rhsPtr)
{
    // Explicit Euler. Initializes new stage.
    a_rhsPtr->setExplicitRHS(a_kvelE, a_kqE, a_vel0, a_p0, a_q0, a_time0, a_refluxDt);
    SDC::copy(a_vel, a_vel0);
    SDC::copy(a_q  , a_q0);
    SDC::plus(a_vel, a_dt, a_kvelE);
    SDC::plus(a_q  , a_dt, a_kqE);

    // Approximate projector. Initializes new pressure.
    SDC::copy(a_p, a_p0);
    a_rhsPtr->projectPredict(a_vel, a_p, a_time0, a_dt);

    // Implicit Euler.
    SDC::copy(a_kvelI, a_vel);
    SDC::copy(a_kqI, a_q);
    a_rhsPtr->solveImplicit(a_vel, a_q, a_dt, a_time0, a_refluxDt);
    SDC::axby(-1.0 / a_dt, a_kvelI, 1.0 / a_dt, a_vel);
    SDC::axby(-1.0 / a_dt, a_kqI, 1.0 / a_dt, a_q);

    // Projection correction. Completes pressure update.
    a_rhsPtr->projectCorrect(a_vel, a_p, a_time0 + a_dt, a_dt);

    nanCheck(a_vel);
    nanCheck(a_p);
    nanCheck(a_q);

    nanCheck(a_kvelE);
    nanCheck(a_kqE);

    nanCheck(a_kvelI);
    nanCheck(a_kqI);
}


// -----------------------------------------------------------------------------
void
SDC::FEadvance(LevelData<FluxBox>&   a_vel,
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

    nanCheck(a_vel);
    nanCheck(a_p  );
    nanCheck(a_q  );

    this->setOldQ(a_vel, a_p, a_q, a_oldTime);
    m_dt = a_dt;


    // --- Begin stage ---
    pout() << Format::fixed << "FE timestep..."
           << "\n"
           << Format::indent() << flush;

    // Evaluate forces.
    constexpr size_t m = 0;
    a_rhsPtr->setExplicitRHS(m_kvelE[m], m_kqE[m], a_vel, a_p, a_q, a_oldTime, a_dt);
    a_rhsPtr->setImplicitRHS(m_kvelI[m+1], m_kqI[m+1], a_vel, a_p, a_q, a_oldTime, a_dt);

    nanCheck(m_kvelE[m]);
    nanCheck(m_kqE[m]);
    nanCheck(m_kvelI[m+1]);
    nanCheck(m_kqI[m+1]);

    // Update q (Ri)
    // Q += dt*(kQE + kQI)
    this->plus(a_vel, a_dt, m_kvelE[m]);
    this->plus(a_vel, a_dt, m_kvelI[m+1]);
    this->plus(a_q, a_dt, m_kqE[m]);
    this->plus(a_q, a_dt, m_kqI[m+1]);

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
    m_finalForcesReady = false;
    m_interpDataReady  = false;
    // ...also, we will not call postStep. I doubt you want to prepare
    // plots using this final state.
}


// -----------------------------------------------------------------------------
Real
SDC::velTimeInterp(LevelData<FluxBox>& a_vel,
                   const Real          a_time,
                   int                 a_srcComp,
                   int                 a_destComp,
                   int                 a_numComp) const
{
    return this->velLinearInterp(a_vel, a_time, a_srcComp, a_destComp, a_numComp);
}


// -----------------------------------------------------------------------------
Real
SDC::pTimeInterp(LevelData<FArrayBox>& a_p,
                 const Real            a_time,
                 int                   a_srcComp,
                 int                   a_destComp,
                 int                   a_numComp) const
{
    return this->pLinearInterp(a_p, a_time, a_srcComp, a_destComp, a_numComp);
}


// -----------------------------------------------------------------------------
Real
SDC::qTimeInterp(LevelData<FArrayBox>& a_q,
                 const Real            a_time,
                 int                   a_srcComp,
                 int                   a_destComp,
                 int                   a_numComp) const
{
    return this->qLinearInterp(a_q, a_time, a_srcComp, a_destComp, a_numComp);
}


// -----------------------------------------------------------------------------
Real
SDC::computeTheta(const Real a_time) const
{
    const Real theta = (a_time - m_oldTime) / m_dt;
    if (abs(theta - 0.0) < timeEps) return 0.0;
    if (abs(theta - 1.0) < timeEps) return 1.0;
    return theta;
}


// -----------------------------------------------------------------------------
Real
SDC::velLinearInterp(LevelData<FluxBox>& a_vel,
                     const Real          a_time,
                     int                 a_srcComp,
                     int                 a_destComp,
                     int                 a_numComp) const
{
    CH_assert(a_numComp <= m_velNumComps);
    CH_assert(a_vel.getBoxes().compatible(m_grids));
    CH_assert(a_vel.getBoxes().physDomain().size() ==
              m_grids.physDomain().size());

    // DataIterator dit = m_grids.dataIterator();
    // const LevelData<FluxBox>& oldVel = m_vel[0];
    // const LevelData<FluxBox>& newVel = m_vel[RKC::numStages - 1];

    // if (a_numComp == -1) a_numComp = m_velNumComps - a_srcComp;

    // // If we are at m_oldTime, just copy the old data.
    // CH_assert(m_Q0Ready);
    // if (abs(a_time - m_oldTime) < timeEps) {
    //     for (dit.reset(); dit.ok(); ++dit) {
    //         for (int dir = 0; dir < SpaceDim; ++dir) {
    //             a_vel[dit][dir].copy(
    //                 oldVel[dit][dir], a_srcComp, a_destComp, a_numComp);
    //         }
    //     }
    //     return 0.0;
    // }

    // // Linear interpolation...
    // CH_assert(m_interpDataReady);

    // Compute theta and if nearby, snap into 0 or 1 values to prevent
    // extrapolation when not expected.
    const Real theta = this->computeTheta(a_time);

    // if (abs(theta) < timeEps) {
    //     // A bit redundant...
    //     for (dit.reset(); dit.ok(); ++dit) {
    //         for (int dir = 0; dir < SpaceDim; ++dir) {
    //             a_vel[dit][dir].copy(
    //                 oldVel[dit][dir], a_srcComp, a_destComp, a_numComp);
    //         }
    //     }
    //     return 0.0;
    // } else if (abs(theta - 1.0) < timeEps) {
    //     for (dit.reset(); dit.ok(); ++dit) {
    //         for (int dir = 0; dir < SpaceDim; ++dir) {
    //             a_vel[dit][dir].copy(
    //                 newVel[dit][dir], a_srcComp, a_destComp, a_numComp);
    //         }
    //     }
    //     return 1.0;
    // }

    // // Perform interpolation.
    // for (dit.reset(); dit.ok(); ++dit) {
    //     for (int dir = 0; dir < SpaceDim; ++dir) {
    //         a_vel[dit][dir].setVal(0.0);

    //         a_vel[dit][dir].plus(oldVel[dit][dir],
    //                              1.0 - theta,
    //                              a_srcComp,
    //                              a_destComp,
    //                              a_numComp);

    //         a_vel[dit][dir].plus(newVel[dit][dir],
    //                              theta,
    //                              a_srcComp,
    //                              a_destComp,
    //                              a_numComp);
    //     }
    // }

    CH_assert(0.0 <= theta && theta <= 1.0);
    return theta;
}


// -----------------------------------------------------------------------------
Real
SDC::pLinearInterp(LevelData<FArrayBox>& a_p,
                   const Real            a_time,
                   int                   a_srcComp,
                   int                   a_destComp,
                   int                   a_numComp) const
{
    CH_assert(a_numComp <= m_pNumComps);
    CH_assert(a_p.getBoxes().compatible(m_grids));
    CH_assert(a_p.getBoxes().physDomain().size() ==
              m_grids.physDomain().size());

    // DataIterator dit = m_grids.dataIterator();
    // const LevelData<FArrayBox>& oldP = m_p[0];
    // const LevelData<FArrayBox>& newP = m_p[RKC::numStages - 1];

    // if (a_numComp == -1) a_numComp = m_pNumComps - a_srcComp;

    // // If we are at m_oldTime, just copy the old data.
    // CH_assert(m_Q0Ready);
    // if (abs(a_time - m_oldTime) < timeEps) {
    //     for (dit.reset(); dit.ok(); ++dit) {
    //         a_p[dit].copy(oldP[dit], a_srcComp, a_destComp, a_numComp);
    //     }
    //     return 0.0;
    // }

    // // Linear interpolation...
    // CH_assert(m_interpDataReady);

    // Compute theta and if nearby, snap into 0 or 1 values to prevent
    // extrapolation when not expected.
    const Real theta = this->computeTheta(a_time);

    // if (abs(theta) < timeEps) {
    //     // A bit redundant...
    //     for (dit.reset(); dit.ok(); ++dit) {
    //         a_p[dit].copy(oldP[dit], a_srcComp, a_destComp, a_numComp);
    //     }
    //     return 0.0;
    // } else if (abs(theta - 1.0) < timeEps) {
    //     for (dit.reset(); dit.ok(); ++dit) {
    //         a_p[dit].copy(newP[dit], a_srcComp, a_destComp, a_numComp);
    //     }
    //     return 1.0;
    // }

    // // Perform interpolation.
    // for (dit.reset(); dit.ok(); ++dit) {
    //     a_p[dit].setVal(0.0);
    //     a_p[dit].plus(oldP[dit], 1.0-theta, a_srcComp, a_destComp, a_numComp);
    //     a_p[dit].plus(newP[dit], theta, a_srcComp, a_destComp, a_numComp);
    // }

    CH_assert(0.0 <= theta && theta <= 1.0);
    return theta;
}


// -----------------------------------------------------------------------------
Real
SDC::qLinearInterp(LevelData<FArrayBox>& a_q,
                   const Real            a_time,
                   int                   a_srcComp,
                   int                   a_destComp,
                   int                   a_numComp) const
{
    CH_assert(a_numComp <= m_qNumComps);
    CH_assert(a_q.getBoxes().compatible(m_grids));
    CH_assert(a_q.getBoxes().physDomain().size() ==
              m_grids.physDomain().size());

    // DataIterator dit = m_grids.dataIterator();
    // const LevelData<FArrayBox>& oldQ = m_q[0];
    // const LevelData<FArrayBox>& newQ = m_q[RKC::numStages - 1];

    // if (a_numComp == -1) a_numComp = m_qNumComps - a_srcComp;

    // // If we are at m_oldTime, just copy the old data.
    // CH_assert(m_Q0Ready);
    // if (abs(a_time - m_oldTime) < timeEps) {
    //     for (dit.reset(); dit.ok(); ++dit) {
    //         a_q[dit].copy(oldQ[dit], a_srcComp, a_destComp, a_numComp);
    //     }
    //     return 0.0;
    // }

    // // Linear interpolation...
    // CH_assert(m_interpDataReady);

    // Compute theta and if nearby, snap into 0 or 1 values to prevent
    // extrapolation when not expected.
    const Real theta = this->computeTheta(a_time);

    // if (abs(theta) < timeEps) {
    //     // A bit redundant...
    //     for (dit.reset(); dit.ok(); ++dit) {
    //         a_q[dit].copy(oldQ[dit], a_srcComp, a_destComp, a_numComp);
    //     }
    //     return 0.0;
    // } else if (abs(theta - 1.0) < timeEps) {
    //     for (dit.reset(); dit.ok(); ++dit) {
    //         a_q[dit].copy(newQ[dit], a_srcComp, a_destComp, a_numComp);
    //     }
    //     return 1.0;
    // }

    // // Perform interpolation.
    // for (dit.reset(); dit.ok(); ++dit) {
    //     a_q[dit].setVal(0.0);
    //     a_q[dit].plus(oldQ[dit], 1.0-theta, a_srcComp, a_destComp, a_numComp);
    //     a_q[dit].plus(newQ[dit], theta, a_srcComp, a_destComp, a_numComp);
    // }

    CH_assert(0.0 <= theta && theta <= 1.0);
    return theta;
}


// -----------------------------------------------------------------------------
void
SDC::setOldQ(const LevelData<FluxBox>&   a_velOld,
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
    for (size_t m = 0; m < NumNodes; ++m) {
        debugInitLevel(m_vel[m]);
        debugInitLevel(m_p[m]);
        debugInitLevel(m_q[m]);

        debugInitLevel(m_kvelE[m]);
        debugInitLevel(m_kqE[m]);

        debugInitLevel(m_kvelI[m]);
        debugInitLevel(m_kqI[m]);

        debugInitLevel(m_kvelOld[m]);
        debugInitLevel(m_kqOld[m]);
    }

    m_oldTime = quietNAN;
    m_newTime = quietNAN;
    m_dt      = quietNAN;
#endif

    for (DataIterator dit(m_grids); dit.ok(); ++dit) {
        for (int dir = 0; dir < SpaceDim; ++dir) {
            m_vel[0][dit][dir].copy(a_velOld[dit][dir]);
        }
        m_p[0][dit].copy(a_pOld[dit]);
        m_q[0][dit].copy(a_qOld[dit]);
    }
    m_oldTime = a_oldTime;

    m_Q0Ready          = true;
    m_finalForcesReady = false;
    m_interpDataReady  = false;
}


// -----------------------------------------------------------------------------
void
SDC::copy(LevelData<FluxBox>&       a_x,
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
void
SDC::copy(LevelData<FArrayBox>&       a_x,
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
void
SDC::plus(LevelData<FluxBox>&       a_x,
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
void
SDC::plus(LevelData<FArrayBox>&       a_x,
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
void
SDC::axby(const Real                a_a,
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
void
SDC::axby(const Real                  a_a,
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
