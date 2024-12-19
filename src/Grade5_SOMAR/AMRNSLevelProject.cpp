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
#include "AMRNSLevel.H"
#include "Analysis.H"
#include "Integral.H"
#include "MiscUtils.H"
#include <chrono>


// TEMP!!!
#include "LevelLepticSolver.H"
#include "memusage.H"


#ifndef NDEBUG
    // Debug mode
#   define nanCheck(x) checkForValidNAN(x)
#else
    // Release mode
#   define nanCheck(x)
#endif

#define WRITE_DIAGNOSTIC_HDF5(d, lg) \
    {                                \
    }                                \
    while (0)

/*
 #define WRITE_DIAGNOSTIC_HDF5(d, lg)         \
     {                                        \
         static int idx = 0;                  \
         char       fn[80];                   \
         sprintf(fn, #d "_%04d.hdf5", idx++); \
         IO::writeHDF5(fn, d, lg)             \
     }                                        \
     while (0)
*/


// -----------------------------------------------------------------------------
void
AMRNSLevel::projectPredict(LevelData<FluxBox>&   a_vel,
                           LevelData<FArrayBox>& a_p,
                           const Real            a_time,
                           const Real            a_projDt)
{
    BEGIN_FLOWCHART();

    const ProblemContext* ctx = ProblemContext::getInstance();
    if (!ctx->proj.doLevelProj) return;

    // Gather references / set parameters
    const bool               writeToPout      = (s_verbosity >= 1);
    const int                normType         = ctx->proj.normType;
    const DisjointBoxLayout& grids            = m_levGeoPtr->getBoxes();
    const Real               pressureTime     = a_time;  // - 0.5 * a_projDt;
    constexpr bool           unprojectCrseVel = false;

    // Prepare data on the coarser level.
    LevelData<FArrayBox>* crsePPtr       = nullptr;
    LevelData<FluxBox>*   crseVelStarPtr = nullptr;
    if (m_level > 0) {
        const AMRNSLevel*        crsePtr   = this->crseNSPtr();
        const DisjointBoxLayout& crseGrids = crsePtr->getBoxes();

        crsePPtr = new LevelData<FArrayBox>(crseGrids, 1, IntVect::Unit);
        crsePtr->fillPressure(*crsePPtr, pressureTime, true);

        constexpr bool setBCs = !unprojectCrseVel;
        crseVelStarPtr = new LevelData<FluxBox>(crseGrids, 1, IntVect::Unit);
        crsePtr->fillVelocity(*crseVelStarPtr, a_time, setBCs);

        if (unprojectCrseVel) {
            LevelData<FArrayBox>* crseCrsePPtr = nullptr; // TEMP!!!
            if (m_level > 1) {
                const AMRNSLevel*        crseCrsePtr   = crsePtr->crseNSPtr();
                const DisjointBoxLayout& crseCrseGrids = crseCrsePtr->getBoxes();

                crseCrsePPtr = new LevelData<FArrayBox>(crseCrseGrids, 1);
                crseCrsePtr->fillPressure(*crseCrsePPtr, pressureTime, false);
            }

            LevelData<FluxBox> crseGradP(crseGrids, 1);
            crsePtr->m_projOpPtr->levelGradient(
                crseGradP, *crsePPtr, crseCrsePPtr, pressureTime, false, false);

            for (DataIterator dit(crseGrids); dit.ok(); ++dit) {
                for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                    auto&       crseVelStarFAB = (*crseVelStarPtr)[dit][velComp];
                    const auto& crseGradPFAB   = crseGradP[dit][velComp];

                    crseVelStarFAB.plus(crseGradPFAB, a_projDt);
                }
            }

            BCTools::extrapAllGhosts(*crseVelStarPtr, 2);
            crseVelStarPtr->exchange();
        }
    }

    // Prepare pressure on this level.
    this->setPressureBC(a_p, pressureTime, false, crsePPtr);

    // Prepare advecting velocity on this level.
    this->setVelBC(a_vel, a_time, false, crseVelStarPtr);
    this->sendToAdvectingVelocity(a_vel, a_vel);
    nanCheck(a_vel);

    // Free memory
    delete crseVelStarPtr;
    crseVelStarPtr = nullptr;

    // Compute initial divergence and diagnostics.
    LevelData<FArrayBox> divVel(grids, 1);
    CH_assert(m_projOpPtr);
    m_projOpPtr->levelDivergence(divVel, a_vel);
    nanCheck(divVel);

    Real vol;
    const Real initDivNorm = m_projOpPtr->norm(divVel, normType);
    const Real initDivInt  = Integral::sum(vol, divVel, *m_levGeoPtr);
    const Real initBdryInt = Integral::bdrySum(a_vel, *m_levGeoPtr);

    // Compute Grad[p].
    LevelData<FluxBox> gradP(grids, 1);
    m_projOpPtr->levelGradient(gradP, a_p, crsePPtr, pressureTime, false, false);
    nanCheck(gradP);

    // Free memory
    delete crsePPtr;
    crsePPtr = nullptr;

    // Project.
    for (DataIterator dit(grids); dit.ok(); ++dit) {
        D_TERM(
        a_vel[dit][0].plus(gradP[dit][0], -a_projDt);,
        a_vel[dit][1].plus(gradP[dit][1], -a_projDt);,
        a_vel[dit][2].plus(gradP[dit][2], -a_projDt);)
    }
    nanCheck(a_vel);

    // // Set projected vel BCs.
    // this->sendToCartesianVelocity(a_vel, a_vel);
    // this->setVelBC(a_vel, a_time, false);
    // this->sendToAdvectingVelocity(a_vel, a_vel);
    // nanCheck(a_vel);

    // Compute final divergence and diagnostics.
    m_projOpPtr->levelDivergence(divVel, a_vel);
    nanCheck(divVel);

    const Real finalDivNorm = m_projOpPtr->norm(divVel, normType);
    const Real finalDivInt  = Integral::sum(divVel, *m_levGeoPtr);
    const Real finalBdryInt = Integral::bdrySum(a_vel, *m_levGeoPtr);

    if (finalDivNorm > initDivNorm) {
        // Lagged pressure failed. Use one FMG iter to drive divVel down.
        const auto saveOpts = m_levelProjSolverPtr->getOptions();
        m_levelProjSolverPtr->modifyOptionsExceptMaxDepth(LevelProjSolver::getQuickAndDirtyOptions());
        this->projectCorrect(a_vel, a_p, a_time, a_projDt);
        m_levelProjSolverPtr->modifyOptionsExceptMaxDepth(saveOpts);

        m_projOpPtr->levelDivergence(divVel, a_vel);
        nanCheck(divVel);
        const Real finalDivNorm2 = m_projOpPtr->norm(divVel, normType);

        if (writeToPout) {
            pout() << "MAC lagged projection:"
                   << Format::pushFlags
                   << Format::indent()
                   << Format::scientific;
            pout() << "\npre  |Div[vel]|_" << normType << " = " << initDivNorm
                   << "\npost |Div[vel]|_" << normType << " = " << finalDivNorm2
                   << '\n';
            pout() << Format::unindent
                   << Format::popFlags
                   << std::flush;
        }
    } else {
        // Write diagnostics.
        if (writeToPout) {
            bool writeIntegrals = false;
            writeIntegrals |= (initBdryInt >= 1.0e-6);
            writeIntegrals |= (initDivInt >= 1.0e-6);
            writeIntegrals |= (finalBdryInt >= 1.0e-6);
            writeIntegrals |= (finalDivInt >= 1.0e-6);
            writeIntegrals |= (finalDivNorm >= initDivNorm / 10.0);

            pout() << "MAC lagged projection:"
                   << Format::pushFlags
                   << Format::indent()
                   << Format::scientific;

            // if (writeIntegrals) {
            //     pout() << "\npre  Integral[vel.n dS]  = " << initBdryInt
            //            << "\npre  Integral[DivVel dV] = " << initDivVelInt
            //            << "\npost Integral[vel.n dS]  = " << finalBdryInt
            //            << "\npost Integral[DivVel dV] = " << finalDivVelInt;
            // }

            pout() << "\npre  |Div[vel]|_" << normType << " = " << initDivNorm
                   << "\npost |Div[vel]|_" << normType << " = " << finalDivNorm
                   << "\n";

            pout() << Format::unindent
                   << Format::popFlags
                   << std::flush;
        }
    }


    // // Set projected vel BCs.
    // this->sendToCartesianVelocity(a_vel, a_vel);
    // this->setVelBC(a_vel, a_time, false);
    // this->sendToAdvectingVelocity(a_vel, a_vel);
    // nanCheck(a_vel);

    this->sendToCartesianVelocity(a_vel, a_vel);
    this->setVelBC(a_vel, a_time, false);
}


// -----------------------------------------------------------------------------
void
AMRNSLevel::projectCorrect(LevelData<FluxBox>&   a_vel,
                           LevelData<FArrayBox>& a_p,
                           const Real            a_time,
                           const Real            a_projDt,
                           const Real            a_gammaDt)
{
    BEGIN_FLOWCHART();

    const ProblemContext* ctx = ProblemContext::getInstance();

    // ====== Projection ======
    // NaN checks in debug mode no matter what.
    nanCheck(a_vel);
    nanCheck(a_p);

    // Set velocity BCs no matter what.
    this->setVelBC(a_vel, a_time, false);
    nanCheck(a_vel);

    // Are we skipping the projection? (Not typical. Used when debugging.)
    if (!ctx->proj.doLevelProj) return;

    // Gather references / set parameters.
    const bool               setPhiToZero  = true;
    const bool               writeToPout   = (s_verbosity >= 1);
    const int                normType      = ctx->proj.normType;
    const DisjointBoxLayout& grids         = m_levGeoPtr->getBoxes();
    DataIterator             dit           = grids.dataIterator();

    // Begin indentation...
    if (writeToPout) {
        pout() << "MAC projection:"
               << Format::pushFlags
               << Format::indent()
               << endl;
    }

    // Allocate workspace
    LevelData<FArrayBox> divVel(grids, 1);
    LevelData<FArrayBox> phi(grids, 1, IntVect::Unit);
    LevelData<FluxBox>   gradPhi(grids, 1);

    // Div[vel] can only be computed with the advecting velocity.
    this->sendToAdvectingVelocity(a_vel, a_vel);

    // Compute initial divergence and diagnostics.
    m_projOpPtr->levelDivergence(divVel, a_vel);
    nanCheck(divVel);

    // Real vol;
    const Real initDivNorm = m_projOpPtr->norm(divVel, normType);
    // const Real initDivInt  = Integral::sum(vol, divVel, *m_levGeoPtr);
    // const Real initBdryInt = Integral::bdrySum(a_vel, *m_levGeoPtr);

    // Project...
    {
        if (writeToPout) {
            pout() << "MG solver:" << Format::indent() << endl;
        }

        if (m_level == 0) {
            pout() << Format::scientific;
            m_projOpPtr->setToZero(phi);
            m_projOpPtr->levelEquationIsConsistent(
                phi, nullptr, divVel, a_time, true, true);
        }

        // Solve L[phi] = Div[vel*].
        const auto start = std::chrono::high_resolution_clock::now();
        m_levelProjSolverPtr->solve(phi, nullptr, divVel, a_time, true, setPhiToZero);
        const auto finish = std::chrono::high_resolution_clock::now();

        // Write solve time.
        if (writeToPout) {
            std::chrono::duration<double> elapsed = finish - start;
            pout() << "Solve time: " << elapsed.count() << " s"
                   << Format::unindent << endl;
        }

        // Compute Grad[p].
        nanCheck(phi);
        m_projOpPtr->levelGradient(gradPhi, phi, nullptr, a_time, true, true);
        nanCheck(gradPhi);

        // Project a_vel and update a_p.
        for (dit.reset(); dit.ok(); ++dit) {
            D_TERM(
            a_vel[dit][0].plus(gradPhi[dit][0], -1.0);,
            a_vel[dit][1].plus(gradPhi[dit][1], -1.0);,
            a_vel[dit][2].plus(gradPhi[dit][2], -1.0);)

            a_p[dit].plus(phi[dit], 1.0 / a_projDt);
            // if (!RealCmp::isZero(a_gammaDt) && ctx->rhs.doViscousForcing &&
            //     ctx->rhs.doImplicitDiffusion) {
            //     a_p[dit].plus(divVel[dit], a_gammaDt * ctx->rhs.nu / a_projDt);
            // }
        }
        nanCheck(a_vel);
        nanCheck(a_p);
    }

    // Compute final divergence and diagnostics.
    m_projOpPtr->levelDivergence(divVel, a_vel);
    nanCheck(divVel);

    const Real finalDivNorm = m_projOpPtr->norm(divVel, normType);
    // const Real finalDivInt  = Integral::sum(divVel, *m_levGeoPtr);
    // const Real finalBdryInt = Integral::bdrySum(a_vel, *m_levGeoPtr);

    // Write diagnostics and restore indentation state.
    if (writeToPout) {
        pout() << Format::scientific << "pre  |Div[vel]|_" << normType << " = "
               << initDivNorm << "\npost |Div[vel]|_" << normType << " = "
               << finalDivNorm << std::endl
               << Format::unindent
               << Format::popFlags
               << std::flush;
    }

    // Restore Cartesian velocity.
    this->sendToCartesianVelocity(a_vel, a_vel);
}


// -----------------------------------------------------------------------------
void
AMRNSLevel::projectDownToThis(const bool a_sync)
{
#if 0
    // This is code that I used to test the line relaxation method.
    // Will delete in a future push.
    const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();

    LevelData<FArrayBox> phi(grids, 1, IntVect::Unit);
    LevelData<FArrayBox> sol(grids, 1, IntVect::Unit);
    LevelData<FArrayBox> rhs(grids, 1);
    for (DataIterator dit(grids); dit.ok(); ++dit) {
        const Box valid = grids[dit];

        const RealVect L = m_levGeoPtr->getDomainLength();
        const RealVect kk(D_DECL(16.0, 16.0, 1.0));
        // const RealVect kk(D_DECL(16.0, 1.0, 1.0));

        const RealVect k = kk * 2.0 * Pi / L;

        FArrayBox xyzFAB(valid, SpaceDim);
        m_levGeoPtr->fill_physCoor(xyzFAB);

        for (BoxIterator bit(valid); bit.ok(); ++bit) {
            const IntVect& cc = bit();
            const Real z = xyzFAB(cc,2) / L[SpaceDim - 1];

            // sol[dit](cc) = D_TERM(
            //                cos(k[0] * xyzFAB(cc,0)),
            //              * cos(k[1] * xyzFAB(cc,1)),
            //              * cos(k[2] * xyzFAB(cc,2)));
            sol[dit](cc) = D_TERM(
                           cos(k[0] * xyzFAB(cc,0)),
                         * cos(k[1] * xyzFAB(cc,1)),
                         * 16.0 * z * z * (z + 1.0) * (z + 1.0));
        }
    }

    // std::shared_ptr<const Elliptic::LepticOperator<LevelData<FArrayBox>>> lepticOpPtr(
    //     new Elliptic::PoissonOp(*m_levGeoPtr, this->getCrseGridsPtr(), 1));
    // const auto* opPtr = &lepticOpPtr->getLevelOp();
    auto opPtr = m_projOpPtr;

    constexpr bool homogPhysBCs = true;
    constexpr bool homogCFIBCs  = true;
    constexpr bool setPhiToZero = true;
    // opPtr->assignLocal(phi, sol);

    opPtr->applyOp(rhs, sol, nullptr, 0.0, homogPhysBCs, homogCFIBCs);


    // auto& solver = *m_levelProjSolverPtr;

    Elliptic::LevelLepticSolver solver;
    solver.define(opPtr);

    // Elliptic::LepticSolver<LevelData<FArrayBox>> solver;
    // solver.define(lepticOpPtr);
    // // solver.options().convergenceMetric = -1.0;
    // solver.options().minIters = 10;
    // solver.options().maxIters = 100;
    // solver.options().normType = 2;
    // solver.options().verbosity = 5;

    // Elliptic::MGSolver<LevelData<FArrayBox>> solver;
    // {
    //     auto opts = solver.getOptions();
    //     opts.bottomOptions.absTol = 1.0e-12;
    //     opts.bottomOptions.relTol = 1.0e-12;
    //     // opts.bottomOptions.verbosity = 5;
    //     opts.bottomOptions.normType = 2;
    //     opts.bottomOptions.numSmoothPrecond = 0;
    //     // opts.verbosity = 5;
    //     opts.absTol = 1.0e-12;
    //     opts.relTol = 1.0e-12;
    //     opts.numSmoothDown = 2;
    //     opts.numSmoothUp = 2;
    //     opts.numSmoothPrecond = 2;
    //     opts.numSmoothBottom = 2;
    //     // opts.numSmoothUpFMG = 16;
    //     opts.maxDepth = 0;
    //     opts.numCycles = 1;
    //     opts.normType = 2;
    //     solver.define(*opPtr, opts);
    // }

    // Elliptic::BiCGStabSolver<LevelData<FArrayBox>> solver;
    // solver.define(opPtr);
    // solver.options().numSmoothPrecond = 1;
    // solver.options().verbosity = 2;
    // solver.options().normType = 2;

    // Elliptic::RelaxSolver<LevelData<FArrayBox>> solver;
    // solver.define(opPtr);
    // solver.options().normType = 2;
    // solver.options().itersPerNormCheck = 1;
    // solver.options().minIters = 10;
    // solver.options().maxIters = 100;
    // solver.options().verbosity = 5;
    // // solver.options().convergenceMetric = 1.0;


    // opPtr->preCond(phi, rhs, 0.0, solver.getOptions().numSmoothPrecond);

    {
        pout() << Format::scientific;
        m_projOpPtr->setToZero(phi);
        m_projOpPtr->levelEquationIsConsistent(
            phi, nullptr, rhs, 0.0, true, true);
    }


    const auto start = std::chrono::high_resolution_clock::now();
    solver.solve(phi, nullptr, rhs, 0.0, homogPhysBCs, setPhiToZero);
    const auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    pout() << "Solve time: " << elapsed.count() << " s" << endl;
    IO::tout(0) << "Solve time: " << elapsed.count() << " s" << endl;


    for (DataIterator dit(grids); dit.ok(); ++dit) {
        const Box valid = grids[dit];
        phi[dit] -= sol[dit];
    }
    const Real errNorm0 = Analysis::pNorm(phi, 0);
    const Real errNorm2 = Analysis::pNorm(phi, 2);
    IO::tout(0) << "|error|_oo = " << errNorm0 << '\n';
    IO::tout(0) << "|error|_2  = " << errNorm2 << '\n';
    pout() << "|error|_oo = " << errNorm0 << '\n';
    pout() << "|error|_2  = " << errNorm2 << '\n';
    pout() << std::endl;

    IO::writeHDF5("phi.hdf5", phi, *m_levGeoPtr, 0.0, 0.0);
    barrier();
    CH_verify(false);
#endif


#if 0
    // This is code that I used to test the horizontal solver.
    // Will delete in a future push.

    std::shared_ptr<const Elliptic::LepticOperator<LevelData<FArrayBox>>> opPtr(
        new Elliptic::PoissonOp(*m_levGeoPtr, this->getCrseGridsPtr(), 1));
    const auto  hOpPtr   = opPtr->createHorizontalOperator();

    const auto& op       = *opPtr;
    const auto& levelOp  = opPtr->getLevelOp();
    const auto& hOp      = *hOpPtr;
    const auto& hLevelOp = hOp.getLevelOp();
    const auto& hMGOpPtr = hOp.getMGOp();

    auto hLevelOpPtr =
        std::shared_ptr<const Elliptic::LevelOperator<LevelData<FArrayBox>>>(
            &hLevelOp, &noDelete2);


    const DisjointBoxLayout& grids  = m_levGeoPtr->getBoxes();
    const ProblemDomain&     domain = grids.physDomain();
    const RealVect           L      = m_levGeoPtr->getDomainLength();

    DisjointBoxLayout hGrids;
    {
        bool hPeriodic[SpaceDim];
        D_TERM(hPeriodic[SpaceDim - 1] = false;,
            hPeriodic[0] = domain.isPeriodic(0);,
            hPeriodic[1] = domain.isPeriodic(1);)

        ProblemDomain horizDomain(
            Subspace::flattenBox(domain.domainBox(), SpaceDim - 1), hPeriodic);

        LayoutTools::RangeTransform rangeTransform(SpaceDim - 1, 0);
        hGrids.deepCopy(grids);
        hGrids.transform(rangeTransform);
        hGrids.setDomain(horizDomain);
        hGrids.close();
        CH_verify(hGrids.compatible(grids));
    }
    const ProblemDomain& hDomain = hGrids.physDomain();
    const RealVect hL = [&L]() {
        RealVect val = L;
        val[SpaceDim - 1] = 1.0;
        return val;
    }();

    LevelData<FArrayBox> hPhi(hGrids, 1, IntVect::Unit);
    LevelData<FArrayBox> hSol(hGrids, 1, IntVect::Unit);
    LevelData<FArrayBox> hRhs(hGrids, 1);
    LevelData<FArrayBox> hLPhi(hGrids, 1);

    for (DataIterator dit(hGrids); dit.ok(); ++dit) {
        const Box hValid = hGrids[dit];

        const RealVect kk(D_DECL(16.0, 16.0, 0.0));
        const RealVect k = kk * 2.0 * Pi / L;
        const Real lambda = -k.dotProduct(k);

        FArrayBox hxyzFAB(hValid, SpaceDim);
        m_levGeoPtr->fill_physCoor(hxyzFAB);

        for (BoxIterator hbit(hValid); hbit.ok(); ++hbit) {
            const IntVect& hcc = hbit();
            const Real x = hxyzFAB(hcc, 0);
            const Real y = hxyzFAB(hcc, 1);

            hSol[dit](hcc) = D_TERM(,
                             cos(k[0] * x),
                           * cos(k[1] * y));

            hRhs[dit](hcc) = lambda * hSol[dit](hcc);
        }
    }

    constexpr bool homogPhysBCs = true;
    constexpr bool homogCFIBCs  = true;
    constexpr bool setPhiToZero = true;
    // opPtr->assignLocal(phi, sol);

    hLevelOp.applyOp(hLPhi, hSol, nullptr, 0.0, homogPhysBCs, homogCFIBCs);
    // hLevelOp.setToZero(hPhi);
    // hLevelOp.assignLocal(hPhi, hSol);

    // // Check if applyOp is doing as expected.
    // hLevelOp.incr(hRhs, hLPhi, -1.0);
    // {
    //     const Real norm0 = Analysis::pNorm(hRhs, 0);
    //     IO::tout(0) << "L[phi] inf-norm = " << norm0 << endl;
    // }

    // Elliptic::LepticSolver<LevelData<FArrayBox>> solver;
    // solver.define(opPtr);
    // solver.options().minIters = 10;
    // solver.options().maxIters = 100;
    // solver.options().normType = 2;
    // solver.options().verbosity = 5;

    Elliptic::MGSolver<LevelData<FArrayBox>> solver;
    {
        auto opts = solver.getOptions();
        opts.bottomOptions.absTol = 1.0e-12;
        opts.bottomOptions.relTol = 1.0e-12;
        // opts.bottomOptions.verbosity = 5;
        opts.bottomOptions.normType = 2;
        opts.bottomOptions.numSmoothPrecond = 0;
        // opts.verbosity = 5;
        opts.absTol = 1.0e-12;
        opts.relTol = 1.0e-12;
        opts.numSmoothDown = 4;
        opts.numSmoothUp = 4;
        opts.numSmoothPrecond = 0;
        opts.numSmoothBottom = 0;
        // opts.numSmoothUpFMG = 16;
        // opts.maxDepth = -1;
        // opts.numCycles = 1;
        opts.normType = 2;
        solver.define(hMGOpPtr, opts);
    }

    // Elliptic::BiCGStabSolver<LevelData<FArrayBox>> solver;
    // solver.define(hLevelOpPtr);
    // solver.options().absTol = 1.0e-12;
    // solver.options().relTol = 1.0e-12;
    // solver.options().numSmoothPrecond = 10;
    // solver.options().verbosity = 4;
    // solver.options().normType = 2;

    // Elliptic::RelaxSolver<LevelData<FArrayBox>> solver;
    // solver.define(hLevelOpPtr);
    // solver.options().absTol = 1.0e-12;
    // solver.options().relTol = 1.0e-12;
    // solver.options().normType = 2;
    // solver.options().itersPerNormCheck = 1;
    // solver.options().minIters = 10;
    // solver.options().maxIters = 100;
    // solver.options().verbosity = 5;
    // // solver.options().convergenceMetric = 1.0;


    // opPtr->preCond(phi, rhs, 0.0, solver.getOptions().numSmoothPrecond);


    const auto start = std::chrono::high_resolution_clock::now();
    {
        solver.solve(hPhi, nullptr, hLPhi, 0.0, homogPhysBCs, setPhiToZero);
    }
    const auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    pout() << "Solve time: " << elapsed.count() << " s" << endl;
    IO::tout(0) << "Solve time: " << elapsed.count() << " s" << endl;

    hLevelOp.incr(hPhi, hSol, -1.0);
    const Real errNorm0 = Analysis::pNorm(hPhi, 0);
    const Real errNorm2 = Analysis::pNorm(hPhi, 2);
    IO::tout(0) << "|error|_oo = " << errNorm0 << '\n';
    IO::tout(0) << "|error|_2  = " << errNorm2 << '\n';
    pout() << "|error|_oo = " << errNorm0 << '\n';
    pout() << "|error|_2  = " << errNorm2 << '\n';
    pout() << std::endl;

    LevelGeometry hLevGeo(hDomain, hL, nullptr, &m_levGeoPtr->getGeoSource());
    hLevGeo.createMetricCache(hGrids);

    IO::writeHDF5("phi.hdf5", hPhi, hLevGeo, 0.0, 0.0);
    barrier();
    CH_verify(false);
#endif




    BEGIN_FLOWCHART();

    // Gather references / set parameters.
    const ProblemContext* ctx          = ProblemContext::getInstance();
    const bool            writeToPout  = (s_verbosity >= 1);
    const int             normType     = ctx->proj.normType;

    const size_t lbase = (m_level > 0)? (m_level - 1): 0;
    const size_t lmin = m_level;
    const size_t lmax = this->finestNSPtr()->m_level;

    const Real syncDt = this->m_dt; // lmin's dt.

    // Begin intentation...
    if (writeToPout) {
        pout() << "AMR proj on levels " << lmin << " to " << lmax << ":"
               << Format::indent() << Format::pushFlags << endl;
    }

    // Gather velocity references, create workspace.
    Vector<LevelData<FluxBox>*>   amrVel;
    Vector<LevelData<FArrayBox>*> amrPhi;
    Vector<LevelData<FluxBox>*>   amrGradPhi;
    Vector<LevelData<FArrayBox>*> amrRhs;

    this->allocateAndDefine(amrPhi, 1, IntVect::Unit, lbase, lmax);
    this->allocateAndAliasVel(amrVel, m_time, Interval(), lmin, lmax);
    this->allocateAndDefine(amrGradPhi, 1, IntVect::Zero, lmin, lmax);
    this->allocateAndDefine(amrRhs, 1, IntVect::Zero, lmin, lmax);

    // Initial / base level phi is zero.
    setValLevels(amrPhi, lbase, lmax, 0.0);

    // Base level velocity needs to be interpolated in time.
    if (lbase < lmin) {
        const AMRNSLevel*        basePtr   = this->getLevel(lbase);
        const DisjointBoxLayout& baseGrids = basePtr->getBoxes();
        const int                numComps  = amrVel[lmin]->nComp();
        const IntVect&           ghostVect = amrVel[lmin]->ghostVect();

        amrVel[lbase] = new LevelData<FluxBox>(baseGrids, numComps, ghostVect);
        basePtr->fillVelocity(*amrVel[lbase], m_time, false);
    }

    // Send to advecting velocity
    for (size_t lev = lbase; lev <= lmax; ++lev) {
        const AMRNSLevel*   levPtr = this->getLevel(lev);
        LevelData<FluxBox>& vel    = *amrVel[lev];

        levPtr->sendToAdvectingVelocity(vel, vel);
    }

    // Compute initial divergence.
    for (size_t lev = lmin; lev <= lmax; ++lev) {
        const AMRNSLevel* levPtr = this->getLevel(lev);
        const auto        opPtr  = levPtr->m_projOpPtr;

        const LevelData<FluxBox>* fineVelPtr = nullptr;
        if (lev < lmax) fineVelPtr = amrVel[lev + 1];

        opPtr->compDivergence(*amrRhs[lev], *amrVel[lev], fineVelPtr);
    }
    const Real initDivNorm = Analysis::pNorm(amrRhs, normType);
    // this->m_levelProjectorOpPtr->AMRNorm(amrDivVel, normType);

    // Check eq. constistency.
    this->averageDown(amrRhs, lmin, false);
    {
        const AMRNSLevel* levPtr  = this->getLevel(lmin);
        auto              op      = levPtr->m_projOpPtr;
        auto&             phi     = *amrPhi[lmin];
        auto&             gradPhi = *amrGradPhi[lmin];
        auto&             rhs     = *amrRhs[lmin];
        auto&             levGeo  = *levPtr->m_levGeoPtr;

        pout() << Format::scientific;
        op->setToZero(phi);

        // const bool consistent = op->levelEquationIsConsistent(
        //     phi, nullptr, rhs, m_time, true, true);
        // if (!consistent) {
        //     pout() << "level " << lmin << " equation is NOT consistent.";
        // }

        op->compGradient(gradPhi, phi, nullptr, m_time, true, true);
        const Real bdrySum = Integral::bdrySum(gradPhi, levGeo);
        const Real rhsSum  = Integral::sum(rhs, levGeo, false);
        const bool consistent =
            (abs(bdrySum - rhsSum) < m_amrProjSolverPtr->getOptions().relTol);

        if (writeToPout || !consistent) {
            pout() << "bdrySum = " << bdrySum
                << "\nrhsSum  = " << rhsSum
                << endl;
        }
        if (!consistent) {
            std::ostringstream msg;
            msg << "level " << lmin
                << " equation is not consistent BEFORE solve.";
            pout() << msg.str() << endl;
            // IO::tout(0) << Format::yellow << msg.str() << Format::none << endl;
        }
    }

    // Solve L[phi] = Div[vel*].
    {
        const bool setPhitoZero = true;
        const bool useHomogBCs  = true;

        auto start = std::chrono::high_resolution_clock::now();
        m_amrProjSolverPtr->solve(amrPhi,
                                  makeConstVector(amrRhs),
                                  m_time,
                                  useHomogBCs,
                                  setPhitoZero);
        auto finish = std::chrono::high_resolution_clock::now();

        // Write solve time.
        if (writeToPout) {
            std::chrono::duration<double> elapsed = finish - start;
            pout() << "Solve time: " << elapsed.count() << " s" << endl;
        }
    }

    // Compute Grad[phi] and project.
    for (size_t lev = lmin; lev <= lmax; ++lev) {
        // Gather references.
        const AMRNSLevel*     levPtr  = this->getLevel(lev);
        const auto            opPtr   = levPtr->m_projOpPtr;
        LevelData<FluxBox>&   gradPhi = *amrGradPhi[lev];
        LevelData<FArrayBox>& phi     = *amrPhi[lev];
        LevelData<FluxBox>&   vel     = *amrVel[lev];
        DataIterator          dit     = vel.dataIterator();

        const LevelData<FArrayBox>* crsePhiPtr = nullptr;
        if (lev > 0) crsePhiPtr = amrPhi[lev - 1];

        // Grad[phi]
        opPtr->compGradient(gradPhi, phi, crsePhiPtr, m_time, true, false);

        // Project.
        for (dit.reset(); dit.ok(); ++dit) {
            D_TERM(
            vel[dit][0].plus(gradPhi[dit][0], -1.0);,
            vel[dit][1].plus(gradPhi[dit][1], -1.0);,
            vel[dit][2].plus(gradPhi[dit][2], -1.0);)
        }

        if (a_sync) {
            LevelData<FArrayBox>& pressure = levPtr->m_statePtr->p;
            for (dit.reset(); dit.ok(); ++dit) {
                pressure[dit].plus(phi[dit], 1.0 / syncDt);
            }
        }
    }

    // Compute final divergence.
    for (size_t lev = lmin; lev <= lmax; ++lev) {
        const AMRNSLevel* levPtr = this->getLevel(lev);
        const auto        opPtr  = levPtr->m_projOpPtr;

        LevelData<FluxBox>* fineVelPtr = nullptr;
        if (lev < lmax) fineVelPtr = amrVel[lev + 1];

        opPtr->compDivergence(*amrRhs[lev], *amrVel[lev], fineVelPtr);
    }
    const Real finalDivNorm = Analysis::pNorm(amrRhs, normType);

    // Check eq. constistency.
    this->averageDown(amrRhs, lmin, false);
    {
        const AMRNSLevel* levPtr  = this->getLevel(lmin);
        auto              op      = levPtr->m_projOpPtr;
        auto&             phi     = *amrPhi[lmin];
        auto&             gradPhi = *amrGradPhi[lmin];
        auto&             rhs     = *amrRhs[lmin];
        auto&             levGeo  = *levPtr->m_levGeoPtr;

        pout() << Format::scientific;
        // const bool consistent = op->levelEquationIsConsistent(
        //     phi, nullptr, rhs, m_time, true, true);
        // if (!consistent) {
        //     pout() << "level " << lmin << " equation is NOT consistent.";
        // }

        op->compGradient(gradPhi, phi, nullptr, m_time, true, true);
        const Real bdrySum = Integral::bdrySum(gradPhi, levGeo);
        const Real rhsSum = Integral::sum(rhs, levGeo, false);
        pout() << "bdrySum = " << bdrySum
               << "\nrhsSum  = " << rhsSum
               << endl;
        if (abs(bdrySum - rhsSum) > m_amrProjSolverPtr->getOptions().relTol) {
            std::ostringstream msg;
            msg << "level " << lmin
                << " equation is not consistent AFTER solve.";
            pout() << msg.str() << endl;
            // IO::tout(0) << Format::yellow << msg.str() << Format::none << endl;
        }
    }

    // Send amrVel back to Cartesian velocity.
    for (size_t lev = lbase; lev <= lmax; ++lev) {
        const AMRNSLevel*   levPtr = this->getLevel(lev);
        LevelData<FluxBox>& vel    = *amrVel[lev];

        levPtr->sendToCartesianVelocity(vel, vel);
    }

    // Write diagnostics.
    if (writeToPout) {
        pout() << Format::scientific << "pre  |Div[vel]|_" << normType << " = "
               << initDivNorm << "\npost |Div[vel]|_" << normType << " = "
               << finalDivNorm << endl;
    }

    // Restore indentation state.
    if (writeToPout) {
        pout() << Format::popFlags << Format::unindent << std::flush;
    }

    // Free memory
    this->deallocate(amrRhs);
    this->deallocate(amrGradPhi);
    this->deallocate(amrPhi);
    this->deallocate(amrVel);
}
