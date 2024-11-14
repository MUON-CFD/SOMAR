#include "AMRNSLevel.H"
#include "AMRNSLevelF_F.H"
#include "CFInterpF_F.H"
#include "Subspace.H"
#include "lapack.H"
#include "Comm.H"
#include "ParmParse.H"
#include "MiscUtils.H"
#include "Debug.H"



// -----------------------------------------------------------------------------
void
AMRNSLevel::setStratification(std::vector<Real>& a_vTbar,
                              std::vector<Real>& a_vSbar,
                              std::vector<Real>& a_vz,
                              const Real         a_zmin,
                              const Real         a_zmax) const
{
    BEGIN_FLOWCHART();
    pout() << "Stratification parameters:\n" << Format::indent() << std::flush;

    const int stratType = queryParm<int>("strat", "type", 0);

    if (stratType == 0) {
        // No stratification.
        // Setting min and max values to 0 should interp to 0 everywhere.
        a_vz.resize(2);
        a_vTbar.resize(2);
        a_vSbar.resize(2);

        a_vz[0] = a_zmin;
        a_vTbar[0] = 0.0;
        a_vSbar[0] = 0.0;

        a_vz[1] = a_zmax;
        a_vTbar[1] = 0.0;
        a_vSbar[1] = 0.0;

    } else if (stratType == 1) {
        // Linear stratification
        const Real L    = a_zmax - a_zmin;
        const Real T0   = getParm<Real>("strat", "T0");
        const Real S0   = getParm<Real>("strat", "S0");
        const Real dTdz = getParm<Real>("strat", "dTdz");
        const Real dSdz = getParm<Real>("strat", "dSdz");

        a_vz.resize(2);
        a_vTbar.resize(2);
        a_vSbar.resize(2);

        a_vz[0] = a_zmin;
        a_vTbar[0] = T0;
        a_vSbar[0] = S0;

        a_vz[1] = a_zmax;
        a_vTbar[1] = T0 + L * dTdz;
        a_vSbar[1] = S0 + L * dSdz;

        pout() << "type = 1 (LINEAR)\n";
        pout() << "T0 = " << T0 << '\n';
        pout() << "dTdz = " << dTdz << '\n';
        pout() << "S0 = " << T0 << '\n';
        pout() << "dSdz = " << dTdz << '\n';

    } else if (stratType == 2) {
        // Tanh stratification
        const Real T0 = getParm<Real>("strat", "T0");
        const Real dT = getParm<Real>("strat", "dT");
        const Real S0 = getParm<Real>("strat", "S0");
        const Real dS = getParm<Real>("strat", "dS");
        const Real z0    = getParm<Real>("strat", "z0");
        const Real sigma = std::max(getParm<Real>("strat", "sigma"),
                                    0.1 * this->getMinAMRDx(SpaceDim - 1));

        Real dz = m_levGeoPtr->getDXi(SpaceDim - 1);
        Box vertBox = Subspace::verticalDataBox(this->getDomain());
        vertBox.grow(SpaceDim - 1, 1);

        // Refine until we have at least 10 cells within deltaz.
        constexpr int r = 2;
        const IntVect riv = IntVect::Unit + (r - 1) * BASISV(SpaceDim - 1);
        while (dz > 0.1 * sigma) {
            vertBox.refine(riv);
            dz /= Real(r);
        }

        const int kmin = vertBox.smallEnd(SpaceDim - 1);
        const int kmax = vertBox.  bigEnd(SpaceDim - 1);
        CH_assert((Real(kmin) + 0.5) * dz <= a_zmin);
        CH_assert((Real(kmax) + 0.5) * dz >= a_zmax);

        a_vz.reserve(vertBox.size(SpaceDim - 1));
        a_vTbar.reserve(vertBox.size(SpaceDim - 1));
        a_vSbar.reserve(vertBox.size(SpaceDim - 1));
        Real z, tanhArg;

        for (int k = kmin; k <= kmax; ++k) {
            z = (Real(k) + 0.5) * dz;
            a_vz.push_back(z);

            tanhArg = tanh((z - z0) / sigma);
            a_vTbar.push_back(T0 - dT * tanhArg);
            a_vSbar.push_back(S0 - dS * tanhArg);
        }

        pout() << "type = " << stratType << " (TANH)\n";
        pout() << "T0    = " << T0 << '\n';
        pout() << "dT    = " << dT << '\n';
        pout() << "S0    = " << S0 << '\n';
        pout() << "dS    = " << dS << '\n';
        pout() << "z0    = " << z0 << '\n';
        pout() << "sigma = " << sigma
               << " (possibly adjusted to 0.1*min(dz) avoid zero).";

    } else {
        MAYDAYERROR(
            "strat.type = "
            << stratType
            << " not recognized. Choose a pre-defined type or override "
               "AMRNSLevel::setStratification in your custom physics class.");
    }

    pout() << Format::unindent << std::endl;
}


// -----------------------------------------------------------------------------
void
AMRNSLevel::setStratificationMembers()
{
    if (m_bbarPtr) return;

    const int             zdir = SpaceDim - 1;
    const IntVect         ez   = BASISV(zdir);
    const ProblemContext* ctx  = ProblemContext::getInstance();

    if (m_level == 0) {
        m_hasStrat = (queryParm<int>("strat", "type", 0) != 0);
    } else {
        m_hasStrat = this->coarsestNSPtr()->m_hasStrat;
    }

    // Do we need to solve the splines?
    if (!m_bbarSplinePtr && m_hasStrat) {
        if (m_level > 0) {
            // Just use splines from coarser level.
            const AMRNSLevel* crsePtr = this->coarsestNSPtr();
            m_bbarSplinePtr = crsePtr->m_bbarSplinePtr;
            m_TbarSplinePtr = crsePtr->m_TbarSplinePtr;
            m_SbarSplinePtr = crsePtr->m_SbarSplinePtr;

        } else {

            // Get Tbar, Sbar, and zbar.
            std::vector<Real> vz(0);
            std::vector<Real> vTbar(0);
            std::vector<Real> vSbar(0);
            {
                Box vertBox = Subspace::verticalDataBox(this->getDomainBox());
                vertBox.surroundingNodes(SpaceDim - 1);

                FArrayBox zFAB(vertBox, 1);
                m_levGeoPtr->fill_physCoor(zFAB, 0, zdir);
                const Real zmin = zFAB.min();
                const Real zmax = zFAB.max();

                this->setStratification(vTbar, vSbar, vz, zmin, zmax);

                if (vz.size() < 2) {
                    MAYDAYERROR(
                        "Your setStratification function should set at least 2 "
                        "points so that SOMAR can interpolate from zmin = "
                        << zmin << " to zmax = " << zmax << ".");
                }

                // Sort the values so that z is a increasing function.
                const auto perm = sortPermutation(vz);
                applyPermutation(vz, perm);
                applyPermutation(vTbar, perm);
                applyPermutation(vSbar, perm);

                // No duplicates allowed.
                if (std::unique(vz.begin(), vz.end()) != vz.end()) {
                    MAYDAYERROR(
                        "Your setStratification function must supply unique z "
                        "values. SOMAR found multiple entries with the same z "
                        "value.");
                }

                // Provided data must span domain.
                if (vz.front() > zmin) {
                    MAYDAYERROR(
                        "Your setStratification function should set at least one "
                        "point at or below zmin = "
                        << zmin
                        << " so that SOMAR can properly interpolate throughout the "
                        "domain.");
                }

                if (vz.back() < zmax) {
                    MAYDAYERROR(
                        "Your setStratification function should set at least one "
                        "point at or above zmax = "
                        << zmax
                        << " so that SOMAR can properly interpolate throughout the "
                        "domain.");
                }
            }
            const size_t Nz = vz.size();

            // Compute bbar
            std::vector<Real> vbbar(Nz);
            {
                Box vertBox(IntVect::Zero, ez * Nz);

                FArrayBox bbarFAB(vertBox, 1);
                FArrayBox TbarFAB(vertBox, 1);
                FArrayBox SbarFAB(vertBox, 1);
                FArrayBox zFAB(vertBox, 1);

                IntVect iv = vertBox.smallEnd();
                for (size_t idx = 0; idx < Nz; ++idx, ++iv[zdir]) {
                    zFAB(iv)    = vz[idx];
                    TbarFAB(iv) = vTbar[idx];
                    SbarFAB(iv) = vSbar[idx];
                }

                this->equationOfState(bbarFAB, TbarFAB, SbarFAB, zFAB);

                iv = vertBox.smallEnd();
                for (size_t idx = 0; idx < Nz; ++idx, ++iv[zdir]) {
                    vbbar[idx] = bbarFAB(iv);
                }
            }

            // Set up *bar splines
            m_bbarSplinePtr.reset(new CubicSpline);
            m_TbarSplinePtr.reset(new CubicSpline);
            m_SbarSplinePtr.reset(new CubicSpline);

            m_bbarSplinePtr->solve(vbbar, vz);
            m_TbarSplinePtr->solve(vTbar, vz);
            m_SbarSplinePtr->solve(vSbar, vz);

        } // if not / if is level 0.
    } // end if setting up splines.

    // Splines are ready. Use them...
    const Box    vertBox = Subspace::verticalDataBox(m_levGeoPtr->getDomain());
    const Box    vertBoxGrow = Box(vertBox).grow(zdir, 1);
    const size_t Nz          = vertBox.size(zdir);

    // Compute z.
    FArrayBox zFAB(vertBoxGrow, 1);
    m_levGeoPtr->fill_physCoor(zFAB, 0, zdir);

    // Allocation
    m_bbarPtr.reset(new FArrayBox(vertBoxGrow, 1));
    m_TbarPtr.reset(new FArrayBox(vertBoxGrow, 1));
    m_SbarPtr.reset(new FArrayBox(vertBoxGrow, 1));
    m_NsqPtr .reset(new FArrayBox(vertBoxGrow, 1));
    m_phi1Ptr.reset(new FArrayBox(vertBoxGrow, 1));

    if (!m_hasStrat) {
        m_bbarPtr->setVal(0.0);
        m_TbarPtr->setVal(0.0);
        m_SbarPtr->setVal(0.0);
        m_NsqPtr ->setVal(0.0);
        m_phi1Ptr->setVal(0.0);
        m_c1 = 0.0;

    } else {
        std::vector<Real> vbbar(Nz);
        std::vector<Real> vTbar(Nz);
        std::vector<Real> vSbar(Nz);

        std::vector<Real> vz(Nz);
        IntVect iv = vertBox.smallEnd();
        for (size_t idx = 0; idx < Nz; ++idx, ++iv[zdir]) {
            vz[idx] = zFAB(iv);
        }

        // Interpolate
        const int interpOrder = queryParm<int>("strat", "interpOrder", 1);
        if (interpOrder != 1 && interpOrder != 3) {
            MAYDAYERROR("strat.interpOrder = " << interpOrder << " not recognized."
                        " Must be 1 (linear) or 3 (cubic splines).");
        }

        const bool useLinearInterp = (interpOrder == 1);
        m_bbarSplinePtr->interp(vbbar, vz, useLinearInterp);
        m_TbarSplinePtr->interp(vTbar, vz, useLinearInterp);
        m_SbarSplinePtr->interp(vSbar, vz, useLinearInterp);

        // Send to FABs
        iv = vertBox.smallEnd();
        for (size_t idx = 0; idx < Nz; ++idx, ++iv[zdir]) {
            (*m_bbarPtr)(iv) = vbbar[idx];
            (*m_TbarPtr)(iv) = vTbar[idx];
            (*m_SbarPtr)(iv) = vSbar[idx];
        }

        // Compute Nsq = -db/dz using simple finite differences.
        {
            iv = vertBox.smallEnd();
            size_t idx = 0;

            (*m_NsqPtr)(iv) = -(vbbar[idx + 1] - vbbar[idx]) /
                               (vz[idx + 1] - vz[idx]);
            ++idx;
            ++iv[zdir];

            while (idx < Nz - 1) {
                (*m_NsqPtr)(iv) = -(vbbar[idx + 1] - vbbar[idx - 1]) /
                                   (vz[idx + 1] - vz[idx - 1]);
                ++idx;
                ++iv[zdir];
            }

            (*m_NsqPtr)(iv) = -(vbbar[idx] - vbbar[idx - 1]) /
                               (vz[idx] - vz[idx - 1]);

        }

        // Extrapolation to ghosts.
        for (SideIterator sit; sit.ok(); ++sit) {
            BCTools::extrap(*m_bbarPtr, vertBox, zdir, sit(), 2);
            BCTools::extrap(*m_TbarPtr, vertBox, zdir, sit(), 2);
            BCTools::extrap(*m_SbarPtr, vertBox, zdir, sit(), 2);
            BCTools::extrap(*m_NsqPtr , vertBox, zdir, sit(), 0);
        }
        // Make sure Nsq is positive everywhere inside the domain.
        static bool flowIsUnstable = false;
        flowIsUnstable |= (m_NsqPtr->min(vertBox, 0) <= 0.0);

        if (flowIsUnstable) {
            static bool userWarned = false;
            if (!userWarned) {
                if (Comm::iAmRoot()) {
                    MAYDAYWARNING(
                        "Your stratification choices caused Nsq to be non-positive "
                        "somewhere on level "
                        << m_level
                        << ". If this was intentional, then disregard this warning."
                        " Just know that the vertical structure functions will not "
                        "be available.");
                }
            }
            userWarned = true;

            m_phi1Ptr->setVal(quietNAN);
            m_c1 = quietNAN;
        } else {
            // Compute c1 and phi1.
            this->computeStructureFunctions(m_c1, *m_phi1Ptr, *m_NsqPtr);

            // Constant extrapolation to ghosts.
            // TODO: Do we need ghosts? Should we set homog Diri BCs?
            for (SideIterator sit; sit.ok(); ++sit) {
                BCTools::extrap(*m_phi1Ptr, vertBox, zdir, sit(), 0);
            }
        }
    }

    // Write the stratification data to file.
    {
        bool isPeriodic[CH_SPACEDIM];
        D_TERM(
        isPeriodic[0] = 0;,
        isPeriodic[1] = 0;,
        isPeriodic[2] = 0;)
        isPeriodic[SpaceDim - 1] =
            m_levGeoPtr->getDomain().isPeriodic(SpaceDim - 1);

        ProblemDomain vertDomain(vertBox, isPeriodic);

        std::vector<Box> vbox(1, vertBox);
        std::vector<int> vproc(1, 0);
        DisjointBoxLayout vertGrids(vbox, vproc, vertDomain);

        LevelGeometry vertLG(vertDomain,
                             m_levGeoPtr->getDomainLength(),
                             nullptr,
                             &m_levGeoPtr->getGeoSource());
        vertLG.createMetricCache(vertGrids);

        std::vector<std::string> compNames(6);
        compNames[0] = "z";
        compNames[1] = "bbar";
        compNames[2] = "Tbar";
        compNames[3] = "Sbar";
        compNames[4] = "Nsq";
        compNames[5] = "phi_1";

        LevelData<FArrayBox> stratData(vertGrids, compNames.size(), ez);
        Subspace::horizontalExtrusion(stratData, 0, zFAB, 0, 1);
        Subspace::horizontalExtrusion(stratData, 1, *m_bbarPtr, 0, 1);
        Subspace::horizontalExtrusion(stratData, 2, *m_TbarPtr, 0, 1);
        Subspace::horizontalExtrusion(stratData, 3, *m_SbarPtr, 0, 1);
        Subspace::horizontalExtrusion(stratData, 4, *m_NsqPtr, 0, 1);
        Subspace::horizontalExtrusion(stratData, 5, *m_phi1Ptr, 0, 1);

        char fn[80];
        const char* prefix = ctx->output.plotPrefix.c_str();
        sprintf(fn, "%sstratData_lev%d.hdf5", prefix, m_level);
        IO::writeHDF5(fn, stratData, vertLG, 0.0, 0.0, compNames);
    }

    if (ctx->output.verbosity >= 1 && m_hasStrat) {
        pout() << "Mode-1 wave speed = " << m_c1 << '\n';
    }
}


#define ARRAY1D(a) \
    Box(IntVect(D_DECL((1), (1), (1))), IntVect(D_DECL((a), (1), (1)))), (1)
#define ARRAY2D(a, b) \
    Box(IntVect(D_DECL((1), (1), (1))), IntVect(D_DECL((a), (b), (1)))), (1)

// -----------------------------------------------------------------------------
void
AMRNSLevel::computeStructureFunctions(Real&            a_c1,
                                      FArrayBox&       a_phi1FAB,
                                      const FArrayBox& a_NsqFAB) const
{
    BEGIN_FLOWCHART();

    // Collect data structures
    const RealVect&      dXi    = m_levGeoPtr->getDXi();
    const ProblemDomain& domain = this->getDomain();
    const Box&           domBox = domain.domainBox();
    const IntVect        ez     = BASISV(SpaceDim - 1);

    CH_assert(!domain.isPeriodic(SpaceDim - 1)); // Would need to set different phi BCs.

    if (m_level == 0) {
        // Create data holders
        char                      dflag[1]  = { 'S' };
        static const Real         abstol    = 2.0 * lapack::dlamch_(dflag);
        const int                 N         = domain.size(SpaceDim - 1);
        const Real                dZeta     = dXi[SpaceDim - 1];
        const GeoSourceInterface& geoSource = m_levGeoPtr->getGeoSource();

        const Box valid(domBox.smallEnd(SpaceDim - 1) * ez,
                        domBox.bigEnd(SpaceDim - 1) * ez);
        const Box FCvalid = surroundingNodes(valid, SpaceDim - 1);

        // Fill dXi^i/dz
        FArrayBox jacFCFAB(FCvalid, 1);
        geoSource.fill_dXidx(jacFCFAB, 0, SpaceDim - 1, dXi);

        FArrayBox jacCCFAB(valid, 1);
        geoSource.fill_dXidx(jacCCFAB, 0, SpaceDim - 1, dXi);

        IntVect loIV = valid.smallEnd();

        BaseFab<int>  IFAIL(ARRAY1D(1));      // (1)
        BaseFab<int>  IWORK(ARRAY1D(5 * N));  // (5*N)
        BaseFab<Real> AB(ARRAY2D(2, N));      // (2,N)
        BaseFab<Real> BB(ARRAY2D(1, N));      // (1,N)
        BaseFab<Real> Q(ARRAY2D(N, N));       // (N,N)
        BaseFab<Real> W(ARRAY1D(N));          // (N)
        BaseFab<Real> WORK(ARRAY1D(7 * N));   // (7*N)
        BaseFab<Real> Z(ARRAY2D(N, N));       // (N,N)

        FORT_SOLVEVERTEIGENPROBLEM(CHF_REAL(a_c1),
                                   CHF_FRA1(a_phi1FAB, 0),
                                   CHF_CONST_FRA1(a_NsqFAB, 0),
                                   CHF_CONST_FRA1(jacFCFAB, 0),
                                   CHF_CONST_FRA1(jacCCFAB, 0),
                                   CHF_CONST_INTVECT(loIV),
                                   CHF_CONST_INT(N),
                                   CHF_CONST_REAL(dZeta),
                                   CHF_CONST_REAL(abstol),
                                   CHF_FIA1(IFAIL, 0),
                                   CHF_FIA1(IWORK, 0),
                                   CHF_FRA1(AB, 0),
                                   CHF_FRA1(BB, 0),
                                   CHF_FRA1(Q, 0),
                                   CHF_FRA1(W, 0),
                                   CHF_FRA1(WORK, 0),
                                   CHF_FRA1(Z, 0));

        // Set BCs -- Homogeneous Dirichlet on phi at top and bottom.
        {
            IntVect gv, nv;

            nv = valid.smallEnd();
            gv = nv - ez;
            a_phi1FAB(gv) = -a_phi1FAB(nv);

            nv = valid.bigEnd();
            gv = valid.bigEnd() + ez;
            a_phi1FAB(gv) = -a_phi1FAB(nv);
        }

    } else {
        // On the finer levels, we can just interpolate up.
        const int            zdir       = SpaceDim - 1;
        const AMRNSLevel*    crsePtr    = this->crseNSPtr();
        const FArrayBox&     crsePhi1   = *crsePtr->m_phi1Ptr;
        const ProblemDomain& crseDomain = crsePtr->getDomain();

        IntVect refRatio = IntVect::Unit;
        refRatio[zdir] = this->getCrseRefRatio()[zdir];

        Box crseValid(IntVect::Zero, IntVect::Zero);
        crseValid.setRange(zdir,
                           crseDomain.domainBox().smallEnd(zdir),
                           crseDomain.size(zdir));

        const Box refbox(IntVect::Zero, refRatio - IntVect::Unit);

        // c1 is easy...
        a_c1 = crsePtr->m_c1;

        // Now phi1...
        // Begin with piecewise constant interpolation (order = 0).
        FORT_CFINTERP_UNMAPPEDINTERPCONSTANT(CHF_FRA(a_phi1FAB),
                                             CHF_CONST_FRA(crsePhi1),
                                             CHF_BOX(crseValid),
                                             CHF_CONST_INTVECT(refRatio),
                                             CHF_BOX(refbox));

        // Update to a linear interpolation (the correction).
        if (refRatio[zdir] > 1) {
            FArrayBox slope(crseValid, 1);
            FORT_CFINTERP_INTERPCENTRALSLOPE(CHF_FRA(slope),
                                            CHF_CONST_FRA(crsePhi1),
                                            CHF_BOX(crseValid),
                                            CHF_CONST_INT(zdir));

            FORT_CFINTERP_UNMAPPEDINTERPLINEAR(CHF_FRA(a_phi1FAB),
                                               CHF_CONST_FRA(slope),
                                               CHF_BOX(crseValid),
                                               CHF_CONST_INT(zdir),
                                               CHF_CONST_INTVECT(refRatio),
                                               CHF_BOX(refbox));
        }
    }
}

#undef ARRAY1D
#undef ARRAY2D
