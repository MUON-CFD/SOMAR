#include "AMRNSLevel.H"
#include "AMRNSLevelF_F.H"
#include "CFInterpF_F.H"
#include "Comm.H"
#include "lapack.H"
#include "MiscUtils.H"
#include "Subspace.H"

#include "LDFlubOps.H"
#include <chrono>


// ******************** Fill with time-interpolated state **********************

// -----------------------------------------------------------------------------
// Fills a data holder with the Cartesian-based velocity.
// All ghosts will be filled, even at edges and vertices.
// -----------------------------------------------------------------------------
void
AMRNSLevel::fillVelocity(LevelData<FluxBox>& a_vel,
                         const Real          a_time,
                         const bool          a_setBCs) const
{
    CH_assert(a_time >= 0.0);
    CH_assert(a_vel.ghostVect() >= IntVect::Unit || !a_setBCs);
    CH_assert(a_vel.nComp() == 1); // For now.
    CH_assert(a_vel.getBoxes().physDomain() == this->getDomain());
    debugInitLevel(a_vel);

    // Interpolate state
    if (abs(a_time - m_time) < timeEps) {
        DataIterator dit = a_vel.dataIterator();
        Box overlap;
        for (dit.reset(); dit.ok(); ++dit) {
            for (int dir = 0; dir < SpaceDim; ++dir) {
                a_vel[dit][dir].copy(m_statePtr->vel[dit][dir]);
            }
        }
    } else {
        m_parkPtr->velTimeInterp(a_vel, a_time);
    }

    // Set all BCs, if needed.
    if (a_setBCs) {
        this->setVelBC(a_vel, a_time, false);
    }
}


// -----------------------------------------------------------------------------
// Fills a data holder with the pressure. All ghosts will be filled,
// even at edges and vertices.
// -----------------------------------------------------------------------------
void
AMRNSLevel::fillPressure(LevelData<FArrayBox>& a_p,
                         const Real            a_time,
                         const bool            a_setBCs) const
{
    CH_assert(a_time >= 0.0);
    CH_assert(a_p.nComp() == 1);
    CH_assert(a_p.getBoxes().physDomain() == this->getDomain());
    debugInitLevel(a_p);

    // Interpolate state
    if (abs(a_time - m_time) < timeEps) {
        DataIterator dit = a_p.dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            a_p[dit].copy(m_statePtr->p[dit]);
        }
    } else {
        m_parkPtr->pTimeInterp(a_p, a_time);
    }

    // Set all BCs, if needed.
    if (a_setBCs) {
        this->setPressureBC(a_p, a_time, false);
    }
}


// -----------------------------------------------------------------------------
// Fills a data holder with the temperature. All ghosts will be filled,
// even at edges and vertices.
// -----------------------------------------------------------------------------
void
AMRNSLevel::fillTemperature(LevelData<FArrayBox>& a_state,
                            const Real            a_time,
                            const bool            a_setBCs,
                            const bool            a_removeBackground) const
{
    CH_assert(a_time >= 0.0);
    CH_assert(a_state.nComp() == 1);
    CH_assert(a_state.getBoxes().physDomain() == this->getDomain());
    debugInitLevel(a_state);

    // Interpolate state
    if (abs(a_time - m_time) < timeEps) {
        DataIterator dit = a_state.dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            a_state[dit].copy(m_statePtr->T[dit]);
        }
    } else {
        m_parkPtr->qTimeInterp(a_state, a_time, m_statePtr->TComp, 0, 1);
    }

    // Set all BCs, if needed.
    if (a_setBCs) {
        this->setTemperatureBC(a_state, a_time, false);
    }

    if (a_removeBackground) {
        CH_assert(a_state.ghostVect() <= IntVect::Unit);
        Subspace::addHorizontalExtrusion(a_state, 0, *m_TbarPtr, 0, 1, -1.0);
    }
}


// -----------------------------------------------------------------------------
// Fills a data holder with the salinity. All ghosts will be filled,
// even at edges and vertices.
// -----------------------------------------------------------------------------
void
AMRNSLevel::fillSalinity(LevelData<FArrayBox>& a_state,
                         const Real            a_time,
                         const bool            a_setBCs,
                         const bool            a_removeBackground) const
{
    CH_assert(a_time >= 0.0);
    CH_assert(a_state.nComp() == 1);
    CH_assert(a_state.getBoxes().physDomain() == this->getDomain());
    debugInitLevel(a_state);

    // Interpolate state
    if (abs(a_time - m_time) < timeEps) {
        DataIterator dit = a_state.dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            a_state[dit].copy(m_statePtr->S[dit]);
        }
    } else {
        m_parkPtr->qTimeInterp(a_state, a_time, m_statePtr->SComp, 0, 1);
    }

    // Set all BCs, if needed.
    if (a_setBCs) {
        this->setSalinityBC(a_state, a_time, false);
    }

    if (a_removeBackground) {
        CH_assert(a_state.ghostVect() <= IntVect::Unit);
        Subspace::addHorizontalExtrusion(a_state, 0, *m_SbarPtr, 0, 1, -1.0);
    }
}


// -----------------------------------------------------------------------------
// Fills a data holder with the scalars. All ghosts will be filled,
// even at edges and vertices.
// -----------------------------------------------------------------------------
void
AMRNSLevel::fillScalar(LevelData<FArrayBox>& a_s,
                       const Real            a_time,
                       const bool            a_setBCs) const
{
    CH_assert(a_time >= 0.0);
    CH_assert(a_s.nComp() == m_statePtr->numScalars);
    CH_assert(a_s.getBoxes().physDomain() == this->getDomain());
    debugInitLevel(a_s);

    // Interpolate state
    if (abs(a_time - m_time) < timeEps) {
        DataIterator dit = a_s.dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            a_s[dit].copy(m_statePtr->scalars[dit]);
        }
    } else {
        const int startComp = m_statePtr->scalarsInterval.begin();
        const int numComps = m_statePtr->numScalars;
        m_parkPtr->qTimeInterp(a_s, a_time, startComp, 0, numComps);
    }

    // Set all BCs, if needed.
    if (a_setBCs) {
        this->setScalarBC(a_s, a_time, false);
    }
}


// -----------------------------------------------------------------------------
// Computes advVel^i = J * dXi^i/dx^i * cartVel^i (No sum)
// a_cartVel must have exactly 1 component.
// This function can work in-place.
// -----------------------------------------------------------------------------
void
AMRNSLevel::sendToAdvectingVelocity(LevelData<FluxBox>&       a_advVel,
                                    const LevelData<FluxBox>& a_cartVel) const
{
    // Gather references
    const GeoSourceInterface& geoSrc = m_levGeoPtr->getGeoSource();
    const RealVect&           dXi    = m_levGeoPtr->getDXi();
    const DisjointBoxLayout&  grids  = a_cartVel.getBoxes();
    DataIterator              dit    = a_cartVel.dataIterator();

    // Define a_advVel if needed.
    if (!a_advVel.isDefined()) {
        a_advVel.define(grids, 1, a_cartVel.ghostVect());
    }
    CH_assert(a_advVel.getBoxes() == grids);

    // Copy a_cartVel to a_advVel if we are not working in-place.
    if (&a_advVel != &a_cartVel) {
        CH_assert(a_advVel.getBoxes() == a_cartVel.getBoxes());
        for (dit.reset(); dit.ok(); ++dit) {
            for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                a_advVel[dit][velComp].copy(a_cartVel[dit][velComp]);
            }
        }
    }

    // Convert!
    // Note that J * dXi^0/dx^0 = dx^1/dx^1 * dXi^2/dx^2, etc.
    // We use this identity here because the RHS is easier to compute.
    // It avoids computing dXi^0/dx^0 in the u^0 calculation, etc.
    // This also allows us to avoid asking for the coordinates of cell-centers.
    CH_assert(a_cartVel.nComp() == 1);
    for (dit.reset(); dit.ok(); ++dit) {
        for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
            FArrayBox& velFAB = a_advVel[dit][fcDir];
            FArrayBox  dxdXiFAB(velFAB.box(), 1);

            for (int offset = 1; offset < SpaceDim; ++offset) {
                const int mu = (fcDir + offset) % SpaceDim;
                geoSrc.fill_dxdXi(dxdXiFAB, 0, mu, dXi);
                velFAB.mult(dxdXiFAB, 0, 0, 1);
            }
        } // fcDir
    } // dit
}


// -----------------------------------------------------------------------------
// Undo the sendToAdvectingVelocity operation. Result will be in the
// Cartesian basis. This function can work in-place.
// -----------------------------------------------------------------------------
void
AMRNSLevel::sendToCartesianVelocity(LevelData<FluxBox>&       a_cartVel,
                                    const LevelData<FluxBox>& a_advVel) const
{
    // Gather references
    const GeoSourceInterface& geoSrc = m_levGeoPtr->getGeoSource();
    const RealVect&           dXi    = m_levGeoPtr->getDXi();
    const DisjointBoxLayout&  grids  = a_advVel.getBoxes();
    DataIterator              dit    = a_advVel.dataIterator();

    // Define a_cartVel if needed.
    if (!a_cartVel.isDefined()) {
        a_cartVel.define(grids, a_advVel.nComp(), a_advVel.ghostVect());
    }
    CH_assert(a_cartVel.getBoxes() == grids);

    // Copy a_advVel to a_cartVel if we are not working in-place.
    if (&a_cartVel != &a_advVel) {
        CH_assert(a_cartVel.getBoxes() == a_advVel.getBoxes());
        for (dit.reset(); dit.ok(); ++dit) {
            for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                a_cartVel[dit][velComp].copy(a_advVel[dit][velComp]);
            }
        }
    }

    // Convert, being careful that this is the EXACT inverse of sendToAdvectingVelocity.
    CH_assert(a_advVel.nComp() == 1);
    for (dit.reset(); dit.ok(); ++dit) {
        for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
            FArrayBox& velFAB = a_cartVel[dit][fcDir];
            FArrayBox  dxdXiFAB(velFAB.box(), 1);

            for (int offset = 1; offset < SpaceDim; ++offset) {
                const int mu = (fcDir + offset) % SpaceDim;
                geoSrc.fill_dxdXi(dxdXiFAB, 0, mu, dXi);
                velFAB.divide(dxdXiFAB, 0, 0, 1);
            }
        } // fcDir
    } // dit
}


// -----------------------------------------------------------------------------
void
AMRNSLevel::sendToCartesianVelocity(BoxLayoutData<FluxBox>&       a_cartVel,
                                    const BoxLayoutData<FluxBox>& a_advVel) const
{
    // Gather references
    const GeoSourceInterface& geoSrc = m_levGeoPtr->getGeoSource();
    const RealVect&           dXi    = m_levGeoPtr->getDXi();
    const BoxLayout&          layout = a_advVel.boxLayout();
    DataIterator              dit    = a_advVel.dataIterator();

    // Define a_cartVel if needed.
    if (!a_cartVel.isDefined()) {
        a_cartVel.define(layout, a_advVel.nComp());
    }
    CH_assert(a_cartVel.boxLayout() == layout);

    // Copy a_advVel to a_cartVel if we are not working in-place.
    if (&a_cartVel != &a_advVel) {
        CH_assert(a_cartVel.boxLayout() == a_advVel.boxLayout());
        for (dit.reset(); dit.ok(); ++dit) {
            for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                a_cartVel[dit][velComp].copy(a_advVel[dit][velComp]);
            }
        }
    }

    // Convert, being careful that this is the EXACT inverse of sendToAdvectingVelocity.
    CH_assert(a_advVel.nComp() == 1);
    for (dit.reset(); dit.ok(); ++dit) {
        for (int fcDir = 0; fcDir < SpaceDim; ++fcDir) {
            FArrayBox& velFAB = a_cartVel[dit][fcDir];
            FArrayBox  dxdXiFAB(velFAB.box(), 1);

            for (int offset = 1; offset < SpaceDim; ++offset) {
                const int mu = (fcDir + offset) % SpaceDim;
                geoSrc.fill_dxdXi(dxdXiFAB, 0, mu, dXi);
                velFAB.divide(dxdXiFAB, 0, 0, 1);
            }
        } // fcDir
    } // dit
}


// -----------------------------------------------------------------------------
// The equation of state.
// Fills a data holder with the buoyancy given temperature and salinity.
//
// The FABs will be at this level's index space, but may be defined over a
// patch that is not in this level's DisjointBoxLayout. Therefore, you
// can compare the FAB's boxes to this level's ProblemDomain to figure out
// where the ghosts are, but you should not loop over the grids to find
// an associated DataIndex. For any reasonable EoS, this should not be
// a problem.
//
// T and S must be the total temperature and salinity, not just the deviation
// from the background. Likewise, this function returns the total buoyancy.
//
// By default, this uses a simple linear equation of state:
//  b = -\alpha*g*T + \beta*g*S
// with \alpha = 2.0e-4 K^-1, \beta = 8.0e-4 psu^-1, and g = 9.81 m/s^2.
// -----------------------------------------------------------------------------
void
AMRNSLevel::equationOfState(FArrayBox&       a_bFAB,
                            const FArrayBox& a_TFAB,
                            const FArrayBox& a_SFAB,
                            const FArrayBox& /*a_zFAB*/) const
{
    static const Real g      = 9.81;
    static const Real alpha  = 2.0e-4;
    static const Real beta   = 8.0e-4;

    static const Real alphag = alpha * g;
    static const Real betag  = beta * g;

    const Box& region = a_bFAB.box();

    FORT_DEFAULTLINEAREOS(
        CHF_FRA1(a_bFAB, 0),
        CHF_CONST_FRA1(a_TFAB, 0),
        CHF_CONST_FRA1(a_SFAB, 0),
        CHF_BOX(region),
        CHF_CONST_REAL(alphag),
        CHF_CONST_REAL(betag));
}


// -----------------------------------------------------------------------------
// Fills a data holder with the total energy on a single level.
// Don't expect conservation on the upper levels as the grids change.
// -----------------------------------------------------------------------------
void
AMRNSLevel::computeTotalEnergy(LevelData<FArrayBox>& a_energy,
                               const Real            a_time) const
{
    CH_assert(a_energy.getBoxes().physDomain() == m_levGeoPtr->getDomain());
    CH_assert(a_energy.getBoxes().compatible(m_levGeoPtr->getBoxes()));
    CH_assert(a_energy.nComp() == 1);

    const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();

    // Start clean
    setValLevel(a_energy, 0.0);

    // Add KE
    {
        LevelData<FluxBox> cartVel(grids, 1);
        this->fillVelocity(cartVel, a_time, false);
        this->addKE(a_energy, cartVel);
    }

    // Compute total b.
    LevelData<FArrayBox> b(grids, 1);
    LevelData<FArrayBox> z(grids, 1);
    {
        LevelData<FArrayBox> T(grids, 1);
        this->fillTemperature(T, a_time, 2);

        LevelData<FArrayBox> S(grids, 1);
        this->fillSalinity(S, a_time, 2);

        for (DataIterator dit(grids); dit.ok(); ++dit) {
            m_levGeoPtr->fill_physCoor(z[dit], 0, SpaceDim - 1);
            this->equationOfState(b[dit], T[dit], S[dit], z[dit]);
        }
    }

    // Add gravitational contribution.
    if (m_hasStrat) {
        // Use APE if we are using a background stratification
        // We need to convert b to the perturbation.
        Subspace::addHorizontalExtrusion(b, 0, *m_bbarPtr, 0, 1, -1.0);
        this->addAPE(a_energy, b);

    } else {
        // Plain ol' GPE.
        this->addGPE(a_energy, b, z);
    }
}


// -----------------------------------------------------------------------------
// Adds the kinetic energy to a data holder on a single level.
// -----------------------------------------------------------------------------
void
AMRNSLevel::addKE(LevelData<FArrayBox>&     a_energy,
                  const LevelData<FluxBox>& a_cartVel) const
{
    CH_assert(a_energy.getBoxes().physDomain() == m_levGeoPtr->getDomain());
    CH_assert(a_energy.getBoxes().compatible(m_levGeoPtr->getBoxes()));
    CH_assert(a_energy.nComp() == 1);

    CH_assert(a_cartVel.getBoxes().physDomain() == m_levGeoPtr->getDomain());
    CH_assert(a_cartVel.getBoxes().compatible(m_levGeoPtr->getBoxes()));
    CH_assert(a_cartVel.nComp() == 1);

    const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();
    DataIterator dit = grids.dataIterator();

    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox&     energyFAB = a_energy[dit];
        const FluxBox& uFlub     = a_cartVel[dit];
        const Box&     valid     = grids[dit];

        FORT_ADDKE(
            CHF_FRA1(energyFAB,0),
            CHF_CONST_FRA1(uFlub[0],0),
            CHF_CONST_FRA1(uFlub[1],0),
            CHF_CONST_FRA1(uFlub[SpaceDim - 1],0),
            CHF_BOX(valid));
    }
}


// -----------------------------------------------------------------------------
// Adds the potential energy to a data holder on a single level.
// Use this when there is no background stratification.
// -----------------------------------------------------------------------------
void
AMRNSLevel::addGPE(LevelData<FArrayBox>&       a_energy,
                   const LevelData<FArrayBox>& a_btot,
                   const LevelData<FArrayBox>& a_z [[maybe_unused]]) const
{
    CH_assert(a_energy.getBoxes().physDomain() == m_levGeoPtr->getDomain());
    CH_assert(a_energy.getBoxes().compatible(m_levGeoPtr->getBoxes()));
    CH_assert(a_energy.nComp() == 1);

    CH_assert(a_btot.getBoxes().physDomain() == m_levGeoPtr->getDomain());
    CH_assert(a_btot.getBoxes().compatible(m_levGeoPtr->getBoxes()));
    CH_assert(a_btot.nComp() == 1);

    CH_assert(a_z.getBoxes().physDomain() == m_levGeoPtr->getDomain());
    CH_assert(a_z.getBoxes().compatible(m_levGeoPtr->getBoxes()));
    CH_assert(a_z.nComp() == 1);

    const IntVect            ez     = BASISV(SpaceDim - 1);
    const Box&               domBox = m_problem_domain.domainBox();
    const RealVect&          dXi    = m_levGeoPtr->getDXi();
    const DisjointBoxLayout& grids  = a_energy.getBoxes();
    DataIterator             dit    = a_energy.dataIterator();

    FArrayBox zFAB;
    {
        Box zBox = Subspace::flattenBox(domBox, ez);
        zFAB.define(zBox, 1);
        m_levGeoPtr->getGeoSource().fill_physCoor(zFAB, 0, SpaceDim - 1, dXi);
    }

    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox&       energyFAB = a_energy[dit];
        const FArrayBox& bFAB      = a_btot[dit];
        const Box&       valid     = grids[dit];

        FORT_ADDGPE(
            CHF_FRA1(energyFAB, 0),
            CHF_CONST_FRA1(bFAB, 0),
            CHF_CONST_FRA1(zFAB, 0),
            CHF_BOX(valid));
    }
}


// -----------------------------------------------------------------------------
// Adds the available potential energy to a data holder on a single level.
// Use this when Nsq != 0.0.
// -----------------------------------------------------------------------------
void
AMRNSLevel::addAPE(LevelData<FArrayBox>&       a_energy,
                   const LevelData<FArrayBox>& a_bpert) const
{
    CH_assert(a_energy.getBoxes().physDomain() == m_levGeoPtr->getDomain());
    CH_assert(a_energy.getBoxes().compatible(m_levGeoPtr->getBoxes()));
    CH_assert(a_energy.nComp() == 1);

    CH_assert(a_bpert.getBoxes().physDomain() == m_levGeoPtr->getDomain());
    CH_assert(a_bpert.getBoxes().compatible(m_levGeoPtr->getBoxes()));
    CH_assert(a_bpert.nComp() == 1);

    const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();
    DataIterator dit = grids.dataIterator();

    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox&       energyFAB = a_energy[dit];
        const FArrayBox& bpertFAB  = a_bpert[dit];
        const FArrayBox& NsqFAB    = *m_NsqPtr;
        const Box&       valid     = grids[dit];

        FORT_ADDAPE(
            CHF_FRA1(energyFAB, 0),
            CHF_CONST_FRA1(bpertFAB, 0),
            CHF_CONST_FRA1(NsqFAB, 0),
            CHF_BOX(valid));
    }
}


// -----------------------------------------------------------------------------
void
AMRNSLevel::computePhysicalPressure(
    Vector<LevelData<FArrayBox>*>& a_amrPressure,
    const Real                     a_oldTime,
    const Real                     a_newTime) const
{
    CH_assert(m_level == 0);

    // Gather references / set parameters.
    const bool   writeToPout = (s_verbosity >= 1);
    const size_t lmin        = 0;
    const size_t lmax        = this->finestNSPtr()->m_level;
    const Real   syncDt      = a_newTime - a_oldTime;

    // Begin intentation...
    if (writeToPout) {
        pout() << "AMR pressure computation on levels " << lmin << " to "
               << lmax << ":" << Format::indent() << Format::pushFlags << endl;
    }

    // Allocation
    Vector<LevelData<FluxBox>*>   amrRhs0;
    Vector<LevelData<FArrayBox>*> amrPhi0;
    Vector<LevelData<FArrayBox>*> amrRes;
    Vector<LevelData<FArrayBox>*> amrE;
    this->allocateAndDefine(amrRhs0, 1, IntVect::Zero, lmin, lmax);
    this->allocateAndDefine(amrPhi0, 1, IntVect::Unit, lmin, lmax);
    this->allocateAndDefine(amrRes, 1, IntVect::Zero, lmin, lmax);
    this->allocateAndDefine(amrE, 1, IntVect::Unit, lmin, lmax);

    // Compute amrRhs0 = (NS rhs w/out pressure gradient) - du/dt.
    for (size_t lev = lmin; lev <= lmax; ++lev) {
        const AMRNSLevel*        levPtr   = this->getLevel(lev);
        const DisjointBoxLayout& grids    = levPtr->getBoxes();
        const Real               refluxDt = 0.0;

        LevelData<FluxBox>& rhs0 = *amrRhs0[lev];
        LevelData<FluxBox> FI(grids, 1);

        LevelData<FluxBox> vel(grids, m_velPtr->nComp(), IntVect::Unit);
        LevelData<FArrayBox> p(grids, m_pPtr->nComp(), IntVect::Unit);
        LevelData<FArrayBox> q(grids, m_qPtr->nComp(), IntVect::Unit);

        levPtr->fillVelocity(vel, a_newTime);
        levPtr->fillPressure(p, a_newTime);
        {
            LevelData<FArrayBox> qComp;
            aliasLevelData(qComp, &q, m_statePtr->scalarsInterval);
            levPtr->fillScalar(qComp, a_newTime);

            aliasLevelData(qComp, &q, m_statePtr->TInterval);
            levPtr->fillTemperature(qComp, a_newTime);

            aliasLevelData(qComp, &q, m_statePtr->SInterval);
            levPtr->fillSalinity(qComp, a_newTime);
        }

        // rhs = kvelE + kvelI
        {
            LevelData<FArrayBox> kq(grids, q.nComp());

            AMRNSLevel& castLev = const_cast<AMRNSLevel&>(*levPtr);
            castLev.setExplicitRHS(rhs0, kq, vel, p, q, a_newTime, refluxDt);
            castLev.setImplicitRHS(FI, kq, vel, p, q, a_newTime, refluxDt);

            for (DataIterator dit(grids); dit.ok(); ++dit) {
                rhs0[dit].plus(FI[dit], grids[dit], 0, 0);
            }
        }

        // rhs -= du/dt
        {
            for (DataIterator dit(grids); dit.ok(); ++dit) {
                for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                    rhs0[dit][velComp].plus(vel[dit][velComp], -1.0 / syncDt);
                }
            }

            levPtr->fillVelocity(vel, a_oldTime);
            for (DataIterator dit(grids); dit.ok(); ++dit) {
                for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                    rhs0[dit][velComp].plus(vel[dit][velComp], 1.0 / syncDt);
                }
            }
        }

        // Convert to mapped coordinates. (aka prep for divergence)
        levPtr->sendToAdvectingVelocity(rhs0, rhs0);
    }
    // this->averageDown(amrRhs0, lmin, true);

    // Compute p0 = initial pressure guess with BCs set.
    for (size_t lev = lmin; lev <= lmax; ++lev) {
        const AMRNSLevel*         levPtr = this->getLevel(lev);
        LevelData<FArrayBox>&     phi0   = *amrPhi0[lev];
        const LevelData<FluxBox>& rhs0   = *amrRhs0[lev];

        struct PressNeumBCs : public BCTools::BCFunction {
            PressNeumBCs(const LevelData<FluxBox>& a_bc)
            : m_bc(a_bc)
            {}

            virtual void
            operator()(FArrayBox&            a_alpha,
                       FArrayBox&            a_beta,
                       FArrayBox&            a_bcFAB,
                       const FArrayBox&      /* a_stateFAB */,
                       const FArrayBox&      /* a_xFAB */,
                       const DataIndex&      a_di,
                       const int             a_bdryDir,
                       const Side::LoHiSide& a_side,
                       const Real            /* a_time */,
                       const bool            a_homogBCs) const
            {
                a_alpha.setVal(0.0);
                a_beta.setVal(1.0);
                // Technically wrong, but we won't use the transverse ghosts.
                if (!a_homogBCs) {
                    a_bcFAB.copy(m_bc[a_di][a_bdryDir]);
                    if (a_side == Side::Lo) {
                        a_bcFAB *= -1.0;
                    }
                }
            }
        private:
            const LevelData<FluxBox>& m_bc;
        };

        levPtr->fillPressure(phi0, a_newTime);
        BCTools::applyBC(phi0,
                         false,
                         PressNeumBCs(rhs0),
                         a_newTime,
                         levPtr->getLevGeo(),
                         levPtr->m_allPhysBdryIter);
    }

    // Compute res = D[rhs0 - G[phi0]].
    // Do this in two steps to help refluxing.
    //  1. rhs0 = rhs0 - G[phi0]
    for (size_t lev = lmin; lev <= lmax; ++lev) {
        const AMRNSLevel*        levPtr = this->getLevel(lev);
        const DisjointBoxLayout& grids  = levPtr->getBoxes();

        const auto& fd = *(levPtr->m_finiteDiffPtr);

        using LDFlubTraits = Elliptic::StateTraits<LevelData<FluxBox>>;
        const Elliptic::StateOps<LevelData<FluxBox>, LDFlubTraits> ops;

        LevelData<FluxBox>&         rhs0 = *(amrRhs0[lev]);
        const LevelData<FArrayBox>& phi0 = *(amrPhi0[lev]);
        LevelData<FluxBox>          Gphi0(grids, 1);

        if (lev > 0) {
            fd.levelGradientMAC(Gphi0, phi0);
        } else {
            fd.compGradientMAC(Gphi0, phi0);
        }
        ops.incr(rhs0, Gphi0, -1.0);
    }
    // this->averageDown(amrRhs0, lmin, true);

    //  2. res = D[rhs0].
    for (size_t lev = lmin; lev <= lmax; ++lev) {
        const AMRNSLevel* levPtr = this->getLevel(lev);
        const auto&       fd     = *(levPtr->m_finiteDiffPtr);

        using LDFlubTraits = Elliptic::StateTraits<LevelData<FluxBox>>;
        const Elliptic::StateOps<LevelData<FluxBox>, LDFlubTraits> ops;

        auto&       res  = *(amrRes[lev]);
        const auto& rhs0 = *(amrRhs0[lev]);

        if (lev < lmax) {
            const auto& fineRhs0 = *(amrRhs0[lev + 1]);
            fd.compDivergenceMAC(res, rhs0, fineRhs0);
        } else {
            fd.levelDivergenceMAC(res, rhs0);
        }
    }

    // ----


    // Check eq. constistency.
    this->averageDown(amrRes, lmin, false);
    {
        const AMRNSLevel* levPtr  = this->getLevel(lmin);
        auto              op      = levPtr->m_projOpPtr;
        auto&             e       = *amrE[lmin];
        auto&             gradPhi = *amrRhs0[lmin];   // Use as temp space
        auto&             levGeo  = *levPtr->m_levGeoPtr;

        pout() << Format::scientific;
        op->setToZero(e);

        // const bool consistent = op->levelEquationIsConsistent(
        //     phi, nullptr, rhs, newTime, true, true);
        // if (!consistent) {
        //     pout() << "level " << lmin << " equation is NOT consistent.";
        // }

        op->compGradient(gradPhi, e, nullptr, a_newTime, true, true);
        const Real bdrySum = Integral::bdrySum(gradPhi, levGeo);
        const Real rhsSum  = Integral::sum(e, levGeo, false);
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

    // Solve L[e] = Div[res].
    {
        const bool setPhitoZero = true;
        const bool useHomogBCs  = true;

        auto start = std::chrono::high_resolution_clock::now();
        m_amrProjSolverPtr->solve(amrE,
                                  makeConstVector(amrRes),
                                  a_newTime,
                                  useHomogBCs,
                                  setPhitoZero);
        auto finish = std::chrono::high_resolution_clock::now();

        // Write solve time.
        if (writeToPout) {
            std::chrono::duration<double> elapsed = finish - start;
            pout() << "Solve time: " << elapsed.count() << " s" << endl;
        }
    }

    // Update initial guess.
    for (size_t lev = lmin; lev <= lmax; ++lev) {
        const AMRNSLevel*        levPtr   = this->getLevel(lev);
        const DisjointBoxLayout& grids    = levPtr->getBoxes();

        LevelData<FArrayBox>&       phi  = *a_amrPressure[lev];
        const LevelData<FArrayBox>& phi0 = *amrPhi0[lev];
        const LevelData<FArrayBox>& e    = *amrE[lev];

        for (DataIterator dit(grids); dit.ok(); ++dit) {
            phi[dit].copy(phi0[dit]);
            phi[dit].plus(e[dit]);
        }
    }

    // Restore indentation state.
    if (writeToPout) {
        pout() << Format::popFlags << Format::unindent << std::flush;
    }

    // Free memory
    this->deallocate(amrE);
    this->deallocate(amrRes);
    this->deallocate(amrPhi0);
    this->deallocate(amrRhs0);
}
