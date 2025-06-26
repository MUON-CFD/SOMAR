#include "AMRNSLevel.H"
#include "ScalarBC.H"
#include "Subspace.H"
#include "VelBC.H"


// **************************** Set physical BCs *******************************

// -----------------------------------------------------------------------------
BCTools::BCFunction*
AMRNSLevel::createVelPhysBC(const int a_dir, const Side::LoHiSide& a_side) const
{
    const auto* ctx   = ProblemContext::getInstance();
    const auto  btype = ctx->rhs.velBCType[a_dir][int(a_side)];

    if (m_problem_domain.isPeriodic(a_dir)) {
        return nullptr;

    } else if (btype == RHSParameters::VelBCType::NO_SLIP) {
        return new VelBC::NoSlip;

    } else if (btype == RHSParameters::VelBCType::FREE_SLIP) {
        return new VelBC::FreeSlip;

    } else if (btype == RHSParameters::VelBCType::TIDAL) {
        return new VelBC::Tidal;

    } else if (btype == RHSParameters::VelBCType::OUTFLOW) {
        return new VelBC::Outflow;

    } else if (btype == RHSParameters::VelBCType::CUSTOM) {
        MAYDAYERROR(
            "For CUSTOM BCs, you must override AMRNSLevel::createVelBC.");
        return nullptr;

    } else {
        MAYDAYERROR("velBCType["
                    << a_dir << "]["
                    << ((a_side == Side::Lo) ? "Side::Lo" : "Side::Hi")
                    << "] = " << btype << " has no default functionality.");
        return nullptr;
    }
}


// -----------------------------------------------------------------------------
BCTools::BCFunction*
AMRNSLevel::createPressurePhysBC(const int             /* a_dir */,
                                 const Side::LoHiSide& /* a_side */) const
{
    return new BCTools::HomogNeumBC;
}


// -----------------------------------------------------------------------------
BCTools::BCFunction*
AMRNSLevel::createTemperaturePhysBC(const int             a_dir,
                                    const Side::LoHiSide& a_side) const
{
    const auto* ctx   = ProblemContext::getInstance();
    const auto  btype = ctx->rhs.tempBCType[a_dir][int(a_side)];

    if (m_problem_domain.isPeriodic(a_dir)) {
        return nullptr;

    } else if (btype == RHSParameters::TempBCType::ZERO_NEUM) {
        return new ScalarBC::NormalGradient;

    } else if (btype == RHSParameters::TempBCType::ZERO_NEUM_ON_PERT) {
        // Heat flux = background value.
        // Homogeneous Neumann BCs on Tpert.
        // Neumann BCs on T: Grad[T].n = Grad[Tbar].n.
        CH_assert(m_TbarPtr);
        return new ScalarBC::NormalGradient(*m_TbarPtr, *m_levGeoPtr);

    } else if (btype == RHSParameters::TempBCType::TIDAL) {
        // Tidal inflow/outflow.
        UNDEFINED_FUNCTION();  // I'm still working on this.
        return new ScalarBC::Tidal;

    } else if (btype == RHSParameters::TempBCType::OUTFLOW) {
        // Homogeneous Neumann BCs on Tpert.
        // Neumann BCs on T: Grad[T].n = Grad[Tbar].n.
        CH_assert(m_TbarPtr);
        return new ScalarBC::NormalGradient(*m_TbarPtr, *m_levGeoPtr);

    } else if (btype == RHSParameters::TempBCType::CUSTOM) {
        MAYDAYERROR(
            "For CUSTOM BCs, you must override "
            "AMRNSLevel::createTemperatureBC.");
        return nullptr;

    } else if (btype == RHSParameters::TempBCType::ZERO_DIRI_ON_PERT) {
        // Temp = background value.
        // Homogeneous Dirichlet BCs on Tpert.
        // Dirichlet BCs on T: T = Tbar.
        CH_assert(m_TbarPtr);
        return new ScalarBC::HomogDiriOnPert(m_TbarPtr, *m_levGeoPtr);

    } else {
        MAYDAYERROR("tempBCType["
                    << a_dir << "]["
                    << ((a_side == Side::Lo) ? "Side::Lo" : "Side::Hi")
                    << "] = " << btype << " has no default functionality.");
        return nullptr;
    }
}


// -----------------------------------------------------------------------------
BCTools::BCFunction*
AMRNSLevel::createSalinityPhysBC(const int             a_dir,
                                 const Side::LoHiSide& a_side) const
{
    const auto* ctx   = ProblemContext::getInstance();
    const auto  btype = ctx->rhs.salinityBCType[a_dir][int(a_side)];

    if (m_problem_domain.isPeriodic(a_dir)) {
        return nullptr;

    } else if (btype == RHSParameters::SalinityBCType::ZERO_NEUM) {
        // Isohalines normal to wall.
        // Homogeneous Neumann BCs on S: Grad[S].n = 0.
        return new ScalarBC::NormalGradient;

    } else if (btype == RHSParameters::SalinityBCType::ZERO_NEUM_ON_PERT) {
        // Salt flux = background value.
        // Neumann BCs on S: Grad[S].n = Grad[Sbar].n.
        CH_assert(m_SbarPtr);
        return new ScalarBC::NormalGradient(*m_SbarPtr, *m_levGeoPtr);

    } else if (btype == RHSParameters::SalinityBCType::TIDAL) {
        // Tidal inflow/outflow.
        UNDEFINED_FUNCTION();  // I'm still working on this.
        return new ScalarBC::Tidal;

    } else if (btype == RHSParameters::SalinityBCType::OUTFLOW) {
        // Homogeneous Neumann BCs on Spert.
        // Neumann BCs on T: Grad[S].n = Grad[Sbar].n.
        CH_assert(m_SbarPtr);
        return new ScalarBC::NormalGradient(*m_SbarPtr, *m_levGeoPtr);

    } else if (btype == RHSParameters::SalinityBCType::ZERO_DIRI_ON_PERT) {
        // Salinity = background value.
        // Homogeneous Dirichlet BCs on Spert.
        // Dirichlet BCs on S: S = Sbar.
        CH_assert(m_SbarPtr);
        return new ScalarBC::HomogDiriOnPert(m_SbarPtr, *m_levGeoPtr);

    } else if (btype == RHSParameters::SalinityBCType::CUSTOM) {
        MAYDAYERROR(
            "For CUSTOM BCs, you must override AMRNSLevel::createSalinityBC.");
        return nullptr;

    } else {
        MAYDAYERROR("salinityBCType["
                    << a_dir << "]["
                    << ((a_side == Side::Lo) ? "Side::Lo" : "Side::Hi")
                    << "] = " << btype << " has no default functionality.");
        return nullptr;
    }
}


// -----------------------------------------------------------------------------
BCTools::BCFunction*
AMRNSLevel::createScalarsPhysBC(const int             /* a_dir */,
                                const Side::LoHiSide& /* a_side */) const
{
    if (this->numScalars() > 0) {
        const ProblemDomain& domain = this->getDomain();

        bool isTotallyPeriodic = true;
        D_TERM(isTotallyPeriodic &= domain.isPeriodic(0);
            , isTotallyPeriodic &= domain.isPeriodic(1);
            , isTotallyPeriodic &= domain.isPeriodic(2);)

        if (!isTotallyPeriodic) {
            MAYDAYERROR(
                "You must override AMRNSLevel::createScalarsBC in your physics "
                "class if you plan to use a user-defined scalar in a non-periodic "
                "domain.");
        }
    }

    return nullptr;
}


// -----------------------------------------------------------------------------
void
AMRNSLevel::setVelPhysBC(LevelData<FluxBox>& a_vel,
                         const Real          a_time,
                         const bool          a_homogBCs,
                         const SideArray&    a_doSides) const
{
    if (a_vel.getBoxes() == this->getBoxes()) {
        for (int bdryDir = 0; bdryDir < SpaceDim; ++bdryDir) {
            for (SideIterator sit; sit.ok(); ++sit) {
                const int iside = int(sit());
                if (a_doSides[bdryDir][iside] == 0) continue;

                for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                    BCTools::applyBC(a_vel,
                                    velComp,
                                    a_homogBCs,
                                    *this->velPhysBC(),
                                    a_time,
                                    *m_levGeoPtr,
                                    m_singlePhysBdryIter[bdryDir][iside]);
                }
            }
        }
    } else {
        for (int bdryDir = 0; bdryDir < SpaceDim; ++bdryDir) {
            for (SideIterator sit; sit.ok(); ++sit) {
                const int iside = int(sit());
                if (a_doSides[bdryDir][iside] == 0) continue;

                PhysBdryIter it(a_vel.getBoxes(), a_doSides);
                for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                    BCTools::applyBC(a_vel,
                                     velComp,
                                     a_homogBCs,
                                     *this->velPhysBC(),
                                     a_time,
                                     *m_levGeoPtr, // This does not need to be compatible.
                                     it);
                }
            }
        }
    }
}


// -----------------------------------------------------------------------------
void
AMRNSLevel::setPressurePhysBC(LevelData<FArrayBox>& a_p,
                              const Real            a_time,
                              const bool            a_homogBCs,
                              const SideArray&      a_doSides) const
{
    for (int bdryDir = 0; bdryDir < SpaceDim; ++bdryDir) {
        for (SideIterator sit; sit.ok(); ++sit) {
            const int iside = int(sit());
            if (a_doSides[bdryDir][iside] == 0) continue;
            BCTools::applyBC(a_p,
                             a_homogBCs,
                             *this->pressurePhysBC(),
                             a_time,
                             *m_levGeoPtr,
                             m_singlePhysBdryIter[bdryDir][iside]);
        }
    }
}


// -----------------------------------------------------------------------------
void
AMRNSLevel::setTemperaturePhysBC(LevelData<FArrayBox>& a_T,
                                 const Real            a_time,
                                 const bool            a_homogBCs,
                                 const SideArray&      a_doSides) const
{
    if (a_T.getBoxes() == this->getBoxes()) {
        for (int bdryDir = 0; bdryDir < SpaceDim; ++bdryDir) {
            for (SideIterator sit; sit.ok(); ++sit) {
                const int iside = int(sit());
                if (a_doSides[bdryDir][iside] == 0) continue;

                BCTools::applyBC(a_T,
                                 a_homogBCs,
                                 *this->temperaturePhysBC(),
                                 a_time,
                                 *m_levGeoPtr,
                                 m_singlePhysBdryIter[bdryDir][iside]);
            }
        }
    } else {
        for (int bdryDir = 0; bdryDir < SpaceDim; ++bdryDir) {
            for (SideIterator sit; sit.ok(); ++sit) {
                const int iside = int(sit());
                if (a_doSides[bdryDir][iside] == 0) continue;

                PhysBdryIter it(a_T.getBoxes(), a_doSides);
                BCTools::applyBC(a_T,
                                 a_homogBCs,
                                 *this->temperaturePhysBC(),
                                 a_time,
                                 *m_levGeoPtr,
                                 it);
            }
        }
    }
}


// -----------------------------------------------------------------------------
// Set physical BCs on salinity.
// -----------------------------------------------------------------------------
void
AMRNSLevel::setSalinityPhysBC(LevelData<FArrayBox>& a_S,
                              const Real            a_time,
                              const bool            a_homogBCs,
                              const SideArray&      a_doSides) const
{
    if (a_S.getBoxes() == this->getBoxes()) {
        for (int bdryDir = 0; bdryDir < SpaceDim; ++bdryDir) {
            for (SideIterator sit; sit.ok(); ++sit) {
                const int iside = int(sit());
                if (a_doSides[bdryDir][iside] == 0) continue;

                BCTools::applyBC(a_S,
                                 a_homogBCs,
                                 *this->salinityPhysBC(),
                                 a_time,
                                 *m_levGeoPtr,
                                 m_singlePhysBdryIter[bdryDir][iside]);
            }
        }
    } else {
        for (int bdryDir = 0; bdryDir < SpaceDim; ++bdryDir) {
            for (SideIterator sit; sit.ok(); ++sit) {
                const int iside = int(sit());
                if (a_doSides[bdryDir][iside] == 0) continue;

                PhysBdryIter it(a_S.getBoxes(), a_doSides);
                BCTools::applyBC(a_S,
                                 a_homogBCs,
                                 *this->salinityPhysBC(),
                                 a_time,
                                 *m_levGeoPtr,
                                 it);
            }
        }
    }
}


// -----------------------------------------------------------------------------
// Set physical BCs on all scalars.
// Default: Throws an error if domain is not totally periodic.
// -----------------------------------------------------------------------------
void
AMRNSLevel::setScalarsPhysBC(LevelData<FArrayBox>& a_scalars,
                             const Real            a_time,
                             const bool            a_homogBCs,
                             const SideArray&      a_doSides) const
{
    for (int bdryDir = 0; bdryDir < SpaceDim; ++bdryDir) {
        if (m_problem_domain.isPeriodic(bdryDir)) continue;

        for (int iside = 0; iside < 2; ++iside) {
            if (a_doSides[bdryDir][iside] == 0) continue;

            PhysBdryIter& it = m_singlePhysBdryIter[bdryDir][iside];
            BCTools::applyBC(a_scalars,
                             a_homogBCs,
                             *this->scalarsPhysBC(),
                             a_time,
                             *m_levGeoPtr,
                             it);

        } // iside
    } // bdryDir
}


// ****************************** Set all BCs *********************************

// -----------------------------------------------------------------------------
// Resets all BCs on all state variables. This will involve several
// exchanges and CFI interpolations. Use sparingly.
// NOTE: Unlike the other setBC functions, this assumes vel is in the
//  Cartesian basis!
// -----------------------------------------------------------------------------
void
AMRNSLevel::setBC(State& a_state, const Real a_time) const
{
    CH_assert(a_state.grids.compatible(this->getBoxes()));
    CH_assert(a_state.grids == this->getBoxes());

    this->setVelBC(a_state.vel, a_time, false);
    this->setPressureBC(a_state.p, a_time, false);
    this->setTemperatureBC(a_state.T, a_time, false);
    this->setSalinityBC(a_state.S, a_time, false);
    if (this->numScalars() > 0) {
        this->setScalarBC(a_state.scalars, a_time, false);
    }
}


// -----------------------------------------------------------------------------
// Set all BCs on velocity.
// This assumes a_vel is in the curvilinear (mapped) basis.
// -----------------------------------------------------------------------------
void
AMRNSLevel::setVelBC(LevelData<FluxBox>&       a_vel,
                     const Real                a_time,
                     const bool                a_homogBCs,
                     const LevelData<FluxBox>* a_crseVelPtr) const
{
    const DisjointBoxLayout& grids = this->getBoxes();
    DataIterator             dit   = grids.dataIterator();

    const bool exCopiersAreCached =
        a_vel.ghostVect() == m_statePtr->vel.ghostVect() &&
        a_vel.getBoxes() == m_statePtr->vel.getBoxes();

    const bool hasGhosts = (a_vel.ghostVect().sum() > 0);

    debugCheckValidFaceOverlap(a_vel);

    // Without this, ghosts at some boundary vertices are not filled.
    // The stencil for IB needs multiple ghost layers. We need either this
    // extrap or a biased IB stencil. I'm going with this extrap for now.
    BCTools::extrapAllGhosts(a_vel, 2);

    // Set CFI BCs. This MUST happen first!
    if (m_level > 0) {
        // Gather needed structures
        AMRNSLevel*              crsePtr   = this->crseNSPtr();
        const DisjointBoxLayout& crseGrids = crsePtr->getBoxes();

        // Time interpolation
        LevelData<FluxBox> crseVel;
        if (a_homogBCs) {
            crseVel.define(crseGrids, 1, IntVect::Unit);
            setValLevel(crseVel, 0.0);
        } else {
            if (a_crseVelPtr) {
                // Alias the user's crse vel.
                CH_assert(a_crseVelPtr->nComp() == a_vel.nComp());
                aliasLevelData(crseVel,
                               const_cast<LevelData<FluxBox>*>(a_crseVelPtr),
                               a_crseVelPtr->interval());

            } else if (abs(m_time - crsePtr->m_time) < timeEps) {
                // Alias the coarse level's vel.
                CH_assert(crsePtr->m_velPtr->nComp() == a_vel.nComp());
                aliasLevelData(crseVel,
                               const_cast<LevelData<FluxBox>*>(crsePtr->m_velPtr),
                               crsePtr->m_velPtr->interval());

            } else {
                // Interpolate the coarse level's vel.
                crseVel.define(crseGrids, 1, IntVect::Unit);
                crsePtr->m_parkPtr->velTimeInterp(crseVel, a_time);
            }

            // The interpolation procedure requires an advecting velocity.
            crsePtr->setVelBC(crseVel, a_time, a_homogBCs);
            crsePtr->sendToAdvectingVelocity(crseVel, crseVel);
        }  // end if !a_homogBCs

        // Space interpolation.
        this->sendToAdvectingVelocity(a_vel, a_vel);
        debugCheckValidFaceOverlap(a_vel);
        debugCheckValidFaceOverlap(crseVel);
        m_cfInterpPtr->interpGhostsAtCFI(a_vel, crseVel, a_homogBCs);
        this->sendToCartesianVelocity(a_vel, a_vel);

        // Restore the user's crseVel, if necessary.
        if (!a_homogBCs && a_crseVelPtr) {
            crsePtr->sendToCartesianVelocity(crseVel, crseVel);
        }
    }

    // Set physical BCs directly at faces and at ghosts.
    this->setVelPhysBC(a_vel, a_time, a_homogBCs);

    // Do exchanges.
    if (hasGhosts) {
        if (exCopiersAreCached) {
            a_vel.exchangeBegin(m_statePtr->velExCopier);
            a_vel.exchangeEnd();

            a_vel.exchangeBegin(m_statePtr->velExCornerCopier1);
            a_vel.exchangeEnd();

            a_vel.exchangeBegin(m_statePtr->velExCornerCopier2);
            a_vel.exchangeEnd();

        } else {
            a_vel.exchange();
        }
    }

    // Fill ghosts at corners of domain. This must happen last!
    BCTools::extrapDomainCorners(a_vel, 2);

#ifndef NDEBUG
    for (dit.reset(); dit.ok(); ++dit) {
        for (int velComp = 0; velComp < SpaceDim; ++velComp) {
            checkForNAN(a_vel[dit][velComp], a_vel[dit][velComp].box());
        }
    }
#endif
}


// -----------------------------------------------------------------------------
// Set all BCs on the pressure.
// -----------------------------------------------------------------------------
void
AMRNSLevel::setPressureBC(LevelData<FArrayBox>&       a_p,
                          const Real                  a_time,
                          const bool                  a_homogBCs,
                          const LevelData<FArrayBox>* a_crsePPtr) const
{
    const bool exCopiersAreCached =
        a_p.ghostVect() == m_statePtr->p.ghostVect() &&
        a_p.getBoxes() == m_statePtr->p.getBoxes();

    // Filling all of the ghosts properly in mapped coordinates is a nightmare.
    // Neum BCs need all nearby cells to be filled, including corner ghosts,
    // while the exchange of corner ghosts requires the physical BCs to be set.
    // It's a chicken and egg problem that requires iteration. Rather than do
    // that, I'm gonna cheat. Most of these extrapolated ghosts will be
    // properly overwritten. Some will stick around -- namely, the edge and
    // vertex ghosts that live at the intersection of the CFI and physical
    // boundary.
    BCTools::extrapAllGhosts(a_p, 2);

    // Set CFI BCs.
    // Do this first since it shouldn't depend on this level's ghosts.
    if (m_level > 0) {
        // Gather coarse and fine data structures.
        AMRNSLevel*              crsePtr   = this->crseNSPtr();
        const DisjointBoxLayout& crseGrids = crsePtr->getBoxes();

        if (a_homogBCs) {
            m_cfInterpPtr->homogInterpAtCFI(a_p);
        } else {
            // Gather the coarse level data.
            LevelData<FArrayBox> crseP;
            if (a_crsePPtr) {
                // Alias the user's crse p.
                CH_assert(a_crsePPtr->nComp() == a_p.nComp());
                aliasLevelData(crseP,
                               const_cast<LevelData<FArrayBox>*>(a_crsePPtr),
                               a_crsePPtr->interval());

            } else if (abs(m_time - crsePtr->m_time) < timeEps) {
                // Alias the coarse level's p.
                CH_assert(crsePtr->m_pPtr->nComp() == a_p.nComp());
                aliasLevelData(crseP,
                               const_cast<LevelData<FArrayBox>*>(crsePtr->m_pPtr),
                               crsePtr->m_statePtr->pInterval);

            } else {
                // Interpolate the coarse level's p.
                crseP.define(crseGrids, 1);
                crsePtr->m_parkPtr->pTimeInterp(crseP, a_time, m_statePtr->pComp, 0, 1);
            }

            m_cfInterpPtr->interpAtCFI(a_p, crseP);
        }
    }  // m_level > 0

    // Do exchange #1.
    // This must be done before we potentially call BCTools::neum.
    if (exCopiersAreCached) {
        a_p.exchange(m_statePtr->pExCopier);
    } else {
        Copier cp;
        State::pDefineExCopier(cp, a_p);
        a_p.exchange(cp);
    }

    // Set physical BCs.
    // This is where we would potentially call BCTools::neum.
    this->setPressurePhysBC(a_p, a_time, a_homogBCs);

    // Do exchange #2.
    // This must be done after all other ghosts are filled.
    if (exCopiersAreCached) {
        a_p.exchange(m_statePtr->pExCornerCopier);
    } else {
        CornerCopier ccp;
        State::pDefineExCornerCopier(ccp, a_p);
        a_p.exchange(ccp);
    }

    // Fill ghosts at corners of domain. This must happen last!
    BCTools::extrapDomainCorners(a_p, 2);
}


// -----------------------------------------------------------------------------
// Set all BCs on temperature.
// -----------------------------------------------------------------------------
void
AMRNSLevel::setTemperatureBC(LevelData<FArrayBox>&       a_T,
                             const Real                  a_time,
                             const bool                  a_homogBCs,
                             const LevelData<FArrayBox>* a_crseTPtr) const
{
    const bool exCopiersAreCached =
        a_T.ghostVect() == m_statePtr->T.ghostVect() &&
        a_T.getBoxes() == m_statePtr->T.getBoxes();

    const bool hasGhosts = (a_T.ghostVect().sum() > 0);

    // Diffusive terms need edges too.
    // Filling all of the ghosts properly in mapped coordinates is a nightmare.
    // Neum BCs need all nearby cells to be filled, including corner ghosts,
    // while the exchange of corner ghosts requires the physical BCs to be set.
    // It's a chicken and egg problem that requires iteration. Rather than do
    // that, I'm gonna cheat. Most of these extrapolated ghosts will be
    // properly overwritten. Some will stick around -- namely, the edge and
    // vertex ghosts that live at the intersection of the CFI and physical
    // boundary.
    Subspace::addHorizontalExtrusion(a_T, 0, *m_TbarPtr, 0, 1, -1.0);
    BCTools::extrapAllGhosts(a_T, 2);
    Subspace::addHorizontalExtrusion(a_T, 0, *m_TbarPtr, 0, 1, 1.0);

    // Set CFI BCs.
    // Do this first since it shouldn't depend on this level's ghosts.
    if (m_level > 0) {
        // Gather coarse and fine data structures.
        AMRNSLevel*              crsePtr   = this->crseNSPtr();
        const DisjointBoxLayout& crseGrids = crsePtr->getBoxes();

        if (a_homogBCs) {
            m_cfInterpPtr->homogInterpAtCFI(a_T);
        } else {
            // Gather the coarse level data.
            LevelData<FArrayBox> crseT;
            if (a_crseTPtr) {
                // Alias the user's crse T.
                CH_assert(a_crseTPtr->nComp() == a_T.nComp());
                aliasLevelData(crseT,
                               const_cast<LevelData<FArrayBox>*>(a_crseTPtr),
                               a_crseTPtr->interval());

            } else if (abs(m_time - crsePtr->m_time) < timeEps) {
                // Alias the coarse level's T.
                aliasLevelData(crseT,
                               const_cast<LevelData<FArrayBox>*>(crsePtr->m_qPtr),
                               crsePtr->m_statePtr->TInterval);

            } else {
                // Interpolate the coarse level's T.
                crseT.define(crseGrids, 1);
                crsePtr->m_parkPtr->qTimeInterp(crseT, a_time, m_statePtr->TComp, 0, 1);
            }

            // Interpolate to fine CFI ghosts.
            m_cfInterpPtr->interpAtCFI(a_T, crseT);
        }
    }

    // Do exchange #1 on perturbation.
    // This must be done before we potentially apply Neumann BCs.
    Subspace::addHorizontalExtrusion(a_T, 0, *m_TbarPtr, 0, 1, -1.0);
    if (hasGhosts) {
        if (exCopiersAreCached) {
            a_T.exchange(m_statePtr->qExCopier);
        } else {
            Copier cp;
            State::qDefineExCopier(cp, a_T);
            a_T.exchange(cp);
        }
    }
    Subspace::addHorizontalExtrusion(a_T, 0, *m_TbarPtr, 0, 1, 1.0);

    // Set physical BCs.
    // This is where we would potentially apply Neumann BCs.
    this->setTemperaturePhysBC(a_T, a_time, a_homogBCs);

    // Do exchange #2 on perturbation.
    // This must be done after all other ghosts are filled.
    Subspace::addHorizontalExtrusion(a_T, 0, *m_TbarPtr, 0, 1, -1.0);
    if (hasGhosts) {
        if (exCopiersAreCached) {
            a_T.exchange(m_statePtr->qExCornerCopier);
        } else {
            CornerCopier ccp;
            State::qDefineExCornerCopier(ccp, a_T);
            a_T.exchange(ccp);
        }
    }

    // Fill ghosts at corners of domain. This must happen last!
    BCTools::extrapDomainCorners(a_T, 2);
    Subspace::addHorizontalExtrusion(a_T, 0, *m_TbarPtr, 0, 1, 1.0);
}


// -----------------------------------------------------------------------------
// Set all BCs on salinity.
// -----------------------------------------------------------------------------
void
AMRNSLevel::setSalinityBC(LevelData<FArrayBox>&       a_S,
                          const Real                  a_time,
                          const bool                  a_homogBCs,
                          const LevelData<FArrayBox>* a_crseSPtr) const
{
    const bool exCopiersAreCached =
        a_S.ghostVect() == m_statePtr->S.ghostVect() &&
        a_S.getBoxes() == m_statePtr->S.getBoxes();

    const bool hasGhosts = (a_S.ghostVect().sum() > 0);

    // Diffusive terms need edges too.
    // Filling all of the ghosts properly in mapped coordinates is a nightmare.
    // Neum BCs need all nearby cells to be filled, including corner ghosts,
    // while the exchange of corner ghosts requires the physical BCs to be set.
    // It's a chicken and egg problem that requires iteration. Rather than do
    // that, I'm gonna cheat. Most of these extrapolated ghosts will be
    // properly overwritten. Some will stick around -- namely, the edge and
    // vertex ghosts that live at the intersection of the CFI and physical
    // boundary.
    Subspace::addHorizontalExtrusion(a_S, 0, *m_SbarPtr, 0, 1, -1.0);
    BCTools::extrapAllGhosts(a_S, 2);
    Subspace::addHorizontalExtrusion(a_S, 0, *m_SbarPtr, 0, 1, 1.0);

    // Set CFI BCs.
    // Do this first since it shouldn't depend on this level's ghosts.
    if (m_level > 0) {
        // Gather coarse and fine data structures.
        AMRNSLevel*              crsePtr   = this->crseNSPtr();
        const DisjointBoxLayout& crseGrids = crsePtr->getBoxes();

        if (a_homogBCs) {
            m_cfInterpPtr->homogInterpAtCFI(a_S);
        } else {
            // Gather the coarse level data.
            LevelData<FArrayBox> crseS;
            if (a_crseSPtr) {
                // Alias the user's crse S.
                CH_assert(a_crseSPtr->nComp() == a_S.nComp());
                aliasLevelData(crseS,
                               const_cast<LevelData<FArrayBox>*>(a_crseSPtr),
                               a_crseSPtr->interval());

            } else if (abs(m_time - crsePtr->m_time) < timeEps) {
                // Alias the coarse level's S.
                aliasLevelData(crseS,
                               const_cast<LevelData<FArrayBox>*>(crsePtr->m_qPtr),
                               crsePtr->m_statePtr->SInterval);

            } else {
                // Interpolate the coarse level's S.
                crseS.define(crseGrids, 1);
                crsePtr->m_parkPtr->qTimeInterp(crseS, a_time, m_statePtr->SComp, 0, 1);
            }

            // Interpolate to fine CFI ghosts.
            m_cfInterpPtr->interpAtCFI(a_S, crseS);
        }
    }

    // Do exchange #1 on perturbation.
    // This must be done before we potentially apply Neumann BCs.
    Subspace::addHorizontalExtrusion(a_S, 0, *m_SbarPtr, 0, 1, -1.0);
    if (hasGhosts) {
        if (exCopiersAreCached) {
            a_S.exchange(m_statePtr->qExCopier);
        } else {
            Copier cp;
            State::qDefineExCopier(cp, a_S);
            a_S.exchange(cp);
        }
    }
    Subspace::addHorizontalExtrusion(a_S, 0, *m_SbarPtr, 0, 1, 1.0);

    // Set physical BCs.
    // This is where we would potentially apply Neumann BCs.
    this->setSalinityPhysBC(a_S, a_time, a_homogBCs);

    // Do exchange #2 on perturbation.
    // This must be done after all other ghosts are filled.
    Subspace::addHorizontalExtrusion(a_S, 0, *m_SbarPtr, 0, 1, -1.0);
    if (hasGhosts) {
        if (exCopiersAreCached) {
            a_S.exchange(m_statePtr->qExCornerCopier);
        } else {
            CornerCopier ccp;
            State::qDefineExCornerCopier(ccp, a_S);
            a_S.exchange(ccp);
        }
    }

    // Fill ghosts at corners of domain. This must happen last!
    BCTools::extrapDomainCorners(a_S, 2);
    Subspace::addHorizontalExtrusion(a_S, 0, *m_SbarPtr, 0, 1, 1.0);
}


// -----------------------------------------------------------------------------
// Set all BCs on all scalars.
// -----------------------------------------------------------------------------
void
AMRNSLevel::setScalarBC(LevelData<FArrayBox>&       a_s,
                        const Real                  a_time,
                        const bool                  a_homogBCs,
                        const LevelData<FArrayBox>* a_crseSPtr) const
{
    CH_assert(a_s.nComp() == this->numScalars());

    const bool exCopiersAreCached =
        a_s.ghostVect() == m_statePtr->scalars.ghostVect() &&
        a_s.getBoxes() == m_statePtr->scalars.getBoxes();

    const bool hasGhosts = (a_s.ghostVect().sum() > 0);

    // Filling all of the ghosts properly in mapped coordinates is a nightmare.
    // Neum BCs need all nearby cells to be filled, including corner ghosts,
    // while the exchange of corner ghosts requires the physical BCs to be set.
    // It's a chicken and egg problem that requires iteration. Rather than do
    // that, I'm gonna cheat. Most of these extrapolated ghosts will be
    // properly overwritten. Some will stick around -- namely, the edge and
    // vertex ghosts that live at the intersection of the CFI and physical
    // boundary.
    BCTools::extrapAllGhosts(a_s, 2);

    // Set CFI BCs.
    // Do this first since it shouldn't depend on this level's ghosts.
    if (m_level > 0) {
        // Gather coarse and fine data structures.
        AMRNSLevel*              crsePtr   = this->crseNSPtr();
        const DisjointBoxLayout& crseGrids = crsePtr->getBoxes();
        const int startComp = m_statePtr->scalarsInterval.begin();
        const int numComps  = m_statePtr->numScalars;

        if (a_homogBCs) {
            m_cfInterpPtr->homogInterpAtCFI(a_s);
        } else {
            // Gather the coarse level data.
            LevelData<FArrayBox> crseS;
            if (a_crseSPtr) {
                // Alias the user's crse s.
                CH_assert(a_crseSPtr->nComp() == a_s.nComp());
                aliasLevelData(crseS,
                               const_cast<LevelData<FArrayBox>*>(a_crseSPtr),
                               a_crseSPtr->interval());

            } else if (abs(m_time - crsePtr->m_time) < timeEps) {
                // Alias the coarse level's s.
                aliasLevelData(crseS,
                               const_cast<LevelData<FArrayBox>*>(crsePtr->m_qPtr),
                               crsePtr->m_statePtr->scalarsInterval);

            } else {
                // Interpolate the coarse level's s.
                crseS.define(crseGrids, numComps);
                crsePtr->m_parkPtr->qTimeInterp(crseS, a_time, startComp, 0, numComps);
            }

            // Interpolate to fine CFI ghosts.
            m_cfInterpPtr->interpAtCFI(a_s, crseS);
        }
    }

    // Do exchange #1.
    // This must be done before we potentially call BCTools::neum.
    if (hasGhosts) {
        if (exCopiersAreCached) {
            a_s.exchange(m_statePtr->qExCopier);
        } else {
            Copier cp;
            State::qDefineExCopier(cp, a_s);
            a_s.exchange(cp);
        }
    }

    // Set physical BCs
    this->setScalarsPhysBC(a_s, a_time, a_homogBCs);

    // Do exchange #2.
    // This must be done after all other ghosts are filled.
    if (hasGhosts) {
        if (exCopiersAreCached) {
            a_s.exchange(m_statePtr->qExCornerCopier);
        } else {
            CornerCopier ccp;
            State::qDefineExCornerCopier(ccp, a_s);
            a_s.exchange(ccp);
        }
    }

    // Fill ghosts at corners of domain. This must happen last!
    BCTools::extrapDomainCorners(a_s, 2);
}
