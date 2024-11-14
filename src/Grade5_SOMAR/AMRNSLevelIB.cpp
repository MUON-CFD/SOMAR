#include "AMRNSLevel.H"
// #include "IBMFCOp.H"
#include "MiscUtils.H"
#include "Debug.H"

// -----------------------------------------------------------------------------
void
AMRNSLevel::constructIB(const Real a_time)
{
    // Default = no-op.
}

// -----------------------------------------------------------------------------
void
AMRNSLevel::moveIB(const Real a_time)
{
    // Default = no-op.
}

// -----------------------------------------------------------------------------
void
AMRNSLevel::applyIBForcing(LevelData<FluxBox>& a_vel,
                           const Real          a_time,
                           const Real          a_dt,
                           const bool          a_makeFast)
{
    const auto& ctxib = ProblemContext::getInstance()->ib;

    // Do we have an IB?
    if (!ctxib.doIB) return;

    // Sanity checks
    CH_assert(m_IBPtr);
    debugCheckValidFaceOverlap(a_vel);

    // 1. Move the IB, if needed.
    if (const RealVect movingIBVel = this->computeIBVel(a_time, a_dt);
        movingIBVel.vectorLengthSq() > smallReal) {

        const IB* crseIBPtr =
            (m_level > 0) ? &(this->crseNSPtr()->getIB()) : nullptr;
        m_IBPtr->move(movingIBVel, a_dt, crseIBPtr);
    }

    // Gather amrVel
    Vector<LevelData<FluxBox>*> amrVel;
    this->allocateAndDefine(amrVel, 1, IntVect::Unit, 0, m_level);
    for (int lev = 0; lev < m_level; ++lev) {
        const AMRNSLevel* levPtr   = this->getLevel(lev);
        constexpr bool    setBCs   = false;
        constexpr bool    homogBCs = false;

        levPtr->fillVelocity(*amrVel[lev], a_time, setBCs);
        if (lev > 0) {
            levPtr->setVelBC(*amrVel[lev], a_time, homogBCs, amrVel[lev - 1]);
        } else {
            levPtr->setVelBC(*amrVel[lev], a_time, homogBCs);
        }
    }
    {
        const int      lev      = m_level;
        constexpr bool homogBCs = false;

        debugInitLevel(*amrVel[lev]);
        a_vel.copyTo(*amrVel[lev]);
        if (lev > 0) {
            this->setVelBC(*amrVel[lev], a_time, homogBCs, amrVel[lev - 1]);
        } else {
            this->setVelBC(*amrVel[lev], a_time, homogBCs);
        }
    }
    for (int lev = 0; lev <= m_level; ++lev) {
        this->sendToAdvectingVelocity(*amrVel[lev], *amrVel[lev]);
    }

    // 3. Send the velocity to the IB's grids.
    LevelData<FluxBox> velIB;
    m_IBPtr->sendToIBGrids(velIB, amrVel);
    this->sendToCartesianVelocity(velIB, velIB);

    // 4. Apply the boundary forcing to the velocity.
    const size_t numIters = ProblemContext::getInstance()->ib.maxExplicitIters;
    for (size_t iter = 0; iter < numIters; ++iter) {
        // Set ghosts
        constexpr bool homogBCs = false;
        this->setVelPhysBC(velIB, a_time, homogBCs);
        velIB.exchange();
        // BCTools::extrapDomainCorners(velIB, 2); // TODO: Do we need this?

        // Apply
        m_IBPtr->applyBoundaryForce(velIB, ctxib.forceScale, a_makeFast);
    }

    // 5. Send forced velocity back to user's grids.
    m_IBPtr->sendToUserGrids(a_vel, velIB);
    if (m_level > 0) {
        this->sendToCartesianVelocity(*amrVel[m_level - 1], *amrVel[m_level - 1]);
        this->setVelBC(a_vel, a_time, false, amrVel[m_level - 1]);
    } else {
        this->setVelBC(a_vel, a_time, false);
    }
    checkForValidNAN(a_vel);

    // Free memory
    this->deallocate(amrVel);

    // // TEMP!!!
    // IO::writeHDF5("test.hdf5", a_vel, *m_levGeoPtr, m_dt, a_time);
    // barrier();
    // CH_verify(false);
}
