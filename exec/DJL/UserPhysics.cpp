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
#include "UserPhysics.H"
#include "BilinearInterp.H"
#include "Convert.H"
#include "FourthOrder.H"
#include "IO.H"
#include "Subspace.H"
#include "Comm.H"


// Static members
UserPhysics::SrcData UserPhysics::s_src;


//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
UserPhysics::UserPhysics() : AMRNSLevel::AMRNSLevel()
{
    UserPhysics::setStaticMembers();
}


//------------------------------------------------------------------------------
// Virtual destructor for good measure.
//------------------------------------------------------------------------------
UserPhysics::~UserPhysics()
{
}


// -----------------------------------------------------------------------------
// Returns the names of the scalars. This will be used to put scalar names
// in the plot files.
// -----------------------------------------------------------------------------
std::string
UserPhysics::getScalarName(const int a_comp) const
{
    switch (a_comp) {
        // If you have scalars, add cases here.
        default:
            MayDay::Error("UserPhysics::getScalarName: a_comp out of range.");
    }

    return "Undefined";
}


// -----------------------------------------------------------------------------
void
UserPhysics::setStaticMembers()
{
    if (s_src.isDefined()) return;

    // Open the input file.
    const std::string filename = "DJL_CC_1024x1024_APE205.txt";
    std::ifstream inFile(filename);
    if (!inFile.is_open()) {
        MAYDAYERROR("Could not open file " << filename << ".");
    }

    // L, H, nx, nz, A, c
    IO::readASCII(s_src.L, inFile);
    IO::readASCII(s_src.H, inFile);
    IO::readASCII(s_src.nx, inFile);
    IO::readASCII(s_src.nz, inFile);
    IO::readASCII(s_src.A, inFile);
    IO::readASCII(s_src.c, inFile);

    // x
    s_src.x = std::vector<Real>(s_src.nx, quietNAN);
    IO::readASCII(s_src.x, inFile);

    // z
    s_src.z = std::vector<Real>(s_src.nz, quietNAN);
    IO::readASCII(s_src.z, inFile);

    // CC eta
    const Box srcEtaBox(
        IntVect::Zero,
        (s_src.nx - 1) * BASISV(0) + (s_src.nz - 1) * BASISV(SpaceDim - 1));
    s_src.etaFAB.define(srcEtaBox, 1);
    IO::readASCII(s_src.etaFAB, 0, srcEtaBox, inFile);

    // We are done with the file.
    inFile.close();

    POUT(s_src.L);
    POUT(s_src.H);
    POUT(s_src.A);
    POUT(s_src.c);
}


// -----------------------------------------------------------------------------
std::vector<Real>
UserPhysics::getLayer_eta(const Real a_z)
{
    if (a_z <= s_src.z.front() || s_src.z.back() <= a_z) {
        return std::vector<Real>(s_src.nx, 0.0);
    }

    // Find the layers surrounding a_z.
    const size_t srck = std::lower_bound(s_src.z.begin(), s_src.z.end(), a_z) -
                        s_src.z.begin() - 1;

    CH_verify(srck < s_src.z.size() - 1);
    CH_verify(s_src.z[srck] < a_z);
    CH_verify(a_z < s_src.z[srck + 1]);

    std::vector<Real> eta(s_src.nx, quietNAN);

    if (srck == 0 || srck == s_src.nz - 2) {
        // Linear interp
        const Real cl = (s_src.z[srck + 1] - a_z) / (s_src.z[srck + 1] - s_src.z[srck]);
        const Real cr = 1.0 - cl;

        IntVect minIV = (srck + 0) * BASISV(SpaceDim - 1);
        IntVect maxIV = (srck + 1) * BASISV(SpaceDim - 1);
        for (size_t i = 0; i < eta.size(); ++i) {
            eta[i] = cl * s_src.etaFAB(minIV) + cr * s_src.etaFAB(maxIV);
            ++minIV[0];
            ++maxIV[0];
        }
    } else {
        // Cubic interp
        const Real zeta = (a_z - 0.5 * (s_src.z[srck + 1] + s_src.z[srck]));

        IntVect ivll = (srck - 1) * BASISV(SpaceDim - 1);
        IntVect ivl  = (srck + 0) * BASISV(SpaceDim - 1);
        IntVect ivr  = (srck + 1) * BASISV(SpaceDim - 1);
        IntVect ivrr = (srck + 2) * BASISV(SpaceDim - 1);

        for (size_t i = 0; i < eta.size(); ++i) {
            const Real etall = s_src.etaFAB(ivll);
            const Real etal  = s_src.etaFAB(ivl);
            const Real etar  = s_src.etaFAB(ivr);
            const Real etarr = s_src.etaFAB(ivrr);

            const Real a0 = (-etall +  9.0 * etal +  9.0 * etar - etarr) / 16.0;
            const Real a1 = ( etall - 27.0 * etal + 27.0 * etar - etarr) / 48.0;
            const Real a2 = ( etall -        etal -        etar + etarr) / 16.0;
            const Real a3 = (-etall +  3.0 * etal -  3.0 * etar + etarr) / 48.0;

            eta[i] = a0 + zeta * (a1 + zeta * (a2 + zeta * a3));

            ++ivll[0];
            ++ivl[0];
            ++ivr[0];
            ++ivrr[0];
        }
    }

    return eta;
}


// -----------------------------------------------------------------------------
std::vector<Real>
UserPhysics::getLayer_eta_z(const Real a_z)
{
    if (a_z <= s_src.z.front() || s_src.z.back() <= a_z) {
        return std::vector<Real>(s_src.nx, 0.0);
    }

    // Find the layers surrounding a_z.
    const size_t srck = std::lower_bound(s_src.z.begin(), s_src.z.end(), a_z) -
                        s_src.z.begin() - 1;

    CH_verify(srck < s_src.z.size() - 1);
    CH_verify(s_src.z[srck] < a_z);
    CH_verify(a_z < s_src.z[srck + 1]);

    std::vector<Real> eta_z(s_src.nx, quietNAN);

    if (srck == 0 || srck == s_src.nz - 2) {
        // Deriv of linear interp
        IntVect minIV = srck * BASISV(SpaceDim - 1);
        IntVect maxIV = (srck + 1) * BASISV(SpaceDim - 1);
        const Real scale = 1.0 / (s_src.z[srck + 1] - s_src.z[srck]);

        for (size_t i = 0; i < eta_z.size(); ++i) {
            eta_z[i] = (s_src.etaFAB(maxIV) - s_src.etaFAB(minIV)) * scale;
            ++minIV[0];
            ++maxIV[0];
        }
    } else {
        // Deriv of cubic interp
        const Real zeta_z = 2.0 / (s_src.z[srck + 1] - s_src.z[srck]);
        const Real zeta   = (a_z - 0.5 * (s_src.z[srck + 1] + s_src.z[srck]));

        IntVect ivll = (srck - 1) * BASISV(SpaceDim - 1);
        IntVect ivl  = (srck + 0) * BASISV(SpaceDim - 1);
        IntVect ivr  = (srck + 1) * BASISV(SpaceDim - 1);
        IntVect ivrr = (srck + 2) * BASISV(SpaceDim - 1);

        for (size_t i = 0; i < eta_z.size(); ++i) {
            const Real etall = s_src.etaFAB(ivll);
            const Real etal  = s_src.etaFAB(ivl);
            const Real etar  = s_src.etaFAB(ivr);
            const Real etarr = s_src.etaFAB(ivrr);

            // const Real a0 = (-etall +  9.0 * etal +  9.0 * etar - etarr) / 16.0;
            const Real a1 = ( etall - 27.0 * etal + 27.0 * etar - etarr) / 48.0;
            const Real a2 = ( etall -        etal -        etar + etarr) / 16.0;
            const Real a3 = (-etall +  3.0 * etal -  3.0 * etar + etarr) / 48.0;

            eta_z[i] = (a1 + zeta * (2.0 * a2 + zeta * 3.0 * a3)) * zeta_z;

            ++ivll[0];
            ++ivl[0];
            ++ivr[0];
            ++ivrr[0];
        }
    }

    return eta_z;
}


// -----------------------------------------------------------------------------
void
UserPhysics::computeVelSoln(LevelData<FluxBox>& a_vel, const Real a_time) const
{
    const auto&     grids  = m_levGeoPtr->getBoxes();
    const auto&     domBox = m_levGeoPtr->getDomainBox();
    const Real      minX   = s_src.x.front();
    const Real      maxX   = s_src.x.back();

    if constexpr (SpaceDim == 2) {

        IntVect iv = domBox.smallEnd();
        for (size_t k = 0; k < (size_t)domBox.size(SpaceDim-1); ++k, ++iv[SpaceDim - 1]) {

            // Create spline of this layer's eta
            const Real z = m_levGeoPtr->getCellX(SpaceDim - 1, iv);

            CubicSpline eta_z;
            eta_z.solve(UserPhysics::getLayer_eta_z(z), s_src.x, 0.0, 0.0);

            for (DataIterator dit(grids); dit.ok(); ++dit) {
                const Box ccValid = grids[dit];
                if (iv[SpaceDim - 1] < ccValid.smallEnd(SpaceDim - 1)) continue;
                if (ccValid.bigEnd(SpaceDim - 1) < iv[SpaceDim - 1]) continue;

                // Compute u = cos(theta) * c * eta_z
                iv[0] = ccValid.smallEnd(0);
                size_t i = iv[0] - domBox.smallEnd(0);
                for (; iv[0] <= ccValid.bigEnd(0) + 1; ++iv[0], ++i) {
                    const Real x      = m_levGeoPtr->getFaceX(0, iv, 0) - s_xOffset - s_src.c * a_time;
                    const Real interp = (minX <= x && x <= maxX) ? eta_z.interp(x) : 0.0;

                    a_vel[dit][0](iv) = s_src.c * interp;
                } // i

                if constexpr (SpaceDim == 3) {
                    a_vel[dit][1].setVal(0.0);
                }
            } // dit
        } // k

        // Loop over each nodal vertical layer to compute w.
        iv = domBox.smallEnd();
        for (size_t k = 0; k < (size_t)domBox.size(SpaceDim-1) + 1; ++k, ++iv[SpaceDim - 1]) {
            // Create spline of this layer's eta
            const Real z = m_levGeoPtr->getFaceX(SpaceDim - 1, iv, SpaceDim - 1);

            CubicSpline eta;
            eta.solve(UserPhysics::getLayer_eta(z), s_src.x, 0.0, 0.0);

            for (DataIterator dit(grids); dit.ok(); ++dit) {
                const Box fcValid = grids[dit].surroundingNodes(SpaceDim - 1);
                if (iv[SpaceDim - 1] < fcValid.smallEnd(SpaceDim - 1)) continue;
                if (fcValid.bigEnd(SpaceDim - 1) < iv[SpaceDim - 1]) continue;

                // Compute w = -c * eta_x
                iv[0] = fcValid.smallEnd(0);
                size_t i = iv[0] - domBox.smallEnd(0);
                for (; iv[0] <= fcValid.bigEnd(0); ++iv[0], ++i) {
                    const Real x      = m_levGeoPtr->getFaceX(0, iv, SpaceDim - 1) - s_xOffset - s_src.c * a_time;
                    const Real interp = (minX <= x && x <= maxX) ? eta.interpFirstDeriv(x) : 0.0;

                    a_vel[dit][SpaceDim - 1](iv) = -s_src.c * interp;
                } // i
            } // dit
        } // k

    } else {

        const Real  cosTheta = cos(s_rotAngle);
        const Real  sinTheta = sin(s_rotAngle);

        IntVect iv = domBox.smallEnd();
        for (size_t k = 0; k < (size_t)domBox.size(SpaceDim-1); ++k, ++iv[SpaceDim - 1]) {

            // Create spline of this layer's eta
            const Real z = m_levGeoPtr->getCellX(SpaceDim - 1, iv);

            CubicSpline eta;
            eta.solve(UserPhysics::getLayer_eta(z), s_src.x, 0.0, 0.0);

            CubicSpline eta_z;
            eta_z.solve(UserPhysics::getLayer_eta_z(z), s_src.x, 0.0, 0.0);

            for (DataIterator dit(grids); dit.ok(); ++dit) {
                const Box ccValid = grids[dit];
                if (iv[SpaceDim - 1] < ccValid.smallEnd(SpaceDim - 1)) continue;
                if (ccValid.bigEnd(SpaceDim - 1) < iv[SpaceDim - 1]) continue;

                // Compute u = cos(theta) * c * eta_z
                iv[1] = ccValid.smallEnd(1);
                int j = iv[1] - domBox.smallEnd(1);
                for (; iv[1] <= ccValid.bigEnd(1); ++iv[1], ++j) {

                    iv[0] = ccValid.smallEnd(0);
                    size_t i = iv[0] - domBox.smallEnd(0);
                    for (; iv[0] <= ccValid.bigEnd(0) + 1; ++iv[0], ++i) {

                        const Real x = m_levGeoPtr->getFaceX(0, iv, 0) - s_xOffset;
                        const Real y = m_levGeoPtr->getFaceX(1, iv, 0) - s_yOffset;

                        const Real xprime =  x * cosTheta + y * sinTheta;
                        const Real yprime = -x * sinTheta + y * cosTheta;

                        const Real env    = UserPhysics::extrusionEnvelope(yprime);
                        const Real interp =
                            (minX <= xprime && xprime <= maxX)
                                ? eta_z.interp(xprime)
                                : 0.0;

                        a_vel[dit][0](iv) = env * cosTheta * s_src.c * interp;
                    }  // i
                } // j

                // Compute v = -sin(theta) * c * eta_z
                iv[1] = ccValid.smallEnd(1);
                j = iv[1] - domBox.smallEnd(1);
                for (; iv[1] <= ccValid.bigEnd(1) + 1; ++iv[1], ++j) {

                    iv[0] = ccValid.smallEnd(0);
                    size_t i = iv[0] - domBox.smallEnd(0);
                    for (; iv[0] <= ccValid.bigEnd(0); ++iv[0], ++i) {

                        const Real x = m_levGeoPtr->getFaceX(0, iv, 1) - s_xOffset;
                        const Real y = m_levGeoPtr->getFaceX(1, iv, 1) - s_yOffset;

                        const Real xprime =  x * cosTheta + y * sinTheta;
                        const Real yprime = -x * sinTheta + y * cosTheta;

                        const Real env = UserPhysics::extrusionEnvelope(yprime);
                        const Real interp =
                            (minX <= xprime && xprime <= maxX)
                                ? eta_z.interp(xprime)
                                : 0.0;

                        a_vel[dit][1](iv) = env * sinTheta * s_src.c * interp;
                    }  // i
                } // j
            } // dit
        } // k

        // Loop over each nodal vertical layer to compute w.
        iv = domBox.smallEnd();
        for (size_t k = 0; k < (size_t)domBox.size(SpaceDim-1) + 1; ++k, ++iv[SpaceDim - 1]) {

            // Create spline of this layer's eta
            const Real z = m_levGeoPtr->getFaceX(SpaceDim - 1, iv, SpaceDim - 1);

            CubicSpline eta;
            eta.solve(UserPhysics::getLayer_eta(z), s_src.x, 0.0, 0.0);

            for (DataIterator dit(grids); dit.ok(); ++dit) {
                const Box fcValid = grids[dit].surroundingNodes(SpaceDim - 1);
                if (iv[SpaceDim - 1] < fcValid.smallEnd(SpaceDim - 1)) continue;
                if (fcValid.bigEnd(SpaceDim - 1) < iv[SpaceDim - 1]) continue;

                // Compute w = -c * eta_x
                iv[1] = fcValid.smallEnd(1);
                int j = iv[1] - domBox.smallEnd(1);
                for (; iv[1] <= fcValid.bigEnd(1); ++iv[1], ++j) {

                    iv[0] = fcValid.smallEnd(0);
                    size_t i = iv[0] - domBox.smallEnd(0);
                    for (; iv[0] <= fcValid.bigEnd(0); ++iv[0], ++i) {

                        const Real x = m_levGeoPtr->getFaceX(0, iv, SpaceDim - 1) - s_xOffset;
                        const Real y = m_levGeoPtr->getFaceX(1, iv, SpaceDim - 1) - s_yOffset;

                        const Real xprime =  x * cosTheta + y * sinTheta;
                        const Real yprime = -x * sinTheta + y * cosTheta;

                        const Real env = UserPhysics::extrusionEnvelope(yprime);
                        const Real interp =
                            (minX <= xprime && xprime <= maxX)
                                ? eta.interpFirstDeriv(xprime)
                                : 0.0;

                        a_vel[dit][SpaceDim - 1](iv) = -env * s_src.c * interp;
                    }  // i
                } // j
            } // dit
        } // k
    }
}


// -----------------------------------------------------------------------------
void
UserPhysics::computeBSoln(LevelData<FArrayBox>& a_b, const Real a_time) const
{
    const auto&     grids  = m_levGeoPtr->getBoxes();
    const auto&     domBox = m_levGeoPtr->getDomainBox();
    const Real      minX   = s_src.x.front();
    const Real      maxX   = s_src.x.back();

    // Loop over each zonal vertical layer to compute T = b.
    IntVect iv = domBox.smallEnd();
    for (size_t k = 0; k < (size_t)domBox.size(SpaceDim-1); ++k, ++iv[SpaceDim - 1]) {

        // Create spline of this layer's eta
        const Real z = m_levGeoPtr->getCellX(SpaceDim - 1, iv);

        CubicSpline eta;
        eta.solve(UserPhysics::getLayer_eta(z), s_src.x, 0.0, 0.0);

        CubicSpline eta_z;
        eta_z.solve(UserPhysics::getLayer_eta_z(z), s_src.x, 0.0, 0.0);

        for (DataIterator dit(grids); dit.ok(); ++dit) {
            const Box ccValid = grids[dit];
            if (iv[SpaceDim - 1] < ccValid.smallEnd(SpaceDim - 1)) continue;
            if (ccValid.bigEnd(SpaceDim - 1) < iv[SpaceDim - 1]) continue;

            // Compute T = bbar(z-eta(x'))
#if CH_SPACEDIM == 2
            iv[0] = ccValid.smallEnd(0);
            size_t i = iv[0] - domBox.smallEnd(0);
            for (; iv[0] <= ccValid.bigEnd(0); ++iv[0], ++i) {
                const Real x      = m_levGeoPtr->getCellX(0, iv) - s_xOffset - s_src.c * a_time;
                const Real etaVal = (minX <= x && x <= maxX) ? eta.interp(x) : 0.0;

                a_b[dit](iv) = m_TbarSplinePtr->interp(z - etaVal);
            } // i
#else
            const Real cosTheta = cos(s_rotAngle);
            const Real sinTheta = sin(s_rotAngle);

            iv[1] = ccValid.smallEnd(1);
            int j = iv[1] - domBox.smallEnd(1);
            for (; iv[1] <= ccValid.bigEnd(1); ++iv[1], ++j) {

                iv[0] = ccValid.smallEnd(0);
                size_t i = iv[0] - domBox.smallEnd(0);
                for (; iv[0] <= ccValid.bigEnd(0); ++iv[0], ++i) {

                    const Real x = m_levGeoPtr->getCellX(0, iv) - s_xOffset;
                    const Real y = m_levGeoPtr->getCellX(1, iv) - s_yOffset;

                    const Real xprime =  x * cosTheta + y * sinTheta;
                    const Real yprime = -x * sinTheta + y * cosTheta;

                    const Real env    = UserPhysics::extrusionEnvelope(yprime);
                    const Real Tbar   = (*m_TbarPtr)(iv * BASISV(SpaceDim - 1));
                    const Real etaVal = (minX <= xprime && xprime <= maxX) ? eta.interp(xprime) : 0.0;
                    const Real interp = m_TbarSplinePtr->interp(z - etaVal);

                    a_b[dit](iv) = env * interp + (1.0 - env) * Tbar;
                }  // i
            } // j
#endif
        } // dit
    } // k
}


// -----------------------------------------------------------------------------
void
UserPhysics::setICs(State& a_state)
{
    BEGIN_FLOWCHART();

    // In this simple example, remember that there is no temperature.
    // T is just a proxy for b. That is why we fill T with b values.
    // In other words, we set a_state.T with b_total.
    this->computeVelSoln(a_state.vel, 0.0);
    this->computeBSoln(a_state.T, 0.0);

    this->setBC(a_state, 0.0);
    FourthOrder::nodalToAvg(a_state.vel);
    for (DataIterator dit(a_state.grids); dit.ok(); ++dit) {
        FourthOrder::nodalToAvg(a_state.T[dit], a_state.T[dit].box());
    }

    LayoutTools::averageOverlappingValidFaces(a_state.vel);
    if (m_level > 0) {
        this->getLevel(0)->averageDownToThis(m_level);
    }
    debugCheckValidFaceOverlap(a_state.vel);
}


// -----------------------------------------------------------------------------
void
UserPhysics::tagCells(IntVectSet& a_tags)
{
    // AMRNSLevel::tagCells(a_tags);

    const ProblemContext* ctx       = ProblemContext::getInstance();
    const int             growTags  = ctx->amr.growTags;
    const int             numLevels = ctx->amr.numLevels;
    const auto&           grids     = m_levGeoPtr->getBoxes();
    const auto            ez        = BASISV(SpaceDim - 1);
    const auto            dZeta     = m_levGeoPtr->getDXi(SpaceDim - 1);

    // Add tags based on projected b'.
    if (1) {
        LevelData<FArrayBox> Tpert(grids, 1);
        this->fillTemperature(Tpert, m_time, false, true);

        IntVectSet newTags;
        for (DataIterator dit(grids); dit.ok(); ++dit) {
            const auto&   bpertFAB = Tpert[dit];
            const auto&   phi1FAB   = *m_phi1Ptr;
            const auto&   JFAB      = m_levGeoPtr->getCCJ()[dit];
            const Box     valid     = grids[dit];
            const Box     flatValid = Subspace::flattenBox(valid, SpaceDim - 1);
            const IntVect iv_zmin   = valid.smallEnd(SpaceDim - 1) * ez;
            const IntVect iv_zmax   = valid.bigEnd(SpaceDim - 1) * ez;

            for (BoxIterator flatbit(flatValid); flatbit.ok(); ++flatbit) {
                const Box intBox(flatbit() + iv_zmin, flatbit() + iv_zmax);

                Real bproj = 0.0;
                for (BoxIterator bit(intBox); bit.ok(); ++bit) {
                    const IntVect& iv = bit();
                    bproj += bpertFAB(iv) * phi1FAB(iv * ez) * JFAB(iv) * dZeta;
                }
                bproj = std::abs(bproj);

                if (bproj >= s_bProjTagTol) {
                    newTags |= intBox;
                }
            }
        }

        if (!newTags.isEmpty()) {
            if (m_level == numLevels - 2) {
                newTags.grow(growTags);
            }
            newTags &= m_levGeoPtr->getDomainBox();
            a_tags |= newTags;
        }
    }

    // Add tags based on p.
    if (0) {
    // if (!RealCmp::isZero(m_time) && !RealCmp::isZero(s_pTagTol)) {
        const auto& domBox = grids.physDomain().domainBox();

        const IntVect iv_zmin = domBox.smallEnd(SpaceDim - 1) * ez;
        const IntVect iv_zmax = domBox.bigEnd(SpaceDim - 1) * ez;
        const IntVect zOffset = ez * (domBox.smallEnd(SpaceDim - 1) +
                                      0.80 * domBox.size(SpaceDim - 1));

        IntVectSet newTags;
        for (DataIterator dit(grids); dit.ok(); ++dit) {
            const auto&   dataFAB   = (*m_pPtr)[dit];
            const Box     valid     = grids[dit];
            const Box     flatValid = Subspace::flattenBox(valid, SpaceDim - 1);

            for (BoxIterator flatbit(flatValid); flatbit.ok(); ++flatbit) {
                const Box intBox(flatbit() + iv_zmin, flatbit() + iv_zmax);
                const IntVect iv = flatbit() + zOffset;

                if (dataFAB(iv) > s_pTagTol) {
                    newTags |= intBox;
                }
            }
        }

        if (!newTags.isEmpty()) {
            if (m_level == numLevels - 2) {
                newTags.grow(growTags);
            }
            newTags &= m_levGeoPtr->getDomainBox();
            a_tags |= newTags;
        }
    }

    // Write diagnostic info.
    {
        int numTags = a_tags.numPts();
        Comm::reduce(numTags, MPI_SUM);
        IO::tout(0) << "Level " << m_level << " tagging:\n";
        IO::tout(0) << "\tratio = "
                    << ((Real)numTags /
                        (Real)m_levGeoPtr->getDomainBox().numPts())
                    << '\n';
        IO::tout(0) << std::flush;
    }
}
