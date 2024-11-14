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
#include "LockExchangePhysics.H"
#include "BCTools.H"
#include "SetValLevel.H"
#include "BoxIterator.H"
#include "Debug.H"
#include "Zipper.H"
#include "FETools.H"
#include "IB.H"
#include "DEM.H"
#include "Comm.H"

#include "Convert.H"
#include "AnisotropicRefinementTools.H"

#include "VTK.H"


//------------------------------------------------------------------------------
LockExchangePhysics::LockExchangePhysics()
: AMRNSLevel::AMRNSLevel()
{
}


//------------------------------------------------------------------------------
LockExchangePhysics::~LockExchangePhysics()
{
}


//------------------------------------------------------------------------------
std::string
LockExchangePhysics::getScalarName(const int a_comp) const
{
    CH_assert(a_comp >= 0);

    switch (a_comp) {
        case X_TRACER: return "x_Tracer";
        case Y_TRACER: return "y_Tracer";
        case Z_TRACER: return "z_Tracer";
        case LAMBDA:   return "lambda";
    }

    MAYDAYERROR("LockExchangePhysics::getScalarName: a_comp out of range.");
    return "Undefined";
}


//------------------------------------------------------------------------------
void
LockExchangePhysics::setICs(State& a_state)
{
    const DisjointBoxLayout& grids = a_state.grids;
    DataIterator dit = grids.dataIterator();

    // Start clean.
    setValLevel(a_state.vel, 0.0);
    setValLevel(a_state.p, 0.0);
    setValLevel(a_state.q, 0.0);

    // tracers
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox tracers(Interval(X_TRACER, Z_TRACER), a_state.scalars[dit]);
        m_levGeoPtr->fill_physCoor(tracers);
    }

    //lambda
    for (auto [s] : zip(a_state.scalars))
        s.setVal(1.0,LAMBDA);



    // for (dit.reset(); dit.ok(); ++dit) {
    //     a_state.scalars[dit].setVal(1.0, LAMBDA);
    // }

    // b = T
    LevelData<FArrayBox>& b = a_state.T;

    const Real xhalf  = 0.0;   // Interface location
    const Real deltax = 0.05;  // Interface width

    // const Real pertA = 0.0;
    const Real pertA = ((SpaceDim > 2)? 0.025: 0.0);
    const Real pertK = 2.0 * Pi / m_levGeoPtr->getDomainLength(1);

    const Real bmin = 0.0;
    const Real bmax = 1.0;

    for (auto [ccRegion,bFAB] : zip(grids,b))
    {
    // for (dit.reset(); dit.ok(); ++dit) {
    //     FArrayBox& bFAB = b[dit];
    //     const Box& ccRegion = bFAB.box();
        BoxIterator bit(ccRegion);
        Real x, y, ifx, frac;

        // Compute Cartesian cell coordinates
        FArrayBox posFAB(ccRegion, SpaceDim);
        m_levGeoPtr->fill_physCoor(posFAB);

        // Loop over ccRegion and set bFAB.
        for (bit.reset(); bit.ok(); ++bit) {
            const IntVect& cc = bit();
            x   = posFAB(cc, 0);
            y   = posFAB(cc, 1);
            ifx = xhalf + pertA * sin(pertK * y);

            frac = 0.5 * tanh((x - ifx) / deltax) + 0.5;
            bFAB(cc) = bmin + (bmax - bmin) * frac;
        }
    }
}



//------------------------------------------------------------------------------
BCTools::BCFunction*
LockExchangePhysics::createScalarsPhysBC(const int             a_dir,
                                         const Side::LoHiSide& a_side) const
{
    struct ScalarBC: public BCTools::BCFunction
    {
        // The BC will be set according to alpha * s + beta * ds/dn = B.
        // Note that ds/dn is an outward normal derivative!
        //
        // We want s = {x,y,z,1}. So alpha = 1, beta = 0, B = {x,y,z,1}.
        virtual void
        operator()(FArrayBox&            a_alpha,
                   FArrayBox&            a_beta,
                   FArrayBox&            a_bcFAB,
                   const FArrayBox&      /*a_stateFAB*/,
                   const FArrayBox&      a_xFAB,
                   const DataIndex&      /*a_di*/,
                   const int             /*a_bdryDir*/,
                   const Side::LoHiSide& /*a_side*/,
                   const Real            /*a_time*/,
                   const bool            a_homogBCs) const override
        {
            a_alpha.setVal(1.0);
            a_beta.setVal(0.0);
            if (!a_homogBCs) {
                a_bcFAB.copy(a_xFAB, 0, 0, SpaceDim);
                a_bcFAB.setVal(1.0, SpaceDim);
            }
        }
    };

    return new ScalarBC;
}


// // -----------------------------------------------------------------------------
// void
// LockExchangePhysics::tagCellsInit(IntVectSet& a_tags)
// {
//     a_tags.makeEmpty();
//     if (m_IBPtr) {
//         m_IBPtr->addLocalIBStencil(a_tags);
//     }
// }


// // -----------------------------------------------------------------------------
// void
// LockExchangePhysics::tagCells(IntVectSet& a_tags)
// {
// //     CH_verify(s_hasIB);
// //     if (!m_IBPtr->hasGrids()) {
// //         m_IBPtr->regrid(*m_levGeoPtr);
// //     }
// //     CH_verify(m_IBPtr->hasGrids());
// //     CH_verify(m_IBPtr->hasStencils());

//     a_tags.makeEmpty();
//     if (m_IBPtr) {
//         m_IBPtr->addLocalIBStencil(a_tags);
//     }
//     a_tags &= this->getDomainBox();

// //     const int growTags = ProblemContext::getInstance()->amr.growTags;

// //     D_TERM(
// //     a_tags.grow(0, growTags * m_IBPtr->supportGhosts()[0]);,
// //     a_tags.grow(1, growTags * m_IBPtr->supportGhosts()[1]);,
// //     a_tags.grow(2, growTags * m_IBPtr->supportGhosts()[2]);)


// //     // const Box&               domBox = m_levGeoPtr->getDomainBox();
// //     // const RealVect&          L      = m_levGeoPtr->getDomainLength();
// //     // const RealVect&          dXi    = m_levGeoPtr->getDXi();
// //     // // const Real               ivol   = 1.0 / dXi.product();  // Should be J*dx...
// //     // // const RealVect           dA     = dXi.product() / dXi;
// //     // const DisjointBoxLayout& grids  = m_levGeoPtr->getBoxes();
// //     // DataIterator             dit    = grids.dataIterator();

// //     // // Define the IB.
// //     // RealVect xmin;
// //     // xmin[0] = -5.0;
// //     // xmin[1] = 0.0;
// //     // xmin[SpaceDim - 1] = 0.0;

// //     // RealVect xmax;
// //     // xmax[0] = xmin[0] + 2.0;
// //     // xmax[1] = L[1];
// //     // xmax[SpaceDim - 1] = xmin[SpaceDim - 1] + 0.5 * L[SpaceDim - 1];

// //     // Box ibBox(domBox);
// //     // for (int d = 0; d < SpaceDim; ++d) {
// //     //     const IntVect ed = BASISV(d);

// //     //     Box xBox = Subspace::flattenBox(ibBox, ed);
// //     //     xBox.surroundingNodes(d);

// //     //     FArrayBox xFAB(xBox, 1);
// //     //     m_levGeoPtr->getGeoSource().fill_physCoor(xFAB, 0, d, dXi);

// //     //     xBox.enclosedCells(d);
// //     //     BoxIterator bit(xBox);
// //     //     CH_verify(!xBox.isEmpty());

// //     //     Real xl, xr;
// //     //     int imin = domBox.smallEnd(d);
// //     //     int imax = domBox.bigEnd(d);
// //     //     for (bit.reset(); bit.ok(); ++bit) {
// //     //         const IntVect& cc = bit();
// //     //         xl = xFAB(cc);
// //     //         xr = xFAB(cc + ed);
// //     //         if ((xl <= xmin[d]) && (xmin[d] < xr)) {
// //     //             imin = cc[d];
// //     //         }
// //     //         if ((xl < xmax[d]) && (xmax[d] <= xr)) {
// //     //             imax = cc[d];
// //     //         }
// //     //     }

// //     //     ibBox.setRange(d, imin, imax - imin + 1);
// //     //     CH_verify(!ibBox.isEmpty());
// //     // }

// //     // ibBox.grow(8);
// //     // ibBox &= domBox;

// //     // AMRNSLevel::tagCells(a_tags);
// //     // a_tags |= ibBox;
// }


//------------------------------------------------------------------------------
void
LockExchangePhysics::constructIB(Real a_time)
{
    const ProblemContext* ctx = ProblemContext::getInstance();
    if (!ctx->ib.doIB) return;

    // BEGIN_FLOWCHART();

    // constexpr int mode = 0; // Ellipse
    constexpr int mode = 1; // Gaussian bump.
    // constexpr int mode = 2; // 3D scaled Kaena Ridge.
    // constexpr int mode = ((SpaceDim == 2) ? 1 : 2);

    if constexpr (mode == 0) { // Ellipse
        CH_verify(SpaceDim == 2);

        const size_t numElems = 2000;
        // const RealVect radius{  2.0, 0.5 };
        // const RealVect center{ -4.0, 1.0 };

        const RealVect radius{D_DECL( 1.0, 1.0, quietNAN)};
        const RealVect center{D_DECL(-3.0, 0.0, quietNAN)};

        const auto [nodes, elems] =
            FiniteElements::Utils::createEllipse(radius, center, numElems);

        if (m_level == 0) {
            m_IBPtr.reset(new IB(nodes, elems, *m_levGeoPtr));
        } else {
            const IB& crseIB = this->crseNSPtr()->getIB();
            m_IBPtr.reset(new IB(nodes, elems, *m_levGeoPtr, &crseIB));
        }

        m_IBPtr->setNoSlipBCs();
        // m_IBPtr->setFreeSlipBCs();

    } else if constexpr (mode == 1) { // Gaussian bump
        const ProblemContext* ctx = ProblemContext::getInstance();
        // const Real L = ctx->base.L[0];
        const Real H = ctx->base.L[SpaceDim - 1];

        // Create a bump.
        const size_t Nx    = 8000;
        const Real   x0    = -3.0;
        const Real   A     = H / 2.0;
        const Real   sigma = 0.5;

        // The bump is resolved within x = [x0-5*sigma, x+5*sigma].
        const Real xmin = x0 - 5.0 * sigma;
        const Real xmax = x0 + 5.0 * sigma;
        const Real dx   = (xmax - xmin) / Real(Nx - 3);
        IO::tout(0) << "FE dx = " << dx << endl;

        std::vector<Real> x(Nx), z(Nx);

        // Create the gaussian bump.
        // The first and last points are far beyond the domain extents.
        x[0] = xmin - 0.5*dx;
        z[0] = -m_levGeoPtr->getDXi(SpaceDim - 1) * 0.5;

        x[1] = xmin;
        z[1] = 0.0;

        for (size_t i = 2; i < Nx - 2; ++i) {
            x[i] = xmin + Real(i - 1) * dx;

            const Real arg = (x[i] - x0) / sigma;
            z[i] = A * exp(-0.5 * arg * arg);
        }

        x[Nx - 2] = xmax;
        z[Nx - 2] = 0.0;

        x[Nx - 1] = xmax + 0.5*dx;
        z[Nx - 1] = -m_levGeoPtr->getDXi(SpaceDim - 1) * 0.5;

        // Create finite element representation.
        const auto [nodes, elems] =
            FiniteElements::Utils::createUnstructuredLines(x, z, false);

        if (m_level == 0) {
            m_IBPtr.reset(new IB(nodes, elems, *m_levGeoPtr));
        } else {
            const IB& crseIB = this->crseNSPtr()->getIB();
            m_IBPtr.reset(new IB(nodes, elems, *m_levGeoPtr, &crseIB));
        }

        m_IBPtr->setNoSlipBCs();
        // m_IBPtr->setFreeSlipBCs();

    } else if constexpr (mode == 2) { // 3D scaled Kaena Ridge
        // Load and modify DEM.
        // DEM dem("DEM/GoM1.asc");
        // DEM dem("DEM/KaenaRidge.asc");
        // DEM dem("DEM/KaenaRidge2.asc");
        // DEM dem("DEM/LuzonStrait1.asc");
        DEM dem("DEM/MontereyCanyon1.asc");
        // DEM dem("DEM/MontereyCanyon2.asc");
        // dem.diffuse(0.8, 10);
        // dem.rescaleXY({-8.5, -6.0}, {3.5, 6.0});
        // dem.rescaleZ(0.0, 1.3);

        // dem.diffuse(0.8, 1);
        // dem.diffuse(0.8, 10);
        constexpr Real margin = 0.0;
        dem.rescaleXY({-8.5 - margin, -6.0 - margin },
                      { 3.5 + margin,  6.0 + margin });
        dem.rescaleZ(0.0, 1.3);
        // dem.refine({2,2});
        // dem.coarsen({2,2});

        // Create IB representation.
        if (m_level == 0) {
            const auto domainBox     = grow(m_levGeoPtr->getDomainBox(), 4); // BUG ????
            const auto domainExtents = m_levGeoPtr->getBoxExtents(domainBox);
            const auto feDomainPtr   = dem.createFEDomain(domainExtents);
            m_IBPtr.reset(new IB(feDomainPtr, *m_levGeoPtr));

            // m_IBPtr.reset(new IB(dem, *m_levGeoPtr));

            dem.writeHDF5("IB_3d.hdf5");
        } else {
            const IB& crseIB = this->crseNSPtr()->getIB();
            m_IBPtr.reset(new IB(dem, *m_levGeoPtr, &crseIB));
        }

        // Set BCs.
        m_IBPtr->setNoSlipBCs();
        // m_IBPtr->setFreeSlipBCs();

    } else {
        MAYDAYERROR("mode not valid.");
    }


    // Write IB to file for visualization.
    if (m_level == 0) {
        const std::string baseFileName = SpaceDim == 2 ? "IB_2d" : "IB_3d";
        FiniteElements::VTK::write(
            "./hdf5_output", baseFileName, m_IBPtr->getFELayout());
    }
}


// // -----------------------------------------------------------------------------
// void
// LockExchangePhysics::writePlotHeader(HeaderData&        a_header,
//                                      const std::string& a_filename) const
// {
//     BEGIN_FLOWCHART();

//     if (m_time > 0.0) {
//         CH_verify(m_level == 0);

//         const Real oldTime = m_parkPtr->oldTime();
//         const Real newTime = m_parkPtr->newTime();

//         Vector<LevelData<FArrayBox>*> amrPressure;
//         this->allocateAndDefine(amrPressure, 1, IntVect::Unit);
//         this->computePhysicalPressure(amrPressure, oldTime, newTime);

//         Vector<LevelData<FArrayBox>*> amrData;
//         this->allocateAndDefine(amrData, 2 + SpaceDim, IntVect::Unit);
//         for (size_t l = 0; l < amrPressure.size(); ++l) {
//             const auto& grids = amrPressure[l]->getBoxes();
//             for (DataIterator dit(grids); dit.ok(); ++dit) {
//                 FArrayBox&       dataFAB     = (*amrData[l])[dit];
//                 const FArrayBox& pressureFAB = (*amrPressure[l])[dit];
//                 const FArrayBox& phiFAB      = (getLevel(l)->getState().p)[dit];
//                 const FluxBox&   velFlub     = (getLevel(l)->getState().vel)[dit];

//                 dataFAB.copy(pressureFAB, 0, 0, 1);
//                 dataFAB.copy(phiFAB, 0, 1, 1);
//                 for (int velComp = 0; velComp < SpaceDim; ++velComp) {
//                     dataFAB.setVal(0.0, 2 + velComp);
//                     FABAlgebra::CCaddFC(dataFAB, 2, grids[dit], velFlub[velComp], 0);
//                 }
//             }
//         }

//         this->deallocate(amrData);
//         this->deallocate(amrPressure);
//     }

//     AMRNSLevel::writePlotHeader(a_header, a_filename);
// }
