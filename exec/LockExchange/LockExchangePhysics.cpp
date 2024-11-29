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
 *  https://github.com/MUON-CFD/somar.
 ******************************************************************************/
#include "LockExchangePhysics.H"
#include "BCTools.H"
#include "SetValLevel.H"
#include "BoxIterator.H"
#include "Debug.H"
#include "Zipper.H"
#include "Comm.H"

#include "Convert.H"
#include "AnisotropicRefinementTools.H"


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
    for (auto [s] : zip(a_state.scalars)) {
        s.setVal(1.0,LAMBDA);
    }

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

