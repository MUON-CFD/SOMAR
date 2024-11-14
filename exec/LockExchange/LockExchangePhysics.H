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
#ifndef ___LockExchangePhysics_H__INCLUDED___
#define ___LockExchangePhysics_H__INCLUDED___

#include "AMRNSLevel.H"


class LockExchangePhysics: public AMRNSLevel
{
public:
    // Constructor
    LockExchangePhysics();

    // Virtual destructor for good measure.
    virtual
    ~LockExchangePhysics();

    // Custom scalars ----------------------------------------------------------
    // Custom scalars will automatically be advected and, if requested,
    // diffused. If your needs are more specialized than that, then it
    // is up to you to define a forcing function (see below).

    // This tells the Navier-Stokes class how many additional
    // (temperature and salinity not counted) scalars to create. This is inline
    // so that is is just as fast as reading from a const int.
    virtual int
    numScalars() const override
    {
        return SpaceDim + 1;
    }

    // Scalar identifiers. This is not absolutely necessary, but it makes some
    // bits of code easier to read. I set Y_TRACER in 2D to avoid "variable
    // not declared" compilation errors.
#if CH_SPACEDIM == 2
    enum { Y_TRACER = -999, X_TRACER = 0, Z_TRACER = 1, LAMBDA = 2 };
#else
    enum { X_TRACER = 0, Y_TRACER = 1, Z_TRACER = 2, LAMBDA = 3 };
#endif

    // Returns the names of the scalars. This will be used to put scalar names
    // in the plot files.
    virtual std::string
    getScalarName(const int a_comp) const override;


    // ICs ---------------------------------------------------------------------
    // Set the ICs on all fields at once. This includes the velocity, buoyancy,
    //  eddy viscosity and diffusivity, and all scalars. Note that I did not
    //  say pressure. Upon entry, the pressure will already be set to it's
    //  initial guess of 0. If you have a better estimate, feel free to provide
    //  it. We will use it as an initial guess for the projector.
    // By the time this function is called, *m_levGeoPtr will be ready for use.
    virtual void
    setICs(State& a_state) override;

    // Returns maximum stable time step at t = 0.
    virtual Real
    computeInitialDt() override {return 0.005;}


    // BCs ---------------------------------------------------------------------
    // Make these functions very fast. They will be called many times
    // per timestep, potenially by AMRMG solvers.

    virtual BCTools::BCFunction*
    createScalarsPhysBC(const int             a_dir,
                        const Side::LoHiSide& a_side) const override;

    // // Set physical BCs on all scalars.
    // virtual void
    // setScalarPhysBC(LevelData<FArrayBox>& a_state,
    //                 const Real            a_time,
    //                 const bool            a_homogBCs) const override;

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
    // T and S must be the total temperature and salinity, not just the
    // deviation from the background. Likewise, this function returns the total
    // buoyancy.
    //
    // By default, this uses a simple linear equation of state:
    //  b = -\alpha*g*T + \beta*g*S
    // with \alpha = 2.0e-4 K^-1, \beta = 8.0e-4 psu^-1, and g = 9.81 m/s^2.
    virtual void
    equationOfState(FArrayBox&       a_bFAB,
                    const FArrayBox& a_TFAB,
                    const FArrayBox& /*a_SFAB*/,
                    const FArrayBox& /*a_zFAB*/) const override
    {
        a_bFAB.copy(a_TFAB);
    }

    // //
    // virtual void
    // tagCellsInit(IntVectSet& a_tags);

    // //
    // virtual void
    // tagCells(IntVectSet& a_tags);

    ///
    virtual void
    constructIB(Real a_time = 0.0) override;

    // ///
    // virtual void
    // writePlotHeader(HeaderData& a_header, const std::string& a_filename) const;
};

#endif //!___LockExchangePhysics_H__INCLUDED___
