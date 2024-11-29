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
#ifndef ___LevelLepticSolver_H__INCLUDED___
#define ___LevelLepticSolver_H__INCLUDED___

#include <memory>
#include "LevelSolver.H"
#include "LepticOperator.H"
#include "CFRegion.H"
#include "CornerCopier.H"
#include "BoundaryData.H"
#include "MGSolver.H"

namespace Elliptic {


// -----------------------------------------------------------------------------
// Solves L[phi] = Div[u*] on a single level using the leptic iterative method.
// See: http://dx.doi.org/10.1016/j.jcp.2011.06.022
//
// You are looking at the prototype of the leptic method. This solver has
// been working well for us, but the class itself is not 100% user-friendly.
// There are a lot of parameters to set -- many of which are hard coded. In the
// next version of this code, we hope to create a more versatile solver that is
// a bit more automatic. That is, if you just need to solve equations without
// getting under the hood, we should be able to provide a leptic solver with
// less knobs to turn. -ES
// -----------------------------------------------------------------------------
class LevelLepticSolver : public LevelSolver<LevelData<FArrayBox>>
{
public:
    typedef LevelData<FArrayBox> StateType;

    /// Knobs for you to turn.
    struct Options {
        Options() {}

        Real absTol             = 1.0e-15;
        Real relTol             = 1.0e-15;
        int  maxOrder           = 10;
        int  normType           = 2;
        int  verbosity          = 1;
        Real hang               = 0.01;
        int  maxDivergingOrders = 2;

        MGSolver<LevelData<FArrayBox>>::Options horizOptions;
    };

    /// Constructor -- leaves object unusable.
    LevelLepticSolver();

    /// Destructor
    virtual ~LevelLepticSolver();

    /// Deallocates all memory and brings this solver back to an undefined
    /// state.
    virtual void
    undefine();

    /// \brief Full define -- leave object in a usable state.
    ///  This is an override of the pure virtual LinearSolver function.
    /// \param a_linearOpPtr Must be castable to a LepticOperator.
    virtual void
    define(std::shared_ptr<const LevelOperator<StateType>> a_linearOpPtr) override;

    /// \brief Full define + set options.
    /// \param a_linearOpPtr Must be castable to a LepticOperator.
    /// \param a_opts        The options for the leptic solver.
    virtual void
    define(std::shared_ptr<const LevelOperator<StateType>> a_linearOpPtr,
           const Options&                                  a_opts);

    /// \brief Returns this solver's settings. If called before this object is
    /// defined, these will be the default settings.
    /// \note There is no corresponding setOptions() function becuase the solver
    /// needs to know the options during define() time.
    inline virtual Options
    getOptions() const {
        return m_options;
    }

    /// @brief This is useful when you want to change tolerances or verbosities
    ///  on the fly. Be careful, not all changes will be respected.
    /// @param a_options
    inline virtual void
    setOptions(const Options& a_options)
    {
        m_options = a_options;
        m_horizSolverPtr->setOptions(a_options.horizOptions);
    }

    /// @brief Sets the verbosities of the main, horizontal MG, and horizontal
    /// MG bottom solvers.
    /// @param a_verbs Components are 0 = main, 1 = horiz MG, 2 = horiz MG bottom.
    inline virtual void
    setVerbosities(const std::array<int, 3> a_verbs) {
        m_options.verbosity = a_verbs[0];
        m_options.horizOptions.verbosity = a_verbs[1];
        m_options.horizOptions.bottomOptions.verbosity = a_verbs[2];
    }

    inline virtual std::array<int, 3>
    getVerbosities() const {
        return { m_options.verbosity,
                 m_options.horizOptions.verbosity,
                 m_options.horizOptions.bottomOptions.verbosity };
    }

    /// @brief Returns the residual norms after a solve.
    /// @warning If called before a solve, it may return an empty vector or old
    /// results.
    inline virtual const std::vector<Real>&
    getResNorms() const
    {
        return m_resNorms;
    }

    /// \brief Solve L[phi] = rhs.
    virtual SolverStatus
    solve(StateType&       a_phi,
          const StateType* a_crsePhiPtr,
          const StateType& a_rhs,
          const Real       a_time,
          const bool       a_useHomogBCs,
          const bool       a_setPhiToZero,
          const Real       a_convergenceMetric = -1.0) const override;

protected:
    /// \brief Initializes all parameters with default values.
    virtual void
    setDefaultOptions();

    /// Computes the excess function,
    /// excess = hiNeumBC - loNeumBC - Integral[rhs].
    virtual void
    computeVerticalExcess(BoxLayoutData<FArrayBox>&   a_excess,
                          const LevelData<FArrayBox>& a_rhs,
                          const BoundaryData<Real>&   a_bdryData) const;

    /// \brief Computes the vertical solutions.
    /// \details
    ///  vertRhs will be temporarily altered, but then restored.
    ///  No BCs will be set on vertPhi and this is a non-blocking function.
    virtual void
    verticalLineSolver(LevelData<FArrayBox>&     a_vertPhi,
                       LevelData<FArrayBox>&     a_vertRhs,
                       const BoundaryData<Real>& a_bdryData) const;

    /// Computes the horizontal solutions.
    virtual SolverStatus
    horizontalSolver(LevelData<FArrayBox>&       a_phi,
                     const LevelData<FArrayBox>& a_rhs,
                     const Real                  a_time) const;


    /// Adds a horizontal correction to the vertical solution.
    virtual void
    addHorizontalCorrection(LevelData<FArrayBox>&       a_vertPhi,
                            const LevelData<FArrayBox>& a_horizCor) const;

    virtual void
    addFlatToHoriz(LevelData<FArrayBox>&           a_horiz,
                   const BoxLayoutData<FArrayBox>& a_flat) const;

    virtual void
    copyHorizToFlat(BoxLayoutData<FArrayBox>&  a_flat,
                   const LevelData<FArrayBox>& a_horiz) const;

    /// Brings average[a_phi] to zero.
    static void
    setZeroAvg(LevelData<FArrayBox>& a_phi);


    // Member variables...
    bool                      m_isDefined;
    Options                   m_options;
    mutable std::vector<Real> m_resNorms;

    // Original grid stuff
    RealVect          m_L;
    RealVect          m_dXi;
    RealVect          m_dXiCrse;
    ProblemDomain     m_domain;
    DisjointBoxLayout m_origGrids;

    // Vertical grid stuff (orig grids but not split in vertical)
    DisjointBoxLayout                                     m_grids;
    std::shared_ptr<LevelData<FluxBox>>                   m_JgupPtr;
    Copier                                                m_origToVertCopier;
    Copier                                                m_vertToOrigCopier;
    std::shared_ptr<LevelOperator<LevelData<FArrayBox>>>  m_vertOpPtr;

    // Flat grid stuff
    BoxLayout         m_flatGrids;        // Compatible with m_grids. Lo end is same as domain.
    Vector<DataIndex> m_flatDI;           // Boxes that span the domain...
    Vector<DataIndex> m_flatDIComplement; // ...and the rest.
    BoxLayout         m_shiftedFlatGrids; // Compatible with m_flatGrids. Lo end is 0.

    // Horizontal grid stuff
    bool                                              m_doHorizSolve;
    bool                                              m_horizRemoveAvg;
    ProblemDomain                                     m_horizDomain;
    DisjointBoxLayout                                 m_horizGrids;
    Copier                                            m_shiftedFlatToHorizCopier;
    Copier                                            m_horizToShiftedFlatCopier;
    std::shared_ptr<MGOperator<LevelData<FArrayBox>>> m_horizMGOpPtr;
    std::shared_ptr<MGSolver<LevelData<FArrayBox>>>   m_horizSolverPtr;

#if CH_SPACEDIM == 2
    static constexpr IntVect s_hmask = IntVect(1, 0);
    static constexpr IntVect s_vmask = IntVect(0, 1);
#else
    static constexpr IntVect s_hmask = IntVect(1, 1, 0);
    static constexpr IntVect s_vmask = IntVect(0, 0, 1);
#endif
};


};  // namespace Elliptic
#endif  //!___LevelLepticSolver_H__INCLUDED___
