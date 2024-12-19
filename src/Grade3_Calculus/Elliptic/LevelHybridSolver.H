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
#ifndef ___LevelHybridSolver_H__INCLUDED___
#define ___LevelHybridSolver_H__INCLUDED___

#include "LevelLepticSolver.H"

namespace Elliptic {



class LevelHybridSolver
{
public:
    typedef LevelData<FArrayBox> StateType;

    // -------------------------------------------------------------------------
    /// \name Option setters / getters
    /// \{

    enum class SolveMode {
        Undefined,
        MG,
        Leptic,
        Leptic_MG,
    };

    // Knobs for you to turn.
    struct Options {
        Options() {}

        Real absTol         = -1.0;
        Real relTol         = -1.0;
        int  maxSolverSwaps = -1;
        int  normType       = -1;
        int  verbosity      = -1;
        Real hang           = -1.0;

        LevelLepticSolver::Options   lepticOptions;
        MGSolver<StateType>::Options mgOptions;
    };

    /// @brief Creates an Options struct set to defaults as defined by the input
    /// file.
    static Options
    getDefaultOptions();

    /// @brief Creates an Options struct set to silently perform a single FMG iter.
    static Options
    getQuickAndDirtyOptions();

    /// Returns the current solver's options.
    inline virtual const Options&
    getOptions() const {
        return m_options;
    }

    /// @brief This allows you to change absTol, relTol, etc on-the-fly.
    /// This is useful when you want to use an already-defined solver to do a
    /// quick and dirty estimate. Just raise the tols or set maxIters to 1.
    ///
    /// @param a_opt The new options. All MG maxDepth changes will be ignored.
    ///              Your changes to lepticOptions.verbosity will be ignored as
    ///              well.
    ///
    /// @details The reason we can't change maxDepth after the solver is defined
    /// is because we'd need to reallocate and redefine all of the solvers in
    /// the MG hierarchy - an expensive proposition if all you want is a quick
    /// and dirty solution. It's probably faster (and definitely easier) to just
    /// use the MG hierarchy that's available.
    ///
    /// @warning This can only be called on an already-defined solver.
    virtual void
    modifyOptionsExceptMaxDepth(const Options& a_opt);
    /// \}


    // -------------------------------------------------------------------------
    /// \name Constructors / destructors
    /// \{

    /// Constructor -- leaves object unusable.
    LevelHybridSolver();

    /// Destructor
    ~LevelHybridSolver();

    /// \brief Full define + set options.
    /// \param a_mgOpPtr Must be castable to a LepticOperator.
    /// \param a_opts
    void
    define(std::shared_ptr<const MGOperator<StateType>> a_mgOpPtr,
           const Options&                               a_opts);

    /// @brief Full define -- leave object in a usable state.
    /// @param a_lepticSolverPtr
    /// @param a_mgSolverPtr
    /// \param a_mgOpPtr Must be castable to a LepticOperator.
    /// \param a_opts    Only the hybrid solver options will be used,
    ///                  not the leptic or MG options.
    void
    define(std::shared_ptr<LevelLepticSolver>           a_lepticSolverPtr,
           std::shared_ptr<MGSolver<StateType>>         a_mgSolverPtr,
           std::shared_ptr<const MGOperator<StateType>> a_mgOpPtr,
           const Options&                               a_opts);

    /// Simple check.
    inline virtual bool
    isDefined() const {
        return m_isDefined;
    }

    /// Free memory. Leaves object unusable.
    void
    clear();


    // -------------------------------------------------------------------------
    /// \name Solvers
    /// \{

    /// \brief Solve L[phi] = rhs.
    SolverStatus
    solve(StateType&       a_phi,
          const StateType* a_crsePhiPtr,
          const StateType& a_rhs,
          const Real       a_time,
          const bool       a_useHomogBCs,
          const bool       a_setPhiToZero,
          const Real       a_convergenceMetric = -1.0) const;

    /// \brief Solve residual eq L[cor] = res, assuming homog BCs.
    SolverStatus
    solveResidualEq(StateType&       a_cor,
                    const StateType& a_res,
                    const Real       a_time,
                    const Real       a_convergenceMetric = -1.0) const;
    /// \}


    /// Returns a copy of our exit status.
    inline SolverStatus
    getSolverStatus() const
    {
        CH_assert(m_isDefined);
        return m_solverStatus;
    }

protected:
    /// @brief
    /// @param a_mgOpPtr
    static SolveMode
    computeSolveMode(std::shared_ptr<const MGOperator<StateType>>& a_mgOpPtr,
                     const Options&                                a_options);

    bool                      m_isDefined;
    Options                   m_options;

    mutable std::vector<Real> m_resNorms;
    mutable SolverStatus      m_solverStatus;

    std::shared_ptr<const MGOperator<StateType>> m_mgOpPtr;
    SolveMode                                    m_solveMode;
    std::shared_ptr<LevelLepticSolver>           m_lepticSolverPtr;
    std::shared_ptr<MGSolver<StateType>>         m_mgSolverPtr;
};


};  // namespace Elliptic
#endif //!___LevelHybridSolver_H__INCLUDED___
