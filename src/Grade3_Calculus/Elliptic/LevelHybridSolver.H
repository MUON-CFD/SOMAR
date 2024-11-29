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

    enum class SolveMode {
        Undefined,
        MG,
        Leptic,
        Leptic_MG,
    };

    // Knobs for you to turn.
    struct Options {
        Options() {}

        Real absTol         = 1.0e-15;
        Real relTol         = 1.0e-15;
        int  maxSolverSwaps = 10;
        int  normType       = 2;
        int  verbosity      = 1;
        Real hang           = 0.01;

        LevelLepticSolver::Options   lepticOptions;
        MGSolver<StateType>::Options mgOptions;
    };

    /// Constructor -- leaves object unusable.
    LevelHybridSolver();

    /// Destructor
    ~LevelHybridSolver();

    /// Deallocates all memory and brings this solver back to an undefined
    /// state.
    void
    undefine();

    /// \brief Full define -- leave object in a usable state.
    /// \param a_mgOpPtr Must be castable to a LepticOperator.
    void
    define(std::shared_ptr<const MGOperator<StateType>> a_mgOpPtr);

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

    /// \brief Initializes all parameters with default values.
    void
    setDefaultOptions();

    /// \brief Returns this solver's settings. If called before this object is
    /// defined, these will be the default settings.
    /// \note There is no corresponding setOptions() function becuase the solver
    /// needs to know the options during define() time.
    inline Options
    getOptions() const {
        return m_options;
    }


    /// @brief Only some settings will be affected by using this function
    /// instead of redefining the whole solver.
    /// @param a_options
    inline void
    setOptions(const Options& a_options)
    {
        m_options = a_options;
        if (m_mgSolverPtr) {
            m_mgSolverPtr->setOptions(a_options.mgOptions);
        }
        if (m_lepticSolverPtr) {
            m_lepticSolverPtr->setOptions(a_options.lepticOptions);
        }
    }

    /// Returns a copy of our exit status.
    inline SolverStatus
    getSolverStatus() const
    {
        CH_assert(m_isDefined);
        return m_solverStatus;
    }

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
protected:
    /// @brief
    /// @param a_mgOpPtr
    static std::vector<SolveMode>
    computeSolveModes(std::shared_ptr<const MGOperator<StateType>>& a_mgOpPtr,
                      const Options&                                a_options);

    bool                      m_isDefined;
    Options                   m_options;

    mutable std::vector<Real> m_resNorms;
    mutable SolverStatus      m_solverStatus;

    std::shared_ptr<const MGOperator<StateType>> m_mgOpPtr;
    std::vector<SolveMode>                       m_solveModeOrder;
    std::shared_ptr<LevelLepticSolver>           m_lepticSolverPtr;
    std::shared_ptr<MGSolver<StateType>>         m_mgSolverPtr;
};


};  // namespace Elliptic
#endif //!___LevelHybridSolver_H__INCLUDED___
