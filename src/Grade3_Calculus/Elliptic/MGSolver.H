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
#ifndef ___ElllipticMGSolver_H__INCLUDED___
#define ___ElllipticMGSolver_H__INCLUDED___

#include "MGOperator.H"
#include "LevelSolver.H"

namespace Elliptic {


/**
 * @interface MGSolver
 * @brief     The multi-grid elliptic solver.
 * @date      2019/11/26
 * @author    Edward.Santilli@jefferson.edu
 * @note
 *  Currently, this only works when solving L[phi]=rhs where phi and rhs are
 *  LevelData<FArrayBox>es.
 */
template <class T>
class MGSolver
{
public:
    typedef T                                  StateType;
    typedef typename StateTraits<T>::GridsType GridsType;
    typedef MGOperator<T>                      MGOpType;
    typedef BiCGStabSolver<T>                  BottomSolverType;

    /// Knobs for you to turn.
    struct Options {
        /// A bug in the C++11 standard requires this,
        /// even though it does nothing.
        Options() {}

        /// Stop if |res| < absTol
        Real absTol            = 1.0e-12;
        /// Stop if |res| / |init res| < relTol
        Real relTol            = 1.0e-12;
        Real convergenceMetric = -1.0;
        int  numSmoothDown     = 4;
        int  numSmoothUp       = 4;
        int  numSmoothBottom   = 2;
        int  numSmoothPrecond  = 2;
        int  prolongOrder      = 1;
        int  prolongOrderFMG   = 3;
        int  numSmoothUpFMG    = 0;
        /// -1 = coarsen as much as possible
        int  maxDepth          = -1;
        int  numCycles         = 1;
        int  maxIters          = 10;
        Real hang              = 0.01;
        /**
         * @brief 0 = inf-norm, 1 = 1-norm, etc...
         * @detail
         *  I suggest using the 2-norm (defualt) so that norms across different
         *  MG depths are comparable.
         */
        int  normType          = 2;
        int  verbosity         = 4;

        typename BottomSolverType::Options bottomOptions;
    };

    /// Default constructor. Leaves object unusable.
    MGSolver();

    /// Destructor
    virtual ~MGSolver();

    /// @brief Full define. Leaves object useable.
    /// @param a_topOpPtr        Must be well-defined at depth = 0.
    /// @param a_opt             The MG options.
    /// @param a_mgRefSchedule   You can provide your own refRatio at each depth.
    ///                          Otherwise, SemicoarseningStrategy will be used.
    /// @param a_useBottomSolver Default is to use BiCGStab. Otherwise, only
    ///                          relaxation will be used at the bottom.
    virtual void
    define(std::shared_ptr<const MGOpType> a_topOpPtr,
           Options                         a_opt,
           Vector<IntVect>                 a_mgRefSchedule = Vector<IntVect>(),
           const bool                      a_useBottomSolver = true);

    /// @brief Full define. Leaves object useable.
    /// @param a_topOpPtr        Must be well-defined at depth = 0.
    /// @param a_opt             The MG options.
    /// @param a_mgRefSchedule   You can provide your own refRatio at each depth.
    ///                          Otherwise, SemicoarseningStrategy will be used.
    /// @param a_useBottomSolver Default is to use BiCGStab. Otherwise, only
    ///                          relaxation will be used at the bottom.
    virtual void
    define(const MGOpType& a_topOp,
           Options         a_opt,
           Vector<IntVect> a_mgRefSchedule   = Vector<IntVect>(),
           const bool      a_useBottomSolver = true);

    /// @brief Full define. Leaves object useable. This version does not set
    ///        the MG options, just in case you already did that.
    /// @param a_topOpPtr        Must be well-defined at depth = 0.
    /// @param a_mgRefSchedule   You can provide your own refRatio at each depth.
    ///                          Otherwise, SemicoarseningStrategy will be used.
    /// @param a_useBottomSolver Default is to use BiCGStab. Otherwise, only
    ///                          relaxation will be used at the bottom.
    virtual void
    define(std::shared_ptr<const MGOpType> a_topOpPtr,
           Vector<IntVect>                 a_mgRefSchedule = Vector<IntVect>(),
           const bool                      a_useBottomSolver = true);

    /// @brief Full define. Leaves object useable. This version does not set
    ///        the MG options, just in case you already did that.
    /// @param a_topOpPtr        Must be well-defined at depth = 0.
    /// @param a_mgRefSchedule   You can provide your own refRatio at each depth.
    ///                          Otherwise, SemicoarseningStrategy will be used.
    /// @param a_useBottomSolver Default is to use BiCGStab. Otherwise, only
    ///                          relaxation will be used at the bottom.
    virtual void
    define(const MGOpType& a_topOp,
           Vector<IntVect> a_mgRefSchedule   = Vector<IntVect>(),
           const bool      a_useBottomSolver = true);

    /// Simple check.
    virtual bool
    isDefined() const;

    /// Free's memory. Leaves object unusable.
    virtual void
    clear();

    /// Returns a copy of our exit status.
    virtual SolverStatus
    getSolverStatus() const
    {
        return m_solverStatus;
    }

    /// \brief Initializes all parameters with default values.
    virtual void
    setDefaultOptions();

    /// Returns the current solver's options.
    virtual const Options&
    getOptions() const;

    /// I'm giving you access to this so that you can change absTol, relTol, etc
    /// on-the-fly. This is useful when you want to use an already-defined
    /// MGSolver to do a quick and dirty estimate. Just increase the tolerances
    /// and set maxIters to 1.
    /// NOTE: You cannot use this to change maxDepth.
    virtual void
    setOptions(const Options& a_opt);

    /**
     *
     * @brief      Solves L[phi] = rhs using the most appropriate methods.
     *
     * @param[inout] a_phi                On entry, the initial guess.
     *                                    On exit, the solution.
     * @param[in]    a_crsePhiPtr         Used to set CFI BCs, if needed.
     *                                    Can be null if not needed.
     * @param[in]    a_rhs
     * @param[in]    a_time               Used to set BCs.
     * @param[in]    a_useHomogBCs
     * @param[in]    a_setPhiToZero       Should we set initial guess to zero?
     * @param[in]    a_convergenceMetric  If supplied, this will be used as the
     *                                    initial residual when computing
     *                                    relative residuals. This setting takes
     *                                    precedence over
     *                                    Options::convergenceMetric.
     *
     * @return     The SolverStatus object.
     *             This tells the caller if we converged.
     */
    virtual SolverStatus
    solve(StateType&       a_phi,
          const StateType* a_crsePhiPtr,
          const StateType& a_rhs,
          const Real       a_time,
          const bool       a_useHomogBCs       = false,
          const bool       a_setPhiToZero      = false,
          const Real       a_convergenceMetric = -1.0);

    /**
     * @brief      Solves L[phi] = rhs using V-Cycles.
     *
     * @param[inout] a_phi                On entry, the initial guess.
     *                                    On exit, the solution.
     * @param[in]    a_crsePhiPtr         Used to set CFI BCs, if needed.
     *                                    Can be null if not needed.
     * @param[in]    a_rhs
     * @param[in]    a_time               Used to set BCs.
     * @param[in]    a_useHomogBCs
     * @param[in]    a_setPhiToZero       Should we set initial guess to zero?
     * @param[in]    a_convergenceMetric  If supplied, this will be used as the
     *                                    initial residual when computing
     *                                    relative residuals. This setting takes
     *                                    precedence over
     *                                    Options::convergenceMetric.
     *
     * @return     The SolverStatus object.
     *             This tells the caller if we converged.
     */
    virtual SolverStatus
    vCycle(StateType&       a_phi,
           const StateType* a_crsePhiPtr,
           const StateType& a_rhs,
           const Real       a_time,
           const bool       a_useHomogBCs       = false,
           const bool       a_setPhiToZero      = false,
           const Real       a_convergenceMetric = -1.0);

    /**
     * @brief      Solves L[phi] = rhs using Full Multigrid.
     *
     * @param[inout] a_phi                On entry, the initial guess.
     *                                    On exit, the solution.
     * @param[in]    a_crsePhiPtr         Used to set CFI BCs, if needed.
     *                                    Can be null if not needed.
     * @param[in]    a_rhs
     * @param[in]    a_time               Used to set BCs.
     * @param[in]    a_useHomogBCs
     * @param[in]    a_setPhiToZero       Should we set initial guess to zero?
     * @param[in]    a_convergenceMetric  If supplied, this will be used as the
     *                                    initial residual when computing
     *                                    relative residuals. This setting takes
     *                                    precedence over
     *                                    Options::convergenceMetric.
     *
     * @return     The SolverStatus object.
     *             This tells the caller if we converged.
     */
    virtual SolverStatus
    fmg(StateType&       a_phi,
        const StateType* a_crsePhiPtr,
        const StateType& a_rhs,
        const Real       a_time,
        const bool       a_useHomogBCs       = false,
        const bool       a_setPhiToZero      = false,
        const Real       a_convergenceMetric = -1.0);

    /// Solves L[cor] = res with homogBCs at a specified depth.
    virtual void
    vCycle_residualEq(StateType&       a_cor,
                      const StateType& a_res,
                      const Real       a_time,
                      const int        a_depth);

    /// Solves L[cor] = res with homogBCs at a specified depth.
    virtual void
    fmg_residualEq(StateType&       a_cor,
                   const StateType& a_res,
                   const Real       a_time,
                   const int        a_depth);

protected:
    /// pout formatting utility.
    virtual void
    indent(const int a_depth = -1) const;

    /// pout formatting utility.
    virtual void
    unindent(const int a_depth = -1) const;

    /// References our exit status.
    virtual SolverStatus&
    getSolverStatusRef() const
    {
        return m_solverStatus;
    }

    // Member variables...
    Vector<IntVect>                   m_refSchedule;
    Vector<RealVect>                  m_dXi;

    // Use smart pointers here because we own them and we want them all deleted.
    Vector<shared_ptr<MGOpType>>      m_opPtrs;
    unique_ptr<BottomSolverType>      m_bottomSolverPtr;

    Options                           m_opt;
    mutable SolverStatus              m_solverStatus;
};


}; // namespace Elliptic

#define H955d3c5f3254f48107821d4c7c238980
#include "MGSolverI.H"
#undef H955d3c5f3254f48107821d4c7c238980

#endif //!___ElllipticMGSolver_H__INCLUDED___


