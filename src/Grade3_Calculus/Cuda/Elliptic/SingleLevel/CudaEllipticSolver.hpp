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
#pragma once

#include "FunctionTimer.hpp"
#include "isHilbert.hpp"
#include "isLinearOperatorOverHilbert.hpp"
#include <cusp/krylov/bicgstab.h>
#include <cusp/krylov/gmres.h>
#include <cusp/precond/aggregation/smoothed_aggregation.h>
#include <cusp/monitor.h>
namespace Elliptic
{
namespace SingleLevel
{
/**
*  The generic EllipticSolver interface.
*/
template <typename HilbertSpace, typename LinearOperator>
class EllipticSolver
{
  static_assert(Elliptic::traits::isLinearOperatorOverHilbert<
                LinearOperator, HilbertSpace>::value);

public:
  EllipticSolver();

  virtual ~EllipticSolver();

  virtual void define(const std::shared_ptr<LinearOperator> a_opPtr);

  virtual const LinearOperator &getOp() const { return *m_opPtr; } /*!< returns the LinearOperator associated with this solver*/

/// The generic interface for a solve.
  virtual int solve(HilbertSpace &a_phi, const HilbertSpace &a_rhs,
                    const Real a_time, const bool a_useHomogBCs,
                    const bool a_setPhiToZero) const = 0;

  // I'm doing this without an enum so that child classes can add to
  // the list if needed.
  static constexpr int UNDEFINED = -1;
  static constexpr int DIVERGED = 0;
  static constexpr int CONVERGED = 1;
  static constexpr int SINGULAR = 2;
  static constexpr int MAXITERS = 3;
  static constexpr int HANG = 4;

  /// Returns the exit status.
  virtual int getExitStatus() const { return m_exitStatus; }

  /// Convert the exit status to a readable message.
  virtual const char *getExitMessage(const int a_exitStatus) const;

protected:
  /// Set the exit status.
  virtual void setExitStatus(const int a_exitStatus) const
  {
    m_exitStatus = a_exitStatus;
  }

  std::shared_ptr<LinearOperator> m_opPtr;
  mutable int m_exitStatus = UNDEFINED;

private:
  EllipticSolver &operator=(const EllipticSolver &that){};
};

/**
*This does nothing! It's useful for debugging V-Cycles.
*/
template <class T, class L>
class CudaNullSolver : public EllipticSolver<T, L>
{
public:
  CudaNullSolver() {}

  virtual ~CudaNullSolver() {}

  virtual int solve(T &a_phi, const T &a_rhs, const Real a_time,
                    const bool a_useHomogBCs, const bool a_setPhiToZero) const
  {
    return this->getExitStatus();
  }
};

/**
* BiCGStab linear solver.
* Copied almost directly from Chombo.
*/
template <class T, class L>
class BiCGStabSolver : public EllipticSolver<T, L>
{

public:
  /// Knobs for you to turn.
  struct Options
  {
    Options() {} // A bug in the C++11 standard requires this.

    Real absTol = 1.0e-6; //!< Stop if |res| < absTol.
    Real relTol = 1.0e-6; //!< Stop if |res| / |init res| < relTol.
    Real small = 1.0e-30;
    Real hang = 1.0e-7;
    Real convergenceMetric = -1.0; //!< If < 0, use initial residual norm.
    int maxIters = 8000;
    int maxRestarts = 5;
    int normType = 2; //!< We only use l2 norm here.
    int verbosity = 3;
    int numSmoothPrecond = 2;
  };

  /// Constructor. You must still call define() to supply an operator.
  BiCGStabSolver(Options a_opt = Options());

  /// Destructor
  virtual ~BiCGStabSolver();

  /// Retrieves the solver's options.
  virtual const Options &getOptions() const { return m_opt; }

  /// Replace this solver's options.
  virtual void setOptions(const Options &a_opt) { m_opt = a_opt; }

  /// This allows you to manipulate individual options.
  virtual Options &options() { return m_opt; }

  /// Solve
  virtual int solve(T &a_phi, const T &a_rhs, const Real a_time,
                    const bool a_useHomogBCs, const bool a_setPhiToZero) const;

protected:
  Options m_opt;

private:
  BiCGStabSolver &operator=(const BiCGStabSolver &that){};
};

/**
 * This uses the cusp BiCGStab solver. It inherits from BiCGStabSolver.
 * Note that it really only makes sense if using the global solver, since
 * the cusp solver does not know how to exchange across mpi solvers.
 */
template <class Hilbert, class LinearOperator>
class CuspBiCGStabSolver : public BiCGStabSolver<Hilbert, LinearOperator>
{
  private:
  //typedef typename  cusp::precond::aggregation::smoothed_aggregation<int, Real, cusp::device_memory> PreCond;
  typedef typename  cusp::precond::diagonal<Real, cusp::device_memory> PreCond; //!< There
public:
  CuspBiCGStabSolver() : BiCGStabSolver<Hilbert, LinearOperator>() {}

  virtual ~CuspBiCGStabSolver()
  {
    delete m_Precond;
  }
  /**
   * In the definition, we need to create the Preconditioner.
   * There are several possibilities. Here we use the standard diagonal
   * preconditioner.
   */
  virtual void define(const std::shared_ptr<LinearOperator> a_OpPtr)
  {
    CH_assert(a_OpPtr.get());

    this->m_opPtr = a_OpPtr;
    m_Precond = new PreCond(*(a_OpPtr->getOperatorMatrixPtr()));


  }
  /**
   * Solves the system. We use the monitor class from cusp to control convergence.
   * The only return values are 1 if converged or 0 if not.
   */
  virtual int solve(Hilbert &a_phi, const Hilbert &a_rhs, const Real a_time, const bool a_useHomogBCs, const bool a_setPhiToZero) const
  {


    cusp::monitor<Real> monitor(a_rhs(), this->m_opt.maxIters, this->m_opt.relTol, this->m_opt.absTol, false);
    cusp::krylov::bicgstab(*(this->m_opPtr->getOperatorMatrixPtr()), a_phi(), a_rhs(), monitor, *m_Precond);
    //cusp::krylov::gmres(*(this->m_opPtr->getOperatorMatrixPtr()), a_phi(), RHS(), 50, monitor, *m_Precond);
    //monitor.print();
    return monitor.converged() ? 1:0;
  }

private:
  CuspBiCGStabSolver &operator=(const CuspBiCGStabSolver &that){};

  PreCond *m_Precond;


};
#define H9929fb0cc38edce433e3b0842ce1c397
#include "CudaEllipticSolverImpl.hpp"

#undef H9929fb0cc38edce433e3b0842ce1c397
} // namespace SingleLevel
} // namespace Elliptic
