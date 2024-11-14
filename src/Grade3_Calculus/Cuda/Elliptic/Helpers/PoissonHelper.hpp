/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2019
 *    Jefferson University and
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
#pragma once
#include <chrono>

#include "CudaMGSolver.hpp"
#include "LevelDataHilbert.hpp"
#include "LevelGeometry.H"
#include <cusp/eigen/spectral_radius.h>
// #include <cusp/eigen/lobpcg.h>
// #include <cusp/eigen/lanczos.h>

#include "AMRNSLevel.H"
namespace Elliptic
{
namespace Helpers
{

template <typename Container, typename LinearOperator, typename BC>
class PoissonHelper
{
public:
  typedef Elliptic::SingleLevel::BiCGStabSolver<Container, LinearOperator>
      SingleLevelSolverType;
  typedef Elliptic::MultiGrid::MGSolver<Container, LinearOperator,
                                        SingleLevelSolverType, BC>
      SolverType;

  PoissonHelper(const LevelGeometry *thisLevelGeometry, const int level)
      : m_levGeoPtr(thisLevelGeometry), m_isDefined(true), m_level(level)
  {
    this->setupPoissonSolver();
  }
  virtual ~PoissonHelper() { this->cleanupPoissonSolver(); }

  SolverType *getSolverPtr()
  {
    if (!m_isDefined)
      MayDay::Error(
          "AMRNSLevel::PoissonHelpers. Solver pointer requested "
          "before solver is defined.");
    return &m_solver;
  }

  void project(LevelData<FluxBox> &vel, LevelData<FArrayBox> &phi,
               LevelData<FluxBox> &gradPhi) 
  {
    Container RHS, LHS;

    this->setupContainers(RHS, LHS);
    this->m_solver.getOpPtr()->divergence(vel, RHS);

    // the quality of convergence is really sensitive to how close to 0 the
    // average is.
    int counter = 0;
    while (abs(RHS.removeAverage(1. / phi.getBoxes().numCells())) > 1e-20 &&
           counter++ < 4)
    {
    }

    this->LDToContainer(phi, LHS);
    this->solve(LHS, RHS);
    this->ContainerToLD(LHS, phi);
    phi.exchange();
    this->getSolverPtr()->getOpPtr()->gradient(LHS, gradPhi);
    DataIterator dit = this->m_levGeoPtr->getBoxes().dataIterator();
    for (dit.reset(); dit.ok(); ++dit)
    {
      D_TERM(vel[dit][0].plus(gradPhi[dit][0], -1.0);
             , vel[dit][1].plus(gradPhi[dit][1], -1.0);
             , vel[dit][2].plus(gradPhi[dit][2], -1.0);)
    }
  }
  auto solve(Container &LHS, const Container &RHS)
      -> decltype(std::declval<Container>().norm()) const
  {
    auto start = std::chrono::high_resolution_clock::now();
    auto res = this->m_solver.solve(LHS, RHS, 0.0, true);
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    this->m_MFlops = ((double)this->m_OperationCount / elapsed.count() / 1.e6) * res.size();
    return res.back();
  }

  void LDToContainer(const LevelData<FArrayBox> &LD, Container &C) const
  {
    this->m_solver.getOpPtr()->LDToContainer(LD, C);
  }

  void ContainerToLD(const Container &C, LevelData<FArrayBox> &LD) const
  {
    this->m_solver.getOpPtr()->ContainerToLD(C, LD);
  }

  auto testSolver(Container &ERR)
      -> decltype(std::declval<Container>().norm()) const
  {
    Container RHS, LHS;
    this->setupContainers(RHS, LHS, ERR);
    LHS.fillRandom(-1., 1.);
    this->m_solver.getOpPtr()->applyOp(RHS, LHS, 0.0, true);
    LHS.setToZero();

    auto res = this->solve(LHS, RHS);
    this->m_solver.getOpPtr()->residual(ERR, LHS, RHS, 0.0, true);
    ERR.scale(1. / RHS.norm());
    return res;
  }
  auto testSolver() -> decltype(std::declval<Container>().norm()) const
  {
    Container RHS, LHS, ERR;
    this->setupContainers(RHS, LHS, ERR);
    LHS.fillRandom(-1., 1.);
    this->m_solver.getOpPtr()->applyOp(RHS, LHS, 0.0, true);
    LHS.setToZero();

    return this->solve(LHS, RHS);

    // this->m_solver.getOpPtr()->residual(ERR, LHS, RHS, 0.0, true);
    // ERR.scale(1. / RHS.norm());
    // return ERR.norm();
  }
  auto testProjector() 
  {
    Container RHS, LHS, ERR;
    this->setupContainers(RHS, LHS, ERR);
    LevelData<FArrayBox> div, phi;
    LevelData<FluxBox> vel(this->m_levGeoPtr->getBoxes(), 1);
    LevelData<FluxBox> gradPhi(this->m_levGeoPtr->getBoxes(), 1);
    // first we generate a random velocity field applying the gradient operator to
    // a random potential field
    LHS.fillRandom(-1., 1.);
    this->setupContainers(div, phi);
    this->m_solver.getOpPtr()->gradient(LHS, gradPhi);
    DataIterator dit = this->m_levGeoPtr->getBoxes().dataIterator();
    for (dit.reset(); dit.ok(); ++dit)
    {
      D_TERM(vel[dit][0].copy(gradPhi[dit][0]);
             , vel[dit][1].copy(gradPhi[dit][1]);
             , vel[dit][2].copy(gradPhi[dit][2]);)
    }
    // then we project 
    this->m_solver.getOpPtr()->divergence(vel, RHS);
    

    LHS.setToZero();
    auto res=this->solve(LHS, RHS);
    this->m_solver.getOpPtr()->gradient(LHS, gradPhi);

    for (dit.reset(); dit.ok(); ++dit)
    {
      D_TERM(vel[dit][0].plus(gradPhi[dit][0], -1.0);
             , vel[dit][1].plus(gradPhi[dit][1], -1.0);
             , vel[dit][2].plus(gradPhi[dit][2], -1.0);)
    }
    // and finally we calculate the residual 
    this->m_solver.getOpPtr()->divergence(vel, ERR);
    auto res_divgrad=ERR.norm()/RHS.norm();
    IO::tout(0) << " Residual of solver " << res << ". Residual of divgrad " << res_divgrad << endl;
    return res_divgrad;
  }

  void findLowestMode(const std::shared_ptr<LinearOperator> L, LevelData<FArrayBox> & EigenVector)
  { // uses Python to calculate the first two smallest eigenpairs. 
    static int c=0;
    std::vector<int> cols, row_offsets;
    std::vector<Real> Data, eig;
    auto &Matrix=*(L->getOperatorMatrixPtr());
    cols.resize(Matrix.column_indices.size());
    row_offsets.resize(Matrix.row_offsets.size());
    Data.resize(Matrix.values.size());
    eig.resize(row_offsets.size()-1);

    thrust::copy(Matrix.column_indices.begin(), Matrix.column_indices.end(), cols.begin());
    thrust::copy(Matrix.row_offsets.begin(), Matrix.row_offsets.end(), row_offsets.begin());
    thrust::copy(Matrix.values.begin(), Matrix.values.end(), Data.begin());

    Python.PythonFunction("ElliKitWrapper", "LowestMode", cols, row_offsets, Data, eig);
    Container EIG(eig.size());
    thrust::copy(eig.begin(), eig.end(), EIG().begin());
    
    L->ContainerToLD(EIG, EigenVector);
    
  }
  // void findLowestMode()
  // {
  //   cusp::array1d<Real, cusp::device_memory> S(5,0.);
  //   cusp::array2d<Real, cusp::device_memory, cusp::column_major> V;

  //   // initialize Lanczos Options
  //   cusp::eigen::lanczos_options<Real> options;

  //   options.tol=1e-6;
  //   options.maxIter =100;
  //   options.verbose=true;
  //   options.computeEigVecs = true;
  //   options.reorth = cusp::eigen::Full;
  //   // compute eigen pairs
  //   cusp::eigen::lanczos(*(this->m_solver.getOpPtr()->getOperatorMatrixPtr()), S, V, options);
  //   for (int i=0; i<S.size(); ++i)
  //   IO::tout(0) << "Largest EigenValue # " << i << " is " << S[i] << endl;

  // }
  template <typename T_,
            typename std::enable_if<
                std::is_same<T_, Elliptic::Helpers::NodeData<
                                     Real, typename T_::MyMemory>>::value,
                int>::type = 0>
  void setupContainers(T_ &RHS, T_ &LHS, T_ &ERR) const
  {
    RHS.resize(this->m_levGeoPtr->getBoxes().numCells());
    LHS.resize(this->m_levGeoPtr->getBoxes().numCells());
    ERR.resize(this->m_levGeoPtr->getBoxes().numCells());
  }

  template <typename T_, typename S_ = typename std::enable_if<
                             std::is_polymorphic<T_>::value, T_>::type>
  void setupContainers(T_ &RHS, T_ &LHS, T_ &ERR) const
  {
    if (RHS.isDefined())
      RHS.clear();
    if (LHS.isDefined())
      LHS.clear();
    if (ERR.isDefined())
      ERR.clear();
    // Grids for the lhs are defined peeling off the first layer next to a
    // physical boundary
    const DisjointBoxLayout &rhsGrids = this->m_levGeoPtr->getBoxes();
    DisjointBoxLayout lhsGrids;
    lhsGrids.deepCopy(rhsGrids);
    lhsGrids.close();

    ResizeGrids<BC> RG(rhsGrids, this->m_levGeoPtr->getDomain());
    lhsGrids.transform(RG); //
    RHS.define(rhsGrids, 1);
    LHS.define(lhsGrids, 1, IntVect::Unit);
    ERR.define(rhsGrids, 1);
  }
  template <typename T_,
            typename std::enable_if<
                std::is_same<T_, Elliptic::Helpers::NodeData<
                                     Real, typename T_::MyMemory>>::value,
                int>::type = 0>
  void setupContainers(T_ &RHS, T_ &LHS) const
  {
    RHS.resize(this->m_levGeoPtr->getBoxes().numCells());
    LHS.resize(this->m_levGeoPtr->getBoxes().numCells());
  }
  template <typename T_, typename S_ = typename std::enable_if<
                             std::is_polymorphic<T_>::value, T_>::type>
  void setupContainers(T_ &RHS, T_ &LHS) const
  {
    if (RHS.isDefined())
      RHS.clear();
    if (LHS.isDefined())
      LHS.clear();

    const DisjointBoxLayout &rhsGrids = this->m_levGeoPtr->getBoxes();
    DisjointBoxLayout lhsGrids;
    lhsGrids.deepCopy(rhsGrids);
    lhsGrids.close();

    ResizeGrids<BC> RG(rhsGrids, this->m_levGeoPtr->getDomain());
    lhsGrids.transform(RG); //
    RHS.define(rhsGrids, 1);
    LHS.define(lhsGrids, 1, IntVect::Unit);
  }

  double getMFlops()
  {
    return this->m_MFlops;
  }

private:
  PoissonHelper() : m_isDefined(false) {}
  PoissonHelper &operator=(const PoissonHelper &that){};
  Vector<LevelGeometry *> m_levGeoPtrsPoisson;
  const LevelGeometry *m_levGeoPtr;
  bool m_isDefined = false;
  int m_level = 0;

public:
  SolverType m_solver;
  private:
  mutable unsigned int m_OperationCount = 0;
  mutable float m_MFlops;

  void setupPoissonSolver()
  {
    // Create structures to host data and solver
    static int c=0;
    // Tailor solver options
    typename SolverType::Options options;

    options.maxDepth = -1;
    options.absTol = 1.0e-12;
    options.relTol = 1.0e-12;
    options.maxIters = 10;
    options.numSmoothDown = 4;
    options.numSmoothUp = 4;
    options.numSmoothPrecond = 2;
    options.numSmoothBottom = 2;
    options.verbosity = 0;

    const LevelGeometry &levGeo = *(m_levGeoPtr);
    BC phiBC(&levGeo);

    // definition of the hyerarchy of operators;
    std::vector<std::shared_ptr<LinearOperator>> OpsPtrs;

    auto refSchedule =
        Elliptic::Helpers::SemicoarseningStrategy(levGeo.getDomainLength())
            .createMGRefSchedule(levGeo.getBoxes(),
                                 LinearOperator::minBoxSize(),
                                 options.maxDepth);

    // Create the levGeo hierarchy. We will own all of these.
    this->m_levGeoPtrsPoisson = levGeo.createMGLevGeos(refSchedule);
    std::vector<IntVect> refScheduleVector(refSchedule.size());
    for (int i = 0; i < m_levGeoPtrsPoisson.size(); ++i)
    {
      refScheduleVector.at(i) = refSchedule[i];
      OpsPtrs.emplace_back(
          new LinearOperator(this->m_levGeoPtrsPoisson[i], &phiBC, this->m_level));
      this->m_OperationCount +=
          (14 * (options.numSmoothDown + options.numSmoothUp + 2) +
           1) *
          this->m_levGeoPtrsPoisson[i]->getBoxes().numCells();
    }
    // cross link the operators in the hierarchy
    for (int i = 0; i < OpsPtrs.size() - 1; ++i)
    {
      OpsPtrs[i]->setCoarserOpPtr(OpsPtrs[i + 1].get());
      OpsPtrs[i + 1]->setFinerOpPtr(OpsPtrs[i].get());
    }
    // post linking definitions. Used when the operator al level n
    // needs to have access to info provided by the operator at 
    // higher or lower levels. Thuss we do it here 
    // after the cross linking. 
    OpsPtrs[0]->upHierarchyPostDefinitions();
    //
    OpsPtrs.back()->downHierarchyPostDefinitions();
    // define the slover
    m_solver.define(OpsPtrs, refScheduleVector, options);

    // If required calculate the lowest eigenpairs and plot them
    // LevelData<FArrayBox> Eig;
    // this->findLowestMode(OpsPtrs[0], Eig);
    // char fn[80];
    // sprintf(fn, "eig_%01d.hdf5", c++);
    // IO::writeHDF5(fn, Eig , levGeo);
  }

  void cleanupPoissonSolver()
  {
    m_levGeoPtr->destroyMGLevGeos(this->m_levGeoPtrsPoisson);
  }
};

} // namespace Helpers
} // namespace Elliptic