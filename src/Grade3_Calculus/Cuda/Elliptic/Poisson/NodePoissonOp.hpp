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
#pragma once
#include "ElliKit.hpp"
#include "FArrayBox.H"
#include "LevelGeometry.H"
#include <cusp/eigen/spectral_radius.h>
#include <cusp/elementwise.h>
#include "cuda_profiler_api.h"
#include "FunctorsLibrary.hpp"
namespace Elliptic
{
namespace Poisson
{
template <typename Memory, typename BCType,
          typename Container = Elliptic::Helpers::NodeData<Real, Memory>>
class NodePoissonOp
{
protected:
  void fillInvDiags(){};

  // Member variables
  const bool m_areOperatorsSet;
  const BCType *m_bcPtr;
  const bool m_isDefined;
  const LevelGeometry *const m_levGeoPtr;
  IntVect m_GhostVector;
  IntVect m_RefRatio;
  LayoutData<Elliptic::PythonInterface::ElliKit<int, Real, cusp::host_memory,
                                                cusp::host_memory, BCType>>
      m_MatOp;
  typedef typename Elliptic::PythonInterface::ElliKit<int, Real, Memory, Memory,
                                                      BCType>::Matrix Matrix;
  typedef
      typename Elliptic::PythonInterface::ElliKit<int, Real, cusp::host_memory,
                                                  cusp::host_memory,
                                                  BCType>::Matrix HostMatrix;
  Matrix m_GlobalLaplacian, m_GlobalRestrictor, m_GlobalProlongator;

  HostMatrix m_EnvFieldToGlobalfield;
  typename Elliptic::PythonInterface::ElliKit<
      int, Real, Memory, Memory, BCType>::Relaxator *m_GlobalRelaxerPtr;
  typename Elliptic::PythonInterface::ElliKit<int, Real, Memory, Memory,
                                              BCType>::Preconditioner
      *m_GlobalPreconditionerPtr;
  cusp::array1d<bool, cusp::host_memory> m_mask;
  int m_Nrows;
#ifdef USE_CUSPARSE
  std::vector<cudaStream_t> m_Streams;
#else
  std::vector<int> m_Streams;
#endif
  Copier m_cp;
  NodePoissonOp *m_finerOpPtr = nullptr;
  NodePoissonOp *m_coarserOpPtr = nullptr;
  mutable decltype(m_EnvFieldToGlobalfield.row_offsets) m_GLIndPtr;

public:
  // Full Constructor.

  NodePoissonOp(const LevelGeometry *a_levGeoPtr, const BCType *a_bcPtr, const int level,
                bool setOperators = true, IntVect GhostVector = IntVect::Unit);

  // Destructor for good measure.
  ~NodePoissonOp();

  void extrapAndExchange(Container &a_state) const
  {
    // BCTools::extrapAllGhosts(a_state, 3);
    return;
  }

  // Computes lhs = J*L[phi].
  void applyOp(Container &a_lhs, const Container &a_phi, const Real a_time,
               const bool a_homogBCs) const;

  // Sets phi to a preconditioned value.
  // Typically, this will be phi = rhs / invDiags, then relaxed.
  void preCond(Container &a_phi, const Container &a_rhs, const Real a_time,
               const int a_relaxIters) const;

  // Apply relaxation to the homogeneous eq. L[phi] = rhs.
  // Do not set a_phi to zero first!
  // a_omega is the over (> 1.0) or under (< 1.0) relaxation parameter.
  void relax(Container &a_phi, const Container &a_rhs, const Real a_time,
             const int a_iters, const Real a_omega) const;

  void upHierarchyPostDefinitions();

  void downHierarchyPostDefinitions();

  IntVect getGhostVector() const { return m_GhostVector; }
   const LevelGeometry  &getLevGeo() const  { return *(this->m_levGeoPtr);}
  bool isDefined() const { return m_isDefined; }

  void createCoarsened(Container &a_crseCorrection,
                       const Container &a_fineCorrection,
                       Container &a_crseResidual,
                       const Container &a_fineResidual) const
  {
    int m = 1;
    for (int d = 0; d < SpaceDim; ++d)
    {
      m *= this->m_RefRatio[d];
    }
    int nrows = this->m_Nrows / m;
    a_crseCorrection().resize(nrows);
    a_crseResidual().resize(nrows);
  }

  void divergence(const LevelData<FluxBox> &a_vel, Container &DIV) const
  {
    Elliptic::Helpers::LevelDataHilbert<CudaFArrayBox<cusp::host_memory>> accumulator;
    accumulator.define(this->m_levGeoPtr->getBoxes(), 1);
    DataIterator dit = accumulator.dataIterator();
    for (dit.reset(); dit.ok(); ++dit)
    {
      CudaFArrayBox<cusp::host_memory> dummy(accumulator[dit].box(), 1);

      accumulator[dit].setVal(0.0);
      switch (SpaceDim)
      {
      case 3:
        m_MatOp[dit].ApplyOperator("Divergence z", a_vel[dit][2], dummy);
        accumulator[dit] += dummy;
      case 2:
        m_MatOp[dit].ApplyOperator("Divergence y", a_vel[dit][1], dummy);
        accumulator[dit] += dummy;
      default:
        m_MatOp[dit].ApplyOperator("Divergence x", a_vel[dit][0], dummy);
        accumulator[dit] += dummy;
      }
    }

    this->LDToContainer(accumulator, DIV);
  }
  void gradient(const Container &a_phi, LevelData<FluxBox> &a_gradPhi) const
  {
    const DisjointBoxLayout &grids = this->m_levGeoPtr->getBoxes();
    // remember, the LD Phi is defined on a grid which includes the overlapping
    // ghost points hence the following
    DisjointBoxLayout lhsGrids;
    lhsGrids.deepCopy(grids);
    lhsGrids.close();
    ResizeGrids<BCType> RG(grids, this->m_levGeoPtr->getDomain());
    lhsGrids.transform(RG); //
    Elliptic::Helpers::LevelDataHilbert<CudaFArrayBox<cusp::host_memory>> Phi;
    Phi.define(lhsGrids, 1, IntVect::Unit);

    DataIterator dit = Phi.dataIterator();
    this->ContainerToLD(a_phi, Phi);

    for (dit.reset(); dit.ok(); ++dit)
    {
      switch (SpaceDim)
      {
      case 3:
      {
        m_MatOp[dit].ApplyOperator("Gradient z", Phi[dit], a_gradPhi[dit][2]);
      }
      case 2:
      {
        m_MatOp[dit].ApplyOperator("Gradient y", Phi[dit], a_gradPhi[dit][1]);
      }
      default:
      {
        m_MatOp[dit].ApplyOperator("Gradient x", Phi[dit], a_gradPhi[dit][0]);
      }
      }
    }
  }

  void residual(Container &a_res, const Container &a_phi,
                const Container &a_rhs, const Real a_time,
                const bool a_homogBCs) const
  {
    this->applyOp(a_res, a_phi, a_time, a_homogBCs);

    a_res.axpy(a_rhs, a_res, 1.0, -1.0);
  }
  void restrictResidual(Container &a_crseRes, const Container &a_fineRes,
                        const Real a_time, const IntVect &a_refRatio) const
  {
    this->restrict(a_fineRes, a_crseRes);
  }
  void prolongIncrement(Container &a_finePhi, Container &a_crseCor,
                        const Real a_time, const IntVect &a_refRatio,
                        const int a_interpOrder) const
  {
    this->prolong(a_crseCor, a_finePhi);
  }

  double RelaxCoefficient()
  {
    // // this formula is taken from the Wikipedia article on SOR.
    // decltype(this->m_GlobalLaplacian) Diag, DinvA, R, Id;
    // cusp::array1d<Real, cusp::device_memory> diagonal;
    // // R=Laplacian^2;
    // cusp::copy(this->m_GlobalLaplacian, R);
    // thrust::replace_copy_if(R.values.begin(),
    //                         R.values.end(),
    //                         R.column_indices.begin(),
    //                         R.values.begin(),
    //                         Functors::equal_functor(this->m_Nrows),
    //                         0.0);
    // Id.resize(this->m_Nrows-1, this->m_Nrows, R.row_offsets[this->m_Nrows]);
    // thrust::copy(R.values.begin(), R.values.begin()+R.row_offsets[this->m_Nrows], Id.values.begin());
    // thrust::copy(R.column_indices.begin(), R.column_indices.begin()+R.row_offsets[this->m_Nrows], Id.column_indices.begin());
    // cusp::transpose(Id,R);
    // Id.resize(this->m_Nrows-1, this->m_Nrows-1, R.row_offsets[this->m_Nrows]);
    // thrust::copy(R.values.begin(), R.values.begin()+R.row_offsets[this->m_Nrows], Id.values.begin());
    // thrust::copy(R.column_indices.begin(), R.column_indices.begin()+R.row_offsets[this->m_Nrows], Id.column_indices.begin());
    // cusp::transpose(Id,R);

    // auto nu=cusp::eigen::estimate_rho_Dinv_A(R);
    // IO::tout(0) <<"Spectral radius of DinvL minus zero entry = " << nu << endl;
    // auto rho=cusp::eigen::estimate_rho_Dinv_A(this->m_GlobalLaplacian);
    // IO::tout(0) <<"Spectral radius of DinvL = " << rho << endl;

    // cusp::extract_diagonal(R, diagonal);
    // Diag.resize(this->m_Nrows, this->m_Nrows, this->m_Nrows);

    // thrust::transform(diagonal.begin(), diagonal.end(), Diag.values.begin(), Functors::invert_functor<Real>(1.));
    // thrust::sequence(Diag.column_indices.begin(), Diag.column_indices.end(), 0);
    // thrust::sequence(Diag.row_offsets.begin(), Diag.row_offsets.end(), 0);
    // cusp::multiply(Diag, R, DinvA);
    // thrust::fill(Diag.values.begin(), Diag.values.end(), 1.);
    // cusp::subtract(Diag, DinvA, R);

    // auto mu = cusp::eigen::estimate_spectral_radius(R, 20);
    // IO::tout(0) << "Spectral radius " << mu << endl;
    // IO::tout(0) << " optimal omega " << 1 + pow(mu / (1.0 + sqrt(1.0 - pow(mu, 2))), 2) << endl;
    return 0.9; //1+pow(mu / (1.0 + sqrt(1.0 - pow(mu, 2))), 2);
  }


  const NodePoissonOp *getCoarserOpPtr() const { return m_coarserOpPtr; }
  const NodePoissonOp *getFinerOpPtr() const { return m_finerOpPtr; }

  void setCoarserOpPtr(NodePoissonOp *Ptr) { m_coarserOpPtr = Ptr; }
  void setFinerOpPtr(NodePoissonOp *Ptr) { m_finerOpPtr = Ptr; }

  auto getOperatorMatrixPtr() { return &(this->m_GlobalLaplacian); }
  void LDToContainer(const LevelData<FArrayBox> &a_phi,
                     Container &WholeField) const;
  void ContainerToLD(const Container &WholeField,
                     LevelData<FArrayBox> &a_lhs) const;
  static int numComps() { return 1; }
  static IntVect ghostVect() { return IntVect::Unit; }
  static IntVect minBoxSize() { return 4 * IntVect::Unit; }
  static void eliminateZerosRowsAndColumns(
      typename Elliptic::PythonInterface::ElliKit<int, Real, Memory, Memory,
                                                  BCType>::Matrix &Op);
  int m_level = 0;

private:
  template <typename M1, typename M2>
  void makeCopyOfMatrix(const M1 &m1, M2 &m2)
  {
    m2.resize(m1.num_rows, m1.num_cols, m1.num_entries);
    thrust::copy(m1.values.begin(), m1.values.end(), m2.values.begin());
    thrust::copy(m1.row_offsets.begin(), m1.row_offsets.end(),
                 m2.row_offsets.begin());
    thrust::copy(m1.column_indices.begin(), m1.column_indices.end(),
                 m2.column_indices.begin());
  }

  void prolong(const Container &a_x, Container &a_y) const;
  void restrict(const Container &a_x, Container &a_y) const;
  NodePoissonOp() {} // disallow empty constructor and copy constructor
  NodePoissonOp(const NodePoissonOp &that) {}
  NodePoissonOp &operator=(const NodePoissonOp &that) {}
};
#include "NodePoissonOpImpl.hpp"
} // namespace Poisson
} // namespace Elliptic
