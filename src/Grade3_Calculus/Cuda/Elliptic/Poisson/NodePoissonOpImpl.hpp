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

// -----------------------------------------------------------------------------
// Constructor.
// -----------------------------------------------------------------------------

template <typename M, typename BCType, typename FABC>
NodePoissonOp<M, BCType, FABC>::NodePoissonOp(const LevelGeometry *a_levGeoPtr,
                                              const BCType *a_bcPtr,
                                              const int level,
                                              bool setOperators,
                                              IntVect GhostVector)
    : m_levGeoPtr(a_levGeoPtr), m_bcPtr(a_bcPtr), m_GhostVector(GhostVector),
      m_isDefined(true), m_Nrows(0), m_areOperatorsSet(setOperators),
      m_GlobalLaplacian(
          a_levGeoPtr->getBoxes().physDomain().domainBox().numPts(),
          a_levGeoPtr->getBoxes().physDomain().domainBox().numPts(), 0)

{
  CH_assert(a_levGeoPtr);
  CH_assert(a_bcPtr);
  // define copier
  const DisjointBoxLayout &grids = a_levGeoPtr->getBoxes();

  if (!setOperators)
    return;

  const DisjointBoxLayout &fineGrids = this->m_levGeoPtr->getBoxes();
  if (this->m_levGeoPtr->getCoarserPtr() != nullptr) {
    const DisjointBoxLayout &crseGrids =
        this->m_levGeoPtr->getCoarserPtr()->getBoxes();
    CH_assert(fineGrids.compatible(crseGrids));

    const ProblemDomain &fineDomain = fineGrids.physDomain();
    const ProblemDomain &crseDomain = crseGrids.physDomain();

    this->m_RefRatio = fineDomain.size() / crseDomain.size();
    CH_assert(this->m_RefRatio.product() > 1);
  } else {
    this->m_RefRatio = IntVect::Unit;
  }

  // Ellikit requires the physical grids. Internally, it will add ghost layers
  // to the internal boundaries.
  this->m_MatOp.define(grids); // define the LayoutData for the operators

  // Allocate/alias the workspaces.

  DataIterator dit = this->m_MatOp.dataIterator();

  HostMatrix dummy1, dummy2, DomainToValiGrid, GlobalLaplacian;
  GlobalLaplacian.resize(a_levGeoPtr->getBoxes().physDomain().domainBox().numPts(),
                         a_levGeoPtr->getBoxes().physDomain().domainBox().numPts(), 0);
  m_Nrows = 0;


  this->m_Streams.resize(grids.numBoxes(procID()));
  int c = 0;
  for (dit.reset(); dit.ok(); ++dit) {

#ifdef USE_CUSPARSE
    cudaStreamCreate(&m_Streams[c]);
#else
    m_Streams.at(c) = 0;
#endif
    // EllipticBCType *BCPtr = new BCTools::ProjNeumBC; // hardcoded for
    // projection

    this->m_MatOp[dit].define(a_levGeoPtr->getDXi(), a_bcPtr, grids[dit], grids,
                              1., grids.physDomain(), this->m_RefRatio, this->m_level,
                              this->m_Streams.at(c++));
    // delete BCPtr;

    auto &Laplacian = this->m_MatOp[dit].retrieveOperatorMatrix("Laplacian");
    auto &GridToDomain =
        this->m_MatOp[dit].retrieveOperatorMatrix("GridToDomain");
    auto &DomainToGrid =
        this->m_MatOp[dit].retrieveOperatorMatrix("DomainToGrid");

    // work on the GlobalLaplacian

    cusp::multiply(Laplacian, DomainToGrid, dummy1);
    cusp::multiply(GridToDomain, dummy1, dummy2);
    cusp::add(GlobalLaplacian, dummy2, GlobalLaplacian);
    m_MatOp[dit].deleteOperator("Laplacian");
    m_MatOp[dit].deleteOperator("LaplacianOG");
    m_MatOp[dit].deleteOperator("LaplacianNG");

    m_Nrows += GridToDomain.num_cols;
  }
  // we define the GlobalIndexPoiter of the row offsets.
  this->m_GLIndPtr = GlobalLaplacian.row_offsets;
  // and the mask to eliminate zero rows/columns
  this->m_mask.resize(this->m_GLIndPtr.size());
  this->m_mask[0] = true;
  thrust::transform(this->m_GLIndPtr.begin() + 1, this->m_GLIndPtr.end(),
                    this->m_GLIndPtr.begin(), this->m_mask.begin() + 1,
                    thrust::greater<Real>());
  // define the EnvelopingFieldToGlobalFieldMatrix
  cusp::array1d<int, cusp::host_memory> column_indices(
      a_levGeoPtr->getBoxes().physDomain().domainBox().numPts());

  thrust::sequence(column_indices.begin(), column_indices.end());

  m_EnvFieldToGlobalfield.resize(
      m_Nrows, a_levGeoPtr->getBoxes().physDomain().domainBox().numPts(),
      m_Nrows);
  thrust::sequence(m_EnvFieldToGlobalfield.row_offsets.begin(),
                   m_EnvFieldToGlobalfield.row_offsets.end());
  thrust::fill(m_EnvFieldToGlobalfield.values.begin(),
               m_EnvFieldToGlobalfield.values.end(), 1.0);
  thrust::copy_if(
      column_indices.begin(), column_indices.end(), this->m_mask.begin() + 1,
      m_EnvFieldToGlobalfield.column_indices.begin(), thrust::identity<bool>());

  this->makeCopyOfMatrix(GlobalLaplacian, this->m_GlobalLaplacian);
  this->eliminateZerosRowsAndColumns(this->m_GlobalLaplacian);

  m_GlobalRelaxerPtr = new typename Elliptic::PythonInterface::ElliKit<
      int, Real, M, M, BCType>::SORRelax(&(this->m_GlobalLaplacian), this->RelaxCoefficient());
  //  m_GlobalRelaxerPtr = new typename Elliptic::PythonInterface::ElliKit<
  //     int, Real, M, M, BCType>::PolynomialRelax(&(this->m_GlobalLaplacian));

  m_GlobalPreconditionerPtr =
      new typename Elliptic::PythonInterface::ElliKit<int, Real, M, M, BCType>::
          JacobiPreconditioner(&(this->m_GlobalLaplacian));
}

template <typename M, typename BCType, typename Container>
void NodePoissonOp<M, BCType, Container>::eliminateZerosRowsAndColumns(
    typename Elliptic::PythonInterface::ElliKit<int, Real, M, M, BCType>::Matrix
        &Op) {
  auto GLIndPtr = Op.row_offsets;
  cusp::array1d<bool, M> mask;
  mask.resize(GLIndPtr.size());
  mask[0] = true;
  thrust::transform(GLIndPtr.begin() + 1, GLIndPtr.end(), GLIndPtr.begin(),
                    mask.begin() + 1, thrust::greater<Real>());
  int Nrows = thrust::count(mask.begin(), mask.end(), true) - 1;
  Op.row_offsets.resize(Nrows + 1);
  Op.num_rows = Nrows;
  thrust::copy_if(GLIndPtr.begin(), GLIndPtr.end(), mask.begin(),
                  Op.row_offsets.begin(), thrust::identity<bool>());
  typename Elliptic::PythonInterface::ElliKit<int, Real, M, M, BCType>::Matrix
      dummy1;
  cusp::transpose(Op, dummy1);
  GLIndPtr = dummy1.row_offsets;
  mask.resize(GLIndPtr.size());
  mask[0] = true;
  thrust::transform(GLIndPtr.begin() + 1, GLIndPtr.end(), GLIndPtr.begin(),
                    mask.begin() + 1, thrust::greater<Real>());
  Nrows = thrust::count(mask.begin(), mask.end(), true) - 1;
  dummy1.row_offsets.resize(Nrows + 1);
  dummy1.num_rows = Nrows;
  thrust::copy_if(GLIndPtr.begin(), GLIndPtr.end(), mask.begin(),
                  dummy1.row_offsets.begin(), thrust::identity<bool>());

  cusp::transpose(dummy1, Op);
}

//------------------------------------------------------------------------------
// Restriction and Prolongation operations require access to the whole
// hierarchy. They are defined here. upHierarchy should be started from
// the finest level and will recurse itself all the way up.
// downHierarchy should be started from the coarsest and will move all the way
// down.
//------------------------------------------------------------------------------
template <typename M, typename BCType, typename Container>
void NodePoissonOp<M, BCType, Container>::upHierarchyPostDefinitions() {
  if (!this->m_areOperatorsSet)
    return;
  auto coarserPtr = this->m_coarserOpPtr;
  if (coarserPtr != nullptr) {

    const auto &grids = this->m_levGeoPtr->getBoxes();
    typedef decltype(this->m_EnvFieldToGlobalfield) Matrix;
    Matrix DomainToValidGrid, dummy, dummy1, GlobalRestrictor;

    int Nrows = coarserPtr->m_levGeoPtr->getDomainBox().numPts();
    int Ncolumns = this->m_levGeoPtr->getDomainBox().numPts();
    GlobalRestrictor.resize(Nrows, Ncolumns, 0);
    auto dit = grids.dataIterator();

    for (dit.reset(); dit.ok(); ++dit) {

      auto &Restrictor =
          this->m_MatOp[dit].retrieveOperatorMatrix("Restriction");
      auto &GridToDomain =
          this->m_MatOp[dit].retrieveOperatorMatrix("GridToDomain");
      cusp::transpose(GridToDomain, DomainToValidGrid);
      cusp::multiply(Restrictor, DomainToValidGrid, dummy);
      this->m_MatOp[dit].deleteOperator("Restriction");
      this->m_MatOp[dit].deleteOperator("Prolongation");
      auto &CoarserGridToDomain =
          coarserPtr->m_MatOp[dit].retrieveOperatorMatrix("GridToDomain");
      cusp::multiply(CoarserGridToDomain, dummy, dummy1);
      cusp::add(GlobalRestrictor, dummy1, GlobalRestrictor);
    }
    //
    this->makeCopyOfMatrix(GlobalRestrictor, this->m_GlobalRestrictor);
    this->eliminateZerosRowsAndColumns(this->m_GlobalRestrictor);

    coarserPtr->upHierarchyPostDefinitions();
  }
}
template <typename M, typename BCType, typename Container>
void NodePoissonOp<M, BCType, Container>::downHierarchyPostDefinitions() {
  if (!this->m_areOperatorsSet)
    return;
  auto finerPtr = this->m_finerOpPtr;

  typedef decltype(this->m_GlobalLaplacian) Matrix;

  if (finerPtr != nullptr) {
    cusp::transpose(
        finerPtr->m_GlobalRestrictor,
        this->m_GlobalProlongator); // note that this is not normalized.
    finerPtr->downHierarchyPostDefinitions();
  }
}
// -----------------------------------------------------------------------------
// Destructor for good reasons!
// -----------------------------------------------------------------------------
template <typename M, typename BCType, typename FABC>
NodePoissonOp<M, BCType, FABC>::~NodePoissonOp() {

  m_Streams.resize(0);
  m_bcPtr = nullptr;
  delete m_GlobalRelaxerPtr;
  m_GlobalRelaxerPtr = nullptr;
  delete m_GlobalPreconditionerPtr;
  m_GlobalPreconditionerPtr = nullptr;
}

// -----------------------------------------------------------------------------
// prolongIncrement with a_interpOrder > 0 requires ghosts be filled.
// You do not need to worry about the edge or vertex ghosts.
// -----------------------------------------------------------------------------
// template <typename M, typename BCType, typename Container>
// void NodePoissonOp<M, BCType, Container>::fillGhostsHomog(Container &a_state,
//                                                           const Real a_time)
//                                                           const
// {
//     //m_bcPtr->fillGhosts(a_state, a_time, *(this->m_levGeoPtr), true);
//     return;
// }
// -----------------------------------------------------------------------------
// packs an entire LayoutData into a single entity and vicecersa
// -----------------------------------------------------------------------------
template <typename M, typename BCType, typename Container>
void NodePoissonOp<M, BCType, Container>::LDToContainer(
    const LevelData<FArrayBox> &a_phi,
    Container &WholeField) const { //

  auto Nrows = this->m_levGeoPtr->getDomainBox().numPts();
  HostMatrix EnvelopingField;
  EnvelopingField.resize(Nrows, 1, 0);

  HostMatrix dummy;
  DataIterator dit = a_phi.dataIterator();
  auto &grids = this->m_levGeoPtr->getBoxes();
  for (dit.reset(); dit.ok(); ++dit) {
    CudaFArrayBox<M> CopyOfPhiOnValid(grids[dit], 1);
    CopyOfPhiOnValid.copyHetero(a_phi[dit]);

    HostMatrix LocalField(CopyOfPhiOnValid.m_view.size(), 1,
                      CopyOfPhiOnValid.m_view.size());
    LocalField.values = CopyOfPhiOnValid.m_view;
    LocalField.column_indices =
        cusp::array1d<int, M>(CopyOfPhiOnValid.m_view.size(), 0);
    LocalField.row_offsets =
        cusp::array1d<int, M>(CopyOfPhiOnValid.m_view.size() + 1);
    thrust::sequence(LocalField.row_offsets.begin(),
                     LocalField.row_offsets.end(), 0, 1);
    auto &GridToDomain = m_MatOp[dit].retrieveOperatorMatrix("GridToDomain");

    cusp::multiply(GridToDomain, LocalField, dummy);

    cusp::add(EnvelopingField, dummy, EnvelopingField);
  }
  // now restrict EnvelopingField field to nonzero entries
  cusp::multiply(this->m_EnvFieldToGlobalfield, EnvelopingField, dummy);
  WholeField().resize(this->m_Nrows);
  cusp::convert(dummy, WholeField());
}
//
template <typename M, typename BCType, typename Container>
void NodePoissonOp<M, BCType, Container>::ContainerToLD(
    const Container &WholeField, LevelData<FArrayBox> &a_lhs) const {
  //
  // repack the data

  HostMatrix EnvelopingField;
  EnvelopingField.resize(this->m_Nrows, 1, this->m_Nrows);
  thrust::copy(WholeField().begin(), WholeField().end(), EnvelopingField.values.begin());
  thrust::fill(EnvelopingField.column_indices.begin(),
               EnvelopingField.column_indices.end(), 0);
  thrust::sequence(EnvelopingField.row_offsets.begin(),
                   EnvelopingField.row_offsets.begin() + this->m_Nrows + 1);

  DataIterator dit = a_lhs.dataIterator();

  for (dit.reset(); dit.ok(); ++dit) {
    auto &DomainToGrid = m_MatOp[dit].retrieveOperator("DomainToGrid");
    CudaFArrayBox<cusp::host_memory> LF(DomainToGrid.getMyOutputGrid(), 1);
    Container LocalField;
    HostMatrix dummy;
    HostMatrix dummy2;
    HostMatrix LocalFieldMat;
    cusp::transpose(this->m_EnvFieldToGlobalfield, dummy2);
    cusp::multiply(*(DomainToGrid.getMatrixPtr()), dummy2, dummy);
    cusp::multiply(dummy, EnvelopingField, LocalFieldMat);
    cusp::convert(LocalFieldMat, LocalField());
    thrust::copy(LocalField().begin(), LocalField().end(), LF.m_view.begin());

    a_lhs[dit].copyHetero<CudaFArrayBox<cusp::host_memory>>(LF);
  }
}
// -----------------------------------------------------------------------------
// Computes lhs = J*L[phi].
// -----------------------------------------------------------------------------
template <typename M, typename BCType, typename Container>
void NodePoissonOp<M, BCType, Container>::applyOp(Container &a_lhs,
                                                  const Container &a_phi,
                                                  const Real a_time,
                                                  const bool a_homogBCs) const {
  cusp::multiply(this->m_GlobalLaplacian, a_phi(), a_lhs());
}
//
template <typename M, typename BCType, typename Container>
void NodePoissonOp<M, BCType, Container>::restrict(const Container &a_x,
                                                   Container &a_y) const {
  cusp::multiply(this->m_GlobalRestrictor, a_x(), a_y());
}
//
template <typename M, typename BCType, typename Container>
void NodePoissonOp<M, BCType, Container>::prolong(const Container &a_x,
                                                  Container &a_y) const {
  Container dummy;
  dummy.create(a_y);
  cusp::multiply(this->m_coarserOpPtr->m_GlobalProlongator, a_x(), dummy());
  int m = 1;

  for (int d = 0; d < SpaceDim; ++d) {
    m *= this->m_RefRatio[d];
  }
  a_y.incr(dummy, (Real)m);
}

// -----------------------------------------------------------------------------
// Sets phi to a preconditioned value.
// Typically, this will be phi = rhs / invDiags, then relaxed.
// -----------------------------------------------------------------------------
template <typename M, typename BCType, typename Container>
void NodePoissonOp<M, BCType, Container>::preCond(
    Container &a_phi, const Container &a_rhs, const Real a_time,
    const int a_relaxIters) const {
  (*m_GlobalPreconditionerPtr)(a_rhs(), a_phi());

  const Real omega = 1.0;
  this->relax(a_phi, a_rhs, a_time, a_relaxIters, omega);
}

// -----------------------------------------------------------------------------
// Apply relaxation to the homogeneous eq. L[phi] = rhs.
// Do not set a_phi to zero first!
// a_omega is the over (> 1.0) or under (< 1.0) relaxation parameter.
// -----------------------------------------------------------------------------
template <typename M, typename BCType, typename Container>
void NodePoissonOp<M, BCType, Container>::relax(Container &a_phi,
                                                const Container &a_rhs,
                                                const Real a_time,
                                                const int a_iters,
                                                const Real a_omega) const {
  for (int i = 0; i < a_iters; ++i) {
    (*(m_GlobalRelaxerPtr))(a_rhs(), a_phi());
  }
}
