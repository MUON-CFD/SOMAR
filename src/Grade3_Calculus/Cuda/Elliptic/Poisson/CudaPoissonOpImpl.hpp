
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
// -----------------------------------------------------------------------------
// Constructor.
// The one of the BC ptrs can be nullptr, but not both.
// If you supply them, ws1 and ws3 need 1 comp and ws2 needs SpaceDim comps.
// -----------------------------------------------------------------------------

template <typename FABC, typename M, typename BCType>
PoissonOp<FABC, M, BCType>::PoissonOp(const LevelGeometry *a_levGeoPtr,
                                      const BCType *a_bcPtr,
                                      const int level,
                                      bool setOperators,
                                      IntVect GhostVector)
    : m_levGeoPtr(a_levGeoPtr), m_bcPtr(a_bcPtr), m_GhostVector(GhostVector),
      m_isDefined(true),m_level(level),
      m_GlobalLaplacian(
          a_levGeoPtr->getBoxes().physDomain().domainBox().numPts(),
          a_levGeoPtr->getBoxes().physDomain().domainBox().numPts(), 0)
{
  CH_assert(a_levGeoPtr);
  CH_assert(a_bcPtr);
  // define copier
  const DisjointBoxLayout &grids = a_levGeoPtr->getBoxes();
  m_cp.exchangeDefine(grids, IntVect::Unit);

  if (!setOperators)
    return;

  RealVect RR =
      (a_levGeoPtr->getCoarserPtr() != nullptr)
          ? a_levGeoPtr->getCoarserPtr()->getDXi() / (a_levGeoPtr->getDXi())
          : RealVect::Unit;

  IntVect RefRatio(D_DECL6((int)RR[0], (int)RR[1], (int)RR[2], 0, 0, 0));
  // Ellikit requires the physical grids. Internally, it will add ghost layers
  // to the internal boundaries.
  this->m_MatOp.define(grids); // define the LayoutData for the operators

  // Allocate/alias the workspaces.

  this->m_ws3.define(grids, 1);

  DataIterator dit = this->m_MatOp.dataIterator();

  m_Streams.resize(grids.numBoxes(procID()));
  int c = 0;
  for (dit.reset(); dit.ok(); ++dit)
  {
#ifdef USE_CUSPARSE
    cudaStreamCreate(&m_Streams[c]);
#else
    m_Streams.at(c) = 0;
#endif
    Real constexpr omega = 1.0; // use this for SOR relaxator
    // EllipticBCType *BCPtr = new BCTools::ProjNeumBC; // hardcoded for
    // projection
    std::vector<std::string> Operators{"Laplacian", "LaplacianDMinusOne", "LaplacianNG",
                                       "LaplacianOG", "LaplacianNGDMinusOne", "LaplacianOGDMinusOne",
                                       "Restriction"};
    m_MatOp[dit].operatorsNeeded(Operators);
    m_MatOp[dit].define(a_levGeoPtr->getDXi(), a_bcPtr, grids[dit], grids,
                        omega, grids.physDomain(), RefRatio, m_level, m_Streams.at(c++));
    // delete BCPtr;
  }
}

// -----------------------------------------------------------------------------
// Destructor for good measure.
// -----------------------------------------------------------------------------
template <typename FABC, typename M, typename BCType>
PoissonOp<FABC, M, BCType>::~PoissonOp()
{
  m_ws3.clear();
  m_Streams.resize(0);
  m_bcPtr = nullptr;
}

// -----------------------------------------------------------------------------
// prolongIncrement with a_interpOrder > 0 requires ghosts be filled.
// You do not need to worry about the edge or vertex ghosts.
// -----------------------------------------------------------------------------
template <typename FABC, typename M, typename BCType>
void PoissonOp<FABC, M, BCType>::fillGhostsHomog(
    Elliptic::Helpers::LevelDataHilbert<FABC> &a_state, const Real a_time) const
{
  // m_bcPtr->fillGhosts(a_state, a_time, *(this->m_levGeoPtr), true);
  return;
}

// -----------------------------------------------------------------------------
// Computes lhs = J*L[phi].
// -----------------------------------------------------------------------------
template <typename FABC, typename M, typename BCType>
void PoissonOp<FABC, M, BCType>::applyOp(Elliptic::Helpers::LevelDataHilbert<FABC> &a_lhs,
                                         Elliptic::Helpers::LevelDataHilbert<FABC> &a_phi,
                                         const Real a_time,
                                         const bool a_homogBCs) const
{
  CH_assert(a_lhs.getBoxes().physDomain() == this->m_levGeoPtr->getDomain());
  CH_assert(a_phi.getBoxes().physDomain() == this->m_levGeoPtr->getDomain());
  CH_assert(a_lhs.getBoxes() == this->m_levGeoPtr->getBoxes());

  DataIterator dit = a_phi.dataIterator();
  const DisjointBoxLayout &grids = a_lhs.getBoxes();

  // Util::FuncTimer("extrap and exchange ", [&]() -> void {
  // cudaProfilerStart() ;
  this->extrapAndExchange(a_phi);
  // cudaProfilerStop();
  // });

  for (dit.reset(); dit.ok(); ++dit)
  {
    m_MatOp[dit].ApplyOperator("Laplacian", a_phi[dit], a_lhs[dit]);
  }

  // If space needs to be saved the formulation below is equivalent,
  // and does not require the storage of the Laplacian

  //         //Util::FuncTimer("LapOnG ", [&]() -> void {
  //         m_MatOp[dit].ApplyOperator("LaplacianOG", a_phi[dit], a_lhs[dit]);
  //         //});
  //     CudaFArrayBox<M> phiInt(grids[dit], 1);
  //     phiInt.copy(a_phi[dit]);

  //     //Util::FuncTimer("LapNG", [&]() -> void {
  //         m_MatOp[dit].ApplyOperator("LaplacianNG", phiInt, m_ws3[dit]);
  //         //});
  //     thrust::transform(m_ws3[dit].begin(),
  //                       m_ws3[dit].end(),
  //                       a_lhs[dit].begin(),
  //                       a_lhs[dit].begin(),
  //                       thrust::plus<Real>());
  // }
}

template <typename FABC, typename M, typename BCType>
void PoissonOp<FABC, M, BCType>::restrict(
    const Elliptic::Helpers::LevelDataHilbert<FABC> &a_x, Elliptic::Helpers::LevelDataHilbert<FABC> &a_y) const
{

  DataIterator dit = a_x.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
  {
    m_MatOp[dit].ApplyOperator("Restriction", a_x[dit], a_y[dit]);
  }
}
template <typename FABC, typename M, typename BCType>
void PoissonOp<FABC, M, BCType>::prolong(
    const Elliptic::Helpers::LevelDataHilbert<FABC> &a_x, Elliptic::Helpers::LevelDataHilbert<FABC> &a_y) const
{
  assert(0);
}

// -----------------------------------------------------------------------------
// Sets phi to a preconditioned value.
// Typically, this will be phi = rhs / invDiags, then relaxed.
// -----------------------------------------------------------------------------
template <typename FABC, typename M, typename BCType>
void PoissonOp<FABC, M, BCType>::preCond(
    Elliptic::Helpers::LevelDataHilbert<FABC> &a_phi, const Elliptic::Helpers::LevelDataHilbert<FABC> &a_rhs,
    const Real a_time, const int a_relaxIters) const
{
  const DisjointBoxLayout &grids = this->m_levGeoPtr->getBoxes();
  DataIterator dit = grids.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
  {
    m_MatOp[dit].Precond("LaplacianNG", a_rhs[dit], a_phi[dit]);
  }

  const Real omega = 1.0;
  this->relax(a_phi, a_rhs, a_time, a_relaxIters, omega);
}

// -----------------------------------------------------------------------------
// Apply relaxation to the homogeneous eq. L[phi] = rhs.
// Do not set a_phi to zero first!
// a_omega is the over (> 1.0) or under (< 1.0) relaxation parameter.
// -----------------------------------------------------------------------------
template <typename FABC, typename M, typename BCType>
void PoissonOp<FABC, M, BCType>::relax(Elliptic::Helpers::LevelDataHilbert<FABC> &a_phi,
                                       const Elliptic::Helpers::LevelDataHilbert<FABC> &a_rhs,
                                       const Real a_time, const int a_iters,
                                       const Real a_omega) const
{
  DataIterator dit = a_phi.dataIterator();
  const DisjointBoxLayout &grids = a_rhs.getBoxes();
  for (int i = 0; i < a_iters; ++i)
  {
    this->extrapAndExchange(a_phi);
    for (dit.reset(); dit.ok(); ++dit)
    {

      int Nphi;
      auto phiIt = a_phi[dit].subVolumeIt(grids[dit], 0, 1, Nphi);
      FABC dummy(grids[dit], 1);

      m_MatOp[dit].ApplyOperator("LaplacianOG", a_phi[dit], m_ws3[dit]);
      // dummy=Lphi
      thrust::transform(a_rhs[dit].begin(), a_rhs[dit].end(),
                        m_ws3[dit].begin(), dummy.begin(),
                        thrust::minus<Real>());
      // now dummy=rhs-LapOG(phi)
      FABC phiInt(grids[dit], 1);
      phiInt.copy(a_phi[dit]);

      m_MatOp[dit].Relax("LaplacianNG", dummy.m_varr, phiInt.m_varr);
      thrust::copy(phiInt.begin(), phiInt.end(), phiIt);
    }
  }
  this->extrapAndExchange(a_phi);
}

// #define USE_RELAXATION_METRICS

#ifdef USE_RELAXATION_METRICS
#define relaxationMetrics(p, r, t, i)                           \
  do                                                            \
  {                                                             \
    Elliptic::Helpers::LevelDataHilbert<FArrayBox<M>> localRes; \
    this->create(localRes, p);                                  \
    this->residual(localRes, p, r, t, true);                    \
    Real n = this->norm(localRes, 2);                           \
    pout() << "Iter " << (i) << " res: " << n << "\n";          \
  } while (0)
#else
#define relaxationMetrics(c, r, t, i)
#endif

template <typename FABC, typename M, typename BCType>
void PoissonOp<FABC, M, BCType>::divergence(
    const LevelData<CudaFluxBox<M>> &a_vel,
    Elliptic::Helpers::LevelDataHilbert<FABC> &a_div) const
{
  DataIterator dit = a_div.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    a_div[dit].setVal(0.);
    FABC dummy(a_div[dit].box(), 1);
    switch (SpaceDim)
    {
    case 3:
      m_MatOp[dit].ApplyOperator("Divergence z", a_vel[dit][2], dummy);
      a_div[dit] += dummy;
    case 2:
      m_MatOp[dit].ApplyOperator("Divergence y", a_vel[dit][1], dummy);
      a_div[dit] += dummy;
    default:
      m_MatOp[dit].ApplyOperator("Divergence x", a_vel[dit][0], dummy);
      a_div[dit] += dummy;
    }
  }
}
template <typename FABC, typename M, typename BCType>
void PoissonOp<FABC, M, BCType>::gradient(
    Elliptic::Helpers::LevelDataHilbert<FABC> &a_phi,
    LevelData<CudaFluxBox<M>> &a_gradPhi) const
{
  const DisjointBoxLayout &grids = this->m_levGeoPtr->getBoxes();
  this->extrapAndExchange(a_phi);
  DataIterator dit = a_phi.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    switch (SpaceDim)
    {
    case 3:
    {
      m_MatOp[dit].ApplyOperator("Gradient z", a_phi[dit], a_gradPhi[dit][2]);
    }
    case 2:
    {
      m_MatOp[dit].ApplyOperator("Gradient y", a_phi[dit], a_gradPhi[dit][1]);
    }
    default:
    {
      m_MatOp[dit].ApplyOperator("Gradient x", a_phi[dit], a_gradPhi[dit][0]);
    }
    }
  }
}

// Restrict to coarser MG depth:
// a_crseRes = I[h->2h](a_fineRes).
// Do not compute the residual, that was done for us. Just coarsen.
// a_crse can be overwritten. This op is at the fine level.
template <typename FABC, typename M, typename BCType>
void PoissonOp<FABC, M, BCType>::restrictResidual(
    Elliptic::Helpers::LevelDataHilbert<FABC> &a_crseRes, const Elliptic::Helpers::LevelDataHilbert<FABC> &a_fineRes,
    const Real a_time, const IntVect &a_refRatio
    // const  typename Elliptic::Poisson::PoissonOp<FABC, M, BCType>
    // *a_OpPtr
    ) const
{
  CH_TIME("FullWeightingPS::restrictResidual");

  // Tripping an assert in C is easier to debug than an invalid memory access
  // in fortran!
  CH_assert(a_crseRes.nComp() == a_fineRes.nComp());
  this->restrict(a_fineRes, a_crseRes);
}
template <typename FABC, typename M, typename BCType>
void PoissonOp<FABC, M, BCType>::prolongIncrement(
    Elliptic::Helpers::LevelDataHilbert<FABC> &a_finePhi, Elliptic::Helpers::LevelDataHilbert<FABC> &a_crseCor,
    const Real a_time, const IntVect &a_refRatio,
    // const PoissonOp<FABC, M, BCType>  *a_crseOpPtr,
    const int a_interpOrder) const
{
  BUG(" This one works only for 1 component FABs.");
  CH_TIME("ZeroAvgConstInterpPS::prolongIncrement");

  // Gather grids, domains, refinement ratios...
  // const DisjointBoxLayout &fineGrids = a_finePhi.getBoxes();
  // const DisjointBoxLayout &crseGrids = a_crseCor.getBoxes();
  const DisjointBoxLayout &fineGrids = this->m_levGeoPtr->getBoxes();
  const DisjointBoxLayout &crseGrids =
      this->m_levGeoPtr->getCoarserPtr()->getBoxes();
  CH_assert(fineGrids.compatible(crseGrids));

  const ProblemDomain &fineDomain = fineGrids.physDomain();
  const ProblemDomain &crseDomain = crseGrids.physDomain();

  const IntVect mgRefRatio = fineDomain.size() / crseDomain.size();
  CH_assert(mgRefRatio.product() > 1);

  // These will accumulate averaging data.
  Real localSum = 0.0;
  Real localVol = 0.0;
  const Real dxProduct = this->m_levGeoPtr->getDXi().product();
  // a_crseOpPtr->prolong(a_crseCor, a_finePhi);
  DataIterator dit = fineGrids.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {

    // Create references for convenience
    auto &fineFAB = a_finePhi[dit];
    const auto &crseFAB = a_crseCor[dit];
    const Box &fineValid = fineGrids[dit];
    const Box &crseValid = crseGrids[dit];
    CH_assert(fineFAB.nComp() ==
              1);             // this things works only for pressure (one component)
    FABC VolEl(fineValid, 1); // used to accumulate results

    // const HostFArrayBox &JinvFAB = (*m_CCJinvPtr)[dit];

    int Nfine;
    auto fineIt = fineFAB.subVolumeIt(fineValid, 0, 1, Nfine);
    int Ncrse;
    auto coarseIt = crseFAB.subVolumeIt(crseValid, 0, 1, Ncrse);
    CH_assert(Nfine == mgRefRatio.product() * Ncrse);

    fineToCoarse_functor P(fineValid.size(), crseValid.size(), mgRefRatio);

    auto permutator = thrust::make_transform_iterator(
        thrust::counting_iterator<unsigned long>(0), P);
    auto coarsePermIt = thrust::make_permutation_iterator(coarseIt, permutator);
    thrust::transform(fineIt, fineIt + Nfine, coarsePermIt, fineIt,
                      thrust::plus<Real>());

    // also compute the local sum and volume

    VolEl.setVal(dxProduct);
    // VolEl.divide(JinvFAB, fineValid, 0, 0);
    localVol = thrust::reduce(VolEl.begin(), VolEl.end());

    VolEl.mult(fineFAB, fineValid, 0, 0);
    localSum = thrust::reduce(VolEl.begin(), VolEl.end());
  }

  // Compute global sum (this is where the MPI communication happens)
#ifdef CH_MPI
  Real globalSum = 0.0;
  int result = MPI_Allreduce(&localSum, &globalSum, 1, MPI_CH_REAL, MPI_SUM,
                             Chombo_MPI::comm);

  if (result != MPI_SUCCESS)
  {
    MayDay::Error("Sorry, but I had a communication error in "
                  "ZeroAvgConstInterpPS::prolongIncrement");
  }

  Real globalVol = 0.0;
  result = MPI_Allreduce(&localVol, &globalVol, 1, MPI_CH_REAL, MPI_SUM,
                         Chombo_MPI::comm);

  if (result != MPI_SUCCESS)
  {
    MayDay::Error("Sorry, but I had a communication error in "
                  "ZeroAvgConstInterpPS::prolongIncrement");
  }

#else
  Real globalSum = localSum;
  Real globalVol = localVol;
#endif

  // Remove the average from phi.
  Real avgPhi = globalSum / globalVol;

  for (dit.reset(); dit.ok(); ++dit)
  {
    a_finePhi[dit] -= avgPhi;
  }
}
