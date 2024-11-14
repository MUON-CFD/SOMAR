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

#include "ElliKit.hpp"
//#include "CudaMGEllipticOp.H"
#include "BCTools.H"
#include "EllipticOp.H"
#include "FArrayBox.H"
#include "FiniteDiff.H"
#include "FunctionTimer.hpp"
#include "LevelDataHilbert.hpp"
#include "LevelGeometry.H"
#include "ResizeGrids.hpp"
#include "cuda_profiler_api.h"
namespace Elliptic
{
namespace Poisson
{
template <typename FABContainer, typename Memory, typename BCType>
class PoissonOp
{
public:
  // Full Constructor.

  PoissonOp(const LevelGeometry *a_levGeoPtr, const BCType *a_bcPtr, const int level,
                bool setOperators = true, IntVect GhostVector = IntVect::Unit);

  // Destructor for good measure.
  ~PoissonOp();



  void extrapAndExchange(Elliptic::Helpers::LevelDataHilbert<FABContainer> &a_state) const
  {
    // BCTools::extrapAllGhosts(a_state, 3);

    a_state.exchange(this->m_cp);
  }
  // prolongIncrement with a_interpOrder > 0 requires ghosts be filled.
  // You do not need to worry about the edge or vertex ghosts.
  void fillGhostsHomog(Elliptic::Helpers::LevelDataHilbert<FABContainer> &a_state,
                       const Real a_time) const;
  // use Matrix vector mult for the restriction
  void restrict(const Elliptic::Helpers::LevelDataHilbert<FABContainer> &a_x,
                Elliptic::Helpers::LevelDataHilbert<FABContainer> &a_y) const;
  void prolong(const Elliptic::Helpers::LevelDataHilbert<FABContainer> &a_x,
               Elliptic::Helpers::LevelDataHilbert<FABContainer> &a_y) const;
  // Computes lhs = J*L[phi].
  void applyOp(Elliptic::Helpers::LevelDataHilbert<FABContainer> &a_lhs,
               Elliptic::Helpers::LevelDataHilbert<FABContainer> &a_phi, const Real a_time,
               const bool a_homogBCs) const;

  // Sets phi to a preconditioned value.
  // Typically, this will be phi = rhs / invDiags, then relaxed.
  void preCond(Elliptic::Helpers::LevelDataHilbert<FABContainer> &a_phi,
               const Elliptic::Helpers::LevelDataHilbert<FABContainer> &a_rhs, const Real a_time,
               const int a_relaxIters) const;

  // Apply relaxation to the homogeneous eq. L[phi] = rhs.
  // Do not set a_phi to zero first!
  // a_omega is the over (> 1.0) or under (< 1.0) relaxation parameter.
  void relax(Elliptic::Helpers::LevelDataHilbert<FABContainer> &a_phi,
             const Elliptic::Helpers::LevelDataHilbert<FABContainer> &a_rhs, const Real a_time,
             const int a_iters, const Real a_omega) const;

  void upHierarchyPostDefinitions() {}

  void downHierarchyPostDefinitions() {}
  IntVect getGhostVector() const { return m_GhostVector; }
  bool isDefined() const { return m_isDefined; }
  void
  createCoarsened(Elliptic::Helpers::LevelDataHilbert<FABContainer> &a_crseCorrection,
                  const Elliptic::Helpers::LevelDataHilbert<FABContainer> &a_fineCorrection,
                  Elliptic::Helpers::LevelDataHilbert<FABContainer> &a_crseResidual,
                  const Elliptic::Helpers::LevelDataHilbert<FABContainer> &a_fineResidual) const
  {
    const LevelGeometry *coarseLG = m_levGeoPtr->getCoarserPtr();
    const DisjointBoxLayout &coarseGrids = coarseLG->getBoxes();
    DisjointBoxLayout coarseGridsCorr;
    coarseGridsCorr.deepCopy(coarseGrids);
    coarseGridsCorr.close();
    ResizeGrids<BCType> RG(coarseGrids, coarseLG->getDomain());
    coarseGridsCorr.transform(RG);
    this->createCoarsened(a_crseCorrection, a_fineCorrection, coarseGridsCorr,
                          IntVect::Unit);
    this->createCoarsened(a_crseResidual, a_fineResidual, coarseGrids,
                          IntVect::Unit);
  }
  void createCoarsened(Elliptic::Helpers::LevelDataHilbert<FABContainer> &a_crse,
                       const Elliptic::Helpers::LevelDataHilbert<FABContainer> &a_fine,
                       const DisjointBoxLayout &a_crseGrids,
                       const IntVect a_refRatio) const
  {
    a_crse.define(a_crseGrids, a_fine.nComp(), a_fine.ghostVect());
  }
  void divergence(const LevelData<CudaFluxBox<Memory>> &a_vel,
                  Elliptic::Helpers::LevelDataHilbert<FABContainer> &a_div) const;
  void gradient(Elliptic::Helpers::LevelDataHilbert<FABContainer> &a_phi,
                LevelData<CudaFluxBox<Memory>> &a_gradPhi) const;
  void residual(Elliptic::Helpers::LevelDataHilbert<FABContainer> &a_lhs,
                Elliptic::Helpers::LevelDataHilbert<FABContainer> &a_phi,
                const Elliptic::Helpers::LevelDataHilbert<FABContainer> &a_rhs, const Real a_time,
                const bool a_homogBCs) const
  {
    this->applyOp(a_lhs, a_phi, a_time, a_homogBCs);
    a_lhs.axpy(a_rhs, a_lhs, 1.0, -1.0);
  }
  void
  restrictResidual(Elliptic::Helpers::LevelDataHilbert<FABContainer> &a_crseRes,
                   const Elliptic::Helpers::LevelDataHilbert<FABContainer> &a_fineRes,
                   const Real a_time, const IntVect &a_refRatio
                   // const PoissonOp<FABContainer, Memory, BCType> *a_OpPtr
                   ) const;
  void prolongIncrement(
      Elliptic::Helpers::LevelDataHilbert<FABContainer> &a_finePhi,
      Elliptic::Helpers::LevelDataHilbert<FABContainer> &a_crseCor, const Real a_time,
      const IntVect &a_refRatio,
      // const PoissonOp<FABContainer, Memory, BCType> *a_crseOpPtr,
      const int a_interpOrder) const;

  const PoissonOp *getCoarserOpPtr() const { return m_coarserOpPtr; }
  const PoissonOp *getFinerOpPtr() const { return m_finerOpPtr; }

  void setCoarserOpPtr(PoissonOp *Ptr) { m_coarserOpPtr = Ptr; }
  void setFinerOpPtr(PoissonOp *Ptr) { m_finerOpPtr = Ptr; }
  // -----------------------------------------------------------------------------
  // packs an entire LayoutData into a single entity and vicecersa
  // -----------------------------------------------------------------------------

  void LDToContainer(const LevelData<FArrayBox> &LD,
                     Elliptic::Helpers::LevelDataHilbert<FABContainer> &Container) const
  { //
    DataIterator dit = LD.dataIterator();
    for (dit.reset(); dit.ok(); ++dit)
    {
      Container[dit].copyHetero(LD[dit]);
    }
  }
  //

  void ContainerToLD(const Elliptic::Helpers::LevelDataHilbert<FABContainer> &Container,
                     LevelData<FArrayBox> &LD) const
  {
    DataIterator dit = LD.dataIterator();
    for (dit.reset(); dit.ok(); ++dit)
    {
      LD[dit].copyHetero(Container[dit]);
    }
  }
  static int numComps() { return 1; }
  static IntVect ghostVect() { return IntVect::Unit; }
  static IntVect minBoxSize() { return 4 * IntVect::Unit; }

protected:
  // Computes 1.0 / diagonal matrix elements.
  // Comp 0 will have homog BCs rolled in.
  // Comp 1 will not.
  void fillInvDiags(){};

  // Member variables
  const BCType *m_bcPtr;
  const bool m_isDefined;
  mutable Elliptic::Helpers::LevelDataHilbert<CudaFArrayBox<Memory>> m_ws3;
  const LevelGeometry *const m_levGeoPtr;
  IntVect m_GhostVector;

  LayoutData<Elliptic::PythonInterface::ElliKit<
      int, Real, typename FABContainer::MyMemory, Memory, BCType>>
      m_MatOp;
  typename Elliptic::PythonInterface::ElliKit<
      int, Real, typename FABContainer::MyMemory, Memory, BCType>::Operator
      m_GlobalLaplacian;
  cusp::array1d<bool, Memory> m_mask;
  int m_Nrows;
#ifdef USE_CUSPARSE
  std::vector<cudaStream_t> m_Streams;
#else
  std::vector<int> m_Streams;
#endif
  Copier m_cp;
  PoissonOp *m_finerOpPtr = nullptr;
  PoissonOp *m_coarserOpPtr = nullptr;
  int m_level = 0;

private:
  PoissonOp() {} // disallow empty constructor
  PoissonOp(const PoissonOp &source) {}
  PoissonOp &operator=(const PoissonOp &that) {}
};
#include "CudaPoissonOpImpl.hpp"
} // namespace Poisson
} // namespace Elliptic
