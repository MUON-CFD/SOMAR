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
#ifndef ___DiffusiveOp_H__INCLUDED___
#define ___DiffusiveOp_H__INCLUDED___

#include <memory>
#include "MGOperator.H"
#include "LDFABOps.H"
#include "LevelGeometry.H"
#include "BCTools.H"
#include "SOMAR_Constants.H"
#include "CFInterp.H"


namespace Elliptic {


class DiffusiveOp : public MGOperator<LevelData<FArrayBox>>
{
public:
    typedef LevelData<FArrayBox>                          StateType;
    typedef StateTraits<LevelData<FArrayBox>>::GridsType  GridsType;
    typedef StateTraits<LevelData<FArrayBox>>::CopierType CopierType;
    typedef MGOperator<LevelData<FArrayBox>>              MGOpType;

    /// Standard, single level constructor.
    /// Defines a Helmholtz-type operator,
    /// Op[phi] = [alpha*I + beta*D*(kappa + eddyKappa)*G](phi),
    /// where phi is cell-centered.
    DiffusiveOp(const LevelGeometry&                        a_levGeo,
                const DisjointBoxLayout*                    a_crseGridsPtr,
                const Real                                  a_alpha,
                const Real                                  a_beta,
                const std::vector<Real>&                    a_vKappa,
                const LevelData<FArrayBox>&                 a_eddyNu,
                const std::vector<Real>&                    a_vEddyPrandtl,
                const std::shared_ptr<BCTools::BCFunction>& a_bcFuncPtr,
                const int                                   a_numComps);

    ///
    virtual
    ~DiffusiveOp() {}

    /// Modifies the coefficients of this Helmholtz-type operator,
    /// Op[u] = [alpha*I + beta*D*nu*G](u).
    /// This will call cacheMatrixElements for you.
    virtual void
    setAlphaAndBeta(const Real a_alpha, const Real a_beta);

    ///
    virtual inline int
    numComps() const
    {
        CH_assert(m_resPtr);
        return m_resPtr->nComp();
    }

    // -------------------------------------------------------------------------
    /// \name LevelOperator overrides
    /// \{

    /// Sets BCs and fills all ghosts of the Cartesian (not stretched) a_phi.
    /// This is needed by the MG and AMR prolongation functions, but we include
    /// it here because it is likely to be used in applyOp. If you are not using
    /// ghosts or if you are not creating an MG or AMR operator, this can be a
    /// simple no-op. a_crsePhiPtr can be nullptr if not needed.
    virtual void
    applyBCs(StateType&       a_phi,
             const StateType* a_crsePhiPtr,
             const Real       a_time,
             const bool       a_homogPhysBCs,
             const bool       a_homogCFIBCs) const override;

    /// Computes lhs = L[phi].
    virtual void
    applyOp(StateType&       a_lhs,
            StateType&       a_phi,
            const StateType* a_crsePhiPtr,
            const Real       a_time,
            const bool       a_homogPhysBCs,
            const bool       a_homogCFIBCs) const override;

    /// If this operator has a null space, a_phi needs to be projected.
    /// In other words, if L[e^1] = 0, and phi = c_1*e^1 + c_2*e^2 + ...,
    /// then this function should remove c_1*e^1.
    /// (Typically, you would just remove the average.)
    virtual void
    removeKernel(StateType& /* a_phi */) const override
    {
        CH_verify(!RealCmp::isZero(m_alpha));
    }

    /// Sets phi to a preconditioned value.
    /// Typically, this will be phi = rhs / invDiags, then relaxed, but you
    /// can ignore a_relaxIters if you have a better idea.
    virtual void
    preCond(StateType&       a_phi,
            const StateType& a_rhs,
            const Real       a_time,
            const int        a_relaxIters) const override;

    /// Applies relaxation to the residual equation.
    /// If you don't plan to use this in preCond, or if you aren't creating
    /// an MG or AMR op, then this function can be a simple no-op.
    virtual void
    relax(StateType&       a_cor,
          const StateType& a_res,
          const Real       a_time,
          const int        a_relaxIters) const override;
    /// \}


    // MGOperator overrides ----------------------------------------------------

    /// Factory method.
    /// Allocate + define a coarsened version of *this.
    /// Deletion is left to the caller.
    ///
    /// a_refRatio specifies how much *coarsening* will be needed to get from
    /// this to the new operator.
    ///
    /// The new operator will only be used with homogeneous BCs.
    virtual MGOperator<StateType>*
    newMGOperator(const IntVect& a_refRatio) const override __malloc;

    /// What is the smallest Box we can use in MG?
    virtual IntVect
    minBoxSize() const override
    {
        // return IntVect(D_DECL(4,4,4));
        return IntVect(D_DECL(2,2,2));
        // return IntVect(D_DECL(1,1,1));
    }

    /// What grids is this op defined over?
    virtual const DisjointBoxLayout&
    getBoxes() const override
    {
        return m_grids;
    }

    /// What is the resolution of this op's grids?
    virtual const RealVect&
    getDXi() const override
    {
        return m_dXi;
    }

    /// Restrict to coarser MG depth: a_crseRes = I[h->2h](a_fineRes).
    /// This op is at the fine level.
    /// a_crseRes and a_fineRes should be defined over compatible grids.
    virtual void
    MGRestrict(StateType&       a_crseRes,
               const StateType& a_fineRes,
               const Real       a_time,
               const IntVect&   a_refRatio,
               const MGOpType&  a_crseOp) const override;

    /// Prolong (interpolate) a_crse to finer MG depth and add to a_fine:
    ///   a_finePhi += I[2h->h](a_crseCor).
    /// This op is at the fine level.
    /// a_crse and a_fine should be defined over compatible grids.
    ///
    /// We do not overwrite a_finePhi! We ADD the correction.
    ///
    /// If a_interpOrder > 1, then all ghosts of a_crseCor must be extrapolated,
    /// including edge and vertex ghosts.
    virtual void
    MGProlong(StateType&      a_finePhi,
              StateType&      a_crseCor,
              const Real      a_time,
              const IntVect&  a_refRatio,
              const MGOpType& a_crseOp,
              const int       a_interpOrder) const override;

protected:
    // Clone a_srcOp.
    // Metric data will be aliased.
    DiffusiveOp(const DiffusiveOp& a_srcOp);

    // Defines this as a coarsened version of a_srcOp.
    // a_crseGrids must be compatible with a_srcOp::m_grids.
    // This constructor is used in MGSolvers to create coarsened ops.
    // BCs will be homogeneous.
    DiffusiveOp(const DiffusiveOp&       a_srcOp,
                const DisjointBoxLayout& a_crseGrids,
                const IntVect&           a_refRatio);

    virtual void
    jacobi_relax(StateType&       a_cor,
                 const StateType& a_res,
                 const Real       a_time,
                 const int        a_relaxIters) const;

    virtual void
    gsrb_relax(StateType&       a_cor,
               const StateType& a_res,
               const Real       a_time,
               const int        a_relaxIters) const;

    // Member variables --------------------------------------------------------
    Real                      m_alpha;
    Real                      m_beta;

    ProblemDomain             m_domain;
    DisjointBoxLayout         m_grids;
    RealVect                  m_dXi;
    const GeoSourceInterface& m_geoSrc;

    const DisjointBoxLayout*  m_crseAMRGridsPtr;
    RealVect                  m_amrCrseDXi;

    std::shared_ptr<LevelData<FArrayBox>> m_Jptr;
    std::shared_ptr<LevelData<FluxBox>>   m_JgupPtr;

    std::shared_ptr<LevelData<FArrayBox>> m_invDiagsPtr;

    // BC stuff...
    std::shared_ptr<BCTools::BCFunction>  m_bcFuncPtr;
    mutable PhysBdryIter                  m_physBdryIter;
    mutable CFIIter                       m_cfiIter;
    std::shared_ptr<CFInterp>             m_cfInterpPtr;
    Copier                                m_exCopier;

    mutable std::shared_ptr<LevelData<FArrayBox>> m_resPtr;
};


}; // end namespace Elliptic

#endif //!___DiffusiveOp_H__INCLUDED___
