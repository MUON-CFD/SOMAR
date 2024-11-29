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
 *  https://github.com/somarhub.
 ******************************************************************************/
#ifndef ___EllipticPoissonOp_H__INCLUDED___
#define ___EllipticPoissonOp_H__INCLUDED___

#include "AMRMGOperator.H"
#include "AnisotropicFluxRegister.H"
#include "BCTools.H"
#include "BdryIter.H"
#include "CFInterp.H"
#include "CornerCopier.H"
#include "LDFABOps.H"
#include "LevelGeometry.H"
#include "ProjectorParameters.H"

// #include "LepticOperator.H"
#include "LevelLepticSolver.H"

namespace Elliptic
{


/**
 * @class     PoissonOp
 * @brief     Defines the operator alpha*I + beta*L for the AMRMGSolver.
 * @author    Edward.Santilli@jefferson.edu
 * @warning   Homog Neum BCs are hardcoded for now.
 */
class PoissonOp : public AMRMGOperator<LevelData<FArrayBox>>
                , public LepticOperator
{
public:
    typedef LevelData<FArrayBox>                          StateType;
    typedef StateTraits<LevelData<FArrayBox>>::GridsType  GridsType;
    typedef StateTraits<LevelData<FArrayBox>>::CopierType CopierType;
    typedef MGOperator<LevelData<FArrayBox>>              MGOpType;
    typedef AMRMGOperator<LevelData<FArrayBox>>           AMRMGOpType;

    /// \brief Constructor for op at an AMR level, not in the MG depths.
    PoissonOp(const LevelGeometry&                       a_levGeo,
              const DisjointBoxLayout&                   a_fineGrids,
              const DisjointBoxLayout&                   a_crseGrids,
              const int                                  a_numComps,
              const std::shared_ptr<BCTools::BCFunction> a_bcFuncPtr,
              const Real                                 a_alpha   = 0.0,
              const Real                                 a_beta    = 1.0,
              const LevelData<FluxBox>*                  a_JgupPtr = nullptr);

    /// \brief Constructor for a LepticOperator at an AMR level, not in the MG
    ///  depths.
    /// \note This has very few features. It sets up an op with homogeneous
    ///  Neumann BCs.
    PoissonOp(const LevelGeometry&     a_levGeo,
              const DisjointBoxLayout& a_crseGrids,
              const int                a_numComps);

    /// Destructor for good measure.
    virtual ~PoissonOp();

    /// Use this to create a Helmholtz-type operator,
    /// Op[phi] = [alpha*J*I + beta*J*Laplacian](phi).
    /// This will call setupInvDiags and set m_hasNullSpace.
    virtual void
    setAlphaAndBeta(const Real a_alpha, const Real a_beta);


    // LepticOperator overrides ------------------------------------------------
    virtual const RealVect&
    getDomainLength() const override
    {
        return m_L;
    }

    // Also needed by MGOperator.
    // virtual const RealVect&
    // getDXi() const override
    // {
    //     return m_dXi;
    // }

    virtual const RealVect&
    getDXiCrse() const override
    {
        return m_amrCrseDXi;
    }

    virtual const LevelData<FluxBox>&
    getFCJgup() const override
    {
        return m_Jgup;
    }

    virtual IntVect
    getBlockFactor() const override
    {
        IntVect blockFactor = IntVect(D_DECL(4, 4, 4));
        for (int d = 0; d < SpaceDim; ++d) {
            if (m_activeDirs[d] == 0) blockFactor[d] = 1;
        }
        return blockFactor;
    }

    virtual std::shared_ptr<LevelOperator<LevelData<FArrayBox>>>
    createLevelOperator(const DisjointBoxLayout& a_grids) const override;

    virtual std::shared_ptr<MGOperator<LevelData<FArrayBox>>>
    createHorizontalMGOperator(const DisjointBoxLayout& a_hgrids) const override;


    // LinearOperator overrides ------------------------------------------------

    /// Sets BCs and fills all ghosts of a_phi.
    /// This is needed by the MG and AMR prolongation functions, but we include
    /// it here because it is likely to be used in applyOp. If you are not using
    /// ghosts or if you are not creating an MG or AMR operator, this can be a
    /// simple no-op. a_crsePhiPtr can be nullptr if not needed.
    virtual inline void
    applyBCs(LevelData<FArrayBox>&       a_phi,
             const LevelData<FArrayBox>* a_crsePhiPtr,
             const Real                  a_time,
             const bool                  a_homogPhysBCs,
             const bool                  a_homogCFIBCs) const
    {
        this->applyBCsBegin(a_phi, a_crsePhiPtr, a_time, a_homogPhysBCs, a_homogCFIBCs);
        this->applyBCsEnd(a_phi);
    }

    /// Starts async exchange and sets BCs.
    virtual void
    applyBCsBegin(LevelData<FArrayBox>&       a_phi,
                  const LevelData<FArrayBox>* a_crsePhiPtr,
                  const Real                  a_time,
                  const bool                  a_homogPhysBCs,
                  const bool                  a_homogCFIBCs) const;

    /// Ends async exchange. You MUST call applyBCsBegin first.
    virtual inline void
    applyBCsEnd(LevelData<FArrayBox>& a_phi) const
    {
        a_phi.exchangeEnd();
    }

    /// Computes lhs = L[phi], assuming the BCs have already been set.
    virtual void
    applyOpNoBCs(LevelData<FArrayBox>&       a_lhs,
                 const LevelData<FArrayBox>& a_phi) const;


    /// Computes lhs = L[phi].
    virtual void
    applyOp(LevelData<FArrayBox>&       a_lhs,
            LevelData<FArrayBox>&       a_phi,
            const LevelData<FArrayBox>* a_crsePhiPtr,
            const Real                  a_time,
            const bool                  a_homogPhysBCs,
            const bool                  a_homogCFIBCs) const
    {
        this->applyBCs(
            a_phi, a_crsePhiPtr, a_time, a_homogPhysBCs, a_homogCFIBCs);
        this->applyOpNoBCs(a_lhs, a_phi);
    }

    /// If this operator has a null space, a_phi needs to be projected.
    /// In other words, if L[e^1] = 0, and phi = c_1*e^1 + c_2*e^2 + ...,
    /// then this function should remove c_1*e^1.
    /// (Typically, you would just remove the average.)
    virtual void
    removeKernel(LevelData<FArrayBox>& a_phi) const;

    /// Checks if L[phi] = rhs is a solvable problem.
    /// For this op, we check if Sum[rhs] = 0.
    virtual bool
    levelEquationIsConsistent(const StateType& a_phi,
                              const StateType* a_crsePhiPtr,
                              const StateType& a_rhs,
                              const Real       a_time,
                              const bool       a_homogPhysBCs,
                              const bool       a_homogCFIBCs) const;

    /// Sets phi to a preconditioned value.
    /// Typically, this will be phi = rhs / invDiags, then relaxed, but you
    /// can ignore a_relaxIters if you have a better idea.
    virtual void
    preCond(LevelData<FArrayBox>&       a_phi,
            const LevelData<FArrayBox>& a_rhs,
            const Real                  a_time,
            const int                   a_relaxIters) const;

    /// Applies relaxation to the residual equation.
    /// If you don't plan to use this in preCond, or if you aren't creating
    /// an MG or AMR op, then this function can be a simple no-op.
    virtual void
    relax(LevelData<FArrayBox>&       a_cor,
          const LevelData<FArrayBox>& a_res,
          const Real                  a_time,
          const int                   a_iters) const;


    // MGOperator overrides ----------------------------------------------------

    /// Factory method.
    /// Allocate + define a coarsened version of *this.
    /// Deletion is left to the caller.
    ///
    /// a_refRatio specifies how much *coarsening* will be needed to get from
    /// this to the new operator.
    ///
    /// The new operator will only be used with homogeneous BCs.
    virtual MGOperator<LevelData<FArrayBox>>*
    newMGOperator(const IntVect& a_refRatio) const __malloc;

    /// What is the smallest Box we can use in MG?
    virtual IntVect
    minBoxSize() const
    {
        // If you change this, be sure to set to 1 in inactive dirs.
        return IntVect::Unit;
    }

    /// What grids is this op defined over?
    virtual const DisjointBoxLayout&
    getBoxes() const
    {
        return m_grids;
    }

    /// What is the resolution of this op's grids?
    virtual const RealVect&
    getDXi() const
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
               const MGOpType&  a_crseOp) const;

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
              const int       a_interpOrder) const;

    // AMRMGOperator overrides -------------------------------------------------

    /// Apply the AMR operator, including coarse-fine matching.
    /// The CFI BCs are always inhomogeneous in this function.
    virtual void
    AMROperator(LevelData<FArrayBox>&       a_LofPhi,
                LevelData<FArrayBox>&       a_phiFine,  // You may set ghosts.
                LevelData<FArrayBox>&       a_phi,      // You may set ghosts.
                const LevelData<FArrayBox>& a_phiCoarse,
                const IntVect&              a_fineRefRatio,
                const IntVect&              a_crseRefRatio,
                const Real                  a_time,
                const bool                  a_homogPhysBCs,
                const AMRMGOpType&          a_finerOp) const;

    /// Apply the AMR operator, including coarse-fine matching.
    /// The CFI BCs are always inhomogeneous in this function.
    /// assume no finer AMR level.
    virtual void
    AMROperatorNF(LevelData<FArrayBox>&       a_LofPhi,
                  LevelData<FArrayBox>&       a_phi,  // You may set ghosts.
                  const LevelData<FArrayBox>& a_phiCoarse,
                  const IntVect&              a_crseRefRatio,
                  const Real                  a_time,
                  const bool                  a_homogPhysBCs) const;

    /// Apply the AMR operator, including coarse-fine matching.
    /// Assume no coarser AMR level.
    virtual void
    AMROperatorNC(LevelData<FArrayBox>& a_LofPhi,
                  LevelData<FArrayBox>& a_phiFine,  // You may set ghosts.
                  LevelData<FArrayBox>& a_phi,      // You may set ghosts.
                  const IntVect&        a_fineRefRatio,
                  const Real            a_time,
                  const bool            a_homogPhysBCs,
                  const AMRMGOpType&    a_finerOp) const;

    /// \brief
    ///  Compute the p-Norm of a_x only on where not covered by a_fineRes.
    ///  If a_p = 0, then we compute an inf-norm.
    /// \details
    ///  This function must perform all MPI communication.
    ///  We want the norm on this level to be comparable to the norms on
    ///  other levels, so we choose powScale = dx*dy*dz and compute the norm as
    ///   |a_res|_p = ( Sum |a_res|^p * powScale ) ^(1/p).
    ///
    ///  Sorry if the name is confusing. It means we want the norm of a single
    ///  level, but that level may be covered by another level in an AMR solve.
    ///  If a_fineResPtr is nullptr, then we just need plain ol' norm().
    virtual Real
    AMRNormLevel(const StateType& a_res,
                 const StateType* a_fineResPtr,
                 const IntVect&   a_refRatio,
                 const int        a_p) const;

    // Miscellaneous methods ---------------------------------------------------

    /// Basic accessor
    virtual const ProblemDomain&
    getDomain() const
    {
        return m_domain;
    }

    /// If L[a_phi] = a divergence of fluxes, then refluxing will be needed.
    /// This function provides the fluxes needed by the refluxing scheme.
    /// a_phi's ghosts must be filled prior to call.
    virtual void
    getFlux(FArrayBox&       a_flux,
            const FArrayBox& a_phi,
            const Box&       a_fluxBox,
            const DataIndex& a_di,
            const int        a_fluxDir) const;

    /// Used by the AMR ops with a finer level.
    virtual void
    reflux(LevelData<FArrayBox>&       a_res,
           LevelData<FArrayBox>&       a_finePhi,  // Ghosts will be filled.
           const LevelData<FArrayBox>& a_phi,
           const IntVect&              a_fineRefRatio,
           const AMRMGOpType&          a_finerOp) const;

    /// This version is useful when you know the fluxes and don't need to
    /// recompute them.
    virtual void
    reflux(LevelData<FArrayBox>&     a_div,
           const LevelData<FluxBox>& a_flux,
           const LevelData<FluxBox>& a_fineFlux) const;

    /// Set phi BCs and compute Grad[phi].
    /// Since we need to set CFI BCs, levelGrad = compGrad.
    /// Both are included for completeness.
    virtual void
    levelGradient(LevelData<FluxBox>&         a_gradPhi,
                  LevelData<FArrayBox>&       a_phi,
                  const LevelData<FArrayBox>* a_crsePhiPtr,
                  const Real                  a_time,
                  const bool                  a_homogPhysBCs,
                  const bool                  a_homogCFIBCs) const;

    /// Set phi BCs and compute composite Grad[phi].
    virtual void
    compGradient(LevelData<FluxBox>&         a_gradPhi,
                 LevelData<FArrayBox>&       a_phi,
                 const LevelData<FArrayBox>* a_crsePhiPtr,
                 const Real                  a_time,
                 const bool                  a_homogPhysBCs,
                 const bool                  a_homogCFIBCs) const;

    /// Divergence of fluxes without flux corrections from the finer level.
    virtual void
    levelDivergence(LevelData<FArrayBox>&     a_div,
                    const LevelData<FluxBox>& a_flux) const;

    /// Composite divergence of fluxes with refluxing.
    virtual void
    compDivergence(LevelData<FArrayBox>&     a_div,
                   const LevelData<FluxBox>& a_flux,
                   const LevelData<FluxBox>* a_fineFluxPtr) const;

    // // LepticOperator overrides ------------------------------------------------
    // virtual Real
    // domainHeight() const
    // {
    //     return m_H;
    // }

    // virtual std::shared_ptr<const LepticOperator<StateType>>
    // createHorizontalOperator() const;

    // /// @brief Creates a flattened version of a data holder. This should be
    // /// called by a horizontal-only (flattened) LepticOperator.
    // /// @param a_hState The resulting horizontal (flattened) data holder.
    // /// @param a_state  The full-dimension (not flattened) data holder.
    // virtual void
    // createHorizontalState(StateType& a_hState, const StateType& a_state) const;

    // /// @brief Compute the vertical integral of a_Jrhs and place it in a_hInt.
    // /// @param a_hInt The horizontal (flat) data holder for the result.
    // /// @param a_Jrhs The full-dimension (not flat) data to integrate.
    // /// @details
    // ///  This can be called by any LepticOperator, flattened or not flattened.
    // virtual void
    // verticalIntegral(StateType& a_hInt, const StateType& a_Jrhs) const;

    // /// @brief Adds a_hInc(x,y)*a_scale to a_phi(x,y,z) for all z in the domain.
    // /// @param a_phi   The full-dimension (not flat) data to be incremented.
    // /// @param a_hInc  The horizontal (flat) data that serves as the increment.
    // /// @param a_scale A scaling factor for the increment.
    // virtual void
    // addVerticalExtrusion(StateType&       a_phi,
    //                      const StateType& a_hInc,
    //                      const Real       a_scale) const;

    // virtual void
    // solveVerticalODEs(StateType&       a_fphi,
    //                   const StateType& a_fJrhs,
    //                   const StateType& a_hUpperBC) const;

protected:
    /// Clone a_srcOp.
    /// Metric data will be aliased.
    PoissonOp(const PoissonOp& a_srcOp);

    /// Clones a_srcOp, but used a new set of grids.
    PoissonOp(const PoissonOp&         a_srcOp,
              const DisjointBoxLayout& a_newGrids);

    /// Defines this as a coarsened version of a_srcOp.
    /// a_crseGrids must be compatible with a_srcOp::m_grids.
    /// This constructor is used in MGSolvers to create coarsened ops.
    /// BCs will be homogeneous.
    PoissonOp(const PoissonOp&         a_srcOp,
              const DisjointBoxLayout& a_crseGrids,
              const IntVect&           a_refRatio);

    /// \brief Constructor for a horizontal-only op that can be used in a
    ///  LepticSolver at an AMR level, not in the MG depths.
    /// \note This has very few features. It sets up an op with homogeneous
    ///  Neumann BCs.
    PoissonOp(const GeoSourceInterface&                  a_geoSrc,
              const DisjointBoxLayout&                   a_hGrids,
              const DisjointBoxLayout&                   a_hCrseGrids,
              const int                                  a_numComps,
              const RealVect&                            a_L,
              const RealVect&                            a_dXi,
              const std::shared_ptr<BCTools::BCFunction> a_bcFuncPtr,
              const int                                  a_relaxMethod);

    /// Recomputes the Poisson op matrix elements. The diagonal is stored as
    /// 1/diag to speed up relaxation schemes.
    /// Must be called after setAlphaAndBeta, but before setting m_hasNullSpace.
    virtual void
    cacheMatrixElements();

    /// Checks if Op[ones] = zero.
    /// Must be called after setAlphaAndBeta and after setupInvDiags.
    /// This will not set m_hasNullSpace. Do that yourself after the call.
    virtual bool
    checkForNullSpace() const;

    // Relaxation methods ------------------------------------------------------
    /// Point-relaxation scheme used in Martin's thesis.
    /// cor^{i+1} = cor^{i} + lambda * (L[cor^{i}] - rhs)
    ///           = phi^{i} - lambda * res
    /// lambda = omega * min(dx)^2 / 4
    virtual void
    point_relax(LevelData<FArrayBox>&       a_cor,
                const LevelData<FArrayBox>& a_res,
                const Real                  a_time,
                const int                   a_iters) const;

    /// Jacobi relaxation.
    virtual void
    jacobi_relax(LevelData<FArrayBox>&       a_phi,
                 const LevelData<FArrayBox>& a_rhs,
                 const Real                  a_time,
                 const int                   a_iters) const;

    /// Red-Black Jacobi relaxation.
    virtual void
    jacobiRB_relax(LevelData<FArrayBox>&       a_phi,
                   const LevelData<FArrayBox>& a_rhs,
                   const Real                  a_time,
                   const int                   a_iters) const;

    /// Gauss-Seidel relaxation.
    virtual void
    gs_relax(LevelData<FArrayBox>&       a_phi,
             const LevelData<FArrayBox>& a_rhs,
             const Real                  a_time,
             const int                   a_iters) const;

    /// Red-Black Gauss-Seidel relaxation.
    virtual void
    gsrb_relax(LevelData<FArrayBox>&       a_phi,
               const LevelData<FArrayBox>& a_rhs,
               const Real                  a_time,
               const int                   a_iters) const;

    /// Red-Black Gauss-Seidel relaxation using a vertical line solver.
    virtual void
    vertLineGSRB_relax(LevelData<FArrayBox>&       a_phi,
                       const LevelData<FArrayBox>& a_rhs,
                       const Real                  a_time,
                       const int                   a_iters) const;


    // Member variables --------------------------------------------------------
#if CH_SPACEDIM == 2
    static constexpr IntVect s_hunit = IntVect(1, 0);
    static constexpr IntVect s_vunit = IntVect(0, 1);
#else
    static constexpr IntVect s_hunit = IntVect(1, 1, 0);
    static constexpr IntVect s_vunit = IntVect(0, 0, 1);
#endif

    const DisjointBoxLayout   m_crseAMRGrids;
    const GeoSourceInterface& m_geoSrc;
    const int                 m_whoMadeMe;  // Just for debugging.

    IntVect           m_activeDirs;
    ProblemDomain     m_domain;
    DisjointBoxLayout m_grids;
    RealVect          m_dXi;
    int               m_numComps;
    int               m_relaxMethod;

    // J, Jgup cache
    LevelData<FluxBox>   m_Jgup;
    LevelData<FArrayBox> m_J;

    // AMR stuff...
    RealVect                        m_amrCrseDXi;
    mutable AnisotropicFluxRegister m_fluxReg;

    // BC stuff...
    std::shared_ptr<BCTools::BCFunction> m_bcFuncPtr;
    mutable PhysBdryIter                 m_physBdryIter;
    mutable CFIIter                      m_cfiIter;
    CFInterp                             m_cfInterp;
    Copier                               m_exCopier;

    // Op = J * (alpha * I + beta * Laplacian)
    Real m_alpha;
    Real m_beta;

    // Cached matrix elements.
    // Warning:
    //  m_M is scaled by 1/J. (A true Laplacian)
    //  1/m_Dinv is the diag of J*(alpha*I + beta*Laplacian).
    FArrayBox            m_M[SpaceDim];  // comps: 0 = lower, 1 = upper
    LevelData<FArrayBox> m_Dinv;

    CornerCopier m_HOProlongCornerCopier;

    bool m_hasNullSpace;

    // For line relaxation
    LevelData<FArrayBox> m_vertTriDiagsLoBCs;
    LevelData<FArrayBox> m_vertTriDiagsHiBCs;

    // Leptic-specific stuff
    RealVect m_L;
};


};  // namespace Elliptic
#endif  //!___EllipticPoissonOp_H__INCLUDED___
