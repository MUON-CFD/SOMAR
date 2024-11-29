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
#ifndef ___FiniteDiff_H__INCLUDED___
#define ___FiniteDiff_H__INCLUDED___

#include "FArrayBox.H"
#include "LevelData.H"
#include "LevelGeometry.H"
#include "StaggeredFluxLD.H"
#include "Debug.H"


class FiniteDiff
{
public:
    // Non-static declarations...

    // Default constructor
    FiniteDiff();

    // Full constructor
    FiniteDiff(const LevelGeometry& a_levGeo);

    // Full virtual constructor
    virtual void
    define(const LevelGeometry& a_levGeo);

    // Disallow copy and assignment
    FiniteDiff& operator=(const FiniteDiff&) = delete;
    FiniteDiff(const FiniteDiff&) = delete;

    // Destructor
    virtual
    ~FiniteDiff();

    // "Destructor"
    virtual void
    clear();

    // Basic accessors
    inline virtual bool
    isDefined() const {return (m_levGeoPtr != nullptr);}

    inline virtual const LevelGeometry&
    getLevGeo() const {return *m_levGeoPtr;}

    /// @brief Computes the staggered derivative of a scalar onto a_pdBox.
    ///        This function does not fill ghosts or bias the stencil.
    /// @param a_pd       The dest holder.
    /// @param a_pdComp   The comp of a_pd that will be filled.
    /// @param a_pdBox    Where to compute the derivs.
    ///                   Must have same centering as a_pd
    /// @param a_phi      The source data. Ghosts must be filled if needed.
    ///                   Must have the opposite centering of a_pd in a_derivDir.
    /// @param a_phiComp  The comp of a_phi that will be differentiated.
    /// @param a_derivDir The differencing direction.
    /// @param a_dXi      The grid spacing in a_derivDir.
    /// @param a_addToPD  If true, we will add to a_pd rather than overwrite it.
    static void
    partialD(FArrayBox&       a_pd,
             const int        a_pdComp,
             const Box&       a_pdBox,
             const FArrayBox& a_phi,
             const int        a_phiComp,
             const int        a_derivDir,
             const Real       a_dXi,
             const bool       a_addToPD = false);

    /// @brief Computes the colocated derivative of a scalar onto a_pdBox.
    ///        This function avoids ghosts by biasing the stencil.
    ///        The resulting slopes are second order.
    /// @param a_pd          The dest holder.
    /// @param a_pdComp      The comp of a_pd that will be filled.
    /// @param a_phi         The source data.
    /// @param a_phiComp     The comp of a_phi that will be differentiated.
    /// @param a_validPhiBox We will only reference src data from this box.
    /// @param a_derivDir    The differencing direction.
    /// @param a_dXi         The grid spacing in a_derivDir.
    static void
    computeSlopes(FArrayBox&       a_pd,
                  const int        a_pdComp,
                  const FArrayBox& a_phi,
                  const int        a_phiComp,
                  const Box&       a_validPhiBox,
                  const int        a_derivDir,
                  const Real       a_dXi);

    /// @brief Computes the second differences of a_phi and adds it to a_pd.
    ///        This stays within a_validPhibox by using biased stencils.
    /// @param a_pd          The dest holder.
    /// @param a_pdComp      The comp of a_pd that will be filled.
    /// @param a_phi         The source data.
    /// @param a_phiComp     The comp of a_phi that will be differentiated.
    /// @param a_validPhiBox We will only reference src data from this box.
    /// @param a_derivDir    The differencing direction.
    /// @param a_scale       For a true second derivative, set a_scale to 1/dx^2.
    static void
    addSecondDifference(FArrayBox&       a_pd,
                        const int        a_pdComp,
                        const FArrayBox& a_phi,
                        const int        a_phiComp,
                        const Box&       a_validPhiBox,
                        const int        a_derivDir,
                        const Real       a_scale);

    /// @brief Computes the Laplacian of a_phi and adds it to a_pd.
    ///        This stays within a_validPhibox by using biased stencils.
    /// @param a_pd          The dest holder.
    /// @param a_pdComp      The comp of a_pd that will be filled.
    /// @param a_phi         The source data.
    /// @param a_phiComp     The comp of a_phi that will be differentiated.
    /// @param a_validPhiBox We will only reference src data from this box.
    /// @param a_scale       For a true second derivative, set a_scale to 1/dx^2 in each dir.
    ///                      We do not compute in dirs with a_scale[dir] = 0.
    static void
    addLaplacian(FArrayBox&       a_pd,
                 const int        a_pdComp,
                 const FArrayBox& a_phi,
                 const int        a_phiComp,
                 const Box&       a_validPhiBox,
                 const RealVect&  a_scale);


    /// \{
    /// \name Level operators

    // Computes the FC, single-level gradient of a CC scalar.
    //
    // This function assumes:
    //   1. All BCs have been set.
    //   2. *this, a_gradPhi, and a_phi are defined over the same grids.
    //   3. a_gradPhi and a_phi have the same number of components.
    virtual void
    levelGradientMAC(LevelData<FluxBox>&         a_gradPhi,
                     const LevelData<FArrayBox>& a_phi) const;

    // Computes the CC, single-level divergence of FC fluxes.
    //
    // This function assumes:
    //   1. All BCs have been set, even on the finer data (if it exists.)
    //   2. The fluxes are scaled by J.
    //   3. a_div and a_flux are defined over the same grids.
    //
    // a_scaleByJinv: A true divergence will scale the result by 1/J, but
    //    there are some algorithms that avoid this for efficiency reasons.
    //
    // a_scale: These will be multiplied into each partial derivative.
    //    For example, if you wish to compute a horizontal divergence in 3D,
    //    set this to (1.0, 1.0, 0.0).
    virtual void
    levelDivergenceMAC(LevelData<FArrayBox>&     a_div,
                       const LevelData<FluxBox>& a_flux,
                       const bool                a_scaleByJinv = true,
                       const RealVect&           a_scale = RealVect::Unit) const;

    virtual void
    levelDivergenceMAC_4thOrder(LevelData<FArrayBox>&     a_div,
                                const LevelData<FluxBox>& a_flux,
                                const bool                a_scaleByJinv = true,
                                const RealVect& a_scale = RealVect::Unit) const;

    // Take the gradient of staggered velocity vectors.
    // a_grad[i][a] = \partial_i u^a,
    // where i is curvilinear-based and a is Cartesian-based.
    virtual void
    levelVectorGradient(StaggeredFluxLD&          a_grad,
                        const LevelData<FluxBox>& a_cartVect) const;

    /// Compute the staggered divergence of fluxes.
    /// a_div[a] = Sum_i \partial_i F^{ia},
    /// where i is curvilinear-based and a is Cartesian-based.
    /// If a_addToDiv == true,  we will compute a_div += scale * divergence.
    /// If a_addToDiv == false, we will compute a_div  = scale * divergence.
    virtual void
    levelVectorDivergence(LevelData<FluxBox>&    a_div,
                          const StaggeredFluxLD& a_flux,
                          const bool             a_scaleByJinv = true,
                          const bool             a_addToDiv    = false,
                          const Real             a_scale       = 1.0) const;

    virtual void
    levelVectorDivergence_4thOrder(LevelData<FluxBox>&    a_div,
                                   const StaggeredFluxLD& a_flux,
                                   const bool             a_scaleByJinv = true,
                                   const bool             a_addToDiv    = false,
                                   const Real             a_scale = 1.0) const;
    /// \}


    /// \{
    /// \name Composite operators

    // Computes the FC gradient of a CC scalar.
    //
    // This function assumes:
    //   1. All BCs have been set.
    //   2. *this, a_gradPhi, and a_phi are defined over the same grids.
    //   3. a_gradPhi and a_phi have the same number of components.
    virtual void
    compGradientMAC(LevelData<FluxBox>&         a_gradPhi,
                    const LevelData<FArrayBox>& a_phi) const;

    // Computes the CC divergence of FC fluxes.
    //
    // This function assumes:
    //   1. All BCs have been set, even on the finer data (if it exists.)
    //   2. The fluxes are scaled by J.
    //   3. a_div and a_flux are defined over the same grids.
    //
    // a_scaleByJinv: A true divergence will scale the result by 1/J, but
    //    there are some algorithms that avoid this for efficiency reasons.
    //
    // a_scale: These will be multiplied into each partial derivative.
    //    For example, if you wish to compute a horizontal divergence in 3D,
    //    set this to (1.0, 1.0, 0.0).
    virtual void
    compDivergenceMAC(LevelData<FArrayBox>&     a_div,
                      const LevelData<FluxBox>& a_flux,
                      const LevelData<FluxBox>& a_fineFlux,
                      const bool                a_scaleByJinv = true,
                      const RealVect&           a_scale = RealVect::Unit) const;

    /// \}

protected:
    // Member variables.
    const LevelGeometry* m_levGeoPtr;
};



#endif //!___FiniteDiff_H__INCLUDED___
