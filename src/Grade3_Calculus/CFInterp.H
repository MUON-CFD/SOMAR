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
#ifndef ___CFInterp_H__INCLUDED___
#define ___CFInterp_H__INCLUDED___

#include "LevelData.H"
#include "RealVect.H"
#include "Tuple.H"
#include "LevelGeometry.H"
#include "BdryIter.H"


// #define ALLOW_DIVFREEINTERP


/**
 * \class   CFInterp
 * \brief   Interpolates data across grids at different levels of refinement.
 * \date    2019/11/02
 * \author  Edward.Santilli@jefferson.edu
 */
class CFInterp
{
private:
    // Sets the default parameter values.
    static const bool useLimiter = true;

public:
    /// \name Construction / destruction
    /*\{*/

    /**
     * \brief Weak constructor. Leaves object unusable.
     */
    CFInterp();

    /**
     * \brief      Full constructor. Avoids use of a_levGeoPtr->getCoarserPtr().
     *
     * \param[in]  a_levGeoPtr  Contains all of the fine level's geometric info.
     * \param[in]  a_userCrseGrids  The coarse level structure.
     */
    CFInterp(const LevelGeometry&     a_levGeo,
             const DisjointBoxLayout& a_userCrseGrids);

    /**
     * \brief      Destructor
     */
    virtual ~CFInterp();

    /**
     * \brief      Full constructor. Avoids use of a_levGeoPtr->getCoarserPtr().
     *
     * \param[in]  a_levGeoPtr  Contains all of the fine level's geometric info.
     * \param[in]  a_userCrseGrids  The coarse level structure.
     */
    virtual void
    define(const LevelGeometry&     a_levGeo,
           const DisjointBoxLayout& a_userCrseGrids);

    /**
     * \brief      Full constructor. Does not use levGeos.
     *
     * \param[in]  a_grids              The fine level's grids.
     * \param[in]  a_dXi                The fine level's grid spacing.
     * \param[in]  a_userCrseGrids      The coarse level's grids.
     */
    virtual void
    define(const DisjointBoxLayout& a_grids,
           const RealVect&          a_dXi,
           const DisjointBoxLayout& a_userCrseGrids);

    /**
     * \brief      Frees memory. Called by the destructor.
     */
    virtual void
    undefine();

    /*\}*/

    /// \name Basic accessors
    /*\{*/

    /**
     * \brief      Determines if this object is ready for use.
     *
     * \return     True if defined, False otherwise.
     */
    virtual inline bool
    isDefined() const;

    /**
     * \brief      Gets the ProblemDomain.
     *
     * \return     The ProblemDomain.
     */
    virtual inline const ProblemDomain&
    getDomain() const;

    /**
     * \brief      Gets the ProblemDomain's Box.
     *
     * \return     The ProblemDomain's Box.
     */
    virtual inline const Box&
    getDomainBox() const;

    /**
     * \brief      Gets this (fine) level's grids.
     *
     * \return     The DisjointBoxLayout.
     */
    virtual inline const DisjointBoxLayout&
    getBoxes() const;

    /**
     * \brief      Gets the user's coarse level's grids.
     *
     * \return     The DisjointBoxLayout.
     */
    virtual inline const DisjointBoxLayout&
    getCrseBoxes() const;

    /**
     * \brief      Gets the refinement ratio between this and the coarser level.
     *
     * \return     The refinement ratio.
     */
    virtual inline const IntVect&
    getRefRatio() const;

    /*\}*/

    /// \name CC interpolators
    /*\{*/

    /**
     * \brief      Interpolate at CF interface and extrapolate to first ghost
     *             layer.
     *
     * \param      a_fine     The fine level data.
     * \param[in]  a_crse     The coarse level data.
     * \param[in]  a_isHomog  If true, then a_crse can be undefined.
     *
     * \warning
     *  This interpolator does not need J-weighting. Do NOT send in J*data.
     */
    virtual void
    interpAtCFI(LevelData<FArrayBox>&       a_fine,
                const LevelData<FArrayBox>& a_crse,
                const bool                  a_isHomog = false) const;

    /**
     * \brief      Interpolate at CF interface and extrapolate to first ghost
     *             layer.
     *
     * \param      a_fine  The fine level data.
     *
     * \details
     *  Interpolate at CF interface and extrapolate to first ghost layer.
     *  Uses zeros as coarse level data.
     * \warning
     *  This interpolator does not need J-weighting. Do NOT send in J*fine.
     */
    virtual void
    homogInterpAtCFI(LevelData<FArrayBox>& a_fine) const;

    /**
     * \brief      Interpolate at CF interface and extrapolate to first ghost
     *             layer. This is a static version useful for MG.
     *
     * \param      a_fine     The data whose ghosts need fillin'.
     * \param[in]  a_fineDXi  a_fine's resolution.
     * \param[in]  a_crseDXi  The resolution of the homogeneous data (zeros).
     * \param      a_cfiIter  A pre-defined CFI iterator.
     *
     * \details
     *  Extrapolate at CF interface to the first ghost layer. This is a static
     *  version useful for MG. When using in an AMRMG solver, for example when
     *  solving Poisson's eq on AMR level 1, a_crseDXi is the dXi on AMR level
     *  0, MG depth 0 and a_fineDXi is the dXi on AMR level 1, MG depth d. Yes,
     *  this means a_fineDXi may actually be coarser than a_crseDXi.
     *
     * \warning
     *  This interpolator does not need J-weighting. Do NOT send in J*fine.
     */
    static void
    homogInterpAtCFI(LevelData<FArrayBox>& a_fine,
                     const RealVect&       a_fineDXi,
                     const RealVect&       a_crseDXi,
                     CFIIter&              a_cfiIter);

    /**
     * \brief      Interpolates data from coarse grids to fine grids.
     *
     * \param[out] a_fine           The fine level data.
     * \param[in]  a_crse           The coarse level data.
     * \param[in]  a_interpOrder    0 = constant (injection), 1 = linear
     * \param[in]  a_doSlopeLimiter Use Chombo's multidimensional limiter?
     *
     * \warning
     *  This interpolator does not need J-weighting. Do NOT send in J*crse.
     */
    virtual void
    refine(LevelData<FArrayBox>&       a_fine,
           const LevelData<FArrayBox>& a_crse,
           const int                   a_interpOrder    = 1,
           const bool                  a_doSlopeLimiter = useLimiter) const;

    /**
     * \brief      Interpolates data from coarse grids to fine grids.
     *
     * \param[out] a_fine         The fine level data.
     * \param[in]  a_crse         The coarse level data.
     * \param[in]  a_interpOrder  Can be 0 (injection) or 1 (linear).
     *
     * \details
     *  a_crse's grids must be a coarsened version of a_fine's grids.
     * \warning
     *  This interpolator does not need J-weighting. Do NOT send in J*crse.
     */
    virtual void
    localRefine(LevelData<FArrayBox>&       a_fine,
                const LevelData<FArrayBox>& a_crse,
                const int                   a_interpOrder    = 1,
                const bool                  a_doSlopeLimiter = useLimiter) const;

    /**
     * \brief      Interpolates data from coarse grids to fine grids.
     *
     * \param[out] a_fineFAB        The fine level's data.
     * \param[in]  a_crseFAB        The coarse level's data.
     * \param[in]  a_crseInterpBox  The region that will be refined.
     * \param[in]  a_interpOrder    Can be 0 (injection) or 1 (linear).
     * \param[in]  a_doSlopeLimiter Use Chombo's multi-dimensional limiter?
     *
     * \details
     *  a_fineFAB's Box must be a refined version of a_crseFAB's Box.
     * \warning
     *  This interpolator does not need J-weighting. Do NOT send in J*data.
     */
    virtual void
    localRefine(FArrayBox&       a_fineFAB,
                const FArrayBox& a_crseFAB,
                const Box&       a_crseInterpBox,
                const int        a_interpOrder    = 1,
                const bool       a_doSlopeLimiter = useLimiter) const;


    ///
    static void
    localRefine(FArrayBox&       a_fineFAB,
                const FArrayBox& a_crseFAB,
                const Box&       a_crseInterpBox,
                const Box&       a_crseCentralBox,
                const IntVect&   a_refRatio,
                const int        a_interpOrder    = 1,
                const bool       a_doSlopeLimiter = useLimiter);

    /**
     * \brief      Interpolates data from fine grids to coarse grids.
     *
     * \param[out] a_crse          The coarse level data.
     * \param[in]  a_fine          The fine level data.
     * \param[in]  a_doHarmonicAvg Set true to do harmonic averaging.
     * \param[in]  a_fineJPtr      If provided, this will do a J-weighted avg.
     */
    virtual void
    coarsen(LevelData<FArrayBox>&       a_crse,
            const LevelData<FArrayBox>& a_fine,
            const bool                  a_doHarmonicAvg,
            const LevelData<FArrayBox>* a_fineJPtr) const;

    /**
     * \brief      Interpolates data from fine grids to compatible coarse grids.
     *
     * \param[out] a_localCrse     The coarse data on fine-compatible grids.
     * \param[in]  a_fine          The fine level data.
     * \param[in]  a_doHarmonicAvg Set true to do harmonic averaging.
     * \param[in]  a_fineJPtr      If provided, this will do a J-weighted avg.
     *
     * \details
     *   a_localCrse's grids must be a coarsened version of a_fine's grids.
     */
    static void
    localCoarsen(LevelData<FArrayBox>&       a_localCrse,
                 const LevelData<FArrayBox>& a_fine,
                 const bool                  a_doHarmonicAvg,
                 const LevelData<FArrayBox>* a_fineJPtr);

    /**
     * \brief      Interpolates data from a fine patch to a coarse patch.
     *
     * \param[out] a_crseFAB       The CC fine level data.
     * \param[in]  a_fineFAB       The CC coarse level data.
     * \param[in]  a_crseInterpBox The data's destination.
     * \param[in]  a_refRatio      The refinement ratio between crse and fine.
     * \param[in]  a_doHarmonicAvg Set true to do harmonic averaging.
     * \param[in]  a_fineJFABPtr   If provided, this will do a J-weighted avg.
     */
    static void
    localCoarsen(FArrayBox&       a_crseFAB,
                 const FArrayBox& a_fineFAB,
                 const Box&       a_crseInterpBox,
                 const IntVect&   a_refRatio,
                 const bool       a_doHarmonicAvg,
                 const FArrayBox* a_fineJFABPtr);

    /*\}*/

    /// \name FC interpolators
    /*\{*/

    /**
     * \brief      Interpolate ghost faces at CF interface.
     *
     * \param[inout] a_fineAdvVel   The fine level, advecting vel, J*u^i.
     * \param[in]    a_crseAdvVel   The crse level, advecting vel, J*u^i.
     * \param[in]    a_isHomog      If true, then a_crse can be undefined.
     * \param[in]    a_fineLevGeo   For accessing metric info.
     *
     * \warning
     *  Be sure to set all other ghosts AFTER calling this function!
     */
    virtual void
    interpGhostsAtCFI(LevelData<FluxBox>&       a_fineAdvVel,
                      const LevelData<FluxBox>& a_crseAdvVel,
                      const bool                a_isHomog) const;

#ifdef ALLOW_DIVFREEINTERP
    /**
     * \brief      Interpolates divergence-free velocity data from the coarse
     *             level to the fine level.
     *
     * \param[out] a_fine     The fine level data.
     * \param[in]  a_crse     The coarse level data.
     * \param[in]  a_fineParentPtr If provided, the scheme will use this data
     *                             at the fine parent faces. This can be
     *                             defined on different grids than a_fine.
     *                             Parents will be copied wherever possible.
     *
     * \pre        The coarse level data must be divergence-free.
     * \pre        a_fine and a_crse MUST be scaled by J.
     *
     * \details
     *  Interpolates velocity data from the coarse level to the fine level.
     *  This will fill all valid fine-level faces, but will not fill ghosts.
     *  The interpolated data will be divergence-free to machine precision.
     *  Also, the vorticity in the xy-plane and upwelling/downwelling will be
     *  approximately preserved.
     *
     *  If you supply a_fineDataPtr, and it overlaps with a_fine's grids,
     *  then its data will be deemed more accurate than an interpolation of
     *  a_crse. In the end, it will be copied to a_fine in the overlap
     *  (using the copier, if provided).
     */
    virtual void
    divFreeInterp(LevelData<FluxBox>&       a_fine,
                  const LevelData<FluxBox>& a_crse,
                  const LevelData<FluxBox>* a_fineDataPtr = nullptr,
                  const Copier*             a_fineDataCopierPtr = nullptr) const;
#endif // ALLOW_DIVFREEINTERP

    /**
     * \brief      Interpolate at CF interface and extrapolate to first ghost
     *             layer. This is a static version useful for MG.
     *
     * \param      a_fine     The data whose ghosts need fillin'.
     * \param[in]  a_fineDXi  a_fine's resolution.
     * \param[in]  a_crseDXi  The resolution of the homogeneous data (zeros).
     * \param      a_cfRegion Tells us the locations of the fine FC ghosts.
     *
     * \details
     *  Extrapolate at CF interface to the first ghost layer. This is a static
     *  version useful for MG. When using in an AMRMG solver, for example when
     *  solving Poisson's eq on AMR level 1, a_crseDXi is the dXi on AMR level
     *  0, MG depth 0 and a_fineDXi is the dXi on AMR level 1, MG depth d. Yes,
     *  this meand a_fineDXi may actually be coarser than a_crseDXi.
     *
     * \warning
     *  This interpolator does not need J-weighting. Do NOT send in J*fine.
     *  Also, this function will not fill the FC data directly at the CFI.
     */
    static void
    homogInterpAtCFI(LevelData<FluxBox>&                   a_fine,
                     const RealVect&                       a_fineDXi,
                     const RealVect&                       a_crseDXi,
                     const std::array<CFRegion, SpaceDim>& a_cfRegion);

    /**
     * \brief      Interpolates FC data from coarse grids to fine grids.
     *
     * \param[out] a_fineAdvVel     The fine level advecting velocity.
     * \param[in]  a_crseAdvVel     The coarse level advecting velocity.
     */
    virtual void
    refine(LevelData<FluxBox>&       a_fineAdvVel,
           const LevelData<FluxBox>& a_crseAdvVel) const;

    /**
     * \brief      Interpolates FC data from coarse grids to fine grids.
     *
     * \param[out] a_fineAdvVel     The fine level advecting velocity.
     * \param[in]  a_crseAdvVel     The coarse level advecting velocity.
     *
     * \note
     *  a_crseAdvVel's grids must be a coarsened version of a_fineAdvVel's grids.
     */
    virtual void
    localRefine(LevelData<FluxBox>&       a_fineAdvVel,
                const LevelData<FluxBox>& a_crseAdvVel) const;

    /**
     * \brief      Interpolates FC data from a coarse patch to a fine patch.
     *
     * \param[out] a_fineAdvVelFB     The fine level's data.
     * \param[in]  a_crseAdvVelFB     The coarse level's data.
     * \param[in]  a_ccCrseInterpBox  The crse region to be interpolated.
     *
     * \note
     *  a_fineAdvVelFB's Box must be a refined version of a_crseAdvVelFB's Box.
     */
    virtual void
    localRefine(FluxBox&       a_fineAdvVelFB,
                const FluxBox& a_crseAdvVelFB,
                const Box&     a_ccCrseeInterpBox) const;

    /**
     * \brief      Interpolates FC data from a coarse patch to a fine patch.
     *
     * \param[out] a_fineAdvVelFAB    The fine level's data.
     * \param[in]  a_crseAdvVelFAB    The coarse level's data.
     * \param[in]  a_fcCrseInterpBox  The coarse region to be interpolated.
     * \param[in]  a_fcDir            The FABs and Boxes must be FC in this dir.
     * \param[in]  a_refRatio         The refinement ratio.
     *
     * \note
     * a_fineAdvVelFAB's Box must be a refined version of a_crseAdvVelFAB's Box.
     */
    static void
    localRefineFace(FArrayBox&       a_fineAdvVelFAB,
                    const FArrayBox& a_crseAdvVelFAB,
                    const Box&       a_fcCrseInterpBox,
                    const int        a_fcDir,
                    const IntVect&   a_refRatio);

    /**
     * \brief    Restrict data to coarser grids by averaging along parent faces.
     *
     * \param[out] a_crse     The coarse level data.
     * \param[in]  a_fine     The fine level data.
     * \param[in]  a_fineJPtr The area weights at each fine face.
     *
     * \details
     *  Restrict to compatible coarser grids by averaging along parent faces.
     *  Data that is not at a parent face is not considered, which makes
     *  this similar to injection in the face-normal directions.
     */
    virtual void
    coarsen(LevelData<FluxBox>&       a_crse,
            const LevelData<FluxBox>& a_fine,
            const LevelData<FluxBox>* a_fineJPtr) const;

    /**
     * \brief    Restrict data to coarser grids by averaging along parent faces.
     *
     * \param[out] a_localCrse  The coarse data on fine-compatible grids.
     * \param[in]  a_fine       The fine level data.
     * \param[in]  a_fineJPtr   The area weights at each fine face. Can be NULL.
     *
     * \details
     *  Restrict to coarser grids by averaging along parent faces.
     *  a_localCrse's grids must be a coarsened version of a_fine's grids.
     *  Data that is not at a parent face is not considered, which makes
     *  this similar to injection in the face-normal directions.
     */
    static void
    localCoarsen(LevelData<FluxBox>&       a_localCrse,
                 const LevelData<FluxBox>& a_fine,
                 const LevelData<FluxBox>* a_fineJPtr);

    /**
     * \brief  Restrict data to a coarser patch by averaging along parent faces.
     *
     * \param[out] a_crseFB       The coarse data.
     * \param[in]  a_fineFB       The fine data.
     * \param[in]  a_crseCCValid  The region that will be coarsened.
     * \param[in]  a_refRatio     The refinement between crse and fine.
     * \param[in]  a_fineJPtr   The area weights at each fine face. Can be NULL.
     *
     * \details
     *  Restrict to a coarser grid by averaging along parent faces.
     *  Data that is not at a parent face is not considered, which makes this
     *  similar to injection in the face-normal directions. This assumes
     *  a_crseFB is a coarsened version of a_fineFB.
     */
    static void
    localCoarsen(FluxBox&       a_crseFB,
                 const FluxBox& a_fineFB,
                 const Box&     a_crseCCValid,
                 const IntVect& a_refRatio,
                 const FluxBox* a_fineJPtr);

    /**
     * \brief  Restrict data to a coarser patch by averaging along parent faces.
     *
     * \param[out] a_crseFAB      The coarse data.
     * \param[in]  a_fineFAB      The fine data.
     * \param[in]  a_crseFCValid  The region that will be coarsened.
     * \param[in]  a_fcDir        The FC-normal direction.
     * \param[in]  a_refRatio     The refinement between crse and fine.
     * \param[in]  a_fineJPtr   The area weights at each fine face. Can be NULL.
     *
     * \details
     *  Restrict to a coarser grid by averaging along parent faces.
     *  Data that is not at a parent face is not considered, which makes this
     *  similar to injection in the face-normal directions. This assumes
     *  a_crseFB is a coarsened version of a_fineFB.
     */
    static void
    localCoarsenFace(FArrayBox&       a_crseFAB,
                     const FArrayBox& a_fineFAB,
                     const Box&       a_crseFCValid,
                     const int        a_fcDir,
                     const IntVect&   a_refRatio,
                     const FArrayBox* a_fineJPtr);

    /*\}*/

    /// \name EC interpolators
    /// \{

    /**
     * \brief  Restrict data to a coarser patch by averaging along parent edges.
     *
     * \param[out] a_crseFAB      The coarse data.
     * \param[in]  a_fineFAB      The fine data.
     * \param[in]  a_crseFCValid  The region that will be coarsened.
     * \param[in]  a_ecDir        The EC-normal direction.
     * \param[in]  a_refRatio     The refinement between crse and fine.
     * \param[in]  a_fineJPtr   The area weights at each fine face. Can be NULL.
     *
     * \details
     *  Restrict to a coarser grid by averaging along parent edges.
     *  Data that is not at a parent edge is not considered, which makes this
     *  similar to injection in the edge-ortho directions. This assumes
     *  a_crseFB is a coarsened version of a_fineFB.
     */
    static void
    localCoarsenEdge(FArrayBox&       a_crse,
                     const FArrayBox& a_fine,
                     const Box&       a_crseECValid,
                     const int        a_ecDir,
                     const IntVect&   a_refRatio,
                     const FArrayBox* a_fineJPtr);
    /// \}

    /// \name NC interpolators
    /// \{

    /**
     * \brief  Restrict data to a coarser patch by injection from parent nodes.
     *
     * \param[out] a_crseFAB      The coarse data.
     * \param[in]  a_fineFAB      The fine data.
     * \param[in]  a_crseNCValid  The region that will be coarsened.
     * \param[in]  a_refRatio     The refinement between crse and fine.
     *
     * \details
     *  Restrict data to a coarser patch by injection from parent nodes.
     *  Data that is not at a parent node is not considered.
     *  This assumes a_crseFB is a coarsened version of a_fineFB.
     */
    static void
    localCoarsenNode(FArrayBox&       a_crse,
                     const FArrayBox& a_fine,
                     const Box&       a_crseNCValid,
                     const IntVect&   a_refRatio);
    /// \}


    /**
     * @brief      Ensures all coarse data is identical at overlapping faces.
     *
     * @param      a_crse  The coarse level data.
     *
     * \details
     *  Static utility.
     *  After averaging down, we need to ensure consistency under the CFI. If a
     *  coarse level exchange boundary coincides with the CFI, then the two
     *  adjacent coarse patches will contain different data in their overlap --
     *  one patch will have newly averaged (valid) data while the other patch
     *  will have old, unaveraged (invalid) data. This function copies valid
     *  data to the invalid data on the adjacent patch.
     */
    static void
    validateAtCoarseCFI(LevelData<FluxBox>&      a_crse,
                        const DisjointBoxLayout& a_fineGrids,
                        const IntVect&           a_refRatio) __deprecated;

protected:

    /**
     * \brief      Prepares a_userCrseData to be interpolated.
     *
     * \param[out] a_localCrseData  Coarse data on fine-level compatible grids.
     * \param[in]  a_userCrseData   Coarse data as provided by the caller.
     *
     * \pre
     *  All ghost data must be filled prior to call.
     *
     * \details
     *  All valid and ghost data will be copied.
     *  Result will be on grids that are compatible with the fine level's grids.
     */
    virtual void
    localizeCrseData(LevelData<FArrayBox>&       a_localCrseData,
                     const LevelData<FArrayBox>& a_userCrseData) const;

    /**
     * \brief      Prepares a_userCrseData to be interpolated.
     *
     * \param[out] a_localCrseData  Coarse data on fine-level compatible grids.
     * \param[in]  a_userCrseData   Coarse data as provided by the caller.
     *
     * \pre
     *  All ghost data must be filled prior to call.
     *
     * \details
     *  All valid and ghost data will be copied.
     *  Result will be on grids that are compatible with the fine level's grids.
     */
    virtual void
    localizeCrseData(LevelData<FluxBox>&       a_localCrseData,
                     const LevelData<FluxBox>& a_userCrseData) const;

    /**
     * \brief      1D, conservative, quadratic interpolation for use along
     *             parent cell edges in 2D.
     *
     * \param[out] a_fineFAB    Fine-level interpolated data
     * \param[in]  a_crseFAB    Coarse-level data to be interpolated.
     * \param[in]  a_crseBox    The region to be interpolated.
     * \param[in]  a_interpDir  Interpolation direction.
     *                          Tangent to the parent edges.
     * \param[in]  a_refRatio   The refinement ratio
     *
     * \pre
     *  This function assumes all crseFAB ghosts have been filled.
     *
     * \details
     *  - crseBox must have the same centering as the FABs.
     *  - crseBox and the FABs must be cell-centered in interpDir.
     *  - All fine data above crseBox will be filled.
     *  - crseFAB must be large enough for the stencil (one ghost layer).
     *
     *  This stencil comes from setting u = F'(x) with ul, uc, and ur matched.
     *  This is not quite the same as the 3D version!
     *  Result has an O(dx^2) error.
     */
    static void
    interpAlongParentEdgesIn2D(FArrayBox&       a_fineFAB,
                               const FArrayBox& a_crseFAB,
                               const Box&       a_crseBox,
                               const int        a_interpDir,
                               const IntVect&   a_refRatio);

    /**
     * \brief      2D, conservative, quadratic interpolation for use along
     *             parent cell faces in 3D.
     *
     * \param[out] a_fineFAB    Fine-level interpolated data
     * \param[in]  a_crseFAB    Coarse-level data to be interpolated.
     * \param[in]  a_crseBox    The region to be interpolated.
     * \param[in]  a_fcDir      The dir normal to the faces.
     * \param[in]  a_refRatio   The refinement ratio
     *
     * \pre
     *  This function assumes all crseFAB ghosts have been filled.
     *
     * \details
     *  - crseBox must have the same centering as the FABs.
     *  - crseBox and the FABs must be face-centered in fcDir.
     *  - Interpolation occurs in all directions other than fcDir.
     *  - All fine data at parent faces (directly above crseBox) will be filled.
     *  - None of the fine data interior to the parent faces will be altered.
     *  - crseFAB must be large enough for the stencil (one ghost layer).
     *
     *  This stencil comes from setting u(x,y) = F_{x,y} (sub means partials),
     *  setting F = axy + bx^2y + cxy^2 + dx^3y + ex^2y^2 + fxy^3
     *  noting that only terms with x and y would survive differentiation,
     *  then matching uc, the 2 first derivs, and the 3 second derivs.
     */
    static void
    interpAlongParentFacesIn3D(FArrayBox&       a_fineFAB,
                               const FArrayBox& a_crseFAB,
                               const Box&       a_crseBox,
                               const int        a_fcDir,
                               const IntVect&   a_refRatio);

#ifdef ALLOW_DIVFREEINTERP
    /**
     * \brief      Divergence-free CF interpolator.
     * \see        divFreeInterp
     * \details    a_crse's grids must be a coarsened version of a_fine's grids.
     *             a_fineParentPtr can be defined on any set of grids.
     */
    virtual void
    localDivFreeInterp(LevelData<FluxBox>&       a_fine,
                       const LevelData<FluxBox>& a_crse,
                       const LevelData<FluxBox>* a_fineDataPtr = nullptr,
                       const Copier* a_fineDataCopierPtr       = nullptr) const;

    /**
     * \brief      Divergence-free CF interpolator.
     * \param[out] a_fineAdvVelFB     The fine level, advecting vel, J*u^i.
     * \param[in]  a_crseAdvVelFB     The crse level, advecting vel, J*u^i.
     * \param[in]  a_ccCrseInterpBox  The region to be interpolated.
     * \see        divFreeInterp
     * \details    a_fineFB's Box must be a refined version of a_crseFB's Box.
     *             This performs step 1 of 2 in the process -- interpolates
     *             (u,v) to the parent faces and in 3D, interpolates w to all
     *             faces.
     */
    virtual void
    localDivFreeInterp_Step1(FluxBox&       a_fineAdvVelFB,
                             const FluxBox& a_crseAdvVelFB,
                             const Box&     a_ccCrseInterpBox) const;

    /**
     * \brief        Divergence-free CF interpolator.
     * \param[inout] a_fineAdvVelFB     The fine-level, advecting vel, J*u^i.
     * \param[in]    a_crseAdvVelFB     The crse-level, advecting vel, J*u^i.
     * \param[in]    a_ccCrseInterpBox  The region to be interpolated.
     * \param[in]    a_crseJgupFB       The crse-level J*g^{ii}.
     * \see          divFreeInterp
     * \details      a_fineFB's Box must be a refined version of a_crseFB's Box.
     *               This performs step 2 of 2 in the process -- interpolation
     *               of (u,v) from the parent faces to the parent interiors.
     */
    virtual void
    localDivFreeInterp_Step2(FluxBox&       a_fineAdvVelFB,
                             const FluxBox& a_crseAdvVelFB,
                             const Box&     a_ccCrseInterpBox/*,
                             const FluxBox& a_crseJgupFB*/) const;

    /**
     * \brief      Computes the curl of velocity in xy-plane.
     *
     * \param[out] a_curlFAB   NC in 2D, EC in 3D.
     * \param[in]  a_curlBox   Where do you want the curl? Same NC/EC in 2D/3D.
     * \param[in]  a_covelFB   The FC covariant velocity, u_i.
     * \param[in]  a_dXi       The grid spacing.
     *
     * \details
     *  This returns J*curl[vel] in the mapped basis, which technically isn't
     *  the curl, but the components of d(one-form version of vel).
     */
    virtual void
    horizNodeCurl(FArrayBox&      a_curlFAB,
                  const Box&      a_curlBox,
                  const FluxBox&  a_covelFB,
                  const RealVect& a_dXi) const;

    /**
     * \brief      Bi-linear interpolation used on horiz curl.
     * \details
     *  Used to interpolate the horizontal curl.
     *  In 3D, this does NOT interpolate in the z dir, it only fills a_fineFAB
     *  at the bottom z-centered edge of each parent cell. It is up to you to
     *  interpolate to the remainiing z-centered edges.
     */
    virtual void
    horizNodeInterp(FArrayBox&       a_fineFAB,
                    const FArrayBox& a_crseFAB,
                    Box              a_crseBox,
                    const IntVect&   a_refRatio) const;

    /**
     * \brief      Vertical, explicit, cubic interpolation for 3D only.
     *             FABs must be FC in vertical dir.
     * \details
     *  Used to interpolate the vertical velocity.
     *  This will call interpAlongParentFacesIn3D first, then interpolate
     *  to faces between.
     */
    virtual void
    vertVelInterp(FArrayBox&       a_fineFAB,
                  const FArrayBox& a_crseFAB,
                  const Box&       a_ccCrseInterpBox) const;
#endif // ALLOW_DIVFREEINTERP

    /// Is there a CFI that requires interpolation?
    inline virtual bool
    hasCFI() const
    {
        CH_assert(m_isDefined);
        return !m_cfiIter.isEmpty();
    }

    bool m_isDefined;

    /// \name Fine-level stuff
    /*\{*/
    mutable CFIIter   m_cfiIter;
    ProblemDomain     m_domain;
    DisjointBoxLayout m_grids;
    IntVect           m_refRatio;
    RealVect          m_fineDXi;
    RealVect          m_crseDXi;
    /*\}*/

    /// \name Coarsened fine-level stuff. These are compatible with m_grids.
    /*\{*/
    DisjointBoxLayout m_crseGrids;
    IntVect           m_crseGhosts;
    Copier            m_crseExCopier;
    Copier            m_crseExCornerCopier;
    /*\}*/


    /// \name User's coarse grid stuff
    /*\{*/

    /// User's coarse-level grids.
    DisjointBoxLayout m_userCrseGrids;

    /// Copiers for m_userCrseGridsPtr <-> m_crseGrids operations.
    Copier m_crseToUserCopier;

    /*\}*/


    /**
     * \name Matrix inversion tools
     *
     * We need to solve
     * \f$[\nabla\cdot \vec{u}^f = 0;
     * \nabla\times\vec{u}^f = \nabla\times\vec{u}^c]\f$
     * for the fine velocity \f$\vec{u}^f\f$. So, we write in matrix form
     * \f$M\vec{u}^f = \vec{b}\f$ and invert \f$\vec{u}^f = M^{-1}\vec{b}\f$.
     * The matrix \f$M^{-1}\f$ only depends on this level's refRatio and grid
     * spacing, so we can set it and forget it in the constructor.
     *
     * Each row of \f$M\f$ corresponds to an equation -- a divergence condition
     * at each CC and a curl condition at each NC. Each column of \f$M\f$
     * corresponds to the unknown vel comps at each FC.
     * To help identify the rows and columns, we will create BaseFab<int>s
     * that map each unknown location to the column index. Similarly, each
     * equation is computed at a location, so we map each of these location
     * to a row index.
     *
     * When the matrix gets inverted, the row and column index BaseFabs will
     * be flipped to correctly describe \f$M^{-1}\f$.
     *
     * These BaseFabs are actually Tuples. Component 0 = u/div and 1 = v/curl.
     *
     * I am using BaseFabs because they are actually pretty convenient
     * for this type of job.
    */
    /*\{*/
    Tuple<BaseFab<int>, 2> m_rowIdx;
    Tuple<BaseFab<int>, 2> m_colIdx;
    BaseFab<Real>          m_Minv;   // The cached matrix inverse.
    BaseFab<Real>          m_DCt;    // Transpose of div/curl matrix.

#ifdef ALLOW_DIVFREEINTERP
    // [working] 0 = original,
    // [2D only] 1 = standard least squares (LS) on div and curl conditions.
    // [not yet] 2 = LS to minimize div, constrained by curl.
    // [not yet] 3 = LS on div and quad interp on vel. (Vanella, et al)
    static constexpr int s_mode = 1;

    // Sets the unknown velocity indices.
    virtual void defineRowIdx(const IntVect& a_refRatio);

    // Mode 0
    virtual void defineColIdx_Original(const IntVect& a_refRatio);
    virtual void invertMatrix_Original(const IntVect& a_refRatio);

    // Mode 1
    virtual void defineColIdx_StdLS(const IntVect& a_refRatio);
    virtual void invertMatrix_StdLS(const IntVect& a_refRatio);
#endif // ALLOW_DIVFREEINTERP
    /*\}*/
};


// -----------------------------------------------------------------------------
// Is this object ready for use?
// -----------------------------------------------------------------------------
bool
CFInterp::isDefined() const
{
    return m_isDefined;
}


// -----------------------------------------------------------------------------
// Basic accessor
// -----------------------------------------------------------------------------
const ProblemDomain&
CFInterp::getDomain() const
{
    return m_domain;
}


// -----------------------------------------------------------------------------
// Basic accessor
// -----------------------------------------------------------------------------
const Box&
CFInterp::getDomainBox() const
{
    return m_domain.domainBox();
}


// -----------------------------------------------------------------------------
// Basic accessor
// -----------------------------------------------------------------------------
const DisjointBoxLayout&
CFInterp::getBoxes() const
{
    return m_grids;
}


// -----------------------------------------------------------------------------
// Basic accessor
// -----------------------------------------------------------------------------
const DisjointBoxLayout&
CFInterp::getCrseBoxes() const
{
    return m_userCrseGrids;
}


// -----------------------------------------------------------------------------
// Basic accessor
// -----------------------------------------------------------------------------
const IntVect&
CFInterp::getRefRatio() const
{
    return m_refRatio;
}


#endif  //!___CFInterp_H__INCLUDED___
