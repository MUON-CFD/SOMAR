/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2019
 *    Jefferson (Philadelphia University + Thomas Jefferson University) and
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
#ifndef ___FluxRegisterFace_H__INCLUDED___
#define ___FluxRegisterFace_H__INCLUDED___

#include "LevelGeometry.H"
#include "StaggeredFluxLD.H"
#include <array>


/**
 * \class   FluxRegisterFace
 * \brief   A FluxRegister for data defined on a StaggeredFluxLD.
 * \date    2019/12/13
 * \author  Edward.Santilli@jefferson.edu
 */
class FluxRegisterFace
{
    // This object is not copyable.
    FluxRegisterFace(const FluxRegisterFace&) = delete;
    FluxRegisterFace& operator=(const FluxRegisterFace&) = delete;

public:

    /// \name Construction / destruction
    /// \{

    /// Weak constructor. Leaves object undefined.
    FluxRegisterFace ();

    /// Full constructor. Leaves object in a usable state.
    FluxRegisterFace(const LevelGeometry* a_fineLevGeoPtr,
                     const LevelGeometry* a_crseLevGeoPtr,
                     const int            a_numComps           = 1,  // Remove
                     const bool           a_doNormalFluxes     = true,
                     const bool           a_doTransverseFluxes = true);

    /// Full virtual constructor. Leaves object in a usable state.
    virtual void
    define(const LevelGeometry* a_fineLevGeoPtr,
           const LevelGeometry* a_crseLevGeoPtr,
           const int            a_numComps           = 1, // Remove
           const bool           a_doNormalFluxes     = true,
           const bool           a_doTransverseFluxes = true);

    /// Virtual destructor.
    virtual ~FluxRegisterFace();

    /// Free memory. Leaves object undefined.
    virtual void
    clear();

    /// \}

    /// \name Basic accessors and inquiries
    /// \{

    /// Returns if this object is ready for use.
    virtual bool
    isDefined() const {return m_isDefined;}

    /// Returns the fine-level grids on which this object was defined.
    virtual const DisjointBoxLayout&
    getFineBoxes() const {return m_fineGrids;}

    /// Returns the coarse-level grids on which this object was defined.
    virtual const DisjointBoxLayout&
    getCrseBoxes() const {return m_crseGrids;}

    /// \}


    /// \name Main refluxing functions
    /// \{

    /// Initialize the flux register.
    virtual void
    setToZero();

    /// \}

    /**
     * @brief      Increments the register with (scaled) data from a_crseFlux.
     *
     * @param      a_crseFluxFAB   The edge-centered coarse fluxes.
     * @param[in]  a_scale         Typically, this is dt / crse dx in a_cfiDir.
     * @param[in]  a_crseDI
     * @param[in]  a_srcInterval   The a_crseFAB comps to take data from.
     * @param[in]  a_destInterval  The register comps to increment.
     * @param[in]  a_fcDir         The FC dir of the velocity comp being refluxed.
     * @param[in]  a_divDir        The boundary-normal direction / derivative
     *                             direction in the divergence.
     * @param[in]  a_cfiSide       Is this the Hi or Lo side of the _fine_ patch?
     *
     * @details
     *  a_crseFAB must contain the coarse fluxes for m_crseGrids[a_crseDI].
     *  This operation is local.
     *  incrementCoarse does register -= scale * crseFlux
     *  for both Lo and Hi side.
     */
    virtual void
    incrementCoarse(const FArrayBox&      a_crseFluxFAB,
                    const Real            a_scale,
                    const DataIndex&      a_crseDI,
                    const Interval&       a_srcInterval,
                    const Interval&       a_destInterval,
                    const int             a_fcDir,
                    const int             a_divDir,
                    const Side::LoHiSide& a_cfiSide);

    /**
     * @brief      Increments the register with (scaled) data from a_fineFlux.
     *
     * @param      a_fineFluxFAB   The edge-centered fine fluxes.
     * @param[in]  a_scale         Typically, this is dt / fine dx in a_cfiDir.
     * @param[in]  a_fineDI
     * @param[in]  a_srcInterval   The a_fineFAB comps to take data from.
     * @param[in]  a_destInterval  The register comps to increment.
     * @param[in]  a_fcDir         The FC dir of the velocity comp being refluxed.
     * @param[in]  a_divDir        The boundary-normal direction / derivative
     *                             direction in the divergence.
     * @param[in]  a_cfiSide       Is this the Hi or Lo side of the _fine_ patch?
     *
     * @details
     *  a_fineFAB must contain the fine fluxes for m_fineGrids[a_fineDI].
     *  This operation is local.
     *  incrementFine does register +=  sign(side)*incrScale*Ave(fineFlux)
     *  for both Lo and Hi side.
     */
    virtual void
    incrementFine(const FArrayBox&      a_fineFluxFAB,
                  const Real            a_scale,
                  const DataIndex&      a_fineDI,
                  const Interval&       a_srcInterval,
                  const Interval&       a_destInterval,
                  const int             a_fcDir,
                  const int             a_divDir,
                  const Side::LoHiSide& a_cfiSide);

    /// Reflux does:
    ///  deltaF = fine flux - crse flux
    ///  face divergence -= iside * scale * deltaF
    virtual void
    reflux(LevelData<FluxBox>& a_crseDiv,
           const Real          a_scale = 1.0);

protected:

    // Member variables...
    bool m_isDefined;
    bool m_isNoOp;

    bool m_doNormalFluxes;
    bool m_doTransverseFluxes;

    IntVect m_refRatio;

    DisjointBoxLayout m_fineGrids;
    DisjointBoxLayout m_coarsenedFineGrids;
    DisjointBoxLayout m_crseGrids;

    int m_numComps;


    // Find and store the refluxed velocity locations...

    // Format: T[velComp][cfiDir][side][di]
    template <typename T>
    using Spread = std::array<std::array<std::array<T, 2>, SpaceDim>, SpaceDim>;

    // Crse vel locations to be refluxed.
    //
    // Defined on m_crseGrids, but side is the CFI cide. Careful, this is
    // confusing. The high CFI side could be the low side of m_crseGrids[di]!
    //
    // Format: m_crseLocations[cfiDir][fcDir][side][di], which equals
    //         m_crseLocations[advectingDir][advectedDir][side][di]
    Spread<LayoutData<Vector<Box>>> m_crseLocations;

    /// Compute crse vel locations to be refluxed with normal fluxes.
    static void
    normalLocations(LayoutData<Vector<Box>>& a_locations,
                    const DisjointBoxLayout& a_coarsenedFineGrids,
                    const DisjointBoxLayout& a_crseGrids,
                    const size_t             a_bdryDir,
                    const Side::LoHiSide     a_cfiSide);

    /// Compute crse vel locations to be refluxed with transverse fluxes.
    static void
    transverseLocations(LayoutData<Vector<Box>>& a_locations,
                        const DisjointBoxLayout& a_coarsenedFineGrids,
                        const DisjointBoxLayout& a_crseGrids,
                        const size_t             a_bdryDir,
                        const size_t             a_fcDir,
                        const Side::LoHiSide     a_cfiSide);

    // Data storage...
    LevelData<FluxBox>                    m_coarsenedFineFlux;
    LevelData<FluxBox>                    m_crseFlux;
    std::array<StaggeredCopier, SpaceDim> m_reverseCopier;
};


#endif //!___FluxRegisterFace_H__INCLUDED___
