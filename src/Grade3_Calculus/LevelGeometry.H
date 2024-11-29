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
#ifndef __LEVELGEOMETRY_H__INCLUDED__
#define __LEVELGEOMETRY_H__INCLUDED__

#include <array>
#include "LevelData.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "CFRegion.H"
#include "GeoSourceInterface.H"
#include "GNUC_Extensions.H"

// -----------------------------------------------------------------------------
// A container for the metric, connection, Jacobian, etc...
// This class completely describes the geometry of the problem's domain.
// -----------------------------------------------------------------------------
class LevelGeometry
{
public:
    // Construction / destruction ----------------------------------------------

     // Full constructor (calls define)
    LevelGeometry(const ProblemDomain&      a_domain,
                  const RealVect&           a_domainLength,
                  LevelGeometry*            a_crseLevGeoPtr = nullptr,
                  const GeoSourceInterface* a_geoSourcePtr  = nullptr);

    // Puts this Geometry object in a usable state. Called by the full
    // constructor. a_crseLevGeoPtr links this with the coarser level, which
    // sets the refRatios. The finer pointer will still be nullptr. It is
    // expected that these defines will be called from coarsest to finest. If
    // a_geoSourcePtr is nullptr, we will attempt to take a pointer from the
    // coarser LevelGeometry.
    void
    define(const ProblemDomain&      a_domain,
           const RealVect&           a_domainLength,
           LevelGeometry*            a_crseLevGeoPtr = nullptr,
           const GeoSourceInterface* a_geoSourcePtr  = nullptr);

    // Destructor
    ~LevelGeometry ();


    // Multigrid construction --------------------------------------------------

    // Creates a new version of this levGeo, but coarsened by a_crseRefRatio.
    // It is up to you to delete the resulting object when you are through.
    // If a_crseRefRatio == (1,1,1), then the cached metric will be aliased.
    // The result will NOT be linked to coarser/finer levGeos.
    LevelGeometry*
    newCoarsenedLG(const IntVect& a_crseRefRatio) const;


    // Basic accessors ---------------------------------------------------------

    // Returns the ProblemDomain
    inline const ProblemDomain& getDomain () const;

    // Returns the ProblemDomain's Box.
    inline const Box& getDomainBox () const;

    // Returns the domain's physical extents
    inline const RealVect& getDomainLength () const;
    inline Real getDomainLength (int a_dir) const;

    // Returns this level's grid spacing
    inline const RealVect& getDXi () const;
    inline Real getDXi (const int a_dir) const;


    // GeoSourceInterface accessors --------------------------------------------

    // Returns a pointer to the GeoSourceInterface. Originally, this function
    // was supplied for use in MG functions that need to compute metric data
    // with a coarser dXi than what is in the metric cache. But, there is no
    // harm in providing access to everyone.
    inline const GeoSourceInterface& getGeoSource () const;

    // Returns current coordinate map name
    inline const char* getCoorMapName () const;

    // Returns if this metric is constant in space
    inline bool isUniform () const;


    // Coarse/fine level accessors ---------------------------------------------
    // Try to minimize your usage of these functions, they may be phased out.
    //
    // Over the years, these versions caused headaches. For example, when I need
    // a levGeo for a PoissonOp at level l, depth > 0, should it be linked to
    // the other levGeos in the hierarchy? If so, the AMR or MG heirarchy?
    // There are benefits and headaches either way. Another example is when
    // you have AMR levels 0 through 5, but want to do something with only
    // levels 0 through 2. So you sever the link between the levGeos at levels
    // 2 and 3, which may unintentionally mess up other parts of the code.
    //
    // The best thing to do is let each levGeo know only about its own level and
    // ask AMRNSLevel, AMRMGSolver, etc for hierarchical info when necessary.
    // Many of these functions can be found in AMRNSLevel. Use those versions.
    //
    // I'm leaving these here for now. The IO and Integral classes use these
    // functions and I don't want to get into re-coding those classes right now.

    // Return a pointer to the coarser LevelGeometry
    inline const LevelGeometry* getCoarserPtr () const;

    // Set/Return a pointer to the finer LevelGeometry
    inline const LevelGeometry* getFinerPtr () const;

    // Returns the coarsest LevGeo in the hierarchy.
    const LevelGeometry* getCoarsestPtr() const;

    // Returns all levGeos over the AMR hierarchy
    Vector<const LevelGeometry*> getAMRLevGeos () const;

    // Returns all fineRefRatios over the AMR hierarchy
    Vector<IntVect> getAMRRefRatios () const;

    // Returns the LevelGeometry that is at a_domBox's index space or nullptr
    // if none found.
    const LevelGeometry* getLevGeoPtr (const Box& a_domBox) const;

    // Returns all DisjointBoxLayouts over the AMR hierarchy
    Vector<DisjointBoxLayout> getAMRGrids () const;

    // Returns the refinement ratio to the next coarser level
    inline const IntVect& getCrseRefRatio () const;

    // Returns the refinement ratio to the next coarser level
    inline int getCrseRefRatio (const int a_dir) const;

    // Returns the refinement ratio to the next finer level
    inline const IntVect& getFineRefRatio () const;

    // Returns the refinement ratio to the next finer level
    inline int getFineRefRatio (const int a_dir) const;

    // Returns a pointer to the coarser grids. nullptr if this is level 0.
    inline DisjointBoxLayout* getCrseGridsPtr () const;

    // Returns a pointer to the finer grids. nullptr if this is the finest level.
    inline DisjointBoxLayout* getFineGridsPtr () const;

    // Returns the number of levels attached through getCoarserPtr/getFinerPtr.
    // You can call this from any level.
    int getNumLevels() const;


    // Metric cache ------------------------------------------------------------

    // Evaluates commonly used metric data over the specified grids.
    // Before calling this function, this LevelGeometry has no knowledge of this
    // level's grids.
    void createMetricCache (const DisjointBoxLayout& a_newGrids);

    // This copies an already established metric cache to a new set of grids.
    // The new grids can have a new decomposition and set of processor
    // assignments, but can not contain cells in new locations. For that, you
    // will need to call createMetricCache() to recompute the metric quantities.
    void redistributeCache (const DisjointBoxLayout& a_newGrids);

// //     // Regrids the current geometry data with flat (vertically averaged) metrics.
// //     void regridVertAvg (const DisjointBoxLayout& a_newGrids,
// //                         const Box&               a_fullDomainBox);

    // Clears the metric cache and removes all knowledge of this level's grids.
    void clearMetricCache ();

    // Tests if this LevelGeometry has been fit with grids
    inline bool hasGrids () const;

    // Returns the current DisjointBoxLayout
    inline const DisjointBoxLayout& getBoxes () const;

    // Returns the CFRegion
    inline const CFRegion& getCFRegion() const;

    // Metric cache accessors
    inline const LevelData<FArrayBox>& getCCJ () const;
    inline const LevelData<FArrayBox>& getCCJinv () const;
    inline const LevelData<FluxBox>& getFCJgup () const;

    Vector<const LevelData<FArrayBox>*> getAMRCCJ() const;
    Vector<const LevelData<FArrayBox>*> getAMRCCJinv() const;
    Vector<const LevelData<FluxBox>*> getAMRFCJgup() const;

    // Fill functions ----------------------------------------------------------

    // Fills a FAB with displacements from Xi to physical locations.
    // For use with VisIt's displace operator.
    void fill_displacement (FArrayBox& a_dest) const;

    // Fills a mapped box with physical locations
    // a_dest must have SpaceDim comps.
    void fill_physCoor (FArrayBox& a_dest) const;

    // Fills a mapped box with physical locations
    // Single component version.
    void fill_physCoor(FArrayBox& a_dest,
                       const int  a_destComp,
                       const int  a_coorComp) const;

//     // Ensures that metric data is the average of finer data.
//     void averageMetricsDown ();


    // Coordinate accessors ----------------------------------------------------

    /// Retrieves the coordinates of the entire size of a Box.
    /// The Box will determine if the results are cell or face centered.
    std::vector<Real>
    getX(const int a_coorDir, const Box& a_bx) const;

    /// Retrieves the cell coordinates.
    inline Real
    getCellX(const int a_coorDir, const IntVect& a_iv) const;

    /// Retrieves the cell coordinates.
    inline RealVect
    getCellX(const IntVect& a_iv) const;

    /// Retrieves the face coordinates.
    inline Real
    getFaceX(const int a_coorDir, const IntVect& a_iv, const int a_fcDir) const;

    /// Retrieves the face coordinates.
    inline RealVect
    getFaceX(const IntVect& a_iv, const int a_fcDir) const;

    /// Retrieves the lower-left coordinates of a Box.
    inline RealVect
    getBoxLoX(const Box& a_box) const;

    /// Retrieves the upper-right coordinates of a Box.
    inline RealVect
    getBoxHiX(const Box& a_box) const;

    /// Retrieves the lower-left and upper-right coordinates of a Box.
    inline std::array<RealVect, 2>
    getBoxExtents(const Box& a_box) const;

    // Cell/face size accessors ------------------------------------------------

    /// Computes the cell width at a_cc.
    inline Real
    getCellDx(const IntVect& a_cc, const int a_dir) const;

    /// Computes the cell width at a_cc.
    inline RealVect
    getCellDx(const IntVect& a_iv) const;

    /// Similar to getCellDx, but for FC locations.
    /// Returns the average of adjacent getCellDx's.
    inline Real
    getFaceDx(const int      a_coorDir,
              const IntVect& a_iv,
              const int      a_fcDir) const;

    /// Similar to getCellDx, but for FC locations.
    /// Returns the average of adjacent getCellDx's.
    inline RealVect
    getFaceDx(const IntVect& a_iv, const int a_fcDir) const;

    /// Returns the volume of a cell, dx*dy*dz.
    inline Real
    getCellVolume(const IntVect& a_iv) const;

    /// Similar to getCellVolume, but centered at a face.
    /// NOTE: This does not return an area, it returns a volume!
    inline Real
    getFaceVolume(const IntVect& a_iv, const int a_fcDir) const;

    /// Returns the smallest cell size on this level.
    /// This searches over the whole domain, not just the extant grids.
    Real getMinCellDx(const int a_coorDir) const;

    /// Returns the largest cell size on this level.
    /// This searches over the whole domain, not just the extant grids.
    Real getMaxCellDx(const int a_coorDir) const;

    // Search functions --------------------------------------------------------

    /// Finds the indices of the cell that contains a_X.
    /// Returns an estimate if a_X cannot be found within m_xNodeFAB's limits.
    int
    findContainingCellIdxDir(const Real a_X, const int a_dir) const;

    /// Finds the indices of the cell that contains a_X.
    /// Returns an estimate if a_X cannot be found within m_xNodeFAB's limits.
    IntVect
    findContainingCellIdx(const RealVect& a_X) const;


    // Debugging functions -----------------------------------------------------

    // Dumps debugging info to pout()
#ifndef NDEBUG
    void dump () const;
#else
    inline void dump () const {;}
#endif


    // Utilities ---------------------------------------------------------------

    // Sends a vector from the mapped basis to the Cartesian basis.
    // uCart^{a} = [dx^{a} / dXi^{a}] * uMapped^{a}
    void sendToCartesianBasis (FluxBox& a_vectFlub) const;

    // Sends a vector from the Cartesian basis to the mapped basis.
    // uMapped^{a} = [dXi^{a} / dx^{a}] * uCart^{a}
    void sendToMappedBasis (FluxBox& a_vectFlub) const;

    // Multiplies every element by J
    void multByJ (LevelData<FArrayBox>& a_data) const;
    void multByJ (FArrayBox& a_data, const DataIndex& a_di, int a_comp = -1) const;
    void multByJ (LevelData<FluxBox>& a_data) const;
    void multByJ (FluxBox& a_data, const DataIndex& a_di, int a_comp = -1) const;

    // Divides every element by J
    void divByJ (LevelData<FArrayBox>& a_data) const;
    void divByJ (FArrayBox& a_data, const DataIndex& a_di, int a_comp = -1) const;
    void divByJ (LevelData<FluxBox>& a_data) const;
    void divByJ (FluxBox& a_data, const DataIndex& a_di, int a_comp = -1) const;

    // Contracts a CC contra-vector with gdn, making it covariant.
    // Single comp, single grid version.
    void makeCovariant (FArrayBox&       a_coVect,
                        const int        a_coVectComp,
                        const FArrayBox& a_contraVect,
                        const int        a_contraVectComp,
                        const Box&       a_region,
                        const DataIndex& a_di,
                        const int        a_mu) const;

    // Contracts a CC contra-vector with gdn, making it covariant.
    // Single grid version.
    void makeCovariant (FArrayBox&       a_coVect,
                        const FArrayBox& a_contraVect,
                        const Box&       a_region,
                        const DataIndex& a_di) const;

    // Contracts a CC covariant one-form with gup, making it a contravariant vector.
    // Single comp, single grid version.
    void makeContravariant (FArrayBox&       a_contraVect,
                            const int        a_contraVectComp,
                            const FArrayBox& a_coVect,
                            const int        a_coVectComp,
                            const Box&       a_region,
                            const DataIndex& a_di,
                            const int        a_mu) const;

    // Contracts a CC covariant one-form with gup, making it a contravariant vector.
    // Single grid version.
    void makeContravariant (FArrayBox&       a_contraVect,
                            const FArrayBox& a_coVect,
                            const Box&       a_region,
                            const int        a_vectComp,
                            const DataIndex& a_di) const;

    // Compute the magnitude of a vector
    void contractVectors (FArrayBox&       a_mag,
                          const FArrayBox& a_vec1,
                          const FArrayBox& a_vec2,
                          const DataIndex& a_di,
                          const bool       a_isCartesian = false) const;

protected:
    // Basic construction parameters
    ProblemDomain             m_domain;
    RealVect                  m_domainLength;  // The domain's physical extents
    RealVect                  m_dXi;           // This level's grid spacing
    const GeoSourceInterface* m_geoSourcePtr;  // The metric source functions

    // Links to coarse / fine levels
    IntVect        m_crseRefRatio; // Refinement ratio to the next coarser level
    IntVect        m_fineRefRatio; // Refinement ratio to the next finer level
    LevelGeometry* m_coarserPtr;   // Pointer to coarser LevelGeometry
    LevelGeometry* m_finerPtr;     // Pointer to finer LevelGeometry.

    // Metric cache data
    bool                 m_hasGrids;
    DisjointBoxLayout    m_grids;
    CFRegion             m_cfRegion;
    LevelData<FArrayBox> m_CCJCache;
    LevelData<FArrayBox> m_CCJinvCache;
    LevelData<FluxBox>   m_FCJgupCache;

    std::array<FArrayBox, CH_SPACEDIM> m_xCellFAB;
    std::array<FArrayBox, CH_SPACEDIM> m_xNodeFAB;

private:
    // Weak construction not allowed
    LevelGeometry () = delete;

    // Copy constructor not allowed
    LevelGeometry (const LevelGeometry&) = delete;

    // Assignment operator is not allowed
    void operator= (const LevelGeometry&) = delete;
};



// -----------------------------------------------------------------------------
// Implementation of inline functions
// -----------------------------------------------------------------------------

// Returns the current ProblemDomain
const ProblemDomain&
LevelGeometry::getDomain () const
{
    return m_domain;
}


// Returns the ProblemDomain's Box.
const Box&
LevelGeometry::getDomainBox () const
{
    return m_domain.domainBox();
}


// Returns the domain's physical extents
const RealVect&
LevelGeometry::getDomainLength () const
{
    return m_domainLength;
}


// Returns the domain's physical extents in one direction
Real
LevelGeometry::getDomainLength (int a_dir) const
{
    CH_assert(0 <= a_dir && a_dir < SpaceDim);
    return m_domainLength[a_dir];
}


// Returns this level's grid spacing
const RealVect&
LevelGeometry::getDXi () const
{
    return m_dXi;
}


// Returns this level's grid spacing
Real
LevelGeometry::getDXi (const int a_dir) const
{
    CH_assert(0 <= a_dir && a_dir < SpaceDim);
    return m_dXi[a_dir];
}


// Returns a pointer to the GeoSourceInterface. Originally, this function
// was supplied for use in MG functions that need to compute metric data
// with a coarser dXi than what is in the metric cache. But, there is no
// harm in providing access to everyone.
const GeoSourceInterface&
LevelGeometry::getGeoSource () const
{
    return *m_geoSourcePtr;
}


// Returns the coordinate map's name
const char*
LevelGeometry::getCoorMapName () const
{
    return m_geoSourcePtr->getCoorMapName();
}


// Returns whether or not this metric is uniform
bool
LevelGeometry::isUniform () const
{
    return m_geoSourcePtr->isUniform();
}


// Return a pointer to the coarser LevelGeometry
const LevelGeometry*
LevelGeometry::getCoarserPtr () const
{
    return m_coarserPtr;
}


// Return a pointer to the finer LevelGeometry
const LevelGeometry*
LevelGeometry::getFinerPtr () const
{
    return m_finerPtr;
}


// Returns the refinement ratio to the next coarser level
const IntVect&
LevelGeometry::getCrseRefRatio () const
{
    return m_crseRefRatio;
}


// Returns the refinement ratio to the next coarser level
int
LevelGeometry::getCrseRefRatio (const int a_dir) const
{
    CH_assert(0 <= a_dir && a_dir < SpaceDim);
    return m_crseRefRatio[a_dir];
}


// Returns the refinement ratio to the next finer level
const IntVect&
LevelGeometry::getFineRefRatio () const
{
    return m_fineRefRatio;
}


// Returns the refinement ratio to the next finer level
int
LevelGeometry::getFineRefRatio (const int a_dir) const
{
    CH_assert(0 <= a_dir && a_dir < SpaceDim);
    return m_fineRefRatio[a_dir];
}


// Returns a pointer to the coarser grids. nullptr if this is level 0.
DisjointBoxLayout*
LevelGeometry::getCrseGridsPtr () const
{
    return (m_coarserPtr? &(m_coarserPtr->m_grids): nullptr);
}


// Returns a pointer to the finer grids. nullptr if this is the finest level.
DisjointBoxLayout*
LevelGeometry::getFineGridsPtr () const
{
    return (m_finerPtr? &(m_finerPtr->m_grids): nullptr);
}


// Tests if this LevelGeometry has been fit with grids
bool
LevelGeometry::hasGrids () const
{
    return m_hasGrids;
}


// Returns the current DisjointBoxLayout
const DisjointBoxLayout&
LevelGeometry::getBoxes () const
{
    CH_assert(m_hasGrids);
    return m_grids;
}


// Returns the CFRegion
const CFRegion&
LevelGeometry::getCFRegion() const
{
    CH_assert(m_hasGrids);
    return m_cfRegion;
}


// Returns the current CC J field
const LevelData<FArrayBox>&
LevelGeometry::getCCJ () const {
    CH_assert(m_hasGrids);
    return m_CCJCache;
}


// Returns the current CC 1/J field
const LevelData<FArrayBox>&
LevelGeometry::getCCJinv () const {
    CH_assert(m_hasGrids);
    return m_CCJinvCache;
}


// Returns the current FC contravariant metric tensor * J
const LevelData<FluxBox>&
LevelGeometry::getFCJgup () const {
    CH_assert(m_hasGrids);
    return m_FCJgupCache;
}


// -----------------------------------------------------------------------------
Real
LevelGeometry::getCellX(const int a_coorDir, const IntVect& a_iv) const
{
    return m_xCellFAB[a_coorDir](a_iv * BASISV(a_coorDir));
}


// -----------------------------------------------------------------------------
RealVect
LevelGeometry::getCellX(const IntVect& a_iv) const
{
    return RealVect(D_DECL(this->getCellX(0, a_iv),
                           this->getCellX(1, a_iv),
                           this->getCellX(2, a_iv)));
}


// -----------------------------------------------------------------------------
Real
LevelGeometry::getFaceX(const int      a_coorDir,
                        const IntVect& a_iv,
                        const int      a_fcDir) const
{
    if (a_coorDir == a_fcDir) {
        return m_xNodeFAB[a_coorDir](a_iv * BASISV(a_coorDir));
    } else {
        return m_xCellFAB[a_coorDir](a_iv * BASISV(a_coorDir));
    }
}


// -----------------------------------------------------------------------------
RealVect
LevelGeometry::getFaceX(const IntVect& a_iv, const int a_fcDir) const
{
    return RealVect(D_DECL(this->getFaceX(0, a_iv, a_fcDir),
                           this->getFaceX(1, a_iv, a_fcDir),
                           this->getFaceX(2, a_iv, a_fcDir)));
}


// -----------------------------------------------------------------------------
RealVect
LevelGeometry::getBoxLoX(const Box& a_box) const
{
    return RealVect(D_DECL(this->getFaceX(0, a_box.smallEnd(), 0),
                           this->getFaceX(1, a_box.smallEnd(), 1),
                           this->getFaceX(2, a_box.smallEnd(), 2)));
}


// -----------------------------------------------------------------------------
RealVect
LevelGeometry::getBoxHiX(const Box& a_box) const
{
    return RealVect(D_DECL(this->getFaceX(0, a_box.bigEnd() + BASISV(0), 0),
                           this->getFaceX(1, a_box.bigEnd() + BASISV(1), 1),
                           this->getFaceX(2, a_box.bigEnd() + BASISV(2), 2)));
}

// -----------------------------------------------------------------------------
std::array<RealVect, 2>
LevelGeometry::getBoxExtents(const Box& a_box) const
{
    return {this->getBoxLoX(a_box), this->getBoxHiX(a_box)};
}

// -----------------------------------------------------------------------------
Real
LevelGeometry::getCellDx(const IntVect& a_cc, const int a_dir) const
{
    const IntVect& e  = BASISV(a_dir);
    const IntVect  iv = a_cc * e;

    CH_assert(m_xNodeFAB[a_dir].box().contains(iv + e));
    CH_assert(m_xNodeFAB[a_dir].box().contains(iv));

    return m_xNodeFAB[a_dir](iv + e) - m_xNodeFAB[a_dir](iv);
}


// -----------------------------------------------------------------------------
RealVect
LevelGeometry::getCellDx(const IntVect& a_iv) const
{
    return RealVect(D_DECL(this->getCellDx(a_iv, 0),
                           this->getCellDx(a_iv, 1),
                           this->getCellDx(a_iv, 2)));
}


// -----------------------------------------------------------------------------
Real
LevelGeometry::getFaceDx(const int      a_coorDir,
                         const IntVect& a_iv,
                         const int      a_fcDir) const
{
    const IntVect& e        = BASISV(a_coorDir);
    const IntVect  maskedIV = a_iv * e;

    if (a_coorDir == a_fcDir) {
        CH_assert(m_xCellFAB[a_coorDir].box().contains(maskedIV));
        CH_assert(m_xCellFAB[a_coorDir].box().contains(maskedIV - e));

        return m_xCellFAB[a_coorDir](maskedIV) -
               m_xCellFAB[a_coorDir](maskedIV - e);

    } else {
        CH_assert(m_xNodeFAB[a_coorDir].box().contains(maskedIV + e));
        CH_assert(m_xNodeFAB[a_coorDir].box().contains(maskedIV));

        return m_xNodeFAB[a_coorDir](maskedIV + e) -
               m_xNodeFAB[a_coorDir](maskedIV);
    }

}


// -----------------------------------------------------------------------------
RealVect
LevelGeometry::getFaceDx(const IntVect& a_iv, const int a_fcDir) const
{
    return RealVect(D_DECL(this->getFaceDx(0, a_iv, a_fcDir),
                           this->getFaceDx(1, a_iv, a_fcDir),
                           this->getFaceDx(2, a_iv, a_fcDir)));
}


// -----------------------------------------------------------------------------
Real
LevelGeometry::getCellVolume(const IntVect& a_iv) const
{
    return this->getCellDx(a_iv).product();
}


// -----------------------------------------------------------------------------
Real
LevelGeometry::getFaceVolume(const IntVect& a_iv, const int a_fcDir) const
{
    return this->getFaceDx(a_iv, a_fcDir).product();
}

#endif //!__LEVELGEOMETRY_H__INCLUDED__
