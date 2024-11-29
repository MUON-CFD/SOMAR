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
#ifndef ___StaggeredCopier_H__INCLUDED___
#define ___StaggeredCopier_H__INCLUDED___

#include "Copier.H"


class StaggeredCopier : public Copier
{
public:
    /// Constructor
    StaggeredCopier();

    /// Throws an error.
    virtual void
    define(const DisjointBoxLayout& a_level,
           const BoxLayout&         a_dest,
           bool                     a_exchange = false,
           IntVect                  a_shift    = IntVect::Zero) override final;

    /// Throws an error.
    virtual void
    define(const DisjointBoxLayout& a_level,
           const BoxLayout&         a_dest,
           const ProblemDomain&     a_domain,
           bool                     a_exchange = false,
           IntVect                  a_shift    = IntVect::Zero) override final;

    /// Throws an error.
    virtual void
    define(const DisjointBoxLayout& a_level,
           const BoxLayout&         a_dest,
           const IntVect&           a_destGhost,
           bool                     a_exchange = false,
           IntVect                  a_shift    = IntVect::Zero) override final;

    /// Throws an error.
    virtual void
    define(const BoxLayout&     a_level,
           const BoxLayout&     a_dest,
           const ProblemDomain& a_domain,
           const IntVect&       a_destGhost,
           bool                 a_exchange = false,
           IntVect              a_shift    = IntVect::Zero) override final;

    /// Throws an error.
    virtual void
    exchangeDefine(const DisjointBoxLayout& a_grids,
                   const IntVect&           a_ghost) override final;

    /// Standard copier for FC data. Only copies valid data.
    virtual void
    define(const BoxLayout&     a_srcLayout,
           const BoxLayout&     a_destLayout,
           const ProblemDomain& a_domain,
           const IntVect&       a_destGhost,
           const int            a_fcDir,
           const IntVect        a_shift = IntVect::Zero);

    /// Throws an error.
    virtual void
    ghostDefine(const DisjointBoxLayout& a_srcGrids,
                const DisjointBoxLayout& a_destGrids,
                const ProblemDomain&     a_domain,
                const IntVect&           a_srcGhost) override final;

    /// Just in case you want to revert.
    virtual void
    standardGhostDefine(const DisjointBoxLayout& a_srcGrids,
                        const DisjointBoxLayout& a_destGrids,
                        const ProblemDomain&     a_domain,
                        const IntVect&           a_srcGhost)
    {
        Copier::define(a_destGrids, a_srcGrids,  a_domain, a_srcGhost);
        this->reverse();
        this->sort();
    }

    /// .
    virtual void
    staggeredGhostDefine(const DisjointBoxLayout& a_srcGrids,
                         const DisjointBoxLayout& a_destGrids,
                         const ProblemDomain&     a_domain,
                         const IntVect&           a_srcGhost,
                         const int                a_fcDir)
    {
        this->define(a_destGrids, a_srcGrids,  a_domain, a_srcGhost, a_fcDir);
        this->reverse();
        this->sort();
    }

    /// \brief Performs the computation necessary for moving from ghost+valid
    ///        cells to valid.
    ///
    /// \param[in] a_srcGrids
    /// \param[in] a_destGrids
    /// \param[in] a_domain
    /// \param[in] a_srcGhost
    /// \param[in] a_fcDir
    ///
    /// \details
    /// Instead of copying from the valid of a_src to the valid+ghost of
    /// a_dest, this routine defines a Copier that moves data from the
    /// ghost+valid region of a_src to the valid region of a_dest.  This is used
    /// by FluxRegisterFace.
    void
    ghostDefine(const DisjointBoxLayout& a_srcGrids,
                const DisjointBoxLayout& a_destGrids,
                const ProblemDomain&     a_domain,
                const IntVect&           a_srcGhost,
                const int                a_fcDir);

    /// \brief Standard exchange Copier for a face-centered (FC) LevelData.
    ///
    /// \param[in] a_grids
    /// \param[in] a_ghostVect
    /// \param[in] a_fcDir
    /// \param[in] a_doValidCorners Exchange the corner faces after the adjacent
    ///                             ghosts are validated.
    ///
    /// \details
    ///  Suppose a 2D FAB is FC in the x-direction. We can identify the faces:
    ///  ┌─────┬─────┬─────┬─────┬─────┬─────┐
    ///  c     t     t     t     t     t     c
    ///  ├─────┼─────┼─────┼─────┼─────┼─────┤     V = valid face
    ///  n     V     V     V     V     V     n     n = normal ghost
    ///  ├─────┼─────┼─────┼─────┼─────┼─────┤     t = transverse ghost
    ///  n     V     V     V     V     V     n     c = corner ghost
    ///  ├─────┼─────┼─────┼─────┼─────┼─────┤
    ///  n     V     V     V     V     V     n
    ///  ├─────┼─────┼─────┼─────┼─────┼─────┤
    ///  n     V     V     V     V     V     n
    ///  ├─────┼─────┼─────┼─────┼─────┼─────┤
    ///  c     t     t     t     t     t     c
    ///  └─────┴─────┴─────┴─────┴─────┴─────┘
    ///  First, this function replaces n and t ghosts with valid data from
    ///  neighboring FABs. This "validates" the n and t faces. Afterwards, if
    ///  a_doValidCorners == true, this function will also replace c ghosts with
    ///  neighboring, valid data.
    ///
    ///  In 2D, the c ghosts are vertices. In 3D, they are edges and vertices.
    ///
    /// \note If you use a standard cell-centered Copier on a LevelData<FluxBox>
    /// overlapping valid faces get swapped or copied. This is not a problem iff
    /// these overlapping valid faces contain identical data, but this is not
    /// always the case. For example, if we average a fine level down to a
    /// coarse level and a fine-level exchange boundary coincides with a
    /// coarse-level exchange boundary. In this situation, overlapping
    /// coarse-level faces contain different data and calling exchange before
    /// fixing the problem would be bad. Which faces would end up with valid
    /// data? The answer depends on the order the faces were copied, which is
    /// random. So, to avoid this and other pitfalls, use this, face-centered
    /// copier, which leaves the overlapping valid faces alone!
    virtual void
    defineValidExchange(const DisjointBoxLayout& a_grids,
                        const IntVect&           a_ghostVect,
                        const int                a_fcDir,
                        const bool               a_doValidCorners = true);

    /// \brief An exchange CornerCopier for face-centered data.
    ///
    /// \param[in] a_grids
    /// \param[in] a_ghostVect
    /// \param[in] a_fcDir
    ///
    /// \details
    ///  Assuming you already used a copier defined via defineValidExchange,
    ///  filled the CFI ghosts, and filled the physical boundary ghosts,
    ///  the ghosts that overlie valid regions have been filled (aka validated).
    ///  We will call these the "valid" ghosts. Next, you may want to take these
    ///  "valid" ghosts and send them to overlying corner ghosts (edges and
    ///  vertices). That's what this function sets up.
    ///
    /// \note
    ///  "Valid" ghosts = normal ghosts
    ///  corner ghosts  = edge ghosts in 3D, vertex ghosts in 2D.
    ///
    /// \warning You must be sure that all other ghosts are filled (validated)
    ///  before performing this type of exchange. That includes CFI ghosts,
    ///  physical boundary ghosts, and exchange ghosts filled by a standard
    ///  defineValidExchange copier. Otherwise, the "valid" data will not be
    ///  "valid!"
    virtual void
    defineInvalidCornerExchange1(const DisjointBoxLayout& a_grids,
                                 const IntVect&           a_ghostVect,
                                 const int                a_fcDir);

    /// \brief The second phase of the defineInvalidCornerExchange process.
    ///
    /// \param[in] a_grids
    /// \param[in] a_ghostVect
    /// \param[in] a_fcDir
    ///
    /// \note
    ///  "Valid" ghosts = edge ghosts
    ///  corner ghosts  = vertex ghosts in 3D
    ///  This is a no-op in 2D.
    ///
    /// \warning You must exchange with a Copier defind via
    ///  defineInvalidCornerExchange1 first!
    virtual void
    defineInvalidCornerExchange2(const DisjointBoxLayout& a_grids,
                                 const IntVect&           a_ghostVect,
                                 const int                a_fcDir);
    ///
    inline virtual IntVect
    getGhostVect() const
    {
        return m_ghostVect;
    }

    ///
    inline virtual int
    getFCDir() const
    {
        return m_fcDir;
    }

protected:
    IntVect m_ghostVect;
    int     m_fcDir;
};


#endif //!___StaggeredCopier_H__INCLUDED___
