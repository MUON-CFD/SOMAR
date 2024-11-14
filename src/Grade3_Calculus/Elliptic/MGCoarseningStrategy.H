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
#ifndef ___EllipticMGCoarseningStrategy_H__INCLUDED___
#define ___EllipticMGCoarseningStrategy_H__INCLUDED___

#include <algorithm>
#include <vector>
#include "IntVect.H"
#include "RealVect.H"
#include "DisjointBoxLayout.H"


namespace Elliptic {


/**
 * @interface MGCoarseningStrategy
 * @brief     The interface for all coarsening strategies. MGSolver will use
 *            this to choose the amount of refinement at each MG depth (the MG
 *            refinement schedule).
 * @date      2019/11/26
 * @author    Edward.Santilli@jefferson.edu
 */
class MGCoarseningStrategy
{
public:
    virtual ~MGCoarseningStrategy() {}

    /**
     * @brief      Returns the MG refinement schedule.
     *
     * @param[in]  a_grids       The grids to be caorsened. (The finest MG
     *                           level.)
     * @param[in]  a_minBoxSize  The smallest a Box can be. Stop coarsening if
     *                           we find a Box with this size, even if we are
     *                           not at a_maxDepth.
     * @param[in]  a_maxDepth    Do not coarsen more than this. the finest MG
     *                           level is depth = 0 and we count up from there.
     *                           If a_depth = -1, then we will coarsen as much
     *                           as possible (until we reach a Box with
     *                           a_minBoxSize).
     *
     * @return     The MG grids at all depths. The Vector index is the depth and
     *             starts at zero. So, the number of coarsenings will be the
     *             size of the return vector.
     */
    virtual Vector<IntVect>
    createMGRefSchedule(const DisjointBoxLayout& a_grids,
                        const IntVect&           a_minBoxSize = IntVect::Unit,
                        const int                a_maxDepth   = -1) = 0;
};


/**
 * @class     SemicoarseningStrategy
 * @brief     This will create a set of refRatios that push the grids towards
 *            isotropy (dx = dy = dz) as quickly as possible.
 * @date      2019/11/26
 * @author    Edward.Santilli@jefferson.edu
 */
class SemicoarseningStrategy : public MGCoarseningStrategy
{
public:
    /**
     * @brief            Full constructor
     *
     * @param[in]  a_L   The length of the computational domain.
     */
    SemicoarseningStrategy(const RealVect& a_L,
                           const int       a_verbosity = 0);

    ///
    virtual Vector<IntVect>
    createMGRefSchedule(const DisjointBoxLayout& a_grids,
                        const IntVect&           a_minBoxSize = IntVect::Unit,
                        const int                a_maxDepth   = -1);

protected:
    /// Returns dx*dy*dz / min(dx, dy, dz) as a measure of the anisotropy.
    /// 1.0 is isotropic. Larger values represent more anisotropy.
    inline static Real
    measureAnisotropy(const RealVect& a_dx) {
        const Real minDx =
            std::min(std::min(a_dx[0], a_dx[1]), a_dx[SpaceDim - 1]);
        return (a_dx / minDx).product();
    }

    /// @brief Decides the best MG refRatio to promote fast convergence.
    ///        Returns IntVect::Unit if no coarsening is possible.
    /// @param a_dx     The current grid's cell spacing.
    /// @param a_grids  The current grids to be coarsened.
    static IntVect
    computeNextRefRatio(const RealVect&          a_dx,
                        const DisjointBoxLayout& a_grids);

    RealVect m_L;
    int      m_verbosity;
    static const std::vector<IntVect> s_refList;
};


/**
 * @class     HorizCoarseningStrategy
 * @brief     This will create a set of refRatios that push the grids towards
 *            isotropy in the horizontal (dx = dy) as quickly as possible.
 *            Useful when using a line relaxation technique.
 * @date      2024/05/27
 * @author    Edward.Santilli@jefferson.edu
 */
class HorizCoarseningStrategy : public MGCoarseningStrategy
{
public:
    /**
     * @brief            Full constructor
     *
     * @param[in]  a_L   The length of the computational domain.
     */
    HorizCoarseningStrategy(const RealVect& a_L,
                            const bool      a_doVertCoarsening,
                            const int       a_verbosity = 0);

    ///
    virtual Vector<IntVect>
    createMGRefSchedule(const DisjointBoxLayout& a_grids,
                        const IntVect&           a_minBoxSize = IntVect::Unit,
                        const int                a_maxDepth   = -1);

protected:
    /// Returns dx*dy / min(dx, dy) as a measure of the anisotropy.
    /// 1.0 is isotropic. Larger values represent more anisotropy.
    inline static Real
    measureAnisotropy(const RealVect& a_dx) {
        const Real minDx = std::min(a_dx[0], a_dx[SpaceDim - 2]);
        return D_TERM(, (a_dx[0] / minDx), * (a_dx[1] / minDx));
    }

    /// @brief Decides the best MG refRatio to promote fast convergence.
    ///        Returns IntVect::Unit if no coarsening is possible.
    /// @param a_dx     The current grid's cell spacing.
    /// @param a_grids  The current grids to be coarsened.
    static IntVect
    computeNextRefRatio(const RealVect&          a_dx,
                        const DisjointBoxLayout& a_grids);

    RealVect m_L;
    bool     m_doVertCoarsening;
    int      m_verbosity;
    static const std::vector<IntVect> s_refList;
};


/**
 * @class     MiniVCycleStrategy
 * @brief     In AMRMGSolvers, you might need to coarsen from l to l-1 where
 *            the refinement ratio is > 2 in some or all of the directions.
 *            this is bad for MG. For this purpose, you must perform a mini
 *            V-cycle before coarsening to l-1. The mini V-Cycle will coarsen
 *            in such a way that it tries to go to isotropy as fast as possible
 *            and it terminates at a depth that is one coarsening (by 2) away
 *            from l-1. In short, we need to eliminate all of the error modes
 *            not captured by the grid at l-1.
 * @date      2020/05/25
 * @author    Edward.Santilli@jefferson.edu
 */
class MiniVCycleStrategy : public MGCoarseningStrategy
{
public:
    /**
     * @brief            Full constructor
     *
     * @param[in]  a_L   The length of the computational domain.
     */
    MiniVCycleStrategy(const RealVect& a_L,
                       const IntVect&  a_amrRefRatio,
                       const int       a_verbosity = 0);

    /// a_depth will be ignored. We coarsen as much as is needed.
    virtual Vector<IntVect>
    createMGRefSchedule(const DisjointBoxLayout& a_grids,
                        const IntVect&           a_minBoxSize = IntVect::Unit,
                        const int                a_maxDepth   = -1);

protected:
    RealVect m_L;
    IntVect  m_amrRefRatio;
    int      m_verbosity;
};


}; // namespace Elliptic
#endif //!___EllipticMGCoarseningStrategy_H__INCLUDED___
