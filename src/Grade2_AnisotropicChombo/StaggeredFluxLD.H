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
#ifndef ___StaggeredFluxLD_H__INCLUDED___
#define ___StaggeredFluxLD_H__INCLUDED___

#include <array>
#include "FArrayBox.H"
#include "LevelData.H"
#include "CornerCopier.H"


/**
 * \class   StaggeredFluxLD
 * \brief   A LevelData-like object to hold fluxes of staggered data.
 * \date    2019/12/13
 * \author  Edward.Santilli@jefferson.edu
 * \details
 *  Suppose \f$\vec{u}\f$ is staggered (that is, a LevelData<FluxBox> with one
 *  component). When we update
 *  \f$\frac{du^{\alpha}}{dt} = \nabla_i(u^i u^{\alpha})\f$,
 *  we need storage for the fluxes \f$u^i u^{\alpha}\f$. These fluxes are
 *  cell-centered for \f$ i=\alpha \f$ and edge-centered for \f$ i\neq \alpha\f$
 *  (node centered in both \f$ i \f$ and \f$\alpha\f$).
 *
 *  This class is essentially a 2D array of LevelData's, each defined with the
 *  correct centering, to hold such data and properly perform the exchanges.
 */
class StaggeredFluxLD:
    public std::array<std::array<LevelData<FArrayBox>, CH_SPACEDIM>, CH_SPACEDIM>
{
public:
    /// Default constructor. Leaves object undefined.
    StaggeredFluxLD();

    /// Full constructor. Leaves object in a usable state.
    /// If a_noDiags == true, then (*this)[i][i] will be left undefined.
    StaggeredFluxLD(const DisjointBoxLayout& a_grids,
                    const bool               a_noDiags = false);

    /// Full virtual constructor. Leaves object in a usable state.
    /// If a_noDiags == true, then (*this)[i][i] will be left undefined.
    virtual void
    define(const DisjointBoxLayout& a_grids,
           const bool               a_noDiags = false);

    /// Free memory. Leaves object undefined.
    virtual void
    clear();

    /// Returns if this object is ready for use.
    virtual bool
    isDefined() const {return m_isDefined;}

    /// Returns the grids on which this object was defined.
    virtual const DisjointBoxLayout&
    getBoxes () const {return m_grids;}

    /// MPI data exchange between neighboring boxes. This is blocking.
    virtual void
    exchange();

    /// Sets all data to a single value.
    virtual void
    setVal(const Real a_val);

protected:
    bool                                  m_isDefined;
    bool                                  m_noDiags;
    DisjointBoxLayout                     m_grids;
    std::array<Copier, CH_SPACEDIM>       m_exCopier;
    std::array<CornerCopier, CH_SPACEDIM> m_exCornerCopier;
};


#endif //!___StaggeredFluxLD_H__INCLUDED___
