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
#ifndef __StretchedMap_H__INCLUDED__
#define __StretchedMap_H__INCLUDED__

#include "GeoSourceInterface.H"


/**
 * @class  StretchedMap
 * @brief  Capable of refining/stretching the grids at the boundaries or in the
 * domain center.
 * @author Edward.Santilli@jefferson.edu
 */
class StretchedMap : public GeoSourceInterface
{
public:
    /// Constructor
    /// Don't get too crazy with the amplitude. If it is too large, the grids
    /// will become kinked and your simulation will crash.
    /// A positve ampl will resolve the middle of the domain.
    /// A negative ampl will resolve the ends of the domain.
    StretchedMap(const RealVect& a_xmin,
                 const RealVect& a_xmax,
                 const RealVect& a_ampl);

    /// Destructor
    virtual ~StretchedMap();

    /// 1. Must return the name of the coordinate mapping
    inline virtual const char*
    getCoorMapName() const
    {
        return "StretchedMap";
    }

    /// 2. Must return whether or not this metric is uniform
    inline virtual bool
    isUniform() const
    {
        return false;
    }

    /// 3. Interplates the function x^mu(xi^mu).
    virtual void
    interp(Vector<Real>& a_x, const Vector<Real>& a_xi, const int a_mu) const;

protected:
    const RealVect m_xmin;
    const RealVect m_xmax;
    const RealVect m_ampl;
    const RealVect m_k;
};


#endif //!__StretchedMap_H__INCLUDED__

