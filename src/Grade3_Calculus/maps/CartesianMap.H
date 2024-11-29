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
#ifndef __CARTESIANMAP_H__INCLUDED__
#define __CARTESIANMAP_H__INCLUDED__

#include "GeoSourceInterface.H"


/**
 * @class  CartesianMap
 * @brief  The Cartestian coordinate system's metric source functions.
 * @author Edward.Santilli@jefferson.edu
 */
class CartesianMap: public GeoSourceInterface
{
public:
    /// Destructor
    virtual ~CartesianMap ();


    /// \name Mandatory overrides...
    /// \{

    /// 1. Must return the name of the coordinate mapping
    virtual const char*
    getCoorMapName() const override;

    /// 2. Must return whether or not this metric is uniform
    virtual bool
    isUniform() const override;

    /// 3. Interplates the function x^mu(xi^mu).
    virtual void
    interp(Vector<Real>&       a_x,
           const Vector<Real>& a_xi,
           const int           a_mu) const override;
    /// \}


    /// \name Optional overrides...
    /// \{

    /// Computes the physical coordinates of the cells or nodes of a box.
    virtual void
    fill_physCoor(Vector<Real>& a_x,
                  const int     a_mu,
                  const Real    a_dXi,
                  const Box&    a_box) const override;

    /// Fills a mapped box with Cartesian locations.
    virtual void
    fill_physCoor(FArrayBox&      a_dest,
                  const int       a_destComp,
                  const int       a_mu,
                  const RealVect& a_dXi) const override;

    /// Fills a mapped box with Cartesian locations.
    /// a_dest must have SpaceDim comps.
    virtual void
    fill_physCoor(FArrayBox&      a_dest,
                  const RealVect& a_dXi,
                  const RealVect  a_scale = RealVect::Unit) const override;

    /// Fills an FArrayBox with the Jacobian matrix elements d[x^mu] / d[xi^nu].
    virtual void
    fill_dxdXi(FArrayBox&      a_dest,
               const int       a_destComp,
               const int       a_mu,
               const RealVect& a_dXi,
               const Real      a_scale = 1.0) const override;

    /// Fills an FArrayBox with J = det[Jacobian]
    virtual void
    fill_J(FArrayBox&      a_dest,
           const int       a_destComp,
           const RealVect& a_dXi,
           const Real      a_scale = 1.0) const override;

    /// Fills an FArrayBox with the inverse Jacobian matrix elements
    /// d[xi^mu] / d[x^nu].
    virtual void
    fill_dXidx(FArrayBox&      a_dest,
               const int       a_destComp,
               const int       a_mu,
               const RealVect& a_dXi,
               const Real      a_scale = 1.0) const override;

    /// Fills an FArrayBox with 1/J
    virtual void
    fill_Jinv(FArrayBox&      a_dest,
              const int       a_destComp,
              const RealVect& a_dXi,
              const Real      a_scale = 1.0) const override;

    /// Fills an FArrayBox with the covariant metric elements
    virtual void
    fill_gdn(FArrayBox&      a_dest,
             const int       a_destComp,
             const int       a_mu,
             const RealVect& a_dXi,
             const Real      a_scale = 1.0) const override;

    /// Fills an FArrayBox with the contravariant metric elements
    virtual void
    fill_gup(FArrayBox&      a_dest,
             const int       a_destComp,
             const int       a_mu,
             const RealVect& a_dXi,
             const Real      a_scale = 1.0) const override;

    /// Fills an FArrayBox with J * gup
    virtual void
    fill_Jgup(FArrayBox&      a_dest,
              const int       a_destComp,
              const int       a_mu,
              const RealVect& a_dXi,
              const Real      a_scale = 1.0) const override;
    /// \}

};


#endif //!__CARTESIANMAP_H__INCLUDED__

