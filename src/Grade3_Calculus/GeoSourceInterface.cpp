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
#include "GeoSourceInterface.H"
#include "GeoSourceInterfaceF_F.H"
#include "FABAlgebra.H"
#include "FiniteDiff.H"
#include "Convert.H"
#include "Debug.H"



// -----------------------------------------------------------------------------
// Default destructor
// -----------------------------------------------------------------------------
GeoSourceInterface::~GeoSourceInterface ()
{;}


// -----------------------------------------------------------------------------
// Computes the physical coordinates of the cells or nodes of a box.
// -----------------------------------------------------------------------------
void
GeoSourceInterface::fill_physCoor(Vector<Real>& a_x,
                                  const int     a_mu,
                                  const Real    a_dXi,
                                  const Box&    a_box) const
{
    const size_t Nx = a_box.size(a_mu);

    // Compute xi at each box location.
    Vector<Real> vxi(Nx);
    {
        // Compute initial xi.
        Real xi = a_box.smallEnd(a_mu) * a_dXi;
        if (a_box.type(a_mu) == 0) {
            xi += 0.5 * a_dXi;
        }

        // Step along, computing xi at each location.
        size_t idx = 0;
        while (idx < Nx) {
            vxi[idx] = xi;
            ++idx;
            xi += a_dXi;
        }
    }

    // Resize a_x, if needed.
    if (a_x.size() != Nx) {
        a_x.resize(Nx);
    }

    // Convert to Cartesian locations.
    this->interp(a_x, vxi, a_mu);
}


// -----------------------------------------------------------------------------
// Fills a mapped box with Cartesian locations.
// -----------------------------------------------------------------------------
void
GeoSourceInterface::fill_physCoor(FArrayBox&      a_dest,
                                  const int       a_destComp,
                                  const int       a_mu,
                                  const RealVect& a_dXi) const
{
    const Box&     destBox     = a_dest.box();
    const IntVect& smallEnd    = destBox.smallEnd();
    const IntVect& bigEnd      = destBox.bigEnd();

    Vector<Real> vx;
    this->fill_physCoor(vx, a_mu, a_dXi[a_mu], destBox);

    // Send to a_dest.
    IntVect iv;
    int idx;
    if (a_mu == 0) {
#if CH_SPACEDIM > 2
        for (iv[2] = smallEnd[2]; iv[2] <= bigEnd[2]; ++iv[2]) {
#endif
            for (iv[1] = smallEnd[1]; iv[1] <= bigEnd[1]; ++iv[1]) {
                for (iv[0] = smallEnd[0], idx = 0; iv[0] <= bigEnd[0]; ++iv[0], ++idx) {
                    a_dest(iv, a_destComp) = vx[idx];
                }
            }
#if CH_SPACEDIM > 2
        }
#endif

    } else if (a_mu == 1) {
#if CH_SPACEDIM > 2
        for (iv[2] = smallEnd[2]; iv[2] <= bigEnd[2]; ++iv[2]) {
#endif
            for (iv[1] = smallEnd[1], idx = 0; iv[1] <= bigEnd[1]; ++iv[1], ++idx) {
                Real y = vx[idx];
                for (iv[0] = smallEnd[0]; iv[0] <= bigEnd[0]; ++iv[0]) {
                    a_dest(iv, a_destComp) = y;
                }
            }
#if CH_SPACEDIM > 2
        }
    } else if (a_mu == 2) {
        for (iv[2] = smallEnd[2], idx = 0; iv[2] <= bigEnd[2]; ++iv[2], ++idx) {
            Real z = vx[idx];
            for (iv[1] = smallEnd[1]; iv[1] <= bigEnd[1]; ++iv[1]) {
                for (iv[0] = smallEnd[0]; iv[0] <= bigEnd[0]; ++iv[0]) {
                    a_dest(iv, a_destComp) = z;
                }
            }
        }
#endif
    } else {
        MAYDAYERROR("a_mu = " << a_mu << ". Cannot be larger than SpaceDim = "
                              << SpaceDim);
    } // if dir == 0, 1, or 2.
}


// -----------------------------------------------------------------------------
// Fills a mapped box with Cartesian locations (a_dest must have SpaceDim comps)
// -----------------------------------------------------------------------------
void GeoSourceInterface::fill_physCoor (FArrayBox&      a_dest,
                                        const RealVect& a_dXi,
                                        const RealVect  a_scale) const
{
    CH_TIME("GeoSourceInterface::fill_physCoor (all comps)");

    for (int dir = 0; dir < CH_SPACEDIM; ++dir) {
        this->fill_physCoor(a_dest,
                            dir,
                            dir,
                            a_dXi);

        if (a_scale[dir] != 1.0) {
            a_dest.mult(a_scale[dir], dir, 1);
        }
    }
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with the Jacobian matrix elements d[x^mu] / d[Xi^mu].
// This is a speed bottleneck!!!
// -----------------------------------------------------------------------------
void GeoSourceInterface::fill_dxdXi (FArrayBox&      a_dest,
                                     const int       a_destComp,
                                     const int       a_mu,
                                     const RealVect& a_dXi,
                                     const Real      a_scale) const
{
    CH_TIME("GeoSourceInterface::fill_dxdXi");

    // Sanity checks
    CH_assert(a_dest.interval().contains(a_destComp));
    CH_assert((0 <= a_mu) && (a_mu < SpaceDim));

    // Stagger the xmu FAB wrt the dest FAB.
    const Box& destBox = a_dest.box();
    const int destBoxType = destBox.type()[a_mu];

    Box xmuBox = destBox;
    if (destBoxType == IndexType::CELL) {
        xmuBox.surroundingNodes(a_mu);
    } else {
        xmuBox.grow(a_mu, 1);
        xmuBox.enclosedCells(a_mu);
    }

    // Fill a FAB with the mapping x^{mu}(Xi)
    FArrayBox xmu(xmuBox, 1);
    this->fill_physCoor(xmu, 0, a_mu, a_dXi);

    // Differentiate the mapping function and apply the scaling
    const Real scaledDXi = a_dXi[a_mu] / a_scale;
    FiniteDiff::partialD(a_dest, a_destComp, destBox, xmu, 0, a_mu, scaledDXi);
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with J = det[Jacobian]
// This is a speed bottleneck!!!
// TODO: try calling dxdXi and using its results. This way, the speed
// bottleneck may be circumvented by rewriting one function instead of two.
// -----------------------------------------------------------------------------
void GeoSourceInterface::fill_J (FArrayBox&      a_dest,
                                 const int       a_destComp,
                                 const RealVect& a_dXi,
                                 const Real      a_scale) const
{
    CH_TIME("GeoSourceInterface::fill_J");

    // Sanity check
    CH_assert(a_dest.interval().contains(a_destComp));

    const Box destBox = a_dest.box();
    FArrayBox dxdXiFAB(destBox, 1);

    a_dest.setVal(a_scale, a_destComp);
    for (int mu = 0; mu < SpaceDim; ++mu) {
        this->fill_dxdXi(dxdXiFAB, 0, mu, a_dXi);
        a_dest.mult(dxdXiFAB, 0, a_destComp, 1);
    }
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with the inverse Jacobian matrix elements d[xi^mu] / d[x^mu].
// -----------------------------------------------------------------------------
void GeoSourceInterface::fill_dXidx (FArrayBox&       a_dest,
                                     const int        a_destComp,
                                     const int        a_mu,
                                     const RealVect&  a_dXi,
                                     const Real       a_scale) const
{
    CH_TIME("GeoSourceInterface::fill_dXidx");

    // Sanity checks
    CH_assert(a_dest.interval().contains(a_destComp));
    CH_assert((0 <= a_mu) && (a_mu < SpaceDim));

    // Compute a_scale / dx/dXi
    this->fill_dxdXi(a_dest, a_destComp, a_mu, a_dXi);
    a_dest.invert(a_scale, a_destComp);
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with 1/J
// -----------------------------------------------------------------------------
void GeoSourceInterface::fill_Jinv (FArrayBox&       a_dest,
                                    const int        a_destComp,
                                    const RealVect&  a_dXi,
                                    const Real       a_scale) const
{
    CH_TIME("GeoSourceInterface::fill_Jinv");

    // Sanity check
    CH_assert(a_dest.interval().contains(a_destComp));

    // Calculate or copy J
    this->fill_J(a_dest, a_destComp, a_dXi);

    // Invert and apply scaling
    a_dest.invert(a_scale, a_destComp);
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with the covariant metric elements
// gdn_{mu,mu} = dx^{mu}/dXi^{mu} * dx^{mu}/dXi^{mu}
// -----------------------------------------------------------------------------
void GeoSourceInterface::fill_gdn (FArrayBox&      a_dest,
                                   const int       a_destComp,
                                   const int       a_mu,
                                   const RealVect& a_dXi,
                                   const Real      a_scale) const
{
    CH_TIME("GeoSourceInterface::fill_gdn");

    // Sanity checks
    CH_assert(a_dest.interval().contains(a_destComp));
    CH_assert((0 <= a_mu) && (a_mu < SpaceDim));

    // Calculate dx^{mu}/dXi^{mu}
    this->fill_dxdXi(a_dest, a_destComp, a_mu, a_dXi);

    // Square it.
    a_dest.mult(a_dest, a_destComp, a_destComp, 1);

    // Apply the scaling
    if (a_scale != 1.0) {
        a_dest.mult(a_scale, a_destComp, 1);
    }
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with the contravariant metric elements
// gup^{mu,mu} = dXi^{mu}/dx^{mu} * dXi^{mu}/dx^{mu}
// -----------------------------------------------------------------------------
void GeoSourceInterface::fill_gup (FArrayBox&       a_dest,
                                   const int        a_destComp,
                                   const int        a_mu,
                                   const RealVect&  a_dXi,
                                   const Real       a_scale) const
{
    CH_TIME("GeoSourceInterface::fill_gup");

    // Sanity checks
    CH_assert(a_dest.interval().contains(a_destComp));
    CH_assert((0 <= a_mu) && (a_mu < SpaceDim));

    // Calculate gdn
    this->fill_gdn(a_dest, a_destComp, a_mu, a_dXi);

    // Invert and apply scaling
    a_dest.invert(a_scale, a_destComp);
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with detJ * gup
// TODO: Make this faster by requesting detJ.
// -----------------------------------------------------------------------------
void GeoSourceInterface::fill_Jgup (FArrayBox&       a_dest,
                                    const int        a_destComp,
                                    const int        a_mu,
                                    const RealVect&  a_dXi,
                                    const Real       a_scale) const
{
    CH_TIME("GeoSourceInterface::fill_Jgup");

    // Sanity checks
    CH_assert(a_dest.interval().contains(a_destComp));
    CH_assert((0 <= a_mu) && (a_mu < SpaceDim));

    FArrayBox dxdXi(a_dest.box(), 1);

    a_dest.setVal(a_scale, a_destComp);
    for (int i = 0; i < SpaceDim; ++i) {
        this->fill_dxdXi(dxdXi, 0, i, a_dXi);
        if (i != a_mu) {
            a_dest.mult(dxdXi, 0, a_destComp, 1);
        } else {
            a_dest.divide(dxdXi, 0, a_destComp, 1);
        }
    }
}

