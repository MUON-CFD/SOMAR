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
#ifndef ___Analysis_H__INCLUDED___
#define ___Analysis_H__INCLUDED___

#include "LevelData.H"
#include "RealVect.H"
#include "FluxBox.H"

// Tools to perform convergence analysis - max, min, and p-norms.
class Analysis
{
public:
    // -------------------------------------------------------------------------
    // max
    // -------------------------------------------------------------------------

    // Computes the maximum over each cell/node in a_region.
    // Warnings:
    //  1. Does not take the abs value!
    //  2. Behavior is undefined if a_data does not contain a_region.
    static Real
    max (const FArrayBox& a_data,
         const Box&       a_region,
         const Interval   a_interv = Interval(0,0));

    // Computes the maximum over each cell/node in a_layout. If you want to
    // include the ghosts of a_data, simply include those ghosts in the layout.
    // This is the most general single-level version. It can ignore cells/nodes
    // that are covered by a finer BoxLayout.
    // Warnings:
    //  1. Does not take the abs value!
    //  2. Behavior is undefined if a_data does not contain a_layout.
    //  3. a_data must be compatible with a_layout.
    static Real
    max (const LevelData<FArrayBox>& a_data,
         const BoxLayout&            a_layout,
         const Interval&             a_interv,
         const BoxLayout*            a_fineLayoutPtr,
         const IntVect&              a_refRatio);

    // Computes the maximum over each cell/node in a_layout. If you want to
    // include the ghosts of a_data, simply include those ghosts in the layout.
    // Warnings:
    //  1. Does not take the abs value!
    //  2. Behavior is undefined if a_data does not contain a_layout.
    //  3. a_data must be compatible with a_layout.
    static Real
    max (const LevelData<FArrayBox>& a_data,
         const BoxLayout&            a_layout,
         const Interval              a_interv = Interval(0,0));

    // Computes the maximum over each cell/node in the valid cells of a_data.
    // Warning: Does not take the abs value!
    static Real
    max (const LevelData<FArrayBox>& a_data,
         const Interval              a_interv = Interval(0,0));


    // Computes the maximum over each valid cell/node in a_data.
    // Warning: Does not take the abs value!
    static Real
    max (const Vector<LevelData<FArrayBox>*>& a_data,
         const int                            a_lmin,
         const int                            a_lmax,
         const Interval&                      a_interv,
         const Vector<IntVect>&               a_refRatio);

    // -------------------------------------------------------------------------
    // min
    // -------------------------------------------------------------------------

    // Computes the minimum over each cell/node in a_region.
    // Warnings:
    //  1. Does not take the abs value!
    //  2. Behavior is undefined if a_data does not contain a_region.
    static Real
    min (const FArrayBox& a_data,
         const Box&       a_region,
         const Interval   a_interv = Interval(0,0));

    // Computes the minimum over each cell/node in a_layout. If you want to
    // include the ghosts of a_data, simply include those ghosts in the layout.
    // This is the most general single-level version. It can ignore cells/nodes
    // that are covered by a finer BoxLayout.
    // Warnings:
    //  1. Does not take the abs value!
    //  2. Behavior is undefined if a_data does not contain a_layout.
    //  3. a_data must be compatible with a_layout.
    static Real
    min (const LevelData<FArrayBox>& a_data,
         const BoxLayout&            a_layout,
         const Interval&             a_interv,
         const BoxLayout*            a_fineLayoutPtr,
         const IntVect&              a_refRatio);

    // Computes the minimum over each cell/node in a_layout. If you want to
    // include the ghosts of a_data, simply include those ghosts in the layout.
    // Warnings:
    //  1. Does not take the abs value!
    //  2. Behavior is undefined if a_data does not contain a_layout.
    //  3. a_data must be compatible with a_layout.
    static Real
    min (const LevelData<FArrayBox>& a_data,
         const BoxLayout&            a_layout,
         const Interval              a_interv = Interval(0,0));

    // Computes the minimum over each cell/node in the valid cells of a_data.
    // Warning: Does not take the abs value!
    static Real
    min (const LevelData<FArrayBox>& a_data,
         const Interval              a_interv = Interval(0,0));


    // Computes the minimum over each valid cell/node in a_data.
    // Warning: Does not take the abs value!
    static Real
    min (const Vector<LevelData<FArrayBox>*>& a_data,
         const int                            a_lmin,
         const int                            a_lmax,
         const Interval&                      a_interv,
         const Vector<IntVect>&               a_refRatio);


    // -------------------------------------------------------------------------
    // p-norm
    // -------------------------------------------------------------------------
protected:
    // Computes a portion of the the p-norm over a_ccRegion using the formula
    //   a_normAccum += Int[|a_data|^p dV]
    //   a_volAccum  += Int[dV]
    // where Int = the integral over a_ccRegion and dV = J * dXi*dEta*dZeta.
    // After some simplification, this function actually uses
    //   a_normAccum += Sum[|a_data|^p * J]
    //   a_volAccum  += Sum[J]
    // where Sum = the sum over a_ccRegion.
    //
    // For the infinity-norm, use a_p = 0;
    //
    // So that the integral (sum) doesn't improperly scale the boundary cells/nodes,
    // this function will convert a_data to cell-centered before summing.
    //
    // This is not a public function and was not created to be particularly
    // user-friendly. It is instead a very general, reliable function used
    // by the public, specialized versions.
    //
    // Warnings:
    //  1. a_ccJ must be cell-centered, as the name suggests.
    //  2. Behavior is undefined if a_data or a_ccJ do not contain a_ccRegion.
    static void
    pNorm_local_Accumulator (Real&            a_normAccum,
                             Real&            a_volAccum,
                             const FArrayBox& a_data,
                             const Interval&  a_dataInterv,
                             const FArrayBox& a_ccJ,
                             const Box&       a_ccRegion,
                             const int        a_p);

    // Computes a portion of the the p-norm over a_layout using the formula
    //   a_normAccum += Int[|a_data|^p dV]
    //   a_volAccum  += Int[dV]
    // where Int = the integral over a_layout and dV = J * dXi*dEta*dZeta.
    // After some simplification, this function actually uses
    //   a_normAccum += Sum[|a_data|^p * J]
    //   a_volAccum  += Sum[J]
    // where Sum = the sum over a_layout.
    //
    // For the infinity-norm, use a_p = 0;
    //
    // So that the integral (sum) doesn't improperly scale the boundary cells/nodes,
    // this function will convert a_data to cell-centered before summing.
    //
    // This is not a public function and was not created to be particularly
    // user-friendly. It is instead a very general, reliable function used
    // by the public, specialized versions.
    //
    // Warnings:
    //  1. a_ccJ must be cell-centered, as the name suggests.
    //  2. Behavior is undefined if a_data or a_ccJ do not contain a_layout.
    static void
    pNorm_local_Accumulator (Real&                       a_normAccum,
                             Real&                       a_volAccum,
                             const LevelData<FArrayBox>& a_data,
                             const Interval&             a_dataInterv,
                             const LevelData<FArrayBox>& a_ccJ,
                             const BoxLayout&            a_layout,
                             const BoxLayout*            a_fineLayoutPtr,
                             const IntVect&              a_refRatio,
                             const int                   a_p);

public:
    // Computes the p-norm over a_ccRegion using the formula
    //   ||a_data||_p = { Int[|a_data|^p dV] / Int[dV] } ^ (1/p)
    // where Int = the integral over a_ccRegion and dV = J * dXi*dEta*dZeta.
    // After some simplification, this function actually uses
    //   ||a_data||_p = { Sum[|a_data|^p * J] / Sum[J] } ^ (1/p)
    // where Sum = the sum over a_ccRegion.
    //
    // So that the integral (sum) doesn't improperly scale the boundary cells/nodes,
    // this function will convert a_data to cell-centered before summing.
    //
    // If a_p is negative, we will take the |a_p| norm without dividing by J or
    // taking the root. This is useful for computing norms over several patches
    // and is mainly provided for internal use.
    //
    // If you supply a_volPtr, then Sum[J] will be added to it.
    //
    // Warnings:
    //  1. a_ccJ must be cell-centered, as the name suggests.
    //  2. Behavior is undefined if a_data and a_ccJ does not contain a_ccRegion.
    static Real
    pNorm (const FArrayBox& a_data,
           const Interval&  a_dataInterv,
           const FArrayBox& a_ccJ,
           const Box&       a_ccRegion,
           const int        a_p);

    // Computes the p-norm over a_layout using the formula
    //   ||a_data||_p = { Int[|a_data|^p dV] / Int[dV] } ^ (1/p)
    // where Int = the integral over a_layout and dV = J * dXi*dEta*dZeta.
    // After some simplification, this function actually uses
    //   ||a_data||_p = { Sum[|a_data|^p * J] / Sum[J] } ^ (1/p)
    // where Sum = the sum over a_layout.
    //
    // This is the most general single-level version. It can ignore cells/nodes
    // that are covered by a finer BoxLayout.
    //
    // If you want to include the ghosts of a_data, simply include those ghosts
    // in the layout.
    //
    // So that the integral (sum) doesn't improperly scale the boundary cells/nodes,
    // this function will convert a_data to cell-centered before summing.
    //
    // Warnings:
    //  1. a_ccJ must be cell-centered, as the name suggests.
    //  2. Behavior is undefined if a_data and a_ccJ does not contain a_layout.
    static Real
    pNorm (const LevelData<FArrayBox>& a_data,
           const Interval&             a_dataInterv,
           const LevelData<FArrayBox>& a_ccJ,
           const BoxLayout&            a_layout,
           const BoxLayout*            a_fineLayoutPtr,
           const IntVect&              a_refRatio,
           const int                   a_p);

    // Computes the p-norm over a_data's grids using the formula
    //   ||a_data||_p = { Int[|a_data|^p dV] / Int[dV] } ^ (1/p)
    // where Int = the integral over a_data's grids and dV = J * dXi*dEta*dZeta.
    // After some simplification, this function actually uses
    //   ||a_data||_p = { Sum[|a_data|^p * J] / Sum[J] } ^ (1/p)
    // where Sum = the sum over a_data's grids.
    //
    // So that the integral (sum) doesn't improperly scale the boundary cells/nodes,
    // this function will convert a_data to cell-centered before summing.
    //
    // Warnings:
    //  1. a_ccJ must be cell-centered, as the name suggests.
    //  2. Behavior is undefined if a_data and a_ccJ are defined over different grids.
    static Real
    pNorm (const LevelData<FArrayBox>& a_data,
           const Interval&             a_dataInterv,
           const LevelData<FArrayBox>& a_ccJ,
           const int                   a_p);

    // Computes the p-norm over a_data's grids using the formula
    //   ||a_data||_p = { Int[|a_data|^p dV] / Int[dV] } ^ (1/p)
    // where Int = the integral over a_data's grids and dV = J * dXi*dEta*dZeta.
    // After some simplification, this function actually uses
    //   ||a_data||_p = { Sum[|a_data|^p * J] / Sum[J] } ^ (1/p)
    // where Sum = the sum over a_data's grids.
    //
    // So that the integral (sum) doesn't improperly scale the boundary cells/nodes,
    // this function will convert a_data to cell-centered before summing.
    //
    // Warnings:
    //  1. a_ccJ must be cell-centered, as the name suggests.
    //  2. Behavior is undefined if a_data and a_ccJ are defined over different grids.
    static Real
    pNorm (const LevelData<FArrayBox>& a_data,
           const LevelData<FArrayBox>& a_ccJ,
           const int                   a_p);

    // Computes the p-norm over a_data's grids using the formula
    //   ||a_data||_p = { Int[|a_data|^p dV] / Int[dV] } ^ (1/p)
    // where Int = the integral over a_data's grids and dV = dXi*dEta*dZeta.
    // After some simplification, this function actually uses
    //   ||a_data||_p = { Sum[|a_data|^p] / Number of cells } ^ (1/p)
    // where Sum = the sum over a_data's grids.
    //
    // So that the integral (sum) doesn't improperly scale the boundary cells/nodes,
    // this function will convert a_data to cell-centered before summing.
    //
    // This version is supposed to be convenient, so it will ignore NULL levels.
    static Real
    pNorm (const LevelData<FArrayBox>& a_data,
           const int                   a_p);

    // Computes the p-norm over all valid cells using the formula
    //   ||a_data||_p = { Int[|a_data|^p dV] / Int[dV] } ^ (1/p)
    // where Int = the integral over all valid cells and dV = J * dXi*dEta*dZeta.
    // After some simplification, this function actually uses
    //   ||a_data||_p = { Sum[|a_data|^p * J] / Sum[J] } ^ (1/p)
    // where Sum = the sum over a_data's valid regions.
    //
    // So that the integral (sum) doesn't improperly scale the boundary cells/nodes,
    // this function will convert a_data to cell-centered before summing.
    //
    // Warnings:
    //  1. a_ccJ must be cell-centered, as the name suggests.
    //  2. Behavior is undefined if a_data and a_ccJ are defined over different grids.
    static Real
    pNorm (const Vector<LevelData<FArrayBox>*>& a_data,
           const Interval&                      a_dataInterv,
           const Vector<LevelData<FArrayBox>*>& a_ccJ,
           const int                            a_p);

    // Computes the p-norm over all valid cells using the formula
    //   ||a_data||_p = { Int[|a_data|^p dV] / Int[dV] } ^ (1/p)
    // where Int = the integral over all valid cells and dV = J * dXi*dEta*dZeta.
    // After some simplification, this function actually uses
    //   ||a_data||_p = { Sum[|a_data|^p * J] / Sum[J] } ^ (1/p)
    // where Sum = the sum over a_data's valid regions
    // and J = 1 on all levels.
    //
    // So that the integral (sum) doesn't improperly scale the boundary cells/nodes,
    // this function will convert a_data to cell-centered before summing.
    static Real
    pNorm (const Vector<LevelData<FArrayBox>*>& a_data,
           const Interval&                      a_dataInterv,
           const int                            a_p);

    // Computes the p-norm over all valid cells using the formula
    //   ||a_data||_p = { Int[|a_data|^p dV] / Int[dV] } ^ (1/p)
    // where Int = the integral over all valid cells and dV = J * dXi*dEta*dZeta.
    // After some simplification, this function actually uses
    //   ||a_data||_p = { Sum[|a_data|^p * J] / Sum[J] } ^ (1/p)
    // where Sum = the sum over a_data's valid regions
    // and J = 1 on all levels.
    //
    // So that the integral (sum) doesn't improperly scale the boundary cells/nodes,
    // this function will convert a_data to cell-centered before summing.
    static Real
    pNorm (const Vector<LevelData<FArrayBox>*>& a_data,
           const int                            a_p);

    // -------------------------------------------------------------------------
    // p-norm (FC versions...)
    // -------------------------------------------------------------------------

    // FC AMR version.
    static RealVect
    pNorm (const Vector<LevelData<FluxBox>*>&   a_data,
           const Interval&                      a_dataInterv,
           const Vector<LevelData<FArrayBox>*>& a_ccJ,
           const int                            a_p);

    // FC single-level version
    static RealVect
    pNorm (const LevelData<FluxBox>&   a_data,
           const Interval&             a_dataInterv,
           const LevelData<FArrayBox>& a_ccJ,
           const BoxLayout&            a_layout,
           const BoxLayout*            a_fineLayoutPtr,
           const IntVect&              a_refRatio,
           const int                   a_p);

    // FC single-level version.
    static RealVect
    pNorm (const LevelData<FluxBox>& a_data,
           const int                 a_p);
};


#endif //!___Analysis_H__INCLUDED___
