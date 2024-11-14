#include "Analysis.H"
#include "AnalysisF_F.H"
#include "AnisotropicRefinementTools.H"
#include "Convert.H"
#include "Comm.H"
#include "SOMAR_Constants.H"
#include "Debug.H"
#include <limits>



// -----------------------------------------------------------------------------
// Computes the maximum over each cell/node in a_region.
// Warnings:
//  1. Does not take the abs value!
//  2. Behavior is undefined if a_data does not contain a_region.
// -----------------------------------------------------------------------------
Real
Analysis::max (const FArrayBox& a_data,
               const Box&       a_region,
               const Interval   a_interv)
{
    CH_assert(a_data.box().sameType(a_region));
    CH_assert(a_data.box().contains(a_region));

    int comp = a_interv.begin();
    Real tmpVal;
    Real retVal = a_data.max(a_region, comp);
    ++comp;

    while (comp <= a_interv.end()) {
        tmpVal = a_data.max(a_region, comp);
        if (tmpVal > retVal) {
            retVal = tmpVal;
        }
        ++comp;
    }

    return retVal;
}


// -----------------------------------------------------------------------------
// Computes the maximum over each cell/node in a_layout. If you want to
// include the ghosts of a_data, simply include those ghosts in the layout.
// This is the most general single-level version. It can ignore cells/nodes
// that are covered by a finer BoxLayout.
// Warnings:
//  1. Does not take the abs value!
//  2. Behavior is undefined if a_data does not contain a_layout.
//  3. a_data must be compatible with a_layout.
// -----------------------------------------------------------------------------
Real
Analysis::max (const LevelData<FArrayBox>& a_data,
               const BoxLayout&            a_layout,
               const Interval&             a_interv,
               const BoxLayout*            a_fineLayoutPtr,
               const IntVect&              a_refRatio)
{
    CH_assert(a_data.getBoxes().compatible(a_layout));

    DataIterator dit = a_data.dataIterator();
    const int startComp = a_interv.begin();
    const int numComp = a_interv.size();

    // Coarsen the fine grids.
    BoxLayout invalidLayout;
    if (a_fineLayoutPtr != nullptr) {
        coarsen(invalidLayout, *a_fineLayoutPtr, a_refRatio);
    }

    // Begin with the smallest number representable by a Real.
    Real maxVal = minReal;
    Real tmpVal;

    for (dit.reset(); dit.ok(); ++dit) {
        const FArrayBox& dataFAB = a_data[dit];

        // Figure out where we will be computing the local max
        Box region = a_layout[dit];
        region.convert(dataFAB.box().type());

        // Do we need to eliminate invalid regions?
        if (a_fineLayoutPtr == nullptr) {
            // No finer grids, just compute the local max
            tmpVal = Analysis::max(dataFAB, region, a_interv);
            if (tmpVal > maxVal) {
                maxVal = tmpVal;
            }

        } else {
            // Finer grids exist...
            // First, make a copy of the data.
            FArrayBox validDataFAB(region, startComp);
            validDataFAB.copy(dataFAB, region, startComp, region, 0, numComp);

            // Then eliminate invalid data.
            LayoutIterator lit = invalidLayout.layoutIterator();
            for (lit.reset(); lit.ok(); ++lit) {
                Box invalidBox = invalidLayout[lit];
                validDataFAB.setVal(minReal, invalidBox, 0, numComp);
            }

            // Now compute the local max
            tmpVal = Analysis::max(validDataFAB, region, validDataFAB.interval());
            if (tmpVal > maxVal) {
                maxVal = tmpVal;
            }
        }
    }

    // Compute global max (this is where the MPI communication happens)
    Comm::reduce(maxVal, MPI_MAX);

    return maxVal;
}


// -----------------------------------------------------------------------------
// Computes the maximum over each cell/node in a_layout. If you want to
// include the ghosts of a_data, simply include those ghosts in the layout.
// Warnings:
//  1. Does not take the abs value!
//  2. Behavior is undefined if a_data does not contain a_layout.
//  3. a_data must be compatible with a_layout.
// -----------------------------------------------------------------------------
Real
Analysis::max (const LevelData<FArrayBox>& a_data,
               const BoxLayout&            a_layout,
               const Interval              a_interv)
{
    return Analysis::max(a_data, a_layout, a_interv,
                         nullptr, IntVect::Unit);
}


// -----------------------------------------------------------------------------
// Computes the maximum over each cell/node in the valid cells of a_data.
// Warning: Does not take the abs value!
// -----------------------------------------------------------------------------
Real
Analysis::max (const LevelData<FArrayBox>& a_data,
               const Interval              a_interv)
{
    return Analysis::max(a_data, a_data.getBoxes(), a_interv,
                         nullptr, IntVect::Unit);
}


// -----------------------------------------------------------------------------
// Computes the maximum over each valid cell/node in a_data.
// Warning: Does not take the abs value!
// -----------------------------------------------------------------------------
Real
Analysis::max (const Vector<LevelData<FArrayBox>*>& a_data,
               const int                            a_lmin,
               const int                            a_lmax,
               const Interval&                      a_interv,
               const Vector<IntVect>&               a_refRatio)
{
    // Sanity checks
    CH_assert(0 <= a_lmin);
    CH_assert(a_lmin <= a_lmax);
    CH_assert(a_lmax < static_cast<int>(a_data.size()));

    Real maxVal = minReal;
    Real tmpVal;

    for (int lev = a_lmin; lev < a_lmax; ++lev) {
        tmpVal = Analysis::max(*a_data[lev],
                               a_data[lev]->getBoxes(),
                               a_interv,
                               &(a_data[lev+1]->getBoxes()),
                               a_refRatio[lev+1]);

        if (tmpVal > maxVal) {
            maxVal = tmpVal;
        }
    }

    tmpVal = Analysis::max(*a_data[a_lmax],
                           a_data[a_lmax]->getBoxes(),
                           a_interv,
                           nullptr,
                           IntVect::Unit);

    if (tmpVal > maxVal) {
        maxVal = tmpVal;
    }

    return maxVal;
}


// -----------------------------------------------------------------------------
// Computes the minimum over each cell/node in a_region.
// Warnings:
//  1. Does not take the abs value!
//  2. Behavior is undefined if a_data does not contain a_region.
// -----------------------------------------------------------------------------
Real
Analysis::min (const FArrayBox& a_data,
               const Box&       a_region,
               const Interval   a_interv)
{
    CH_assert(a_data.box().sameType(a_region));
    CH_assert(a_data.box().contains(a_region));

    int comp = a_interv.begin();
    Real tmpVal;
    Real retVal = a_data.min(a_region, comp);
    ++comp;

    while (comp <= a_interv.end()) {
        tmpVal = a_data.min(a_region, comp);
        if (tmpVal < retVal) {
            retVal = tmpVal;
        }
        ++comp;
    }

    return retVal;
}


// -----------------------------------------------------------------------------
// Computes the minimum over each cell/node in a_layout. If you want to
// include the ghosts of a_data, simply include those ghosts in the layout.
// This is the most general single-level version. It can ignore cells/nodes
// that are covered by a finer BoxLayout.
// Warnings:
//  1. Does not take the abs value!
//  2. Behavior is undefined if a_data does not contain a_layout.
//  3. a_data must be compatible with a_layout.
// -----------------------------------------------------------------------------
Real
Analysis::min (const LevelData<FArrayBox>& a_data,
               const BoxLayout&            a_layout,
               const Interval&             a_interv,
               const BoxLayout*            a_fineLayoutPtr,
               const IntVect&              a_refRatio)
{
    TEST();
    CH_assert(a_data.getBoxes().compatible(a_layout));

    DataIterator dit = a_data.dataIterator();
    const int startComp = a_interv.begin();
    const int numComp = a_interv.size();

    // Coarsen the fine grids.
    BoxLayout invalidLayout;
    if (a_fineLayoutPtr != nullptr) {
        coarsen(invalidLayout, *a_fineLayoutPtr, a_refRatio);
    }

    // Begin with the largest number representable by a Real.
    Real minVal = maxReal;
    Real tmpVal;

    for (dit.reset(); dit.ok(); ++dit) {
        const FArrayBox& dataFAB = a_data[dit];

        // Figure out where we will be computing the local min
        Box region = a_layout[dit];
        region.convert(dataFAB.box().type());

        // Do we need to eliminate invalid regions?
        if (a_fineLayoutPtr == nullptr) {
            // No finer grids, just compute the local min
            tmpVal = Analysis::min(dataFAB, region, a_interv);
            if (tmpVal < minVal) {
                minVal = tmpVal;
            }

        } else {
            // Finer grids exist...
            // First, make a copy of the data.
            FArrayBox validDataFAB(region, startComp);
            validDataFAB.copy(dataFAB, region, startComp, region, 0, numComp);

            // Then eliminate invalid data.
            LayoutIterator lit = invalidLayout.layoutIterator();
            for (lit.reset(); lit.ok(); ++lit) {
                Box invalidBox = invalidLayout[lit];
                validDataFAB.setVal(maxReal, invalidBox, 0, numComp);
            }

            // Now compute the local min
            tmpVal = Analysis::min(validDataFAB, region, validDataFAB.interval());
            if (tmpVal < minVal) {
                minVal = tmpVal;
            }
        }
    }

    // Compute global min (this is where the MPI communication happens)
    Comm::reduce(minVal, MPI_MIN);

    return minVal;
}


// -----------------------------------------------------------------------------
// Computes the minimum over each cell/node in a_layout. If you want to
// include the ghosts of a_data, simply include those ghosts in the layout.
// Warnings:
//  1. Does not take the abs value!
//  2. Behavior is undefined if a_data does not contain a_layout.
//  3. a_data must be compatible with a_layout.
// -----------------------------------------------------------------------------
Real
Analysis::min (const LevelData<FArrayBox>& a_data,
               const BoxLayout&            a_layout,
               const Interval              a_interv)
{
    return Analysis::min(a_data, a_layout, a_interv,
                         nullptr, IntVect::Unit);
}


// -----------------------------------------------------------------------------
// Computes the minimum over each cell/node in the valid cells of a_data.
// Warning: Does not take the abs value!
// -----------------------------------------------------------------------------
Real
Analysis::min (const LevelData<FArrayBox>& a_data,
               const Interval              a_interv)
{
    return Analysis::min(a_data, a_data.getBoxes(), a_interv,
                         nullptr, IntVect::Unit);
}


// -----------------------------------------------------------------------------
// Computes the minimum over each valid cell/node in a_data.
// Warning: Does not take the abs value!
// -----------------------------------------------------------------------------
Real
Analysis::min (const Vector<LevelData<FArrayBox>*>& a_data,
               const int                            a_lmin,
               const int                            a_lmax,
               const Interval&                      a_interv,
               const Vector<IntVect>&               a_refRatio)
{
    // Sanity checks
    CH_assert(0 <= a_lmin);
    CH_assert(a_lmin <= a_lmax);
    CH_assert(a_lmax < static_cast<int>(a_data.size()));

    Real minVal = maxReal;
    Real tmpVal;

    for (int lev = a_lmin; lev < a_lmax; ++lev) {
        tmpVal = Analysis::min(*a_data[lev],
                               a_data[lev]->getBoxes(),
                               a_interv,
                               &(a_data[lev+1]->getBoxes()),
                               a_refRatio[lev+1]);

        if (tmpVal < minVal) {
            minVal = tmpVal;
        }
    }

    tmpVal = Analysis::min(*a_data[a_lmax],
                           a_data[a_lmax]->getBoxes(),
                           a_interv,
                           nullptr,
                           IntVect::Unit);

    if (tmpVal < minVal) {
        minVal = tmpVal;
    }

    return minVal;
}


// -----------------------------------------------------------------------------
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
// -----------------------------------------------------------------------------
void
Analysis::pNorm_local_Accumulator (Real&            a_normAccum,
                                   Real&            a_volAccum,
                                   const FArrayBox& a_data,
                                   const Interval&  a_dataInterv,
                                   const FArrayBox& a_ccJ,
                                   const Box&       a_ccRegion,
                                   const int        a_p)
{
    // Sanity checks
    CH_assert(a_ccJ.box().type() == IntVect::Zero);
    CH_assert(a_ccRegion.type() == IntVect::Zero);

    CH_assert(a_ccJ.box().contains(a_ccRegion));
    CH_assert(enclosedCells(a_data.box()).contains(a_ccRegion));

    CH_assert(a_p >= 0);

    const int startComp = a_dataInterv.begin();
    const int numComps = a_dataInterv.size();

    // Convert to CC if needed.
    FArrayBox ccData;
    if (a_data.box().type() != IntVect::Zero) {
        // First, compute the region that contains data.
        Box dataBox = a_ccRegion;
        dataBox.convert(a_data.box().type());

        // Recenter J and compute J*data over that region.
        FArrayBox Jdata(dataBox, numComps);
        Convert::Simple(Jdata, dataBox, a_ccJ);
        Jdata.mult(a_data, dataBox, startComp, 0, numComps);

        // Convert J*data to CC.
        ccData.define(a_ccRegion, numComps);
        Convert::Simple(ccData, a_ccRegion, Jdata);

        // Remove the J scaling.
        ccData.divide(a_ccJ);

    } else {
        ccData.define(a_dataInterv, (FArrayBox&)a_data);
    }

    // Accumulate Sum[|a_data|^p] and Sum[J].
    for (int comp = 0; comp < numComps; ++comp) {
        FORT_ANALYSIS_PNORMACCUMULATOR (
            CHF_REAL(a_normAccum),
            CHF_CONST_FRA1(ccData, comp),
            CHF_CONST_FRA1(a_ccJ, 0),
            CHF_BOX(a_ccRegion),
            CHF_CONST_INT(a_p));
    }
    a_volAccum += a_ccJ.sum(a_ccRegion, 0);
}


// -----------------------------------------------------------------------------
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
// -----------------------------------------------------------------------------
void
Analysis::pNorm_local_Accumulator (Real&                       a_normAccum,
                                   Real&                       a_volAccum,
                                   const LevelData<FArrayBox>& a_data,
                                   const Interval&             a_dataInterv,
                                   const LevelData<FArrayBox>& a_ccJ,
                                   const BoxLayout&            a_layout,
                                   const BoxLayout*            a_fineLayoutPtr,
                                   const IntVect&              a_refRatio,
                                   const int                   a_p)
{
    CH_assert(a_data.getBoxes().compatible(a_layout));
    CH_assert(a_ccJ .getBoxes().compatible(a_layout));

    DataIterator dit = a_data.dataIterator();
    const int startComp = a_dataInterv.begin();
    const int numComps = a_dataInterv.size();

    // Coarsen the fine grids.
    BoxLayout invalidLayout;
    if (a_fineLayoutPtr != nullptr) {
        coarsen(invalidLayout, *a_fineLayoutPtr, a_refRatio);
    }

    for (dit.reset(); dit.ok(); ++dit) {
        const FArrayBox& dataFAB = a_data[dit];
        const FArrayBox& ccJFAB = a_ccJ[dit];
        const Box& ccRegion = a_layout[dit];

        // Do we need to eliminate invalid regions?
        if (a_fineLayoutPtr == nullptr) {
            // No finer grids, just accumulate.
            Analysis::pNorm_local_Accumulator(a_normAccum,
                                              a_volAccum,
                                              dataFAB,
                                              a_dataInterv,
                                              ccJFAB,
                                              ccRegion,
                                              a_p);
        } else {
            // Finer grids exist...

            // First, figure out where the data lives.
            Box region = ccRegion;
            region.convert(dataFAB.box().type());
            CH_assert(dataFAB.box().contains(region));

            // Next, make a copy of the data.
            FArrayBox validDataFAB(region, numComps);
            validDataFAB.copy(dataFAB, region, startComp, region, 0, numComps);

            // Then eliminate invalid data.
            LayoutIterator lit = invalidLayout.layoutIterator();
            for (lit.reset(); lit.ok(); ++lit) {
                Box invalidBox = invalidLayout[lit];
                invalidBox.convert(region.type());
                invalidBox &= region;

                validDataFAB.setVal(0.0, invalidBox, 0, numComps);
            }

            // Now compute the local norm
            Analysis::pNorm_local_Accumulator(a_normAccum,
                                              a_volAccum,
                                              validDataFAB,
                                              validDataFAB.interval(),
                                              ccJFAB,
                                              ccRegion,
                                              a_p);
        }
    }
}


// -----------------------------------------------------------------------------
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
//  2. Behavior is undefined if a_data or a_ccJ do not contain a_ccRegion.
// -----------------------------------------------------------------------------
Real
Analysis::pNorm (const FArrayBox& a_data,
                 const Interval&  a_dataInterv,
                 const FArrayBox& a_ccJ,
                 const Box&       a_ccRegion,
                 const int        a_p)
{
    Real pNorm = 0.0;
    Real vol = 0.0;
    Analysis::pNorm_local_Accumulator(pNorm,
                                      vol,
                                      a_data,
                                      a_dataInterv,
                                      a_ccJ,
                                      a_ccRegion,
                                      a_p);

    if (a_p >= 1) {
        // Divide by Sum[J]
        pNorm /= vol;

        // Take the pth-root.
        if (a_p >= 2) {
            pNorm = pow(pNorm, 1.0 / a_p);
        }
    }

    return pNorm;
}


// -----------------------------------------------------------------------------
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
// -----------------------------------------------------------------------------
Real
Analysis::pNorm (const LevelData<FArrayBox>& a_data,
                 const Interval&             a_dataInterv,
                 const LevelData<FArrayBox>& a_ccJ,
                 const BoxLayout&            a_layout,
                 const BoxLayout*            a_fineLayoutPtr,
                 const IntVect&              a_refRatio,
                 const int                   a_p)
{
    // Accumulate from 0.
    Real pNorm = 0.0;
    Real vol = 0.0;
    Analysis::pNorm_local_Accumulator(pNorm,
                                      vol,
                                      a_data,
                                      a_dataInterv,
                                      a_ccJ,
                                      a_layout,
                                      a_fineLayoutPtr,
                                      a_refRatio,
                                      a_p);
    // Reduce.
    if (a_p == 0) {
        Comm::reduce(pNorm, MPI_MAX);
        return pNorm;
    } else {
        Comm::reduce(pNorm, MPI_SUM);
        Comm::reduce(vol  , MPI_SUM);
    }

    // Divide by Sum[J].
    pNorm /= vol;

    // Take the pth-root.
    if (a_p >= 2) {
        pNorm = pow(pNorm, 1.0 / a_p);
    }

    return pNorm;
}


// -----------------------------------------------------------------------------
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
// -----------------------------------------------------------------------------
Real
Analysis::pNorm (const LevelData<FArrayBox>& a_data,
                 const Interval&             a_dataInterv,
                 const LevelData<FArrayBox>& a_ccJ,
                 const int                   a_p)
{
    return Analysis::pNorm(a_data,
                           a_dataInterv,
                           a_ccJ,
                           a_data.getBoxes(),
                           nullptr,
                           IntVect::Unit,
                           a_p);
}


// -----------------------------------------------------------------------------
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
// -----------------------------------------------------------------------------
Real
Analysis::pNorm (const LevelData<FArrayBox>& a_data,
                 const LevelData<FArrayBox>& a_ccJ,
                 const int                   a_p)
{
    return Analysis::pNorm(a_data,
                           a_data.interval(),
                           a_ccJ,
                           a_data.getBoxes(),
                           nullptr,
                           IntVect::Unit,
                           a_p);
}


// -----------------------------------------------------------------------------
// Computes the p-norm over a_data's grids using the formula
//   ||a_data||_p = { Int[|a_data|^p dV] / Int[dV] } ^ (1/p)
// where Int = the integral over a_data's grids and dV = dXi*dEta*dZeta.
// After some simplification, this function actually uses
//   ||a_data||_p = { Sum[|a_data|^p] / Number of cells } ^ (1/p)
// where Sum = the sum over a_data's grids.
//
// So that the integral (sum) doesn't improperly scale the boundary cells/nodes,
// this function will convert a_data to cell-centered before summing.
// -----------------------------------------------------------------------------
Real
Analysis::pNorm (const LevelData<FArrayBox>& a_data,
                 const int                   a_p)
{
    const DisjointBoxLayout& grids = a_data.getBoxes();
    DataIterator dit = grids.dataIterator();

    LevelData<FArrayBox> ccJ(grids, 1, a_data.ghostVect() + IntVect::Unit);
    for (dit.reset(); dit.ok(); ++dit) {
        ccJ[dit].setVal(1.0);
    }

    return Analysis::pNorm(a_data,
                           a_data.interval(),
                           ccJ,
                           a_data.getBoxes(),
                           nullptr,
                           IntVect::Unit,
                           a_p);
}


// -----------------------------------------------------------------------------
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
// -----------------------------------------------------------------------------
Real
Analysis::pNorm (const Vector<LevelData<FArrayBox>*>& a_data,
                 const Interval&                      a_dataInterv,
                 const Vector<LevelData<FArrayBox>*>& a_ccJ,
                 const int                            a_p)
{
    const int numLevels = a_data.size();
    CH_assert(static_cast<int>(a_ccJ.size()) >= numLevels);

    // Compute integrals on each level and add up local results.
    Real pNorm = 0.0;
    Real vol = 0.0;
    for (int l = 0; l < numLevels; ++l) {
        // Sanity checks
        CH_assert(a_data[l] != nullptr);
        CH_assert(a_ccJ[l] != nullptr);
        CH_assert(a_data[l]->getBoxes().compatible(a_ccJ[l]->getBoxes()));

        // This level's grids
        const DisjointBoxLayout& grids = a_data[l]->getBoxes();

        // Fine level's info.
        const DisjointBoxLayout* fineGridsPtr = nullptr;
        IntVect refRatio = IntVect::Unit;
        if (l < numLevels-1) {
            fineGridsPtr = &(a_data[l+1]->getBoxes());

            const Box& crseDomBox = grids.physDomain().domainBox();
            const Box& fineDomBox = fineGridsPtr->physDomain().domainBox();
            refRatio = calculateRefinementRatio(crseDomBox, fineDomBox);
        }

        // Accumulate!
        pNorm_local_Accumulator(pNorm,
                                vol,
                                *a_data[l],
                                a_dataInterv,
                                *a_ccJ[l],
                                grids,
                                fineGridsPtr,
                                refRatio,
                                a_p);
    }

    // Reduce.
    if (a_p == 0) {
        Comm::reduce(pNorm, MPI_MAX);
        return pNorm;
    } else {
        Comm::reduce(pNorm, MPI_SUM);
        Comm::reduce(vol  , MPI_SUM);
    }

    // Divide by Sum[J].
    pNorm /= vol;

    // Take the pth-root.
    if (a_p >= 2) {
        pNorm = pow(pNorm, 1.0 / a_p);
    }

    return pNorm;
}


// -----------------------------------------------------------------------------
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
// -----------------------------------------------------------------------------
Real
Analysis::pNorm (const Vector<LevelData<FArrayBox>*>& a_data,
                 const Interval&                      a_dataInterv,
                 const int                            a_p)
{
    // Create ccJ = 1.
    const int numLevels = a_data.size();
    Vector<LevelData<FArrayBox>*> ccJ(numLevels, nullptr);
    for (int l = 0; l < numLevels; ++l) {
        const DisjointBoxLayout& grids = a_data[l]->getBoxes();
        DataIterator dit = grids.dataIterator();

        ccJ[l] = new LevelData<FArrayBox>(grids, 1);

        for (dit.reset(); dit.ok(); ++dit) {
            (*ccJ[l])[dit].setVal(1.0);
        }
    }

    // Compute norm.
    const Real pNorm = Analysis::pNorm(a_data, a_dataInterv, ccJ, a_p);

    // Free memory.
    for (int l = 0; l < numLevels; ++l) {
        delete ccJ[l];
    }

    return pNorm;
}


// -----------------------------------------------------------------------------
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
//
// This version is supposed to be convenient, so it will ignore nullptr levels.
// -----------------------------------------------------------------------------
Real
Analysis::pNorm (const Vector<LevelData<FArrayBox>*>& a_data,
                 const int                            a_p)
{
    // Trim the nullptrs.
    Vector<LevelData<FArrayBox>*> tv(0);
    for (unsigned int idx = 0; idx < a_data.size(); ++idx) {
        if (a_data[idx] != nullptr) {
            tv.push_back(a_data[idx]);
        }
    }

    // Compute norm of non-nullptr levels.
    return pNorm(tv, tv[0]->interval(), a_p);;
}


// FC versions...

// -----------------------------------------------------------------------------
// FC AMR version.
// -----------------------------------------------------------------------------
RealVect
Analysis::pNorm (const Vector<LevelData<FluxBox>*>&   a_data,
                 const Interval&                      a_dataInterv,
                 const Vector<LevelData<FArrayBox>*>& a_ccJ,
                 const int                            a_p)
{
    const int numComps = a_dataInterv.size();
    const int numLevels = a_data.size();
    CH_assert(static_cast<int>(a_ccJ.size()) >= numLevels);

    Vector<LevelData<FArrayBox>*> ccData(numLevels, nullptr);
    for (int l = 0; l < numLevels; ++l) {
        const DisjointBoxLayout& grids = a_data[l]->getBoxes();
        CH_assert(a_ccJ[l] != nullptr);
        CH_assert(a_ccJ[l]->getBoxes().compatible(grids));

        ccData[l] = new LevelData<FArrayBox>(grids, numComps);
    }

    RealVect pNorm = RealVect::Zero;
    for (int dir = 0; dir < SpaceDim; ++dir) {
        for (int l = 0; l < numLevels; ++l) {
            const DisjointBoxLayout& grids = a_data[l]->getBoxes();
            DataIterator dit = grids.dataIterator();

            for (dit.reset(); dit.ok(); ++dit) {
                Convert::Simple((*ccData[l])[dit],
                                Interval(0, numComps-1),
                                grids[dit],
                                (*a_data[l])[dit][dir],
                                a_dataInterv);
            }
        } // end loop over levels (l)

        pNorm[dir] = Analysis::pNorm(ccData,
                                     Interval(0, numComps-1),
                                     a_ccJ,
                                     a_p);
    } // end loop over directions (dir)

    return pNorm;
}


// -----------------------------------------------------------------------------
// Computes the p-norm over a_layout using the formula
//   ||a_data||_p = { Int[|a_data|^p dV] / Int[dV] } ^ (1/p)
// where Int = the integral over a_layout and dV = J * dXi*dEta*dZeta.
// After some simplification, this function actually uses
//   ||a_data||_p = { Sum[|a_data|^p * J] / Sum[J] } ^ (1/p)
// where Sum = the sum over a_layout.
//
// This function will perform the p-norm on all SpaceDim components
// separately and return the results in a RealVect.
//
// This is the most general single-level version for FC data. It can ignore
// cells/nodes that are covered by a finer BoxLayout.
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
// -----------------------------------------------------------------------------
RealVect
Analysis::pNorm (const LevelData<FluxBox>&   a_data,
                 const Interval&             /*a_dataInterv*/,
                 const LevelData<FArrayBox>& a_ccJ,
                 const BoxLayout&            a_layout,
                 const BoxLayout*            a_fineLayoutPtr,
                 const IntVect&              a_refRatio,
                 const int                   a_p)
{
    RealVect pNormVect;

    for (int dir = 0; dir < SpaceDim; ++dir) {
        LevelData<FArrayBox> ccData(a_data.getBoxes(), a_data.nComp(), a_data.ghostVect());

        DataIterator dit = a_data.dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            Convert::Simple(ccData[dit], a_data[dit][dir]);
        }

        pNormVect[dir] = Analysis::pNorm(ccData,
                                         ccData.interval(),
                                         a_ccJ,
                                         a_layout,
                                         a_fineLayoutPtr,
                                         a_refRatio,
                                         a_p);
    }

    return pNormVect;
}


// -----------------------------------------------------------------------------
// Computes the p-norm over a_data's grids using the formula
//   ||a_data||_p = { Int[|a_data|^p dV] / Int[dV] } ^ (1/p)
// where Int = the integral over a_data's grids and dV = dXi*dEta*dZeta.
// After some simplification, this function actually uses
//   ||a_data||_p = { Sum[|a_data|^p] / Number of cells } ^ (1/p)
// where Sum = the sum over a_data's grids.
//
// This function will perform the p-norm on all SpaceDim components
// separately and return the results in a RealVect.
//
// So that the integral (sum) doesn't improperly scale the boundary cells/nodes,
// this function will convert a_data to cell-centered before summing.
// -----------------------------------------------------------------------------
RealVect
Analysis::pNorm (const LevelData<FluxBox>& a_data,
                 const int                 a_p)
{
    const DisjointBoxLayout& grids = a_data.getBoxes();
    DataIterator dit = grids.dataIterator();

    LevelData<FArrayBox> ccJ(grids, 1, IntVect::Zero);
    for (dit.reset(); dit.ok(); ++dit) {
        ccJ[dit].setVal(1.0);
    }

    return Analysis::pNorm(a_data,
                           a_data.interval(),
                           ccJ,
                           a_data.getBoxes(),
                           nullptr,
                           IntVect::Unit,
                           a_p);
}
