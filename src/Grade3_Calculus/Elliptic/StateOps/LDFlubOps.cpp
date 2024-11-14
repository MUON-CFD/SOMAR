#include "LDFlubOps.H"
#include "Comm.H"
#include "FABAlgebra.H"
#include "AnisotropicRefinementTools.H"
#include "Debug.H"


#ifndef NDEBUG
#   define nanCheck(x) checkForValidNAN(x)
#else
#   define nanCheck(x)
#endif


namespace Elliptic {

using StateType  = LevelData<FluxBox>;
using TraitsType = StateTraits<LevelData<FluxBox>>;
using CopierType = StateTraits<LevelData<FluxBox>>::CopierType;
using GridsType  = StateTraits<LevelData<FluxBox>>::GridsType;


// -----------------------------------------------------------------------------
// Create data holder a_lhs that mirrors a_rhs. You do not need to copy
// the data of a_rhs, just make a holder the same size.
// -----------------------------------------------------------------------------
void
StateOps<StateType, TraitsType>::create(StateType&       a_lhs,
                                        const StateType& a_rhs) const
{
    a_lhs.define(a_rhs.getBoxes(), a_rhs.nComp(), a_rhs.ghostVect());
    debugInitLevel(a_lhs);
}


// -----------------------------------------------------------------------------
// Opposite of create -- perform any operations required before lhs goes
// out of scope. In general, the only time this needs to be defined in
// a derived class is if the create() function called new. Otherwise, the
// default no-op function is sufficient.
// -----------------------------------------------------------------------------
void
StateOps<StateType, TraitsType>::clear(StateType& a_lhs) const
{
    a_lhs.clear();
}


// -----------------------------------------------------------------------------
// Set a_lhs equal to a_rhs.
// If StateType does not have copiers, just send in a_copierPtr = nullptr.
// -----------------------------------------------------------------------------
void
StateOps<StateType, TraitsType>::assign(StateType&        a_lhs,
                                        const StateType&  a_rhs,
                                        const CopierType* a_copierPtr) const
{
    CH_assert(a_lhs.nComp() == a_rhs.nComp());

    if (a_rhs.getBoxes().compatible(a_lhs.getBoxes())) {
        this->assignLocal(a_lhs, a_rhs);
    } else {
        if (a_copierPtr) {
            a_rhs.copyTo(a_lhs, *a_copierPtr);
        } else {
            a_rhs.copyTo(a_lhs);
        }
    }
}


// -----------------------------------------------------------------------------
// Set a_lhs equal to a_rhs. Non-blocking version.
// a_lhs and a_rhs MUST be defined on compatible grids.
// -----------------------------------------------------------------------------
void
StateOps<StateType, TraitsType>::assignLocal(StateType&       a_lhs,
                                             const StateType& a_rhs) const
{
    CH_assert(a_lhs.getBoxes().compatible(a_rhs.getBoxes()));

    DataIterator dit = a_rhs.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        D_TERM(
        a_lhs[dit][0].copy(a_rhs[dit][0]);,
        a_lhs[dit][1].copy(a_rhs[dit][1]);,
        a_lhs[dit][2].copy(a_rhs[dit][2]);)
    }
}


// -----------------------------------------------------------------------------
// Compute and return the dot product of a_1 and a_2. In most contexts,
// this means return the sum over all data points of a_1*a_2.
// This does not use the metric tensor, so a_1 and a_2 need to be
// Cartesian-based.
// -----------------------------------------------------------------------------
Real
StateOps<StateType, TraitsType>::dotProduct(const StateType& a_1,
                                            const StateType& a_2) const
{
    nanCheck(a_1);
    nanCheck(a_2);

    CH_assert(a_1.getBoxes() == a_2.getBoxes());
    CH_assert(a_1.nComp() == a_2.nComp());

    const DisjointBoxLayout& grids = a_1.getBoxes();
    DataIterator             dit   = grids.dataIterator();
    Real                     val   = 0.0;

    for (dit.reset(); dit.ok(); ++dit) {
        for (int dir = 0; dir < SpaceDim; ++dir) {
            const FArrayBox& a1FAB = a_1[dit][dir];
            const FArrayBox& a2FAB = a_2[dit][dir];

            const Box validLo  = bdryLo(grids[dit], dir);
            const Box validMid = grids[dit].grow(dir, -1).surroundingNodes(dir);
            const Box validHi  = bdryHi(grids[dit], dir);

            val += a1FAB.dotProduct(a2FAB, validLo) * 0.5;
            val += a1FAB.dotProduct(a2FAB, validMid);
            val += a1FAB.dotProduct(a2FAB, validHi) * 0.5;
        }
    }

    Comm::reduce(val, MPI_SUM);

    return val;
}


// -----------------------------------------------------------------------------
// \brief
//  Compute the p-Norm of a_x. If a_p = 0, then we compute an inf-norm.
// \details
//  This function must perform all MPI communication.
//  If you want the norm on this level to be comparable to the norms on
//  other levels, then a_powScale must be chosen carefully. For example,
//  with cell-averaged data, a_powScale = dx * dy * dz would work.
//  It is used like this: ( Sum[|a_x|^p] * a_powScale )^(1/p) for p > 0
//  and not used for p = 0 (the inf-norm).
// -----------------------------------------------------------------------------
Real
StateOps<StateType, TraitsType>::norm(const StateType& a_x,
                                      const int        a_p,
                                      const Real       a_powScale) const
{
    nanCheck(a_x);

    const int                numComps = a_x.nComp();
    const DisjointBoxLayout& grids    = a_x.getBoxes();
    DataIterator             dit      = a_x.dataIterator();
    Real                     retVal   = 0.0;
    Real                     boxVal;

    for (dit.reset(); dit.ok(); ++dit) {
        for (int d = 0; d < SpaceDim; ++d) {
            Box fcValid = grids[dit];
            fcValid.surroundingNodes(d);

            boxVal = a_x[dit][d].norm(fcValid, a_p, 0, numComps);

            if (a_p == 0) {
                retVal = max(retVal, boxVal);
            } else {
                retVal += pow(boxVal, a_p);
            }
        }
    }

    if (a_p == 0) {
        Comm::reduce(retVal, MPI_MAX);
    } else {
        Comm::reduce(retVal, MPI_SUM);
        retVal = pow(retVal * a_powScale, 1.0 / Real(a_p));
    }

    return retVal;
}


// -----------------------------------------------------------------------------
// Increment by scaled amount (a_lhs += a_scale*a_x).
// -----------------------------------------------------------------------------
void
StateOps<StateType, TraitsType>::incr(StateType&       a_lhs,
                                      const StateType& a_x,
                                      Real             a_scale) const
{
    nanCheck(a_lhs);
    nanCheck(a_x);

    CH_assert(a_x.getBoxes() == a_lhs.getBoxes());
    CH_assert(a_x.nComp() == a_lhs.nComp());

    DataIterator dit = a_lhs.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        D_TERM(
        a_lhs[dit][0].plus(a_x[dit][0], a_scale);,
        a_lhs[dit][1].plus(a_x[dit][1], a_scale);,
        a_lhs[dit][2].plus(a_x[dit][2], a_scale);)
    }

    nanCheck(a_lhs);
}


// -----------------------------------------------------------------------------
// Set input to a scaled sum (a_lhs = a_a*a_x + a_b*a_y).
// -----------------------------------------------------------------------------
void
StateOps<StateType, TraitsType>::axby(StateType&       a_lhs,
                                      const StateType& a_x,
                                      const StateType& a_y,
                                      Real             a_a,
                                      Real             a_b) const
{
    nanCheck(a_x);
    nanCheck(a_y);

    CH_assert(a_x.getBoxes() == a_lhs.getBoxes());
    CH_assert(a_y.getBoxes() == a_lhs.getBoxes());
    CH_assert(a_x.nComp() == a_lhs.nComp());
    CH_assert(a_y.nComp() == a_lhs.nComp());

    debugInitLevel(a_lhs);

    DataIterator dit = a_lhs.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        D_TERM(
        FABAlgebra::axby(a_lhs[dit][0], a_x[dit][0], a_y[dit][0], a_a, a_b);,
        FABAlgebra::axby(a_lhs[dit][1], a_x[dit][1], a_y[dit][1], a_a, a_b);,
        FABAlgebra::axby(a_lhs[dit][2], a_x[dit][2], a_y[dit][2], a_a, a_b);)
    }

    nanCheck(a_lhs);
}


// -----------------------------------------------------------------------------
// Multiply the input by a given scale (a_lhs *= a_scale).
// -----------------------------------------------------------------------------
void
StateOps<StateType, TraitsType>::scale(StateType&  a_lhs,
                                       const Real& a_scale) const
{
    nanCheck(a_lhs);

    DataIterator dit = a_lhs.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        a_lhs[dit] *= a_scale;
    }

    nanCheck(a_lhs);
}


// -----------------------------------------------------------------------------
// Set a_lhs to zero.
// -----------------------------------------------------------------------------
void
StateOps<StateType, TraitsType>::setToZero(StateType& a_lhs) const
{
    DataIterator dit = a_lhs.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        a_lhs[dit].setVal(0.0);
    }
}


// -----------------------------------------------------------------------------
// Create a coarsened version of a_fine.
// You do not need to fill a_crse with data, just define it properly.
// -----------------------------------------------------------------------------
void
StateOps<StateType, TraitsType>::createCoarsened(
    StateType&       a_crse,
    const StateType& a_fine,
    const IntVect&   a_refRatio) const
{
    const GridsType& fineGrids = a_fine.getBoxes();

    GridsType crseGrids;
    ::coarsen(crseGrids, fineGrids, a_refRatio);

    a_crse.define(crseGrids, a_fine.nComp(), a_fine.ghostVect());
    debugInitLevel(a_crse);
}


// -----------------------------------------------------------------------------
// Create a copier for a_rhs.copyTo(a_lhs).
// If your StateType does not support copiers, leave the function empty.
// -----------------------------------------------------------------------------
void
StateOps<StateType, TraitsType>::buildCopier(CopierType&      a_copier,
                                             const StateType& a_lhs,
                                             const StateType& a_rhs) const
{
    const DisjointBoxLayout& srcGrids  = a_rhs.getBoxes();
    const DisjointBoxLayout& destGrids = a_lhs.getBoxes();
    const ProblemDomain&     domain    = destGrids.physDomain();
    CH_assert(domain == srcGrids.physDomain());

    for (int d = 0; d < SpaceDim; ++d) {
        a_copier[d].define(srcGrids, destGrids, domain, a_lhs.ghostVect(), d);
    }
}


// -----------------------------------------------------------------------------
// Create a copier for a_lhs --> a_rhs.
// If your StateType does not support copiers, leave the function empty.
// -----------------------------------------------------------------------------
void
StateOps<StateType, TraitsType>::buildReverseCopier(
    CopierType&       a_reverseCopier,
    const CopierType& a_copier,
    const StateType&  /*a_lhs*/,
    const StateType&  /*a_rhs*/) const
{
    a_reverseCopier = a_copier;
    for (int d = 0; d < SpaceDim; ++d) {
        a_reverseCopier[d].reverse();
    }
}


// -----------------------------------------------------------------------------
// Set lhs = 0 where covered by rhs.
// rhs can be clobbered, it is temporary storage.
// If you don't have a copier for any reason, just send in nullptr.
// -----------------------------------------------------------------------------
void
StateOps<StateType, TraitsType>::zeroCovered(
    StateType& a_lhs, StateType& a_rhs, const CopierType* a_copierPtr) const
{
    this->setToZero(a_rhs);
    this->assign(a_lhs, a_rhs, a_copierPtr);
}


// -----------------------------------------------------------------------------
// Compute the max norm. Don't worry about MPI communication, that is left
// to the caller. This function should be non-blocking.
// -----------------------------------------------------------------------------
Real
StateOps<StateType, TraitsType>::localMaxNorm(const StateType& a_x) const
{
    nanCheck(a_x);

    const int                p        = 0;
    const int                numComps = a_x.nComp();
    const DisjointBoxLayout& grids    = a_x.getBoxes();
    DataIterator             dit      = a_x.dataIterator();
    Real                     retVal   = 0.0;
    Real                     boxVal;
    Box                      fcValid;

    for (dit.reset(); dit.ok(); ++dit) {
        for (int d = 0; d < SpaceDim; ++d) {
            fcValid = grids[dit];
            fcValid.surroundingNodes(d);

            boxVal = a_x[dit][d].norm(fcValid, p, 0, numComps);
            retVal = max(retVal, boxVal);
        }
    }

    return retVal;
}


}; // namespace Elliptic
