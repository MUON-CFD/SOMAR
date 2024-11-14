#include "FABAlgebra.H"
#include "FABAlgebraF_F.H"


// -----------------------------------------------------------------------------
// Calculates dest = dest + p0 * p1
// -----------------------------------------------------------------------------
void
FABAlgebra::AddProd2 (FArrayBox&       a_dest,
                      const int        a_destComp,
                      const FArrayBox& a_p0,
                      const int        a_p0Comp,
                      const FArrayBox& a_p1,
                      const int        a_p1Comp,
                      const Box&       a_destBox)
{
    // Sanity checks
    CH_assert(a_dest.interval().contains(a_destComp));
    CH_assert(a_p0.interval().contains(a_p0Comp));
    CH_assert(a_p1.interval().contains(a_p1Comp));

    CH_assert(a_dest.box().type() == a_destBox.type());
    CH_assert(a_p0.box().type() == a_destBox.type());
    CH_assert(a_p1.box().type() == a_destBox.type());

    CH_assert(a_dest.contains(a_destBox));
    CH_assert(a_p0.contains(a_destBox));
    CH_assert(a_p1.contains(a_destBox));

    // Do it!
    FORT_ADDPROD2(
        CHF_FRA1(a_dest, a_destComp),
        CHF_CONST_FRA1(a_p0, a_p0Comp),
        CHF_CONST_FRA1(a_p1, a_p1Comp),
        CHF_BOX(a_destBox));
}


// -----------------------------------------------------------------------------
// Calculates dest = dest - p0 * p1
// -----------------------------------------------------------------------------
void
FABAlgebra::SubProd2 (FArrayBox&       a_dest,
                      const int        a_destComp,
                      const FArrayBox& a_p0,
                      const int        a_p0Comp,
                      const FArrayBox& a_p1,
                      const int        a_p1Comp,
                      const Box&       a_destBox)
{
    // Sanity checks
    CH_assert(a_dest.interval().contains(a_destComp));
    CH_assert(a_p0.interval().contains(a_p0Comp));
    CH_assert(a_p1.interval().contains(a_p1Comp));

    CH_assert(a_dest.box().type() == a_destBox.type());
    CH_assert(a_p0.box().type() == a_destBox.type());
    CH_assert(a_p1.box().type() == a_destBox.type());

    CH_assert(a_dest.contains(a_destBox));
    CH_assert(a_p0.contains(a_destBox));
    CH_assert(a_p1.contains(a_destBox));

    // Do it!
    FORT_SUBPROD2(
        CHF_FRA1(a_dest, a_destComp),
        CHF_CONST_FRA1(a_p0, a_p0Comp),
        CHF_CONST_FRA1(a_p1, a_p1Comp),
        CHF_BOX(a_destBox));
}


// -----------------------------------------------------------------------------
// Just like FArrayBox::axby, except this works on the intersection of
// all the FABs.
// -----------------------------------------------------------------------------
void
FABAlgebra::axby(FArrayBox&       a_dest,
                 const FArrayBox& a_x,
                 const FArrayBox& a_y,
                 const Real       a_a,
                 const Real       a_b)
{
    // Sanity checks
    CH_assert(a_x.nComp() == a_dest.nComp());
    CH_assert(a_y.nComp() == a_dest.nComp());

    // Compute intersection of FABs.
    Box region = a_dest.box();
    region &= a_x.box();
    region &= a_y.box();

    // Do it!
    FORT_AXBY(
        CHF_FRA(a_dest),
        CHF_CONST_FRA(a_x),
        CHF_CONST_FRA(a_y),
        CHF_CONST_REAL(a_a),
        CHF_CONST_REAL(a_b),
        CHF_BOX(region));

}


// -----------------------------------------------------------------------------
// Calculates dest = dest^pow
// -----------------------------------------------------------------------------
void
FABAlgebra::pow(FArrayBox& a_dest,
                const Real a_pow)
{
    Box destBox = a_dest.box();
    FORT_POWFAB(
        CHF_FRA(a_dest),
        CHF_BOX(destBox),
        CHF_CONST_REAL(a_pow));
}


// -----------------------------------------------------------------------------
// Calculates dest = dest^pow
// -----------------------------------------------------------------------------
void
FABAlgebra::pow(FArrayBox& a_dest,
                const int  a_destComp,
                const Box& a_destBox,
                const Real a_pow)
{
    // Sanity checks
    CH_assert(a_dest.interval().contains(a_destComp));
    CH_assert(a_dest.box().type() == a_destBox.type());
    CH_assert(a_dest.contains(a_destBox));

    // Do it!
    FArrayBox destAlias(Interval(a_destComp, a_destComp), a_dest);
    FORT_POWFAB(
        CHF_FRA(destAlias),
        CHF_BOX(a_destBox),
        CHF_CONST_REAL(a_pow));
}


// -----------------------------------------------------------------------------
// Computes EC = EC * (alpha + beta*Av[CC])
// -----------------------------------------------------------------------------
void
FABAlgebra::ECmultCC(FArrayBox&       a_ecFAB,
                     const int        a_ecFABComp,
                     const Box&       a_ecBox,
                     const FArrayBox& a_ccFAB,
                     const int        a_ccFABComp,
                     const Real       a_alpha,
                     const Real       a_beta)
{
    // Find the EC dirs.
    int ecDir0 = 0;
    for (; ecDir0 < SpaceDim; ++ecDir0) {
        if (a_ecBox.type(ecDir0) == IndexType::NODE) break;
    }
    CH_assert(ecDir0 < SpaceDim - 1);

    int ecDir1 = ecDir0 + 1;
    for (; ecDir1 < SpaceDim; ++ecDir1) {
        if (a_ecBox.type(ecDir1) == IndexType::NODE) break;
    }
    CH_assert(ecDir1 < SpaceDim);


    // Sanity checks
    CH_assert(a_ecFAB.interval().contains(a_ecFABComp));
    CH_assert(a_ccFAB.interval().contains(a_ccFABComp));
    CH_assert(a_ecFAB.box().contains(a_ecBox));
    CH_assert(a_ccFAB.box().contains(
        enclosedCells( grow(a_ecBox, BASISV(ecDir0) + BASISV(ecDir1)) )
    ));

    // Do it!
    FORT_ECMULTCC (
        CHF_FRA1(a_ecFAB, a_ecFABComp),
        CHF_CONST_FRA1(a_ccFAB, a_ccFABComp),
        CHF_CONST_REAL(a_alpha),
        CHF_CONST_REAL(a_beta),
        CHF_BOX(a_ecBox),
        CHF_CONST_INT(ecDir0),
        CHF_CONST_INT(ecDir1));
}


// -----------------------------------------------------------------------------
// Computes FC = FC * (alpha + beta*Av[CC])
// -----------------------------------------------------------------------------
void
FABAlgebra::FCmultCC(FArrayBox&       a_fcFAB,
                     const int        a_fcFABComp,
                     const Box&       a_fcBox,
                     const FArrayBox& a_ccFAB,
                     const int        a_ccFABComp,
                     const Real       a_alpha,
                     const Real       a_beta)
{
    // Find the FC dir.
    int fcDir = 0;
    for (; fcDir < SpaceDim; ++fcDir) {
        if (a_fcBox.type(fcDir) == IndexType::NODE) break;
    }
    CH_assert(fcDir != SpaceDim);

    // Sanity checks
    CH_assert(a_fcFAB.interval().contains(a_fcFABComp));
    CH_assert(a_ccFAB.interval().contains(a_ccFABComp));
    CH_assert(a_fcFAB.box().contains(a_fcBox));
    CH_assert(a_ccFAB.box().contains(
        enclosedCells(grow(a_fcBox, BASISV(fcDir)), fcDir)));

    // Do it!
    FORT_FCMULTCC (
        CHF_FRA1(a_fcFAB, a_fcFABComp),
        CHF_CONST_FRA1(a_ccFAB, a_ccFABComp),
        CHF_CONST_REAL(a_alpha),
        CHF_CONST_REAL(a_beta),
        CHF_BOX(a_fcBox),
        CHF_CONST_INT(fcDir));
}


// -----------------------------------------------------------------------------
// Computes CC1 = CC1 * (alpha + beta*CC2)
// -----------------------------------------------------------------------------
void
FABAlgebra::CCmultCC(FArrayBox&       a_cc1FAB,
                     const int        a_cc1FABComp,
                     const Box&       a_region,
                     const FArrayBox& a_cc2FAB,
                     const int        a_cc2FABComp,
                     const Real       a_alpha,
                     const Real       a_beta)
{
    // Sanity checks
    CH_assert(a_cc1FAB.interval().contains(a_cc1FABComp));
    CH_assert(a_cc2FAB.interval().contains(a_cc2FABComp));
    CH_assert(a_cc1FAB.box().contains(a_region));
    CH_assert(a_cc2FAB.box().contains(a_region));

    // Do it!
    FORT_CCMULTCC (
        CHF_FRA1(a_cc1FAB, a_cc1FABComp),
        CHF_CONST_FRA1(a_cc2FAB, a_cc2FABComp),
        CHF_CONST_REAL(a_alpha),
        CHF_CONST_REAL(a_beta),
        CHF_BOX(a_region));
}


// -----------------------------------------------------------------------------
// Computes FC = FC + scale * Av[CC].
// -----------------------------------------------------------------------------
void
FABAlgebra::FCaddCC(FArrayBox&       a_fcFAB,
                    const int        a_fcFABComp,
                    const Box&       a_fcBox,
                    const FArrayBox& a_ccFAB,
                    const int        a_ccFABComp,
                    const Real       a_scale)
{
    // Find the FC dir.
    int fcDir = 0;
    for (; fcDir < SpaceDim; ++fcDir) {
        if (a_fcBox.type(fcDir) == IndexType::NODE) break;
    }
    CH_assert(fcDir != SpaceDim);

    // Sanity checks
    CH_assert(a_fcFAB.interval().contains(a_fcFABComp));
    CH_assert(a_ccFAB.interval().contains(a_ccFABComp));
    CH_assert(a_fcFAB.box().contains(a_fcBox));
    CH_assert(a_ccFAB.box().contains(
        enclosedCells(grow(a_fcBox, BASISV(fcDir)), fcDir)));

    // Do it!
    FORT_FCADDCC (
        CHF_FRA1(a_fcFAB, a_fcFABComp),
        CHF_CONST_FRA1(a_ccFAB, a_ccFABComp),
        CHF_CONST_REAL(a_scale),
        CHF_BOX(a_fcBox),
        CHF_CONST_INT(fcDir));
}


// -----------------------------------------------------------------------------
// Computes CC = CC + scale * Av[FC].
// -----------------------------------------------------------------------------
void
FABAlgebra::CCaddFC(FArrayBox&       a_ccFAB,
                    const int        a_ccFABComp,
                    const Box&       a_ccBox,
                    const FArrayBox& a_fcFAB,
                    const int        a_fcFABComp,
                    const Real       a_scale)
{
    // Find the FC dir.
    int fcDir = 0;
    for (; fcDir < SpaceDim; ++fcDir) {
        if (a_fcFAB.box().type(fcDir) == IndexType::NODE) break;
    }
    CH_assert(fcDir != SpaceDim);

    // Sanity checks
    CH_assert(a_ccFAB.interval().contains(a_ccFABComp));
    CH_assert(a_fcFAB.interval().contains(a_fcFABComp));
    CH_assert(a_ccFAB.box().contains(a_ccBox));
    CH_assert(a_fcFAB.box().contains(surroundingNodes(a_ccBox, fcDir)));

    // Do it!
    FORT_CCADDFC (
        CHF_FRA1(a_ccFAB, a_ccFABComp),
        CHF_CONST_FRA1(a_fcFAB, a_fcFABComp),
        CHF_CONST_REAL(a_scale),
        CHF_BOX(a_ccBox),
        CHF_CONST_INT(fcDir));
}
