#include "ScalarBC.H"
#include "ScalarBCF_F.H"

#include "ProblemContext.H"
#include "Subspace.H"

namespace ScalarBC {


// -----------------------------------------------------------------------------
// Constructor -- no background stratification
NormalGradient::NormalGradient()
: m_bcVal({ 0, 0 })
{
}

// Constructor -- with a background stratification
NormalGradient::NormalGradient(const FArrayBox&     a_bgStateFAB,
                               const LevelGeometry& a_levGeo)
{
    // Gather references, compute domain extents, etc.
    const IntVect ev     = BASISV(SpaceDim - 1);
    const Real    dZeta  = a_levGeo.getDXi(SpaceDim - 1);
    const Box&    domBox = a_levGeo.getDomainBox();
    const IntVect loIV   = domBox.smallEnd() * ev;
    const IntVect hiIV   = domBox.bigEnd() * ev;

    // Sanity checks.
    CH_assert(a_bgStateFAB.nComp() == 1);
    CH_assert(a_bgStateFAB.box().type() == IntVect::Zero);
    CH_assert(a_bgStateFAB.box().smallEnd() * (IntVect::Unit - ev) == IntVect::Zero);
    CH_assert(a_bgStateFAB.box().bigEnd() * (IntVect::Unit - ev) == IntVect::Zero);

    // Lower BC value = outward normal derivative of bgState
    //                = downward derivative of bgState
    //                = -d(bgState)/dZeta
    m_bcVal[0] = -(a_bgStateFAB(loIV) - a_bgStateFAB(loIV - ev)) / dZeta;

    // Upper BC value = outward normal derivative of bgState
    //                = upward derivative of bgState
    //                = d(bgState)/dZeta
    m_bcVal[1] = (a_bgStateFAB(hiIV + ev) - a_bgStateFAB(hiIV)) / dZeta;
}

// The BCFunction override.
void
NormalGradient::operator()(FArrayBox&            a_alpha,
                           FArrayBox&            a_beta,
                           FArrayBox&            a_bcFAB,
                           const FArrayBox&      /*a_stateFAB*/,
                           const FArrayBox&      /*a_xFAB*/,
                           const DataIndex&      /*a_di*/,
                           const int             a_bdryDir,
                           const Side::LoHiSide& a_side,
                           const Real            /*a_time*/,
                           const bool            a_homogBCs) const
{
    a_alpha.setVal(0.0);
    a_beta.setVal(1.0);
    if (!a_homogBCs) {
        if (a_bdryDir < SpaceDim - 1) {
            a_bcFAB.setVal(0.0);
        } else {
            a_bcFAB.setVal(m_bcVal[int(a_side)]);
        }
    }
}


// -----------------------------------------------------------------------------
// Constructor -- no background stratification
Tidal::Tidal()
{
    const ProblemContext* ctx = ProblemContext::getInstance();
    m_tidalOmega     = ctx->rhs.tidalOmega;
    m_tidalInitPhase = ctx->rhs.tidalInitPhase;
}

// The BCFunction override.
void
Tidal::operator()(FArrayBox&            a_alpha,
                  FArrayBox&            a_beta,
                  FArrayBox&            a_bcFAB,
                  const FArrayBox&      /*a_stateFAB*/,
                  const FArrayBox&      /*a_xFAB*/,
                  const DataIndex&      /*a_di*/,
                  const int             a_bdryDir,
                  const Side::LoHiSide& a_side,
                  const Real            a_time,
                  const bool            a_homogBCs) const
{
    bool isInflow;
    {
        const Real omega     = m_tidalOmega[a_bdryDir];
        const Real initPhase = m_tidalInitPhase[a_bdryDir];

        const Real phase         = omega * a_time + initPhase;
        const int  tidalFlowSign = (sin(phase) >= 0.0) ? (+1) : (-1);

        isInflow = false;
        isInflow |= (tidalFlowSign == +1 && a_side == Side::Lo);
        isInflow |= (tidalFlowSign == -1 && a_side == Side::Hi);
    }

    if (isInflow) {
        // Set equal to inflow value = bgState.
        a_alpha.setVal(1.0);
        a_beta.setVal(0.0);
        if (!a_homogBCs) a_bcFAB.setVal(0.0); // totalState = bgState. Override if needed.

    } else {
        // Outflow. Just set Grad[state].n = 0.
        a_alpha.setVal(0.0);
        a_beta.setVal(1.0);
        if (!a_homogBCs) a_bcFAB.setVal(0.0);
    }
}


// -----------------------------------------------------------------------------
// Constructor
BackgroundScalarWrapper::BackgroundScalarWrapper(
    const std::shared_ptr<BCTools::BCFunction>& a_baseBCPtr,
    const std::shared_ptr<FArrayBox>&           a_bgStateFABPtr,
    const std::shared_ptr<FArrayBox>&           a_dzFABPtr)
: m_baseBCPtr(a_baseBCPtr)
, m_bgStateFABPtr(a_bgStateFABPtr)
, m_dzFABPtr(a_dzFABPtr)
{
    CH_verify(m_baseBCPtr);
}

// The BCFunction override.
void
BackgroundScalarWrapper::operator()(FArrayBox&            a_alpha,
                                    FArrayBox&            a_beta,
                                    FArrayBox&            a_bcFAB,
                                    const FArrayBox&      a_stateFAB,
                                    const FArrayBox&      a_xFAB,
                                    const DataIndex&      a_di,
                                    const int             a_bdryDir,
                                    const Side::LoHiSide& a_side,
                                    const Real            a_time,
                                    const bool            a_homogBCs) const
{
    m_baseBCPtr->operator()(a_alpha,
                            a_beta,
                            a_bcFAB,
                            a_stateFAB,
                            a_xFAB,
                            a_di,
                            a_bdryDir,
                            a_side,
                            a_time,
                            a_homogBCs);

    if (m_bgStateFABPtr && !a_homogBCs) {
        const FArrayBox& phiBarFAB = *m_bgStateFABPtr;
        const FArrayBox& dzFAB     = *m_dzFABPtr;
        const int        isign     = sign(a_side);

        Box ghostBox = a_bcFAB.box();
        ghostBox.shiftHalf(a_bdryDir, isign);

        IntVect shiftIV = IntVect::Zero;
        if (a_side == Side::Lo) shiftIV[a_bdryDir] = 1;

        FORT_SCALARBC_SETBACKGROUNDRHS(
            CHF_FRA_SHIFT(a_bcFAB, shiftIV),
            CHF_CONST_FRA(phiBarFAB),
            CHF_CONST_FRA_SHIFT(a_alpha, shiftIV),
            CHF_CONST_FRA_SHIFT(a_beta, shiftIV),
            CHF_CONST_FRA1_SHIFT(dzFAB, 0, shiftIV),
            CHF_BOX(ghostBox),
            CHF_CONST_INT(a_bdryDir),
            CHF_CONST_INT(isign));
    }
}


};  // end namespace ScalarBC
