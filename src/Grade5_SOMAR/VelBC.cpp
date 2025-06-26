#include "VelBC.H"

#include "ProblemContext.H"
namespace VelBC
{

// For all of these functions, the BCs are set via
//   alpha * phi + beta * dphi/dn = BC.
// We need to tell the BC code what the alpha, beta, and BC values are at
// each point on the boundary.


// -----------------------------------------------------------------------------
void
NoSlip::operator()(FArrayBox& a_alpha,
                   FArrayBox& a_beta,
                   FArrayBox& a_bcFAB,
                   const FArrayBox& /*a_stateFAB*/,
                   const FArrayBox& /*a_xFAB*/,
                   const DataIndex& /*a_di*/,
                   const int /*a_bdryDir*/,
                   const Side::LoHiSide& /*a_side*/,
                   const Real /*a_time*/,
                   const bool a_homogBCs) const
{
    a_alpha.setVal(1.0);
    a_beta.setVal(0.0);
    if (!a_homogBCs) a_bcFAB.setVal(0.0);
}


// -----------------------------------------------------------------------------
void
FreeSlip::operator()(FArrayBox&       a_alpha,
                     FArrayBox&       a_beta,
                     FArrayBox&       a_bcFAB,
                     const FArrayBox& a_stateFAB,
                     const FArrayBox& /*a_xFAB*/,
                     const DataIndex& /*a_di*/,
                     const int a_bdryDir,
                     const Side::LoHiSide& /*a_side*/,
                     const Real /*a_time*/,
                     const bool a_homogBCs) const
{
    const int velComp = BCTools::getNodalDir(a_stateFAB.box());
    if (velComp == a_bdryDir) {
        // Normal component -- no flow through wall.
        a_alpha.setVal(1.0);
        a_beta.setVal(0.0);
        if (!a_homogBCs) a_bcFAB.setVal(0.0);
    } else {
        // Transverse component -- no stress normal to wall.
        a_alpha.setVal(0.0);
        a_beta.setVal(1.0);
        if (!a_homogBCs) a_bcFAB.setVal(0.0);
    }
}


// -----------------------------------------------------------------------------
Tidal::Tidal()
{
    const ProblemContext* ctx = ProblemContext::getInstance();
    m_tidalU0                 = ctx->rhs.tidalU0;
    m_tidalOmega              = ctx->rhs.tidalOmega;
    m_tidalInitPhase          = ctx->rhs.tidalInitPhase;
}

void
Tidal::operator()(FArrayBox&       a_alpha,
                  FArrayBox&       a_beta,
                  FArrayBox&       a_bcFAB,
                  const FArrayBox& a_stateFAB,
                  const FArrayBox& /*a_xFAB*/,
                  const DataIndex& /*a_di*/,
                  const int             a_bdryDir,
                  const Side::LoHiSide& a_side,
                  const Real            a_time,
                  const bool            a_homogBCs) const
{
    const int velComp = BCTools::getNodalDir(a_stateFAB.box());

    Real tidalU;
    {
        const Real omega     = m_tidalOmega[velComp];
        const Real initPhase = m_tidalInitPhase[velComp];
        const Real tidalU0   = m_tidalU0[velComp];

        tidalU = tidalU0 * sin(omega * a_time + initPhase);
    }

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
        // Set equal to inflow velocity.
        a_alpha.setVal(1.0);
        a_beta.setVal(0.0);
        if (!a_homogBCs) a_bcFAB.setVal(tidalU);

    } else {
        if (velComp == a_bdryDir) {
            // Set equal to outflow velocity in boundary-normal direction...
            a_alpha.setVal(1.0);
            a_beta.setVal(0.0);
            if (!a_homogBCs) a_bcFAB.setVal(tidalU);

        } else {
            // ...and no-stress in boundary-transverse directions.
            a_alpha.setVal(0.0);
            a_beta.setVal(1.0);
            if (!a_homogBCs) a_bcFAB.setVal(0.0);
        }
    }
}


// -----------------------------------------------------------------------------
Inflow::Inflow(const RealVect& a_vel)
: m_type(UNIFORM_VALUE)
, m_vel(a_vel)
{
}

Inflow::Inflow(std::function<RealVect(const RealVect&)> a_func)
: m_type(FUNCTION)
, m_func(a_func)
{
}

void
Inflow::operator()(FArrayBox&       a_alpha,
                   FArrayBox&       a_beta,
                   FArrayBox&       a_bcFAB,
                   const FArrayBox& a_stateFAB,
                   const FArrayBox& a_xFAB,
                   const DataIndex& /*a_di*/,
                   const int             /* a_bdryDir */,
                   const Side::LoHiSide& /* a_side */,
                   const Real            /* a_time */,
                   const bool            a_homogBCs) const
{
    a_alpha.setVal(1.0);
    a_beta.setVal(0.0);

    if (a_homogBCs) return;

    const int velComp = BCTools::getNodalDir(a_stateFAB.box());

    if (m_type == UNIFORM_VALUE) {
        a_bcFAB.setVal(m_vel[velComp]);
    } else {
        for (BoxIterator bit(a_bcFAB.box()); bit.ok(); ++bit) {
            const IntVect& iv = bit();
            const RealVect x(D_DECL(a_xFAB(iv,0), a_xFAB(iv,1), a_xFAB(iv,2)));
            a_bcFAB(iv) = m_func(x)[velComp];
        }
    }
}


// -----------------------------------------------------------------------------
void
Outflow::operator()(FArrayBox& a_alpha,
                    FArrayBox& a_beta,
                    FArrayBox& a_bcFAB,
                    const FArrayBox& /*a_stateFAB*/,
                    const FArrayBox& /*a_xFAB*/,
                    const DataIndex& /*a_di*/,
                    const int /*a_bdryDir*/,
                    const Side::LoHiSide& /*a_side*/,
                    const Real /*a_time*/,
                    const bool a_homogBCs) const
{
    a_alpha.setVal(0.0);
    a_beta.setVal(1.0);
    if (!a_homogBCs) a_bcFAB.setVal(0.0);
}


};  // end namespace VelBC
