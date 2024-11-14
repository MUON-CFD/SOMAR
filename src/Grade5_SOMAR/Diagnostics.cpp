#include "Diagnostics.H"
#include "SetValLevel.H"
#include "FABAlgebra.H"
#include "FiniteDiff.H"
#include "Convert.H"
#include "IB.H"
#include "StaggeredFluxLD.H"
#include "ProblemContext.H"
namespace Diagnostics{
/// Calculates a_N2
void N2(LevelData<FArrayBox> & a_N2,
            const LevelData<FArrayBox> & a_b,
            const LevelGeometry* const a_levGeoPtr)
{
    const DisjointBoxLayout& grids = a_levGeoPtr->getBoxes();
    DataIterator             dit   = grids.dataIterator();

    LevelData<FluxBox> gradB(grids, 1);
    FiniteDiff(*a_levGeoPtr).levelGradientMAC(gradB, a_b);

    setValLevel(a_N2, 0.0);
    for (dit.reset(); dit.ok(); ++dit) {
        const Box ccValid = grids[dit];
        FABAlgebra::CCaddFC(
            a_N2[dit], 0, ccValid, gradB[dit][SpaceDim - 1], 0, -1.0);
    }
}

void
S2(LevelData<FArrayBox>&      a_S2,
   const LevelData<FluxBox>&  a_vel,
   const LevelGeometry* const a_levGeoPtr)
{
    const DisjointBoxLayout&  grids   = a_levGeoPtr->getBoxes();
    DataIterator              dit     = grids.dataIterator();
    const LevelData<FluxBox>& cartVel = a_vel;

    LevelData<FArrayBox> V(grids, SpaceDim, IntVect::Unit);
    Convert::FacesToCells(V, cartVel);

    LevelData<FArrayBox> U(grids, 1, IntVect::Unit);
    V.copyTo(Interval(0, 0), U, Interval(0, 0));
    U.exchange();

    LevelData<FluxBox> gradU(grids, 1);
    FiniteDiff(*a_levGeoPtr).levelGradientMAC(gradU, U);

    setValLevel(a_S2, 0.0);
    for (dit.reset(); dit.ok(); ++dit) {
        const Box ccValid = grids[dit];
        FABAlgebra::CCaddFC(
            a_S2[dit], 0, ccValid, gradU[dit][SpaceDim - 1], 0, 1.);
        a_S2[dit] *= a_S2[dit];
        a_S2[dit] += 1e-5;
    }
}

void Ri(LevelData<FArrayBox>& a_Ri,
        const LevelData<FluxBox> & a_vel,
        const LevelData<FArrayBox> & a_b ,
        const LevelGeometry* const a_levGeoPtr)
{
    const DisjointBoxLayout& grids=a_Ri.getBoxes();

    LevelData<FArrayBox> s2(grids,1);
    S2(s2,a_vel, a_levGeoPtr);

    N2(a_Ri, a_b, a_levGeoPtr);



    DataIterator             dit   = grids.dataIterator();
    for (dit.reset(); dit.ok(); ++dit)
    {
        a_Ri[dit]/=s2[dit];
    }
}

// Pressure Force on obstacle
// -----------------------------------------------------------------------------

RealVect
calculatePressureForce(const LevelData<FArrayBox>& a_pressure, const LevelData<FArrayBox>& a_buoyancy, 
                            const RealVect& a_accel, const IB& a_IB) 
{
#if 0
    if (!(a_IB.m_CCStencilsDefined))
        MayDay::Error("To interpolate pressure on Lagrangian Markers EC Stencils must be defined");
    RealVect P = RealVect::Zero;
    const Real h = [&]() {
        ProblemContext const * __nowarn_unused ctx =
            ProblemContext::getInstance();
        const RealVect L  = ctx->base.L;
        const IntVect  nx = ctx->base.nx;
        Real h = L[0] / nx[0] + L[1] / nx[1];
        if (SpaceDim == 3) h += L[2] / nx[2];
        return h;
    }();
    Real Area = 0.;


    for (auto [_p, _b, _markers, _mask] : zip(a_pressure, a_buoyancy, a_IB.markers(), a_IB.mask())) {
        if (!_mask.go()) continue;
        
        for (const auto& m : _markers) {
            if (!m.interior) continue;
            // Compute Pressure at the probe marker
            Real Pe = 0.0;
            const auto& stencil = m.probeStencils[0]; //m.stencils[0];
            
            for (size_t e = 0; e < stencil.iv.size(); ++e) {
                const Real pe = _p(stencil.iv[e]);
                Pe += pe * stencil.coeff[e];
            } 

            // buoyancy
            Real b = 0.0;
            for (auto [_iv, _coeff] : zip(stencil.iv, stencil.coeff))
            {
                b += _b(_iv) * _coeff;
            }
            // extrapolated pressure at lagrangian marker (as in Vanella and Balaras, JCP 2009)
            const Real Pl = Pe + m.NVECT.dotProduct(-a_accel-b*BASISV(SpaceDim-1)) * h;

            P += m.NVECT * ( Pl * m.AREA); // note the sign of Pl.
            Area += m.AREA;
            
        } // m
        
    } // dit
    for (int dir=0; dir < SpaceDim; ++dir){
        Comm::reduce(P[dir], MPI_SUM);
    }
    
    return P;
    #else
    return RealVect::Zero;
    #endif
    }

// =============================================================================
// Stress force on obstacle
// -----------------------------------------------------------------------------
RealVect
calculateStressForce(const StaggeredFluxLD& a_S, const IB& a_IB)  
{
#if 0 
    if (!(a_IB.m_ECStencilsDefined))
        MayDay::Error("To interpolate Stress on Lagrangian Markers EC Stencils must be defined");
    RealVect Stress = RealVect::Zero;

    for (DataIterator dit(a_IB.grids()); dit.ok(); ++dit) {

        if (!a_IB.m_mask[dit].go()) continue;
        std::array<std::array<Real,SpaceDim>,SpaceDim> St;

        for (const auto& m : a_IB.markers()[dit]) {
            if (!m.interior) continue;
            // Compute stress at this Lagrangian marker.
            for (int j=0; j<SpaceDim; ++j){
                for (int i=0; i<SpaceDim; ++i){
                    St[j][i]=0.;
                    const FArrayBox& SFAB = a_S[j][i][dit];
                    if (i==j){
                        const auto&      stencil    = m.stencils[0];
                        for (size_t e = 0; e < stencil.iv.size(); ++e) {

                            const Real pe = SFAB(stencil.iv[e]);
                            St[j][i] += pe * stencil.coeff[e];

                        }  // diagonal
                    }else{

                        const IndexType& Idx  = SFAB.box().ixType();
                        const StencilData&     stencil = m.stencils[Idx.hash()];
                        for (size_t e = 0; e < stencil.iv.size(); ++e) {

                            const Real pe = SFAB(stencil.iv[e]);
                            St[j][i] += pe * stencil.coeff[e];

                        }  // off diagonal
                    }

                }
            }
            for (int j=0; j<SpaceDim; ++j){
                for (int i=0; i<SpaceDim; ++i){

                Stress[j] += St[j][i] * m.NVECT[i]  * m.AREA;
                }
            }
        } // m
    } // dit

    for (int dir=0; dir < SpaceDim; ++dir){
        Comm::reduce(Stress[dir], MPI_SUM);
    }
    return Stress;
#else
    return RealVect::Zero;
#endif

}


};
