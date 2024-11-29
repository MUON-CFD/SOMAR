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
#include "TrilinearInterp.H"
#include "TrilinearInterpF_F.H"

// -----------------------------------------------------------------------------
// Simple trilinear interpolation utility.
//
// f_{ABCDEFGH}(u,v,w) = f_{A}*(1-u)*(1-v)*(1-w)
//                     + f_{B}*   u *(1-v)*(1-w)
//                     + f_{C}*(1-u)*(1-v)*   w
//                     + f_{D}*   u *(1-v)*   w
//                     + f_{E}*(1-u)*   v *(1-w)
//                     + f_{F}*   u *   v *(1-w)
//                     + f_{G}*(1-u)*   v *   w
//                     + f_{H}*   u *   v *   w
// where (u,v,w) are the coordinate in the middle of the box and
// the nodes are labeled as A, B, C, E, F, G, and H in this order:
//     Front    Back
//     C---D    G---H
//     |   |    |   |
//     A---B    E---F
// -----------------------------------------------------------------------------
void TrilinearInterp3D (FArrayBox&          a_fInterp,
                        const FArrayBox&    a_xInterp,
                        const FArrayBox&    a_yInterp,
                        const FArrayBox&    a_zInterp,
                        const Box&          a_interpBox,
                        const int           a_xdir,
                        const int           a_ydir,
                        const int           a_zdir,
                        const Vector<Real>& a_x,
                        const Vector<Real>& a_y,
                        const Vector<Real>& a_z,
                        const FArrayBox&    a_f)
{
#ifndef NDEBUG
    {
        // Check centerings.
        CH_assert(a_fInterp.box().type() == a_interpBox.type());
        CH_assert(a_xInterp.box().type() == a_interpBox.type());
        CH_assert(a_yInterp.box().type() == a_interpBox.type());
        CH_assert(a_zInterp.box().type() == a_interpBox.type());

        // Check FAB regions.
        CH_assert(a_fInterp.box().contains(a_interpBox));
        CH_assert(a_xInterp.box().contains(a_interpBox));
        CH_assert(a_yInterp.box().contains(a_interpBox));
        CH_assert(a_zInterp.box().contains(a_interpBox));

        CH_assert(a_f.box().size(a_xdir) == static_cast<int>(a_x.size()));
        CH_assert(a_f.box().size(a_ydir) == static_cast<int>(a_y.size()));
        CH_assert(a_f.box().size(a_zdir) == static_cast<int>(a_z.size()));

        // Check number of comps.
        CH_assert(a_xInterp.nComp() == 1);
        CH_assert(a_yInterp.nComp() == 1);
        CH_assert(a_zInterp.nComp() == 1);
        CH_assert(a_f      .nComp() == a_fInterp.nComp());

        // Check dirs
        CH_assert(0 <= a_xdir);
        CH_assert(a_xdir < SpaceDim);
        CH_assert(0 <= a_ydir);
        CH_assert(a_ydir < SpaceDim);
        CH_assert(0 <= a_zdir);
        CH_assert(a_zdir < SpaceDim);
        CH_assert(a_xdir != a_ydir);
        CH_assert(a_ydir != a_zdir);
    }
#endif

        // The fortran function requires the vector indices to
        // coincide with the constraint indices.
        IntVect shift = IntVect::Zero;
        shift[a_xdir] = a_f.box().smallEnd(a_xdir);
        shift[a_ydir] = a_f.box().smallEnd(a_ydir);
        shift[a_zdir] = a_f.box().smallEnd(a_zdir);

        // Interpolate
        FORT_TRILINEARINTERP3DF (
            CHF_FRA(a_fInterp),
            CHF_CONST_FRA1(a_xInterp,0),
            CHF_CONST_FRA1(a_yInterp,0),
            CHF_CONST_FRA1(a_zInterp,0),
            CHF_BOX(a_interpBox),
            CHF_CONST_INT(a_xdir),
            CHF_CONST_INT(a_ydir),
            CHF_CONST_INT(a_zdir),
            CHF_CONST_VR(a_x),
            CHF_CONST_VR(a_y),
            CHF_CONST_VR(a_z),
            CHF_CONST_FRA_SHIFT(a_f,shift));
}


#if CH_SPACEDIM == 3
#ifndef NDEBUG
// -----------------------------------------------------------------------------
// A test of TrilinearInterp3D. By refining inBox by 2, we should get an
// error drop of 4. This returns the inf-norm of the error.
// -----------------------------------------------------------------------------
#include <random>
#include "AnisotropicMeshRefine.H"
#include "Comm.H"
#include "IO.H"
// #include "Analysis.H"
#include "CartesianMap.H"
#include "BoxIterator.H"

Real
TrilinearInterp3D_Test(
    std::function<Real(Real, Real, Real)> a_func[[maybe_unused]],
    const IntVect&                        a_ref[[maybe_unused]])
{
#if 0
    // auto lambda = [](Real x, Real y, Real z) -> Real {return x*x*y*y*z*z;};

    RealVect inLoL(0.0, 0.0, -10.0);
    RealVect inHiL(1.0, 3.0,  10.0);
    IntVect  inLoIV(0, 0, 0);
    IntVect  inHiIV(9, 9, 9);
    Box inBox(inLoIV, inHiIV);
    inBox.refine(a_ref);
    RealVect inDx = (inHiL - inLoL) / RealVect(inBox.size());

    RealVect outLoL(0.0, 0.2, -9.0);
    RealVect outHiL(1.0, 2.8, 9.0);
    IntVect outLoIV(0, 0, 0);
    IntVect outHiIV(127, 127, 63);
    RealVect outDx = (outHiL - outLoL) / RealVect(outHiIV - outLoIV + IntVect::Unit);
    IntVect maxBoxSize(32, 32, 32);
    int blockFactor = 8;

    Box outDomBox(outLoIV, outHiIV);
    const bool isPeriodic[] = {true, false, false};
    ProblemDomain outDomain(outDomBox, isPeriodic);
    DisjointBoxLayout outGrids;
    {
        Vector<Box> vbox;
        AnisotropicMeshRefine::domainSplit(
            outDomain, vbox, maxBoxSize, blockFactor);
        outGrids.defineAndLoadBalance(vbox, nullptr, outDomain);
    }
    DataIterator dit = outGrids.dataIterator();
    GeoSourceInterface* geoSrcPtr = new CartesianMap;
    LevelGeometry levGeo(outDomain, outHiL - outLoL, nullptr, geoSrcPtr);
    levGeo.createMetricCache(outGrids);


    // Locate input points.
    std::mt19937 gen;
    gen.seed(1);
    std::uniform_real_distribution<Real> dis0(-0.25*inDx[0], +0.25*inDx[0]);
    std::uniform_real_distribution<Real> dis1(-0.25*inDx[1], +0.25*inDx[1]);
    std::uniform_real_distribution<Real> dis2(-0.25*inDx[2], +0.25*inDx[2]);

    Vector<Real> inX(inBox.size(0));
    for(int i = 0; i < inX.size(); ++i) {
        inX[i] = inLoL[0] + (Real(i) + 0.5) * inDx[0] + dis0(gen);
    }

    Vector<Real> inY(inBox.size(1));
    for(int j = 0; j < inY.size(); ++j) {
        inY[j] = inLoL[1] + (Real(j) + 0.5) * inDx[1] + dis1(gen);
    }

    Vector<Real> inZ(inBox.size(2));
    for(int k = 0; k < inZ.size(); ++k) {
        inZ[k] = inLoL[2] + (Real(k) + 0.5) * inDx[2] + dis2(gen);
    }

    // Create input data.
    FArrayBox inDataFAB(inBox, 1);
    Real x, y, z;
    for (BoxIterator bit(inBox); bit.ok(); ++bit) {
        const IntVect& cc = bit();
        x = inX[cc[0]];
        y = inY[cc[1]];
        z = inZ[cc[2]];

        inDataFAB(cc) = a_func(x, y, z);
    }

    // Locate the output points.
    LevelData<FArrayBox> outX(outGrids, 1);
    LevelData<FArrayBox> outY(outGrids, 1);
    LevelData<FArrayBox> outZ(outGrids, 1);
    for (dit.reset(); dit.ok(); ++dit) {
        const Box& valid = outGrids[dit];

        for (BoxIterator bit(valid); bit.ok(); ++bit) {
            const IntVect& cc = bit();
            outX[dit](cc) = outLoL[0] + (Real(cc[0]) + 0.5) * outDx[0];
            outY[dit](cc) = outLoL[1] + (Real(cc[1]) + 0.5) * outDx[1];
            outZ[dit](cc) = outLoL[2] + (Real(cc[2]) + 0.5) * outDx[2];
        }
    }

    // Interpolate
    LevelData<FArrayBox> outData(outGrids, 1);
    for (dit.reset(); dit.ok(); ++dit) {
        const Box& valid = outGrids[dit];

        TrilinearInterp3D(outData[dit],
                          outX[dit],
                          outY[dit],
                          outZ[dit],
                          valid,
                          0, // xdir
                          1, // ydir
                          SpaceDim-1, // zdir
                          inX,
                          inY,
                          inZ,
                          inDataFAB);
    }

    // Check results.
    LevelData<FArrayBox> soln(outGrids, 1);
    for (dit.reset(); dit.ok(); ++dit) {
        const Box& valid = outGrids[dit];
        Real x, y, z, f;

        for (BoxIterator bit(valid); bit.ok(); ++bit) {
            const IntVect& cc = bit();
            x = outX[dit](cc);
            y = outY[dit](cc);
            z = outZ[dit](cc);

            soln[dit](cc) = a_func(x, y, z);
        }
    }
    IO::writeHDF5("outData_000.hdf5", soln, levGeo);    // What to expect
    IO::writeHDF5("outData_001.hdf5", outData, levGeo); // What we have

    for (dit.reset(); dit.ok(); ++dit) {
        outData[dit].plus(soln[dit], -1.0);
    }
    IO::writeHDF5("outData_002.hdf5", outData, levGeo); // Error

    Real norm = Analysis::pNorm(outData, 0);
    IO::tout(0) << "inf-norm = " << norm << endl;

    delete geoSrcPtr;
    geoSrcPtr = nullptr;

    return norm;
#else
    return quietNAN;
#endif // 0
}
#endif // !NDEBUG
#endif // CH_SPACEDIM == 3
