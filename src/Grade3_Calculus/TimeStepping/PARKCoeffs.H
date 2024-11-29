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
#ifndef ___PARKCoeffs_H__INCLUDED___
#define ___PARKCoeffs_H__INCLUDED___

#include "SOMAR_Constants.H"

namespace RK {
    static constexpr Real Unsure   = 1.0;
    static constexpr Real Unstable = 0.0;
    // static constexpr Real Inf      = std::numeric_limits<Real>::max();
};  // namespace RK

// -----------------------------------------------------------------------------
// Forward Euler
struct ForwardEuler_Coeffs
{
    // Number of stages
    static constexpr int numStages = 1;

    // Time coefficients for each stage
    static constexpr Real c[numStages] = {0.} ;

    // Stage coefficients - row # = stage #.
    static constexpr Real aE[numStages][numStages] = {
  {0.}
};
    static constexpr Real aI[numStages][numStages] = {
  {0.}
};

    // Assembly coefficients
    static constexpr Real bE[numStages]  =
  {1.};
    static constexpr Real bI[numStages]  =
  {0.};

    // Embedded scheme
    static constexpr bool hasEmbeddedScheme = false;
    static constexpr Real bEhat[numStages]  =
  {quietNAN};
    static constexpr Real bIhat[numStages]  =
  {quietNAN};
    // Dense output coefficients
    static constexpr bool hasDenseOutput = false;
    static constexpr int numDenseCoeffs = 1;
    static constexpr Real bstar[numDenseCoeffs][numStages] = {
  {quietNAN}
};

    // Stability limits
    static constexpr Real ERKStabilityRe = 2.0;
    static constexpr Real ERKStabilityIm = RK::Unstable;
};


// -----------------------------------------------------------------------------
// Midpoint Method
struct Midpoint_Coeffs
{
    // 0.5 = Midpoint, 1.0 = Heun
    static constexpr Real alpha = 0.5;

    // Number of stages
    static constexpr int numStages = 2;

    // Time coefficients for each stage
    static constexpr Real c[numStages] = {0., Midpoint_Coeffs::alpha} ;

    // Stage coefficients - row # = stage #.
    static constexpr Real aE[numStages][numStages] = {
  {               0., 0. },
  {Midpoint_Coeffs::alpha, 0. }
};
    static constexpr Real aI[numStages][numStages] = {
  {0., 0.},
  {0., 0.}
};

    // Assembly coefficients
    static constexpr Real bE[numStages]  =
  {1. - 1. / (2. * Midpoint_Coeffs::alpha), 1. / (2. * Midpoint_Coeffs::alpha)};
    static constexpr Real bI[numStages] =
  { 0., 0.};

    // Embedded scheme
    static constexpr bool hasEmbeddedScheme = false;
    static constexpr Real bEhat[numStages] =
  {quietNAN, quietNAN};
    static constexpr Real bIhat[numStages] =
  {quietNAN, quietNAN};

    // Dense output coefficients
    static constexpr bool hasDenseOutput = false;
    static constexpr int numDenseCoeffs = 1;
    static constexpr Real bstar[numDenseCoeffs][numStages] = {
  {quietNAN, quietNAN}
};

    // Stability limits
    static constexpr Real ERKStabilityRe = RK::Unsure;
    static constexpr Real ERKStabilityIm = RK::Unsure;
};


// -----------------------------------------------------------------------------
// 2nd order accurate
// 2 and 3 register versions. 3 register requires slightly fewer flops.
// L-stable implicit part
// Has embedded scheme
// SSP -- well suited for hyperbolic systems
// Stage order 1
struct IMEXRKCB2_Coeffs
{
    // Number of stages
    static constexpr int numStages = 3;

    // Time coefficients for each stage
    static constexpr Real c[numStages] = {0.,  2./5.,  1.};

    // Stage coefficients - row # = stage #.
    static constexpr Real aE[numStages][numStages] = {
  {   0.,    0.,    0.},
  {2./5.,    0.,    0.},
  {   0.,    1.,    0.}
};
    static constexpr Real aI[numStages][numStages] = {
  {   0.,    0.,    0.},
  {   0., 2./5.,    0.},
  {   0., 5./6., 1./6.}
};

    // Assembly coefficients
    static constexpr Real bE[numStages]=
  {   0., 5./6., 1./6.};
    static constexpr Real bI[numStages]=
  {   0., 5./6., 1./6.};

    // Embedded scheme
    static constexpr bool hasEmbeddedScheme = true;
    static constexpr Real bEhat[numStages] =
  {   0., 4./5., 1./5.};
    static constexpr Real bIhat[numStages] =
  {   0., 4./5., 1./5.};

    // Dense output coefficients
    static constexpr bool hasDenseOutput = false;
    static constexpr int numDenseCoeffs = 1;
    static constexpr Real bstar[numDenseCoeffs][numStages] = {
  {quietNAN}
};

    // Stability limits
    static constexpr Real ERKStabilityRe = 5.81;
    static constexpr Real ERKStabilityIm = RK::Unsure;
};


// -----------------------------------------------------------------------------
// 3rd order accurate
// 2 and 3 register versions. 3 register requires slightly fewer flops.
// Strongly A-stable implicit part
// ESDIRK -- same implicit operator on all implicit solves.
// Stage order 1
struct IMEXRKCB3b_Coeffs
{
    // Number of stages
    static constexpr int numStages = 4;

    // Time coefficients for each stage
    static constexpr Real c[numStages]= {0., 0.788675134594812882254574390250978727823800875635063438009,  0.211324865405187117745425609749021272176199124364936561990,  0.788675134594812882254574390250978727823800875635063438009};

    // Stage coefficients - row # = stage #.
    static constexpr Real aE[numStages][numStages] = {
  {   0.,    0.,    0.,    0.},
  {   0.788675134594812882254574390250978727823800875635063438009,    0.,    0.,    0.},
  {   0.,    0.211324865405187117745425609749021272176199124364936561990,    0.,    0.},
  {   0.,    0.,    0.788675134594812882254574390250978727823800875635063438009,    0.}
};
    static constexpr Real aI[numStages][numStages] = {
  {   0.,    0.,    0.,    0.},
  {   0.,    0.788675134594812882254574390250978727823800875635063438009,    0.,    0.},
  {   0.,    -0.57735026918962576450914878050195745564760175127012687601,    0.788675134594812882254574390250978727823800875635063438009,    0.},
  {   0.,    0.,    0.,    0.788675134594812882254574390250978727823800875635063438009}
};

    // Assembly coefficients
    static constexpr Real bE[numStages] =
  {   0.,  0., 0.5, 0.5};
    static constexpr Real bI[numStages]  =
  {   0.,  0., 0.5, 0.5};

    // Dense output coefficients
    static constexpr bool hasDenseOutput = false;
    static constexpr int numDenseCoeffs = 1;
    static constexpr Real bstar[numDenseCoeffs][numStages] = {
  {quietNAN}
};

    // Stability limits
    static constexpr Real ERKStabilityRe = 2.21;
    static constexpr Real ERKStabilityIm = 2.6; // Roughly.
};


// -----------------------------------------------------------------------------
// 3rd order accurate
// 2 and 3 register versions. 3 register requires slightly fewer flops.
// L-stable implicit part
// Has embedded scheme
// SSP with c = 0.7027915 -- well suited for hyperbolic systems
// Stage order 1
struct IMEXRKCB3c_Coeffs
{
    // Number of stages
    static constexpr int numStages = 4;

    // Time coefficients for each stage
    static constexpr Real c[numStages] = {0., 3375509829940./4525919076317.,  272778623835./1039454778728.,  1.};

    // Stage coefficients - row # = stage #.
    static constexpr Real aE[numStages][numStages] = {
  {   0.,    0.,    0.,    0.},
  {   3375509829940./4525919076317.,    0.,    0.,    0.},
  {   0.,    272778623835./1039454778728.,    0.,    0.},
  {   0.,    673488652607./2334033219546.,    1660544566939./2334033219546.,    0.}
};
    static constexpr Real aI[numStages][numStages]= {
  {   0.,    0.,    0.,    0.},
  {   0.,    3375509829940./4525919076317.,    0.,    0.},
  {   0.,    -11712383888607531889907./32694570495602105556248.,    566138307881./912153721139.,    0.},
  {   0.,    673488652607./2334033219546.,    493801219040./853653026979.,    184814777513./1389668723319.}
};

    // Assembly coefficients
    static constexpr Real bE[numStages] =
  {   0., 673488652607./2334033219546., 493801219040./853653026979., 184814777513./1389668723319.};
    static constexpr Real bI[numStages] =
  {   0., 673488652607./2334033219546., 493801219040./853653026979., 184814777513./1389668723319.};

    // Embedded scheme
    static constexpr bool hasEmbeddedScheme = true;
    static constexpr Real bEhat[numStages] =
  { 449556814708./1155810555193., 0., 210901428686./ 1400818478499., 480175564215./1042748212601. };
    static constexpr Real bIhat[numStages]  =
  { 0., 366319659506./1093160237145., 270096253287./480244073137., 104228367309./1017021570740. };

    // Dense output coefficients
    static constexpr bool hasDenseOutput = false;
    static constexpr int numDenseCoeffs = 1;
    static constexpr Real bstar[numDenseCoeffs][numStages] = {
  {quietNAN}
};

    // Stability limits
    static constexpr Real ERKStabilityRe = 6.0;
    static constexpr Real ERKStabilityIm = 2.75; // Roughly.
};


// -----------------------------------------------------------------------------
// 3rd order accurate
// 2 and 3 register versions. 3 register requires slightly fewer flops.
// L-stable implicit part
// Stage order 1
// Maximized accuracy of explicit part.
struct IMEXRKCB3e_Coeffs
{
    // Number of stages
    static constexpr int numStages = 4;

    // Time coefficients for each stage
    static constexpr Real c[numStages] = {0., 1./3.,  1.,  1.};

    // Stage coefficients - row # = stage #.
    static constexpr Real aE[numStages][numStages] = {
  {   0.,     0.,      0.,      0.},
  {   1./3.,  0.,      0.,      0.},
  {   0.,     1.,      0.,      0.},
  {   0.,     3./4.,   1./4.,   0.}
};
    static constexpr Real aI[numStages][numStages] = {
  {   0.,     0.,      0.,      0.},
  {   0.,     1./3.,   0.,      0.},
  {   0.,     1./2.,   1./2.,   0.},
  {   0.,     3./4.,   -1./4.,  1./2.}
};

    // Assembly coefficients
    static constexpr Real bE[numStages] =
  {   0., 3./4., -1./4., 1./2.};
    static constexpr Real bI[numStages] =
  {   0., 3./4., -1./4., 1./2.};

    // Dense output coefficients
    static constexpr bool hasDenseOutput = false;
    static constexpr int numDenseCoeffs = 1;
    static constexpr Real bstar[numDenseCoeffs][numStages] = {
  {quietNAN}
};

    // Stability limits
    static constexpr Real ERKStabilityRe = 2.79;
    static constexpr Real ERKStabilityIm = 2.8; // Roughly.
};


// -----------------------------------------------------------------------------
// 3rd order accurate
// 3 and 4 register versions. 4 register requires significantly fewer flops.
// L-stable implicit part
// Has embedded scheme
// Stage order 2 -- well suited for stiff ODEs
struct IMEXRKCB3f_Coeffs
{
    // Number of stages
    static constexpr int numStages = 4;

    // Time coefficients for each stage
    static constexpr Real c[numStages]= {0., 49./50.,  1./25.,  1.};

    // Stage coefficients - row # = stage #.
    static constexpr Real aE[numStages][numStages] = {
  {   0.,    0.,    0.,    0.},
  {   49./50.,    0.,    0.,    0.},
  {   13244205847./647648310246.,    13419997131./686433909488.,    0.,    0.},
  {   -2179897048956./603118880443.,    231677526244./1085522130027.,    3007879347537./683461566472.,    0.}
};
    static constexpr Real aI[numStages][numStages] = {
  {   0.,    0.,    0.,    0.},
  {   49./100.,    49./100.,    0.,    0.},
  {   -785157464198./1093480182337.,    -30736234873./978681420651.,    983779726483./1246172347126.,    0.},
  {   -2179897048956./603118880443.,    99189146040./891495457793.,    6064140186914./1415701440113.,    146791865627./668377518349.}
};

    // Assembly coefficients
    static constexpr Real bE[numStages] =
  {   -2179897048956./603118880443., 99189146040./891495457793., 6064140186914./1415701440113., 146791865627./668377518349.};
    static constexpr Real bI[numStages] =
  {   -2179897048956./603118880443., 99189146040./891495457793., 6064140186914./1415701440113., 146791865627./668377518349.};

    // Embedded scheme
    static constexpr bool hasEmbeddedScheme = true;
    static constexpr Real bEhat[numStages]=
  { 0., 0., 25./48., 23./48. };
    static constexpr Real bIhat[numStages]=
  { 0., 337712514207./759004992869., 311412265155./608745789881., 52826596233./1214539205236. };

    // Dense output coefficients
    static constexpr bool hasDenseOutput = false;
    static constexpr int numDenseCoeffs = 1;
    static constexpr Real bstar[numDenseCoeffs][numStages] = {
  {quietNAN}
};

    // Stability limits
    static constexpr Real ERKStabilityRe = 6.0;
    static constexpr Real ERKStabilityIm = 2.75; // Roughly.
};


// -----------------------------------------------------------------------------
// 2nd order accurate
// 2 and 3 register versions. 3 register requires slightly fewer flops.
// A-stable implicit part
// Stage order 1
struct CN_RKW3_Coeffs
{
    // Number of stages
    static constexpr int numStages = 4;

    // Time coefficients for each stage
    static constexpr Real c[numStages]= { 0.,  8./15.,  2./3.,  1.};

    // Stage coefficients - row # = stage #.
    static constexpr Real aE[numStages][numStages] = {
  {    0.,     0.,    0.,  0.},
  {8./15.,     0.,    0.,  0.},
  { 1./4., 5./12.,    0.,  0.},
  { 1./4.,     0., 3./4.,  0.}
};
    static constexpr Real aI[numStages][numStages]  = {
  {    0.,     0.,     0.,    0.},
  {4./15., 4./15.,     0.,    0.},
  {4./15.,  1./3., 1./15.,    0.},
  {4./15.,  1./3., 7./30., 1./6.}
};

    // Assembly coefficients
    static constexpr Real bE[numStages] =
  {1./4.,  0.,  3./4.,  0.};
    static constexpr Real bI[numStages]=
  {4./15.,  1./3.,  7./30.,  1./6.};

    // Embedded scheme
    static constexpr bool hasEmbeddedScheme = false;
    static constexpr Real bEhat[numStages] =
  {quietNAN, quietNAN, quietNAN, quietNAN};
    static constexpr Real bIhat[numStages] =
  {quietNAN, quietNAN, quietNAN, quietNAN};

    // Dense ouput coefficients
    static constexpr bool hasDenseOutput = true;
    static constexpr int numDenseCoeffs = 2;
    static constexpr Real bstar[numDenseCoeffs][numStages] = {
  {    1.,    0.,    0.,    0.},
  {-3./4.,    0., 3./4.,    0.},
};
    static constexpr Real bIstar[numDenseCoeffs][numStages]= {
  {      1.,      0.,      0.,     0.},
  {-11./15.,   1./3.,  7./30.,  1./6.}
};

    // Stability limits
    static constexpr Real ERKStabilityRe = 2.51; // Roughly
    static constexpr Real ERKStabilityIm = 1.73205080757;
};


// -----------------------------------------------------------------------------
// Classical RK4
struct RK4_Coeffs
{
    // Number of stages
    static constexpr int numStages = 4;

    // Time coefficients for each stage
    static constexpr Real c[numStages]= { 0.,  0.5,  0.5,  1.};

    // Stage coefficients - row # = stage #.
    static constexpr Real aE[numStages][numStages] = {
  {  0.,  0.,  0.,  0.},
  { 0.5,  0.,  0.,  0.},
  {  0., 0.5,  0.,  0.},
  {  0.,  0.,  1.,  0.}
};
    static constexpr Real aI[numStages][numStages]= {
  {    0.,     0.,     0.,    0.},
  {    0.,     0.,     0.,    0.},
  {    0.,     0.,     0.,    0.},
  {    0.,     0.,     0.,    0.}
};

    // Assembly coefficients
    static constexpr Real bE[numStages] =
  {1./6., 1./3.,  1./3., 1./6.};
    static constexpr Real bI[numStages] =
  {   0.,    0.,     0.,    0.};

    // Embedded scheme
    static constexpr bool hasEmbeddedScheme = false;
    static constexpr Real bEhat[numStages] =
  {quietNAN, quietNAN, quietNAN, quietNAN};
    static constexpr Real bIhat[numStages]  =
  {quietNAN, quietNAN, quietNAN, quietNAN};

    // Dense output coefficients
    static constexpr bool hasDenseOutput = false;
    static constexpr int numDenseCoeffs = 1;
    static constexpr Real bstar[numDenseCoeffs][numStages]= {
  {quietNAN, quietNAN, quietNAN, quietNAN}
};

    // Stability limits
    static constexpr Real ERKStabilityRe = 2.7; // Roughly
    static constexpr Real ERKStabilityIm = 2.82;
};


// -----------------------------------------------------------------------------
struct ARK3_2_4L_2_SA_Coeffs
{
    // Number of stages
    static constexpr int numStages = 4;

    // Time coefficients for each stage
    static constexpr Real c[numStages] = { 0., 1767732205903./2027836641118., 3./5., 1. };

    // Stage coefficients - row # = stage #.
    static constexpr Real aE[numStages][numStages] = {
  {0., 0., 0., 0.},
  {1767732205903./2027836641118., 0., 0., 0.},
  {5535828885825./10492691773637., 788022342437./10882634858940., 0., 0.}, // !
  {6485989280629./16251701735622., -4246266847089./9704473918619., 10755448449292./10357097424841., 0.}
};
    static constexpr Real aI[numStages][numStages] = {
  {0., 0., 0., 0.},
  {1767732205903./4055673282236., 1767732205903./4055673282236., 0., 0.},
  {2746238789719./10658868560708., -640167445237/6845629431997., 1767732205903./4055673282236., 0.},
  {1471266399579./7840856788654., -4482444167858./7529755066697, 11266239266428./11593286722821, 1767732205903./4055673282236.}
};

    // Assembly coefficients
    static constexpr Real bE[numStages] =
  {1471266399579./7840856788654., -4482444167858./7529755066697., 11266239266428./11593286722821., 1767732205903./4055673282236.};
    static constexpr Real bI[numStages] =
  {1471266399579./7840856788654., -4482444167858./7529755066697., 11266239266428./11593286722821., 1767732205903./4055673282236.};

    // Embedded scheme
    static constexpr bool hasEmbeddedScheme = true;
    static constexpr Real bEhat[numStages]=
  {2756255671327./12835298489170., -10771552573575./22201958757719., 9247589265047./10645013368117., 2193209047091./5459859503100.};
    static constexpr Real bIhat[numStages] =
  {2756255671327./12835298489170., -10771552573575./22201958757719., 9247589265047./10645013368117., 2193209047091./5459859503100.};

    // Coefficients for dense ouput -- 3rd-order interpolation
    static constexpr bool hasDenseOutput = true;
    static constexpr int numDenseCoeffs = 2;
    static constexpr Real bstar[numDenseCoeffs][numStages]= {
  {4655552711362./22874653954995., -18682724506714./9892148508045., 34259539580243./13192909600954., 584795268549./6622622206610.},
  {-215264564351./13552729205753., 17870216137069./13817060693119., -28141676662227./17317692491321., 2508943948391./7218656332882.},
};

    // Stability limits
    static constexpr Real ERKStabilityRe = 3.66;
    static constexpr Real ERKStabilityIm = 2.48;
};


// -----------------------------------------------------------------------------
struct ARK4_3_6L_2_SA_Coeffs
{
    // Number of stages
    static constexpr int numStages = 6;

    // Time coefficients for each stage
    static constexpr Real c[numStages]= { 0., 0.5, 0.332, 0.62, 0.85, 1. };

    // Stage coefficients - row # = stage #.
    static constexpr Real aE[numStages][numStages]= {
  {0., 0., 0., 0., 0., 0.},
  {0.5, 0., 0., 0., 0., 0.},
  {0.221776, 0.110224, 0., 0., 0., 0.},
  {-0.04884659515311857, -0.17772065232640102, 0.8465672474795197, 0., 0., 0.},
  {-0.15541685842491548, -0.3567050098221991, 1.0587258798684427, 0.30339598837867193, 0., 0.},
  { 0.2014243506726763, 0.008742057842904185, 0.15993995707168115, 0.4038290605220775, 0.22606457389066084, 0.}
};
    static constexpr Real aI[numStages][numStages] = {
  {0., 0., 0., 0., 0., 0.},
  {0.25, 0.25, 0., 0., 0., 0.},
  {0.137776, -0.055776, 0.25, 0., 0., 0.},
  {0.14463686602698217, -0.22393190761334475, 0.4492950415863626, 0.25, 0., 0.},
  {0.09825878328356477, -0.5915442428196704, 0.8101210538282996, 0.283164405707806, 0.25, 0.},
  {0.15791629516167136, 0., 0.18675894052400077, 0.6805652953093346, -0.27524053099500667, 0.25}
};

    // Assembly coefficients
    static constexpr Real bE[numStages] =
  {0.15791629516167136, 0., 0.18675894052400077, 0.6805652953093346, -0.27524053099500667, 0.25};
    static constexpr Real bI[numStages]=
  {0.15791629516167136, 0., 0.18675894052400077, 0.6805652953093346, -0.27524053099500667, 0.25};

    // Embedded scheme
    static constexpr bool hasEmbeddedScheme = true;
    static constexpr Real bEhat[numStages]=
  {4586570599./29645900160., 0., 178811875./945068544., 814220225./1159782912., -3700637./11593932., 61727./225920.};
    static constexpr Real bIhat[numStages]=
  {4586570599./29645900160., 0., 178811875./945068544., 814220225./1159782912., -3700637./11593932., 61727./225920.};

    // Coefficients for dense ouput -- 4th-order interpolation
    static constexpr bool hasDenseOutput = true;
    static constexpr int numDenseCoeffs = 3;
    static constexpr Real bstar[numDenseCoeffs][numStages]= {
  {0.961753400252887, 0., 0.787405595186356, -2.74544192086633, 3.70351728061223, -1.70723435518514},
  {-1.76418754019038, 0., -0.774504669155511, 9.64023584441292, -12.544886411271, 5.44334277620397},
  {0.960350435099165, 0., 0.173858014493155, -6.21422862823726, 8.56612859966376, -3.48610842101883}
};

    // Stability limits
    static constexpr Real ERKStabilityRe = 4.23;
    static constexpr Real ERKStabilityIm = 4.00;
};



#endif //!___PARKCoeffs_H__INCLUDED___
