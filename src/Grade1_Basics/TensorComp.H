#ifdef CH_LANG_CC
/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2017
 *    Jefferson (Philadelphia University + Thomas Jefferson University) and
 *    University of North Carolina at Chapel Hill
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
#endif //CH_LANG_CC
#ifndef ___TensorComp_H__INCLUDED___
#define ___TensorComp_H__INCLUDED___

#ifdef CH_LANG_FORT
#   define xxComp2D 0
#   define xyComp2D 1
#   define yxComp2D 2
#   define yyComp2D 3
#   define xxSymComp2D 0
#   define xySymComp2D 1
#   define yxSymComp2D 1
#   define yySymComp2D 2

#   define xxComp3D 0
#   define xyComp3D 1
#   define xzComp3D 2
#   define yxComp3D 3
#   define yyComp3D 4
#   define yzComp3D 5
#   define zxComp3D 6
#   define zyComp3D 7
#   define zzComp3D 8
#   define xxSymComp3D 0
#   define xySymComp3D 1
#   define xzSymComp3D 2
#   define yxSymComp3D 1
#   define yySymComp3D 3
#   define yzSymComp3D 4
#   define zxSymComp3D 2
#   define zySymComp3D 4
#   define zzSymComp3D 5

#   if CH_SPACEDIM == 2
#       define xxComp xxComp2D
#       define xyComp xyComp2D
#       define yxComp yxComp2D
#       define yyComp yyComp2D
#       define xxSymComp xxSymComp2D
#       define xySymComp xySymComp2D
#       define yxSymComp yxSymComp2D
#       define yySymComp yySymComp2D
#   elif CH_SPACEDIM == 3
#       define xxComp xxComp3D
#       define xyComp xyComp3D
#       define xzComp xzComp3D
#       define yxComp yxComp3D
#       define yyComp yyComp3D
#       define yzComp yzComp3D
#       define zxComp zxComp3D
#       define zyComp zyComp3D
#       define zzComp zzComp3D
#       define xxSymComp xxSymComp3D
#       define xySymComp xySymComp3D
#       define xzSymComp xzSymComp3D
#       define yxSymComp yxSymComp3D
#       define yySymComp yySymComp3D
#       define yzSymComp yzSymComp3D
#       define zxSymComp zxSymComp3D
#       define zySymComp zySymComp3D
#       define zzSymComp zzSymComp3D
#   else
#       error Bad CH_SPACEDIM
#   endif
#endif //CH_LANG_FORT


#ifdef CH_LANG_CC
#include "SPACE.H"
#include "LevelData.H"

class TensorComp
{
public:
    // Utility for converting CC tensor indices to a FAB component
    static inline int CC2D (int, int);
    static inline int symCC2D (int, int);
    static inline int symCC2D (int, int, int);

    static inline int CC3D (int, int);
    static inline int symCC3D (int, int);
    static inline int symCC3D (int, int, int);

    static inline int CC (int, int);
    static inline int symCC (int, int);
    static inline int symCC (int, int, int);

    // Computes v^i = Sum_j [M^{ij} * v^j].
    // This is set up so that is can work in-place on v if needed.
    // This will compute on overlap of boxes in vi, Mij, and vj.
    static void
    contractMatrixVector(LevelData<FArrayBox>&       a_vi,
                         const LevelData<FArrayBox>& a_Mij,
                         const LevelData<FArrayBox>& a_vj,
                         const bool                  a_jIsFastest = true);

    // Computes v^i = Sum_j [M^{ij} * v^j].
    // This is set up so that is can work in-place on v if needed.
    static void
    contractMatrixVector(FArrayBox&       a_viFAB,
                         const FArrayBox& a_MijFAB,
                         const FArrayBox& a_vjFAB,
                         const Box&       a_region,
                         const bool       a_jIsFastest = true);

    // Contract two matrices on only one index. This is set up to work
    // in-place if C = A or B. The second indices must be fastest.
    // Notation:
    //   C_ij = Sum_k [A_ik * B_kj]   if a_dofirstBIdx == true (default)
    //   C_ij = Sum_k [A_ik * B_jk]   if a_dofirstBIdx == false
    static void
    contractMatrixMatrix(FArrayBox&       a_CFAB,
                         const FArrayBox& a_AFAB,
                         const FArrayBox& a_BFAB,
                         const Box&       a_region,
                         const bool       a_doFirstBIdx = true);

    // Sends M_ij to M_ji.
    static void
    transpose(LevelData<FArrayBox>& a_M);

    // Sends M_ij to M_ji.
    static void
    transpose(FArrayBox& a_MFAB);

protected:
    static const int s_CC_2Comp_2D[2][2];
    static const int s_symCC_2Comp_2D[2][2];
    static const int s_symCC_3Comp_2D[2][2][2];

    static const int s_CC_2Comp_3D[3][3];
    static const int s_symCC_2Comp_3D[3][3];
    static const int s_symCC_3Comp_3D[3][3][3];
};



// -----------------------------------------------------------------------------
// Implementation of inline functions
// -----------------------------------------------------------------------------

// Utility for converting tensor indices to a FAB component.
// Mimics TensorCFInterp::gradIndex(i,graddir) which returns i*SpaceDim+graddir.
// mode. The assert is good enough to tell you where things went wrong.
int TensorComp::CC2D (int a_1, int a_2) {
    CH_assert((0 <= a_1) && (a_1 < 2));
    CH_assert((0 <= a_2) && (a_2 < 2));
    return s_CC_2Comp_2D[a_1][a_2];
}


// Utility for converting symmetric tensor indices to a FAB component.
// mode. The assert is good enough to tell you where things went wrong.
int TensorComp::symCC2D (int a_1, int a_2) {
    CH_assert((0 <= a_1) && (a_1 < 2));
    CH_assert((0 <= a_2) && (a_2 < 2));
    return s_symCC_2Comp_2D[a_1][a_2];
}


// Utility for converting tensor indices to a FAB component.
// Tensor is symmetric in the last two indices.
// mode. The assert is good enough to tell you where things went wrong.
int TensorComp::symCC2D (int a_1, int a_2, int a_3) {
    CH_assert((0 <= a_1) && (a_1 < 2));
    CH_assert((0 <= a_2) && (a_2 < 2));
    CH_assert((0 <= a_3) && (a_3 < 2));
    return s_symCC_3Comp_2D[a_1][a_2][a_3];
}


// Utility for converting tensor indices to a FAB component.
// Mimics TensorCFInterp::gradIndex(i,graddir) which returns i*SpaceDim+graddir.
// mode. The assert is good enough to tell you where things went wrong.
int TensorComp::CC3D (int a_1, int a_2) {
    CH_assert((0 <= a_1) && (a_1 < 3));
    CH_assert((0 <= a_2) && (a_2 < 3));
    return s_CC_2Comp_3D[a_1][a_2];
}


// Utility for converting symmetric tensor indices to a FAB component.
// mode. The assert is good enough to tell you where things went wrong.
int TensorComp::symCC3D (int a_1, int a_2) {
    CH_assert((0 <= a_1) && (a_1 < 3));
    CH_assert((0 <= a_2) && (a_2 < 3));
    return s_symCC_2Comp_3D[a_1][a_2];
}


// Utility for converting tensor indices to a FAB component.
// Tensor is symmetric in the last two indices.
// mode. The assert is good enough to tell you where things went wrong.
int TensorComp::symCC3D (int a_1, int a_2, int a_3) {
    CH_assert((0 <= a_1) && (a_1 < 3));
    CH_assert((0 <= a_2) && (a_2 < 3));
    CH_assert((0 <= a_3) && (a_3 < 3));
    return s_symCC_3Comp_3D[a_1][a_2][a_3];
}


// Utility for converting tensor indices to a FAB component.
// Mimics TensorCFInterp::gradIndex(i,graddir) which returns i*SpaceDim+graddir.
// mode. The assert is good enough to tell you where things went wrong.
int TensorComp::CC (int a_1, int a_2) {
#if CH_SPACEDIM == 2
    return CC2D(a_1, a_2);
#else
    return CC3D(a_1, a_2);
#endif
}


// Utility for converting symmetric tensor indices to a FAB component.
// mode. The assert is good enough to tell you where things went wrong.
int TensorComp::symCC (int a_1, int a_2) {
#if CH_SPACEDIM == 2
    return symCC2D(a_1, a_2);
#else
    return symCC3D(a_1, a_2);
#endif
}


// Utility for converting tensor indices to a FAB component.
// Tensor is symmetric in the last two indices.
// mode. The assert is good enough to tell you where things went wrong.
int TensorComp::symCC (int a_1, int a_2, int a_3) {
#if CH_SPACEDIM == 2
    return symCC2D(a_1, a_2, a_3);
#else
    return symCC3D(a_1, a_2, a_3);
#endif
}

#endif //CH_LANG_CC
#endif //!___TensorComp_H__INCLUDED___
