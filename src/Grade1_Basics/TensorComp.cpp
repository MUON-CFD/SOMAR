#include "TensorComp.H"
#include "TensorCompF_F.H"

// Define the constants used to look up tensor indices.

// 2D...
const int TensorComp::s_CC_2Comp_2D[2][2]
    = {{0,1},{2,3}};

const int TensorComp::s_symCC_2Comp_2D[2][2]
    = {{0,1},{1,2}};

const int TensorComp::s_symCC_3Comp_2D[2][2][2]
    = {{{0,1},{1,2}},{{3,4},{4,5}}};


// 3D...
const int TensorComp::s_CC_2Comp_3D[3][3]
    = {{0,1,2},{3,4,5},{6,7,8}};

const int TensorComp::s_symCC_2Comp_3D[3][3]
    = {{0,1,2},{1,3,4},{2,4,5}};

const int TensorComp::s_symCC_3Comp_3D[3][3][3]
    = { {{ 0, 1, 2},{ 1, 3, 4},{ 2, 4, 5}},
        {{ 6, 7, 8},{ 7, 9,10},{ 8,10,11}},
        {{12,13,14},{13,15,16},{14,16,17}} };


// -----------------------------------------------------------------------------
// Computes v^i = Sum_j [M^{ij} * v^j].
// This is set up so that is can work in-place on v if needed.
// -----------------------------------------------------------------------------
void
TensorComp::contractMatrixVector(LevelData<FArrayBox>&       a_vi,
                                 const LevelData<FArrayBox>& a_Mij,
                                 const LevelData<FArrayBox>& a_vj,
                                 const bool                  a_jIsFastest)
{
    CH_assert(a_vi.nComp() == SpaceDim);
    CH_assert(a_vj.nComp() == SpaceDim);
    CH_assert(a_Mij.nComp() == SpaceDim*SpaceDim);

    CH_assert(a_vi.getBoxes().compatible(a_vj.getBoxes()));
    CH_assert(a_vi.getBoxes().compatible(a_Mij.getBoxes()));

    const int jIsFastestInt = (a_jIsFastest? 1: 0);
    DataIterator dit = a_vi.dataIterator();

    for (dit.reset(); dit.ok(); ++dit) {
        const Box overlap = a_vi[dit].box() & a_vj[dit].box() & a_Mij[dit].box();
        FORT_TENSORCOMP_CONTRACTMATRIXVECTORCC(
            CHF_FRA(a_vi[dit]),
            CHF_FRA(a_Mij[dit]),
            CHF_FRA(a_vj[dit]),
            CHF_BOX(overlap),
            CHF_CONST_INT(jIsFastestInt));
    }
}


// -----------------------------------------------------------------------------
// Computes v^i = Sum_j [M^{ij} * v^j].
// This is set up so that is can work in-place on v if needed.
// -----------------------------------------------------------------------------
void
TensorComp::contractMatrixVector(FArrayBox&       a_viFAB,
                                 const FArrayBox& a_MijFAB,
                                 const FArrayBox& a_vjFAB,
                                 const Box&       a_region,
                                 const bool       a_jIsFastest)
{
    CH_assert(a_viFAB .nComp() == SpaceDim);
    CH_assert(a_vjFAB .nComp() == SpaceDim);
    CH_assert(a_MijFAB.nComp() == SpaceDim*SpaceDim);

    CH_assert(a_vjFAB .box().sameType(a_viFAB.box()));
    CH_assert(a_MijFAB.box().sameType(a_viFAB.box()));

    CH_assert(a_viFAB .box().contains(a_region));
    CH_assert(a_vjFAB .box().contains(a_region));
    CH_assert(a_MijFAB.box().contains(a_region));

    const int jIsFastestInt = (a_jIsFastest? 1: 0);
    FORT_TENSORCOMP_CONTRACTMATRIXVECTORCC(
        CHF_FRA(a_viFAB),
        CHF_FRA(a_MijFAB),
        CHF_FRA(a_vjFAB),
        CHF_BOX(a_region),
        CHF_CONST_INT(jIsFastestInt));
}


// -----------------------------------------------------------------------------
// Contract two matrices on only one index. This is set up to work
// in-place if C = A or B. The second indices must be fastest.
// Notation:
//   C_ij = Sum_k [A_ik * B_kj]   if a_dofirstBIdx == true (default)
//   C_ij = Sum_k [A_ik * B_jk]   if a_dofirstBIdx == false
// -----------------------------------------------------------------------------
void
TensorComp::contractMatrixMatrix(FArrayBox&       a_CFAB,
                                 const FArrayBox& a_AFAB,
                                 const FArrayBox& a_BFAB,
                                 const Box&       a_region,
                                 const bool       a_doFirstBIdx)
{
    CH_assert(a_CFAB.nComp() == SpaceDim*SpaceDim);
    CH_assert(a_AFAB.nComp() == SpaceDim*SpaceDim);
    CH_assert(a_BFAB.nComp() == SpaceDim*SpaceDim);

    CH_assert(a_AFAB.box().sameType(a_CFAB.box()));
    CH_assert(a_BFAB.box().sameType(a_CFAB.box()));

    CH_assert(a_CFAB.box().contains(a_region));
    CH_assert(a_AFAB.box().contains(a_region));
    CH_assert(a_BFAB.box().contains(a_region));

    const int mode = (a_doFirstBIdx? 0: 1);
    FORT_TENSORCOMP_CONTRACTMATRIXMATRIXCC(
        CHF_FRA(a_CFAB),
        CHF_FRA(a_AFAB),
        CHF_FRA(a_BFAB),
        CHF_BOX(a_region),
        CHF_CONST_INT(mode));
}


// -----------------------------------------------------------------------------
// Sends M_ij to M_ji.
// -----------------------------------------------------------------------------
void
TensorComp::transpose(LevelData<FArrayBox>& a_M)
{
    DataIterator dit = a_M.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        transpose(a_M[dit]);
    }
}


// -----------------------------------------------------------------------------
// Sends M_ij to M_ji.
// -----------------------------------------------------------------------------
void
TensorComp::transpose(FArrayBox& a_MFAB)
{
    CH_assert(a_MFAB.nComp() == SpaceDim*SpaceDim);

    const Box& region = a_MFAB.box();
    FORT_TENSORCOMP_TRANSPOSE(
        CHF_FRA(a_MFAB),
        CHF_BOX(region));
}

