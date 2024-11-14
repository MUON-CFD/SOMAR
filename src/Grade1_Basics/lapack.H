/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2014 University of North Carolina at Chapel Hill
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
#ifndef __lapack_H__INCLUDED__
#define __lapack_H__INCLUDED__

namespace lapack {
extern "C" {
    // Function declarations for some LAPACK routines...

    // DLAMCH determines double precision machine parameters.
    double dlamch_ (char* CMACH);

    // DGBEQU computes row and column scalings intended to equilibrate an
    // M-by-N band matrix A and reduce its condition number.  R returns the
    // row scale factors and C the column scale factors, chosen to try to
    // make the largest element in each row and column of the matrix B with
    // elements B(i,j)=R(i)*A(i,j)*C(j) have absolute value 1.
    void dgbequ_(int* M, int* N, int* KL, int* KU, double* AB, int* LDAB,
                 double* R, double* C, double* ROWCND, double* COLCND, double* AMAX, int* INFO);

    // DSYGV computes all the eigenvalues, and optionally, the eigenvectors
    // of a real generalized symmetric-definite eigenproblem, of the form
    // A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.
    // Here A and B are assumed to be symmetric and B is also
    // positive definite.
    void dsygv_(int* ITYPE, char* JOBZ, char* UPLO, int* N,
                double* A, int* LDA, double* B, int* LDB,
                double* W, double* WORK, int* LWORK, int* INFO);

    // DSYGVD computes all the eigenvalues, and optionally, the eigenvectors
    // of a real generalized symmetric-definite eigenproblem, of the form
    // A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and
    // B are assumed to be symmetric and B is also positive definite.
    // If eigenvectors are desired, it uses a divide and conquer algorithm.

    // The divide and conquer algorithm makes very mild assumptions about
    // floating point arithmetic. It will work on machines with a guard
    // digit in add/subtract, or on those binary machines without guard
    // digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
    // Cray-2. It could conceivably fail on hexadecimal or decimal machines
    // without guard digits, but we know of none.
    void dsygvd_(int* ITYPE, char* JOBZ, char* UPLO, int* N,
                 double* A, int* LDA, double* B, int* LDB,
                 double* W, double* WORK, int* LWORK, int* IWORK, int* LIWORK, int* INFO);

    // DSYGVX computes selected eigenvalues, and optionally, eigenvectors
    // of a real generalized symmetric-definite eigenproblem, of the form
    // A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A
    // and B are assumed to be symmetric and B is also positive definite.
    // Eigenvalues and eigenvectors can be selected by specifying either a
    // range of values or a range of indices for the desired eigenvalues.
    void dsygvx_(int* ITYPE, char* JOBZ, char* UPLO, int* N, double* A, int* LDA, double* B, int* LDB,
                 double* VL, double* VU, int* IL, int* IU, double* ABSTOL, int* M, double* W, double* Z,
                 int* LDZ, double* WORK, int* LWORK, int* IWORK, int* IFAIL, int* INFO);

    // DSPGV computes all the eigenvalues and, optionally, the eigenvectors
    // of a real generalized symmetric-definite eigenproblem, of the form
    // A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.
    // Here A and B are assumed to be symmetric, stored in packed format,
    // and B is also positive definite.
    void dspgv_(int* ITYPE, char* JOBZ, char* UPLO, int* N, double* AP, double* BP, double* W, double* Z,
                int* LDZ, double* WORK, int* INFO);

    // DSPGVD computes all the eigenvalues, and optionally, the eigenvectors
    // of a real generalized symmetric-definite eigenproblem, of the form
    // A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and
    // B are assumed to be symmetric, stored in packed format, and B is also
    // positive definite.
    // If eigenvectors are desired, it uses a divide and conquer algorithm.

    // The divide and conquer algorithm makes very mild assumptions about
    // floating point arithmetic. It will work on machines with a guard
    // digit in add/subtract, or on those binary machines without guard
    // digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
    // Cray-2. It could conceivably fail on hexadecimal or decimal machines
    // without guard digits, but we know of none.
    void dspgvd_(int* ITYPE, char* JOBZ, char* UPLO, int* N, double* AP, double* BP, double* W, double* Z,
                 int* LDZ, double* WORK, int* LWORK, int* IWORK, int* LIWORK, int* INFO);

    // DSPGVX computes selected eigenvalues, and optionally, eigenvectors
    // of a real generalized symmetric-definite eigenproblem, of the form
    // A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A
    // and B are assumed to be symmetric, stored in packed storage, and B
    // is also positive definite.  Eigenvalues and eigenvectors can be
    // selected by specifying either a range of values or a range of indices
    // for the desired eigenvalues.
    void dspgvx_(int* ITYPE, char* JOBZ, char* RANGE, char* UPLO, int* N, double* AP, double* BP, double* VL,
                 double* VU, int* IL, int* IU, double* ABSTOL, int* M, double* W, double* Z, int* LDZ,
                 double* WORK, int* IWORK, int* IFAIL, int* INFO);

    // DSBGV computes all the eigenvalues, and optionally, the eigenvectors
    // of a real generalized symmetric-definite banded eigenproblem, of
    // the form A*x=(lambda)*B*x. Here A and B are assumed to be symmetric
    // and banded, and B is also positive definite.
    void dsbgv_(char* JOBZ, char* UPLO, int* N, int* KA, int* KB, double* AB, int* LDAB, double* BB,
                int* LDBB, double* W, double* Z, int* LDZ, double* WORK, int* INFO);

    // DSBGVD computes all the eigenvalues, and optionally, the eigenvectors
    // of a real generalized symmetric-definite banded eigenproblem, of the
    // form A*x=(lambda)*B*x.  Here A and B are assumed to be symmetric and
    // banded, and B is also positive definite.  If eigenvectors are
    // desired, it uses a divide and conquer algorithm.

    // The divide and conquer algorithm makes very mild assumptions about
    // floating point arithmetic. It will work on machines with a guard
    // digit in add/subtract, or on those binary machines without guard
    // digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
    // Cray-2. It could conceivably fail on hexadecimal or decimal machines
    // without guard digits, but we know of none.
    void dsbgvd_(char* JOBZ, char* UPLO, int* N, int* KA, int* KB, double* AB, int* LDAB,
                 double* BB, int* LDBB, double* W, double* Z, int* LDZ,
                 double* WORK, int* LWORK, int* IWORK, int* LIWORK, int* INFO);

    // DSBGVX computes selected eigenvalues, and optionally, eigenvectors
    // of a real generalized symmetric-definite banded eigenproblem, of
    // the form A*x=(lambda)*B*x.  Here A and B are assumed to be symmetric
    // and banded, and B is also positive definite.  Eigenvalues and
    // eigenvectors can be selected by specifying either all eigenvalues,
    // a range of values or a range of indices for the desired eigenvalues.
    void dsbgvx_ (char* JOBZ, char* RANGE, char* UPLO, int* N, int* KA, int* KB, double* AB, int* LDAB, double* BB,
                  int* LDBB, double* Q, int* LDQ, double* VL, double* VU, int* IL, int* IU, double* ABSTOL, int* M,
                  double* W, double* Z, int* LDZ, double* WORK, int* IWORK, int* IFAIL, int* INFO);

    // DSYSV computes the solution to a real system of linear equations
    //   A * X = B,
    // where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
    // matrices.


    // -------------------------------------------------------------------------
    /// \name Linear solve. Symmetric matrices
    /// \{

    // The diagonal pivoting method is used to factor A as
    //   A = U * D * U**T,  if UPLO = 'U', or
    //   A = L * D * L**T,  if UPLO = 'L',
    // where U (or L) is a product of permutation and unit upper (lower)
    // triangular matrices, and D is symmetric and block diagonal with
    // 1-by-1 and 2-by-2 diagonal blocks.  The factored form of A is then
    // used to solve the system of equations A * X = B.
    void
    dsysv_(char* UPLO, int* N, int* NRHS, double* A, int* LDA, int* IPIV,
           double* B, int* LDB, double* WORK, int* LWORK, int* INFO);

    // DSYSVX uses the diagonal pivoting factorization to compute the
    // solution to a real system of linear equations A * X = B,
    // where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
    // matrices.

    // Error bounds on the solution and a condition estimate are also
    // provided.
    void
    dsysvx_(char* FACT, char* UPLO, int* N, int* NRHS, double* A, int* LDA,
            double* AF, int* LDAF, int* IPIV, double* B, int* LDB,
            double* X, int* LDX, double* RCOND, double* FERR, double* BERR,
            double* WORK, int* LWORK, int* IWORK, int* INFO);

    // ERROR: This function is missing from my version of Lapack!
    // // DSYSVXX uses the diagonal pivoting factorization to compute the
    // // solution to a double precision system of linear equations A * X = B, where A
    // // is an N-by-N symmetric matrix and X and B are N-by-NRHS matrices.

    // // If requested, both normwise and maximum componentwise error bounds
    // // are returned. DSYSVXX will return a solution with a tiny
    // // guaranteed error (O(eps) where eps is the working machine
    // // precision) unless the matrix is very ill-conditioned, in which
    // // case a warning is returned. Relevant condition numbers also are
    // // calculated and returned.

    // // DSYSVXX accepts user-provided factorizations and equilibration
    // // factors; see the definitions of the FACT and EQUED options.
    // // Solving with refinement and using a factorization from a previous
    // // DSYSVXX call will also produce a solution with either O(eps)
    // // errors or warnings, but we cannot make that claim for general
    // // user-provided factorizations and equilibration factors if they
    // // differ from what DSYSVXX would itself produce.
    // void
    // dsysvxx_(char* FACT, char* UPLO, int* N, int* NRHS, double* A, int* LDA,
    //          double* AF, int* LDAF, int* IPIV, char* EQUED, double* S, double* B, int* LDB,
    //          double* X, int* LDX, double* RCOND, double* RPVGRW,
    //          double* BERR, int* N_ERR_BNDS, double* ERR_BNDS_NORM, double* ERR_BNDS_COMP,
    //          int* NPARAMS, double* PARAMS, double* WORK, int* IWORK, int* INFO);
    /// \}


    // -------------------------------------------------------------------------
    /// \name Linear least squares. General matrices.
    /// \{

    /// @brief Solves overdetermined or underdetermined systems for GE matrices.
    /// @param TRANS [in   ] = 'N': the linear system involves A;
    ///                      = 'T': the linear system involves A**T.
    /// @param M     [in   ] The number of rows of the matrix A.  M >= 0.
    /// @param N     [in   ] The number of columns of the matrices A and B. N >= 0.
    /// @param NRHS  [in   ] The number of right hand sides, i.e., the number of
    ///                      columns of the matrices B and X. NRHS >=0.
    /// @param A     [inout] On entry, the M-by-N matrix A.
    ///                      On exit,
    ///                          if M >= N, A is overwritten by details of its QR
    ///                                  factorization as returned by DGEQRF;
    ///                          if M <  N, A is overwritten by details of its LQ
    ///                                  factorization as returned by DGELQF.
    /// @param LDA   [in   ] The leading dimension of the array A. LDA >= max(1,M).
    /// @param B     [inout] On entry, the matrix B of right hand side vectors, stored
    ///                      columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS
    ///                      if TRANS = 'T'.
    ///                      On exit, if INFO = 0, B is overwritten by the solution
    ///                      vectors, stored columnwise:
    ///                      if TRANS = 'N' and m >= n, rows 1 to n of B contain the least
    ///                      squares solution vectors; the residual sum of squares for the
    ///                      solution in each column is given by the sum of squares of
    ///                      elements N+1 to M in that column;
    ///                      if TRANS = 'N' and m < n, rows 1 to N of B contain the
    ///                      minimum norm solution vectors;
    ///                      if TRANS = 'T' and m >= n, rows 1 to M of B contain the
    ///                      minimum norm solution vectors;
    ///                      if TRANS = 'T' and m < n, rows 1 to M of B contain the
    ///                      least squares solution vectors; the residual sum of squares
    ///                      for the solution in each column is given by the sum of
    ///                      squares of elements M+1 to N in that column.
    /// @param LDB   [in   ] The leading dimension of the array B. LDB >= MAX(1,M,N).
    /// @param WORK  [out  ] WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)).
    ///                      On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
    /// @param LWORK [in   ] LWORK is INTEGER. The dimension of the array WORK.
    ///                      LWORK >= max( 1, MN + max( MN, NRHS ) ).
    ///                      For optimal performance,
    ///                      LWORK >= max( 1, MN + max( MN, NRHS )*NB ).
    ///                      where MN = min(M,N) and NB is the optimum block size.
    ///                      If LWORK = -1, then a workspace query is assumed; the routine
    ///                      only calculates the optimal size of the WORK array, returns
    ///                      this value as the first entry of the WORK array, and no error
    ///                      message related to LWORK is issued by XERBLA.
    /// @param INFO  [out  ] = 0:  successful exit
    ///                      < 0:  if INFO = -i, the i-th argument had an illegal value
    ///                      > 0:  if INFO =  i, the i-th diagonal element of the
    ///                              triangular factor of A is zero, so that A does not have
    ///                              full rank; the least squares solution could not be
    ///                              computed.
    /// @details
    /// DGELS solves overdetermined or underdetermined real linear systems
    ///  involving an M-by-N matrix A, or its transpose, using a QR or LQ
    ///  factorization of A.  It is assumed that A has full rank.
    ///
    ///  The following options are provided:
    ///
    ///  1. If TRANS = 'N' and m >= n:  find the least squares solution of
    ///     an overdetermined system, i.e., solve the least squares problem
    ///                  minimize || B - A*X ||.
    ///
    ///  2. If TRANS = 'N' and m < n:  find the minimum norm solution of
    ///     an underdetermined system A * X = B.
    ///
    ///  3. If TRANS = 'T' and m >= n:  find the minimum norm solution of
    ///     an underdetermined system A**T * X = B.
    ///
    ///  4. If TRANS = 'T' and m < n:  find the least squares solution of
    ///     an overdetermined system, i.e., solve the least squares problem
    ///                  minimize || B - A**T * X ||.
    ///
    ///  Several right hand side vectors b and solution vectors x can be
    ///  handled in a single call; they are stored as the columns of the
    ///  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
    ///  matrix X.
    void
    dgels_(char* TRANS, int* M, int* N, int* NRHS,
           double* A, int* LDA, double* B, int* LDB,
           double* WORK, int* LWORK, int* INFO);


    void
    dgelss_(int* M, int* N, int* NRHS,
            double* A, int* LDA, double* B, int* LDB,
            double* S, double* RCOND, int* RANK,
            double* WORK, int* LWORK, int* INFO);



    /// @brief Solves overdetermined or underdetermined systems for OTHER matrices.
    /// @param M     [in   ] The number of rows of the matrix A.  M >= 0.
    /// @param N     [in   ] The number of columns of the matrices A and B. N >= 0.
    /// @param P     [in   ] The number of rows of the matrix B. 0 <= P <= N <= M+P.
    /// @param A     [inout] On entry, the M-by-N matrix A.
    ///                      On exit, the elements on and above the diagonal of the array
    ///                      contain the min(M,N)-by-N upper trapezoidal matrix T.
    /// @param LDA   [in   ] The leading dimension of the array A. LDA >= max(1,M).
    /// @param B     [inout] On entry, the P-by-N matrix B.
    ///                      On exit, the upper triangle of the subarray B(1:P,N-P+1:N)
    ///                      contains the P-by-P upper triangular matrix R.
    /// @param LDB   [in   ] The leading dimension of the array B. LDB >= max(1,P).
    /// @param C     [inout] On entry, C contains the right hand side vector for the
    ///                      least squares part of the LSE problem.
    ///                      On exit, the residual sum of squares for the solution
    ///                      is given by the sum of squares of elements N-P+1 to M of
    ///                      vector C.
    /// @param D     [inout] On entry, D contains the right hand side vector for the
    ///                      constrained equation. On exit, D is destroyed.
    /// @param X     [out  ] On exit, X is the solution of the LSE problem.
    /// @param WORK  [out  ] On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
    /// @param LWORK [in   ] The dimension of the array WORK. LWORK >= max(1,M+N+P).
    ///                      For optimum performance LWORK >= P+min(M,N)+max(M,N)*NB,
    ///                      where NB is an upper bound for the optimal blocksizes for
    ///                      DGEQRF, SGERQF, DORMQR and SORMRQ.
    ///                      If LWORK = -1, then a workspace query is assumed; the routine
    ///                      only calculates the optimal size of the WORK array, returns
    ///                      this value as the first entry of the WORK array, and no error
    ///                      message related to LWORK is issued by XERBLA.
    /// @param INFO  [out  ] = 0:  successful exit.
    ///                      < 0:  if INFO = -i, the i-th argument had an illegal value.
    ///                      = 1:  the upper triangular factor R associated with B in the
    ///                              generalized RQ factorization of the pair (B, A) is
    ///                              singular, so that rank(B) < P; the least squares
    ///                              solution could not be computed.
    ///                      = 2:  the (N-P) by (N-P) part of the upper trapezoidal factor
    ///                              T associated with A in the generalized RQ factorization
    ///                              of the pair (B, A) is singular, so that
    ///                              rank( (A) ) < N; the least squares solution could not
    ///                                  ( (B) )
    ///                              be computed.
    /// @details
    /// DGGLSE solves the linear equality-constrained least squares (LSE)
    /// problem:
    ///
    ///          minimize || c - A*x ||_2   subject to   B*x = d
    ///
    ///  where A is an M-by-N matrix, B is a P-by-N matrix, c is a given
    ///  M-vector, and d is a given P-vector. It is assumed that
    ///  P <= N <= M+P, and
    ///
    ///           rank(B) = P and  rank( (A) ) = N.
    ///                                ( (B) )
    ///
    ///  These conditions ensure that the LSE problem has a unique solution,
    ///  which is obtained using a generalized RQ factorization of the
    ///  matrices (B, A) given by
    ///
    ///     B = (0 R)*Q,   A = Z*T*Q.
    void
    dgglse_(int* M, int* N, int* P, double* A, int* LDA, double* B, int* LDB,
            double* C, double* D, double* X,
            double* WORK, int* LWORK, int* INFO);
    /// \}


    // -------------------------------------------------------------------------
    /// \name Eigensystems.
    /// \{

    /// @brief Computes the eigenvalues and, optionally, the left and/or right
    /// eigenvectors for SY matrices.
    /// @param JOBZ  [in   ] = 'N':  Compute eigenvalues only;
    ///                      = 'V':  Compute eigenvalues and eigenvectors.
    /// @param UPLO  [in   ] = 'U':  Upper triangle of A is stored;
    ///                      = 'L':  Lower triangle of A is stored.
    /// @param N     [in   ] The order of the matrix A.  N >= 0.
    /// @param A     [inout] A is DOUBLE PRECISION array, dimension (LDA, N)
    ///                      On entry, the symmetric matrix A.  If UPLO = 'U', the
    ///                      leading N-by-N upper triangular part of A contains the
    ///                      upper triangular part of the matrix A.  If UPLO = 'L',
    ///                      the leading N-by-N lower triangular part of A contains
    ///                      the lower triangular part of the matrix A.
    ///                      On exit, if JOBZ = 'V', then if INFO = 0, A contains the
    ///                      orthonormal eigenvectors of the matrix A.
    ///                      If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
    ///                      or the upper triangle (if UPLO='U') of A, including the
    ///                      diagonal, is destroyed.
    /// @param LDA   [in   ] The leading dimension of the array A.  LDA >= max(1,N).
    /// @param W     [out  ] W is DOUBLE PRECISION array, dimension (N)
    ///                      If INFO = 0, the eigenvalues in ascending order.
    /// @param WORK  [out  ] WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
    ///                      On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
    /// @param LWORK [in   ] The length of the array WORK.  LWORK >= max(1,3*N-1).
    ///                      For optimal efficiency, LWORK >= (NB+2)*N,
    ///                      where NB is the blocksize for DSYTRD returned by ILAENV.
    ///                      If LWORK = -1, then a workspace query is assumed; the routine
    ///                      only calculates the optimal size of the WORK array, returns
    ///                      this value as the first entry of the WORK array, and no error
    ///                      message related to LWORK is issued by XERBLA.
    /// @param INFO  [out  ] = 0:  successful exit
    ///                      < 0:  if INFO = -i, the i-th argument had an illegal value
    ///                      > 0:  if INFO = i, the algorithm failed to converge; i
    ///                            off-diagonal elements of an intermediate tridiagonal
    ///                            form did not converge to zero.
    void
    dsyev_(char* JOBZ, char* UPLO, int* N, double* A, int* LDA,
           double* W, double* WORK, int* LWORK, int* INFO);

    /// @brief Computes the eigenvalues and, optionally, the left and/or right
    /// eigenvectors for GE matrices
    /// @param JOBVL  [in   ] = 'N':  do not compute the left generalized eigenvectors;
    ///                       = 'V':  compute the left generalized eigenvectors.
    /// @param JOBVR  [in   ] = 'N':  do not compute the right generalized eigenvectors;
    ///                       = 'V':  compute the right generalized eigenvectors.
    /// @param N      [in   ] The order of the matrices A, B, VL, and VR.  N >= 0.
    /// @param A      [inout] A is DOUBLE PRECISION array, dimension (LDA, N)
    ///                       On entry, the matrix A in the pair (A,B).
    ///                       On exit, A has been overwritten.
    /// @param LDA    [in   ] The leading dimension of A.  LDA >= max(1,N).
    /// @param B      [inout] B is DOUBLE PRECISION array, dimension (LDB, N)
    ///                       On entry, the matrix B in the pair (A,B).
    ///                       On exit, B has been overwritten.
    /// @param LDB    [in   ] The leading dimension of B.  LDB >= max(1,N).
    /// @param ALPHAR [out  ] ALPHAR is DOUBLE PRECISION array, dimension (N)
    /// @param ALPHAI [out  ] ALPHAI is DOUBLE PRECISION array, dimension (N)
    /// @param BETA   [out  ] BETA is DOUBLE PRECISION array, dimension (N)
    ///                       On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
    ///                       be the generalized eigenvalues.  If ALPHAI(j) is zero, then
    ///                       the j-th eigenvalue is real; if positive, then the j-th and
    ///                       (j+1)-st eigenvalues are a complex conjugate pair, with
    ///                       ALPHAI(j+1) negative.
    ///                       Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
    ///                       may easily over- or underflow, and BETA(j) may even be zero.
    ///                       Thus, the user should avoid naively computing the ratio
    ///                       alpha/beta.  However, ALPHAR and ALPHAI will be always less
    ///                       than and usually comparable with norm(A) in magnitude, and
    ///                       BETA always less than and usually comparable with norm(B).
    /// @param VL     [out  ] VL is DOUBLE PRECISION array, dimension (LDVL,N)
    ///                       If JOBVL = 'V', the left eigenvectors u(j) are stored one
    ///                       after another in the columns of VL, in the same order as
    ///                       their eigenvalues. If the j-th eigenvalue is real, then
    ///                       u(j) = VL(:,j), the j-th column of VL. If the j-th and
    ///                       (j+1)-th eigenvalues form a complex conjugate pair, then
    ///                       u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1).
    ///                       Each eigenvector is scaled so the largest component has
    ///                       abs(real part)+abs(imag. part)=1.
    ///                       Not referenced if JOBVL = 'N'.
    /// @param LDVL   [in   ] The leading dimension of the matrix VL. LDVL >= 1, and
    ///                       if JOBVL = 'V', LDVL >= N.
    /// @param VR     [out  ] VR is DOUBLE PRECISION array, dimension (LDVR,N)
    ///                       If JOBVR = 'V', the right eigenvectors v(j) are stored one
    ///                       after another in the columns of VR, in the same order as
    ///                       their eigenvalues. If the j-th eigenvalue is real, then
    ///                       v(j) = VR(:,j), the j-th column of VR. If the j-th and
    ///                       (j+1)-th eigenvalues form a complex conjugate pair, then
    ///                       v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1).
    ///                       Each eigenvector is scaled so the largest component has
    ///                       abs(real part)+abs(imag. part)=1.
    ///                       Not referenced if JOBVR = 'N'.
    /// @param LDVR   [in   ] The leading dimension of the matrix VR. LDVR >= 1, and
    ///                       if JOBVR = 'V', LDVR >= N.
    /// @param WORK   [out  ] WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
    ///                       On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
    /// @param LWORK  [in   ] The dimension of the array WORK.  LWORK >= max(1,8*N).
    ///                       For good performance, LWORK must generally be larger.
    ///                       If LWORK = -1, then a workspace query is assumed; the routine
    ///                       only calculates the optimal size of the WORK array, returns
    ///                       this value as the first entry of the WORK array, and no error
    ///                       message related to LWORK is issued by XERBLA.
    /// @param INFO   [out  ] = 0:  successful exit
    ///                       < 0:  if INFO = -i, the i-th argument had an illegal value.
    ///                       = 1,...,N:
    ///                           The QZ iteration failed.  No eigenvectors have been
    ///                           calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)
    ///                           should be correct for j=INFO+1,...,N.
    ///                       > N:  =N+1: other than QZ iteration failed in DHGEQZ.
    ///                           =N+2: error return from DTGEVC.
    void
    dggev_(char* JOBVL, char* JOBVR,
           int* N, double* A, int* LDA, double* B, int* LDB,
           double* ALPHAR, double* ALPHAI, double* BETA,
           double* VL, int* LDVL, double* VR, int* LDVR,
           double* WORK, int* LWORK, int* INFO);
    /// \}



}; // extern "C"
}; // namespace lapack

#endif //!__lapack_H__INCLUDED__
