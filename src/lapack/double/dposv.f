!> \brief <b> DPOSV computes the solution to system of linear equations A * X = B for PO matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DPOSV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dposv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dposv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dposv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DPOSV computes the solution to a real system of linear equations
!>    A * X = B,
!> where A is an N-by-N symmetric positive definite matrix and X and B
!> are N-by-NRHS matrices.
!>
!> The Cholesky decomposition is used to factor A as
!>    A = U**T* U,  if UPLO = 'U', or
!>    A = L * L**T,  if UPLO = 'L',
!> where U is an upper triangular matrix and L is a lower triangular
!> matrix.  The factored form of A is then used to solve the system of
!> equations A * X = B.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored;
!>          = 'L':  Lower triangle of A is stored.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of linear equations, i.e., the order of the
!>          matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrix B.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!>          N-by-N upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading N-by-N lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>
!>          On exit, if INFO = 0, the factor U or L from the Cholesky
!>          factorization A = U**T*U or A = L*L**T.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
!>          On entry, the N-by-NRHS right hand side matrix B.
!>          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, the leading minor of order i of A is not
!>                positive definite, so the factorization could not be
!>                completed, and the solution has not been computed.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date November 2011
!
!> \ingroup doublePOsolve
!
!  =====================================================================
      SUBROUTINE DPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
!
!  -- LAPACK driver routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!     ..
!
!  =====================================================================
!
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DPOTRF, DPOTRS, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPOSV ', -INFO )
         RETURN
      END IF
!
!     Compute the Cholesky factorization A = U**T*U or A = L*L**T.
!
      CALL DPOTRF( UPLO, N, A, LDA, INFO )
      IF( INFO.EQ.0 ) THEN
!
!        Solve the system A*X = B, overwriting B with X.
!
         CALL DPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
!
      END IF
      RETURN
!
!     End of DPOSV
!
      END
