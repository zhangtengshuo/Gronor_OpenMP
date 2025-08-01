!> \brief \b ZHEGST
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZHEGST + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhegst.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhegst.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhegst.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHEGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, ITYPE, LDA, LDB, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZHEGST reduces a complex Hermitian-definite generalized
!> eigenproblem to standard form.
!>
!> If ITYPE = 1, the problem is A*x = lambda*B*x,
!> and A is overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H)
!>
!> If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
!> B*A*x = lambda*x, and A is overwritten by U*A*U**H or L**H*A*L.
!>
!> B must have been previously factorized as U**H*U or L*L**H by ZPOTRF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ITYPE
!> \verbatim
!>          ITYPE is INTEGER
!>          = 1: compute inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H);
!>          = 2 or 3: compute U*A*U**H or L**H*A*L.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored and B is factored as
!>                  U**H*U;
!>          = 'L':  Lower triangle of A is stored and B is factored as
!>                  L*L**H.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices A and B.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
!>          N-by-N upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading N-by-N lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>
!>          On exit, if INFO = 0, the transformed matrix, stored in the
!>          same format as A.
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
!>          B is COMPLEX*16 array, dimension (LDB,N)
!>          The triangular factor from the Cholesky factorization of B,
!>          as returned by ZPOTRF.
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
!> \date September 2012
!
!> \ingroup complex16HEcomputational
!
!  =====================================================================
      SUBROUTINE ZHEGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
!
!  -- LAPACK computational routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, ITYPE, LDA, LDB, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
      COMPLEX*16         CONE, HALF
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ),                           &
     &                   HALF = ( 0.5D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            K, KB, NB
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZHEGS2, ZHEMM, ZHER2K, ZTRMM, ZTRSM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( ITYPE.LT.1 .OR. ITYPE.GT.3 ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHEGST', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 )                                                              &
     &   RETURN
!
!     Determine the block size for this environment.
!
      NB = ILAENV( 1, 'ZHEGST', UPLO, N, -1, -1, -1 )
!
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
!
!        Use unblocked code
!
         CALL ZHEGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
      ELSE
!
!        Use blocked code
!
         IF( ITYPE.EQ.1 ) THEN
            IF( UPPER ) THEN
!
!              Compute inv(U**H)*A*inv(U)
!
               DO 10 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
!
!                 Update the upper triangle of A(k:n,k:n)
!
                  CALL ZHEGS2( ITYPE, UPLO, KB, A( K, K ), LDA,                 &
     &                         B( K, K ), LDB, INFO )
                  IF( K+KB.LE.N ) THEN
                     CALL ZTRSM( 'Left', UPLO, 'Conjugate transpose',           &
     &                           'Non-unit', KB, N-K-KB+1, CONE,                &
     &                           B( K, K ), LDB, A( K, K+KB ), LDA )
                     CALL ZHEMM( 'Left', UPLO, KB, N-K-KB+1, -HALF,             &
     &                           A( K, K ), LDA, B( K, K+KB ), LDB,             &
     &                           CONE, A( K, K+KB ), LDA )
                     CALL ZHER2K( UPLO, 'Conjugate transpose', N-K-KB+1,        &
     &                            KB, -CONE, A( K, K+KB ), LDA,                 &
     &                            B( K, K+KB ), LDB, ONE,                       &
     &                            A( K+KB, K+KB ), LDA )
                     CALL ZHEMM( 'Left', UPLO, KB, N-K-KB+1, -HALF,             &
     &                           A( K, K ), LDA, B( K, K+KB ), LDB,             &
     &                           CONE, A( K, K+KB ), LDA )
                     CALL ZTRSM( 'Right', UPLO, 'No transpose',                 &
     &                           'Non-unit', KB, N-K-KB+1, CONE,                &
     &                           B( K+KB, K+KB ), LDB, A( K, K+KB ),            &
     &                           LDA )
                  END IF
   10          CONTINUE
            ELSE
!
!              Compute inv(L)*A*inv(L**H)
!
               DO 20 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
!
!                 Update the lower triangle of A(k:n,k:n)
!
                  CALL ZHEGS2( ITYPE, UPLO, KB, A( K, K ), LDA,                 &
     &                         B( K, K ), LDB, INFO )
                  IF( K+KB.LE.N ) THEN
                     CALL ZTRSM( 'Right', UPLO, 'Conjugate transpose',          &
     &                           'Non-unit', N-K-KB+1, KB, CONE,                &
     &                           B( K, K ), LDB, A( K+KB, K ), LDA )
                     CALL ZHEMM( 'Right', UPLO, N-K-KB+1, KB, -HALF,            &
     &                           A( K, K ), LDA, B( K+KB, K ), LDB,             &
     &                           CONE, A( K+KB, K ), LDA )
                     CALL ZHER2K( UPLO, 'No transpose', N-K-KB+1, KB,           &
     &                            -CONE, A( K+KB, K ), LDA,                     &
     &                            B( K+KB, K ), LDB, ONE,                       &
     &                            A( K+KB, K+KB ), LDA )
                     CALL ZHEMM( 'Right', UPLO, N-K-KB+1, KB, -HALF,            &
     &                           A( K, K ), LDA, B( K+KB, K ), LDB,             &
     &                           CONE, A( K+KB, K ), LDA )
                     CALL ZTRSM( 'Left', UPLO, 'No transpose',                  &
     &                           'Non-unit', N-K-KB+1, KB, CONE,                &
     &                           B( K+KB, K+KB ), LDB, A( K+KB, K ),            &
     &                           LDA )
                  END IF
   20          CONTINUE
            END IF
         ELSE
            IF( UPPER ) THEN
!
!              Compute U*A*U**H
!
               DO 30 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
!
!                 Update the upper triangle of A(1:k+kb-1,1:k+kb-1)
!
                  CALL ZTRMM( 'Left', UPLO, 'No transpose', 'Non-unit',         &
     &                        K-1, KB, CONE, B, LDB, A( 1, K ), LDA )
                  CALL ZHEMM( 'Right', UPLO, K-1, KB, HALF, A( K, K ),          &
     &                        LDA, B( 1, K ), LDB, CONE, A( 1, K ),             &
     &                        LDA )
                  CALL ZHER2K( UPLO, 'No transpose', K-1, KB, CONE,             &
     &                         A( 1, K ), LDA, B( 1, K ), LDB, ONE, A,          &
     &                         LDA )
                  CALL ZHEMM( 'Right', UPLO, K-1, KB, HALF, A( K, K ),          &
     &                        LDA, B( 1, K ), LDB, CONE, A( 1, K ),             &
     &                        LDA )
                  CALL ZTRMM( 'Right', UPLO, 'Conjugate transpose',             &
     &                        'Non-unit', K-1, KB, CONE, B( K, K ), LDB,        &
     &                        A( 1, K ), LDA )
                  CALL ZHEGS2( ITYPE, UPLO, KB, A( K, K ), LDA,                 &
     &                         B( K, K ), LDB, INFO )
   30          CONTINUE
            ELSE
!
!              Compute L**H*A*L
!
               DO 40 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
!
!                 Update the lower triangle of A(1:k+kb-1,1:k+kb-1)
!
                  CALL ZTRMM( 'Right', UPLO, 'No transpose', 'Non-unit',        &
     &                        KB, K-1, CONE, B, LDB, A( K, 1 ), LDA )
                  CALL ZHEMM( 'Left', UPLO, KB, K-1, HALF, A( K, K ),           &
     &                        LDA, B( K, 1 ), LDB, CONE, A( K, 1 ),             &
     &                        LDA )
                  CALL ZHER2K( UPLO, 'Conjugate transpose', K-1, KB,            &
     &                         CONE, A( K, 1 ), LDA, B( K, 1 ), LDB,            &
     &                         ONE, A, LDA )
                  CALL ZHEMM( 'Left', UPLO, KB, K-1, HALF, A( K, K ),           &
     &                        LDA, B( K, 1 ), LDB, CONE, A( K, 1 ),             &
     &                        LDA )
                  CALL ZTRMM( 'Left', UPLO, 'Conjugate transpose',              &
     &                        'Non-unit', KB, K-1, CONE, B( K, K ), LDB,        &
     &                        A( K, 1 ), LDA )
                  CALL ZHEGS2( ITYPE, UPLO, KB, A( K, K ), LDA,                 &
     &                         B( K, K ), LDB, INFO )
   40          CONTINUE
            END IF
         END IF
      END IF
      RETURN
!
!     End of ZHEGST
!
      END
