!> \brief \b DSPTRS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DSPTRS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsptrs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsptrs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsptrs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       DOUBLE PRECISION   AP( * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSPTRS solves a system of linear equations A*X = B with a real
!> symmetric matrix A stored in packed format using the factorization
!> A = U*D*U**T or A = L*D*L**T computed by DSPTRF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the details of the factorization are stored
!>          as an upper or lower triangular matrix.
!>          = 'U':  Upper triangular, form is A = U*D*U**T;
!>          = 'L':  Lower triangular, form is A = L*D*L**T.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrix B.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] AP
!> \verbatim
!>          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2)
!>          The block diagonal matrix D and the multipliers used to
!>          obtain the factor U or L as computed by DSPTRF, stored as a
!>          packed triangular matrix.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          Details of the interchanges and the block structure of D
!>          as determined by DSPTRF.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
!>          On entry, the right hand side matrix B.
!>          On exit, the solution matrix X.
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
!>          < 0: if INFO = -i, the i-th argument had an illegal value
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
!> \ingroup doubleOTHERcomputational
!
!  =====================================================================
      SUBROUTINE DSPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   AP( * ), B( LDB, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, K, KC, KP
      DOUBLE PRECISION   AK, AKM1, AKM1K, BK, BKM1, DENOM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMV, DGER, DSCAL, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSPTRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 )                                               &
     &   RETURN
!
      IF( UPPER ) THEN
!
!        Solve A*X = B, where A = U*D*U**T.
!
!        First solve U*D*X = B, overwriting B with X.
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         K = N
         KC = N*( N+1 ) / 2 + 1
   10    CONTINUE
!
!        If K < 1, exit from loop.
!
         IF( K.LT.1 )                                                           &
     &      GO TO 30
!
         KC = KC - K
         IF( IPIV( K ).GT.0 ) THEN
!
!           1 x 1 diagonal block
!
!           Interchange rows K and IPIV(K).
!
            KP = IPIV( K )
            IF( KP.NE.K )                                                       &
     &         CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
!
!           Multiply by inv(U(K)), where U(K) is the transformation
!           stored in column K of A.
!
            CALL DGER( K-1, NRHS, -ONE, AP( KC ), 1, B( K, 1 ), LDB,            &
     &                 B( 1, 1 ), LDB )
!
!           Multiply by the inverse of the diagonal block.
!
            CALL DSCAL( NRHS, ONE / AP( KC+K-1 ), B( K, 1 ), LDB )
            K = K - 1
         ELSE
!
!           2 x 2 diagonal block
!
!           Interchange rows K-1 and -IPIV(K).
!
            KP = -IPIV( K )
            IF( KP.NE.K-1 )                                                     &
     &         CALL DSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB )
!
!           Multiply by inv(U(K)), where U(K) is the transformation
!           stored in columns K-1 and K of A.
!
            CALL DGER( K-2, NRHS, -ONE, AP( KC ), 1, B( K, 1 ), LDB,            &
     &                 B( 1, 1 ), LDB )
            CALL DGER( K-2, NRHS, -ONE, AP( KC-( K-1 ) ), 1,                    &
     &                 B( K-1, 1 ), LDB, B( 1, 1 ), LDB )
!
!           Multiply by the inverse of the diagonal block.
!
            AKM1K = AP( KC+K-2 )
            AKM1 = AP( KC-1 ) / AKM1K
            AK = AP( KC+K-1 ) / AKM1K
            DENOM = AKM1*AK - ONE
            DO 20 J = 1, NRHS
               BKM1 = B( K-1, J ) / AKM1K
               BK = B( K, J ) / AKM1K
               B( K-1, J ) = ( AK*BKM1-BK ) / DENOM
               B( K, J ) = ( AKM1*BK-BKM1 ) / DENOM
   20       CONTINUE
            KC = KC - K + 1
            K = K - 2
         END IF
!
         GO TO 10
   30    CONTINUE
!
!        Next solve U**T*X = B, overwriting B with X.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         K = 1
         KC = 1
   40    CONTINUE
!
!        If K > N, exit from loop.
!
         IF( K.GT.N )                                                           &
     &      GO TO 50
!
         IF( IPIV( K ).GT.0 ) THEN
!
!           1 x 1 diagonal block
!
!           Multiply by inv(U**T(K)), where U(K) is the transformation
!           stored in column K of A.
!
            CALL DGEMV( 'Transpose', K-1, NRHS, -ONE, B, LDB, AP( KC ),         &
     &                  1, ONE, B( K, 1 ), LDB )
!
!           Interchange rows K and IPIV(K).
!
            KP = IPIV( K )
            IF( KP.NE.K )                                                       &
     &         CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            KC = KC + K
            K = K + 1
         ELSE
!
!           2 x 2 diagonal block
!
!           Multiply by inv(U**T(K+1)), where U(K+1) is the transformation
!           stored in columns K and K+1 of A.
!
            CALL DGEMV( 'Transpose', K-1, NRHS, -ONE, B, LDB, AP( KC ),         &
     &                  1, ONE, B( K, 1 ), LDB )
            CALL DGEMV( 'Transpose', K-1, NRHS, -ONE, B, LDB,                   &
     &                  AP( KC+K ), 1, ONE, B( K+1, 1 ), LDB )
!
!           Interchange rows K and -IPIV(K).
!
            KP = -IPIV( K )
            IF( KP.NE.K )                                                       &
     &         CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            KC = KC + 2*K + 1
            K = K + 2
         END IF
!
         GO TO 40
   50    CONTINUE
!
      ELSE
!
!        Solve A*X = B, where A = L*D*L**T.
!
!        First solve L*D*X = B, overwriting B with X.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         K = 1
         KC = 1
   60    CONTINUE
!
!        If K > N, exit from loop.
!
         IF( K.GT.N )                                                           &
     &      GO TO 80
!
         IF( IPIV( K ).GT.0 ) THEN
!
!           1 x 1 diagonal block
!
!           Interchange rows K and IPIV(K).
!
            KP = IPIV( K )
            IF( KP.NE.K )                                                       &
     &         CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
!
!           Multiply by inv(L(K)), where L(K) is the transformation
!           stored in column K of A.
!
            IF( K.LT.N )                                                        &
     &         CALL DGER( N-K, NRHS, -ONE, AP( KC+1 ), 1, B( K, 1 ),            &
     &                    LDB, B( K+1, 1 ), LDB )
!
!           Multiply by the inverse of the diagonal block.
!
            CALL DSCAL( NRHS, ONE / AP( KC ), B( K, 1 ), LDB )
            KC = KC + N - K + 1
            K = K + 1
         ELSE
!
!           2 x 2 diagonal block
!
!           Interchange rows K+1 and -IPIV(K).
!
            KP = -IPIV( K )
            IF( KP.NE.K+1 )                                                     &
     &         CALL DSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB )
!
!           Multiply by inv(L(K)), where L(K) is the transformation
!           stored in columns K and K+1 of A.
!
            IF( K.LT.N-1 ) THEN
               CALL DGER( N-K-1, NRHS, -ONE, AP( KC+2 ), 1, B( K, 1 ),          &
     &                    LDB, B( K+2, 1 ), LDB )
               CALL DGER( N-K-1, NRHS, -ONE, AP( KC+N-K+2 ), 1,                 &
     &                    B( K+1, 1 ), LDB, B( K+2, 1 ), LDB )
            END IF
!
!           Multiply by the inverse of the diagonal block.
!
            AKM1K = AP( KC+1 )
            AKM1 = AP( KC ) / AKM1K
            AK = AP( KC+N-K+1 ) / AKM1K
            DENOM = AKM1*AK - ONE
            DO 70 J = 1, NRHS
               BKM1 = B( K, J ) / AKM1K
               BK = B( K+1, J ) / AKM1K
               B( K, J ) = ( AK*BKM1-BK ) / DENOM
               B( K+1, J ) = ( AKM1*BK-BKM1 ) / DENOM
   70       CONTINUE
            KC = KC + 2*( N-K ) + 1
            K = K + 2
         END IF
!
         GO TO 60
   80    CONTINUE
!
!        Next solve L**T*X = B, overwriting B with X.
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         K = N
         KC = N*( N+1 ) / 2 + 1
   90    CONTINUE
!
!        If K < 1, exit from loop.
!
         IF( K.LT.1 )                                                           &
     &      GO TO 100
!
         KC = KC - ( N-K+1 )
         IF( IPIV( K ).GT.0 ) THEN
!
!           1 x 1 diagonal block
!
!           Multiply by inv(L**T(K)), where L(K) is the transformation
!           stored in column K of A.
!
            IF( K.LT.N )                                                        &
     &         CALL DGEMV( 'Transpose', N-K, NRHS, -ONE, B( K+1, 1 ),           &
     &                     LDB, AP( KC+1 ), 1, ONE, B( K, 1 ), LDB )
!
!           Interchange rows K and IPIV(K).
!
            KP = IPIV( K )
            IF( KP.NE.K )                                                       &
     &         CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K = K - 1
         ELSE
!
!           2 x 2 diagonal block
!
!           Multiply by inv(L**T(K-1)), where L(K-1) is the transformation
!           stored in columns K-1 and K of A.
!
            IF( K.LT.N ) THEN
               CALL DGEMV( 'Transpose', N-K, NRHS, -ONE, B( K+1, 1 ),           &
     &                     LDB, AP( KC+1 ), 1, ONE, B( K, 1 ), LDB )
               CALL DGEMV( 'Transpose', N-K, NRHS, -ONE, B( K+1, 1 ),           &
     &                     LDB, AP( KC-( N-K ) ), 1, ONE, B( K-1, 1 ),          &
     &                     LDB )
            END IF
!
!           Interchange rows K and -IPIV(K).
!
            KP = -IPIV( K )
            IF( KP.NE.K )                                                       &
     &         CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            KC = KC - ( N-K+2 )
            K = K - 2
         END IF
!
         GO TO 90
  100    CONTINUE
      END IF
!
      RETURN
!
!     End of DSPTRS
!
      END
