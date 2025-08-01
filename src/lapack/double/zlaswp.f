!> \brief \b ZLASWP performs a series of row interchanges on a general rectangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLASWP + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaswp.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaswp.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaswp.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLASWP( N, A, LDA, K1, K2, IPIV, INCX )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, K1, K2, LDA, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX*16         A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLASWP performs a series of row interchanges on the matrix A.
!> One row interchange is initiated for each of rows K1 through K2 of A.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the matrix of column dimension N to which the row
!>          interchanges will be applied.
!>          On exit, the permuted matrix.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!> \endverbatim
!>
!> \param[in] K1
!> \verbatim
!>          K1 is INTEGER
!>          The first element of IPIV for which a row interchange will
!>          be done.
!> \endverbatim
!>
!> \param[in] K2
!> \verbatim
!>          K2 is INTEGER
!>          The last element of IPIV for which a row interchange will
!>          be done.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (K2*abs(INCX))
!>          The vector of pivot indices.  Only the elements in positions
!>          K1 through K2 of IPIV are accessed.
!>          IPIV(K) = L implies rows K and L are to be interchanged.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The increment between successive values of IPIV.  If IPIV
!>          is negative, the pivots are applied in reverse order.
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
!> \ingroup complex16OTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Modified by
!>   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZLASWP( N, A, LDA, K1, K2, IPIV, INCX )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
!     ..
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      COMPLEX*16         TEMP
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(I) for each of rows K1 through K2.
!
      IF( INCX.GT.0 ) THEN
         IX0 = K1
         I1 = K1
         I2 = K2
         INC = 1
      ELSE IF( INCX.LT.0 ) THEN
         IX0 = 1 + ( 1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      ELSE
         RETURN
      END IF
!
      N32 = ( N / 32 )*32
      IF( N32.NE.0 ) THEN
         DO 30 J = 1, N32, 32
            IX = IX0
            DO 20 I = I1, I2, INC
               IP = IPIV( IX )
               IF( IP.NE.I ) THEN
                  DO 10 K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
   10             CONTINUE
               END IF
               IX = IX + INCX
   20       CONTINUE
   30    CONTINUE
      END IF
      IF( N32.NE.N ) THEN
         N32 = N32 + 1
         IX = IX0
         DO 50 I = I1, I2, INC
            IP = IPIV( IX )
            IF( IP.NE.I ) THEN
               DO 40 K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
   40          CONTINUE
            END IF
            IX = IX + INCX
   50    CONTINUE
      END IF
!
      RETURN
!
!     End of ZLASWP
!
      END
