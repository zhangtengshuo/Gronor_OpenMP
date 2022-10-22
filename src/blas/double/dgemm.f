#ifdef CRAY
#else
      SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!     .. Parameters ..
!     $Id: dgemm.F 19695 2010-10-29 16:51:02Z d3y133 $
      DOUBLE PRECISION ONE, ZERO
      PARAMETER        (ONE=1.0D+0,ZERO=0.0D+0)
      INTEGER          MB, NB, KB
#ifdef CACHE1M
!x86_64 cache=1024K
      PARAMETER        (MB=256,NB=MB,KB=256)
#elif CACHE6M
!ia64 cache=6M
      PARAMETER        (MB=512,NB=MB,KB=512)
#else
!x86 cache=256K
       PARAMETER        (MB=64,NB=MB,KB=64)
#endif
!     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA, BETA
      INTEGER          K, LDA, LDB, LDC, M, N
      CHARACTER        TRANSA, TRANSB
!     .. Array Arguments ..
!
!     Purpose
!     =======
!
!     DGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!     where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X',
!
!     alpha and beta are scalars, and A, B and C are matrices, with op( A )
!     an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!     Parameters
!     ==========
!
!     TRANSA - CHARACTER*1.
!     On entry, TRANSA specifies the form of op( A ) to be used in
!     the matrix multiplication as follows:
!
!     TRANSA = 'N' or 'n',  op( A ) = A.
!
!     TRANSA = 'T' or 't',  op( A ) = A'.
!
!     TRANSA = 'C' or 'c',  op( A ) = A'.
!
!     Unchanged on exit.
!
!     TRANSB - CHARACTER*1.
!     On entry, TRANSB specifies the form of op( B ) to be used in
!     the matrix multiplication as follows:
!
!     TRANSB = 'N' or 'n',  op( B ) = B.
!
!     TRANSB = 'T' or 't',  op( B ) = B'.
!
!     TRANSB = 'C' or 'c',  op( B ) = B'.
!
!     Unchanged on exit.
!
!     M      - INTEGER.
!     On entry,  M  specifies  the number  of rows  of the  matrix
!     op( A )  and of the  matrix  C.  M  must  be at least  zero.
!     Unchanged on exit.
!
!     N      - INTEGER.
!     On entry,  N  specifies the number  of columns of the matrix
!     op( B ) and the number of columns of the matrix C. N must be
!     at least zero.
!     Unchanged on exit.
!
!     K      - INTEGER.
!     On entry,  K  specifies  the number of columns of the matrix
!     op( A ) and the number of rows of the matrix op( B ). K must
!     be at least  zero.
!     Unchanged on exit.
!
!     ALPHA  - DOUBLE PRECISION.
!     On entry, ALPHA specifies the scalar alpha.
!     Unchanged on exit.
!
!     A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
!     k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!     Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!     part of the array  A  must contain the matrix  A,  otherwise
!     the leading  k by m  part of the array  A  must contain  the
!     matrix A.
!     Unchanged on exit.
!
!     LDA    - INTEGER.
!     On entry, LDA specifies the first dimension of A as declared
!     in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!     LDA must be at least  max( 1, m ), otherwise  LDA must be at
!     least  max( 1, k ).
!     Unchanged on exit.
!
!     B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
!     n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!     Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!     part of the array  B  must contain the matrix  B,  otherwise
!     the leading  n by k  part of the array  B  must contain  the
!     matrix B.
!     Unchanged on exit.
!
!     LDB    - INTEGER.
!     On entry, LDB specifies the first dimension of B as declared
!     in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!     LDB must be at least  max( 1, k ), otherwise  LDB must be at
!     least  max( 1, n ).
!     Unchanged on exit.
!
!     BETA   - DOUBLE PRECISION.
!     On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!     supplied as zero then C need not be set on input.
!     Unchanged on exit.
!
!     C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
!     Before entry, the leading  m by n  part of the array  C must
!     contain the matrix  C,  except when  beta  is zero, in which
!     case C need not be set on entry.
!     On exit, the array  C  is overwritten by the  m by n  matrix
!     ( alpha*op( A )*op( B ) + beta*C ).
!
!     LDC    - INTEGER.
!     On entry, LDC specifies the first dimension of C as declared
!     in  the  calling  (sub)  program.   LDC  must  be  at  least
!     max( 1, m ).
!     Unchanged on exit.
!
!
!     Level 3 Blas routine.
!
!     -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!     This code comes from a report entitled:
!     The IBM RISC System/6000 and Linear Algebra Operations, by
!     Jack J. Dongarra, Peter Mayes, and Giuseppe Radicati di Brozolo,
!     University of Tennessee Computer Science Tech Report: CS - 90 - 122.
!
!
      DOUBLE PRECISION A(LDA,*), B(LDB,*), C(LDC,*)
!     .. External Functions ..
      LOGICAL          LSAME
      EXTERNAL         LSAME
!     .. External Subroutines ..
      EXTERNAL         XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC        MAX, MIN
!     .. Local Scalars ..
      DOUBLE PRECISION T11, T12, T21, T22
      INTEGER          I, IDEPTH, II, ILEN, INFO, ISPAN, J, JDEPTH, JJ,         &
     &                 JLEN, JSPAN, L, LL, LSPAN, NCOLA, NROWA, NROWB
      LOGICAL          NOTA, NOTB
!     .. Local Arrays ..
      DOUBLE PRECISION CH(KB,MB), CH1(KB), CH2(KB)
!     .. Executable Statements ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.
!
      NOTA = LSAME(TRANSA,'N')
      NOTB = LSAME(TRANSB,'N')
      IF (NOTA) THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF (NOTB) THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
!
!     Test the input parameters.
!
      INFO = 0
      IF (( .NOT. NOTA) .AND. ( .NOT. LSAME(TRANSA,'C'))                        &
     &    .AND. ( .NOT. LSAME(TRANSA,'T'))) THEN
         INFO = 1
      ELSE IF (( .NOT. NOTB) .AND. ( .NOT. LSAME(TRANSB,'C'))                   &
     &         .AND. ( .NOT. LSAME(TRANSB,'T'))) THEN
         INFO = 2
      ELSE IF (M.LT.0) THEN
         INFO = 3
      ELSE IF (N.LT.0) THEN
         INFO = 4
      ELSE IF (K.LT.0) THEN
         INFO = 5
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
         INFO = 8
      ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
         INFO = 10
      ELSE IF (LDC.LT.MAX(1,M)) THEN
         INFO = 13
      END IF
      IF (INFO.NE.0) THEN
         CALL XERBLA('DGEMM ',INFO)
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (((ALPHA.EQ.ZERO) .OR. (K.EQ.0))          &
     &     .AND. (BETA.EQ.ONE))) RETURN
      IF (BETA.EQ.ZERO) THEN
         DO 40 J = 1, N
            DO 20 I = 1, M
               C(I,J) = ZERO
   20       CONTINUE
   40    CONTINUE
      ELSE
         DO 80 J = 1, N
            DO 60 I = 1, M
               C(I,J) = BETA*C(I,J)
   60       CONTINUE
   80    CONTINUE
      END IF
!
!     And if  alpha.eq.zero.
!
      IF (ALPHA.EQ.ZERO) RETURN
!
!     Start the operations.
!
      IF (NOTB) THEN
         IF (NOTA) THEN
!
!           Form  C := C + alpha*A*B.
!
            DO 380 L = 1, K, KB
               LSPAN = MIN(KB,K-L+1)
               DO 360 I = 1, M, MB
                  IDEPTH = 2
                  ISPAN = MIN(MB,M-I+1)
                  ILEN = IDEPTH*(ISPAN/IDEPTH)
                  DO 120 II = I, I + ISPAN - 1
                     DO 100 LL = L, L + LSPAN - 1
                        CH(LL-L+1,II-I+1) = ALPHA*A(II,LL)
  100                CONTINUE
  120             CONTINUE
                  DO 340 J = 1, N, NB
                     JDEPTH = 2
                     JSPAN = MIN(NB,N-J+1)
                     JLEN = JDEPTH*(JSPAN/JDEPTH)
                     DO 220 JJ = J, J + JLEN - 1, JDEPTH
                        DO 160 II = I, I + ILEN - 1, IDEPTH
                           T11 = ZERO
                           T21 = ZERO
                           T12 = ZERO
                           T22 = ZERO
                           DO 140 LL = L, L + LSPAN - 1
                              T11 = T11 + CH(LL-L+1,II-I+1)*B(LL,JJ)
                              T21 = T21 + CH(LL-L+1,II-I+2)*B(LL,JJ)
                              T12 = T12 + CH(LL-L+1,II-I+1)*B(LL,JJ+1)
                              T22 = T22 + CH(LL-L+1,II-I+2)*B(LL,JJ+1)
  140                      CONTINUE
                           C(II,JJ) = C(II,JJ) + T11
                           C(II+1,JJ) = C(II+1,JJ) + T21
                           C(II,JJ+1) = C(II,JJ+1) + T12
                           C(II+1,JJ+1) = C(II+1,JJ+1) + T22
  160                   CONTINUE
                        IF (ILEN.LT.ISPAN) THEN
                           DO 200 II = I + ILEN, I + ISPAN - 1
                              T11 = ZERO
                              T12 = ZERO
                              DO 180 LL = L, L + LSPAN - 1
                                 T11 = T11 + CH(LL-L+1,II-I+1)*B(LL,JJ)
                                 T12 = T12 + CH(LL-L+1,II-I+1)*B(LL,            &
     &                                 JJ+1)
  180                         CONTINUE
                              C(II,JJ) = C(II,JJ) + T11
                              C(II,JJ+1) = C(II,JJ+1) + T12
  200                      CONTINUE
                        END IF
  220                CONTINUE
                     IF (JLEN.LT.JSPAN) THEN
                        DO 320 JJ = J + JLEN, J + JSPAN - 1
                           DO 260 II = I, I + ILEN - 1, IDEPTH
                              T11 = ZERO
                              T21 = ZERO
                              DO 240 LL = L, L + LSPAN - 1
                                 T11 = T11 + CH(LL-L+1,II-I+1)*B(LL,JJ)
                                 T21 = T21 + CH(LL-L+1,II-I+2)*B(LL,JJ)
  240                         CONTINUE
                              C(II,JJ) = C(II,JJ) + T11
                              C(II+1,JJ) = C(II+1,JJ) + T21
  260                      CONTINUE
                           IF (ILEN.LT.ISPAN) THEN
                              DO 300 II = I + ILEN, I + ISPAN - 1
                                 T11 = ZERO
                                 DO 280 LL = L, L + LSPAN - 1
                                    T11 = T11 + CH(LL-L+1,II-I+1)*B(LL,         &
     &                                    JJ)
  280                            CONTINUE
                                 C(II,JJ) = C(II,JJ) + T11
  300                         CONTINUE
                           END IF
  320                   CONTINUE
                     END IF
  340             CONTINUE
  360          CONTINUE
  380       CONTINUE
         ELSE
!
!           Form  C := C + alpha*A'*B
!
            DO 680 I = 1, M, MB
               IDEPTH = 2
               ISPAN = MIN(MB,M-I+1)
               ILEN = IDEPTH*(ISPAN/IDEPTH)
               DO 660 L = 1, K, KB
                  LSPAN = MIN(KB,K-L+1)
                  DO 420 II = I, I + ISPAN - 1
                     DO 400 LL = L, L + LSPAN - 1
                        CH(LL-L+1,II-I+1) = ALPHA*A(LL,II)
  400                CONTINUE
  420             CONTINUE
                  DO 640 J = 1, N, NB
                     JDEPTH = 2
                     JSPAN = MIN(NB,N-J+1)
                     JLEN = JDEPTH*(JSPAN/JDEPTH)
                     DO 520 JJ = J, J + JLEN - 1, JDEPTH
                        DO 460 II = I, I + ILEN - 1, IDEPTH
                           T11 = ZERO
                           T21 = ZERO
                           T12 = ZERO
                           T22 = ZERO
                           DO 440 LL = L, L + LSPAN - 1
                              T11 = T11 + CH(LL-L+1,II-I+1)*B(LL,JJ)
                              T21 = T21 + CH(LL-L+1,II-I+2)*B(LL,JJ)
                              T12 = T12 + CH(LL-L+1,II-I+1)*B(LL,JJ+1)
                              T22 = T22 + CH(LL-L+1,II-I+2)*B(LL,JJ+1)
  440                      CONTINUE
                           C(II,JJ) = C(II,JJ) + T11
                           C(II+1,JJ) = C(II+1,JJ) + T21
                           C(II,JJ+1) = C(II,JJ+1) + T12
                           C(II+1,JJ+1) = C(II+1,JJ+1) + T22
  460                   CONTINUE
                        IF (ILEN.LT.ISPAN) THEN
                           DO 500 II = I + ILEN, I + ISPAN - 1
                              T11 = ZERO
                              T12 = ZERO
                              DO 480 LL = L, L + LSPAN - 1
                                 T11 = T11 + CH(LL-L+1,II-I+1)*B(LL,JJ)
                                 T12 = T12 + CH(LL-L+1,II-I+1)*B(LL,            &
     &                                 JJ+1)
  480                         CONTINUE
                              C(II,JJ) = C(II,JJ) + T11
                              C(II,JJ+1) = C(II,JJ+1) + T12
  500                      CONTINUE
                        END IF
  520                CONTINUE
                     IF (JLEN.LT.JSPAN) THEN
                        DO 620 JJ = J + JLEN, J + JSPAN - 1
                           DO 560 II = I, I + ILEN - 1, IDEPTH
                              T11 = ZERO
                              T21 = ZERO
                              DO 540 LL = L, L + LSPAN - 1
                                 T11 = T11 + CH(LL-L+1,II-I+1)*B(LL,JJ)
                                 T21 = T21 + CH(LL-L+1,II-I+2)*B(LL,JJ)
  540                         CONTINUE
                              C(II,JJ) = C(II,JJ) + T11
                              C(II+1,JJ) = C(II+1,JJ) + T21
  560                      CONTINUE
                           IF (ILEN.LT.ISPAN) THEN
                              DO 600 II = I + ILEN, I + ISPAN - 1
                                 T11 = ZERO
                                 DO 580 LL = L, L + LSPAN - 1
                                    T11 = T11 + CH(LL-L+1,II-I+1)*B(LL,         &
     &                                    JJ)
  580                            CONTINUE
                                 C(II,JJ) = C(II,JJ) + T11
  600                         CONTINUE
                           END IF
  620                   CONTINUE
                     END IF
  640             CONTINUE
  660          CONTINUE
  680       CONTINUE
         END IF
      ELSE
         IF (NOTA) THEN
!
!           Form  C := C + alpha*A*B'
!
            DO 1000 J = 1, N, NB
               JDEPTH = 2
               JSPAN = MIN(NB,N-J+1)
               JLEN = JDEPTH*(JSPAN/JDEPTH)
               DO 980 L = 1, K, KB
                  LSPAN = MIN(KB,K-L+1)
                  DO 720 JJ = J, J + JSPAN - 1
                     DO 700 LL = L, L + LSPAN - 1
                        CH(LL-L+1,JJ-J+1) = ALPHA*B(JJ,LL)
  700                CONTINUE
  720             CONTINUE
                  DO 960 I = 1, M, MB
                     IDEPTH = 2
                     ISPAN = MIN(MB,M-I+1)
                     ILEN = IDEPTH*(ISPAN/IDEPTH)
                     DO 840 II = I, I + ILEN - 1, IDEPTH
                        DO 740 LL = L, L + LSPAN - 1
                           CH1(LL-L+1) = A(II,LL)
                           CH2(LL-L+1) = A(II+1,LL)
  740                   CONTINUE
                        DO 780 JJ = J, J + JLEN - 1, JDEPTH
                           T11 = ZERO
                           T21 = ZERO
                           T12 = ZERO
                           T22 = ZERO
                           DO 760 LL = L, L + LSPAN - 1
                              T11 = T11 + CH1(LL-L+1)*CH(LL-L+1,JJ-J+1)
                              T21 = T21 + CH2(LL-L+1)*CH(LL-L+1,JJ-J+1)
                              T12 = T12 + CH1(LL-L+1)*CH(LL-L+1,JJ-J+2)
                              T22 = T22 + CH2(LL-L+1)*CH(LL-L+1,JJ-J+2)
  760                      CONTINUE
                           C(II,JJ) = C(II,JJ) + T11
                           C(II+1,JJ) = C(II+1,JJ) + T21
                           C(II,JJ+1) = C(II,JJ+1) + T12
                           C(II+1,JJ+1) = C(II+1,JJ+1) + T22
  780                   CONTINUE
                        IF (JLEN.LT.JSPAN) THEN
                           DO 820 JJ = J + JLEN, J + JSPAN - 1
                              T11 = ZERO
                              T21 = ZERO
                              DO 800 LL = L, L + LSPAN - 1
                                 T11 = T11 + A(II,LL)*CH(LL-L+1,JJ-J+1)
                                 T21 = T21 + A(II+1,LL)*CH(LL-L+1,              &
     &                                 JJ-J+1)
  800                         CONTINUE
                              C(II,JJ) = C(II,JJ) + T11
                              C(II+1,JJ) = C(II+1,JJ) + T21
  820                      CONTINUE
                        END IF
  840                CONTINUE
                     IF (ILEN.LT.ISPAN) THEN
                        DO 940 II = I + ILEN, I + ISPAN - 1
                           DO 880 JJ = J, J + JLEN - 1, JDEPTH
                              T11 = ZERO
                              T12 = ZERO
                              DO 860 LL = L, L + LSPAN - 1
                                 T11 = T11 + A(II,LL)*CH(LL-L+1,JJ-J+1)
                                 T12 = T12 + A(II,LL)*CH(LL-L+1,JJ-J+2)
  860                         CONTINUE
                              C(II,JJ) = C(II,JJ) + T11
                              C(II,JJ+1) = C(II,JJ+1) + T12
  880                      CONTINUE
                           IF (JLEN.LT.JSPAN) THEN
                              DO 920 JJ = J + JLEN, J + JSPAN - 1
                                 T11 = ZERO
                                 DO 900 LL = L, L + LSPAN - 1
                                    T11 = T11 + A(II,LL)*CH(LL-L+1,             &
     &                                    JJ-J+1)
  900                            CONTINUE
                                 C(II,JJ) = C(II,JJ) + T11
  920                         CONTINUE
                           END IF
  940                   CONTINUE
                     END IF
  960             CONTINUE
  980          CONTINUE
 1000       CONTINUE
         ELSE
!
!           Form  C := C + alpha*A'*B'
!
            DO 1300 J = 1, N, NB
               JDEPTH = 2
               JSPAN = MIN(NB,N-J+1)
               JLEN = JDEPTH*(JSPAN/JDEPTH)
               DO 1280 L = 1, K, KB
                  LSPAN = MIN(KB,K-L+1)
                  DO 1040 JJ = J, J + JSPAN - 1
                     DO 1020 LL = L, L + LSPAN - 1
                        CH(LL-L+1,JJ-J+1) = ALPHA*B(JJ,LL)
 1020                CONTINUE
 1040             CONTINUE
                  DO 1260 I = 1, M, MB
                     IDEPTH = 2
                     ISPAN = MIN(MB,M-I+1)
                     ILEN = IDEPTH*(ISPAN/IDEPTH)
                     DO 1140 II = I, I + ILEN - 1, IDEPTH
                        DO 1080 JJ = J, J + JLEN - 1, JDEPTH
                           T11 = ZERO
                           T21 = ZERO
                           T12 = ZERO
                           T22 = ZERO
                           DO 1060 LL = L, L + LSPAN - 1
                              T11 = T11 + A(LL,II)*CH(LL-L+1,JJ-J+1)
                              T21 = T21 + A(LL,II+1)*CH(LL-L+1,JJ-J+1)
                              T12 = T12 + A(LL,II)*CH(LL-L+1,JJ-J+2)
                              T22 = T22 + A(LL,II+1)*CH(LL-L+1,JJ-J+2)
 1060                      CONTINUE
                           C(II,JJ) = C(II,JJ) + T11
                           C(II+1,JJ) = C(II+1,JJ) + T21
                           C(II,JJ+1) = C(II,JJ+1) + T12
                           C(II+1,JJ+1) = C(II+1,JJ+1) + T22
 1080                   CONTINUE
                        IF (JLEN.LT.JSPAN) THEN
                           DO 1120 JJ = J + JLEN, J + JSPAN - 1
                              T11 = ZERO
                              T21 = ZERO
                              DO 1100 LL = L, L + LSPAN - 1
                                 T11 = T11 + A(LL,II)*CH(LL-L+1,JJ-J+1)
                                 T21 = T21 + A(LL,II+1)*CH(LL-L+1,              &
     &                                 JJ-J+1)
 1100                         CONTINUE
                              C(II,JJ) = C(II,JJ) + T11
                              C(II+1,JJ) = C(II+1,JJ) + T21
 1120                      CONTINUE
                        END IF
 1140                CONTINUE
                     IF (ILEN.LT.ISPAN) THEN
                        DO 1240 II = I + ILEN, I + ISPAN - 1
                           DO 1180 JJ = J, J + JLEN - 1, JDEPTH
                              T11 = ZERO
                              T12 = ZERO
                              DO 1160 LL = L, L + LSPAN - 1
                                 T11 = T11 + A(LL,II)*CH(LL-L+1,JJ-J+1)
                                 T12 = T12 + A(LL,II)*CH(LL-L+1,JJ-J+2)
 1160                         CONTINUE
                              C(II,JJ) = C(II,JJ) + T11
                              C(II,JJ+1) = C(II,JJ+1) + T12
 1180                      CONTINUE
                           IF (JLEN.LT.JSPAN) THEN
                              DO 1220 JJ = J + JLEN, J + JSPAN - 1
                                 T11 = ZERO
                                 DO 1200 LL = L, L + LSPAN - 1
                                    T11 = T11 + A(LL,II)*CH(LL-L+1,             &
     &                                    JJ-J+1)
 1200                            CONTINUE
                                 C(II,JJ) = C(II,JJ) + T11
 1220                         CONTINUE
                           END IF
 1240                   CONTINUE
                     END IF
 1260             CONTINUE
 1280          CONTINUE
 1300       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of DGEMM .
!
      END
#endif
