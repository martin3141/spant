C     ==================================================================
C     QR pre-reduction wrapper for PNNLS (see pnnls.f).
C
C     PNNLS solves the partial NNLS problem directly on the full M by N
C     matrix A, applying a Householder transformation to every column
C     of A each time a variable is moved between the active and passive
C     sets.  When M is much larger than N (the usual case for spant's
C     spectral fitting problems, where M is the number of data points
C     and N the number of basis/spline functions) this does far more
C     work than necessary.
C
C     Because the NNLS active-set iteration only depends on A and B
C     through inner products, and ||A x - b|| = ||Q A x - Q b|| for any
C     orthogonal Q, the M by N problem can be reduced to an N by N
C     triangular one via a QR factorization A = Q R without changing
C     the solution.  PNNLS is then run on R (N by N) and Q'B (first N
C     elements) instead of on the full A and B, which cuts the cost of
C     each Householder update from O(M) to O(N).  The discarded tail of
C     Q'B (elements N+1..M) contributes a fixed amount to the residual
C     norm, which is added back in after PNNLS returns.
C
C     Arguments and their meanings are identical to PNNLS, with the
C     added requirement that M >= N (the caller is expected to enforce
C     this; it is not re-checked here).  MODE is set to 4 if the LAPACK
C     QR factorization step fails, in addition to the failure modes
C     documented in pnnls.f.
C     ==================================================================
      SUBROUTINE QRPNLS (A,MDA,M,N,B,X,RNORM,W,ZZ,INDEX,MODE,K)
C     ------------------------------------------------------------------
      integer MDA, M, N, MODE, K, INDEX(*)
      integer I, J, INFO, LWORK
      double precision A(MDA,*), B(*), X(*), W(*), ZZ(*), RNORM
      double precision TAU(N), WORK(64*N), SM
C     ------------------------------------------------------------------
      LWORK = 64*N

C     A = Q * R ; R is left in the upper triangle of A, and the
C     Householder vectors that implicitly define Q are left below it.
      CALL DGEQRF (M,N,A,MDA,TAU,WORK,LWORK,INFO)
      IF (INFO .ne. 0) then
         MODE = 4
         RETURN
      endif

C     B <- Q**T * B
      CALL DORMQR ('L','T',M,1,N,A,MDA,TAU,B,M,WORK,LWORK,INFO)
      IF (INFO .ne. 0) then
         MODE = 4
         RETURN
      endif

C     The tail of Q**T B (rows N+1..M) is orthogonal to the column
C     space of A and so contributes a fixed amount to the residual
C     norm of any solution; save it before B is overwritten by PNNLS.
      SM = 0.0d0
      DO 20 I = N+1, M
         SM = SM + B(I)*B(I)
 20   continue

C     Zero the strict lower triangle so A's leading N by N block holds
C     R alone (the Householder vectors below it must not be mistaken
C     for data by PNNLS).
      DO 40 J = 1, N
         DO 30 I = J+1, N
            A(I,J) = 0.0d0
 30      continue
 40   continue

      CALL PNNLS (A,MDA,N,N,B,X,RNORM,W,ZZ,INDEX,MODE,K)

      RNORM = sqrt(RNORM*RNORM + SM)

      RETURN
      END
