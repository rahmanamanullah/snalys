C     ****f* solve/NEWRAPH *
C
C   NAME
C     NEWRAPH - Newton-Raphsons method
C
C   DESCRIPTION
C     Simple implementation of Newton-Raphsons one-dimensional
C     equation solver.
C
C   USAGE
C     CALL NEWRAPH(FUNCD, X, DIR, X1, X2, XACC)
C
C   INPUTS
C     FUNCD - to be solved, takes two argument and returns
C             the function value. In the second argument the
C             derative of the function should be on returned.
C     X     - REAL
C     X1    - lower X limit
C     X2    - higher X limit
C     XACC  - requested accuracy.
C
C   SIDE EFFECTS
C     On return X will contain the solution.
C
C   SEE ALSO
C     Numerical Recepies
C
C     ***
      SUBROUTINE NEWRAPH(FUNCD, X, X1, X2, XACC)
      IMPLICIT NONE
      REAL FUNCD
      EXTERNAL FUNCD
      REAL X, X1, X2, XACC

      INTEGER NMAX
      PARAMETER (NMAX=30)
      INTEGER N
      REAL DF,DX,F, RTNEWT

      RTNEWT = .5*(X1 + X2)
      DO 11 N = 1, NMAX
        F = FUNCD(RTNEWT,DF)
        DX=F/DF
        RTNEWT=RTNEWT-DX
        IF((X1-RTNEWT)*(RTNEWT-X2).LT.0.)PAUSE
     $       'rtnewt jumped out of brackets'
        IF(ABS(DX).LT.XACC) RETURN
11    CONTINUE
      PAUSE 'RTNEWT exceeded maximum iterations'
      X = RTNEWT
      END
C
C End NEWRAPH
C
