C     ****h* snalys/powell [1.0]
C
C   NAME
C     Powell minimization -- a minimization
C
C   DESCRIPTION
C     This routine minimizes a function of a arbitrary number
C     of variables. In this case the negative logarithm of the
C     ML-function.
C
C   SEE ALSO
C     "Numerical Recepies". The routine is described in detail at page
C    406 section 10.5.
C
C     ***



C     ****f* snalys/POWELL *
C
C   NAME
C     POWELL -- minimization of multidimensional functions
C
C   DESCRIPTION
C     This routine minimizes a function of a arbitrary number
C     of variables.
C
C   USAGE
C     CALL POWELL(P,XI,N,NP,FTOL,ITER,FRET)
C
C   INPUTS
C     P    - is an array with n elements representing
C            the starting point when the routine is called.
C     XI   - is a matrix whose columns are the set of
C            n directions for which the function is to 
C            be minimized.
C     N    - is the number of variables.
C     NP   - has the same meaning as n in this context.
C     FTOL - the fractional tolerance in the function value
C            such that failure to decrease by more e than this
C            amount on one iteration signals doneness.
C     ITER - INTEGER
C     FRET - REAL
C
C   SIDE EFFECTS
C     P    - is the point of the minimum.
C     ITER - is the number of iterations performed.
C     FRET - is the value of the function in p.
C
C   WARNINGS
C     Warning message will be written if the maximal number of
C     iterations is exceeded.
C
C   TODO
C     Modify it so that it will take the function that is to
C     minimized as an argument.
C
C   USES
C     LOGMLFUNC
C
C     ***
      SUBROUTINE POWELL(P,XI,N,NP,FTOL,ITER,FRET)
      INTEGER ITER,N,NP,NMAX,ITMAX,ENR
      REAL FRET,FTOL,P(NP),XI(NP,NP),LOGMLFUNC
      EXTERNAL LOGMLFUNC
      PARAMETER (NMAX=20,ITMAX=200)
      INTEGER I,IBIG,J
      REAL DEL,FP,FPTT,T,PT(NMAX),PTT(NMAX),XIT(NMAX)
      FRET = LOGMLFUNC( P, P )
      DO 11 J=1,N
        PT(J)=P(J)
11    CONTINUE
      ITER=0
1     ITER=ITER+1
      FP=FRET
      IBIG=0
      DEL=0.
      DO 13 I=1,N
        DO 12 J=1,N
          XIT(J)=XI(J,I)
12      CONTINUE
        FPTT=FRET
        CALL LINMIN(P,XIT,N,FRET)
        IF(ABS(FPTT-FRET).GT.DEL)THEN
          DEL=ABS(FPTT-FRET)
          IBIG=I
        ENDIF
13    CONTINUE
      IF(2.*ABS(FP-FRET).LE.FTOL*(ABS(FP)+ABS(FRET)))RETURN
      IF(ITER.EQ.ITMAX) PAUSE 'powell exceeding maximum iterations'
      DO 14 J=1,N
        PTT(J)=2.*P(J)-PT(J)
        XIT(J)=P(J)-PT(J)
        PT(J)=P(J)
14    CONTINUE
      FPTT=LOGMLFUNC(PTT, PTT)
      IF(FPTT.GE.FP)GOTO 1
      T=2.*(FP-2.*FRET+FPTT)*(FP-FRET-DEL)**2-DEL*(FP-FPTT)**2
      IF(T.GE.0.)GOTO 1
      CALL LINMIN(P,XIT,N,FRE)
      DO 15 J=1,N
        XI(J,IBIG)=XI(J,N)
        XI(J,N)=XIT(J)
15    CONTINUE
      GOTO 1
      END

      SUBROUTINE LINMIN(P,XI,N,FRET)
      INTEGER N,NMAX,ENR,SN_NR
      REAL FRET,P(N),XI(N),TOL
      PARAMETER (NMAX=50,TOL=1.E-4)
CU    USES BRENT,F1DIM,MNBRAK
      INTEGER J,NCOM
      REAL AX,BX,FA,FB,FX,XMIN,XX,PCOM(NMAX),XICOM(NMAX),BRENT
      COMMON /F1COM/ PCOM,XICOM,NCOM
      EXTERNAL F1DIM
      NCOM=N
      DO 11 J=1,N
        PCOM(J)=P(J)
        XICOM(J)=XI(J)
11    CONTINUE
      AX=0.
      XX=1.
      CALL MNBRAK(AX,XX,BX,FA,FX,FB,F1DIM)
      FRET=BRENT(AX,XX,BX,F1DIM,TOL,XMIN)
      DO 12 J=1,N
        XI(J)=XMIN*XI(J)
        P(J)=P(J)+XI(J)
12    CONTINUE
      RETURN
      END

      SUBROUTINE MNBRAK(AX,BX,CX,FA,FB,FC,FUNC)
      REAL AX,BX,CX,FA,FB,FC,FUNC,GOLD,GLIMIT,TINY
      EXTERNAL FUNC
      PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.E-20)
      REAL DUM,FU,Q,R,U,ULIM
      FA=FUNC(AX)
      FB=FUNC(BX)
      IF(FB.GT.FA)THEN
        DUM=AX
        AX=BX
        BX=DUM
        DUM=FB
        FB=FA
        FA=DUM
      ENDIF
      CX=BX+GOLD*(BX-AX)
      FC=FUNC(CX)
1     IF(FB.GE.FC)THEN
        R=(BX-AX)*(FB-FC)
        Q=(BX-CX)*(FB-FA)
        U=BX-((BX-CX)*Q-(BX-AX)*R)/(2.*SIGN(MAX(ABS(Q-R),TINY),Q-R))
        ULIM=BX+GLIMIT*(CX-BX)
        IF((BX-U)*(U-CX).GT.0.)THEN
          FU=FUNC(U)
          IF(FU.LT.FC)THEN
            AX=BX
            FA=FB
            BX=U
            FB=FU
            RETURN
          ELSE IF(FU.GT.FB)THEN
            CX=U
            FC=FU
            RETURN
          ENDIF
          U=CX+GOLD*(CX-BX)
          FU=FUNC(U)
        ELSE IF((CX-U)*(U-ULIM).GT.0.)THEN
          FU=FUNC(U)
          IF(FU.LT.FC)THEN
            BX=CX
            CX=U
            U=CX+GOLD*(CX-BX)
            FB=FC
            FC=FU
            FU=FUNC(U)
          ENDIF
        ELSE IF((U-ULIM)*(ULIM-CX).GE.0.)THEN
          U=ULIM
          FU=FUNC(U)
        ELSE
          U=CX+GOLD*(CX-BX)
          FU=FUNC(U)
        ENDIF
        AX=BX
        BX=CX
        CX=U
        FA=FB
        FB=FC
        FC=FU
        GOTO 1
      ENDIF
      RETURN
      END

      REAL FUNCTION BRENT(AX,BX,CX,F,TOL,XMIN)
      INTEGER ITMAX
      REAL AX,BX,CX,TOL,XMIN,F,CGOLD,ZEPS
      EXTERNAL F
      PARAMETER (ITMAX=100,CGOLD=.3819660,ZEPS=1.0E-10)
      INTEGER ITER
      REAL A,B,D,E,ETEMP,FU,FV,FW,FX,P,Q,R,TOL1,TOL2,U,V,W,X,XM
      A=MIN(AX,CX)
      B=MAX(AX,CX)
      V=BX
      W=V
      X=V
      E=0.
      FX=F(X)
      FV=FX
      FW=FX
      DO 11 ITER=1,ITMAX
        XM=0.5*(A+B)
        TOL1=TOL*ABS(X)+ZEPS
        TOL2=2.*TOL1
        IF(ABS(X-XM).LE.(TOL2-.5*(B-A))) GOTO 3
        IF(ABS(E).GT.TOL1) THEN
          R=(X-W)*(FX-FV)
          Q=(X-V)*(FX-FW)
          P=(X-V)*Q-(X-W)*R
          Q=2.*(Q-R)
          IF(Q.GT.0.) P=-P
          Q=ABS(Q)
          ETEMP=E
          E=D
          IF(ABS(P).GE.ABS(.5*Q*ETEMP).OR.P.LE.Q*(A-X).OR.P.GE.Q*(B-X)) 
     *GOTO 1
          D=P/Q
          U=X+D
          IF(U-A.LT.TOL2 .OR. B-U.LT.TOL2) D=SIGN(TOL1,XM-X)
          GOTO 2
        ENDIF
1       IF(X.GE.XM) THEN
          E=A-X
        ELSE
          E=B-X
        ENDIF
        D=CGOLD*E
2       IF(ABS(D).GE.TOL1) THEN
          U=X+D
        ELSE
          U=X+SIGN(TOL1,D)
        ENDIF
        FU=F(U)
        IF(FU.LE.FX) THEN
          IF(U.GE.X) THEN
            A=X
          ELSE
            B=X
          ENDIF
          V=W
          FV=FW
          W=X
          FW=FX
          X=U
          FX=FU
        ELSE
          IF(U.LT.X) THEN
            A=U
          ELSE
            B=U
          ENDIF
          IF(FU.LE.FW .OR. W.EQ.X) THEN
            V=W
            FV=FW
            W=U
            FW=FU
          ELSE IF(FU.LE.FV .OR. V.EQ.X .OR. V.EQ.W) THEN
            V=U
            FV=FU
          ENDIF
        ENDIF
11    CONTINUE
      PAUSE 'BRENT EXCEED MAXIMUM ITERATIONS'
3     XMIN=X
      BRENT=FX
      RETURN
      END

      FUNCTION F1DIM(X)
      INTEGER NMAX
      REAL F1DIM,LOGMLFUNC,X
      PARAMETER (NMAX=50)
CU    USES FUNC
      INTEGER J,NCOM
      REAL PCOM(NMAX),XICOM(NMAX),XT(NMAX)
      COMMON /F1COM/ PCOM,XICOM,NCOM
      DO 11 J=1,NCOM
        XT(J)=PCOM(J)+X*XICOM(J)
11    CONTINUE
      F1DIM=LOGMLFUNC(XT,XT)
      RETURN
      END

C
C End of powell.f
C
