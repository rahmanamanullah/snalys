C     ****f* solve/RTNEWT *
C
C   NAME
C     RTNEWT - root finding with Newton-Raphson's method
C
C   DESCRIPTION
C     Simple one dimensional root finding in a brackted interval
C     using Newton-Raphsons method.
C
C   USAGE
C     CALL RTNEWT(FUNCD, X1, X2, XACC)
C
C   INPUTS
C     FUNCD - to be solved, takes three argument and returns
C             In the second argument the function value is returned
C             and in the third the derative of the function.
C     X1    - lower X limit
C     X2    - higher X limit
C     XACC  - requested accuracy.
C
C   RESULT
C     The root is returned.
C
C   SEE ALSO
C     Numerical Recepies
C
C     ***
      FUNCTION rtnewt(funcd,x1,x2,xacc)
      INTEGER JMAX
      REAL rtnewt,x1,x2,xacc
      EXTERNAL funcd
      PARAMETER (JMAX=20)
      INTEGER j
      REAL df,dx,f
      rtnewt=.5*(x1+x2)
      do 11 j=1,JMAX
        call funcd(rtnewt,f,df)
        dx=f/df
        rtnewt=rtnewt-dx
        if((x1-rtnewt)*(rtnewt-x2).lt.0.)pause
     *'rtnewt jumped out of brackets'
        if(abs(dx).lt.xacc) return
11    continue
      pause 'rtnewt exceeded maximum iterations'
      END
C
C End of RTNEWT
C
