C     ****i* contourfinder/PYTHAG *
C
C   NAME
C     PYTHAG() -- computes (a**2 + b**2)**0.5
C
C   DESCRIPTION
C     Computes (a**2 + b**2)**0.5 without destructive underflow or
C     overflow.
C
C   AUTHOR
C     William H. Press, Saul A. Teukolsky, William T. Vetterling and
C     Brian P. Flannery
C
C   USAGE
C     PYTHAG( A, B )
C
C   INPUTS
C     A - a REAL value.
C     B - a REAL value.
C
C   RESULT
C     (A**2 + B**2)**0.5 is returned.
C
C   SEE ALSO
C     Numerical Recepies - The Art of Scientific Computing, Second Edition,
C     by William H. Press, Saul A. Teukolsky, William T. Vetterling and
C     Brian P. Flannery, section 2.6 Singular Value Decomposition.
C
C   USED BY
C     SVDCMP()
C
C     ***
      FUNCTION pythag(a,b)
      REAL A,B,PYTHAG
      REAL ABSA,ABSB
      ABSA=ABS(A)
      ABSB=ABS(B)
      IF(ABSA.GT.ABSB)THEN
        PYTHAG=ABSA*SQRT(1.+(ABSB/ABSA)**2)
      ELSE
        IF(ABSB.EQ.0.)THEN
          PYTHAG=0.
        ELSE
          PYTHAG=ABSB*SQRT(1.+(ABSA/ABSB)**2)
        ENDIF
      ENDIF
      RETURN
      END
C     *** End of FUNCTION PYTHAG() ***

C     END OF FILE
