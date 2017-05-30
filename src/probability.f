C     ****h* snoc/probability.f [1.0]
C
C   NAME
C     Probability -- a set of statistics routines
C
C   DESCRIPTION
C     This package contains some distribution function and other
C     statistics related routines.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-03-03
C
C     ***



C     ****f* probability/GAUSS [1.0]
C
C
C   NAME
C     GAUSS -- the normal distributions function
C
C   DESCRIPTION
C     The normal distribution function.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-03-03
C
C   USAGE
C     GAUSS( X , MU , SIGMA, SERR)
C
C   INPUTS
C     MU    - the expectation value
C     X     - the value for which the pdf value is calculated
C     SIGMA - the standard devation
C     SERR  - the measurement error (if none set it to 0)
C
C   RESULTS
C     The value of the p.d.f is returned as DOUBLE PRECISION.
C
C     ***
      DOUBLE PRECISION FUNCTION GAUSS (M_T, BZMAG, SMAG, MERR)
      IMPLICIT NONE
      DOUBLE PRECISION M_T
      REAL BZMAG, SMAG, MERR

      DOUBLE PRECISION PI
      PARAMETER (PI = 3.14159265 )

      GAUSS = 1.0/(DSQRT(2.0*PI)*DBLE(SQRT(SMAG**2 + MERR**2)))*
     $     EXP(-0.5*(DBLE(BZMAG) - M_T)**2/DBLE(SMAG**2 + MERR**2))

      RETURN
      END
C
C End of GAUSS
C


C
C End of probability.f
C
