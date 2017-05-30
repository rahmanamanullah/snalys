C     ****f* snalys/LOGGAUSSLIN *
C
C   NAME
C     LOGGAUSSLIN -- Calculates the logarithm of the GUASSLIN distribution
C
C   DESCRIPTION
C     Calculates the logarithm of a Gaussian distribution with a 10*(m)
C     tail given the input.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-11-13
C
C   USAGE
C     LOGGAUSSLIN( CM , MM , V , B , M0 , ML  )
C
C   INPUTS
C     CM - The calculated magnitude of redshift and cosmology.
C     MM - The measured magnitude.
C     V  - The variance in magnitude.
C     B  - The linear factor.
C     M0 - The shift of the centre point for the Gaussian.
C     ML - The limit where the function can be approximated with
C          only the linear part assuming the expectation value of
C          of the distribution is zero.
C
C   RESULT
C     The logarithm of the distribution function is returned (with
C     exception for the cosmology independent normalization factor).
C
C   SEE ALSO
C     - Snalys Note I
C
C     ***
      REAL FUNCTION LOGGAUSSLIN( CM , MM , V , B , M0 , ML )
      IMPLICIT NONE
      REAL CM, MM, V, B, M0, ML
      
      REAL S, MC, PI, MDIFF
      PARAMETER ( PI = 3.14159265, S = 1.5 , MC = 0. )
      REAL LGELIN, LELIN

C     MDIFF = CM - MM
      MDIFF = MM - CM

C     WRITE (*,*) SQRT(V), V, B, M0

C     V = .16**2
C     B = .02
C     M0 = .03
C     FAC = SQRT(2.0*PI*V)


      IF ( MDIFF .GE. MC ) THEN
C        LOGGAUSSLIN = - .5*( MDIFF - M0 )**2/V - ALOG(FAC)
         LOGGAUSSLIN = - .5*( MDIFF - M0 )**2/V
      ELSE IF ( MDIFF .GE. ML ) THEN
         LOGGAUSSLIN = LGELIN( MDIFF, V, B, M0 ) 
C        LOGGAUSSLIN = ALOG( EXP( - .5*( MDIFF - M0)**2/V )/FAC +
C        LOGGAUSSLIN = ALOG( EXP( - .5*( MDIFF - M0)**2/V ) +
C    $        B*10.0**( S*( MDIFF ) ) )
C    $        B*ABS( MM - CM )*10.0**( S*( MM - CM ) ) )
C        IF ( LOGGAUSSLIN .GT. 0. ) WRITE (*,*) LOGGAUSSLIN,ABS(MM - CM)
      ELSE
         LOGGAUSSLIN = LELIN( MDIFF, V, B, M0 )
C        LOGGAUSSLIN = S*( MDIFF )*ALOG(10.) +
C    $        ALOG( B*ABS( MDIFF ) )
C    $        ALOG( B )
      END IF
      
C     IF ( LOGGAUSSLIN .GT. 0. ) WRITE (*,*) LOGGAUSSLIN, MDIFF, ML

      RETURN
      END
C
C End of LOGGAUSSLIN
C



C     ****f* probability/LGELIN
C
C
C   NAME
C     LGELIN -- the logarithm of the normal distribution function with a
C               linear tail.
C
C   DESCRIPTION
C     This function returns the logarithm of the Gaussian plus a 10^m
C     factor, assuming that the expectation value is zero. Also, setting
C     the normalization constant to one. Note that this is not a pdf since
C     the function goes to infinity as M increases.
C     
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-11-25
C
C   USAGE
C     LGELIN( M, V, B, M0 )
C
C   INPUTS
C     M  - the value for which the function will be calculated.
C     V  - the variance.
C     B  - the weight of the tail.
C     M0 - the displacement from zero of the centre point of the
C          function.
C
C   RESULTS
C     The value of the function in the point M is returned.
C
C     ***
      REAL FUNCTION LGELIN( M, V, B, M0 )
      IMPLICIT NONE
      REAL M, V, B, M0

      REAL S
      PARAMETER ( S = 1.5 )

C
C       version 1
C
C     LGELIN = ALOG( EXP( - .5*(M - M0)**2/V ) + B*10.0**( S*( M ) ) )

C
C       version 2
C
      LGELIN = ALOG( EXP( - .5*(M - M0)**2/V ) +
     $     ABS(M)*B*10.0**( S*( M ) ) )
      RETURN
      END
C
C  End of LGELIN
C


C     ****f* probability/LELIN
C
C
C   NAME
C     LELIN -- the logarithm of the linear tail.
C
C   DESCRIPTION
C     This function returns the logarithm of a 10^M factor, assuming that
C     the expectation value is zero. Also, setting the normalization
C     constant to one. Note that this is not a pdf since the function goes
C     to infinity as M increases.
C     
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-11-25
C
C   USAGE
C     LELIN( M, V, B, M0 )
C
C   INPUTS
C     M  - the value for which the function will be calculated.
C     V  - the variance.
C     B  - the weight of the tail.
C     M0 - the displacement from zero of the centre point of the
C          function.
C
C   RESULTS
C     The value of the function in the point M is returned.
C
C     ***
      REAL FUNCTION LELIN( M, V, B, M0 )
      IMPLICIT NONE
      REAL M, V, B, M0

      REAL S
      PARAMETER ( S = 1.5 )

C
C       version 1
C
C     LELIN = S*M*ALOG(10.) + ALOG(B)

C
C       version 2
C
      LELIN = S*M*ALOG(10.) + ALOG( B*ABS(M) )

      RETURN
      END
C
C  End of LELIN
C



C     ****f* probability/ELIN
C
C
C   NAME
C     ELIN -- the linear tail.
C
C   DESCRIPTION
C     This function returns the 10^M factor, assuming that the
C     expectation value is zero. Also, setting the normalization
C     constant to one. Note that this is not a pdf since the function
C     goes to infinity as M increases.
C     
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-11-25
C
C   USAGE
C     ELIN( M, V, B, M0 )
C
C   INPUTS
C     M  - the value for which the function will be calculated.
C     V  - the variance.
C     B  - the weight of the tail.
C     M0 - the displacement from zero of the centre point of the
C          function.
C
C   RESULTS
C     The value of the function in the point M is returned.
C
C     ***
      REAL FUNCTION ELIN( M, V, B, M0 )
      IMPLICIT NONE
      REAL M, V, B, M0

      REAL S
      PARAMETER ( S = 1.5 )

C
C       version 1
C
C     ELIN = B*10.**(S*M)

C
C       version 2
C
      ELIN = ABS(M)*B*10.**(S*M)
      RETURN
      END
C
C  End of ELIN
C



C     ****f* probability/DELIN
C
C
C   NAME
C     DELIN -- the linear tail.
C
C   DESCRIPTION
C     This function returns the derivative of the 10^M factor, assuming
C     that the expectation value is zero. Also, setting the
C     normalization constant to one. Note that this is not a pdf since
C     the function goes to infinity as M increases.
C     
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-11-25
C
C   USAGE
C     DELIN( M, V, B, M0 )
C
C   INPUTS
C     M  - the value for which the function will be calculated.
C     V  - the variance.
C     B  - the weight of the tail.
C     M0 - the displacement from zero of the centre point of the
C          function.
C
C   RESULTS
C     The derivative of the function in the point M is returned.
C
C     ***
      REAL FUNCTION DELIN( M, V, B, M0 )
      IMPLICIT NONE
      REAL M, V, B, M0

      REAL S
      PARAMETER ( S = 1.5 )

C
C       version 1
C
C     DELIN = S*ALOG(10.)*ELIN( M, V, B, M0 )

C
C       version 2
C
      DELIN = B*(1 + ABS(M)*S*ALOG(10.))*10.**(S*M)
      RETURN
      END
C
C  End of DELIN
C



C     ****f* snalys/DLELIN *
C
C   NAME
C     DLELIN -- The derivative of LELIN.
C
C   DESCRIPTION
C     The derivative of LLIN wit respect to the magnitude.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-11-25
C
C   USAGE
C     DLELIN( M, V, B, M0 )
C
C   INPUTS
C     M  - the value for which the function will be calculated.
C     V  - the variance.
C     B  - the weight of the tail.
C     M0 - the displacement from zero of the centre point of the
C          function.
C
C   RESULTS
C     The value of the derivative of the function in the point M
C     is returned.
C
C     ***
      REAL FUNCTION DLELIN( M, V, B, M0 )
      IMPLICIT NONE
      REAL M, V, B, M0

      REAL S
      PARAMETER ( S = 1.5 )

C 
C       version 1
C
C     DLELIN = S*ALOG(10.)

C
C       version 2
C
      DLELIN = 1./ABS(M) + S*ALOG(10.)
      RETURN
      END
C
C End of DLELIN
C


C
C End of gausslin.f
C
