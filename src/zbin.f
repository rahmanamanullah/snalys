C     ****f* snalys/ZBINFUNC *
C
C   NAME
C     ZBIN -- calculate the chi2 sum for a given cosmology
C
C   DESCRIPTION
C     Calculates the chi2 sum for a given cosmology using the definitions
C     of the red shift bins defined in the include file.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-08-10
C
C   USAGE
C     ZBINFUNC(C)
C
C   INPUTS
C     C - is a vector containing the cosmology.
C
C   RESULT
C     The chi2 sum is returned.
C
C     ***
      REAL FUNCTION ZBINFUNC(C)
      IMPLICIT NONE
      REAL C(*)

      INCLUDE 'cosmo.inc'
      INCLUDE 'zbin.inc'
      LOGICAL UPDATE_GLOBALS
      DOUBLE PRECISION CHI2FUNC
      REAL CALCDL, CALCM
      EXTERNAL CALCDL, CALCM
      REAL EST(MAX_BIN)
      INTEGER N

C  
C        Default return value if something goes wrong.
C  
      ZBINFUNC = 10000.0

      IF ( UPDATE_GLOBALS(C) ) THEN
         DO N = 1, NR_BIN
            EST(N) = REAL( CALCM( CALCDL(0., ZBINS(N)), C(M_SC) ) )
         ENDDO

         ZBINFUNC = REAL(CHI2FUNC( EST, MEASMAG, NR_BIN, SIGMA_MAG ))
      ENDIF
      RETURN
      END
C
C End of ZBINFUNC
C



C     ****f* snalys/CHI2FUNC *
C
C   NAME
C     CHI2FUNC -- calculate the chi2 sum.
C
C   DESCRIPTION
C     Calculates the chi2 sum for a set of measurements and
C     expected values.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-08-10
C
C   USAGE
C     CHI2FUNC(EST,MEAS,NRVAL,SIGMA)
C
C   INPUTS
C     EST   - Vector containing the estimated values
C     MEAS  - Vector containing the measured values.
C     NRVAL - INTEGER defining the number of values.
C     SIGMA - The standard deviation.
C
C   RESULT
C     The chi2 sum is returned.
C
C     ***
      DOUBLE PRECISION FUNCTION CHI2FUNC(EST,MEAS,NRVAL,SIGMA)
      IMPLICIT NONE
      REAL EST(*), MEAS(*), SIGMA
      INTEGER NRVAL

      INTEGER N

      CHI2FUNC = 0.0
      
      DO N = 1, NRVAL
         CHI2FUNC = CHI2FUNC + ( EST(N) - MEAS(N) )**2
      ENDDO
      CHI2FUNC = CHI2FUNC/SIGMA**2
      
      RETURN
      END
C
C End of CHI2FUNC
C

C
C End of zbin.f
C
