C   cosmo.f
C
C
C   rahman@physto.se, 001025
C


C     ****f* snalys/RANDCOSMO *
C
C   NAME
C     RANDCOSMO() -- Generates a randomized cosmology.
C
C   DESCRIPTION
C     Uses the parameter constraints defined in the ini-file and
C     uniformly randomizes the cosmological parameters between
C     those values.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-10-25
C
C   USAGE
C     RANDCOSMO( P )
C
C   INPUTS
C     P - array with the parameters.
C
C   SIDE EFFECTS
C     The randomized cosmology is stored in P.
C
C     ***
      SUBROUTINE RANDCOSMO(P)
      IMPLICIT NONE
      REAL P(*)

      INCLUDE 'cosmo.inc'
      INCLUDE 'operation.inc'
      
      REAL C_END_GET, C_START_GET, RNUNIF, BIGBANG
      INTEGER N

C
C        Uniformly randomizes the cosmological variables between
C        the constraints that are given in the ini-file.
C
      DO N = 1, NR_VAR
         P(N) = (C_END_GET(N) - C_START_GET(N))*RNUNIF(NR_VAR) +
     $        C_START_GET(N)
      ENDDO
C
C        Make sure that we do not start in the non-big bang area.
C
      IF ( .NOT.OPER(FLAT) .AND. O_M.GT.0 .AND. O_X.GT.0 ) THEN
         DO WHILE ( P(O_X).GE.BIGBANG(P(O_M)) .OR.
     $        P(O_X).LE.0.0 .OR. P(O_M).LT.0.0 )
            P(O_M) = RNUNIF(NR_VAR)*(C_END_GET(O_M) -
     $           C_START_GET(O_M)) + C_START_GET(O_M)
            P(O_X) = RNUNIF(NR_VAR)*(C_END_GET(O_X) -
     $           C_START_GET(O_X)) + C_START_GET(O_X)
         ENDDO
         P(O_X) = RNUNIF(NR_VAR)
         P(O_M) = RNUNIF(NR_VAR)
      ENDIF
      END
C
C     End of RANDCOSMO().
C      

C     ****f* snalys/BIGBANG *
C
C   NAME
C     BIGBANG() -- Gives the limiting value in Omega_Lambda
C                  for a Big Bang to be possible.
C
C   DESCRIPTION
C     For a given value of Omega_M the value of Omega_X that corresponds
C     to the limit where a Big Bang model is allowed is calculated.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-03-29
C
C   INPUTS
C     Omega_M - A REAL value, the mass density of the universe.
C
C   RESULT
C     Omega_X at the Big Bang limit is returned.
C
C   SEE ALSO
C     The Cosmological Constant, Sean M. Carroll and William H. Press
C     Annu. Rev. Astron. Astrophys. 1992.30:499-542
C
C     ***
      REAL FUNCTION BIGBANG( OMEGA_M )
      IMPLICIT NONE
      REAL OMEGA_M

C
C       The limit of the function as O_M -> 0.0.
C
      IF ( OMEGA_M.EQ.0.0 ) THEN
         BIGBANG = 1.0
      ELSE IF ( OMEGA_M.LT.0.5D0 ) THEN
         BIGBANG = 4*OMEGA_M*( COSH(LOG( ((1.0 - OMEGA_M)/OMEGA_M) + 
     $        SQRT(((1.0 - OMEGA_M)/OMEGA_M)**2 - 1 ))
     $        /3.0) )**3
      ELSE IF ( OMEGA_M.GT.0.5D0 ) THEN
         BIGBANG = 4*OMEGA_M*( COS(ACOS((1.0 -
     $        OMEGA_M)/OMEGA_M)/3.0) )**3
      ELSE
         BIGBANG = 4*OMEGA_M
      ENDIF

      RETURN
      END
C
C     End of BIGBANG()
C

C
C     End of file cosmo.f
C
