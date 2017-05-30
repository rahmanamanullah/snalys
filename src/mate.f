C     ****f* snalys/MATE *
C
C   NAME
C     MATE -- The dark term in F_INTEGRAND in the case of 
C             inverse power-law quintessence.      
C
C   DESCRIPTION
C     Uses linear interpolation to calculate the dark term from the
C     datafile 'wcube.dat' for any combination of alpha, Omm and z 
C     in the ranges
C      
C     1.0 .LE. alpha .LE. 10.0
C     0.1 .LE.  Omm  .LE.  0.6
C     0.0 .LE.   z   .LT.  2.0
C
C     (values very close to z = 2 might cause problems).      
C
C   AUTHOR
C     Martin Eriksson (mate@physto.se)
C
C   CREATION DATE
C     2001-09-07
C
C   INPUTS
C     ZP    - is a redshift of type DOUBLE PRECISION.
C
C   RESULT
C     The dark term is returned.
C
C     ***      
      DOUBLE PRECISION FUNCTION MATE( ZP )
      IMPLICIT NONE
      DOUBLE PRECISION ZP

      INCLUDE 'cosmo.inc'
      INCLUDE 'mate.inc'
      
      INTEGER N_VAR
      PARAMETER (N_VAR = 3)
      INTEGER NA(N_VAR)
      DOUBLE PRECISION X(N_VAR)
      DOUBLE PRECISION FINT

      DATA NA /GLOB_N_ALPHA, GLOB_N_OMM, GLOB_N_ZP/

C
C        Initialize variable array X
C
      X(1) = GLOB_ALPHA
      X(2) = GLOB_OMM
      X(3) = ZP
      
      MATE = FINT(N_VAR, X, NA, GLOB_WINDEX, GLOB_WCUBE)

      RETURN
      END
C      
C End of MATE
C      


C        ****f* snalys/FINT *
C
C   NAME
C     FINT -- linear interpolation routine
C
C   DESCRIPTION
C     Linear interpolation.
C
C   AUTHOR
C     C. LETERTRE, modified by B. SCHORR
C
C   CREATION DATE
C     1982-07-01
C
C   INPUTS
C
C   RESULT
C
C   ERRORS
C     Error message is written if the number of variables is
C     more than five.
C
C   HISTORY
C     2001-09-04 Slightly modified by Martin Eriksson.
C
C     ***      
      DOUBLE PRECISION FUNCTION FINT(NARG, ARG, NENT, ENT, TABLE)
      IMPLICIT NONE
      INTEGER   NENT(9), NARG
      DOUBLE PRECISION ARG(9), ENT(9), TABLE(9)

      INTEGER   I, INDEX(32), ISTEP, ISHIFT, K, KNOTS, LMAX, LMIN, 
     &     LOCA, LOCB, LOCC, N, NDIM
      DOUBLE PRECISION WEIGHT(32), X, H, ETA

      FINT  =  0.
      IF (NARG .LT. 1  .OR.  NARG .GT. 5)  GOTO 300
      LMAX      =  0
      ISTEP     =  1
      KNOTS     =  1
      INDEX(1)  =  1
      WEIGHT(1) =  1.
      DO 100    N  =  1, NARG
         X     =  ARG(N)
         NDIM  =  NENT(N)
         LOCA  =  LMAX
         LMIN  =  LMAX + 1
         LMAX  =  LMAX + NDIM
         IF (NDIM .GT. 2)  GOTO 10
         IF (NDIM .EQ. 1)  GOTO 100
         H  =  X - ENT(LMIN)
         IF (H .EQ. 0.)  GOTO 90
         ISHIFT  =  ISTEP
         IF (X-ENT(LMIN+1) .EQ. 0.)  GOTO 21
         ISHIFT  =  0
         ETA     =  H / (ENT(LMIN+1) - ENT(LMIN))
         GOTO 30
 10      LOCB  =  LMAX + 1
 11      LOCC  =  (LOCA+LOCB) / 2
         IF (X-ENT(LOCC))  12, 20, 13
 12      LOCB  =  LOCC
         GOTO 14
 13      LOCA  =  LOCC
 14      IF (LOCB-LOCA .GT. 1)  GOTO 11
         LOCA    =  MIN0( MAX0(LOCA,LMIN), LMAX-1 )
         ISHIFT  =  (LOCA - LMIN) * ISTEP
         ETA     =  (X - ENT(LOCA)) / (ENT(LOCA+1) - ENT(LOCA))
         GOTO 30
 20      ISHIFT  =  (LOCC - LMIN) * ISTEP
 21      DO 22  K  =  1, KNOTS
            INDEX(K)  =  INDEX(K) + ISHIFT
 22      CONTINUE
         GOTO 90
 30      DO 31  K  =  1, KNOTS
            INDEX(K)         =  INDEX(K) + ISHIFT
            INDEX(K+KNOTS)   =  INDEX(K) + ISTEP
            WEIGHT(K+KNOTS)  =  WEIGHT(K) * ETA
            WEIGHT(K)        =  WEIGHT(K) - WEIGHT(K+KNOTS)
 31      CONTINUE
         KNOTS  =  2*KNOTS
 90      ISTEP  =  ISTEP * NDIM
 100  CONTINUE
      DO 200    K  =  1, KNOTS
         I  =  INDEX(K)
         FINT  =  FINT + WEIGHT(K) * TABLE(I)
 200  CONTINUE
      RETURN
 300  WRITE(*,*) 'ERROR: Number of variables must be in the range 1-5!'

      RETURN
      END
C
C End of FINT
C      

C
C End of mate.f
C




