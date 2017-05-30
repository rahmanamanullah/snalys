c The davidon minimization routine, stolen from
c Christian Walck, walck@physto.se



C     ****f* snalys/DAVIDO *
C
C   NAME
C     DAVIDO -- minimization by Davidon's variance algorithm (gradient
C               method)
C
C   DESCRIPTION
C     Minimization by Davidon's variance algorithm (gradient method)
C     For quadratic functions the method converges within NPAR
C     iterations to the exact minimum and variance matrix.
C
C   AUTHOR
C     Christian Walck (walck@physto.se)
C
C   CREATION DATE
C     1980-03-04
C
C   USAGE
C     DAVIDO ( NPAR , X , V , FGUSER , EPS )
C
C   INPUTS
C     NPAR   - Number of parameters (see note below)
C     X      - Parameter vector containing start values at
C              entry and position of minimum at return.
C     V      - Variance matrix of estimates (see note below)
C     FGUSER - Function which calculates function value and 
C              its gradients (FUNCTION FGUSER(X,G)).
C              Declare external in calling routine!
C     EPS    - Convergence limit (default=.000001) on twice
C              the estimated excess of the function above
C              its minimum.
C
C   RESULT
C     A LOGICAL variable is returned confiriming that the minimization
C     was successful at the required precision.
C
C   SIDE EFFECTS
C     X - is set to the point where the minimum was found
C     V - contains the variance matrix on return.
C
C   USES
C     FGUSER (user supplied)
C     F110 MXMPY, MXUTY, F121 VDOT, VSUB (GENLIB)
C
C   NOTE
C     The method uses two constants, ALPHA and BETA, specifying bounds
C     on the allowable change that may be made in the variance estimate
C     within one iteration. They should satisfy 0 < ALPHA < 1 < BETA and
C     is set by default to ALPHA=.001 and BETA=10. These constants as
C     well as the maximum number of iterations (by default 10*NPAR) may
C     be altered by including
C                COMMON /CWDAV/ ALPHA , BETA , NITER
C     in the calling program. ALPHA and BETA closer to one will most
C     likely slow down convergence but make the search more cautious.
C
C     Notations: X, XS, G, GS & R are row vectors of dimension NPAR
C                V is a symmetric matrix of dimension NPAR,NPAR
C
C     On entry V could be any positive definite symmetric matrix.
C     Normally a diagonal matrix with overestimated variances of
C     the parameters gives better convergence.
C     
C     NPAR should be given negative if V is initialized, if not
C     V is put to the unity matrix at entry.
C     
C     V is defined as the inverse of the matrix of second deriva-
C     tives, G-1. Note that if a least square fit is performed
C     V = 2 * G-1 and thus, in this case, all elements of V should
C     be multiplied by 2 after the fit.
C     
C     For further details see original paper. 
C
C   SEE ALSO
C     C.    "Variance Algorithm for Minimization" by William C. Davidon
C     [COMP.J. 10 (1968) 406]
C 
C     "Probability and Statistics in Particle Physics" by
C     A.G.Frodesen, O.Skjeggestad and H.Tofte, Universitetsforlaget,
C     Bergen-Oslo-Tromso, pp 363-365
C
C     ***
      LOGICAL FUNCTION DAVIDO ( NPAR , X , V , FGUSER , EPS )
      COMMON /CWDAV/ ALPHA , BETA , NITER
      DIMENSION X(*) , V(*) , XS(50) , G(50) , GS(50) , R(50)
C     DATA ALPHA /0.001/ , BETA /10.0/ , NITER /0/
C      DATA ALPHA /0.5/ , BETA /2.0/ , NITER /200/
      DATA ALPHA /.75/ , BETA /1.5/ , NITER /200/


C  
C      For backward compatibility include entry davidon for CDC and VAX
C  
      ENTRY DAVIDON
C  
C        (0)      Set constants and calculate F and G
C                 Set H = I if not initialized (=> 1st direction along
C                 line of steepest descent).
C  

      DAVIDO = .TRUE.

      NPR = IABS ( NPAR )
      X1 = - BETA / ( BETA - 1.0 )
      X2 = - BETA / ( BETA + 1.0 )
      X3 = - ALPHA / ( 1.0 + ALPHA )
      X4 = ALPHA / ( 1.0 - ALPHA )
      EPSI = 0.000001
      IF (EPS .GT. 0.0) EPSI = EPS
      F = FGUSER ( X , G )
      IF ( NPAR .GE. 0 ) THEN
         DO 50 I = 1, NPR
   50    V((I-1)*NPR+I) = 1.0
      END IF
      MAX = 50 * NPR
      IF ( NITER .GT. 0 ) MAX = NITER
C
      ITER = 0
  100 IF ( ITER .GE. MAX ) THEN
         DAVIDO = .FALSE.
         GO TO 103
      ENDIF
      ITER = ITER + 1

C  
C        (I)      Set XS = X - G.V
C  
      DO 120 I = 1, NPR
         XS(I) = X(I)
         DO 110 J = 1, NPR
  110    XS(I) = XS(I) - V((I-1)*NPR+J) * G(J)
  120 CONTINUE
C  
C        (II)     Calculate F and G at XS
C  
      FS = FGUSER ( XS , GS )
C  
C        (III)    Set R = GS.V and RHO = GS.R'
C                 if RHO < EPSILON go to (VII)
C  
      DO 140 I = 1, NPR
         R(I) = 0.0
         DO 130 J = 1, NPR
  130    R(I) = R(I) + V((I-1)*NPR+J) * GS(J)
  140 CONTINUE
      RHO = 0.0
      DO 150 I = 1, NPR
  150 RHO = RHO + GS(I) * R(I)
      IF ( RHO .LT. EPSI ) GO TO 103
C  
C        (IV)     Set GAMMA = - G.R' / RHO and determine LAMBDA
C  
      GAMMA = 0.0
      DO 160 I = 1, NPR
  160 GAMMA = GAMMA - G(I) * R(I)
      GAMMA = GAMMA / RHO
      IF ( GAMMA.LT.X1 .OR.  GAMMA.GE.X4 ) XLAM = GAMMA / (GAMMA + 1.0)
      IF ( X3.LE.GAMMA .AND. GAMMA.LT.X4)  XLAM = ALPHA
      IF ( X2.LE.GAMMA .AND. GAMMA.LT.X3)  XLAM = - GAMMA/(GAMMA + 1.0)
      IF ( X1.LE.GAMMA .AND. GAMMA.LT.X2)  XLAM = BETA
C
C        (V)      Update V by V = V + ( LAMBDA - 1 ) * R'.R / RHO
C
      DO 101 I = 1 , NPR
      DO 101 J = 1 , NPR
  101 V(NPR*(I-1)+J) = V(NPR*(I-1)+J) + (XLAM-1.0) * R(I) * R(J) / RHO
C  
C        (VI)     If FS >= F go to (I) otherwise set X = XS, F = FS,
C                 G = GS and go to (I)
C  
      IF ( FS .GE. F ) GO TO 100
      F = FS
      DO 102 I = 1 , NPR
      X(I) = XS(I)
  102 G(I) = GS(I)
      GO TO 100
C  
C        (VII)    Minimum found or maximum number of iterations exceeded
C                 return with parameters in X and variance matrix in V.
C                 Modified 940423 to replace parameters only if better minimum.
C  
  103 IF ( FS .LT. F ) THEN
         DO 104 I = 1 , NPR
  104    X(I) = XS(I)
      END IF
      RETURN
      END

C  Double precision version

      LOGICAL FUNCTION DDAVID ( NPAR , X , V , FGUSER , EPS )
C
C     Double precision version of the Davidon minimization routine DAVIDO 
C     (see comments for this routine).
C
C                                                        /Ch.Walck 911211
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE /DCWDAV/
      COMMON /DCWDAV/ ALPHA , BETA , NITER
      DIMENSION X(*) , V(*) , XS(50) , G(50) , GS(50) , R(50)
      DATA ALPHA /0.001D0/ , BETA /10.0D0/ , NITER /0/
C
C        (0)      Set constants and calculate F and G
C                 Set V = I if not initialized (=> 1st direction along
C                 line of steepest descent).
C
      DDAVID = .TRUE.

      NPR  = IABS ( NPAR )
      X1   = -  BETA / ( BETA  - 1.0D0 )
      X2   = -  BETA / ( BETA  + 1.0D0 )
      X3   = - ALPHA / ( 1.0D0 + ALPHA )
      X4   =   ALPHA / ( 1.0D0 - ALPHA )
      EPSI = 1.0D-9
      IF ( EPS .GT. 0.0D0 ) EPSI = EPS
      F    = FGUSER ( X , G )
      IF ( NPAR .GE. 0 ) THEN
         DO 50 I = 1, NPR
   50    V((I-1)*NPR+I) = 1.0D0
      END IF
C     MAX = 10 * NPR
      MAX = 25 * NPR
      IF ( NITER .GT. 0 ) MAX = NITER
C
      ITER = 0
  100 IF ( ITER .GE. MAX ) THEN
         DDAVID = .FALSE.
         GO TO 103
      ENDIF
      ITER = ITER + 1
C
C       (I)      Set XS = X - G.V
C
      DO 120 I = 1, NPR
         XS(I) = X(I)
         DO 110 J = 1, NPR
  110    XS(I) = XS(I) - V((I-1)*NPR+J) * G(J)
  120 CONTINUE
C
C        (II)     Calculate F and G at XS
C
      FS = FGUSER ( XS , GS )
C
C        (III)    Set R = GS.V and RHO = GS.R'
C                 if RHO < EPSILON go to (VII)
C
      DO 140 I = 1, NPR
         R(I) = 0.0D0
         DO 130 J = 1, NPR
  130    R(I) = R(I) + V((I-1)*NPR+J) * GS(J)
  140 CONTINUE
      RHO = 0.0D0
      DO 150 I = 1, NPR
  150 RHO = RHO + GS(I) * R(I)
      IF ( RHO .LT. EPSI ) GO TO 103
C
C        (IV)     Set GAMMA = - G.R' / RHO and determine LAMBDA
C
      GAMMA = 0.0D0
      DO 160 I = 1, NPR
  160 GAMMA = GAMMA - G(I) * R(I)
      GAMMA = GAMMA / RHO
      HLP = GAMMA / ( GAMMA + 1.0D0 )
      IF ( GAMMA.LT.X1 .OR.  GAMMA.GE.X4 ) XLAM = HLP
      IF ( X3.LE.GAMMA .AND. GAMMA.LT.X4 ) XLAM = ALPHA
      IF ( X2.LE.GAMMA .AND. GAMMA.LT.X3 ) XLAM = - HLP
      IF ( X1.LE.GAMMA .AND. GAMMA.LT.X2 ) XLAM = BETA
C
C        (V)      Update V by V = V + ( LAMBDA - 1 ) * R'.R / RHO
C
      DO 101 I = 1 , NPR
      DO 101 J = 1 , NPR
  101 V(NPR*(I-1)+J) = V(NPR*(I-1)+J) + (XLAM-1.0D0) * R(I) * R(J) / RHO
C
C        (VI)     If FS >= F go to (I) otherwise set X = XS, F = FS,
C                 G = GS and go to (I)
C
      IF ( FS .GE. F ) GO TO 100
      F = FS
      DO 102 I = 1 , NPR
      X(I) = XS(I)
  102 G(I) = GS(I)
      GO TO 100
C
C        (VII)    Minimum found or maximum number of iterations exceeded
C                 return with parameters in X and variance matrix in V.
C                 Modified 940423 to replace parameters only if better minimum.
C
  103 IF ( FS .LT. F ) THEN
         DO 104 I = 1 , NPR
  104    X(I) = XS(I)
      END IF
      RETURN
      END
C
C End of davidon.f
C
