C   ml_func.f
C
C   In this file a the negative likelihood function is defined, together
C   with its derivatives... enjoy!
C
C   rahman@physto.se, 000329
C
   
C     ****f* snalys/LOGMLFUNC *
C
C   NAME
C     LOGMLFUNC() -- Calculates the negative logarithm of the ML-function
C                      for a specific cosmology.
C
C   DESCRIPTION
C     The negative logarithm of the ML-function for a
C     specific cosmology is calculated, together with the gradient of
C     the ML-function in that specific point.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-03-29
C
C   USAGE
C     LOGMLFUNC( P, G )
C
C   INPUTS
C     P - is a REAL array of the cosmological parameters that are to
C         be estimated.
C     G - is a REAL array.
C
C   RESULT
C     The value logarithm of the ML-function will be returned.
C
C   SIDE EFFECTS
C     G - is a REAL vector that contains the calculated gradient
C         if N >= 0 in the input.
C
C   USES
C     PGARRO of the PGPLOT library (but this should be moved to
C     some other function).
C
C   TODO
C     - Move the part where the cosmology is changed incase it is
C       non-physical.
C     - The code that takes care of non-physical regions (the
C       non-Big Bang region for instance should be rewritten
C       in a much fancier way).
C
C   HISTORY
C     2001-09-08 Removed the possiblity to draw the Davidon convergence
C                path.
C     2001-09-27 The magnitudes and luminosity distances are all
C                calculated at the same time for each cosmology.
C                This procedure makes the routine about 3% faster.
C  2004-10-27 Added the saving of the values of the nuisance parameters
C                for the minimum chi2 for this parameter combinations, if
C                FITMSC or FITSALP have been set.
C 
C
C     ***
      REAL FUNCTION LOGMLFUNC( P, D )
      IMPLICIT NONE
      REAL P(*), D(*)

      INCLUDE 'snalys.inc'
      INCLUDE 'cosmo.inc'
      INCLUDE 'operation.inc'
      INCLUDE 'supernova.inc'
      INCLUDE 'probability.inc'
      INCLUDE 'error.inc'
      INCLUDE 'precalcs.inc'

      DATA CALCDERIV/.FALSE./

      INTEGER I, N, NMSC, NSALP, STATUS
      LOGICAL UPDATE_GLOBALS
      REAL TOL , Q , XMIN, TMP
      REAL NEGMLSM, NEGMLSALP, NEGMLSMSA, RNUNIF
      REAL C_START_GET, C_STEP_GET, C_NR_GET
      LOGICAL DAVIDO
      EXTERNAL NEGMLSM, NEGMLSALP, NEGMLSMSA, RNUNIF, DAVIDO
      PARAMETER ( TOL = 1.E-7 )
      INTEGER FITNR, NRVAR, GADZ
      PARAMETER ( NRVAR = 2 , GADZ = NRVAR**2 )
      REAL PS(NRVAR) , GS(NRVAR), VAR_MAT(GADZ)


      REAL LOGDIST
      EXTERNAL LOGDIST

      REAL CALCA,BOS_A
      REAL CALCCMB,CMB_R

C
C        Update the global cosmological variables that will be
C        used for the calculations. This function will fail if
C        the parameters have non-physical values.
C
      IF ( UPDATE_GLOBALS(P) ) THEN
C     
C           Calculate the magnitudes and luminosity distances for
C           each supernova.
C     
         FITNR = NRVAR
         IF ( OPER(FITMSC) .OR. OPER(FITSALP) ) THEN
            CALL FOREACHSN( SNZS, SNIX, SNM, SNDL, -SN_NR, SNMS, SNSS,
     $           SNVS, Q, D, -FITNR, STATUS )
            Q = 1.E16

            IF ( STATUS .EQ. 1 ) THEN

C
C                 Check if the nuisence parameters are marked for fitting.
C
               IF ( M_SC .EQ. 0 ) THEN
                  NMSC = 1
               ELSE
                  NMSC = INT(C_NR_GET(M_SC))
               ENDIF
               
               IF ( SALP .EQ. 0 ) THEN
                  NSALP = 1
               ELSE
                  NSALP = INT(C_NR_GET(SALP))
               ENDIF
               
               DO I = 1, NMSC
                  DO N = 1, NSALP
                     IF ( M_SC .EQ. 0 ) THEN
                        PS(1) = GLOB_MSCRIPT
                     ELSE
                        PS(1) = C_START_GET(M_SC) + 
     $                       C_STEP_GET(M_SC)*(I-1.)
                     ENDIF
                     IF ( SALP .EQ. 0 ) THEN
                        PS(2) = GLOB_SALP
                     ELSE
                        PS(2) = C_START_GET(SALP) +
     $                       C_STEP_GET(SALP)*(N-1.)
                     ENDIF

                     TMP = NEGMLSMSA( PS, GS )
                     
                     IF ( TMP .LT. Q ) THEN
                        THIS_SCRIPTM = PS(1)
                        THIS_ALPHA = PS(2)
                        Q = TMP
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF
         ELSE
            IF ( .NOT. CALCDERIV ) FITNR = -NR_VAR
            CALL FOREACHSN( SNZS, SNIX, SNM, SNDL, SN_NR, SNMS, SNSS,
     $           SNVS, Q, D, FITNR, STATUS )


C           DO I = 1, SN_NR
C              WRITE (*,*) LOGDIST( SNM(I), SNMS(I), SNVS(I), 0., 0.,
C    $              0., 1),
C    $              EXP(LOGDIST( SNM(I), SNMS(I), SNVS(I), 0., 0., 0.,
C    $              1)), SNM(I), SNMS(I), SNVS(I)
C
C              WRITE(*,*) SNZS(I),SNM(I),SNMS(I),SNVS(I)
C           ENDDO
         ENDIF

C
C          Set Q to an overwhelmingly large value if FOREACHSN
C          failed.
C     

         IF ( STATUS .NE. 1 ) THEN
            Q = 1.E16
         ENDIF
C     Calculate the A-parameter measured in Baryon Oscillation Measurements
C     for this cosmology
         IF ( OPER(DO_BOS) ) THEN
            BOS_A = CALCA()
            Q = Q + .5*(BOS_A - BOSPRI)**2/DBOSPR**2
c            write(*,*)  bos_a,bospri,q
         ENDIF
C     Calculate the R-parameter measured by CMB
C     for this cosmology
         IF ( OPER(DO_CMB) ) THEN
            CMB_R = CALCCMB()
            Q = Q + .5*(CMB_R - CMBPRI)**2/DCMBPR**2
C            write(*,*)  cmb_r,cmbpri,dcmbpr,q
         ENDIF
C
C     Calculate probability for given OMEGA_M W.R.T PRIOR
      IF ( OPER(DO_OM) ) THEN
        Q = Q + .5*(GLOB_OMM - OMPRI)**2/DOMPR**2
      ENDIF
      ELSE

C        
C
C            If the function takes a non-physical value, reseting
C            these parameters should do the job to bring us back.
C
C            This should be rewritten.
C


         Q = 0.0
         DO I = 1, NR_VAR
            D(I) = 0.0

            IF ( I .EQ. O_M ) THEN
               Q = Q + ( P(I) - .5 )**2
               D(I) = 2. * ( P(I) - .5 )
            ELSE IF ( I .EQ. M_SC ) THEN
               Q = Q + ( P(I) + 3.4 )**2
               D(I) = 2. * ( P(I) + 3.4 )
            ELSE
               Q = Q + P(I)**2
               D(I) = 2. * P(I)
            END IF

            D(I) = 1. * D(I)

         END DO

         Q = 1. * Q  + REAL(SN_NR) * 100.0

C
C     OLD STUFF
C


C                  Q = REAL(SN_NR)*100.0
C
C         IF ( GLOB_OMM .LT. 0.0 ) THEN
C            D(O_M) = -2.0*100.0**ABS(REAL(GLOB_OMM))
C            Q = Q + 100.0*ABS(REAL(GLOB_OMM))**2
C
C         ELSE IF ( GLOB_OMX .LT. -1.5 ) THEN
C            D(O_X) = -2.0*100.0*ABS(REAL(GLOB_OMX))
C            Q = Q + 100.0*ABS(REAL(GLOB_OMX))**2
C         ELSE
C            D(O_X) = 2.0*100.0*ABS(REAL(GLOB_OMX))
C            Q = Q + 100.0*ABS(REAL(GLOB_OMX))**2
C         END IF

C
C         EVEN OLDER STUFF
C
C         DO I = 1, NR_VAR
C            Q  = Q + P(I)**4
C         ENDDO
C         Q = Q*100.0*REAL(SN_NR)
C         TMP = Q**2
C         
C         DO I = 1, NR_VAR
C            IF ( I .EQ. O_M ) THEN
C               IF ( GLOB_OMM.LT.0.0 ) THEN
C                  D(O_M) = -REAL(TMP)
C               ELSE IF ( GLOB_OMM.GT.5.0 ) THEN
C                  D(O_M) = REAL(TMP)
C               ELSE
C                  D(O_M) = 0.0
C               ENDIF
C            ELSE IF ( I .EQ. O_X ) THEN
C               IF ( GLOB_OMM.LT.0.0 ) THEN
C                  D(O_X) = REAL(TMP)
C               ELSE IF ( GLOB_OMX.LE.0.0 ) THEN
C                  D(O_X) = REAL(TMP)
C               ELSE
C                  D(O_X) = -REAL(TMP)
C               ENDIF
C            ELSE
C               D(I) = 0.0
C            ENDIF
C         ENDDO



      END IF

      LOGMLFUNC = Q

      RETURN
      END
C
C     END OF LOGMLFUNC()
C




C     ****f* snalys/CALCM *
C
C   NAME
C     CALCM() -- Calculate the magnitude for a supernova.
C
C   DESCRIPTION
C     Calculate the magnitude for at supernova at a given
C     luminosity distance.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-03-29
C
C   INPUTS
C     DL      - DL is the the luminosity distance (of type REAL).
C     SCRIPTM - is the value for SCRIPTM.
C
C   USES
C     cosmology.f
C
C   RESULT
C     The magnitude is returned.
C
C   HISTORY
C     2001-10-29 Changed so that SCRIPTM is sent to this function instead
C                of beeing in a common block.
C
C     ***
      REAL FUNCTION CALCM( DL, SCRIPTM )
      IMPLICIT NONE
      REAL DL, SCRIPTM
C     
C         Calculate the observed magnitude of the supernova.
C
      CALCM = SCRIPTM + 5.0*ALOG10(DL*3.E5)
      RETURN
      END
C
C  End of CALCM()
C

      



C     ****i* snalys/UPDATE_GLOBALS *
C
C   NAME
C     UPDATE_GLOBALS() -- Update the global cosmological variables.
C
C   DESCRIPTION
C     For the calculation of luminosity distance and the magnitude
C     global variables are used. This is necessary because the
C     functions have to look in a particular way for the general
C     integration routines to be applicable.
C
C     The function returns TRUE if all the parameters and the
C     combination of them are in a physical sensible physical range.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-03-29
C
C   INPUTS
C     P - is a REAL array of the cosmological parameters for
C         which the luminosity distance and the magnitude
C         will be calculated.
C
C   RESULT
C     A boolean is returned, that indicates if the combination of
C     the cosmological parameters are in a sensible physical range.
C
C   USES
C     cosmo.inc, supernova.inc, operation.inc
C
C   HISTORY
C     2001-09-07 Added the update of GLOB_ALPHA.
C     2001-09-27 Added the update of GLOB_OMK
C     2002-01-08 Added the update of GLOB_H0
C
C     ***
      LOGICAL FUNCTION UPDATE_GLOBALS( P )
      IMPLICIT NONE
      REAL P(*)
      INCLUDE 'cosmo.inc'
      INCLUDE 'supernova.inc'
      INCLUDE 'operation.inc'

      INTEGER I, SLASK
      REAL BIGBANG
      EXTERNAL BIGBANG

C
C        Update the necessary global variables.
C

      DO I = 1 , NR_VAR
C
C           Which variable to update.
C   
         SLASK = 9999
         IF ( I .EQ. O_M ) THEN
            SLASK = O_M
         ELSE IF ( I .EQ. O_X ) THEN
            SLASK = O_X
         ELSE IF ( I .EQ. A_X ) THEN
            SLASK = A_X
         ELSE IF ( I .EQ. A_2 ) THEN
            SLASK = A_2
         ELSE IF ( I .EQ. M_SC ) THEN
            SLASK = M_SC
         ELSE IF ( I .EQ. ALP ) THEN
            SLASK = ALP
         ELSE IF ( I .EQ. SALP ) THEN
            SLASK = SALP
         ELSE IF ( I .EQ. H_0 ) THEN
            SLASK = H_0
         ELSE IF ( I .EQ. O_RC) THEN
            SLASK = O_RC
         ELSE IF ( I .EQ. X_DA) THEN
            SLASK = X_DA
         ELSE IF ( I .EQ. B_1) THEN
            SLASK = B_1
         ELSE IF ( I .EQ. B_2) THEN
            SLASK = B_2
         ELSE IF ( I .EQ. B_inv) THEN
            SLASK = B_inv
         ELSE IF ( I .EQ. B_4) THEN
            SLASK = B_4
         ELSE IF ( I .EQ. L_1) THEN
            SLASK = L_1
         ELSE IF ( I .EQ. L_2) THEN
            SLASK = L_2
         ELSE IF ( I .EQ. W_0) THEN
            SLASK = W_0
         ELSE IF ( I .EQ. W_1) THEN
            SLASK = W_1
         ELSE IF ( I .EQ. Q_HP) THEN
            SLASK = Q_HP
         ELSE IF ( I .EQ. Z_X) THEN
            SLASK = Z_X
         ELSE IF ( I.EQ.FV_0 ) THEN
            SLASK =  FV_0 
         ELSE IF (I.EQ.H_0bar ) THEN
            SLASK = H_0bar
         ENDIF

C
C           If i corresponds to one of the relevant global variables
C           update it!
C
         IF ( SLASK.NE.9999 ) THEN
            IF ( SLASK.EQ.O_M ) THEN
               GLOB_OMM = P(SLASK)
            ELSE IF ( SLASK.EQ.O_X ) THEN
               GLOB_OMX = P(SLASK)
            ELSE IF (I.EQ.A_X) THEN
               GLOB_AX = P(SLASK)
            ELSE IF (I.EQ.A_2) THEN
               GLOB_A2 = P(SLASK)
            ELSE IF ( SLASK.EQ.M_SC ) THEN
               GLOB_MSCRIPT = P(SLASK)
            ELSE IF ( SLASK.EQ.ALP ) THEN
               GLOB_ALPHA = P(SLASK)
            ELSE IF ( SLASK.EQ.SALP ) THEN
               GLOB_SALP = P(SLASK)
            ELSE IF ( SLASK.EQ.H_0 ) THEN
               GLOB_H0 = P(SLASK)
            ELSE IF ( SLASK.EQ.O_RC ) THEN
               GLOB_ORC = P(SLASK)
            ELSE IF ( SLASK.EQ.X_DA ) THEN
               GLOB_XDA = P(SLASK)
            ELSE IF ( SLASK.EQ.B_1 ) THEN
               GLOB_B1 = P(SLASK)
            ELSE IF ( SLASK.EQ.B_2 ) THEN
               GLOB_B2 = P(SLASK)
            ELSE IF ( SLASK.EQ.B_inv ) THEN
               GLOB_Binv = P(SLASK)
            ELSE IF ( SLASK.EQ.B_4 ) THEN
               GLOB_B4 = P(SLASK)
            ELSE IF ( SLASK.EQ.L_1 ) THEN
               GLOB_L1 = P(SLASK)
            ELSE IF ( SLASK.EQ.L_2 ) THEN
               GLOB_L2 = P(SLASK)
            ELSE IF ( SLASK.EQ.W_0 ) THEN
               GLOB_W0 = P(SLASK)
            ELSE IF ( SLASK.EQ.W_1 ) THEN
               GLOB_W1 = P(SLASK)
            ELSE IF ( SLASK.EQ.Q_HP ) THEN
               GLOB_QHP = P(SLASK)
            ELSE IF ( SLASK.EQ.Z_X ) THEN
               GLOB_ZX = P(SLASK)
            ELSE IF ( SLASK.EQ.FV_0 ) THEN
               GLOB_FV0 = P(SLASK)
            ELSE IF ( SLASK.EQ.H_0bar ) THEN
               GLOB_H0bar = P(SLASK)
            ENDIF
         ENDIF
      ENDDO

C
C        Set GLOB_OMK depending on the constraint for the fit. Huh...
C        Check if this makes sense.
C

      IF ( OPER(FLAT) ) THEN
         IF ( OPER(XDIM) ) THEN
            GLOB_ORC = (1.D0 - GLOB_OMM)**2/4.D0
         ELSE
            GLOB_OMX = 1.0 - GLOB_OMM
         ENDIF
         GLOB_OMK = 0.0
      ELSE
         IF ( OPER(XDIM) ) THEN
            GLOB_OMK = 1.D0 -(DSQRT(GLOB_ORC) + 
     $           DSQRT(GLOB_ORC+GLOB_OMM))**2
         ELSE
            GLOB_OMK = 1.0 - GLOB_OMM - GLOB_OMX
         ENDIF
      ENDIF

C
C        Check if the cosmological terms are in the physical range.
C
      IF ( REAL(GLOB_OMX) .LE. BIGBANG(REAL(GLOB_OMM)) .AND.
     $     GLOB_OMM .GE. 0.0 .AND.
     $     ( GLOB_OMM*(1.0 + 5.0)**3 + 
     $       GLOB_OMX*(1.0 + 5.0)**(3.0*(1.0 + GLOB_AX + GLOB_A2*5.0)) +
     $       GLOB_OMK*(1.0 + 5.0)**2 ) .GE. 0.0 .AND.
     $     GLOB_OMX .GE. -1.5 ) THEN
         UPDATE_GLOBALS = .TRUE.
      ELSE
         UPDATE_GLOBALS = .FALSE.
      ENDIF

      RETURN
      END
C
C     End of UPDATE_GLOBALS()
C


C     ****f* snalys/CALCDL *
C
C   NAME
C     CALCDL() -- The luminosity distance.
C
C   DESCRIPTION
C     Calculates the luminosity distance between two red shifts.
C     UPDATE_GLOBALS must be run before using this function to
C     define the cosmological parameters.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-12-18
C
C   USAGE
C     CALCDL( Z1 , Z2 )
C
C   INPUTS
C     Z1 - the first red shift
C     Z2 - the second red shift
C
C   RESULT
C     The luminosity distance is returned.
C
C   WARNINGS
C     - If the luminosity distance is negative 5. is returned as
C       the luminosity distance.
C
C   TODO
C     - Update globals should take care of all cases when
C       the luminosity distance might get negative. That
C       check should be removed from this function.
C
C   USES
C     cosmo.inc
C
C     ***
      REAL FUNCTION CALCDL( Z1, Z2 )
      IMPLICIT NONE
      REAL Z1, Z2

      INCLUDE 'cosmo.inc'
      INCLUDE 'error.inc'
      INCLUDE 'operation.inc'

      DOUBLE PRECISION F_INTEGRAND
      REAL TMP
      LOGICAL INTEGRATE
      EXTERNAL F_INTEGRAND, INTEGRATE
      REAL DLTIMESCAPE
      EXTERNAL DLTIMESCAPE

C
C        Make sure that Z2 is greater than Z1
C
      IF ( Z1 .GT. Z2 ) THEN
         TMP = Z2
         Z2 = Z1
         Z1 = TMP
      END IF

C
C        Calculate the integral

      IF (OPER(DO_WILT)) THEN 
        CALCDL = DLTIMESCAPE(Z2)
        print *,z2,calcdl
        RETURN
      ENDIF

C
      IF ( INTEGRATE( F_INTEGRAND, TMP, Z1, Z2 ) ) THEN
         CALCDL = TMP         

C
C           Calculate d_L, depending on the geometry of the universe.
C
         IF ( REAL(GLOB_OMK) .GT. 0. ) THEN ! open universe
            CALCDL = SINH( SQRT(ABS(GLOB_OMK))*CALCDL )*
     $           (1. + Z2)/SQRT(ABS(GLOB_OMK))
         ELSE IF ( REAL(GLOB_OMK) .LT. 0. ) THEN ! closed universe
            CALCDL = SIN( SQRT(ABS(GLOB_OMK))*CALCDL )*
     $           (1. + Z2)/SQRT(ABS(GLOB_OMK))
         ELSE                   ! flat universe
            CALCDL = CALCDL*(1. + Z2)
         ENDIF
      ELSE
         CALCDL = -1.
      END IF
         

C
C        Check that dL is greater than zero.
C
      IF ( CALCDL .LE. 0.0 ) THEN
C         WRITE (STDERR,'(A,F10.7,A,F10.7,A,F5.2)') 
C     $        'warning: dL is less than zero for Om=',
C     $        GLOB_OMM,', Ox=', GLOB_OMX,' and zs=', Z2
         CALCDL = 5.0
      ENDIF


      RETURN
      END
C
C End of CALCDL()
C


C     ****f* snalys/CALCA *
C
C   NAME
C     CALCA() -- The quantity measured in Baryon Oscillation experiments.
C
C   DESCRIPTION
C     Calculates the quanty A in Eisenstein et al (astro-ph/0501171, eqns 4 & 5)
C     UPDATE_GLOBALS must be run before using this function to
C     define the cosmological parameters.
C
C   AUTHOR
C    Ariel Goobar (ariel@physto.se)
C
C   CREATION DATE
C     2005-08-24
C
C   USAGE
C     CALCA()
C
C   RESULT
C    "A" as in Esenstein et al  
C
C   WARNINGS
C     - 
C       
C
C   TODO
C     - 
C
C   USES
C     cosmo.inc
C
C     ***
      REAL FUNCTION CALCA()
      IMPLICIT NONE
      REAL Z1,EZ1
      PARAMETER (Z1 = 0.35)  
      INCLUDE 'cosmo.inc'
      INCLUDE 'error.inc'

      DOUBLE PRECISION F_INTEGRAND
      REAL TMP
      LOGICAL INTEGRATE
      EXTERNAL F_INTEGRAND, INTEGRATE

      EZ1 =  F_INTEGRAND(DBLE(Z1))  
      EZ1 =  EZ1**(1./3.)*SQRT(GLOB_OMM)
      IF ( INTEGRATE( F_INTEGRAND, TMP, 0.,Z1) ) THEN
  
         IF ( REAL(GLOB_OMK) .GT. 0. ) THEN ! open universe
             TMP = SINH( SQRT(ABS(GLOB_OMK))*TMP )
             TMP = TMP/SQRT(ABS(GLOB_OMK))
         ELSE IF ( REAL(GLOB_OMK) .LT. 0. ) THEN ! closed universe
             TMP = SIN( SQRT(ABS(GLOB_OMK))*TMP)
             TMP = TMP/SQRT(ABS(GLOB_OMK))
         ENDIF

         CALCA = EZ1*(TMP/Z1)**(2./3.)
      ELSE
         CALCA = -1.
      END IF

      RETURN
      END


C     ****f* snalys/CALCCMB *
C
C   NAME
C     CALCCMB() -- The R shift CMB parameter
C
C   DESCRIPTION
C     Calculates the quanty R in Bond, Efstathiou & Tegmark, 1997, MNRAS, 291.L33
C     UPDATE_GLOBALS must be run before using this function to
C     define the cosmological parameters.
C
C   AUTHOR
C    Ariel Goobar (ariel@physto.se)
C
C   CREATION DATE
C     2006-10-18
C
C   USAGE
C     CALCCMB()
C
C   RESULT
C    "R" as in Bond et al  
C
C   WARNINGS
C     - 
C       
C
C   TODO
C     - 
C
C   USES
C     cosmo.inc
C
C     ***
      REAL FUNCTION CALCCMB()
      IMPLICIT NONE
      REAL Z1,Z2,CALCDL
      PARAMETER (Z1 = 0., Z2=1089.0)  
      INCLUDE 'cosmo.inc'
      INCLUDE 'error.inc'
      CALCCMB = CALCDL(Z1,Z2)/(1.+Z2)*SQRT(GLOB_OMM)
      RETURN
      END

C     ****f* snalys/FOREACHSN *
C
C   NAME
C     FOREACHSN() -- Cosmology dependent calculations for each SN.
C
C   DESCRIPTION
C     Computes cosmology dependent variables that have to be done for
C     each supernova, like for exampel the luminosity distance, the
C     magnitudes, the derivatives of the negative log-likelihood function.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-11-02
C
C   USAGE
C     CALL FOREACHSN(Z,IDX,M,DL,NR,MM,SS,SNVS,NLNML,GRAD,NRVAR,STATUS)
C
C   INPUTS
C     Z     - is a REAL vector containing the red shifts.
C     IDX   - is an INTEGER vector containging indices of the sorted SNZS
C             vector in ascending order.
C     M     - is a REAL vector of the same size as Z.
C     DL    - is a REAL vector of the same size as Z.
C     NR    - is an INTEGER giving the length of the vectors. If NR is
C             negative only the luminosity distances will be calculated,
C             and not the magnitudes, the chi2 sum or the derivatives.
C     MM    - is a REAL vector containing the real, measured, magnitudes.
C     SS    - is a REAL vector containing the real, measured stretch.
C     SNVS  - is a REAL vector containing the variances of MM
C     NLNML - is of type REAL
C     GRAD  - is a REAL vector of at least length NRVAR
C     NRVAR - is the number of variables to fit. If this value is negative
C             the derivatives will not be calculated.
C
C
C   SIDE EFFECTS
C     DL     - contains the luminosity distance for each supernova.
C     M      - contains the calculated magnitude for each supernova if NR
C              is positive.
C     NLNML  - the calculated negative log-likelihood function.
C     GRAD   - the gradient of the negative log-likelihood function for
C              the given cosmology.
C     STATUS - is set 1 if everything goes well otherwise zero.
C
C   USES
C     operation.inc, error.inc, cosmo.inc, supernova.inc
C
C   TODO
C     - Maybe the input format is a bit stupid...
C     - Move this routine to cosmo-routines.f
C     - Chek other integration possibilities + the option to use
C       the approximation formula for the luminosity distance.
C     - REAL precision is probably more that enough for the
C       integration.
C     - Is it necessary to return the M vector?
C     - Send the cosmological parameters to this function and run
C       UPDATE_GLOBALS from within. This way Script_M does not have
C       to be part of UPDATE_GLOBALS anymore.
C
C   HISTORY
C     2002-01-05 The case of a DO_MULTI was implemented
C     2002-01-11 A minor bug was corrected in the calculation
C                of the gradient for O_M and O_X for a flat
C                universe.
C     2004-04-17 Rewrote the code to handle bad integrations better.
C                If the integration fails for 1 or more redshifts,
C                the status 0 is returned.
C     
C     ***
      SUBROUTINE FOREACHSN(Z,IDX,M,DL,NR,MM,SS,SNVS,NLNML,GRAD,NRVAR,
     $     STATUS)
      IMPLICIT NONE
      REAL Z(*), M(*) , DL(*), MM(*), SS(*), SNVS(*), NLNML, GRAD(*)
      INTEGER IDX(*)
      INTEGER NR, NRVAR, STATUS

      INCLUDE 'supernova.inc'
      INCLUDE 'cosmo.inc'
      INCLUDE 'error.inc'
      INCLUDE 'gausslin.inc'
      INCLUDE 'probability.inc'
      INCLUDE 'operation.inc'

      INTEGER K, N, I
C
C        The arrays I1, I2, ... , are the results of the integrals
C        described in the document 'Probability distributions, likelihood
C        funcitons and its derivatives'
C
      REAL H0, SM, SA, OX, OM, OK, SOK, ORC, LOWZLIM, CALCM,
     $     B1, B2, Binv, B4, W0, W1, QHP, ZX, L1, L2
      LOGICAL INTEGRATE
      EXTERNAL INTEGRATE, CALCM
      INTEGER NRINT
      PARAMETER ( NRINT = 5 )
      REAL INT(NRINT), TMP
      DOUBLE PRECISION F_INTEGRAND, DF_DOMEGA_M, DF_DOMEGA_M2,
     $     DF_DOMEGA_X, DF_DOMEGA_X2, DF_DALPHA_X, DF_DALPHA_2
      EXTERNAL F_INTEGRAND, DF_DOMEGA_M, DF_DOMEGA_M2, DF_DOMEGA_X,
     $     DF_DOMEGA_X2, DF_DALPHA_X, DF_DALPHA_2
      REAL FLATEPS
      PARAMETER ( FLATEPS = 0. )

      REAL LOGDIST
      EXTERNAL LOGDIST
      REAL CALCDL,CALCERR,IMRATIO 
     
      REAL DLTIMESCAPE
      EXTERNAL DLTIMESCAPE
        

C
C        Derivatives
C
      REAL TRIGARG, MAGDIFF
      REAL G(3)
      REAL V, B, M0, ML
      INTEGER DIST
      REAL CONVUNITS
      REAL DSOU,DLIN,DLS,RMOD,REXP,VARREXP
      REAL SQRRAT
      PARAMETER (CONVUNITS = 13044.8324)
      INTEGER J

      real agtemp
C
C        Make some calculations that will be the same for each
C        supernova
C
      N = ABS(NR)                            ! NR might be negative
      OX = GLOB_OMX
      OM = GLOB_OMM
      OK = GLOB_OMK
      ORC = GLOB_ORC
      SM = GLOB_MSCRIPT
      SA = GLOB_SALP
      H0 = GLOB_H0
      B1 = GLOB_B1
      B2 = GLOB_B2
      Binv = GLOB_Binv
      B4 = GLOB_B4
      L1 = GLOB_L1
      L2 = GLOB_L2
      W0 = GLOB_W0
      W1 = GLOB_W1
      QHP = GLOB_QHP
      ZX = GLOB_ZX
      SOK = SQRT(ABS(OK))
      NLNML = 0.
      LOWZLIM = 0.
      STATUS = 1

C
C        Reset the gradient and integral vectors.
C
      IF ( NRVAR .GT. 0 ) THEN
         DO I = 1, NRVAR
            GRAD(I) = 0.
         ENDDO
      ENDIF
      DO I = 1, NRINT
         INT(I) = 0.
      ENDDO


C
C       Do the stuff for each supernova.
C
      IF (OPER(DO_MULTI)) THEN  ! Fits based on multiple image separation
         DO I = 1, N
            DSOU   =CALCDL(0.,SN_DATA(ZS,I))/(1.+SN_DATA(ZS,I))**2
            DLIN   =CALCDL(0.,SN_DATA(ZL,I))/(1.+SN_DATA(ZL,I))**2
            DLS    =CALCDL(SN_DATA(ZL,I),SN_DATA(ZS,I))
            DLS    =         DLS/(1.+SN_DATA(ZS,I))**2
            RMOD   =DLIN*DSOU/DLS*3000./H0
            IMRATIO=SN_DATA(RATIO,I)/SN_DATA(LENSRAT,I)
C            if (imratio.lt.1) print *,'!!!!!', I,IMRATIO,
C     $               SN_DATA(RATIO,I),SN_DATA(LENSRAT,I),' <<<<<<<<<'
            IF (IMRATIO.GT.1) THEN 
            REXP   =2.*SN_DATA(TDEL,I)/SN_DATA(IMSEP,I)**2*
     $           CONVUNITS/(1.+SN_DATA(ZL,I))
            IF (LENSTYPE.eq.'SIS') THEN
               REXP = REXP*(IMRATIO+1.)/(IMRATIO-1.)
            ELSE                ! Compact object
               SQRRAT = SQRT(IMRATIO)
               REXP=REXP*
     $         (SQRRAT+1./SQRRAT+2.)/(SQRRAT-1./SQRRAT+ALOG(IMRATIO))
            ENDIF
C            VARREXP = CALCERR(I) + (REXP*MULTI_EMODEL*IMRATIO**2)**2
            VARREXP = CALCERR(I)
            REXP = REXP*(1. + MULTI_EMODEL*IMRATIO**2)
            NLNML = NLNML - LOGDIST(RMOD, REXP,0.,0.,VARREXP,0.,0.,0.,1)
           ENDIF
         ENDDO
         RETURN                 ! Exit the routine in case of DO_MULTI
      ENDIF                     ! Huhhh... ugly!
      
      DO I = 1 , N
         K = IDX(I)

         IF ( I .GT. 1 ) LOWZLIM = Z(IDX(I-1))
C
C           If this was a nice language I could have used && and one IF
C           clause would have been enough, but this is not a nice language.
C
C
         IF ( STATUS .EQ. 1) THEN
            IF (OPER(DO_WILT)) THEN
             DL(K) = DLTIMESCAPE(Z(K))
            ELSE               
               

            IF ( INTEGRATE( F_INTEGRAND, TMP, LOWZLIM, Z(K) ) ) THEN
               INT(1) = INT(1) + TMP

C
C                 Calculate geometry dependent factors.
C
               TRIGARG = SOK*INT(1)
               IF ( OK .GT. FLATEPS ) THEN ! open universe
                  DL(K) = SINH(TRIGARG)*(1.0 + Z(K))/SOK
                  IF ( NRVAR .GT. 0 ) THEN ! factors for calc derivatives
                     G(3) =  COSH(TRIGARG)
                     G(1) = - SINH(TRIGARG)/SOK**3
                     G(2) = G(3)/ABS(OK)
                  ENDIF
               ELSE IF ( OK .LT. -FLATEPS ) THEN ! closed universe
                  DL(K) = SIN(TRIGARG)*(1.0 + Z(K))/SOK
                  IF ( NRVAR .GT. 0 ) THEN ! factors for calc derivatives
                     G(3) = COS(TRIGARG)
                     G(1) = SIN(TRIGARG)/SOK**3
                     G(2) = - G(3)/ABS(OK)
                  ENDIF
               ELSE             ! flat universe
                  DL(K) = INT(1)*(1.0 + Z(K))
                  IF ( NRVAR .GT. 0 ) THEN ! factors for calc derivatives
                     G(1) = 0.
                     G(2) = 0.
                     G(3) = 1.
                  ENDIF
               ENDIF
            ELSE
               DL(K) = -1.
            END IF
           ENDIF
           
C
C             Something went wrong, probably a bad combinations of
C             parameters that gave a negative dL. We set dL to a
C             a positive value and set status to 0.
C
            IF ( DL(K) .LE. 0.0 ) THEN
               STATUS = 0
               DL(K) = 50.0
            ENDIF

C
C              Calculate magnitudes
C
            IF ( NR .GT. 0 ) THEN
C
C                 Deterimine the type of distribution function that
C                 should be used.
C
               IF ( WHICH_PDF .EQ. P_GAUSS ) THEN
                  DIST = 1
                  V = SNVS(K)
 
               ELSE IF ( WHICH_PDF .EQ. P_GAUSSLIN ) THEN
                  DIST = 2
C     
C                    Parameters imported from gausslin.inc
C
                  V = GAUSSLINPARS(GLPS,K)**2
                  B = GAUSSLINPARS(GLPB,K)
                  M0 = GAUSSLINPARS(GLPM0,K)
                  ML = GAUSSLINPARS(GLMLIM,K)
               END IF
               M(K) = CALCM(DL(K), SM) ! Calculated mag, m(z) 
               NLNML = NLNML - LOGDIST( M(K), MM(K), SA, ! Contribution to nml.
     $              SS(K), V, B, M0, ML, DIST )
               
C            IF ( GLOB_OMM.LT.1E-7 .AND. GLOB_OMX.LT.1E-7 ) THEN
C               WRITE (*,*) 'ML', K, DL(K), SS(K),
C     $              M(K), MM(K), NLNML
C            ENDIF

C
C              Calculate derivatives
C
C     ( IF THIS CODE IS GOING TO BE USED THE INTEGRATION PART MUST
c       BE RE-WRITTEN SINCE THE FORMAT OF 'INTEGRATE' HAS CHANGED )
C
C            IF ( NRVAR .GT. 0 ) THEN
C               MAGDIFF = (M(K) - MM(K))/SNVS(K)
C               IF ( M_SC .GT. 0 ) GRAD(M_SC) = GRAD(M_SC) + MAGDIFF
C               IF ( O_M .GT. 0 ) THEN
C                  IF ( ABS(OK) .LE. FLATEPS ) THEN
C                     INT(2) = INT(2) +
C     $                    INTEGRATE(DF_DOMEGA_M2,LOWZLIM,Z(K))
C                  ELSE
C                     INT(2) = INT(2) +
C     $                    INTEGRATE(DF_DOMEGA_M,LOWZLIM,Z(K))
C                  ENDIF
C                  GRAD(O_M) = GRAD(O_M) + MAGDIFF*
C     $                 (-.5)*(1 + Z(K))*(G(1) + INT(1)*G(2) +
C     $                    INT(2)*G(3))/DL(K)
C               ENDIF
C               IF ( O_X .GT. 0 ) THEN
C                  IF ( ABS(OK) .LE. FLATEPS ) THEN
C                     INT(3) = INT(3) +
C     $                    INTEGRATE(DF_DOMEGA_X2,LOWZLIM,Z(K))
C                  ELSE
C                     INT(3) = INT(3) +
C     $                    INTEGRATE(DF_DOMEGA_X,LOWZLIM,Z(K))
C                  ENDIF
C                  GRAD(O_X) = GRAD(O_X) + MAGDIFF*
C     $                 (-.5)*(1 + Z(K))*(G(1) + INT(1)*G(2) + 
C     $                 INT(3)*G(3))/DL(K)
C               ENDIF
C               IF ( A_X .GT. 0 ) THEN
C                  INT(4) = INT(4) + INTEGRATE(DF_DALPHA_X,LOWZLIM,Z(K))
C                  GRAD(A_X) = GRAD(A_2) + MAGDIFF*
C     $                 (-1.5)*OX*(1. + Z(K))*INT(4)*G(3)/DL(K)
C               ENDIF
C               IF ( A_2 .GT. 0 ) THEN
C                  INT(5) = INT(5) + INTEGRATE(DF_DALPHA_2,LOWZLIM,Z(K))
C                  GRAD(A_2) = GRAD(A_2) + MAGDIFF*
C     $                 (-1.5)*OX*(1. + Z(K))*INT(5)*G(3)/DL(K)
C               ENDIF
C             ENDIF
            ENDIF
         ENDIF

      ENDDO

C
C        Multiply the gradient with the 5/log(10) constant.
C
      IF ( NRVAR .GT. 0 ) THEN
         DO I = 1, NR_VAR
            IF ( I .NE. M_SC ) GRAD(I) = 5.*GRAD(I)/ALOG(10.)
         ENDDO
      ENDIF
  

      END
C
C     End of FOREACHSN()
C




C     ****f* snalys/NEGMLSM *
C
C   NAME
C     NEGMLSM -- Calculates the negative log-likelihood function
C
C   DESCRIPTION
C     Calculates the negative log-likelihood function given
C     an input value for SCRIPTM. Three variables must be sent to this
C     function in a COMMON block, REAL arrays SNDL, SNVS and SNM, and
C     their length SNLENGTH. The reason why these can not be sent as
C     arguments to the function is that most minimization routines will
C     only accept a function with the only argument beeing to variable
C     with respect to which the function should be minimized. Note
C     however that only SNDL and SNVS have to be defined.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-11-14
C
C   USAGE
C     NEGMLSM( SCRIPTM )
C
C   INPUTS
C     SCRIPTM - The REAL value of script M for which the negative
C               log-likelihood function will be calculated.
C
C   USES
C     supernova.inc, probability.inc
C
C   RESULT
C     The negative log-likelihood function is returned.
C
C   SIDE EFFECTS
C     The COMMON REAL array SNM will on return contain the calculated
C     magnitude for SCRIPTM and the given SNDL vector.
C
C   TODO
C     - make it independent on 'supernova.inc'
C     - make it possible to fit ALPHA, and then remove 'cosmo.inc'
C       the way it is done now is really ugly.
C
C     ***
      REAL FUNCTION NEGMLSM(SCRIPTM)
      IMPLICIT NONE
      INCLUDE 'supernova.inc'
      INCLUDE 'precalcs.inc'
      INCLUDE 'probability.inc'
      INCLUDE 'cosmo.inc'
      REAL SCRIPTM, M
      REAL LOGDIST, CALCM
      EXTERNAL LOGDIST, CALCM
      INTEGER I, DIST

      REAL DL, V, B, M0, ML, MM, SS
      
      NEGMLSM = 0

C
C        Check the type of the probability distribution function
C
      DO I = 1, SN_NR
         MM = SNMS(I)
         SS = SNMS(I)
         DL = SNDL(I)

         IF ( WHICH_PDF .EQ. P_GAUSS ) THEN
            V = SNVS(I)
            DIST = 1
         ELSE IF ( WHICH_PDF .EQ. P_GAUSSLIN ) THEN
            V = GAUSSLINPARS(GLPS,I)**2
            B = GAUSSLINPARS(GLPB,I)
            M0 = GAUSSLINPARS(GLPM0,I)
            ML = GAUSSLINPARS(GLMLIM,I)
            DIST = 2
         ELSE
            DIST = 1
         END IF         

         M = CALCM( DL , SCRIPTM )

         NEGMLSM = NEGMLSM - LOGDIST( M, MM, REAL(GLOB_SALP), SS, V, B,
     $        M0, ML, DIST )

      ENDDO
      RETURN
      END
C
C  End of NEGMLSM
C


C     ****f* snalys/NEGMLSALP *
C
C   NAME
C     NEGMLSALP -- Calculates the negative log-likelihood function
C
C   DESCRIPTION
C     Calculates the negative log-likelihood function given
C     an input value for the strech alpha. Three variables must be sent to this
C     function in a COMMON block, REAL arrays SNDL, SNVS and SNM, and
C     their length SNLENGTH. The reason why these can not be sent as
C     arguments to the function is that most minimization routines will
C     only accept a function with the only argument beeing to variable
C     with respect to which the function should be minimized. Note
C     however that only SNDL and SNVS have to be defined.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2003-04-16
C
C   USAGE
C     NEGMLSALP( SALPHA )
C
C   INPUTS
C     SALPHA - The REAL value of the stretch alpha for which the negative
C              log-likelihood function will be calculated.
C
C   USES
C     supernova.inc, probability.inc
C
C   RESULT
C     The neg-log-likelihood function is returned.
C
C   SIDE EFFECTS
C     The COMMON REAL array SNM will on return contain the calculated
C     magnitude for SCRIPTM and the given SNDL vector.
C
C   TODO
C     - make it independent on 'supernova.inc'
C
C     ***
      REAL FUNCTION NEGMLSALP(SALPHA)
      IMPLICIT NONE
      INCLUDE 'supernova.inc'
      INCLUDE 'precalcs.inc'
      INCLUDE 'probability.inc'
      INCLUDE 'cosmo.inc'
      REAL SALPHA, M
      REAL LOGDIST, CALCM
      EXTERNAL LOGDIST, CALCM
      INTEGER I, DIST

      REAL DL, V, B, M0, ML, MM, SS
      
      NEGMLSALP = 0

C
C        Check the type of the probability distribution function
C
      DO I = 1, SN_NR
         MM = SNMS(I)
         SS = SNMS(I)
         DL = SNDL(I)

         IF ( WHICH_PDF .EQ. P_GAUSS ) THEN
            V = SNVS(I)
            DIST = 1
         ELSE IF ( WHICH_PDF .EQ. P_GAUSSLIN ) THEN
            V = GAUSSLINPARS(GLPS,I)**2
            B = GAUSSLINPARS(GLPB,I)
            M0 = GAUSSLINPARS(GLPM0,I)
            ML = GAUSSLINPARS(GLMLIM,I)
            DIST = 2
         ELSE
            DIST = 1
         END IF
         
         M = CALCM( DL , REAL(GLOB_MSCRIPT) )
         NEGMLSALP = NEGMLSALP - LOGDIST( M, MM, SALPHA, SS, V, B,
     $        M0, ML, DIST )
      ENDDO
      RETURN
      END
C
C  End of NEGMLSALP
C


C     ****f* snalys/NEGMLSMSA *
C
C   NAME
C     NEGMLSMSA -- Calculates the negative log-likelihood function
C
C   DESCRIPTION
C     Calculates the negative log-likelihood function and the gradient
C     given an input value for the scriptM and strech alpha. Three
C     variables must be sent to this function in a COMMON block, REAL
C     arrays SNDL, SNVS and SNM, and their length SNLENGTH. The reason
C     why these can not be sent as arguments to the function is that
C     most minimization routines will only accept a function with the
C     only argument beeing to variable with respect to which the
C     function should be minimized. Note however that only SNDL and SNVS
C     have to be defined.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2003-04-16
C
C   USAGE
C     NEGMLSM( P , G )
C
C   INPUTS
C     P - a 2x1 vector where (1,1) is the value of scriptM and (2,1) is the
C         value of stretch alpha.
C     G - same type and dimension as P, contains the gradient in P on return.
C
C   USES
C     supernova.inc, probability.inc
C
C   RESULT
C     The negative log likelihood function is returned
C
C   SIDE EFFECTS
C     The COMMON REAL array SNM will on return contain the calculated
C     magnitude for SCRIPTM and the given SNDL vector.
C
C   TODO
C     - make it independent on 'supernova.inc'
C
C   HISTORY
C     2004-02-27 - Added a check if the variance is zero.
C
C     ***
      REAL FUNCTION NEGMLSMSA( P , G )
      IMPLICIT NONE
      INCLUDE 'supernova.inc'
      INCLUDE 'precalcs.inc'
      INCLUDE 'probability.inc'
      INCLUDE 'cosmo.inc'
      INCLUDE 'error.inc'
C     Jakob J 2006-09-14
C     I had to include operation.inc to be able to acess OPER
C     I hope this will not mess up the code terribly
      INCLUDE 'operation.inc'
      REAL P(*), G(*), M, SCRIPTM, SALPHA
      REAL LOGDIST, CALCM
      EXTERNAL LOGDIST, CALCM
      INTEGER I, DIST

      REAL DL, V, B, M0, ML, MM, SS, TMP
      
C
C     Add redshift dependent error
C
      REAL Z, ZERROR, DHDZ, H2, DH2DZ, HZ, DZ, DARK, P1, P2

      NEGMLSMSA = 0.
      SCRIPTM = P(1)
      SALPHA = P(2)
      G(1) = 0.
      G(2) = 0.

C
C        Check the type of the probability distribution function
C
      DO I = 1, SN_NR
         MM = SNMS(I)
         SS = SNSS(I)
         DL = SNDL(I)
C     Added by Jakob J 2006-09-13
         Z=SNZS(I)
         DZ=SNDZ(I)

         IF ( WHICH_PDF .EQ. P_GAUSS ) THEN
C     Add Redshift dependent error (2006-09-13 Jakob J)
            IF (OPER(Z_ERROR)) THEN
C     The propagated error in the magnitude due to redshift uncertainty 
C     depends on the cosmology and thus which cosmology is fitted
               IF (OPER(DO_LINDER)) THEN
                  P1=1.+GLOB_L1+GLOB_L2
                  P2=GLOB_L2*Z/(1.+Z)
                  DARK=((1.+Z)**(3.*P1))*EXP(-3.*P2)
                  H2=GLOB_OMM*(1.+Z)**3.+GLOB_OMK*(1+Z)**2.
     $                 +GLOB_OMX*DARK
                  
                  IF (H2.GT.0.) THEN
                     HZ=SQRT(H2)
                  ELSE
                     HZ=1.e16
                     WRITE(*,*) 'mlfunc: H2<=0'
                  ENDIF

                  DARK=3.*P1*((1.+Z)**(3.*P1-1.))*EXP(-3.*P2)+
     $                 ((1.+Z)**(3.*P1))*
     $                 (-3.*GLOB_L2/(1.+Z)+
     $                 3.*GLOB_L2*Z/(1.+Z)**2.)*EXP(-3.*P2)
                  DH2DZ=3.*GLOB_OMM*(1.+Z)**2+2.*GLOB_OMK*(1.+Z)+
     $                 GLOB_OMX*DARK
                  DHDZ=DH2DZ/2./HZ
                  
c                  DZ=0.00005

c                  write(*,*) dz
                  ZERROR = (5./ALOG(10.))*(1./(1.+Z)+
     $                 (1.+Z)/DL/HZ*DHDZ)*DZ

                  V = SNVS(I)+ZERROR**2
c                  print *,z,sqrt(v)
               ELSE
                  WRITE(*,*) 'ZERROR will work only in '
     $                 // 'combination with DO_LINDER'
                  V = SNVS(I)
               ENDIF

 
            ELSE
               V = SNVS(I)
            ENDIF
               
            DIST = 1
         ELSE IF ( WHICH_PDF .EQ. P_GAUSSLIN ) THEN
            V = GAUSSLINPARS(GLPS,I)**2
            B = GAUSSLINPARS(GLPB,I)
            M0 = GAUSSLINPARS(GLPM0,I)
            ML = GAUSSLINPARS(GLMLIM,I)
            DIST = 2
         ELSE
            DIST = 1
         END IF
       
         M = SCRIPTM + 5.0*ALOG10(DL*3.E5) - SALPHA*(SS - 1.0)
c         print *, z, DL, M
         IF ( V .GT. 0.0 ) THEN
            NEGMLSMSA = NEGMLSMSA + .5*(M - MM)**2/V
         ELSE
            WRITE (STDERR,'(2A)') 'Warning, variance = 0 for SN, ',
     $           ' - excluded from fit'
         ENDIF
         
C        IF ( GLOB_OMM.LT.1E-7 .AND. GLOB_OMX.LT.-0.9995 ) THEN
C            WRITE (*,*) 'ML', I, SCRIPTM, SALPHA, SS, M
C            WRITE (*,*) 'RESIDS' , .5*(M - MM)**2/V, SALPHA, SCRIPTM,
C     $           SS, MM
C         ENDIF
        
         G(1) = G(1) + (M - MM)/V
         G(2) = G(2) + (M - MM)*(-(SS - 1.0))/V
      ENDDO      
        
C      IF ( GLOB_OMM.LT.1E-7 .AND. GLOB_OMX.LT.-0.9995 ) THEN
C         WRITE (*,*) 'ML', NEGMLSMSA, SCRIPTM, SALPHA
C      ENDIF

      RETURN
      END
C
C  End of NEGMLSMSA
C


C     ****f* snalys/LOGDIST *
C
C   NAME
C     LOGDIST -- Calculates the logarithm of a 
C
C   DESCRIPTION
C     Calculates the logarithm for a specified distribution function,
C     input magnitudes and the variance.
C

C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-11-14
C
C   USAGE
C     LOGDIST( M , MM , ALPHA , S , V , B , M0 , ML , D )
C
C   INPUTS
C     M  - The calculated magnitude.
C     MM - The measured magnitude.
C     ALPHA - The alpha nuisance parameter
C     S  - Stretch
C     V  - The variance in magnitude.
C     B  - Parameter used for D = 2
C     M0 - Parameter used for D = 2
C     ML - Parameter used for D = 2
C     D  - The distribution to be used.
C          1 = Gaussian distribution
C          2 = Gaussian with 10^(s(M - MM)) term for M - MM < 0
C          If D is less than zero, the gradient will also be calculated.
C
C   USES
C     MAXSN, SN_NR
C
C   RESULT
C     The logarithm of the distribution function.
C
C   NOTE
C     If D = 2 is choosen make sure that the COMMON parameters in
C     'gausslin.inc' are defined.
C
C   SEE ALSO
C     Snalys Note I
C
C   HISTORY
C     2003-03-31 Added stretch to the measured 
C
C     ***
      REAL FUNCTION LOGDIST( M, MM, ALPHA, S, V, B, M0, ML, D )
      IMPLICIT NONE
      REAL M, MM, ALPHA, S, V, B, M0, ML
      INTEGER D
      INTEGER I

      REAL LOGGAUSSLIN
      EXTERNAL LOGGAUSSLIN

C      INCLUDE 'cosmo.inc'

C
C        Gaussian distribution
C
      IF ( ABS(D) .EQ. 1 ) THEN
         LOGDIST = - .5*(M - ALPHA*(S - 1.0) - MM)**2/V
C        IF ( GLOB_OMM .LT. 1.E-7 .AND. GLOB_OMX .LT. 1.E-7 ) THEN
C           WRITE (*,*) 'RESIDS ',  LOGDIST, ALPHA,
C     $           GLOB_MSCRIPT, S, MM
C        ENDIF
      ELSE IF ( ABS(D) .EQ. 2 ) THEN
         LOGDIST = LOGGAUSSLIN( M , MM , V , B , M0 , ML )
      END IF

      RETURN
      END
C
C End of LOGDIST
C



C     ****f* snalys/NMLGAUSSLIN *
C
C   NAME
C     NMLGAUSSLIN -- Calculates the negative log-likelihood function
C
C   DESCRIPTION
C     Calculates the negative log-likelihood function for a Gaussian
C     distribution with a 10*(m) tail given an input value for MSCRIPT.
C     Three variables must be sent to this function in a COMMON block,
C     REAL arrays SNDL, SNVS and SNM, and their length SNLENGTH. The reason
C     why these can not be sent as arguments to the function is that
C     most minimization routines will only accept a function with the
C     only argument beeing to variable with respect to which the
C     function should be minimized. Note however that only SNDL and SNVS
C     have to be defined.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-11-13
C
C   USAGE
C     NMLGAUSSLIN( SCRIPTM  )
C
C   INPUTS
C     SCRIPTM - The REAL value of SCRIPTM for which the NMLGAUSSLIN will
C               be calculated.
C
C   USES
C     MAXSN, SN_NR
C
C   RESULT
C     The magnitude is returned.
C
C   SIDE EFFECTS
C     The COMMON REAL array SNM will on return contain the calculated
C     magnitude for MSCRIPT and the given SNDL vector.
C
C   NOTE
C     SNDL is not the luminosity distances but 5.LOG10(dL*3.E5) where
C     dL is the luminosity distance.
C
C   TODO
C     - make it independent on 'supernova.inc'
C
C   SEE ALSO
C     - Snalys Note I
C
C     ***
      REAL FUNCTION NMLGAUSSLIN(SCRIPTM)
      IMPLICIT NONE
      REAL SCRIPTM
      INCLUDE 'supernova.inc'
      INCLUDE 'precalcs.inc'
      INTEGER I

      REAL V, B, M0, ML
      REAL LOGGAUSSLIN, CALCM
      EXTERNAL LOGGAUSSLIN, CALCM
      
      NMLGAUSSLIN = 0
      V = GAUSSLINPARS(GLPS,I)**2
      B = GAUSSLINPARS(GLPB,I)
      M0 = GAUSSLINPARS(GLPM0,I)
      ML = GAUSSLINPARS(GLMLIM,I)
      
      DO I = 1, SN_NR
         SNM(I) = CALCM(SNDL(I), SCRIPTM)
         NMLGAUSSLIN = NMLGAUSSLIN -
     $        LOGGAUSSLIN( SNM(I) , SNMS(I) , V , B , M0 , ML )
      ENDDO      
      RETURN
      END
C
C End of NMLGAUSSLIN
C



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                              C
C        THE INTEGRANDS IN THE NUMERICAL INTEGRATION           C
C                                                              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C     ****i* snalys/F_INTEGRAND *
C
C   NAME
C     F_INTEGRAND() -- The integrand in the expression for the
C                      luminosity distance.
C
C   DESCRIPTION
C     The integrand in the expression for the luminosity distance. This
C     integrand is integrated from 0 to z to receive D_L.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-03-29
C
C   INPUTS
C     ZP - is a redshift of type DOUBLE PRECISION.
C
C   RESULT
C     The result of the integrand for this redshift is returned.
C
C   TODO
C     - Since this function is called very often in the integration
C       as few calculations as possible should be made here, that is
C       all calculations that does not involve z should be done before
C       entering the integration procedure.
C
C   USES
C     cosmo.inc, operation.inc, mate.f
C
C   HISTORY
C     2001-09-07 Made the proper modification so that the Eriksson
C                quintessence stuff can be tested too.
C     2004-01-19 Added the integrand for Jakob
C     2004-02-16 Jakobs integrand was changed
C     2004-04-16 The parameter 'CONVERGED' is set to .FALSE. in case
C                of divergation.
C     2005-05-17 Added the integrand for Edvard (Hannestad parametrization)
C     2005-10-18 modified the DGP xdim integrand for general powers
C    
C     ***
      DOUBLE PRECISION FUNCTION F_INTEGRAND( ZP )
      IMPLICIT NONE
      DOUBLE PRECISION ZP

      INCLUDE 'cosmo.inc'
      INCLUDE 'operation.inc'
      INCLUDE 'integrate.inc'

      DOUBLE PRECISION H_FACTOR, DARK, MATE
      EXTERNAL MATE

      DOUBLE PRECISION XDH,XDF,XDC
      DOUBLE PRECISION XDIMOMK
      EXTERNAL XDIMOMK

      REAL HTIMESCAPE
      EXTERNAL HTIMESCAPE
      
      IF ( OPER(MATE_EQST) ) THEN
         DARK = MATE( ZP )
      ELSEIF ( OPER(DO_HP) ) THEN
         IF (ABS(GLOB_QHP) .LT. 1.D-6) THEN
            DARK = (1.D0 + ZP)**(3.D0*(1.D0 + 
     $           2.D0*GLOB_W0*GLOB_W1/(GLOB_W0 + GLOB_W1)))
         ELSEIF ((ABS(GLOB_W0) .LT. 1.D-6) .OR. 
     $           (ABS(GLOB_W1) .LT. 1.D-6))  THEN
            DARK = (1.D0 + ZP)**3.D0
         ELSE
            DARK = (1.D0 + ZP)**(3.D0*(1.D0 + GLOB_W0))*
     $           ( (GLOB_W0*(1.D0 + ZP)**GLOB_QHP + 
     $           GLOB_W1*(1.D0 + GLOB_ZX)**GLOB_QHP)/
     $           (GLOB_W0 + GLOB_W1*(1.D0 + GLOB_ZX)**GLOB_QHP) )**
     $           (3.D0*(GLOB_W1-GLOB_W0)/GLOB_QHP)
         ENDIF
      ELSEIF ( .NOT. OPER(DO_JACKE) ) THEN
         DARK = H_FACTOR( ZP )*(1.D0 + ZP)**(3.D0*(1.D0 + GLOB_AX +
     $        + ZP*GLOB_A2))
      ELSEIF ( .NOT. OPER(DO_LINDER) ) THEN
         DARK = H_FACTOR( ZP )*(1.D0 + ZP)**(3.D0*(1.D0 + GLOB_AX +
     $        + ZP*GLOB_A2))
      ENDIF
      
      IF ( OPER(XDIM).and.( ABS(GLOB_XDA-1.0d0).LE.1.D-5  )) THEN
         F_INTEGRAND = GLOB_OMK*(1.D0 +ZP)**2 + (DSQRT(GLOB_ORC) + 
     $        DSQRT(GLOB_OMM*(1+ZP)**3+GLOB_ORC))**2
         

      ELSEIF (  OPER(XDIM).and.(abs(GLOB_XDA-1.D0).GT.1.0d-5)) THEN
         IF ( OPER(FLAT) ) THEN
            GLOB_OMK = 0.D0
         ELSE
            GLOB_OMK = XDIMOMK(GLOB_OMM,GLOB_ORC,GLOB_XDA)
c     print *,GLOB_OMM,GLOB_ORC,GLOB_XDA
c       print *,GLOB_OMM+GLOB_OMK+2.d0*DSQRT(GLOB_ORC)*
c     $        (1.D0-GLOB_OMK)**GLOB_XDA/2.D0
         ENDIF
         IF (GLOB_OMK .GT. 1) THEN
            F_INTEGRAND = 0.d0
            CONVERGED = .FALSE.
            RETURN
         ENDIF
         XDF = 1.0D0
         XDC = 2.D0*DSQRT(GLOB_ORC)
         XDH = 1.0D0
         
 10      IF(XDF.LE.1.0D-6) GOTO 30
         IF(XDC*(XDH**GLOB_XDA).GT.XDH*XDH-GLOB_OMM*(1.D0+ZP)**3)THEN
            XDH = XDH+XDF
            GOTO 10
         ELSE
            XDF = XDF/2.0D0
            GOTO 20
      ENDIF


C
C This does not make sense to me!!!!  It will do this whether we are
C interested in x-dimensions or not!
C
 20   IF(XDF.LE.1.0D-6) GOTO 30
      IF(XDC*(XDH**GLOB_XDA).LT.XDH*XDH-GLOB_OMM*(1.D0+ZP)**3)THEN  
         XDH=XDH-XDF
         GOTO 20
      ELSE
         XDF=XDF/2.0D0
         GOTO 10
      ENDIF
 30   F_INTEGRAND=XDH*XDH+GLOB_OMK*(1.D0+ZP)**2
      
      ELSEIF ( OPER(DO_JACKE) ) THEN
         F_INTEGRAND = GLOB_OMM*(1.D0 + ZP)**3 + 1.D0 - GLOB_OMM -
     $        GLOB_B1 - GLOB_B2 - GLOB_Binv - GLOB_B4 + 
     $        GLOB_B1*(1.D0 + ZP) + GLOB_B2*(1.D0 + ZP)**2 +
     $        GLOB_Binv/(1.D0+ZP) + GLOB_B4*(1.D0 + ZP)**4
      ELSEIF ( OPER(DO_LINDER) ) THEN
         DARK=DEXP(-3.D0*GLOB_L2*ZP/(1.D0+ZP))
         DARK=DARK*(1.D0+ZP)**(3.D0*(1.D0+GLOB_L1+GLOB_L2))
         F_INTEGRAND = GLOB_OMM*(1.D0 + ZP)**3 + 
     $        GLOB_OMK*(1.D0 + ZP)**2 + GLOB_OMX*DARK
      ELSEIF ( .NOT. OPER(DO_WILT) ) THEN
         F_INTEGRAND = GLOB_OMM*(1.D0 + ZP)**3 + 
     $        GLOB_OMK*(1.D0 + ZP)**2 + GLOB_OMX*DARK
      ENDIF
      
      IF ( OPER(DO_WILT) ) THEN
         F_INTEGRAND = 1.D0/HTIMESCAPE(ZP)
      ELSEIF ( F_INTEGRAND.GE.0.D0 ) THEN
         F_INTEGRAND = 1.D0/DSQRT(F_INTEGRAND)
      ELSE
         F_INTEGRAND = 0.D0
         CONVERGED = .FALSE.

C         WRITE (STDERR,'(A)') 'error: f_integrand is imaginary'
C         WRITE (STDERR,*) GLOB_OMM, GLOB_OMX, GLOB_AX, GLOB_A2, ZP,
C     $        F_INTEGRAND
      ENDIF

      RETURN
      END
C
C        End of F_integrand()
C
C
C
C     XDIMOMK
C     
C     Author. Malcom Fairbairn 2005-10-19
C     Included into Snalys by Ariel Goobar
C
C     Calculates curvature in Xdim when GLOB_XDA <> 1
C
      DOUBLE PRECISION FUNCTION XDIMOMK(omm,omrc,xda)
      IMPLICIT NONE      

      DOUBLE PRECISION OMK,GLOB_OMM,GLOB_OMR,f,XDA,OMRC,OMM      
      DOUBLE PRECISION GLOB_XDA,DIFF,C,STEP

      omk=0.2D0
      f=1.0d0
      glob_omm = omm
      glob_omr = omrc
      glob_xda = xda
      xdimomk = 0.d0
      IF (glob_xda .GE. 0) THEN 
       c =  1.D0
      ELSE
       c = -1.D0
      ENDIF
 10   CONTINUE
      IF ((glob_omm+omk) .GE. 1.0d0) RETURN
      diff = 1.D0 - glob_omm - omk
      diff = diff-2.d0*dsqrt(glob_omr)*(1.d0-omk)**(glob_xda/2.d0)
      IF (abs(diff) .LT. 1.d-3) GOTO 30
      f=f/1.5d0
      IF (diff .GT. 0 ) THEN      
 12     step = c*f
        IF (omk + step + glob_omm .LT. 1) THEN 
         omk=omk+step
        ELSE
         f = f/2.D0
         goto 12 
        ENDIF
        goto 10
      ELSE        
        goto 20
      ENDIF

 20   IF (f .LE. 1.0d-6) goto 30
      diff = 1.D0 - glob_omm - omk
      diff = diff-2.d0*dsqrt(glob_omr)*(1.d0-omk)**(glob_xda/2.d0) 
      IF (abs(diff) .LT. 1.d-3) GOTO 30
      IF (diff .LT. 0) THEN
         step = c*f
         omk=omk + step 
         goto 20
      ELSE
         f=f/1.5d0
         goto 10
      ENDIF
 30   CONTINUE
      xdimomk=omk 
      RETURN
      END
C
C
C     ****i* snalys/DF_DOMEGA_M *
C
C   NAME
C     DF_DOMEGA_M() -- The integrand in the expression for the
C                      derivative of OMEGA_M.
C
C   DESCRIPTION
C     The integrand in the expression for the calculating the
C     derivative of F with respect to OMEGA_M
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-03-29
C
C   INPUTS
C     ZP - is a redshift of type DOUBLE PRECISION.
C
C   RESULT
C     The result of the integrand for this redshift is returned.
C
C   TODO
C     - Since this function is called very often in the integration
C       as few calculations as possible should be made here, that is
C       all calculations that does not involve z should be done before
C       entering the integration procedure.
C
C     ***
      DOUBLE PRECISION FUNCTION DF_DOMEGA_M( ZP )
      IMPLICIT NONE
      DOUBLE PRECISION ZP

      DOUBLE PRECISION F_INTEGRAND
      
      DF_DOMEGA_M = ( (1.0 + ZP)**3 - (1.0 + ZP)**2 )*
     $     (F_INTEGRAND(ZP))**3

      RETURN
      END
C
C        End of DF_DOMEGA_M()
C



C     ****i* snalys/DF_DOMEGA_M2 *
C
C   NAME
C     DF_DOMEGA_M2() -- The integrand in the expression for the
C                       derivative of OMEGA_M for a flat geometry.
C
C   DESCRIPTION
C     The integrand in the expression for the calculating the
C     derivative of F with respect to OMEGA_M for a flat geometry.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-03-29
C
C   INPUTS
C     ZP - is a redshift of type DOUBLE PRECISION.
C
C   RESULT
C     The result of the integrand for this redshift is returned.
C
C   TODO
C     - Since this function is called very often in the integration
C       as few calculations as possible should be made here, that is
C       all calculations that does not involve z should be done before
C       entering the integration procedure.
C
C     ***
      DOUBLE PRECISION FUNCTION DF_DOMEGA_M2( ZP )
      IMPLICIT NONE
      DOUBLE PRECISION ZP

      DOUBLE PRECISION F_INTEGRAND
      
      DF_DOMEGA_M2 = ( (1.0 + ZP)**3 )*
     $     (F_INTEGRAND(ZP))**3

      RETURN
      END
C
C        End of DF_DOMEGA_M()
C



C     ****i* snalys/DF_DOMEGA_X *
C
C   NAME
C     DF_DOMEGA_X() -- The integrand in the expression for the
C                      derivative of OMEGA_X.
C
C   DESCRIPTION
C     The integrand in the expression for the calculating the
C     derivative of F with respect to OMEGA_X.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-03-29
C
C   INPUTS
C     ZP - is a redshift of type DOUBLE PRECISION.
C
C   RESULT
C     The result of the integrand for this redshift is returned.
C
C   TODO
C     - Since this function is called very often in the integration
C       as few calculations as possible should be made here, that is
C       all calculations that does not involve z should be done before
C       entering the integration procedure.
C
C     ***
      DOUBLE PRECISION FUNCTION DF_DOMEGA_X( ZP )
      IMPLICIT NONE
      DOUBLE PRECISION ZP

      INCLUDE 'cosmo.inc'
      DOUBLE PRECISION F_INTEGRAND, H_FACTOR
      
      DF_DOMEGA_X = ( H_FACTOR(ZP)*
     $     (1.0 + ZP)**(3.0*(1.0 + GLOB_AX + ZP*GLOB_A2)) -
     $     (1.0 + ZP)**2 )*( F_INTEGRAND(ZP) )**3

      RETURN
      END
C
C        End of DF_DOMEGA_X()
C



C     ****i* snalys/DF_DOMEGA_X2 *
C
C   NAME
C     DF_DOMEGA_X2() -- The integrand in the expression for the
C                       derivative of OMEGA_X for a flat geometry.
C
C   DESCRIPTION
C     The integrand in the expression for the calculating the
C     derivative of F with respect to OMEGA_X for a flat geometry.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2002-01-11
C
C   INPUTS
C     ZP - is a redshift of type DOUBLE PRECISION.
C
C   RESULT
C     The result of the integrand for this redshift is returned.
C
C   TODO
C     - Since this function is called very often in the integration
C       as few calculations as possible should be made here, that is
C       all calculations that does not involve z should be done before
C       entering the integration procedure.
C
C     ***
      DOUBLE PRECISION FUNCTION DF_DOMEGA_X2( ZP )
      IMPLICIT NONE
      DOUBLE PRECISION ZP

      INCLUDE 'cosmo.inc'
      DOUBLE PRECISION F_INTEGRAND, H_FACTOR
      
      DF_DOMEGA_X2 = ( H_FACTOR(ZP)*
     $     (1.0 + ZP)**(3.0*(1.0 + GLOB_AX + ZP*GLOB_A2)) )*
     $     ( F_INTEGRAND(ZP) )**3

      RETURN
      END
C
C        End of DF_DOMEGA_X()
C




C     ****i* snalys/DF_DALPHA_X *
C
C   NAME
C     DF_DALPHA_X() -- The integrand in the expression for the
C                      derivative of DF_DALPHA_X.
C
C   DESCRIPTION
C     The integrand in the expression for the calculating the
C     derivative of F with respect to ALPHA_X.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-03-29
C
C   INPUTS
C     ZP - is a redshift of type DOUBLE PRECISION.
C
C   RESULT
C     The result of the integrand for this redshift is returned.
C
C   TODO
C     - Since this function is called very often in the integration
C       as few calculations as possible should be made here, that is
C       all calculations that does not involve z should be done before
C       entering the integration procedure.
C
C     ***
      DOUBLE PRECISION FUNCTION DF_DALPHA_X(ZP)
      IMPLICIT NONE
      DOUBLE PRECISION ZP

      INCLUDE 'cosmo.inc'
      DOUBLE PRECISION F_INTEGRAND, H_FACTOR
      
      DF_DALPHA_X = (1.0 + ZP)**(3.0*(1.0 + GLOB_AX + GLOB_A2*ZP))*
     +     DLOG(1.0 + ZP)*H_FACTOR(ZP)*(F_INTEGRAND(ZP))**3

      RETURN
      END
C
C   End of DF_DALPHA_X()
C




C     ****i* snalys/DF_DALPHA_2 *
C
C   NAME
C     DF_DALPHA_2() -- The integrand in the expression for the
C                      derivative of DF_DALPHA_2.
C
C   DESCRIPTION
C     The integrand in the expression for the calculating the
C     derivative of F with respect to ALPHA_2.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-03-29
C
C   INPUTS
C     ZP - is a redshift of type DOUBLE PRECISION.
C
C   RESULT
C     The result of the integrand for this redshift is returned.
C
C   TODO
C     - Since this function is called very often in the integration
C       as few calculations as possible should be made here, that is
C       all calculations that does not involve z should be done before
C       entering the integration procedure.
C
C     ***
      DOUBLE PRECISION FUNCTION DF_DALPHA_2( ZP )
      IMPLICIT NONE
      DOUBLE PRECISION ZP

      INCLUDE 'cosmo.inc'
      DOUBLE PRECISION F_INTEGRAND, H_FACTOR
      
      DF_DALPHA_2 = (1.0 + ZP)**(3.0*(1.0 + GLOB_AX + GLOB_A2*ZP))*
     +     H_FACTOR(ZP)*(ZP-DLOG(1.0+ZP))*
     +     (F_INTEGRAND(ZP))**3

      RETURN
      END
C
C   End of DF_DALPHA_2()
C



C     ****i* snalys/H_FACTOR *
C
C   NAME
C     H_FACTOR() -- Calculates the h-factor.
C
C   DESCRIPTION
C     When the equation of state is assumed to be z-dependent
C     this factor will no longer be one.
C
C   NOTE
C     This implementation of the h-factor is only valid when
C     the time dependence of the equation of state is linear.
C
C   AUTHOR
C     Martin Goliath (goliath@physto.se)
C
C   CREATION DATE
C     2000-09-22
C
C   INPUTS
C     ZP - is a redshift of type DOUBLE PRECISION.
C
C   RESULT
C     The h-factor is returned.
C
C     ***
      DOUBLE PRECISION FUNCTION H_FACTOR( ZP )
      IMPLICIT NONE
      DOUBLE PRECISION ZP

      INCLUDE 'cosmo.inc'

      H_FACTOR = DEXP(3.D0*GLOB_A2*(ZP - (1.0+ZP)*DLOG(1.0+ZP)))

      RETURN
      END
C
C   End of H_FACTOR
C

C     ****i* snalys/CALCERR *
C
C   NAME
C     CALCERR -- Calculates errors for multiple lensing fitting
C
C   AUTHOR
C     Ariel Goobar (ariel@physto.se)
C
C     ***
      REAL FUNCTION CALCERR(isn)
      IMPLICIT NONE
      INCLUDE 'supernova.inc'
      INTEGER isn
      REAL dRdt,dRdtheta,dRdzl,dRdratio
      REAL ratfac,oneovz,sigt,sigthe,sigzl,sigr
      REAL CONVUNITS,imratio
      PARAMETER (CONVUNITS = 13044.8324)
      sigt   = SN_DATA(ETDEL,isn)/365.    !Convert back to years
      sigthe = SN_DATA(EIMSEP,isn)
      sigzl  = SN_DATA(EZL,isn)
      sigr   = SN_DATA(ERATIO,isn)
*
C      IF (lenstype.eq.'SIS') THEN    ! Must write the code for the other lens types
       imratio  = SN_DATA(RATIO,isn)/SN_DATA(LENSRAT,isn)
       ratfac = (IMRATIO+1.)/(IMRATIO-1.)
       oneovz = 1./(1.+SN_DATA(ZL,isn))
       dRdt = 2./SN_DATA(IMSEP,isn)**2*ratfac*oneovz
       dRdtheta = 2.*SN_DATA(TDEL,isn)*ratfac*oneovz
       dRdtheta =             dRdtheta*(-2./SN_DATA(IMSEP,isn)**3)
       dRdzl    = 2.*SN_DATA(TDEL,isn)/SN_DATA(IMSEP,isn)**2
       dRdzl    =             dRdzl*ratfac*(-1)*oneovz**2
       dRdratio = -4.*SN_DATA(TDEL,isn)/SN_DATA(IMSEP,isn)**2
       dRdratio =        dRdratio*oneovz/(imratio-1.)**2 
       calcerr = (dRdt*sigt)**2+(dRdtheta*sigthe)**2
       calcerr = (calcerr+(dRdzl*sigzl)**2+(dRdratio*sigr)**2)
       calcerr = calcerr*CONVUNITS**2
C      ELSE
C         calculate dRdratio of compact objects...
C      ENDIF 
      END
C
C End of CALCERR
C

C     ****i* snalys/DLTIMESCAPE 
C
C   NAME
C     DLTIMESCAPE -- Calculates DL(z) for "Timescape cosmology"
C                   as in suggested arXiv:0909.0749; D.Wiltshire
C  
C   GLOBAL FIT-PARAMETERS used in wiltfunc: GLOB_fv0,GLOB_H0bar
C   Typical values: fv0   = 0.76
C                   H0bar = 48.
C
C   Comments: The original expressions relate redshift vs time z(t)
C             , thus we resort to a numerical
C             routine zriddr (from Numerical Recipes) to find t(z),
C             which is what we need  to find DL(z)
C
C   AUTHOR
C     Ariel Goobar (ariel@physto.se) 2010-02-05




      REAL FUNCTION DLTIMESCAPE (ZP)

      IMPLICIT NONE
      REAL ZP
      REAL z,x1,x2,t,t0,hofz,ft,wiltfunc,ft0
      EXTERNAL wiltfunc

      REAL zriddr
      EXTERNAL zriddr
      REAL tmp

      common/wiltstuff/z      
      REAL b,fv0,H0bar, conv
      PARAMETER  (conv = 978.56)
 
      include 'cosmo.inc' 

      H0bar = GLOB_H0bar/conv   ! in unitis of Gyr^-1
      fv0   = GLOB_fv0
      x1 = 1.                   !
      x2 = 50.                  ! Range in Gyrs for the Universe to be studied
C
C     Here is where everything happens
     
      z = 0.
      t = zriddr(wiltfunc,x1,x2,0.001) ! Finds t0
      t0 = t
      b = 2.*(1. - fv0)*(2. + fv0)/(9.*fv0*H0bar)     
      tmp = (t**(1./3.) + b**(1./3.))**2
      tmp = tmp/(t**(2./3.)-b**(1./3.)*t**(1./3.)+b**(2./3.))
      tmp = alog(tmp)
      ft = 2.*t**(1./3.) + b**(1./3.)/6.*tmp 
      tmp = (2.*t**(1./3.)-b**(1./3.))/sqrt(3.)/b**(1./3.)
      tmp = atan(tmp)*b**(1./3.)/sqrt(3.)
      ft = ft + tmp
      ft0 = ft
c     print *,'t0=',t0,'ft0=',ft0      
      z  = ZP
      x1 = 1.
      x2 = 35.     
      t= zriddr(wiltfunc,x1,x2,0.001) ! Finds t
      tmp = (t**(1./3.) + b**(1./3.))**2
      tmp = tmp/(t**(2./3.)-b**(1./3.)*t**(1./3.)+b**(2./3.))
      tmp = alog(tmp)
      ft = 2.*t**(1./3.) + b**(1./3.)/6.*tmp 
      tmp = (2.*t**(1./3.)-b**(1./3.))/sqrt(3.)/b**(1./3.)
      tmp = atan(tmp)*b**(1./3.)/sqrt(3.)
      ft = ft + tmp
      DLTIMESCAPE = 306.6*t**(2./3.)*(ft0-ft)*(1.+ z)**2 ! in Mpc
      RETURN
      END

C     

C     ****i* snalys/HTIMESCAPE *
C
C   NAME
C     HTIMESCAPE -- Calculates H(z) for "Timescape cosmology"
C                   as in suggested arXiv:0909.0749; D.Wiltshire
C  
C   GLOBAL FIT-PARAMETERS used in wiltfunc: GLOB_fv0,GLOB_H0bar
C   Typical values: fv0   = 0.76
C                   H0bar = 48.
C
C   Comments: The original expressions relate redshift vs time z(t)
C             and distances vs time H(t), thus we resort to a numerical
C             routine zriddr (from Numerical Recipes) to find t(z),
C             which is what we need  to find H(z)
C
C   AUTHOR
C     Ariel Goobar (ariel@physto.se) 2010-02-01


      REAL FUNCTION HTIMESCAPE (ZP)

      IMPLICIT NONE
      DOUBLE PRECISION ZP
      REAL x1,x2,t,wiltfunc
      EXTERNAL wiltfunc

      REAL zriddr
      EXTERNAL zriddr

      include 'cosmo.inc'

      REAL b,fv0,H0bar, conv,hofz
      
      REAL z
      common/wiltstuff/z
      PARAMETER  (conv = 978.56)
 
      H0bar = GLOB_H0bar/conv   ! in unitis of Gyr^-1
      fv0   = GLOB_fv0

      b = 2.*(1. - fv0)*(2. + fv0)/(9.*fv0*H0bar)     
      z  = REAL(ZP)
      x1 = 1.                   !
      x2 = 35.                  ! Range in Gyrs for the Universe to be studied
      t = zriddr(wiltfunc,x1,x2,0.0001) ! Finds t that matches ZP

      hofz = 3.*(2.*t**2. + 3.*b*t + 2.*b**2.)
      hofz = hofz/t/(2.*t+3.*b)**2.*conv/70. ! Units fixed to give mscript -3.5                 


      HTIMESCAPE = hofz
 
      RETURN
      END
 


      REAL FUNCTION zriddr(func,x1,x2,xacc)
      INTEGER MAXIT
      REAL x1,x2,xacc,func,UNUSED
      PARAMETER (MAXIT=60,UNUSED=-1.11E30)
      EXTERNAL func
CU    USES func
      INTEGER j
      REAL fh,fl,fm,fnew,s,xh,xl,xm,xnew
      fl=func(x1)
      fh=func(x2)
      if((fl.gt.0..and.fh.lt.0.).or.(fl.lt.0..and.fh.gt.0.))then
        xl=x1
        xh=x2
        zriddr=UNUSED
        do 11 j=1,MAXIT
          xm=0.5*(xl+xh)
          fm=func(xm)
          s=sqrt(fm**2-fl*fh)
          if(s.eq.0.)return
          xnew=xm+(xm-xl)*(sign(1.,fl-fh)*fm/s)
          if (abs(xnew-zriddr).le.xacc) return
          zriddr=xnew
          fnew=func(zriddr)
          if (fnew.eq.0.) return
          if(sign(fm,fnew).ne.fm) then
            xl=xm
            fl=fm
            xh=zriddr
            fh=fnew
          else if(sign(fl,fnew).ne.fl) then
            xh=zriddr
            fh=fnew
          else if(sign(fh,fnew).ne.fh) then
            xl=zriddr
            fl=fnew
          else
            pause 'never get here in zriddr'
          endif
          if(abs(xh-xl).le.xacc) return
 11             continue
        pause 'zriddr exceed maximum iterations'
      else if (fl.eq.0.) then
        zriddr=x1
      else if (fh.eq.0.) then
        zriddr=x2
      else
        pause 'root must be bracketed in zriddr'
      endif
      return
      END
 


      REAL FUNCTION wiltfunc(t)
      IMPLICIT NONE
      INCLUDE 'cosmo.inc'
      REAL b,fv0,H0bar,t, conv,hofz,ft,redshift
      REAL tmp


      REAL z

      common/wiltstuff/z
      PARAMETER  (conv = 978.56) 

      H0bar = GLOB_H0bar/conv   ! in unitis of Gyr^-1
      fv0   = GLOB_fv0
      
       
      b = 2.*(1. - fv0)*(2. + fv0)/(9.*fv0*H0bar)
      redshift = 2.**(4./3.)*t**(1./3.)*(t + b)
      redshift=redshift/(fv0**(1./3.)*H0bar*t*(2.*t+3.*b)**(4./3.))- 1.


      wiltfunc = redshift - z                
     
      RETURN
      END



C
C   End of file mlfunc.f
C
