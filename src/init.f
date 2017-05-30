C     ****f* snalys/INITIALIZE *
C
C   NAME
C     INITIALIZE -- handling initialization of parameters files etc.
C
C   DESCRIPTION
C     This is the main init function that handles the initialization of
C     parameters, reading input files, etc.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-10-25
C
C   USAGE
C     INITIALIZE( P_COS, FNAME )
C
C   INPUTS
C     P_COS   - A REAL array
C     FNAME   - CHARACTER string
C
C   RESULT
C     A LOGICAL variable is returned confiriming that the intializtion
C     was successful.
C
C   SIDE EFFECTS
C     P_COS   - if minimization according to Davidon's, Powell's or Genetic's
C               method has been choosen this vector will contain a
C               starting cosmology on return.
C     FNAME   - on return this will contain the name of the output
C               file.
C     COMMON variables from operation.inc, supernova.inc, snalys.inc are
C     set and probably a lot more.
C
C   ERRORS
C     Error message will be written if:
C     - the init file could not be opened.
C     - too many grid steps have been choosen.
C
C   WARNINGS
C     Warning message will be written if:
C     - the cluster option is used on the command line.
C
C   TODO
C
C   HISTORY
C     2001-09-07 Added the option of also estimating the quintessence
C                parameter alpha in the model implemented by Martin
C                Eriksson.
C     2003-03-31 Added the option of also fitting the stretch alpha.
C     2003-07-09 Removed the plotting options
C                The WCUBEFILE is no longer hardcoded
C     2003-09-05 Removed all command line options that could be specified
C                par-file, and removed the default par-file
C     2004-01-15 Added the option of fitting "jacke"-cosmology
C     2004-04-03 Changed the default stretch from 0.0 to 1.0
C     2004-05-17 Added the option of fitting Hannestad parametrization (Edvard)
C     2005-02-16 Added HOSTDM, MIXDM
C     2005-08-08 Added GDM
C     2010-02-01 Added Wiltshare's stuff
C
C     ***
      LOGICAL FUNCTION INITIALIZE(P_COS, FNAME)
      IMPLICIT NONE
      REAL P_COS(*)
      CHARACTER*(*) FNAME

C
C        Common variables and parameters used througout the program.
C
      INCLUDE 'cosmo.inc'
      INCLUDE 'operation.inc'
      INCLUDE 'supernova.inc'
      INCLUDE 'io.inc'
      INCLUDE 'probability.inc'
      INCLUDE 'error.inc'
      INCLUDE 'mate.inc'
      INCLUDE 'snalys.inc'
      INCLUDE 'gausslin.inc'

C
C        Local variables.
C
      CHARACTER*(256) TMP, INIFILE

      CHARACTER*(2) MO_STR, DA_STR

      INTEGER N, I, LENGTH, INIID, IFID, IOS, RMPAR, J, K
      PARAMETER (INIID = 30, IFID = 20 )

      REAL C_START_GET, C_END_GET, C_NR_GET, GETVAL
      REAL C_CORR_GET, ISEED

      LOGICAL READINIT, FOUND(ENR), VALIDDATA
      CHARACTER*9 NAME_NOBLANK 
      INTEGER NBLANK
      REAL IBAND,JBAND

C
C        Set the sigma level to show...
C
      NRSIG = 1

C
C        Set the seed.
C
      CALL RNCLCK( ISEED )
      CALL RNINIT( ISEED )

C
C        Initialize the boolean operation array (what is to be done with
C        the data) and the array that determines what parameters that are
C        to be plotted.
C
      DATA OPER/NR_OP*.TRUE./
      
C
C        Initialize unsigned variables from the common blocks defined in
C        'snalys.inc'
C
      FNAME = ' '
      WCUBEFILE = ' '
      WHICH_PDF = P_GAUSS          ! The default pdf.
      NR_SUPER = -1
      SN_USED = 0
      GAUSS_SIGMA = 0.17
      OPER(QUIET) = .FALSE.        ! Default is not quiet-mode
      OPER(DO_POWELL) = .FALSE.    ! Do not do a powell min by default.
      OPER(DO_GENETIC) = .FALSE.   ! Do not do a genetic min by default.
      OPER(SUB_INTSIG) = .FALSE.   ! Do not subtract the intrinsic dev.
      OPER(MAG_DIST) = .FALSE.     ! Do not produce a magnitude distribution.
      OPER(FLAT) = .FALSE.         ! Make no assumption of the geometry.
      OPER(ADD_LDM) = .FALSE.      ! Do not care about lensing
      OPER(ADD_GDM) = .FALSE.      ! Do not care about grey dust
      OPER(ADD_HOSTDM) = .FALSE.   ! Do not care about host extinction
      OPER(ADD_MIXDM) = .FALSE.    ! Do not care about photon-axion osc.
      OPER(DO_POWELL) = .FALSE.
      OPER(DO_GENETIC) = .FALSE.
      OPER(DO_DAVIDON) = .FALSE.
      OPER(CON_FIND) = .FALSE.
      OPER(ZLIMIT) = .FALSE.
      OPER(MASS_PRIOR) = .FALSE.
      OPER(MATE_EQST) = .FALSE.    ! Do not fit the fancy quintessende mod.
      OPER(FITMSC) = .FALSE.       ! Do not minimize ScriptM
      OPER(FITSALP) = .FALSE.      ! Do not minimize stretch alpha
      OPER(DO_MULTI) = .FALSE.     ! Do not fit multiple images
      OPER(XDIM) = .FALSE.         ! Do not fit X-dimensions
      OPER(FITS) = .FALSE.         ! Do not save output as fits for default.
      OPER(DO_JACKE) = .FALSE.     ! Do not fit strange integral.
      OPER(DO_LINDER) = .FALSE.    ! Do not fit Linder parametrization.
      OPER(DO_WILT) = .FALSE.      ! Do not fit Wiltshare's cosmology
      OPER(DO_HP) = .FALSE.        ! Do not use Hannestad parametrization of EO
      OPER(DO_BOS) = .FALSE.       ! Do not use baryon oscillation prior
      OPER(Z_ERROR) = .FALSE.      ! Do not add redshift dependent error.
      OPER(DO_CMB) = .FALSE.       ! Do not use CMB prior
      OPER(DO_OM) = .FALSE.         ! Do not use Omega_M prior

      OFID = -1

C
C        Reset the cosmology index variables
C
      DATA O_M/-1/,O_X/-1/,A_X/-1/,A_2/-1/,M_SC/-1/,ALP/-1/,SALP/-1/,
     $     O_RC/-1/,B_1/-1/,B_2/-1/,B_inv/-1/,B_4/-1/,W_0/-1/,W_1/-1/,
     $     Q_HP/-1/,Z_X/-1/,X_DA/-1/,L_1/-1/,L_2/-1/

C
C        Read the argument list
C
      INITIALIZE = .FALSE.
      N = 1
      CALL GETARG(N,TMP)
      DO WHILE ( INDEX(TMP,' ').NE.1 )
         IF ( TMP.EQ.'-i' ) THEN
            N = N + 1
            CALL GETARG(N,TMP)
C
C              If a name is given as a second argument, then this is the
C              file name.
C
            IF ((INDEX(TMP,' ').NE.1) .AND. (INDEX(TMP,'-').NE.1)) THEN
               OPEN(INIID,FILE=TMP,IOSTAT=IOS,STATUS='OLD')
               IF (IOS.EQ.0) THEN
                  INITIALIZE = READINIT(INIID) ! Read the init file
               ENDIF
            ENDIF
C
C              If the file was not opened successfully or if there was
C              no argument after '-i' an error message will be written.
C
            IF ( .NOT. INITIALIZE ) THEN
               WRITE (STDERR,'(3A)') 'initialize: unable to open the ',
     $              'parameter file ',TMP
            ENDIF
         ELSEIF ( TMP.EQ.'-o' ) THEN ! The output file
            N = N + 1
            CALL GETARG( N, OUTFILE )
C  
C              If a name is given as a second argument, then the file
C              is opened, otherwise the help dialog is printed.
C  
            IF ( (INDEX( OUTFILE,' ').NE.1) .AND. 
     $           (INDEX(OUTFILE,'-').NE.1 ) ) THEN
               OFID = 40
            ELSE
               WRITE (STDERR,'(2A)') 'initialize: no output file ',
     $              'specified after -o'
               INITIALIZE = .FALSE.
            ENDIF
         ELSEIF ( TMP.EQ.'-q' ) THEN ! Quiet mode
            OPER(QUIET) = .TRUE.
         ELSEIF ( TMP.EQ.'-f' ) THEN ! Output as fits
            OPER(FITS) = .TRUE.
C
C           Make sure that only one data file can be specified
C
         ELSEIF ( INDEX(TMP,'-').NE.1 .AND. INDEX(TMP,' ').NE.1 ) THEN
            FNAME = TMP
         ELSE
            WRITE (STDERR,'(2A)') 'initialize: incorrect argument list',
     $           ' the data file does not exist!'
            INITIALIZE = .FALSE.
         ENDIF

         N = N + 1
         CALL GETARG(N,TMP)     ! Get the N:th argument from the argument list.
      ENDDO

C
C        If no output file was specified use the basename of the
C        data file and add .sna.
C
      IF ( INITIALIZE .AND. OFID.EQ.-1 ) THEN
         N = LEN(FNAME)
         I = 0
         J = 0
         DO WHILE ( INDEX(FNAME(I+1:N),'.').NE.0 )
            I = INDEX(FNAME(I+1:N),'.') + I
            IF ( I.GT.0 ) J = I
         ENDDO
         OFID = 40
         IF ( OPER(FITS) ) THEN
            OUTFILE(1:J+4) = FNAME(1:J) // 'fits'
         ELSE
            OUTFILE(1:J+3) = FNAME(1:J) // 'sna'
         ENDIF
      ENDIF

      IF ( OPER(MATE_EQST) .AND. INDEX(WCUBEFILE,' ').NE.1 ) THEN
         INQUIRE( FILE=WCUBEFILE, EXIST=OPER(MATE_EQST) )
         IF ( OPER(MATE_EQST) ) THEN
C
C              Initialize GLOB_WINDEX with variable values and
C              GLOB_WCUBE with w values
C
            OPEN(66, FILE=WCUBEFILE, STATUS='OLD')
            READ(66, '(4A)')    !Get rid of first row
            DO I=1, GLOB_N_ALPHA
               DO J=1, GLOB_N_OMM
                  DO K=1, GLOB_N_ZP 
                     READ(66, '(4E16.6)') GLOB_WINDEX(I),
     $                    GLOB_WINDEX(GLOB_N_ALPHA + J),
     $                    GLOB_WINDEX(GLOB_N_ALPHA + GLOB_N_OMM + K),
     $                    GLOB_WCUBE(I, J, K)
                  ENDDO
               ENDDO
            ENDDO
            CLOSE(66)
C
C              If the inverse power law quintessence fit is going to
C              be made, we can not fit W0 and W1 at the same time.
C     
            IF ( A_X .GE. 1 ) THEN
               A_X = RMPAR( A_X )
               WRITE (STDERR,'(2A)') 'initialize: W0 can not be ',
     $              'fitted together with power-alpha.'
            ELSEIF ( A_2 .GE. 1 ) THEN
               A_2 = RMPAR( A_2 )
               WRITE (STDERR,'(2A)') 'initialize: W1 can not be ',
     $              'fitted together with power-alpha.'
            ENDIF


C
C              No minimization is possible with Davidon's method for the
C              inverse power-law quintessence model since no derivatives
C              can be calculated.
C
            IF ( OPER(DO_DAVIDON).OR.CALCDERIV ) THEN
               OPER(DO_DAVIDON) = .FALSE.
               CALCDERIV = .FALSE.
               WRITE (STDERR,'(2A)') 'initialize: davidon ',
     $              'minimization can not be done and is disabled'
            ENDIF
         ELSE
            WRITE (STDERR,'(2A)') 'initialize: unable to open ',
     $           'WCUBEFILE, disabling QALPHA fit'
         ENDIF
      ENDIF

C  
C        In case a flat universe is choosen and Omega_X has been 
C        defined as an parameter to be estimated, it must be
C        removed from the parameter list.
C  
      IF ( OPER(FLAT) .AND. O_X.GT.0 ) THEN
         O_X = RMPAR( O_X )
      ENDIF

      IF ( OPER(XDIM) ) THEN
        IF ( O_X.GT.0) O_X = RMPAR( O_X )
        IF ( OPER(FLAT) .AND. O_RC.GT.0 ) THEN
           O_RC = RMPAR( O_RC )   
        ENDIF
      ELSE
         IF ( O_RC.GT.0 ) THEN
            O_RC = RMPAR( O_RC )
         ELSEIF ( X_DA .GT.0) THEN
            X_DA = RMPAR( X_DA )
         ENDIF
      ENDIF

      IF ( .NOT. OPER(DO_JACKE) ) THEN
         IF ( B_1 .GT. 0 ) THEN
            B_1 = RMPAR( B_1 )
         ENDIF
         IF ( B_2 .GT. 0 ) THEN
            B_2 = RMPAR( B_2 )
         ENDIF
         IF ( B_inv .GT. 0 ) THEN
            B_inv = RMPAR( B_inv )
         ENDIF
         IF ( B_4 .GT. 0 ) THEN
            B_4 = RMPAR( B_4 )
         ENDIF
      ENDIF

      IF ( .NOT. OPER(DO_LINDER) ) THEN
         IF ( L_1 .GT. 0 ) THEN
            L_1 = RMPAR( L_1 )
         ENDIF
         IF ( L_2 .GT. 0 ) THEN
            L_2 = RMPAR( L_2 )
         ENDIF
      ENDIF


      IF ( .NOT. OPER(DO_WILT) ) THEN
         IF ( FV_0 .GT. 0 ) THEN
            FV_0 = RMPAR( FV_0 )
         ENDIF
         IF ( H_0BAR .GT. 0 ) THEN
            H_0BAR = RMPAR( H_0BAR )
         ENDIF
      ELSE
         IF (O_M .GT. 0) THEN
            O_M = RMPAR(O_M)
         ENDIF
         IF (O_X .GT. 0) THEN
            O_X = RMPAR(O_X)
         ENDIF
      ENDIF


      IF ( .NOT. OPER(DO_HP) ) THEN
         IF ( W_0 .GT. 0 ) THEN
            W_0 = RMPAR( W_0 )
         ENDIF
         IF ( W_1 .GT. 0 ) THEN
            W_1 = RMPAR( W_1 )
         ENDIF
         IF ( Q_HP .GT. 0 ) THEN
            Q_HP = RMPAR( Q_HP )
         ENDIF
         IF ( Z_X .GT. 0 ) THEN
            Z_X = RMPAR( Z_X )
         ENDIF
      ENDIF

      IF ( .NOT. OPER(DO_MULTI) .AND. H_0 .GT. 0 ) THEN
         H_0 = RMPAR( H_0 )
      ENDIF

C  
C         Does the product of the number of steps exceed the maximal
C         length of the pdf-array produced by the grid search.
C  
      LENGTH = 1
      DO N = 1, NR_VAR
         LENGTH = LENGTH*INT(C_NR_GET(N))
      ENDDO
      IF ( LENGTH.GT.MAX_LEN ) THEN
         WRITE (STDERR,'(2A)') 'initialize: to many steps to run ',
     $        ' a grid search, choose a smaller grid'
         OPER(DO_GRID) = .FALSE.
      ENDIF

C
C        Check if the input filename has been set.
C
      IF ( INDEX(FNAME,' ').NE.1  ) THEN
C  
C          Open the input data file.
C
         OPEN (IFID, FILE=FNAME, STATUS='OLD', IOSTAT=IOS)
C  
C          Check if the file was opened successfully.
C  
         IF (IOS.NE.0) INITIALIZE = .FALSE.
      ELSE
         INITIALIZE = .FALSE.
      ENDIF
      

      IF ( INITIALIZE ) THEN
         SN_NR = 0              ! the number of supernovae

C  
C           Read the data file.
C     
         IF ( OPER(DO_GRID) .OR. OPER(DO_GENETIC) .OR. 
     $       OPER(DO_POWELL) .OR. OPER(DO_DAVIDON) .OR.
     $       OPER(CON_FIND) ) THEN
            IF ( .NOT. OPER(QUIET) )
     $           WRITE (STDOUT,'(A)') 'Reading datafile...'
C           
            DO WHILE (.TRUE.)
C  
C                 If not all of the supernovae are to be considered.
C  
               IF ( ((NR_SUPER.NE.-1).AND.(SN_NR.GE.NR_SUPER))
     $              .OR.(SN_NR.GT.MAX_SN) ) THEN
                  IF ( SN_NR.GT.MAX_SN ) THEN
                     WRITE (STDERR,'(2A)') 'initialize: the maximum ',
     $                    'number of allowed sne read'
                  ENDIF
                  GOTO 25
               ENDIF
               
C  
C                 Reset the FOUND array
C  
               DO N = 1,ENR
                  FOUND(N) = .FALSE.
               ENDDO
C  
C                 Find the beginning of the event.
C  
               DO WHILE ( INDEX(TMP,'event').NE.1 )
                  READ ( IFID,'(A)',END = 25 ) TMP
               ENDDO
               SN_NR = SN_NR + 1
C  
C                 Read each string until the end string and
C                 try to indentify them.
C  
               DO WHILE ( INDEX(TMP,'end').NE.1 )
                  READ ( IFID,'(A)',END = 25 ) TMP
                  DO N = 1,ENR
                     NAME_NOBLANK = ENAMES(N)
                     NBLANK =  INDEX(NAME_NOBLANK,' ')
C                     WRITE (*,*) TMP, NAME_NOBLANK
                     IF (INDEX(TMP,NAME_NOBLANK(1:NBLANK-1)).EQ.1) THEN
                        SN_DATA(N,SN_NR) = GETVAL(TMP)
                        FOUND(N) = .TRUE.
C                        WRITE (*,*) N, ' ',  ENAMES(N), SN_DATA(N,SN_NR)
                     ENDIF
                  ENDDO
               ENDDO


C                 I am not sure what this is, if the redshift is < 0.3,
C                 set the lensing to zero. Why hardcode 1 and 4? What
C                 did I think, or perhaps I am not the one who wrote this?
C
C               IF ( SN_DATA(1,SN_NR) .LT. 0.3 ) THEN
C                  SN_DATA(4,SN_NR) = 0.0
C               ENDIF

C     
C                 Check that that the redshift (zs) and the magnitude
C                 (bzmag) were identified for this supernova. If not,
C                 skip this one and move on to the next in the data file.
C  
C                 If any of the non-vital parameters were not identified
C                 set them to 0.0 in the SN_DATA array.
C   
               VALIDDATA = .TRUE.
               DO N = 1,ENR
                  IF ( .NOT.FOUND(N) ) THEN
                     IF ( OPER(DO_MULTI) ) THEN
                        IF ( N.EQ.TDEL .OR. N.EQ.IMSEP .OR. N.EQ.RATIO
     $                       .OR. N.EQ.ZL .OR. N.EQ.VEL ) THEN
                           VALIDDATA = .FALSE.
                           WRITE (STDERR,'(3A)') 'initialize: midel, ',
     $                          'imsep, secrat and vlmi and zlmi must ',
     $                          'be specified for each SN'
                        ENDIF
                     ELSE
                        IF ( N.EQ.ZS .OR. N.EQ.BZMAG) THEN
                           VALIDDATA = .FALSE.
                           WRITE (STDERR,'(2A)') 'initialize: zs and ',
     $                          'bzmag must be given for each SN'
                        ELSE
                           IF ( N.EQ.STRETCH .AND. SALP.GT.-1 ) THEN
c                              WRITE (STDERR,'(2A)') 'initialize: ',
c     $                          'stretch not specified for SN using 1.0'
                              SN_DATA(N,SN_NR) = 1.0
                           ELSEIF ( N.EQ.LDM .AND. OPER(ADD_LDM) ) THEN
                              WRITE (STDERR,'(2A)') 'initialize: ',
     $                          'linsdm not specified for SN, using 0.0'
                              SN_DATA(N,SN_NR) = 0.0
                           ELSEIF ( N.EQ.GDM .AND. OPER(ADD_GDM) ) THEN
                              WRITE (STDERR,'(2A)') 'initialize: ',
     $                          'greydm not specified for SN, using 0.0'
                              SN_DATA(N,SN_NR) = 0.0
                           ELSEIF ( N.EQ.HOSTDM .AND.
     $                             OPER(ADD_HOSTDM) ) THEN
                              WRITE (STDERR,'(2A)') 'initialize: ',
     $                           'hostdm not specified for SN using 0.0'
                              SN_DATA(N,SN_NR) = 0.0
                           ELSEIF ( N.EQ.DHOSTDM .AND.
     $                             OPER(ADD_HOSTDM) ) THEN
                              WRITE (STDERR,'(3A)') 'initialize: ',
     $                             'dhostdm not specified for SN, ',
     $                             'using 0.0'
                              SN_DATA(N,SN_NR) = 0.0
                           ELSEIF ( N.EQ.MDM .AND.
     $                             OPER(ADD_MIXDM) ) THEN
                              WRITE (STDERR,'(3A)') 'initialize: ',
     $                             'mixdm not specified for SN, ',
     $                             'using 0.0'
                              SN_DATA(N,SN_NR) = 0.0
                           ELSE IF ( N.EQ.ISG .AND. OPER(SUB_INTSIG) ) 
     $                             THEN
                              WRITE (STDERR,'(3A)') 'initialize: ',
     $                             'intsig not specified for SN, ',
     $                             'using 0.0'
                              SN_DATA(N,SN_NR) = 0.0
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
                  IF ( SN_DATA(ZS,SN_NR).LE.0.0 .OR.
     $                 SN_DATA(MERR,SN_NR).LT.0.0 .OR.
     $                 SN_DATA(ZERR,SN_NR).LT.0.0 ) THEN
                     VALIDDATA= .FALSE.
                     WRITE (STDERR,'(A)')
     $                   'initialize: this supernova will be skipped...'
                  ENDIF
               ENDDO
               IF ( .NOT. VALIDDATA ) WRITE (*,*) 'initialize: ',
     $              'SN has not enough input data and will not be used'

C
C                 Only consider SNe that are within a certain
C                 redshift range, if that option is enabled
C
               IF ( OPER(ZLIMIT) .AND. 
     $              ( SN_DATA(ZS,SN_NR).LT.ZRANGE(ZLOW)
     $                .OR. SN_DATA(ZS,SN_NR).GT.ZRANGE(ZHIG) ) ) THEN
                  VALIDDATA = .FALSE.
               ENDIF

C
C                 Only consider SNe that whose SECONDARY images are
C                 within the I-band magnitude limit set by IBANDMAX
C                 and JBANDMAX
C
               IF ( OPER(DO_MULTI) ) THEN
                  IBAND = SN_DATA(BZMAG,SN_NR)+SN_DATA(KCORBI,SN_NR) +
     $                 SN_DATA(LDM,SN_NR) + 
     $                 2.5*ALOG(SN_DATA(RATIO,SN_NR)/
     $                 SN_DATA(LENSRAT,SN_NR))
                  JBAND = SN_DATA(BZMAG,SN_NR)+SN_DATA(KCORBJ,SN_NR) +
     $                 SN_DATA(LDM,SN_NR) + 
     $                 2.5*ALOG(SN_DATA(RATIO,SN_NR)/
     $                 SN_DATA(LENSRAT,SN_NR))
                  
                  IF ( IBAND.GT.IBANDMAX .AND.
     $                 JBAND.GT.JBANDMAX ) THEN
                     VALIDDATA = .FALSE.
                  ENDIF
               ENDIF
C     
C              If some of the data is not in the proper ranges
C              that supernova is excluded.
C  
               IF ( .NOT. VALIDDATA ) SN_NR = SN_NR - 1
            ENDDO

 25         CONTINUE
            IF ( .NOT. OPER(QUIET) ) WRITE (STDOUT,'(A)' ) 'done.'
            CLOSE (IFID)
            SN_USED = SN_NR

C
C              Pre-calculate some values for each supernova that do
C              not depend on cosmology.
C
            CALL NONCOSDEP( SN_USED , GAUSS_SIGMA , COSFP )
         ENDIF
C  
C           Initialize the p_cos
C  
         P_COS(NR_VAR+1) = 1.E32

C  
C           If there is not going to be a grid search, start values for
C           the davidon, powell and genetic minimization routines must be choosen.
C     
         IF ( .NOT.OPER(DO_GRID) ) THEN
            CALL RANDCOSMO(P_COS)
         ENDIF
      ENDIF
      RETURN
      END
C
C   End of INITIALIZE
C



C     ****f* snalys/READINIT *
C
C   NAME
C     READINIT -- reading the init file
C
C   DESCRIPTION
C     Reading the init file where the conditions for the fit is given.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-10-25
C
C   USAGE
C     READINIT( INIID )
C
C   INPUTS
C     INIID   - a valid unit id of type INTEGER to the ini file that
C               is going to be read.
C
C   RESULT
C     A LOGICAL value is returned to confirm that reading of the file
C     was successful.
C
C   SIDE EFFECTS
C     COMMON variables from operation.inc, supernova.inc, snalys.inc,
C     cosmo.inc are set and probably a lot more.
C
C   ERRORS
C     Error messages will be written if any of the parmeters could
C     not be read from the init file.
C
C   HISTORY
C     2001-09-07 Added the option of also estimating the quintessence
C                parameter alpha in the model implemented by Martin
C                Eriksson.
C     2003-03-31 Added the option of fitting the stretch alpha.
C     2003-07-09 Rewrote the format for the parameter file.
C                All keys are given as capitals, and the boolean
C                variables are specified using 1 and 0. Also a
C                function was added to read the REAL value after
C                a specified key.
C                All plotting options and old options were removed.
C                Removed some obselete keywords.
C     2003-07-10 Whenever the nuisance parameters are going to be
C                be fitted they are put in the end of the vector
C                of parameters to be fitted, and actually NR_VAR
C                is decreased, so they will not be part of the grid
C                search in the ordinary manner.
C     2004-01-15 Added the keys realted to "jacke"-cosmology. An alternative
C                parametrization of dark energy using the variables B1 and B2.
C     2004-02-26 Introduced the parameter MAX_LINE_LENGTH
C     2004-05-17 Added the keys realted to the Hannestad parametrization.
C     2004-08-08 Added a condition to only check wether HP and JACKE
C                parameters have been initialized if these paramters
C                are set.
C     2005-02-16 Added HOSTDM, MIXDM
C
C     ***
      LOGICAL FUNCTION READINIT(INIID)
      IMPLICIT NONE
      INTEGER INIID

      INCLUDE 'cosmo.inc'
      INCLUDE 'probability.inc'
      INCLUDE 'operation.inc'
      INCLUDE 'supernova.inc'
      INCLUDE 'error.inc'
      INCLUDE 'mate.inc'
      
      INTEGER LINE_LENGTH, MAX_LINE_LENGTH, N, I, K
      PARAMETER (MAX_LINE_LENGTH = 100)
      CHARACTER*(MAX_LINE_LENGTH) TMP
      CHARACTER*30 TMPLABEL
      REAL TMPVALUE
      REAL C_START_GET, C_NR_GET, C_END_GET, CSTART, CNR, CEND
      LOGICAL RVARVAL, T, GETONEVAL

      READINIT = .TRUE.

C
C        Initialize N
C
      N = 0

C  
C        The program will jump out of this while loop at the end
C        of the file. I can not find any other way to solve this
C        in fortran 77.
C  
      NR_VAR = 0
      DO WHILE (.TRUE.)
         READ (INIID,'(A)', END = 100) TMP ! 100 is the end of the while
C  
C           Search for comment characters (#).
C  
         LINE_LENGTH = INDEX(TMP,'#') - 1
         IF ( LINE_LENGTH .EQ. -1 ) LINE_LENGTH = MAX_LINE_LENGTH
C  
C           Find the label
C           
         I = 1
         CALL FINDNEXTWORD(TMP,I,N,LINE_LENGTH)

C  
C           Read the line if this was not an empty line.
C  
         IF ( I.LT.N ) THEN
            READ ( TMP(I:N),'(A)' ) TMPLABEL
C     
C             Identify the label.
C  
            IF ( INDEX(TMPLABEL, 'W0').GT.0 ) THEN
               T = RVARVAL(TMP,A_X,GLOB_AX,I,N,LINE_LENGTH)
            ELSEIF ( INDEX(TMPLABEL,'W1').GT.0 ) THEN
               T = RVARVAL(TMP,A_2,GLOB_A2,I,N,LINE_LENGTH)
            ELSEIF ( INDEX(TMPLABEL,'QALPHA').GT.0 ) THEN
               T = RVARVAL(TMP,ALP,GLOB_ALPHA,I,N,LINE_LENGTH)
               OPER(MATE_EQST) = .TRUE.
            ELSEIF ( INDEX(TMPLABEL,'ALPHA').GT.0 ) THEN
               T = RVARVAL(TMP,SALP,GLOB_SALP,I,N,LINE_LENGTH)
            ELSEIF ( INDEX(TMPLABEL,'OMEGAX').GT.0 ) THEN
               T = RVARVAL(TMP,O_X,GLOB_OMX,I,N,LINE_LENGTH)
            ELSEIF ( INDEX(TMPLABEL,'OMEGAM').GT.0 ) THEN
               T = RVARVAL(TMP,O_M,GLOB_OMM,I,N,LINE_LENGTH)
            ELSEIF ( INDEX(TMPLABEL,'SCRIPTM').GT.0 ) THEN
               T = RVARVAL(TMP,M_SC,GLOB_MSCRIPT,I,N,
     $              LINE_LENGTH)
            ELSEIF ( INDEX(TMPLABEL,'HUBBLE0').GT.0 ) THEN
               T = RVARVAL(TMP,H_0,GLOB_H0,I,N,
     $              LINE_LENGTH)
            ELSEIF ( INDEX(TMPLABEL,'OMEGARC').GT.0 ) THEN
               T = RVARVAL(TMP,O_RC,GLOB_ORC,I,N,
     $              LINE_LENGTH)
            ELSEIF ( INDEX(TMPLABEL,'XDA').GT.0 ) THEN
               T = RVARVAL(TMP,X_DA,GLOB_XDA,I,N,
     $              LINE_LENGTH)
            ELSEIF ( INDEX(TMPLABEL,'A1').GT.0 ) THEN
               T = RVARVAL(TMP,B_1,GLOB_B1,I,N,LINE_LENGTH)
            ELSEIF ( INDEX(TMPLABEL,'A2').GT.0 ) THEN
               T = RVARVAL(TMP,B_2,GLOB_B2,I,N,LINE_LENGTH)
            ELSEIF ( INDEX(TMPLABEL,'AINV').GT.0 ) THEN
               T = RVARVAL(TMP,B_inv,GLOB_Binv,I,N,LINE_LENGTH)
            ELSEIF ( INDEX(TMPLABEL,'A4').GT.0 ) THEN
               T = RVARVAL(TMP,B_4,GLOB_B4,I,N,LINE_LENGTH)
            ELSEIF ( INDEX(TMPLABEL,'L1').GT.0 ) THEN
               T = RVARVAL(TMP,L_1,GLOB_L1,I,N,LINE_LENGTH)
            ELSEIF ( INDEX(TMPLABEL,'L2').GT.0 ) THEN
               T = RVARVAL(TMP,L_2,GLOB_L2,I,N,LINE_LENGTH)
            ELSEIF ( INDEX(TMPLABEL,'FV0').GT.0 ) THEN
               T = RVARVAL(TMP,FV_0,GLOB_FV0,I,N,LINE_LENGTH)
            ELSEIF ( INDEX(TMPLABEL,'H0BAR').GT.0 ) THEN
               T = RVARVAL(TMP,H_0BAR,GLOB_H0BAR,I,N,LINE_LENGTH)
            ELSEIF ( INDEX(TMPLABEL,'WLT').GT.0 ) THEN
               T = RVARVAL(TMP,W_0,GLOB_W0,I,N,LINE_LENGTH)
            ELSEIF ( INDEX(TMPLABEL,'WET').GT.0 ) THEN
               T = RVARVAL(TMP,W_1,GLOB_W1,I,N,LINE_LENGTH)
            ELSEIF ( INDEX(TMPLABEL,'QHP').GT.0 ) THEN
               T = RVARVAL(TMP,Q_HP,GLOB_QHP,I,N,LINE_LENGTH)
            ELSEIF ( INDEX(TMPLABEL,'ZX').GT.0 ) THEN
               T = RVARVAL(TMP,Z_X,GLOB_ZX,I,N,LINE_LENGTH)
            ELSEIF ( INDEX(TMPLABEL,'ZRANGE').GT.0 ) THEN
               OPER(ZLIMIT) = .TRUE.
               
               K = 0
               DO WHILE ( OPER(ZLIMIT) .AND. K.LT.ZRGPAR )
                  K = K + 1
                  I = N + 1
                  CALL FINDNEXTWORD(TMP,I,N,LINE_LENGTH)
                  
                  ZRANGE(K) = -30000.
                  IF ( I.LE.N ) READ (TMP(I:N),'(F13.6)', END = 15)
     $                 ZRANGE(K)
 15               CONTINUE
                  IF ( ZRANGE(K).EQ.-30000. ) THEN
                     OPER(ZLIMIT) = .FALSE.
                     WRITE (STDERR,'(3A)') 'readinit: error in ', 
     $                    'reading zrange from the ini-file, will not ',
     $                    'use any limits'
                  ENDIF
               ENDDO
               
               IF ( OPER(ZLIMIT) ) THEN
                  IF ( ZRANGE(ZLOW).GT.ZRANGE(ZHIG) ) THEN
                     TMPVALUE = ZRANGE(ZLOW)
                     ZRANGE(ZLOW) = ZRANGE(ZHIG)
                     ZRANGE(ZHIG) = TMPVALUE
                  ENDIF
               ENDIF
            ELSEIF ( INDEX(TMPLABEL,'GAUSSERROR').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(GAUSS_SIGMA,I,TMP,LINE_LENGTH) ) THEN
                  GAUSS_SIGMA = .17
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'gausserror from the parameter file'
               ENDIF
            ELSEIF ( INDEX(TMPLABEL,'USEMERROR').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(TMPVALUE,I,TMP,LINE_LENGTH) ) THEN
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'XDIM from the parameter file'
               ELSE
                  IF ( TMPVALUE.EQ.1. ) THEN
                     OPER(USEMERR) = .TRUE.
                  ELSE
                     OPER(USEMERR) = .FALSE.
                  ENDIF
               ENDIF
            ELSEIF ( INDEX(TMPLABEL,'ETIMEDELAY').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(MULTI_ETDEL,I,TMP,LINE_LENGTH) ) THEN
                  MULTI_ETDEL = .1
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'error in time delay from the parameter file'
               ENDIF
            ELSEIF ( INDEX(TMPLABEL,'EIMSEP').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(MULTI_EIMSEP,I,TMP,LINE_LENGTH) )THEN
                  MULTI_EIMSEP = .1
                  WRITE (STDERR,'(3A)') 'readinit: error in reading ',
     $                 'error in image separation from the parameter ',
     $                 'file'
               ENDIF
            ELSEIF ( INDEX(TMPLABEL,'EMODEL').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(MULTI_EMODEL,I,TMP,LINE_LENGTH) )THEN
                  MULTI_EMODEL = 0.
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $               'error in lens model error from the parameter file'
               ENDIF
            ELSEIF ( INDEX(TMPLABEL,'EIMRAT').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(MULTI_EIMRAT,I,TMP,LINE_LENGTH) )THEN
                  MULTI_EIMRAT = .1
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'error in image ratio from the parameter file'
               ENDIF
            ELSEIF ( INDEX(TMPLABEL,'EZL').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(MULTI_EZL,I,TMP,LINE_LENGTH) )THEN
                  MULTI_EZL = .0
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'error in ZL from the parameter file'
               ENDIF
            ELSEIF ( INDEX(TMPLABEL,'EVEL').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(MULTI_EVEL,I,TMP,LINE_LENGTH) )THEN
                  MULTI_EVEL = .0
                  WRITE (STDERR,'(3A)') 'readinit: error in reading ',
     $                 'error in dispersion velocity from the ',
     $                 'parameter file'
               ENDIF
            ELSEIF ( INDEX(TMPLABEL,'IBANDMAX').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(IBANDMAX,I,TMP,LINE_LENGTH) )THEN
                  IBANDMAX = 99.
                  WRITE (STDERR,'(2A)') 'readinit: error in reading',
     $                 ' ibandmax from the parameter file'
               ENDIF
            ELSEIF ( INDEX(TMPLABEL,'JBANDMAX').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(IBANDMAX,I,TMP,LINE_LENGTH) )THEN
                  IBANDMAX = 99.
                  WRITE (STDERR,'(2A)') 'readinit: error in reading',
     $                 ' ibandmax from the parameter file'
               ENDIF
            ELSEIF ( INDEX(TMPLABEL,'LENSTYPE').GT.0 ) THEN
               I = N + 1
               CALL FINDNEXTWORD(TMP,I,N,LINE_LENGTH)

               LENSTYPE = 'XXX'
               IF ( I.LE.N ) READ (TMP(I:N),'(A3)', END = 24)
     $              LENSTYPE
 24            CONTINUE
               IF (LENSTYPE .EQ. 'XXX') THEN
                  LENSTYPE = 'SIS'
                  WRITE (STDERR,'(2A)') 'readinit: error in reading',
     $                 '  lenstype from the parameter file'
               ENDIF
            ELSEIF ( INDEX(TMPLABEL,'FLAT').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(TMPVALUE,I,TMP,LINE_LENGTH) ) THEN
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'FLAT from the parameter file'
               ELSE
                  IF ( TMPVALUE.EQ.0. ) THEN
                     OPER(FLAT) = .FALSE.
                  ELSE
                     OPER(FLAT) = .TRUE.
                  ENDIF
               ENDIF
            ELSEIF ( INDEX(TMPLABEL,'XDIM').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(TMPVALUE,I,TMP,LINE_LENGTH) ) THEN
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'XDIM from the parameter file'
               ELSE
                  IF ( TMPVALUE.EQ.1. ) THEN
                     OPER(XDIM) = .TRUE.
                  ELSE
                     OPER(XDIM) = .FALSE.
                  ENDIF
               ENDIF
            ELSEIF ( INDEX(TMPLABEL,'DOMULTI').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(TMPVALUE,I,TMP,LINE_LENGTH) ) THEN
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'DOMULTI from the parameter file'
               ELSE
                  IF ( TMPVALUE.EQ.1. ) THEN
                     OPER(DO_MULTI) = .TRUE.
                  ELSE
                     OPER(DO_MULTI) = .FALSE.
                  ENDIF
               ENDIF               
            ELSEIF ( INDEX(TMPLABEL,'DOJACKE').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(TMPVALUE,I,TMP,LINE_LENGTH) ) THEN
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'DOJACKE from the parameter file'
               ELSE
                  IF ( TMPVALUE.EQ.1. ) THEN
                     OPER(DO_JACKE) = .TRUE.
                  ELSE
                     OPER(DO_JACKE) = .FALSE.
                  ENDIF
               ENDIF
               
            ELSEIF ( INDEX(TMPLABEL,'DOLINDER').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(TMPVALUE,I,TMP,LINE_LENGTH) ) THEN
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'DOLINDER from the parameter file'
               ELSE
                  IF ( TMPVALUE.EQ.1. ) THEN
                     OPER(DO_LINDER) = .TRUE.
                  ELSE
                     OPER(DO_LINDER) = .FALSE.
                  ENDIF
               ENDIF   

            ELSEIF ( INDEX(TMPLABEL,'DOWILTSHIRE').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(TMPVALUE,I,TMP,LINE_LENGTH) ) THEN
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'DOWILTSHIRE from the parameter file'
               ELSE
                  IF ( TMPVALUE.EQ.1. ) THEN
                     OPER(DO_WILT) = .TRUE.
                  ELSE
                     OPER(DO_WILT) = .FALSE.
                  ENDIF
               ENDIF   

            ELSEIF ( INDEX(TMPLABEL,'DOBOS').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(TMPVALUE,I,TMP,LINE_LENGTH) ) THEN
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'DOBOS from the parameter file'
               ELSE
                  IF ( TMPVALUE.EQ.1. ) THEN
                     OPER(DO_BOS) = .TRUE.
                  ELSE
                     OPER(DO_BOS) = .FALSE.
                  ENDIF
               ENDIF   
            ELSEIF ( INDEX(TMPLABEL,'BOSPRI').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(BOSPRI,I,TMP,LINE_LENGTH) )THEN
                  BOSPRI = .469
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'error in the bos prior from the parameter file'
               ENDIF
            ELSEIF ( INDEX(TMPLABEL,'DBOSPR').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(DBOSPR,I,TMP,LINE_LENGTH) )THEN
                  DBOSPR = .017
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'error in bos prior uncert. the parameter file'
               ENDIF

            ELSEIF ( INDEX(TMPLABEL,'DOCMB').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(TMPVALUE,I,TMP,LINE_LENGTH) ) THEN
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'DOCMB from the parameter file'
               ELSE
                  IF ( TMPVALUE.EQ.1. ) THEN
                     OPER(DO_CMB) = .TRUE.
                  ELSE
                     OPER(DO_CMB) = .FALSE.
                  ENDIF
               ENDIF   

            ELSEIF ( INDEX(TMPLABEL,'CMBPRI').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(CMBPRI,I,TMP,LINE_LENGTH) )THEN
                  CMBPRI = 1.72
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'error in the CMB prior from the parameter file'
               ENDIF
            ELSEIF ( INDEX(TMPLABEL,'DCMBPR').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(DCMBPR,I,TMP,LINE_LENGTH) )THEN
                  DCMBPR = .031
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'error in CMB prior uncert. the parameter file'
               ENDIF
            ELSEIF ( INDEX(TMPLABEL,'DOOM').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(TMPVALUE,I,TMP,LINE_LENGTH) ) THEN
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'DOOM from the parameter file'
               ELSE
                  IF ( TMPVALUE.EQ.1. ) THEN
                     OPER(DO_OM) = .TRUE.
                  ELSE
                     OPER(DO_OM) = .FALSE.
                  ENDIF
               ENDIF   
            ELSEIF ( INDEX(TMPLABEL,'OMPRI').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(OMPRI,I,TMP,LINE_LENGTH) )THEN
                  CMBPRI = 0.30
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'error in the OM prior from the parameter file'
               ENDIF
            ELSEIF ( INDEX(TMPLABEL,'DOMPR').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(DOMPR,I,TMP,LINE_LENGTH) )THEN
                  DCMBPR = .03
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'error in OM prior uncert. the parameter file'
               ENDIF

            ELSEIF ( INDEX(TMPLABEL,'DOHP').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(TMPVALUE,I,TMP,LINE_LENGTH) ) THEN
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'DOHP from the parameter file'
               ELSE
                  IF ( TMPVALUE.EQ.1. ) THEN
                     OPER(DO_HP) = .TRUE.
                  ELSE
                     OPER(DO_HP) = .FALSE.
                  ENDIF
               ENDIF
            ELSEIF ( INDEX(TMPLABEL,'PDF').GT.0 ) THEN
               I = N + 1
               CALL FINDNEXTWORD(TMP,I,N,LINE_LENGTH)
               IF ( INDEX(TMP(I:N),'GAUSSLIN').GT.0 ) THEN
                  WHICH_PDF = P_GAUSSLIN
                  I = N + 1
                  IF ( .NOT.GETONEVAL(COSFP,I,TMP,LINE_LENGTH) ) THEN
                     COSFP = 0.
                     WRITE (STDERR,'(2A)') 'readinit: error in ',
     $                    'reading GAUSSLIN from the parameter file'
                  ENDIF
               ELSEIF ( INDEX(TMP(I:N),'GAUSS').GT.0 ) THEN
                  WHICH_PDF = P_GAUSS
               ELSE
                  WHICH_PDF = P_GAUSS
               ENDIF
            ELSEIF ( INDEX(TMPLABEL,'WCUBEFILE').GT.0 .AND.
     $               OPER(MATE_EQST) ) THEN
               I = N + 1
               CALL FINDNEXTWORD(TMP,I,N,LINE_LENGTH)
               IF ( I.LE.N ) READ (TMP(I:N),'(A)', END = 27)
     $              WCUBEFILE
 27            CONTINUE
            ELSEIF ( INDEX(TMPLABEL,'DAVIDON').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(TMPVALUE,I,TMP,LINE_LENGTH) ) THEN
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'DAVIDON from the parameter file'
               ELSE
                  IF ( TMPVALUE.EQ.1. ) THEN
                     OPER(DO_DAVIDON) = .TRUE.
                  ELSE
                     OPER(DO_DAVIDON) = .FALSE.
                  ENDIF
               ENDIF
            ELSEIF ( INDEX(TMPLABEL,'POWELL').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(TMPVALUE,I,TMP,LINE_LENGTH) ) THEN
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'POWELL from the parameter file'
               ELSE
                  IF ( TMPVALUE.EQ.1. ) THEN
                     OPER(DO_POWELL) = .TRUE.
                  ELSE
                     OPER(DO_POWELL) = .FALSE.
                  ENDIF
               ENDIF
            ELSEIF ( INDEX(TMPLABEL,'GENETIC').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(TMPVALUE,I,TMP,LINE_LENGTH) ) THEN
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'GENETIC from the parameter file'
               ELSE
                  IF ( TMPVALUE.EQ.1. ) THEN
                     OPER(DO_GENETIC) = .TRUE.
                  ELSE
                     OPER(DO_GENETIC) = .FALSE.
                  ENDIF
               ENDIF
            ELSEIF ( INDEX(TMPLABEL,'GRIDSEARCH').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(TMPVALUE,I,TMP,LINE_LENGTH) ) THEN
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'GRIDSEARCH from the parameter file'
               ELSE
                  IF ( TMPVALUE.EQ.0. ) THEN
                     OPER(DO_GRID) = .FALSE.
                  ELSE
                     OPER(DO_GRID) = .TRUE.
                  ENDIF
               ENDIF
            ELSEIF ( INDEX(TMPLABEL,'SUBINTSIG').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(TMPVALUE,I,TMP,LINE_LENGTH) ) THEN
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'SUBINTSIG from the parameter file'
               ELSE
                  IF ( TMPVALUE.EQ.1. ) THEN
                     OPER(SUB_INTSIG) = .TRUE.
                  ELSE
                     OPER(SUB_INTSIG) = .FALSE.
                  ENDIF
               ENDIF

            ELSEIF ( INDEX(TMPLABEL,'ZERROR').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(TMPVALUE,I,TMP,LINE_LENGTH) ) THEN
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'ZERROR from the parameter file'
               ELSE
                  IF ( TMPVALUE.EQ.1. ) THEN
                     OPER(Z_ERROR) = .TRUE.
                     print *,'z_error = true'
                  ELSE
                     OPER(Z_ERROR) = .FALSE.
                     print *,'z_error = false'
                  ENDIF
               ENDIF

            ELSEIF ( INDEX(TMPLABEL,'LENSEFFECT').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(TMPVALUE,I,TMP,LINE_LENGTH) ) THEN
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'LENSEFFECT from the parameter file'
               ELSE
                  IF ( TMPVALUE.EQ.1. ) THEN
                     OPER(ADD_LDM) = .TRUE.
                  ELSE
                     OPER(ADD_LDM) = .FALSE.
                  ENDIF
               ENDIF
            ELSEIF ( INDEX(TMPLABEL,'GREYDUST').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(TMPVALUE,I,TMP,LINE_LENGTH) ) THEN
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'GREYDUST from the parameter file'
               ELSE
                  IF ( TMPVALUE.EQ.1. ) THEN
                     OPER(ADD_GDM) = .TRUE.
                  ELSE
                     OPER(ADD_GDM) = .FALSE.
                  ENDIF
               ENDIF
            ELSEIF ( INDEX(TMPLABEL,'HOSTEXTINCTION').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(TMPVALUE,I,TMP,LINE_LENGTH) ) THEN
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'HOSTEXTINCTION from the parameter file'
               ELSE
                  IF ( TMPVALUE.EQ.1. ) THEN
                     OPER(ADD_HOSTDM) = .TRUE.
                  ELSE
                     OPER(ADD_HOSTDM) = .FALSE.
                  ENDIF
               ENDIF
            ELSEIF ( INDEX(TMPLABEL,'PAMIXING').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(TMPVALUE,I,TMP,LINE_LENGTH) ) THEN
                  WRITE (STDERR,'(2A)') 'readinit: error in reading ',
     $                 'PAMIXING from the parameter file'
               ELSE
                  IF ( TMPVALUE.EQ.1. ) THEN
                     OPER(ADD_MIXDM) = .TRUE.
                  ELSE
                     OPER(ADD_MIXDM) = .FALSE.
                  ENDIF
               ENDIF
            ELSEIF ( INDEX(TMPLABEL,'NRSNE').GT.0 ) THEN
               I = N + 1
               IF ( .NOT.GETONEVAL(TMPVALUE,I,TMP,LINE_LENGTH) )THEN
                  WRITE (STDERR,'(2A)') 'readinit: error in reading',
     $                 ' ibandmaxfrom the parameter file'
               ELSE
                  NR_SUPER = INT(TMPVALUE)
               ENDIF
            ENDIF
         ENDIF
      ENDDO
 100  CONTINUE

C     
C          Take care of that variables that have not been initialized.
C  
      IF ( O_M.EQ.-1 ) THEN
         O_M = 0
         GLOB_OMM = 0.3
         WRITE (STDERR,'(2A,F4.2)') 'readinit: OMEGAM unspecified, ',
     $        'using: ', GLOB_OMM
      ENDIF
      IF ( O_X.EQ.-1 ) THEN
         O_X = 0
         GLOB_OMX = 0.7
         WRITE (STDERR,'(2A,F4.2)') 'readinit: OMEGAX unspecified, ',
     $        'using: ', GLOB_OMX
      ENDIF
      IF ( A_X.EQ.-1 ) THEN
         A_X = 0
         GLOB_AX = -1.0
         WRITE (STDERR,'(2A,F4.2)') 'readinit: W0 unspecified, ',
     $        'using: ', GLOB_AX
      ENDIF
      IF ( A_2.EQ.-1 ) THEN
         A_2 = 0
         GLOB_A2 = 0.0
         WRITE (STDERR,'(2A,F4.2)') 'readinit: W1 unspecified, ',
     $        'using: ', GLOB_A2
      ENDIF

C
C       Parameters for fitting Alam et al. metamorphosis
C       parametrization model.
C
      IF ( B_1.EQ.-1 .AND. OPER(DO_JACKE) ) THEN
         B_1 = 0
         GLOB_B1 = 0.0
         WRITE (STDERR,'(2A,F4.2)') 'readinit: B1 unspecified, ',
     $        'using: ', GLOB_B1
      ENDIF
      IF ( B_2.EQ.-1 .AND. OPER(DO_JACKE) ) THEN
         B_2 = 0
         GLOB_B2 = 0.0
         WRITE (STDERR,'(2A,F4.2)') 'readinit: B2 unspecified, ',
     $        'using: ', GLOB_B2
      ENDIF
      IF ( B_inv.EQ.-1 .AND. OPER(DO_JACKE) ) THEN
         B_inv = 0
         GLOB_Binv = 0.0
         WRITE (STDERR,'(2A,F4.2)') 'readinit: B_inv unspecified, ',
     $        'using: ', GLOB_Binv
      ENDIF
      IF ( B_4.EQ.-1 .AND. OPER(DO_JACKE) ) THEN
         B_4 = 0
         GLOB_B4 = 0.0
         WRITE (STDERR,'(2A,F4.2)') 'readinit: B4 unspecified, ',
     $        'using: ', GLOB_B4
      ENDIF

C
C     Parameters for fitting Linder parametrization
C
      IF ( L_1.EQ.-1 .AND. OPER(DO_LINDER) ) THEN
         L_1 = 0
         GLOB_L1 = 0.0
         WRITE (STDERR,'(2A,F4.2)') 'readinit: L1 unspecified, ',
     $        'using: ', GLOB_L1
      ENDIF
      IF ( L_2.EQ.-1 .AND. OPER(DO_LINDER) ) THEN
         L_2 = 0
         GLOB_L2 = 0.0
         WRITE (STDERR,'(2A,F4.2)') 'readinit: L2 unspecified, ',
     $        'using: ', GLOB_L2
      ENDIF


C
C     Parameters for fitting Wiltshire's cosmology
C
      IF ( FV_0.EQ.-1 .AND. OPER(DO_WILT) ) THEN
         FV_0 = 0
         GLOB_FV0 = 0.0
         WRITE (STDERR,'(2A,F4.2)') 'readinit: FV0 unspecified, ',
     $        'using: ', GLOB_FV0
      ENDIF
      IF ( H_0BAR.EQ.-1 .AND. OPER(DO_WILT) ) THEN
         H_0BAR = 0
         GLOB_H0BAR = 0.0
         WRITE (STDERR,'(2A,F4.2)') 'readinit: H0BAR unspecified, ',
     $        'using: ', GLOB_H0BAr
      ENDIF


C
C       Parameters for fitting Hannestad parametrization.
C
      IF ( W_0.EQ.-1 .AND. OPER(DO_HP) ) THEN
         W_0 = -1
         GLOB_W0 = -1.0
         WRITE (STDERR,'(2A,F4.2)') 'readinit: WLT unspecified, ',
     $        'using: ', GLOB_W0
      ENDIF
       IF ( W_1.EQ.-1 .AND. OPER(DO_HP) ) THEN
         W_1 = -1
         GLOB_W1 = -1.0
         WRITE (STDERR,'(2A,F4.2)') 'readinit: WET unspecified, ',
     $        'using: ', GLOB_W1
      ENDIF
      IF ( Q_HP.EQ.-1 .AND. OPER(DO_HP) ) THEN
         Q_HP = 0
         GLOB_QHP = 0.
         WRITE (STDERR,'(2A,F4.2)') 'readinit: QHP unspecified, ',
     $        'using: ', GLOB_QHP
      ENDIF
      IF ( Z_X.EQ.-1 .AND. OPER(DO_HP) ) THEN
         Z_X = 1
         GLOB_ZX = 1.0
         WRITE (STDERR,'(2A,F4.2)') 'readinit: ZX unspecified, ',
     $        'using: ', GLOB_ZX
      ENDIF

      IF ( M_SC.EQ.-1 ) THEN
         M_SC = 0
         GLOB_MSCRIPT = MAGIA - 5.0*ALOG10 (85.) + 25.0
         WRITE (STDERR,'(2A,F4.2)') 'readinit: SCRIPTM unspecified, ',
     $        'using: ', GLOB_MSCRIPT
      ELSEIF ( M_SC.GT.0 ) THEN
C
C           ScriptM is treated differently in the grid search and
C           must be put in the end of the parameter list.
C


         IF ( NR_VAR .EQ. O_M ) THEN
            CALL SWITCHPARS( M_SC, O_M )
         ELSEIF ( NR_VAR .EQ. O_X ) THEN
            CALL SWITCHPARS( M_SC, O_X )
         ELSEIF ( NR_VAR .EQ. A_X ) THEN
            CALL SWITCHPARS( M_SC, A_X )
         ELSEIF ( NR_VAR .EQ. A_2 ) THEN
            CALL SWITCHPARS( M_SC, A_2 )
         ELSEIF ( NR_VAR .EQ. ALP ) THEN
            CALL SWITCHPARS( M_SC, ALP )
         ELSEIF ( NR_VAR .EQ. SALP ) THEN
            CALL SWITCHPARS( M_SC, SALP )
         ELSEIF ( NR_VAR .EQ. H_0 ) THEN
            CALL SWITCHPARS( M_SC, H_0 )
         ELSEIF ( NR_VAR .EQ. O_RC ) THEN
            CALL SWITCHPARS( M_SC, O_RC )
         ELSEIF ( NR_VAR .EQ. X_DA ) THEN
            CALL SWITCHPARS( M_SC, X_DA )
         ELSEIF ( NR_VAR .EQ. B_1 ) THEN
            CALL SWITCHPARS( M_SC, B_1 )
         ELSEIF ( NR_VAR .EQ. B_2 ) THEN
            CALL SWITCHPARS( M_SC, B_2 )
          ELSEIF ( NR_VAR .EQ. B_inv ) THEN
            CALL SWITCHPARS( M_SC, B_inv )
         ELSEIF ( NR_VAR .EQ. B_4 ) THEN
            CALL SWITCHPARS( M_SC, B_4 )
         ELSEIF ( NR_VAR .EQ. L_1 ) THEN
            CALL SWITCHPARS( M_SC, L_1 )
         ELSEIF ( NR_VAR .EQ. L_2 ) THEN
            CALL SWITCHPARS( M_SC, L_2 )
         ELSEIF ( NR_VAR .EQ. FV_0 ) THEN
            CALL SWITCHPARS( M_SC, FV_0 )
         ELSEIF ( NR_VAR .EQ. H_0BAR ) THEN
            CALL SWITCHPARS( M_SC, H_0BAR )
         ELSEIF ( NR_VAR .EQ. W_0 ) THEN
            CALL SWITCHPARS( M_SC, W_0 )
         ELSEIF ( NR_VAR .EQ. W_1 ) THEN
            CALL SWITCHPARS( M_SC, W_1 )
         ELSEIF ( NR_VAR .EQ. Q_HP ) THEN
            CALL SWITCHPARS( M_SC, Q_HP )
         ELSEIF ( NR_VAR .EQ. Z_X ) THEN
            CALL SWITCHPARS( M_SC, Z_X )
         ELSEIF ( NR_VAR .EQ. M_SC ) THEN
         ELSE
            WRITE (STDERR,'(2A)') 'readinit: failed to switch the ',
     $           'nuisance parameters'
         ENDIF


         NR_VAR = NR_VAR - 1
         OPER(FITMSC) = .TRUE.
      ENDIF
      IF ( SALP.EQ.-1 ) THEN
         SALP = 0
         GLOB_SALP = 0.0
         WRITE (STDERR,'(2A,F4.2)') 'readinit: ALPHA unspecified ',
     $        'using: ', GLOB_SALP
      ELSEIF ( SALP.GT.0 ) THEN
C
C           Stretch alpha is treated differently in the grid search
C           and must be put in the end of the parameter list.
C
         IF ( NR_VAR .EQ. O_M ) THEN
            CALL SWITCHPARS( SALP, O_M )
         ELSEIF ( NR_VAR .EQ. O_X ) THEN
            CALL SWITCHPARS( SALP, O_X )
         ELSEIF ( NR_VAR .EQ. A_X ) THEN
            CALL SWITCHPARS( SALP, A_X )
         ELSEIF ( NR_VAR .EQ. A_2 ) THEN
            CALL SWITCHPARS( SALP, A_2 )
         ELSEIF ( NR_VAR .EQ. ALP ) THEN
            CALL SWITCHPARS( SALP, ALP )
         ELSEIF ( NR_VAR .EQ. M_SC ) THEN
            CALL SWITCHPARS( SALP, M_SC )
         ELSEIF ( NR_VAR .EQ. H_0 ) THEN
            CALL SWITCHPARS( SALP, H_0 )
         ELSEIF ( NR_VAR .EQ. O_RC ) THEN
            CALL SWITCHPARS( SALP, O_RC )
         ELSEIF ( NR_VAR .EQ. X_DA ) THEN
            CALL SWITCHPARS( SALP, X_DA )
         ELSEIF ( NR_VAR .EQ. B_1 ) THEN
            CALL SWITCHPARS( SALP, B_1 )
         ELSEIF ( NR_VAR .EQ. B_2 ) THEN
            CALL SWITCHPARS( SALP, B_2 )
         ELSEIF ( NR_VAR .EQ. B_inv ) THEN
            CALL SWITCHPARS( SALP, B_inv )
         ELSEIF ( NR_VAR .EQ. B_4 ) THEN
            CALL SWITCHPARS( SALP, B_4 )
         ELSEIF ( NR_VAR .EQ. L_1 ) THEN
            CALL SWITCHPARS( SALP, L_1 )
         ELSEIF ( NR_VAR .EQ. L_2 ) THEN
            CALL SWITCHPARS( SALP, L_2 )
         ELSEIF ( NR_VAR .EQ. FV_0 ) THEN
            CALL SWITCHPARS( SALP, FV_0 )
         ELSEIF ( NR_VAR .EQ. H_0BAR ) THEN
            CALL SWITCHPARS( SALP, H_0BAR )
         ELSEIF ( NR_VAR .EQ. W_0 ) THEN
            CALL SWITCHPARS( SALP, W_0 )
         ELSEIF ( NR_VAR .EQ. W_1 ) THEN
            CALL SWITCHPARS( SALP, W_1 )
         ELSEIF ( NR_VAR .EQ. Q_HP ) THEN
            CALL SWITCHPARS( SALP, Q_HP )
         ELSEIF ( NR_VAR .EQ. Z_X ) THEN
            CALL SWITCHPARS( SALP, Z_X )
         ELSEIF ( NR_VAR .EQ. SALP ) THEN
         ELSE
            WRITE (STDERR,'(2A)') 'readinit: failed to switch the ',
     $           'nuisance parameters'
         ENDIF
         NR_VAR = NR_VAR - 1
         OPER(FITSALP) = .TRUE.
      ENDIF
      IF ( H_0.EQ.-1 .AND. OPER(DO_MULTI) ) THEN
         H_0 = 0
         GLOB_H0 = .65
         WRITE (STDERR,'(2A,F4.2)') 'readinit: HUBBLE0 unspecified ',
     $        'using: ', GLOB_H0
      ENDIF

      RETURN
      END
C
C   End of READINIT
C



C     ****f* snalys/FINDNEXTWORD *
C
C   NAME
C     FINDNEXTWORD -- find the next word in a string that contains space
C
C   DESCRIPTION
C     Finds the next word in a string that contains some white space.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-10-25
C
C   USAGE
C     CALL FINDNEXTWORD( STR , B , E , LENGTH )
C
C   INPUTS
C     STR    - the CHARACTER string to be searched
C     B      - INTEGER that is the position where to start
C              searching for a word in the string.
C     E      - INTEGER
C     LENGTH - the INTEGER length of the whole string.
C
C
C   SIDE EFFECTS
C     B - contains the position of the first character in the word
C     E - contains the position of the last character in the word
C
C   HISTORY
C     2004-02-26 - Changed the conditions from .LE. to .LT. otherwise
C                  the LENGTH+1 element will be accessed in the next
C                  condition check.
C
C     ***
      SUBROUTINE FINDNEXTWORD(STR,B,E,LENGTH)
      IMPLICITNONE
      CHARACTER*(*) STR
      INTEGER B, E, LENGTH

      DO WHILE ( INDEX(STR(B:B),' ').EQ.1 .AND. B.LT.LENGTH )
         B = B + 1
      ENDDO
      E = B
      DO WHILE ( INDEX(STR(E:E),' ').NE.1 .AND. E.LT.LENGTH )
         E = E + 1
      ENDDO
      E = E - 1
      END
C
C   End of FINDNEXTWORD
C


C     ****i* snalys/GETONEVAL *
C
C   NAME
C     GETONEVAL -- read a value after a key
C
C   DESCRIPTION
C     Reads a float value after a string key in a file.
C
C   CREATION DATE
C     2003-07-08
C
C   USAGE
C     GETONEVAL(VAL,N,LINE,LINE_LENGTH)
C
C   INPUTS
C     VAL - Real
C     I - The position in the string where to start searching
C         for the value.
C     LINE - The string
C     LINE_LENGTH - The length of the string

      LOGICAL FUNCTION GETONEVAL(VAL,I,LINE,LINE_LENGTH)
      IMPLICIT NONE
      REAL VAL
      INTEGER I, LINE_LENGTH
      CHARACTER*(*) LINE

      INTEGER N

      GETONEVAL = .TRUE.
      CALL FINDNEXTWORD(LINE,I,N,LINE_LENGTH)

      VAL = -30000.
      IF ( I.LE.N ) READ (LINE(I:N),'(F13.6)', END = 1027) VAL      
 1027 CONTINUE
      IF ( VAL.EQ.-30000. ) GETONEVAL = .FALSE.
      END
C
C     End of GETONEVAL()
C


C     ****i* snalys/RVARVAL *
C
C   NAME
C     RVARVAL -- determines whether a parameter should be fitted
C
C   DESCRIPTION
C     Determines whether a parameter should be fitted or not, and also
C     stores the interval, and the number of grid points for that
C     parameter from the input string which is directly taken from
C     the ini-file format.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-10-25
C
C   USAGE
C     RVARVAL( STR , VNR , GLOB , BSTR , ESTR , LENGTH )
C
C   INPUTS
C     STR    - the CHARACTER string from the ini file
C     VNR    - the INTEGER number of this parameter
C     GLOB   - the global DOUBLE PRECISION variable for this
C              parameter from cosmo.inc.
C     BSTR   - where to start searching in the string.
C     ESTR   - where to stop searching in the string.
C     LENGTH - of the string.
C
C   USES
C     cos_int
C
C   RESULT
C     A LOGICAL value is returned for confirmation.
C
C   SIDE EFFECTS
C     GLOB will be set to the given value if this parameters is not
C     going to be fitted.
C     The COS_INT array is going to be updated with the input data.
C
C   WARNINGS
C     A warning message will be written if the input string is not
C     of the correct type.
C
C     ***
      LOGICAL FUNCTION RVARVAL(STR,VNR,GLOB,BSTR,ESTR,LENGTH)
      IMPLICIT NONE
      CHARACTER*100 STR
      INTEGER VNR, BSTR, ESTR, LENGTH
      DOUBLE PRECISION GLOB

      INCLUDE 'cosmo.inc'
      INCLUDE 'error.inc'

      INTEGER B, E, N, NRNUM
      PARAMETER ( NRNUM = 3 ) 
      REAL TMP(NRNUM), NOVAL
      PARAMETER ( NOVAL = -10000.0 )

      RVARVAL = .TRUE.

C  
C        Try to find the three values defining this parameter.
C  

      E = ESTR
      DO N = 1, NRNUM
         TMP(N) = NOVAL
         
C  
C           Get rid of spaces.
C   
         B = E + 1
         CALL FINDNEXTWORD(STR,B,E,LENGTH)
         
C  
C           Is this a value? Then it should be stored.
C  
         IF ( B.LE.E ) THEN
            IF ( INDEX(STR(B:E),'.').GT.0 ) THEN
               READ ( STR(B:E),'(F13.6)', ERR=80 ) TMP(N)
            ELSE
               READ (STR(B:E),'(I6)', ERR=80) B
               TMP(N) = REAL(B)
            ENDIF
         ENDIF
 80      CONTINUE
      ENDDO
      
C  
C        Determine whether this parameter should be estimated or kept
C        constant.
C  
      IF ( ( TMP(3).EQ.NOVAL .OR. TMP(2).EQ.NOVAL .OR. TMP(3).EQ.1.0 
     $     .OR. TMP(3).EQ.0.0 .OR. TMP(1).EQ.TMP(2) ) .AND. 
     $     TMP(1).NE.NOVAL ) THEN
         GLOB = TMP(1)
         VNR = 0
      ELSEIF ( TMP(1).NE.NOVAL ) THEN
         NR_VAR = NR_VAR + 1
         VNR = NR_VAR
         P_NAMES(VNR) = STR(BSTR:ESTR)
         IF ( TMP(1).GT.TMP(2) ) THEN
            TMP(3) = TMP(1)
            TMP(1) = TMP(2)
            TMP(2) = TMP(3)
         ENDIF
         CALL C_START_PUT(NR_VAR,TMP(1))
         CALL C_END_PUT(NR_VAR,TMP(2))
         CALL C_NR_PUT(NR_VAR,TMP(3))
      ELSE
         WRITE (STDERR,'(3A)') 'rvarval: the parameter ',STR(BSTR:ESTR),
     $        ' was not defined correctly in the parameter file.'
         RVARVAL = .FALSE.
      ENDIF
      
      RETURN
      END
C
C End of RVARVAL
C



C     ****i* snalys/RMPAR *
C
C   NAME
C     RMPAR -- untoggle parameters that already have been toggled
C              for beeing fit
C
C   DESCRIPTION
C     When a parameter is defined to be fit in the ini file, it
C     is added to the list of such variables and the fitting intervals
C     are defined. If for some reason it later turns out that this
C     parameter can not be fit (e.g. flat universe assumption) this
C     function removes the parameter from such lists etc.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-09-07
C
C   USAGE
C     RMPAR( PAR )
C
C   INPUTS
C     PAR    - INTEGER index in the parameter list telling the
C              number of the parameter that should be removed.
C
C   RESULT
C     The INTEGER index of par is returned to 0 if the removal
C     was successful otherwise PAR is returned.
C
C   USES
C     cos_int
C
C   SIDE EFFECTS
C     The global variable for the cosmological parameter from cosmo.inc
C     will be set to the start value in the fitting array COS_INT.
C     This array will of course also be updated accordingly. The COMMON
C     variable NR_VAR will be shortened.
C
C   ERRORS
C     An error message is written if the parameter to be removed does
C     not exist.
C
C   HISTORY
C     2002-01-09 Fixed a bug that was introduced when I added
C                alpha and h0 as new parameters.
C     2003-03-22 Fixed a serious bug in this routine. The cosmo-struct
C                and P_NAMES structures were never updated when a parameter
C                was removed. Suprised that nobody noticed this before!
C
C     ***
      INTEGER FUNCTION RMPAR ( PAR )
      IMPLICIT NONE
      INTEGER PAR

      INCLUDE 'cosmo.inc'
      INCLUDE 'error.inc'

      INTEGER N , TINT
      REAL C_START_GET , C_END_GET , C_NR_GET
      
C
C        Set the global variable for this parameter to
C        C_START_GET.
C
      IF ( PAR .EQ. O_M ) THEN
         GLOB_OMM =  C_START_GET( PAR )
      ELSEIF ( PAR .EQ. O_X ) THEN
         GLOB_OMX =  C_START_GET( PAR )
      ELSEIF ( PAR .EQ. A_X ) THEN
         GLOB_AX =  C_START_GET( PAR )
      ELSEIF ( PAR .EQ. A_2 ) THEN
         GLOB_A2 =  C_START_GET( PAR )         
      ELSEIF ( PAR .EQ. M_SC ) THEN
         GLOB_MSCRIPT =  C_START_GET( PAR )
      ELSEIF ( PAR .EQ. ALP ) THEN
         GLOB_ALPHA = C_START_GET( PAR )
      ELSEIF ( PAR .EQ. SALP ) THEN
         GLOB_SALP = C_START_GET( PAR )
      ELSEIF ( PAR .EQ. H_0 ) THEN
         GLOB_H0 = C_START_GET( PAR )
      ELSEIF ( PAR .EQ. O_RC ) THEN
         GLOB_ORC = C_START_GET( PAR )
      ELSEIF ( PAR .EQ. X_DA ) THEN
         GLOB_XDA = C_START_GET( PAR )
      ELSEIF ( PAR .EQ. B_1 ) THEN
         GLOB_B1 = C_START_GET( PAR )
      ELSEIF ( PAR .EQ. B_2 ) THEN
         GLOB_B2 = C_START_GET( PAR )
      ELSEIF ( PAR .EQ. B_inv ) THEN
         GLOB_Binv = C_START_GET( PAR )
      ELSEIF ( PAR .EQ. B_4 ) THEN
         GLOB_B4 = C_START_GET( PAR )
      ELSEIF ( PAR .EQ. L_1 ) THEN
         GLOB_L1 = C_START_GET( PAR )
      ELSEIF ( PAR .EQ. L_2 ) THEN
         GLOB_L2 = C_START_GET( PAR )
      ELSEIF ( PAR .EQ. FV_0 ) THEN
         GLOB_FV0 = C_START_GET( PAR )
      ELSEIF ( PAR .EQ. H_0BAR ) THEN
         GLOB_H0BAR = C_START_GET( PAR )
      ELSEIF ( PAR .EQ. W_0 ) THEN
         GLOB_W0 = C_START_GET( PAR )
      ELSEIF ( PAR .EQ. W_1 ) THEN
         GLOB_W1 = C_START_GET( PAR )
       ELSEIF ( PAR .EQ. Q_HP ) THEN
         GLOB_QHP = C_START_GET( PAR )
      ELSEIF ( PAR .EQ. Z_X ) THEN
         GLOB_ZX = C_START_GET( PAR )
      ENDIF

      RMPAR = 0
      
      IF ( RMPAR .LT. NR_VAR ) THEN
         DO N = PAR + 1, NR_VAR
            TINT = -1
            IF ( N .EQ. O_M ) THEN
               TINT = O_M
               O_M = O_M - 1
            ELSEIF ( N .EQ. O_X ) THEN
               TINT = O_X
               O_X = O_X - 1
            ELSEIF ( N .EQ. A_X ) THEN
               TINT = A_X
               A_X = A_X - 1
            ELSEIF ( N .EQ. A_2 ) THEN
               TINT = A_2
               A_2 = A_2 - 1
            ELSEIF ( N .EQ. M_SC ) THEN
               TINT = M_SC
               M_SC = M_SC - 1
            ELSEIF ( N .EQ. ALP ) THEN
               TINT = ALP
               ALP = ALP - 1
            ELSEIF ( N .EQ. SALP ) THEN
               TINT = SALP
               ALP = SALP - 1
            ELSEIF ( N .EQ. H_0 ) THEN
               TINT = H_0
               H_0 = H_0 - 1
            ELSEIF ( N .EQ. O_RC ) THEN
               TINT = O_RC
               O_RC = O_RC - 1
            ELSEIF ( N .EQ. X_DA ) THEN
               TINT = X_DA
               X_DA = X_DA - 1
            ELSEIF ( N .EQ. B_1 ) THEN
               TINT = B_1
               B_1 = B_1 - 1
            ELSEIF ( N .EQ. B_2 ) THEN
               TINT = B_2
               B_2 = B_2 - 1
            ELSEIF ( N .EQ. B_inv ) THEN
               TINT = B_inv
               B_inv = B_inv - 1
            ELSEIF ( N .EQ. B_4 ) THEN
               TINT = B_4
               B_4 = B_4 - 1
            ELSEIF ( N .EQ. L_1 ) THEN
               TINT = L_1
               L_1 = L_1 - 1
            ELSEIF ( N .EQ. L_2 ) THEN
               TINT = L_2
               L_2 = L_2 - 1
            ELSEIF ( N .EQ. FV_0 ) THEN
               TINT = FV_0
               FV_0 = FV_0 - 1
            ELSEIF ( N .EQ. H_0BAR ) THEN
               TINT = H_0BAR
               H_0BAR = H_0BAR - 1
            ELSEIF ( N .EQ. W_0 ) THEN
               TINT = W_0
               W_0 = W_0 - 1
            ELSEIF ( N .EQ. W_1 ) THEN
               TINT = W_1
               W_1 = W_1 - 1
            ELSEIF ( N .EQ. Q_HP ) THEN
               TINT = Q_HP
               Q_HP = Q_HP - 1
            ELSEIF ( N .EQ. Z_X ) THEN
               TINT = Z_X
               Z_X = Z_X - 1
            ELSE
               RMPAR = PAR
               WRITE (STDOUT,'(2A)') 'rmpar: this could never ',
     $              'happen...'
            ENDIF
C
C       If the parameter was found we must also update the
C       array structure that keeps the fitting intervals, and
C       parameter names.
C
            IF ( RMPAR .EQ. 0 ) THEN
               CALL C_START_PUT( TINT - 1 , C_START_GET(TINT) )
               CALL C_END_PUT( TINT - 1 , C_END_GET(TINT) )
               CALL C_NR_PUT( TINT - 1 , C_NR_GET(TINT) )
               P_NAMES( TINT - 1 ) = P_NAMES( TINT )
            ENDIF

         ENDDO
      ENDIF
      NR_VAR = NR_VAR - 1
      RETURN
      END
C
C     END OF RMPAR
C



C     ****f* cos_int/SWITCHPARS *
C
C   NAME
C     SWITCHPARS -- switch the order of two parameters in the fit
C                   matrix
C
C   DESCRIPTION
C     When a parameter is defined to be fit in the ini file, it
C     is added to the list of such variables and the fitting intervals
C     are defined. If for some reason it later turns out that this
C     parameter can not be fit (e.g. flat universe assumption) this
C     function removes the parameter from such lists etc.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2003-07-11
C
C   USAGE
C     SWITCHPARS( PAR1, PAR2 )
C
C   INPUTS
C     PAR1   - INTEGER index in the parameter list telling the
C              number of the parameter that should be switched
C              with PAR2.
C     PAR2   - INTEGER index in the parameter list telling the
C              number of the parameter that should be switched
C              with PAR1.
C
C
C   RESULT
C     A BOOLEAN is returned indicating whether the switching was
C     successful or not.
C
C   USES
C     COS_INT
C
C   SIDE EFFECTS
C     * The COS_INT variable will be changed accordingly.
C     * PAR1 and PAR2 will switch values
C
C   ERRORS
C     An error message is written if the parameter to be removed does
C     not exist.
C
C     ***
      SUBROUTINE SWITCHPARS( PAR1 , PAR2 )
      IMPLICIT NONE
      INTEGER PAR1, PAR2

      INCLUDE 'cosmo.inc'
      INCLUDE 'error.inc'

      INTEGER TMPPAR
      REAL C_START_GET , C_END_GET , C_NR_GET,
     $     CSTART, CEND, C_STEP_GET, CSTEP, CNR
      CHARACTER*(10) TMPNAME
      
      IF ( PAR1 .LE. NR_VAR .AND. PAR1 .GT. 0 .AND.
     $     PAR2 .LE. NR_VAR .AND. PAR2 .GT. 0 ) THEN
         CSTART = C_START_GET( PAR1 )
         CNR = C_NR_GET( PAR1 )
         CEND = C_END_GET( PAR1 )
         CSTEP = C_STEP_GET( PAR1 )
         
C         WRITE (*,*) CSTART, CNR, CEND,  C_START_GET( PAR2 ),
C     $        C_NR_GET( PAR2 ), C_END_GET( PAR2 )

         CALL C_START_PUT( PAR1 , C_START_GET( PAR2 ) )
         CALL C_NR_PUT( PAR1 , C_NR_GET( PAR2 ) )
         CALL C_END_PUT( PAR1 , C_END_GET( PAR2 ) )
         CALL C_STEP_PUT( PAR1 , C_STEP_GET( PAR2 ) )
         CALL C_START_PUT( PAR2 , CSTART )
         CALL C_NR_PUT( PAR2 , CNR )
         CALL C_END_PUT( PAR2 , CEND )
         CALL C_STEP_PUT( PAR2, CSTEP )
         
         TMPNAME = P_NAMES(PAR1)
         P_NAMES(PAR1) = P_NAMES(PAR2)
         P_NAMES(PAR2) = TMPNAME

         TMPPAR = PAR1
         PAR1 = PAR2
         PAR2 = TMPPAR
      ENDIF
      
      RETURN
      END
C
C     END OF SWITCHPARS
C
