C     ****f* snalys/ANALYZE *
C
C   NAME
C     ANALYZE() -- fits the best cosmologogy to the SN-data
C
C   DESCRIPTION
C     This is the top function for the analysis, and from here the
C     method for the analysis is controlled.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-10-25
C
C   USAGE
C     ANALYZE( P_COS , OFID , FNAME )
C
C   INPUTS
C     P_COS   - A REAL array.
C     OFID    - the INTEGER ID value of the output file
C     FNAME   - CHARACTER string that contains the name of the output
C               file.
C
C   RESULT
C     A boolean is returned as a confirmation on the result of the fit.
C
C   SIDE EFFECTS
C     P_COS   - a array of type REAL, that will contain the best fit
C               on return
C
C   TODO
C     * Modifications for the 'contourfinder' package whenever it will
C       be implemented.
C
C   HISTORY
C     2001-09-08 Removed the posiblity to interactively draw the path
C                of the Davidon process. In the future it is possible
C                the plotting features will be completely removed. If
C                a debugging tool as this is required it should be
C                enough just to write the points to STDOUT.
C     2010-02-11 Changed to using CPU_TIME for retrieving elapsed time
C
C     ***
      LOGICAL FUNCTION ANALYZE( P_COS, OFID, FNAME )
      IMPLICIT NONE
      REAL P_COS(*)
      INTEGER OFID

      INCLUDE 'cosmo.inc'
      INCLUDE 'supernova.inc'
      INCLUDE 'operation.inc'
      INCLUDE 'error.inc'
      INCLUDE 'snalys.inc'
      INCLUDE 'probability.inc'

      CHARACTER*(*) FNAME

      INTEGER N, I, K, CLEN, IDX(NR_PAR), P_IDX(NR_PAR)
      REAL STIME
      REAL GETVAL
C
C        A temp string
C
      CHARACTER*20 C
C
C        The length of the temp read string.
C
      PARAMETER (CLEN = 20)

      REAL C_STEP_GET, C_START_GET, C_END_GET, IDX_TO_VAL
      REAL LOGMLFUNC, C_NR_GET, ZBINFUNC, CALCM, CALCDL
      EXTERNAL LOGMLFUNC, IDX_TO_VAL, ZBINFUNC, CALCM, CALCDL

      LOGICAL UPDATE_GLOBALS, SUBNETFUNC, SCONFIND, CONFIND, DAVIDO,
     $     PRINTOUT
      EXTERNAL SUBNETFUNC

C
C        f_den is the calculated probability density function (p.d.f) for the
C        nr_var variables, that is it is actually a om_len*ox_len*ax_len
C        array and the "array package" array.f should be used to handle it.
C
      REAL F_DEN(MAX_LEN)

C
C        Most of these variables are used in the powell minimization
C        process.
C        The most probable cosmology (p_cos). var_mat is the correlation
C        matrix for the unknown parameters. xi is the identity matrix.
C
      REAL XI(NR_PAR,NR_PAR),MAGDIST(MAX_SN)
      REAL GRI_COSMO(NR_PAR+1),GEN_COSMO(NR_PAR+1),POW_COSMO(NR_PAR+1),
     $     DAV_COSMO(NR_PAR+1),BIN_COSMO(NR_PAR+1),GEN_START(NR_PAR),
     $     POW_START(NR_PAR),DAV_START(NR_PAR)
      REAL VAR_MAT(NR_PAR**2), FRET, RAND
      INTEGER ITER,GADZ,LEN(NR_PAR),TOT_LEN, NETNR
      PARAMETER (GADZ = NR_PAR**2)
      DATA XI/GADZ*0.0/
      DATA VAR_MAT/GADZ*0.0/

C
C        The gradient vector, error vector, and correlation matrix.
C
      REAL GRAD(NR_PAR), DERR(NR_PAR), CORR(NR_PAR,NR_PAR), QINI,
     $     GERR(2,NR_PAR)

C
C        Retrieve the start time
C
      CALL CPU_TIME(STIME)

C
C        Make xi the identity matrix
C
      DO N = 1, NR_PAR
         XI(N,N) = 0.01
      ENDDO

      ANALYZE = .TRUE.   ! default return value from function ANALYZE


C
C        GRID SEARCH
C
      If ( OPER(DO_GRID) ) THEN
         IF (.NOT.OPER(QUIET) ) WRITE (STDOUT,'(A)')
     $        'Starting Grid Search...'
         TOT_LEN = 1
         DO N = 1, NR_VAR
C
C              The length of the vectors in the array f_den.
C
            LEN(N) = INT (C_NR_GET(N))
            TOT_LEN = TOT_LEN*LEN(N)
         ENDDO
         CALCDERIV = .FALSE.
         CALL GRID_SEARCH( LOGMLFUNC, NR_VAR, TOT_LEN, IDX_TO_VAL, 
     $        P_COS, F_DEN, OPER(QUIET) )
         DO N = 1,NR_PAR + 1
            GRI_COSMO(N) = P_COS(N)
         ENDDO

C
C           Find the standard deviation from the grid search
C
         I = 1
         DO N = 1,NR_VAR
            I = I*INT(C_NR_GET(N))
            LEN(N) = C_NR_GET(N)
            P_IDX(N) = INT((P_COS(N)-C_START_GET(N))/C_STEP_GET(N)) + 1
            GERR(1,N) = REAL(P_IDX(N))
            GERR(2,N) = REAL(P_IDX(N))
         ENDDO
         DO N = 2,I
            IF ( (F_DEN(N-1).GE.(GRI_COSMO(NR_VAR+1) + 
     $           0.5*SIG_Z(NRSIG,1)**2) .AND.
     $           F_DEN(N).LE.(GRI_COSMO(NR_VAR+1) + 
     $           0.5*SIG_Z(NRSIG,1)**2)) .OR.
     $           (F_DEN(N-1).LE.(GRI_COSMO(NR_VAR+1) +
     $           0.5*SIG_Z(NRSIG,1)**2) .AND.
     $            F_DEN(N).GE.(GRI_COSMO(NR_VAR+1) +
     $           0.5*SIG_Z(NRSIG,1)**2))                  ) THEN
               CALL GET_NORMAL_IDX(NR_VAR,N,LEN,IDX)
               DO K = 1,NR_VAR
                  IF ( IDX(K).LT.P_IDX(K).AND.
     $                 IDX(K).LE.INT(GERR(1,K))     ) THEN
                     GERR(1,K) = REAL(IDX(K))
                  ELSE IF ( IDX(K).GT.P_IDX(K).AND.
     $                      IDX(K).GE.INT(GERR(2,K))     ) THEN
                     GERR(2,K) = REAL(IDX(K))
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
         DO N = 1,NR_VAR
            GERR(1,N) = C_START_GET(N) + (GERR(1,N)-1.0)*C_STEP_GET(N)
            GERR(2,N) = C_START_GET(N) + (GERR(2,N)-1.0)*C_STEP_GET(N)
            GERR(1,N) = GRI_COSMO(N) - GERR(1,N)
            GERR(2,N) = GERR(2,N) - GRI_COSMO(N)
         ENDDO
      ENDIF   
     
C
C        POWELL'S METHOD
C
      IF ( OPER(DO_POWELL) ) THEN
         DO N = 1, NR_VAR
            POW_START(N) = P_COS(N)
         ENDDO

         IF ( .NOT. OPER(QUIET) ) THEN
            WRITE (STDOUT,'(A)') 'Starting powell-minimization...'
         ENDIF
         CALL POWELL( P_COS, XI, NR_VAR, NR_PAR, 0.001, ITER, FRET )
         DO N = 1,NR_PAR + 1
            POW_COSMO(N) = P_COS(N)
         ENDDO
         IF ( .NOT. OPER(QUIET) ) THEN
            WRITE (STDOUT,'(A)') 'done.'
         ENDIF
      endif
    
C
C        DAVIDON'S METHOD
C
      IF ( OPER(DO_DAVIDON) ) THEN
         DO N = 1, NR_VAR
            DAV_START(N) = P_COS(N)
         ENDDO
         
         IF (.NOT.OPER(QUIET) ) THEN
            WRITE (STDOUT,'(A)') 'Starting davidon-minimization...'
         ENDIF
         
         I = 1
         OPER(DO_DAVIDON) = .FALSE.
         DO WHILE ( .NOT.OPER(DO_DAVIDON) .AND. I.LE.NR_VAR**2 )
            I = I + 1

C
C              Start values that always will upset Davidon...
C

C           P_COS(O_M) = 2.62089
C           P_COS(O_X) = 1.34048
C           P_COS(M_SC) = -3.57636

            CALCDERIV = .TRUE.
            QINI = LOGMLFUNC( P_COS, GRAD )

C
C              Use this to set the step length
C
            DO N = 1,GADZ
               VAR_MAT(N) = 0.
            END DO

            DO N = 1,NR_VAR
C
C                 Set the step length to 1% of the point value if grad = 0.
C
               IF ( GRAD(N).EQ.0.0 ) THEN
                  VAR_MAT((N-1)*NR_VAR+N) = 0.01*P_COS(N)
C
C                 This may look a bit funny, but it is done to keep the
C                 steplength 0.01 after the davidon routine mess up v.
C
               ELSE
                  VAR_MAT((N-1)*NR_VAR+N) = ABS( 0.01 / GRAD(N) )
               ENDIF
            ENDDO

            OPER(DO_DAVIDON) = DAVIDO ( -NR_VAR, P_COS, VAR_MAT,
     $           LOGMLFUNC, 1.0E-7 )

C
C              If the minimization was unsuccessful, try different
C              start values.
C
            IF ( .NOT. OPER(DO_DAVIDON) ) CALL RANDCOSMO(P_COS)
         ENDDO

         IF ( .NOT. OPER(DO_DAVIDON) ) THEN
            WRITE (STDERR,'(2A)')'warning: maximum nr of iterations in',
     $           ' the davidon minimization'
            WRITE (STDERR,'(A)') '         routine has been reached'
C
C              If the davidon minimizations was unsuccessful do not
C              do anything more related to this way of minimizing
C              (e.g. printing and drawing figures).
C
            OPER(DO_DAVIDON) = .FALSE.
         ENDIF

         DO N = 1, NR_PAR + 1
            DAV_COSMO(N) = P_COS(N)
         ENDDO

C
C           The correlation matrix
C
         DO N = 1,NR_VAR
            DERR(N) = SQRT(MAX(0.0,VAR_MAT(NR_VAR*(N-1)+N)))
            DO I = 1,N
               IF ( DERR(N).GT.0.0 .AND. DERR(I).GT.0.0 ) THEN
                  CORR(N,I) = VAR_MAT(NR_VAR*(N-1)+I)/
     $                 (DERR(N)*DERR(I))
                  CORR(I,N) = CORR(N,I)
               ENDIF
            ENDDO
         ENDDO

         IF ( .NOT.OPER(QUIET) ) WRITE (STDOUT,'(A)') 'done.'
      ENDIF

C      IF ( OPER(CON_FIND) ) THEN
C         NETNR = 1000
C 
C
C           THE CONTOURFIND PACKAGE ( TO BE IMPLIMENTED )
C        
C         IF ( SCONFIND( LOG_ML_GRAD, P_COS, 0.5*SIG_Z(NRSIG,1)**2,
C     $        F_DEN, NETNR, NR_VAR, 0.0001, 1, SUBNETFUNC ) )
C     $        OPER(CON_FIND) = .TRUE.
C         CALCDERIV = .TRUE.
C         IF ( CONFIND( LOG_ML_FUNC, P_COS, 0.5*SIG_Z(NRSIG,1)**2,
C     $        F_DEN, NETNR, NR_VAR, 0.0001, 2) )
C     $        OPER(CON_FIND) = .TRUE.
C         DO N = 1, NETNR/NR_VAR
C            WRITE (*,*) (F_DEN((N-1)*NR_VAR+I), I = 1, NR_VAR)
C         ENDDO
C      ENDIF

C
C        Use the estimation to produce a distribution of the magnitudes.
C
      IF ( OPER(MAG_DIST) ) THEN
c       GLOB_OMM = 0.3
c       GLOB_OMX = 0.7
c       GLOB_AX = -1.0
c       GLOB_A2 = 0.0
c       GLOB_MSCRIPT = -3.39709472656250
C       WRITE (*,*) GLOB_OMM, GLOB_OMX, GLOB_AX
        DO N=1,SN_NR
           MAGDIST(N) = CALCM( CALCDL(0.,SN_DATA(ZS,N)) ,P_COS(M_SC))
     $          - (SN_DATA(BZMAG,N) - SN_DATA(ISG,N))
        ENDDO
      ENDIF

      
C
C        GENETIC'S METHOD
C
      If ( OPER(DO_GENETIC) ) THEN
         IF (.NOT.OPER(QUIET) ) WRITE (STDOUT,'(A)')
     $        'Starting genetic-minimization...'
C        CALCDERIV = .FALSE.
         CALL GENETIC( NR_VAR, P_COS )
         DO N = 1,NR_PAR + 1
            GEN_COSMO(N) = P_COS(N)
         ENDDO
         IF ( .NOT. OPER(QUIET) ) THEN
            WRITE (STDOUT,'(A)') 'done.'
         ENDIF
      ENDIF      
      
      
C
C        If we only vary one or two of the cosmological parameters it
C        is possible to make a graphical illustration of the p.d.f.
C        In these cases, f_den is printed.
C
      PRINTOUT = .FALSE.
      
      DO N = 1, NRMIN
         IF ( OPER(N) ) PRINTOUT = .TRUE.
      ENDDO

      IF ( PRINTOUT ) THEN
         CALL RESULTS( STIME , F_DEN, MAGDIST, BIN_COSMO,
     $        GRI_COSMO, GEN_COSMO, GEN_START, POW_COSMO,
     $        POW_START, DAV_COSMO, DAV_START, GERR,
     $        DERR, CORR )
      ELSE
         ANALYZE = .FALSE.
      ENDIF

      RETURN
      END
C
C     End of ANALYZE
C



C     ****f* snalys/GETVAL *
C
C   NAME
C     GETVAL() -- retrieves a value from a string
C
C   DESCRIPTION
C     If a string consists of a word and a number separated by
C     a arbitrary number of blank spaces, the number is returned.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-10-25
C
C   USAGE
C     GETVAL( C )
C
C   INPUTS
C     C - string of type CHARACTER*20
C
C   RESULT
C     The value is returned of type REAL
C
C     ***
      REAL FUNCTION GETVAL( C )
      CHARACTER*20 C

      INTEGER i
      I = 1
      DO WHILE ( C(I:I).NE.' ' )
         I = I + 1
      ENDDO
      READ (C(I:LEN(C)),'(F20.0)') GETVAL
      RETURN
      END
C
c End of GETVAL()
C


C     ****f* snalys/SUBNETFUNC *
C
C   NAME
C     SUBNETFUNC() -- function to be applied on each subnet
C
C   DESCRIPTION
C     This function defines what is to be done with the points returned
C     for each subnet produced by the CONTOURFINDER package, and is sent
C     as an argument to CONFIND or SCONFIND.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-10-25
C
C   USAGE
C     SUBNETFUNC( P, NET, NRDIM )
C
C   INPUTS
C     P     - an array of type INTEGER, describing the indices in
C             the net that constitute the subnet.
C     NET   - an array of type REAL defining the whole net.
C     NRDIM - INTEGER that gives the dimension of the parameter space.
C
C   RESULT
C     A boolean value is returned, currently always .TRUE.
C
C   SEE ALSO
C     CONTOURFINDER, SCONFIND()
C
C   TODO
C     * Decide what to do with the subnets, once CONTOURFINDER is
C       implemented.
C
C     ***
      LOGICAL FUNCTION SUBNETFUNC( P, NET, NRDIM )
      IMPLICIT NONE
      INTEGER P(*), NRDIM
      REAL NET(*)

      SUBNETFUNC = .TRUE.

C      WRITE (*,*) 'triangle { <',NET(P(1)),',',NET(P(1)+1),
C     $     ',',NET(P(1)+2),'>, <',NET(P(2)),',',NET(P(2)+1),
C     $     ',',NET(P(2)+2),'>, <',NET(P(3)),',',NET(P(3)+1),
C     $     ',',NET(P(3)+2),'> }'

      WRITE (*,*) NET(P(1)),NET(P(1)+1)
      WRITE (*,*) NET(P(2)),NET(P(2)+1)
C     WRITE (*,*) NET(P(3)),NET(P(3)+1)
      RETURN
      END
C
C     END OF SUBNETFUNC
C

C
C End of analyze.f
C
