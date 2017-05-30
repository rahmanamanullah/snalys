C     ****f* snalys/SETUP_PLOT *
C
C   NAME
C     SETUP_PLOT -- open the graphport
C
C   DESCRIPTION
C     Opens the graph port to a specified device, and sets the title
C     for the plot.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-10-25
C
C   USAGE
C     SETUP_PLOT(DEV,PLOT,TITLE)
C
C   INPUTS
C     DEV   - CHARACTER string, the device for the plot.
C     PLOT  - 2D INTEGER vector, specifying what variables
C             to plot.
C     TITLE - the title of the plot.
C
C   RESULT
C     A LOGICAL variable is returned confiriming that the setup
C     was successful.
C
C     ***
      LOGICAL FUNCTION SETUP_PLOT(DEV,PLOT,TITLE)
      IMPLICIT NONE
      CHARACTER*(260) DEV
      CHARACTER*(256) TITLE
      INTEGER PLOT(*)

      INCLUDE 'cosmo.inc'

      CHARACTER*30 LABELS(2)
      INTEGER IER, PGOPEN, N, POLYDIM
      PARAMETER ( POLYDIM = 200 )
      REAL C_START_GET, C_END_GET, BIGBANG, XPTS(POLYDIM), YPTS(POLYDIM)
      EXTERNAL BIGBANG

      SETUP_PLOT = .TRUE.
C     
C           Open the graph port
C        
      IER = PGOPEN(DEV)
      
      IF ( IER.EQ.1 ) THEN
         CALL PGSCR(0,1.,1.,1.) ! Set the background colour.
         CALL PGSCR(1,0.,0.,0.) ! Set the the pen colour.
C     
C     Define the labels
C     
         IF ( PLOT(1).EQ.O_M) THEN
            LABELS(1) = '\\gW\\dM\\d'
         ELSE IF ( PLOT(1).EQ.O_X ) THEN
            IF ( GLOB_AX.EQ.-1.0 ) THEN
               LABELS(1) = '\\gW\\d\\gL\\d'
            ELSE
               LABELS(1) = '\\gW\\dX\\d'
            ENDIF
         ELSE IF ( PLOT(1).EQ.A_X ) THEN
            LABELS(1) = 'w\\d0\\d'
         ELSE IF ( PLOT(1).EQ.A_2 ) THEN
            LABELS(1) = 'w\\d1\\d'
         ENDIF
         
         IF ( PLOT(2).EQ.O_X ) THEN
            IF ( GLOB_AX.EQ.-1.0 ) THEN
               LABELS(2) = '\\gW\\d\\gL \\d'
            ELSE
               LABELS(2) = '\\gW\\dX\\d'
            ENDIF
         ELSE IF ( PLOT(2).EQ.A_X ) THEN
            LABELS(2) = 'w\\d0\\d'
         ELSE IF ( PLOT(2).EQ.A_2 ) THEN
            LABELS(2) = 'w\\d1\\d'
         ENDIF
         
         CALL PGENV(C_START_GET(PLOT(1)),C_END_GET(PLOT(1)),
     +        C_START_GET(PLOT(2)),C_END_GET(PLOT(2)),0,1)
         CALL PGLAB(LABELS(1),LABELS(2),TITLE)
      ELSE
         SETUP_PLOT = .FALSE.
      ENDIF
      
C     
C     Draw the no bigbang area.
C  
      IF ( (PLOT(1).EQ.O_M) .AND. (PLOT(2).EQ.O_X) ) THEN
         IF ( C_END_GET(O_X).GT.1.0 ) THEN
            CALL PGSCI(15)
            XPTS(1) = C_START_GET(O_M)
            XPTS(POLYDIM) = C_END_GET(O_M)
            YPTS(1) = C_END_GET(O_X)
            YPTS(POLYDIM) = C_END_GET(O_X)
            DO N = 2,POLYDIM-1
               XPTS(N) = (C_END_GET(O_M) - C_START_GET(O_M))/
     $              (POLYDIM - 1)*REAL(N-2) + C_START_GET(O_M)
               YPTS(N) = BIGBANG(XPTS(N))
            ENDDO
            
            CALL PGSLS(1)
            CALL PGPOLY (POLYDIM, XPTS, YPTS)
            CALL PGSCI(1)
            CALL PGFUNX(BIGBANG,500,C_START_GET(O_M),
     $           C_END_GET(O_M),1)
            CALL PGMTXT('T',-1.5,0.05,0.0,'No Big Bang')
         ENDIF
      ENDIF

      RETURN
      END
C
C   End of plot.f
C
