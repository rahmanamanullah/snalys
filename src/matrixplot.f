C     ****f* snalys/MATRIXPLOT *
C
C   NAME
C     MATRIXPLOT -- create a contour plot for a matrix.
C
C   DESCRIPTION
C     Plots different contour levels for a matrix.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-10-25
C
C   USAGE
C     MATRIXPLOT( MAT, DIM, C, NRCONT, COLOR, LS, LABEL, PLOT, DISP )
C
C   INPUTS
C     MAT    - the matrix to be plotted
C     DIM    - the size of the matrix
C     C      - the contour levels
C     NRCONT - number of contours levels to be plotted
C     COLOR  - the colour for the contours (see PGPLOT documentation)
C     LS     - the line style (see PGPLOT documentation)
C     LABEL  - the labels for the plot (see PGPLOT documenation)
C     PLOT   - the indices of the two variables that will be plotted.
C     DISP   - do not know what this is but see PGPLOT documentation.
C
C   SEE ALSO
C     The PGPLOT documentation.
C
C     ***
      SUBROUTINE MATRIXPLOT( MAT, DIM, C, NRCONT, COLOR, LS, LABEL,
     $     PLOT, DISP )
      IMPLICIT NONE
      REAL MAT(*), C(*), DISP
      INTEGER DIM(*), NRCONT, COLOR, LS, PLOT(*)
      CHARACTER*(*) LABEL

      INTEGER MAXDIM
      PARAMETER ( MAXDIM = 1000 )

      INTEGER IDX(2), GET_VECTOR_IDX, N, I
      REAL CONT(MAXDIM,MAXDIM), TR(6)
      REAL C_START_GET, C_STEP_GET
C  
C              Convert the matrix mat to a 'real' matrix.
C  
      IF ( DIM(1).LE.MAXDIM .AND. DIM(2).LE.MAXDIM ) THEN
         DO N = 1, DIM(2)
            DO I = 1, DIM(1)
               IDX(1) = I
               IDX(2) = N
               CONT(I,N) = MAT(GET_VECTOR_IDX(2,DIM,IDX))
c             write (*,*) cont(i,n), REAL((i-1))/real(9),
c    $              REAL((n-1))/real(9)
            ENDDO
         ENDDO

C  
C              The plotting transformation matrix.
C              (see the pgplot documentation for pgcons)
C  
         TR(1) = C_START_GET(PLOT(1)) - C_STEP_GET(PLOT(1))
         TR(2) = C_STEP_GET(PLOT(1))
         TR(3) = 0.0
         TR(4) = C_START_GET(PLOT(2)) - C_STEP_GET(PLOT(2))
         TR(5) = 0.0
         TR(6) = C_STEP_GET(PLOT(2))
         
         CALL PGSCI(COLOR)
         CALL PGMTXT('B',DISP,0.95,1.0,LABEL)
         CALL PGSLS(LS)          ! Set the line style.
         CALL PGCONT(CONT,MAXDIM,MAXDIM,1,DIM(1),1,DIM(2),C,
     $        NRCONT,TR)
         CALL PGSCR(0,1.,1.,1.) ! Set the background colour.
         CALL PGSCR(1,0.,0.,0.) ! Set the the pen colour.
      ENDIF
      END
C
C   End of MATRIXPLOT.
C
