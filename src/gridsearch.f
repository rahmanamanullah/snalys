C     ****h* snoc/grid_search.f [1.0]
C
C   NAME
C     Grid Search -- a brute force minimization package.
C   DESCRIPTION
C     This package minimizes a function by calculating its value for
C     different input arguments and returns the arguments that give
C     the smallest function value.
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C   CREATION DATE
C     2000-03-02
C
C     ***


C     ****f* grid_search/GRID_SEARCH [1.0]
C
C
C   NAME
C     GRID_SEARCH -- a brute force minimization routine.
C
C   DESCRIPTION
C     This routine minimizes a function by calculating its value for
C     different input arguments and returns the arguments that give
C     the smallest function value.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-03-02
C
C   USAGE
C     GRID_SEARCH( FUNC, NR_VAR, TOT_LEN, IDX_TO_VAL, PC, RESULT, BEQUIET )
C
C   INPUTS
C     func(c)          - is a function of type DOUBLE PRECISION to be
C                        minimized that takes two REAL vectors as input
C                        arguments.
C     THENRVAR         - the number of variables that func depends on (<100),
C                        i.e. the minimum vector length of c above.
C     TOT_LEN          - the total number of steps/points to be tested.
C     IDX_TO_VAL( I, N, THENRVAR )
C                        - a real function that given a value, n, between 1
C                        and tot_len and a variable number, i, between 1
C                        and THENRVAR returns the value of the variable i
C                        representing n. THENRVAR has the same meaning
C                        as above.
C     PC               - a real vector of at least length THENRVAR+1.
C     RESULT           - a real vector of at least length tot_len.
C     BEQUIET          - a logical variable that controles whether a
C                        progress bar should be printed to standard output.
C   SIDE EFFECTS
C     The vector result contains all function values calculated, and the
C     vector pc will when the upon return consist the minimum point and
C     the element pc(THENRVAR+1) will be the minimum function value.
C
C     The parameters BEST_SCRIPTM and BEST_ALPHA will be set if
C     if the OPER(FITSCM) and OPER(FITSALP) have been set.
C
C   BUGS
C     The program crashes if tot_len is bigger than the dimension of
C     result, or if THENRVAR is bigger than the dimensions of pc or 100.
C
C   USES
C     progressbar.f, cosmo.inc, operation.inc
C
C   HISTORY
C     2004-10-27 - Added the savings of the values for the nuisance
C                  parameters for the best fit.
C
C     ***
      SUBROUTINE GRID_SEARCH( FUNC, THENRVAR, TOT_LEN, IDX_TO_VAL, PC,
     $     RESULT, BEQUIET )
      IMPLICIT NONE
      REAL  IDX_TO_VAL, PC(*), RESULT(*)
      LOGICAL BEQUIET
      EXTERNAL IDX_TO_VAL
      REAL FUNC
      INTEGER TOT_LEN, THENRVAR
      EXTERNAL FUNC

      INTEGER N, I
      REAL C(100), C_NR_GET

      INCLUDE 'cosmo.inc'
      INCLUDE 'operation.inc'

C
C        If not quiet-mode, initialize the progress bar.
C
      IF (.NOT.BEQUIET) CALL PROGRESS(0.0)

C
C        Perform the grid search for all search points.
C
      DO N=1,TOT_LEN
C
C           Calculate the point.
C
         DO I=1,THENRVAR
            C(I) = IDX_TO_VAL( I, N, THENRVAR )
         ENDDO

C
C           Calculate the function value in the point C.
C
         RESULT(N) = REAL( FUNC( C, C ) )

C
C           Store the point, where the minium function value
C           is calculated.
C
         IF ( RESULT(N).LT.PC(THENRVAR+1) ) THEN
            PC(THENRVAR+1) = RESULT(N)
            IF ( OPER(FITMSC) ) THEN
               BEST_SCRIPTM = THIS_SCRIPTM
            ENDIF
            IF ( OPER(FITSALP) ) THEN
               BEST_ALPHA = THIS_ALPHA
            ENDIF
            DO I = 1, THENRVAR
               PC(I) = C(I)
            ENDDO
         ENDIF

C
C           If not in quiet-mode, update the progress bar.
C
         IF (.NOT.BEQUIET) CALL PROGRESS(REAL(N)/REAL(TOT_LEN))
      ENDDO

      END
C
C   End of grid_search()
C

C     ****i* grid_search/progress [1.0]
C
C
C   NAME
C     progress -- prints a progress bar to standard output.
C
C   DESCRIPTION
C     This routine prints a progressbar to standard output, given
C     a percent value representing the current progress status.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-05-25
C
C   USAGE
C     progress(percent)
C
C   INPUTS
C     percent - a real value reprecenting the current progress status
C               in percent.
C
C   HISTORY
C     2001-10-14 Changed so that the progressbar only is updated with
C                there is a change in the printed percentage.
C
C     ***
      SUBROUTINE PROGRESS(PERCENT)
      IMPLICIT NONE
      REAL PERCENT
      INTEGER N, K, BAR_LENGTH, PERCOLD, PERC, NOLD
      COMMON /PROGRESSBAR_OLD/PERCOLD,NOLD
      PARAMETER (BAR_LENGTH = 60)
      CHARACTER EQUAL(BAR_LENGTH),DOT(BAR_LENGTH)
      DATA EQUAL/BAR_LENGTH*'='/,DOT/BAR_LENGTH*'.'/
      DATA PERCOLD/-1/,NOLD/-1/

      PERC = NINT( PERCENT*100 )
      N = NINT( BAR_LENGTH*PERCENT )
      
      IF ( PERCOLD.NE.PERC .OR. N.NE.NOLD .OR. PERCENT.EQ.1.0 ) THEN
         IF (N.EQ.0) THEN
            WRITE(6,'(61A,I3,A)') '|>',(DOT(K),K=1,BAR_LENGTH-1),'|   ',
     $           NINT(PERCENT*100),'%'
         ELSE IF (N.EQ.BAR_LENGTH) THEN
            WRITE (6,'(61A,I3,A)')
     $           '|',(EQUAL(K),K=1,BAR_LENGTH-1),'>|   ',PERC,'%'
         ELSE 
            WRITE (6,'(62A,I3,A)') '|',(EQUAL(K),K=1,N-1),'>',
     $           (DOT(K),K=N+1,BAR_LENGTH),'|   ',PERC,'%'
         ENDIF
         IF ( PERCENT.NE.1.0 )
     $        WRITE (*,'(A)') CHAR(27)//'[A'//CHAR(27)//'[A'
      ENDIF

      NOLD = N
      PERCOLD = PERC
      END
C
C  End of progress()
C

C
C  EOF gridsearch.f
C
