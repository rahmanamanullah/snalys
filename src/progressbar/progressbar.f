C     ****f* progressbar.f/PROGRESS [1.0]
C
C
C   NAME
C     PROGRESS -- prints a progressbar
C
C   DESCRIPTION
C     Prints a progressbar. This routine should be called every time
C     the progressbar should be updated.
C
C   AUTHOR
C     Tomas Oppelstrup (tomaso@nada.kth.se)
C
C   CREATION DATE
C     2000-03-03
C
C   USAGE
C     CALL PROGRESS(PERCENT)
C
C   INPUTS
C     PERCENT - the percentage of the process that has been finished.
C
C     ***
      SUBROUTINE PROGRESS(PERCENT)
      IMPLICIT NONE
      REAL PERCENT
      INTEGER N, K, BAR_LENGTH
      PARAMETER (BAR_LENGTH = 60)
      CHARACTER EQUAL(BAR_LENGTH),DOT(BAR_LENGTH)
      DATA EQUAL/BAR_LENGTH*'='/,DOT/BAR_LENGTH*'.'/

      N = NINT( BAR_LENGTH*PERCENT )
      
      IF (N.EQ.0) THEN
         Write (*,'(61A,I3,A)') '|>',(DOT(K),K=1,BAR_LENGTH-1),'|   ',
     $        NINT(PERCENT*100),'%'
      ELSE IF (N.EQ.BAR_LENGTH) THEN
         WRITE (*,'(61A,I3,A)') '|',(EQUAL(K),K=1,BAR_LENGTH-1),'>|   ',
     $        NINT(PERCENT*100),'%'
      ELSE
         WRITE (*,'(62A,I3,A)') '|',(EQUAL(K),K=1,N-1),'>',
     $        (DOT(K),K=N+1,BAR_LENGTH),'|   ',NINT(PERCENT*100),'%'
      ENDIF
      IF (PERCENT.NE.1.0) WRITE (*,'(A)') CHAR(27)//'[A'//CHAR(27)//'[A'
      END
C
C End of progressbar.f
C
