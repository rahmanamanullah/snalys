C     ****f* snalys/RESFITS *
C
C   NAME
C     RESFITS -- Virtual RESFITS routine
C
C   DESCRIPTION
C     This is a version of the RESFITS routine that will be
C     called if cfitsio is not available.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2003-09-24
C
C   USAGE
C     CALL RESFITS( STIME, F_DEN, GRI, GERR )
C
C   INPUTS
C     STIME     - The start time for the execution (in unix format)
C     F_DEN     - the grid array of type REAL
C     GRI       - the best fit from the GRIDSEARCH.
C     GERR      - the estimated errors from GRIDSEARCH.
C
C     ***
      SUBROUTINE RESFITS( STIME, F_DEN, GRI, GERR )
      IMPLICIT NONE
      INTEGER*4 STIME, TIME
      EXTERNAL TIME
      REAL F_DEN(*), GRI(*), GERR(2,*)

      include 'error.inc'

      WRITE (STDERR,'(2A)') 'resfits: this version was not compiled ',
     $     'with the cfitsio library.'

      END
C
C END OF RESFITS
C
