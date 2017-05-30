C
C MAIN PROGRAM
C
      PROGRAM SNALYS
      IMPLICIT NONE

      CHARACTER*(256) FNAME

      INTEGER N

      LOGICAL OK, ANALYZE, INITIALIZE

      INCLUDE 'cosmo.inc'
      INCLUDE 'operation.inc'
      INCLUDE 'io.inc'
      INCLUDE 'error.inc'

      REAL C_START_GET, C_END_GET
      REAL P_COS(NR_PAR+1)

C  
C        Take care of all the input arguments, read the ini-file, read the
C        input data file.
C  
      OK = INITIALIZE( P_COS, FNAME )
 
      IF ( OK .AND. NR_VAR.LE.NR_PAR ) THEN

C
C           Do the estimations, fits and what ever...
C
         OK = ANALYZE( P_COS, OFID, FNAME )

C
C        Print a help message to standard output.
C
      ELSE
         WRITE (STDOUT,'(2A)') 'snalys [-q -o outfile] -i parfile ',
     $        'datafile'
      Endif
      STOP
      END
C
C End of MAIN PROGRAM
C



C
C End of main.f
C
