C     ****h* snalys/cos_int [0.9a]
C
C   NAME
C     cos_int -- handle the cosmology matrix
C
C   DESCRIPTION
C     The intervals for the grid search and the number of steps
C     are stored in a matrix ("the cosmology matrix"). This module
C     defines a number of selectors and modifiers to handle the matrix.
C
C     The rows of the matrix represent different parameters which and
C     the columns have the following structure:
C
C     1     2     3       4     5
C     begin  end  nrstep  step  count
C     ...    ...  ...     ...   ...
C     ...    ...  ...     ...   ...
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-03-02
C
C   NOTES
C     The cosmology matrix, cos_int, is part of a common block
C     defined in the file 'cosmo.inc'
C
C   BUGS
C     All if-checks, e.g. that the indices are not out of bounds, have
C     been stripped to gain speed. This may (but should not) cause the
C     program to crash if one is not beeing careful.
C
C     ***

C
C     SELECTORS
C     =========
C

C     ****f* snalys/C_START_GET [1.0]
C
C
C   NAME
C     C_START_GET -- returns the start value from the cosmology matrix.
C
C   DESCRIPTION
C     Returns the start value for the grid search, for a specific
C     cosmological parameter.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-03-02
C
C   USAGE
C     C_START_GET( PAR )
C
C   INPUTS
C     PAR - an integer that represents a cosmological parameter.
C
C   RESULT
C     Returns a real value, the start value.
C
C   BUGS
C     If the input integer is out of range the program crashes.
C
C   HISTORY
C     2000-05-25 The if statments were removed to increase the speed
C                of the program.
C
C     ***
      REAL FUNCTION C_START_GET (N)
      IMPLICIT NONE
      INTEGER N
      
      INCLUDE 'cosmo.inc'

      C_START_GET = COS_INT(N,1)
      RETURN
      END
C
C  End of C_START_GET()
C


C     ****f* snalys/C_END_GET [1.0]
C
C
C   NAME
C     C_END_GET -- returns the end value from the cosmology matrix.
C
C   DESCRIPTION
C     Returns the end value for the grid search, for a specific
C     cosmological parameter.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-03-02
C
C   USAGE
C     C_END_GET( PAR )
C
C   INPUTS
C     PAR - an integer that represents a cosmological parameter.
C
C   RESULT
C     Returns a real value, the end value.
C
C   BUGS
C     If the input integer is out of range the program crashes.
C
C   HISTORY
C     2000-05-25 The if statments were removed to increase the speed
C                of the program.
C
C     ***
      REAL FUNCTION C_END_GET (N)
      IMPLICIT NONE
      INTEGER N
      INCLUDE 'cosmo.inc'

      C_END_GET = COS_INT(N,2)
      RETURN
      END
C
C End of C_END_GET()
C


      
C     ****f* snalys/C_NR_GET [1.0]
C
C
C   NAME
C     C_NR_GET -- returns the number of grid steps from the
C                 cosmology matrix.
C   DESCRIPTION
C     Returns the number of grid steps for a specific
C     cosmological parameter.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-03-02
C
C   USAGE
C     C_NR_GET( PAR )
C
C   INPUTS
C     par - an integer that represents a cosmological parameter.
C
C   RESULT
C     Returns a the number of steps as a real value (although it
C     actually is an integer).
C
C   BUGS
C     If the input integer is out of range the program crashes.
C
C   HISTORY
C     2000-05-25 The if statments were removed to increase the speed
C                of the program.
C
C     ***
      REAL FUNCTION C_NR_GET (N)
      IMPLICIT NONE
      INTEGER N
      INCLUDE 'cosmo.inc'
      
      C_NR_GET = COS_INT(N,3)
      RETURN
      END
C
C End of C_NR_GET()
C



C     ****f* snalys/C_STEP_GET [1.0]
C
C
C   NAME
C     C_STEP_GET -- returns the number of grid step length from the
C                   cosmology matrix.
C   DESCRIPTION
C     Returns the grid step length for a specific cosmological
C     parameter.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-03-02
C
C   USAGE
C     C_STEP_GET( PAR )
C
C   INPUTS
C     PAR - an integer that represents a cosmological parameter.
C
C   RESULT
C     Returns a the grid step length as a real value.
C
C   BUGS
C     If the input integer is out of range the program crashes.
C
C   HISTORY
C     2000-05-25 The if statments were removed to increase the speed
C                of the program.
C     
C     ***
      REAL FUNCTION C_STEP_GET (N)
      IMPLICIT NONE
      INTEGER N
      INCLUDE 'cosmo.inc'

      C_STEP_GET = COS_INT(N,4)
      RETURN
      END
C
C End of C_STEP_GET()
C



     
C
C   MUTATORS
C   ========
C

C     ****f* snalys/C_START_PUT [1.0]
C
C
C   NAME
C     C_START_PUT -- sets the start value for the grid search.
C
C   DESCRIPTION
C     This subroutine is used to set the start value of a cosmological
C     parameter for the grid search.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-03-02
C
C   USAGE
C     CALL C_START_PUT( PAR,VALUE )
C
C   INPUTS
C     PAR   - an integer that represents a cosmological parameter.
C     VALUE - the value of set (of type real).
C
C   SIDE EFFECTS
C     The element (n,1) in the cosmology matrix cos_int is set to
C     value.
C
C   WARNINGS
C     If the integer argument is out of array index bounds, a warning
C     message is written.
C
C     ***
      SUBROUTINE C_START_PUT (N, VAL)
      IMPLICIT NONE
      INTEGER N
      REAL VAL
      INCLUDE 'cosmo.inc'

      IF ( (N.GE.1).AND.(N.LE.NR_PAR) ) THEN
         COS_INT(N,1) = VAL
      ELSE
         WRITE (*,'(A)') 'warning: array index out of bounds'
      ENDIF      
      END
C
C End of C_START_PUT
C

C     ****f* snalys/C_END_PUT [1.0]
C
C
C   NAME
C     C_END_PUT -- sets the end value for the grid search.
C
C   DESCRIPTION
C     This subroutine is used to set the end value of a cosmological
C     parameter for the grid search.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-03-02
C
C   USAGE
C     CALL C_END_PUT( PAR,VALUE )
C
C   INPUTS
C     PAR   - An integer that represents a cosmological parameter.
C     VALUE - The value of set (of type real).
C
C   SIDE EFFECTS
C     The element (n,2) in the cosmology matrix cos_int is set to
C     value.
C
C   WARNINGS
C     If the integer argument is out of array index bounds, a warning
C     message is written.
C
C     ***
      SUBROUTINE C_END_PUT (N, VAL)
      IMPLICIT NONE
      INTEGER N
      REAL VAL
      INCLUDE 'cosmo.inc'
      INCLUDE 'error.inc'

      IF ( (N.GE.1).AND.(N.LE.NR_PAR) ) THEN
         COS_INT(N,2) = VAL
      ELSE       
         WRITE (STDERR,'(A)') 'warning: array index out of bounds'
      ENDIF
      END
C
C End of C_END_PUT()
C




C     ****f* snalys/C_NR_PUT [1.0]
C
C
C   NAME
C     C_NR_PUT -- sets the number of grid search steps.
C
C   DESCRIPTION
C     This subroutine is used to set the number of grid search steps
C     for a specific cosmological parameter. The grid step length is
C     also calculated and set.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-03-02
C
C   USAGE
C     CALL C_NR_PUT( PAR,VALUE )
C
C   INPUTS
C     PAR   - an integer that represents a cosmological parameter.
C     VALUE - the value of set (of type real).
C
C   SIDE EFFECTS
C     The element (n,3) in the cosmology matrix cos_int is set to
C     value, and the element (n,4) is also calculated and set.
C
C   WARNINGS
C     A warning message is written if the number of steps value is
C     negative, and the value is changed to 1.
C
C     A warning message is also written when the given integer index
C     is out of bounds.
C
C     ***
      SUBROUTINE C_NR_PUT (N, VAL)
      IMPLICIT NONE
      INTEGER N
      REAL VAL, C_END_GET, C_START_GET, TMP
      INCLUDE 'cosmo.inc'
      INCLUDE 'error.inc'

      IF ( VAL.LT.1.0 ) THEN
         WRITE (STDERR,'(A)') 'warning: the nr of steps must be >=1'
         VAL = 1.0
      ENDIF

      IF ( (N.GE.1).AND.(N.LE.NR_PAR) ) THEN
         COS_INT(N,3) = VAL
         TMP = (C_END_GET(N) - C_START_GET(N))/VAL
         IF (VAL.GT.1) THEN
            CALL C_STEP_PUT (N,(C_END_GET(N)-C_START_GET(N))/(VAL-1.0))
         ELSE
            CALL C_STEP_PUT (N,(C_END_GET(N)-C_START_GET(N))/VAL)
         ENDIF
      ELSE
         WRITE (STDERR,'(A)') 'warning: array index out of bounds'
      ENDIF
      END
C
C End of C_NR_PUT
C







C     ****f* snalys/C_STEP_PUT [1.0]
C
C
C   NAME
C     C_STEP_PUT -- sets the grid search step length.
C
C   DESCRIPTION
C     This subroutine is used to set the grid search step length
C     for a specific cosmological parameter.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-03-02
C
C   USAGE
C     CALL C_STEP_PUT( PAR,VALUE )
C
C   INPUTS
C     PAR   - an integer that represents a cosmological parameter.
C     VALUE - the value of set (of type real).
C
C   SIDE EFFECTS
C     The element (n,4) in the cosmology matrix cos_int is set to
C     value.
C
C   WARNINGS
C     A warning message is written when the given integer index
C     is out of bounds.
C
C     ***
      SUBROUTINE  C_STEP_PUT (N, VAL)
      IMPLICIT NONE
      INTEGER N
      REAL VAL
      INCLUDE 'cosmo.inc'
      INCLUDE 'error.inc'

      IF ( VAL.NE.0 ) THEN
         IF ( (N.GE.1).AND.(N.LE.NR_PAR) ) THEN
            COS_INT(N,4) = VAL
         ELSE
            WRITE (STDERR,'(A)') 'warning: array index out of bounds'
         ENDIF
      ELSE
         WRITE (STDERR,'(A)') 'warning: step length set to zero.'
      ENDIF
      END
C
C End of C_STEP_PUT()
C



C--
C--   THIS SHOULD BE MOVED
C--

      REAL FUNCTION IDX_TO_VAL(I,N,NR_VAR)
      IMPLICIT NONE
      INTEGER I, N, NR_VAR

      INTEGER IDX(100), LEN(100), K
      REAL C_NR_GET, C_START_GET, C_STEP_GET

      DO K = 1, NR_VAR
cc       The length of the vectors in the tensor f_den.
         LEN(K) = INT(C_NR_GET(K))
      ENDDO

      CALL GET_NORMAL_IDX(NR_VAR,N,LEN,IDX)
      IDX_TO_VAL = C_START_GET(I) + REAL(IDX(I)-1)*C_STEP_GET(I)

      END

C-- EOF

