C     ****h* snalys/array [0.9a]
C
C   NAME
C     ARRAY -- package to handle multidimensional arrays
C   DESCRIPTION
C     A package to handle multidimensional ARRAYs in fortran, where
C     the array is stored in one vector, which has the following
C     structure in the case of three indices:
C
C        **1**1**1**1**1**1**1**1**1**1**1**1**1**1**1**1*
C        * *11*11*11  *21*21*21 *31*31*31 *41*41*41      *
C        * *  111  *  *  121  * *  131  * *  141  *      *
C        * *  211  *  *  221  * *  231  * *  241  *  ... *
C        * *  ...  *  *  ...  * *  ...  * *  ...  *      *
C        * *11*11*11  *21*21*21 *31*31*31 *41*41*41      *
C        **1**1**1**1**1**1**1**1**1**1**1**1**1**1**1**1*
C
C        **2**2**2**2**2**2**2**2**2**2**2**2**2**2**2**2*
C        * *12*12*12  *22*22*22 *32*32*32 *42*42*42      *
C        * *  112  *  *  122  * *  132  * *  142  *      *
C        * *  212  *  *  222  * *  232  * *  242  *  ... *
C        * *  ...  *  *  ...  * *  ...  * *  ...  *      *
C        * *12*12*12  *22*22*22 *32*32*32 *42*42*42      *
C        **2**2**2**2**2**2**2**2**2**2**2**2**2**2**2**2*
C
C     that is the order of the storage is
C         111 211 ... 121 221 ... 112 212 ...
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-03-03
C
C     ***

C     
C     SELECTORS
C     =========
C     

C     ****f* array/GET_VAL_ARRAY *
C
C   NAME
C     GET_VAL_ARRAY() -- retrieve the value in the array for an index.
C
C   DESCRIPTION
C     Returns the value in the array for a specific index vector, e.g.
C     (3, 4, 5, 6). The maximum dimension of the array is 20.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-03-03
C
C   USAGE
C      GET_VAL_ARRAY(DIM, ARRAY, LEN, AIDX)
C
C   INPUTS
C     DIM     - INTEGER that specifies the dimension of the array.
C     ARRAY   - The REAL vector that contains the array structure.
C     LEN     - INTEGER vector of at least length DIM, that holds
C               the number of elements in the array for each dimension.
C     AIDX    - INTEGER vector of at least length DIM.
C
C   RESULTS
C     The array value of type REAL is returned.
C
C     ***
      REAL FUNCTION GET_VAL_ARRAY( DIM , ARRAY , LEN , IDX )
      IMPLICIT NONE
      REAL ARRAY(*)
      INTEGER LEN(*), IDX(*), GET_VECTOR_IDX, DIM

      GET_VAL_ARRAY = ARRAY(GET_VECTOR_IDX( DIM , LEN , IDX ))
      RETURN
      END
C
C End of GET_VAL_ARRAY()
C

C     ****i* array/GET_NORMAL_IDX *
C
C   NAME
C     GET_NORMAL_IDX() -- retrieve an index vector
C
C   DESCRIPTION
C     Converts an vector index (e.g. 1535) to an array index vector,
C     e.g. ( 3 , 8 , 1, 4 ). The maximal dimension of the arrays is 20.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-03-03
C
C   USAGE
C      GET_NORMAL_IDX(DIM, VIDX, LEN, AIDX)
C
C   INPUTS
C     DIM     - INTEGER that specifies the dimension of the array.
C     VIDX    - Scalar INTEGER vector index.
C     LEN     - INTEGER vector of at least length DIM, that holds
C               the number of elements in the array for each dimension.
C     AIDX    - INTEGER vector of at least length DIM.
C
C   SIDE EFFECTS
C     AIDX will contain the index vector up on return.
C
C   ERRORS
C     - If one of the elements between 1 and DIM in LEN is 0.
C     - If one of the elements in AIDX turns out to be zero (this can
C       in practice never happen).
C     
C     ***
      SUBROUTINE GET_NORMAL_IDX(DIM, VIDX, LEN, IDX)
      IMPLICIT NONE
      INTEGER VIDX, LEN(*), IDX(*), DIM

      INTEGER MAX_DIM
      PARAMETER ( MAX_DIM = 20 )

      INTEGER N,I,TMP,SIZE(MAX_DIM)
      LOGICAL OK

      OK = .TRUE.
C
C        Calculate the size of each block.
C
      DO N = 1,DIM
         IF ( N.EQ.1 ) THEN
            SIZE(N) = 1
         ELSE
            SIZE(N) = LEN(N-1)*SIZE(N-1)
         ENDIF
         IF (LEN(N).EQ.0) OK = .FALSE.
      ENDDO

      I = VIDX

      IF ( OK ) THEN
         DO N = 1,DIM
            TMP = INT( REAL(I-1)/REAL(SIZE(DIM+1-N)) )
            
            I = I - TMP*(SIZE(DIM+1-N))
            IF (TMP.GE.0) THEN
               IDX(DIM+1-N) = TMP + 1
            ELSE
               WRITE (*,'(2A)') 'error: GET_NORMAL_IDX this can never ',
     $              'happen!'
               IDX(DIM-N+1) = 0
            ENDIF
         ENDDO
      ELSE
         WRITE (*,'(2A)') 'error: GET_NORMAL_IDX illegal tensor ',
     $        'dimension zero'
      ENDIF
      END
C
C End of GET_NORMAL_IDX()
C



C     ****i* array/GET_VECTOR_IDX *
C
C   NAME
C     GET_VECTOR_IDX() -- retrieve an index in the vector
C
C   DESCRIPTION
C     Converts an array index vector, e.g. ( 3 , 8 , 1, 4 ) to a scalar
C     vector (e.g. 1535). The maximal dimension of the arrays is 20.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-03-03
C
C   USAGE
C      GET_VECTOR_IDX(DIM, LEN, AIDX)
C
C   INPUTS
C     DIM     - INTEGER that specifies the dimension of the array.
C     LEN     - INTEGER vector of at least length DIM, that holds
C               the number of elements in the array for each dimension.
C     AIDX    - INTEGER vector of at least length DIM.
C
C   RESULTS
C     The scalar vector index is returned.
C
C     ***
      INTEGER FUNCTION GET_VECTOR_IDX (DIM, LEN, IDX)
      IMPLICIT NONE
      INTEGER LEN(*), IDX(*), DIM

      INTEGER MAX_DIM
      PARAMETER ( MAX_DIM = 20 )

      INTEGER SIZE(MAX_DIM),N

C
C        Calculate the size of each block.
C
      DO N = 1,DIM
         IF ( N.EQ.1 ) THEN
            SIZE(N) = 1
         ELSE
            SIZE(N) = LEN(N-1)*SIZE(N-1)
         ENDIF
      ENDDO

C
C        Find where in the array structure the value should
C        be stored.
C
      GET_VECTOR_IDX = 1
      DO N = 1,DIM
         GET_VECTOR_IDX = GET_VECTOR_IDX + (IDX(N) - 1)*SIZE(N)
      ENDDO
      RETURN
      END
C
C End of GET_VECTOR_IDX()
C

C
C  MODIFIERS
C  =========
C

C     ****f* array/PUT_VAL_ARRAY *
C
C   NAME
C     PUT_VAL_ARRAY() -- stores a value in the array
C
C   DESCRIPTION
C     Stores a value in the array at a specified index. The maximal
C     dimension of the array is 20.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-03-03
C
C   USAGE
C     PUT_VAL_ARRAY(VAL, ARRAY, DIM , LEN , AIDX)
C
C   INPUTS
C     VAL     - REAL value that is to be put in the array.
C     ARRAY   - REAL vector that is the representation of the array.
C     DIM     - INTEGER that specifies the dimension of the array.
C     LEN     - INTEGER vector of at least length DIM, that holds
C               the number of elements in the array for each dimension.
C     AIDX    - INTEGER index vector of at least length DIM.
C
C   SIDE EFFECTS
C     VAL is stored in the array.
C
C     ***
      SUBROUTINE PUT_VAL_ARRAY( VAL , ARRAY , DIM , LEN , IDX )
      IMPLICIT NONE
      REAL ARRAY(*), VAL
      INTEGER DIM, LEN(*), IDX(*)
      INTEGER GET_VECTOR_IDX
C
C        Store the value.
C
      ARRAY(GET_VECTOR_IDX(DIM, LEN, IDX)) = VAL
      END
C
C End of PUT_VAL_ARRAY()
C

C
C End of ARRAY
C
