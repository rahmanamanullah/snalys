C     ****i* contourfinder/NORMVEC [1.0] *
C
C   NAME
C     NORMVEC() -- normalizes a vector.
C
C   DESCRIPTION
C     Normalizes a vector so that the length of it equals one.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-08-20
C
C   USAGE
C     NORMVEC( VEC, START, LENGTH )
C
C   INPUTS
C     VEC    - a vector of type REAL that contains at least N element and
C              is not the null-vector.
C     START  - the offset in VEC where the first element of the vector is
C              stored. This must be less than LENGTH.
C     LENGTH - the length of the vector.
C
C   SIDE EFFECTS
C     The vector VEC is replaced with the normalized vector. If the
C     subroutine was unsuccessful LENGTH will equal -1 in return.
C
C     ***
      SUBROUTINE NORMVEC(VECTOR,START,LENGTH)
      IMPLICIT NONE
      REAL VECTOR(*)
      INTEGER START, LENGTH

      INTEGER N
      REAL NORM, EPS
      PARAMETER ( EPS = 1E-12 )
      
      IF ( START.LE.LENGTH ) THEN
         NORM = 0.0
         DO N = 1, LENGTH
            NORM = NORM + VECTOR(START - 1 + N)**2
         ENDDO
         NORM = SQRT(NORM)

C
C           Check that the given vector is not the null vector.
C
         IF ( NORM.GT.EPS ) THEN
            DO N = 1, LENGTH
               VECTOR(START - 1 + N) = VECTOR(START - 1 + N)/NORM
            ENDDO
         ELSE
            LENGTH = -1
         ENDIF
      ELSE
         LENGTH = -1
      ENDIF
      END
C     *** End of SUBROUTINE NORMVEC() ***

C     END OF FILE
