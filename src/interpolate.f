C     ****h* snalys/interpolate.f [1.0]
C
C   NAME
C     Interpolate -- a set of interpolation related routines.
C
C   DESCRIPTION
C     This package contains a set of routines to perform interpolation.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-10-13
C
C     ***


C     ****i* interpolate/LOCATEINTER
C
C
C   NAME
C     LOCATEINTER -- Locate the position in a vector for interpolation.
C
C   DESCRIPTION
C     Given an array xx(1:n), and given a value x, returns a value j
C     such that x is between xx(j) and xx(j+1). xx(1:n) must be monotonic,
C     either increasing or decreasing. j=0 or j=n is returned to indicate
C     that x is out of range.
C
C   AUTHOR
C     Numerical Recepies (www.nr.com)
C
C   CREATION DATE
C     2001-10-13
C
C   USAGE
C     CALL LOCATEINTER( XX , N , X , J )
C
C   INPUTS
C     XX - REAL array containing the monotonic values.
C     N  - length of XX.
C     X  - the REAL value for which the location in XX will
C          be found.
C     J  - INTEGER
C
C   SIDE EFFECTS
C     On return J will contain the index so that X is located between
C     XX(J) and XX(J+1)
C
C     ***
      SUBROUTINE LOCATEINTER(XX,N,X,J)
      IMPLICIT NONE
      INTEGER J,N
      REAL X,XX(N)
      INTEGER JL,JM,JU
      JL = 0                         ! initalize lower
      JU = N + 1                     ! and upper limits
 10   IF ( JU-JL .GT. 1 ) THEN       ! if we are not yet done,
         JM = (JU + JL)/2            ! compute a midpoint
         IF ( (XX(N).GE.XX(1)) .EQV. (X.GE.XX(JM)) ) THEN
            JL = JM                  ! and replace either the lower limit
         ELSE
            JU = JM                  ! or the upper limit, as appropriate.
         ENDIF
         GOTO 10                     ! Repeat until
      ENDIF                          ! the test condition 10 is satisfied.
      IF ( X .EQ. XX(1) ) THEN
         J = 1
      ELSE IF ( X .EQ. XX(N) ) THEN
         J = N - 1
      ELSE
         J = JL
      ENDIF
      RETURN                         ! and return.
      END
C
C  End of LOCATEINTER
C


C
C     ****f* interpolate/SIMPBIINT
C
C
C   NAME
C     SIMPBIINT -- Simple bilinear interpolation in a grid.
C
C   DESCRIPTION
C     A simple bilinear interpolation in a grid. Given arrays X1A(1:M)
C     and X2A(1:N) of independent variables, and an M by N array of
C     function values YA(1:M,1:N), tabulated at the grid points defined
C     by X1A and X2a; and given values X1 and X2 of the independent
C     variables; this routine returns an interpolated function value Y.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-10-13
C
C   USAGE
C     CALL SIMPBINT( X1A , X2A , YA , M , N , X1 , X2 , Y )
C
C   INPUTS
C     X1A - REAL array of length M containing the monotonic
C           values for the first variable.
C     X2A - REAL array of length N containing the monotonic
C           values for the second variable.
C     YA  - REAL array of dimensions MxN containing the function
C           values.
C     X1  - value of the first variable.
C     X2  - value of the second variable.
C     Y   - REAL
C
C   SIDE EFFECTS
C     On return Y will be the interpolated function value in the point
C     (X1,X2) if the point is within the table boundaries.
C
C     ***
      SUBROUTINE SIMPBIINT(X1A,X2A,YA,M,N,X1,X2,Y)
      IMPLICIT NONE
      INTEGER M,N
      REAL X1,X2,Y,X1A(M),X2A(N),YA(M,N)

      INTEGER J,K
      REAL T,U
C
C        Find the position of the point in the grid.
C
      CALL LOCATEINTER(X1A,M,X1,J)
      CALL LOCATEINTER(X2A,N,X2,K)

C
C        Check that the given X1 and X2 are within table boundaries.
C
      IF ( (J.GT.0) .AND. (J.LT.M) .AND. (K.GT.0) .AND. (K.LT.N) ) THEN

C
C           Both T and U are normalized to lie between 0 and 1.
C      
         T = (X1 - X1A(J))/(X1A(J+1) - X1A(J))
         U = (X2 - X2A(K))/(X2A(K+1) - X2A(K))
      
C
C           Calculate the interpolated value.
C      
         Y = (1 - T)*(1 - U)*YA(J,K) + T*(1 - U)*YA(J+1,K) +
     $        T*U*YA(J+1,K+1) + (1 - T)*U*YA(J,K+1)
      ELSE
         
      ENDIF
      RETURN
      END
C
C End of SIMPBIINT
C
