C     ****h* /solvedirection [1.0] *
C   NAME
C     solvedirection -- solves an equation along a specific line in a
C                       multi-dimensional equation.
C
C   DESCRIPTION
C     The package is used to solve a one-dimensional equation along a
C     specific direction on a multidimensional function.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-08-24
C
C   USED BY
C     contourfinder
C
C     ***

C     ****f* solvedir/SOLVEDIR *
C
C   NAME
C     SOLVEDIR() -- solves a one-dimensional equation along a line
C                   on a multi-dimensional function.
C
C   DESCRIPTION
C     Solves a one-dimensional equation along a specific direction on
C     a multidimensional function.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C
C   CREATION DATE
C     2000-08-23
C
C   USAGE
C     SOLVEDIR( FUNC, VAL, N, D, ACC, P, G, METHOD )
C
C   INPUTS
C     FUNC   - is a REAL function that takes three arguments: FUNC( P, G, N )
C              where P and G are REAL vectors that are at least N elements
C              long (the dimension of the function/number of variables).
C              FUNC should calculate the value of the function in the point P
C              and should also calculate the gradiant in that point and return
C              it in G as a side effect.
C     VAL    - the value for which the equation will be solved.
C     N      - the dimension of the function.
C     D      - the direction along FUNC where the equation will be solved.
C     ACC    - the accuracy that should be used.
C     P      - a REAL vector of length N, that contains the startpoint.
C     G      - a REAL vector of length N.
C     METHOD - an integer that specifies what method should be used to
C              solve the equation:
C                1   Newton-Raphsons method.
C
C   RESULT
C     A boolean value is returned that indicates if the equation was
C     solved properly.
C
C   SIDE EFFECTS
C     The solution is stored in P, and the gradient in that point is
C     stored in G.
C
C   TODO
C     * Impliment more methods.
C     * Maybe some boundary conditions should be implimented.
C
C     ***
      LOGICAL FUNCTION SOLVEDIR(FUNC, VAL, NRDIM, D, ACC, P, G, METHOD)
      IMPLICIT NONE
      REAL FUNC, D(*), G(*), P(*),  VAL, ACC
      INTEGER NRDIM, METHOD
      EXTERNAL FUNC

      LOGICAL NEWRAPDIR, HALFINTDIR
      REAL NORM, EPS
      PARAMETER ( EPS = 1.0E-10 )
      INTEGER N

C
C        Check that the direction is not the null vector.
C
      NORM = 0.0
      DO N = 1, NRDIM
         NORM = NORM + ABS(D(N))
      ENDDO

      IF ( NORM.GT.EPS ) THEN
C
C           Use Newton-Raphsons method.
C
         IF ( METHOD .EQ. 1 ) THEN
            SOLVEDIR = NEWRAPDIR(FUNC, VAL, NRDIM, D, ACC, P, G)
         ELSE IF ( METHOD .EQ. 2 ) THEN
C
C           Use a half-interval method.
C
            SOLVEDIR = HALFINTDIR(FUNC, VAL, NRDIM, D, ACC, P, G)
         ELSE
            SOLVEDIR = .FALSE.
         ENDIF
      ELSE
         SOLVEDIR = .FALSE.
      ENDIF

      RETURN
      END
C     *** End of FUNCTION SOLVEDIR() ***



C     ****i* solvedir/NEWRAPDIR *
C
C   NAME
C     NEWRAPDIR() -- solves a one-dimensional equation along a line
C                    on a multi-dimensional function using Newton-Raphsons
C                    method.
C
C   DESCRIPTION
C     Solves a one-dimensional equation along a specific direction on
C     a multidimensional function using Newton-Raphsons method.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C
C   CREATION DATE
C     2000-08-20
C
C   USAGE
C     NEWRAPDIR( FUNC, VAL, N, D, ACC, P, G )
C
C   INPUTS
C     FUNC   is a REAL function that takes three arguments: FUNC( P, G, N )
C            where P and G are REAL vectors that are at least N elements
C            long (the dimension of the function/number of variables).
C            FUNC should calculate the value of the function in the point P
C            and should also calculate the gradiant in that point and return
C            it in G as a side effect.
C     VAL    is the value for which the equation will be solved.
C     N      is dimension of the function.
C     D      is the direction along FUNC where the equation will be solved.
C     ACC    is the accuracy that should be used.
C     P      is a REAL vector of length N, that contains the startpoint.
C     G      is a REAL vector of length N.
C
C   RESULT
C     A boolean value is returned that indicates if the equation was
C     solved properly.
C
C   SIDE EFFECTS
C     The solution is stored in P, and the gradient in that point is
C     stored in G.
C
C   TODO
C     * Maybe some boundary conditions should be implemented.
C
C   SEE ALSO
C     Any book in covering numerical methods, e.g. Numerical Recepies.
C
C     ***
      LOGICAL FUNCTION NEWRAPDIR( FUNC, VAL, NRDIM, D, ACC, P, G )
      IMPLICIT NONE
      REAL FUNC, D(*), G(*), P(*),  VAL, ACC
      INTEGER NRDIM
      EXTERNAL FUNC
      
      INTEGER N, J, JMAX
      REAL NORM, F, DF, DX, EPS
      PARAMETER ( JMAX = 20, EPS = 1.0E-10 )
      
      NEWRAPDIR = .TRUE.

      P(1) = P(1) + 0.074*D(1)
      P(2) = P(2) + 0.43*D(2)
      P(3) = P(3) + 0.17*D(3)
      

C
C        Solve the equation VAL = FUNC(P) in the direction D.
C
      DX = 1000.0
      J = 0
      write (*,*) 'start'
      DO WHILE ( J.LT.JMAX .AND. ABS(DX).GT.ACC .AND. NEWRAPDIR )
         J = J + 1

         F = FUNC( P, G, NRDIM ) - VAL
         write (*,*) (p(n), n = 1, nrdim), (d(n), n = 1, nrdim),
     $        (g(n), n = 1, nrdim), f, df, dx, val
C
C           Calculate the derivative along the direction D in the point
C           P.
C
         DF = 0.0
         DO N = 1, NRDIM
            DF = DF + D(N)*G(N)
         ENDDO
         
C
C           Calculate the displacement.
C
         IF ( ABS(DF).GT.EPS ) DX = F/DF

c         write (*,*) f, df, dx
c         write (*,*) (p(n), N = 1, nrdim)
C
C           Update the point.
C
         DO N = 1, NRDIM
            P(N) = P(N) - DX*D(N)
         ENDDO
      ENDDO
      write (*,*) 'stop'

C
C        Check if the equation was solved correctly.
C
      IF ( ABS(DX).GT.ACC ) NEWRAPDIR = .FALSE.

      RETURN
      END
C     *** End of FUNCTION NEWRAPDIR() ***


C     ****i* solvedir/HALFINTDIR *
C
C   NAME
C     HALFINTDIR() -- solves a one-dimensional equation along a line
C                     on a multi-dimensional function.
C
C   DESCRIPTION
C     Solves a one-dimensional equation along a specified direction on
C     a multidimensional function using the half-int.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C
C   CREATION DATE
C     2000-10-14
C
C   USAGE
C     HALFINTDIR( FUNC, VAL, N, D, ACC, P, G )
C
C   INPUTS
C     FUNC   is a REAL function that takes three arguments: FUNC( P, G, N )
C            where P and G are REAL vectors that are at least N elements
C            long (the dimension of the function/number of variables).
C            FUNC should calculate the value of the function in the point P
C            and should also calculate the gradiant in that point and return
C            it in G as a side effect.
C     VAL    is the value for which the equation will be solved.
C     N      is dimension of the function.
C     D      is the direction along FUNC where the equation will be solved.
C     ACC    is the accuracy that should be used.
C     P      is a REAL vector of length N, containing the startpoint.
C     G      is a REAL vector of length N.
C
C   RESULT
C     A boolean value is returned that indicates if the equation was
C     solved properly.
C
C   SIDE EFFECTS
C     The solution is stored in P, and the gradient in that point is
C     stored in G.
C
C   TODO
C     * Maybe some boundary conditions should be implemented.
C
C   SEE ALSO
C     Any book in covering numerical methods, e.g. Numerical Recepies.
C
C     ***
      LOGICAL FUNCTION HALFINTDIR( FUNC, VAL, NRDIM, D, ACC, P, G )
      IMPLICIT NONE
      REAL FUNC, D(*), G(*), P(*),  VAL, ACC
      INTEGER NRDIM
      EXTERNAL FUNC
      
      REAL MOVE, INIMOVE, STEP, F1, F2, DIR, EPS, TMP
      PARAMETER ( INIMOVE = 0.01, MOVE = 4.0, EPS = 1.0E-10 )
      INTEGER N, I, MAXINT, MAXDIM
      PARAMETER ( MAXINT = 100, MAXDIM = 20 )
      REAL TMPVEC(MAXDIM), TMPVEC2(MAXDIM), DIST
      
      

      DIR = 1.0
      HALFINTDIR = .TRUE.
C
C        Determine how far from the initial point the routine should move
C        for its first calculation.
C
      STEP = .0
      DO N = 1, NRDIM
         STEP = STEP + ABS( P(N)*D(N) )
      ENDDO
      IF ( STEP .LT. EPS ) STEP = 1.0 ! If P is the origo.
      STEP = INIMOVE*SQRT( STEP )

C
C        Calculate the function value in the initial point.
C
      F1 = FUNC( P, TMPVEC )

C
C        Calculate a second point, displaced STEP from the initial one.
C
      DO N = 1, NRDIM
         G(N) = P(N) + DIR*STEP*D(N)
      ENDDO
      F2 = FUNC( G, TMPVEC )
C
C        Determine in what direction the solution lies.
C
      IF ( (F1 - VAL)*(F1 - F2) .LT. 0.0 ) THEN
         DIR = -1.0
         DO N = 1, NRDIM
            G(N) = P(N) + DIR*STEP*D(N)
         ENDDO
         F2 = FUNC( G, TMPVEC )
      ENDIF

C
C        Now, when the direction from the initial point has been
C        determined it remains to find the boundings for the solution.
C
      I = 0
      DO WHILE ( ( (F1 - VAL)*(F2 - VAL) .GT. 0.0 )
     $     .AND. ( I .LT. 10*MAXINT ) )
         I = I + 1
         
         F1 = F2
         STEP = MOVE*STEP
         DO N = 1, NRDIM
            TMP = G(N)
            G(N) = P(N) + DIR*STEP*D(N)
            P(N) = TMP
         ENDDO
         F2 = FUNC( G, TMPVEC )
      ENDDO
      
C
C        Check if the maximum number of iterations has been
C        exceeded.
C
      IF ( I .LE. 10*MAXINT ) THEN
C
C           Iterate until the solution is found.
C
         I = 0
         DIST = 2*ACC
         DO WHILE ( DIST.GT.ACC .AND. I.LT.MAXINT )
            I = I + 1
C
C              Calculate the function value in the middle of
C              the interval.
C
            TMP = DIST/2.0
            DO N = 1, NRDIM
               TMPVEC2(N) = P(N) + DIR*TMP*D(N)
            ENDDO
            TMP = FUNC( TMPVEC2, TMPVEC )

            IF ( (F1 - VAL)*(TMP - VAL) .LT. 0 ) THEN
               F2 = TMP
               G(N) = TMPVEC(2)
            ELSE
               F1 = TMP
               G(N) = TMPVEC(2)
            ENDIF

C
C              Calculate the length of the interval.
C
            DIST = 0.0
            DO N = 1, NRDIM
               DIST = DIST + (G(N) - P(N))**2
            ENDDO
            DIST = SQRT( DIST )
         
         ENDDO

         IF ( I .GE. MAXINT ) HALFINTDIR = .FALSE.
      ELSE
         HALFINTDIR = .FALSE.
      ENDIF

      RETURN
      END
C     *** End of FUNCTION HALFINTDIR() ***

C     END OF FILE
