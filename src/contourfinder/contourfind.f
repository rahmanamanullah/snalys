C     ****h* snalys/contourfinder [0.9a] *
C   NAME
C     contourfinder -- finds a contour around a minima for a given
C                      function.
C
C   DESCRIPTION
C     The package is used to find a specific contour level around a
C     minimum for a given function. It can not be used to find all
C     contours in a given range in its present state, so occasional
C     islands for example will not be found.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-08-23
C
C   USES
C     solvedir
C
C     ***


C     ****f* contourfinder/CONFIND *
C
C   NAME
C     CONFIND() -- finds a contour level around a specific minimum for a
C                  given function.
C
C   DESCRIPTION
C     The points (the contour) where a n-dimensional function, where n>1,
C     equals a specific value is calculated by CONFIND(). The function
C     will only find the closest contour to the start point.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-08-23
C
C   USAGE
C     CONFIND( FUNC, MINP, CON, NET, NETNR, N, ACC, METHOD )
C
C   INPUTS
C     FUNC   - is a REAL function that takes three arguments: FUNC( P, G, N )
C              where P and G are REAL vectors that are at least N elements
C              long (the dimension of the function/number of variables).
C              FUNC should calculate the value of the function in the point P
C              and should also return the gradient in that point in G.
C     MINP   - the start point for the minimization.
C     CON    - is the REAL value that the function should have on the contour,
C              that is the value for which the equation CON = FUNC( P, G, N) 
C              should be solved. The value CON must of course be greater
C              than the value of the function in the minimum point. The
C              contour corresponding to CON must be closed.
C     NET    - is a REAL vector that contains NETNR elements.
C     NETNR  - is the length of the vector NET.
C     N      - is the dimension of the space, i.e. the number of variables
C              in the function FUNC. This value must be less or equal to
C              MAXDIM().
C     ACC    - is the accuracy that should be used when the contour points
C              are found. In most cases 1E-2 or 1E-3 is enough.
C     METHOD - is an integer that defines the contour finding method that
C              should be used.
C
C   RESULT
C     A boolean value, indicating if the contour was found, is returned.
C
C   SIDE EFFECTS
C     NET contains the contour on return, that is all the points for
C     which the contour equation has been solved. The maximum number
C     of points that can be stored in NET is NET/N, and the points are
C     stored as [x1, x2, x3, ... , xN, y1, y2 , y3, ... , yN, ... ]
C     where x corresponds to the first variable in P and y the second.
C
C   SEE ALSO
C     SCONFIND()
C
C     ***
      LOGICAL FUNCTION CONFIND( FUNC, MINP, CON, NET, NETNR, NRDIM, ACC,
     $     METHOD )
      IMPLICIT NONE
      REAL FUNC
      EXTERNAL FUNC
      REAL NET(*), MINP(*), CON, ACC
      INTEGER NRDIM, NETNR, METHOD

      LOGICAL DUMMYFUNC, CONIFACE
      EXTERNAL DUMMYFUNC

      CONFIND = CONIFACE( FUNC, MINP, CON, NET, NETNR, NRDIM, ACC, 
     $     METHOD, DUMMYFUNC )

      RETURN
      END
C     *** End of FUNCTION CONFIND() ***

C     ****f* contourfinder/SCONFIND *
C
C   NAME
C     SCONFIND() -- finds a contour level around a specific minimum for a
C                   given function.
C
C   DESCRIPTION
C     The (n-1)-dimensional hypersurface where a n-dimensional (n > 1) 
C     function equals a specific value is calculated by SCONFIND(). The
C     function will find the contour closest to the start point.
C
C     The difference between this function and CONFIND() is that this
C     function also takes a function that will be applied to each subnet
C     as an argument.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-08-27
C
C   USAGE
C     SCONFIND( FUNC, MINP, CON, NET, NETNR, N, ACC, METHOD, SUBFUN )
C
C   INPUTS
C     FUNC   - is a REAL function that takes three arguments: FUNC( P, G, N )
C              where P and G are REAL vectors that are at least N elements
C              long (the dimension of the function/number of variables).
C              FUNC should calculate the value of the function in the point P
C              and should also return the gradient in that point in G.
C     MINP   - the start point for the minimization.
C     CON    - is the REAL value that the function should have on the contour,
C              that is the value for which the equation CON = FUNC( P, G, N) 
C              should be solved. The value CON must of course be greater
C              than the value of the function in the minimum point. The
C              contour corresponding to CON must be closed.
C     NET    - is a REAL vector that contains NETNR elements.
C     NETNR  - is the length of the vector NET.
C     N      - is the dimension of the space, i.e. the number of variables
C              in the function FUNC. This value must be less or equal to
C              MAXDIM().
C     ACC    - is the accuracy that should be used when the contour points
C              are found. In most cases 1E-2 or 1E-3 is enough.
C     METHOD - is an integer that defines the contour finding method that
C              should be used.
C     SUBFUN - a LOGICAL function that takes the arguments
C              SUBFUN(P,NET,NRDIM), where P is an INTEGER vector specifying
C              the INDICES in NET where the first coordinates of NRDIM points
C              defining a subnet are stored. Since NET only contains the
C              contour points on return and all information concerning the
C              subnets have been lost this function can be used and will be
C              applied only on the smallest subnets.
C
C   RESULT
C     A boolean value, indicating if the contour was found, is returned.
C
C   SIDE EFFECTS
C     NET contains the contour on return, that is all the points for
C     which the contour equation has been solved. The maximum number
C     of points that can be stored in NET is NET/N, and the points are
C     stored as [x1, x2, x3, ... , xN, y1, y2 , y3, ... , yN, ... ]
C     where x corresponds to the first variable in P and y the second.
C
C   SEE ALSO
C     CONFIND()
C
C     ***
      LOGICAL FUNCTION SCONFIND( FUNC, MINP, CON, NET, NETNR, NRDIM,
     $     ACC, METHOD, SUBFUNC )
      IMPLICIT NONE
      REAL FUNC
      EXTERNAL FUNC
      LOGICAL SUBFUNC
      EXTERNAL SUBFUNC
      REAL NET(*), MINP(*), CON, ACC
      INTEGER NRDIM, NETNR, METHOD

      LOGICAL CONIFACE
      real g(10), tmp
      integer n

      SCONFIND = CONIFACE( FUNC, MINP, CON, NET, NETNR, NRDIM, ACC,
     $     METHOD, SUBFUNC)

      RETURN
      END
C     *** End of FUNCTION SCONFIND() ***

C     ****f* contourfinder/CONIFACE *
C
C   NAME
C     CONIFACE() -- finds a contour level around a specific minimum for a
C                  given function.
C
C   DESCRIPTION
C     The points (the contour) where a n-dimensional function, where n>1,
C     equals a specific value is calculated by CONIFACE(). The contour level
C     closest to the startpoint will only be found.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-08-27
C
C   USAGE
C     CONIFACE( FUNC, MINP, CON, NET, NETNR, N, ACC, METHOD, SUBFUN )
C
C   INPUTS
C     FUNC   - is a REAL function that takes three arguments: FUNC( P, G, N )
C              where P and G are REAL vectors that are at least N elements
C              long (the dimension of the function/number of variables).
C              FUNC should calculate the value of the function in the point P
C              and should also return the gradient in that point in G.
C     MINP   - a start point for the minimization.
C     CON    - is the REAL value that the function should have on the contour,
C              that is the value for which the equation CON = FUNC( P, G, N) 
C              should be solved. The value CON must of course be greater
C              than the value of the function in the minimum point. The
C              contour corresponding to CON must be closed.
C     NET    - is a REAL vector that contains NETNR elements.
C     NETNR  - is the length of the vector NET.
C     N      - is the dimension of the space, i.e. the number of variables
C              in the function FUNC. This value must be less or equal to
C              MAXDIM().
C     ACC    - is the accuracy that should be used when the contour points
C              are found. In most cases 1E-2 or 1E-3 is enough.
C     METHOD - is an integer that defines the contour finding method that
C              should be used.
C     SUBFUN - a function with the arguments SUBFUN( P, NET, NRDIM ), where
C              P is an INTEGER vector specifing the first coordinates of
C              the NRDIM points that defines a subnet. This function will
C              be called for the smallest subnets.
C
C   RESULT
C     A boolean value, indicating if the contour was found, is returned.
C
C   SIDE EFFECTS
C     NET contains the contour on return, that is all the points for
C     which the contour equation has been solved. The maximum number
C     of points that can be stored in NET is NET/N, and the points are
C     stored as [x1, x2, x3, ... , xN, y1, y2 , y3, ... , yN, ... ]
C     where x corresponds to the first variable in P and y the second.
C
C     ***
      LOGICAL FUNCTION CONIFACE( FUNC, MINP, CON, NET, NETNR, NRDIM,
     $     ACC, METHOD, SUBFUNC )
      IMPLICIT NONE
      REAL FUNC
      EXTERNAL FUNC
      LOGICAL SUBFUNC
      EXTERNAL SUBFUNC
      REAL MINP(*), NET(*), CON, ACC
      INTEGER NRDIM, NETNR, METHOD, N, I

      INCLUDE 'confind.inc'
      REAL G(MAXDIM), VARMAT(MAXDIM**2), TMP
      LOGICAL RUBBERNET

C
C        Calculate the gradient in the start point.
C
      TMP = FUNC( MINP, G )
C
C        Set the initial step lengths.
C
      DO N = 1, NRDIM
C
C           Set the step length to 1% of the point value if
C           grad = 0.
C
         IF ( G(N).EQ.0.0 ) THEN
            VARMAT((N-1)*NRDIM+N) = 0.01*MINP(N)
C
C           This may look a bit funny, but it is done to keep the
C           steplength 0.01 after the davidon routine mess up v.
C
         ELSE
            VARMAT((N-1)*NRDIM+N) = ABS(0.01/G(N))
         ENDIF
      ENDDO

C
C        Minimize the function.
C
      CALL DAVIDO (-NRDIM, MINP, VARMAT, FUNC, 1.0E-10)

C
C        Calculate the function value in the min point.
C
      TMP = FUNC( MINP, G )

      CONIFACE = RUBBERNET( FUNC, MINP, VARMAT, TMP + CON, NET, NETNR,
     $     NRDIM, ACC, SUBFUNC )

      RETURN
      END
C     *** End of FUNCTION CONIFACE() ***




C     ****i* contourfinder/PROJCON *
C
C   NAME
C     PROJCON() -- finds a contour around a specific minimum for a
C                  given function, by performing a grid, where the
C                  step length varies.
C
C   DESCRIPTION
C     Given a specific value the points (the contour) around a minimum
C     point where the function equals this value will be calculated.
C     This function uses brute force method and just scans all directions
C     from the minimum, and solve the contour equation along each line.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-08-24
C
C   USAGE
C     PROJCON( FUNC, MINP, VARMAT, CON, NET, NETNR, N, ACC, PLANE )
C
C   INPUTS
C     FUNC  - a REAL function that takes three arguments: FUNC( P, G, N )
C             where P and G are REAL vectors that are at least N elements
C             long (the dimension of the function/number of variables).
C             FUNC should calculate the value of the function in the point P
C             and should also return the gradiant in that point in G.
C     MINP  - a REAL vector not shorter than N, which contains the minima
C             point of the function.
C     VARMAT- the correlation matrix for the variables.
C     CON   - the REAL value that the function should have on the contour,
C             that is the the value for which the equation should be solved.
C             This must of course be greater than the value of the function
C             in the minimum point. The contour corresponding to this value
C             must be closed.
C     NET   - a REAL vector that at least contains NETNR elements.
C     NETNR - the minimum length of NET.
C     N     - dimension of the space, that is the number of variables
C             in the function FUNC.
C     ACC   - the accuracy that should be used when the contours are
C             found. In most cases 1E-2 or 1E-3 is enough.
C     PLANE - an INTEGER vector of length 2, where the plane of the
C             projection is defined.
C
C   SIDE EFFECTS
C     NET contains the contour on return, that is all the points for
C     which the contour equation has been solved. The maximum number
C     of points that can be stored in NET is NET/N
C
C   RESULT
C     The result of the function is a boolean value that represents if
C     the contour was found or not.
C
C   USES
C     RUBBERNET()
C
C     ***
      LOGICAL FUNCTION PROJCON( FUNC, MINP, VARMAT, CON, NET, NETNR,
     $     NRDIM, ACC, P, X, Y )
      IMPLICIT NONE
      REAL FUNC
      EXTERNAL FUNC
      INCLUDE 'confind.inc'
      REAL MINP(*), NET(*), CON, ACC, X(*), Y(*), VARMAT(*)
      INTEGER NRDIM, NETNR, P(2)
      
      INTEGER N, I, LENGTH, NRPOINTS, DIR, OFF, POS, PPOS, PD, NEW,
     $     IDX
      LOGICAL RUBBERNET, DUMMYFUNC
      EXTERNAL DUMMYFUNC
      REAL DX, DR, YDIFF, DY, YMAX, YMIN

      LENGTH = NETNR
      IF ( RUBBERNET(FUNC, MINP, VARMAT, CON, NET, LENGTH, NRDIM, ACC,
     $     DUMMYFUNC) ) THEN

C
C           Extract all the values corresponding to the index PLANE(1)
C
         NRPOINTS = LENGTH/NRDIM - NRDIM + 1
         DO N = 1, NRPOINTS
            X(N) = NET((N-1)*NRDIM + P(1))
         ENDDO
         
C
C           Sort the values. The vector Y will contain the indices of
C           the values in ascending order.
C
         CALL SORTRX( NRPOINTS, X, Y )

C
C           Devide the x-range into equally spaced subintervals of
C           length DX.
C
         OFF = INT(0.1*REAL(NRPOINTS))

         PPOS = 1
         POS = 2
         DIR = 1
         PD = P(1) - P(2)
         IDX = 1
         NEW = 1
         X(1) = Y(1)
         DY = ABS(NET((INT(Y(1))-1)*NRDIM+P(2)) -
     $        NET((INT(Y(2))-1)*NRDIM+P(2)))
         DO WHILE ( POS.GT.1 )
            NEW = POS
C
C              If we run over the index limits.
C
            IF ( POS + OFF .LE. NRPOINTS ) THEN
               PPOS = INT(Y(IDX))
               YMAX = NET((PPOS-1)*NRDIM + P(2)) + DY
               YMIN = NET((PPOS-1)*NRDIM + P(2)) - DY
C
C              Find the point with the smallest y-value if DIR.EQ.1
C              and the one with the largest if DIR.EQ.-1
C
               DO N = POS, POS + OFF
                  YDIFF = ABS( NET((PPOS-1)*NRDIM+P(2)) -
     $                 NET((INT(Y(N))-1)*NRDIM+P(2)) )
                  IF ( DIR .GT. 0 ) THEN
                     IF ( (YDIFF .LT. 5.0*DY)
     $                    .AND. ( NET((INT(Y(N))-1)*NRDIM+P(2)) .LT.
     $                    NET((INT(Y(NEW))-1)*NRDIM+P(2)) )) NEW = N
                  ELSE
                     IF ( (YDIFF .LT. 5.0*DY)
     $                    .AND. ( NET((INT(Y(N))-1)*NRDIM+P(2)) .GT.
     $                    NET((INT(Y(NEW))-1)*NRDIM+P(2)) )) NEW = N
                  ENDIF
               ENDDO
            ENDIF

C
C              Check if a next point was found.
C
            IF ( (ABS( NET((INT(X(IDX))-1)*NRDIM+P(2)) -
     $              NET((INT(Y(NEW))-1)*NRDIM+P(2)) ) .LT. 5.0*DY)
     $           .AND. (N .LT. NRPOINTS) ) THEN
               DY = ABS(NET((INT(X(IDX))-1)*NRDIM+P(1)) -
     $              NET((INT(Y(NEW))-1)*NRDIM+P(1)) )
               IDX = IDX + 1
               X(IDX) = REAL(NEW)
               IF ( POS .GT. 1 ) POS = POS + DIR*OFF
            ELSE
C
C              If no point was found we change the direction in the
C              search for the outer contour.
C
               DIR = -1*DIR
               IF ( DIR .LT. 0 ) POS = POS - OFF
            ENDIF

C
C              If we run over the index bound we have to continue
C              the search in the opposit direction.
C
C            IF ( N .GE. NRPOINTS ) DIR = -1
         ENDDO

C         DO N = 1, IDX
C            WRITE (*,*) (NET((X(N)-1)*NRDIM + I), I = 1, NRDIM)
C         ENDDO
      ELSE
         PROJCON = .FALSE.
      ENDIF
      RETURN
      END
C     *** End of FUNCTION PROJCON ***




C     ****i* contourfinder/RUBBERNET *
C
C   NAME
C     RUBBERNET() -- finds a contour around a specific minimum for a
C                    given function, by using a net, where the density
C                    is continously increased.
C
C   DESCRIPTION
C     Given a specific value the points (the contour) around a minimum
C     point where the function equals this value will be calculated.
C     This function uses the idea of a simple net that is initially
C     anchored to 2*N points where N is the dimension of the function.
C     these points will in the two-dimensional case build up a square,
C     and in the three-dimensional case a cube. The density of the net
C     is then increased by recursively find the contour between the
C     current points.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-08-20
C
C   USAGE
C     RUBBERNET( FUNC, MINP, VARMAT, CON, NET, NETNR, N, ACC, SFUNC )
C
C   INPUTS
C     FUNC  - is a REAL function that takes three arguments: FUNC( P, G, N )
C             where P and G are REAL vectors that are at least N elements
C             long (the dimension of the function/number of variables).
C             FUNC should calculate the value of the function in the point P
C             and should also return the gradiant in that point in G.
C     MINP  - is a REAL vector not shorter than N, which contains the minima
C             point of the function.
C     VARMAT- the correlation matrix for the variables.
C     CON   - is the REAL value that the function should have on the contour,
C             that is the the value for which the equation should be solved.
C             This must of course be greater than the value of the function
C             in the minimum point. The contour corresponding to this value
C             must be closed.
C     NET   - is a REAL vector that at least contains NETNR elements.
C     NETNR - is the minimum length of NET.
C     N     - is dimension of the space, that is the number of variables
C             in the function FUNC.
C     ACC   - is the accuracy that should be used when the contours are
C             found. In most cases 1E-2 or 1E-3 is enough.
C     SFUNC - a function which takes the arguments SFUNC( PIDX, NET, NRDIM ),
C             where PIDX is an INTEGER vector specifying the first coordinate
C             of the NRDIM points that defines a subnet. The function will
C             be applied to the smallest subnets.
C
C   SIDE EFFECTS
C     NET contains the contour on return, that is all the points for
C     which the contour equation has been solved. The maximum number
C     of points that can be stored in NET is NET/N
C
C   RESULT
C     The result of the function is a boolean value that represents if
C     the contour was found or not.
C
C   USES
C     NEXTNETPOINT()
C
C     ***
      LOGICAL FUNCTION RUBBERNET( FUNC, MINP, VARMAT, CON, NET, NETNR,
     $     NRDIM, ACC, SFUNC )
      IMPLICIT NONE
      REAL FUNC
      EXTERNAL FUNC
      LOGICAL SFUNC
      EXTERNAL SFUNC
      INCLUDE 'confind.inc'
      REAL MINP(*), NET(*), CON, ACC, VARMAT(*)
      INTEGER NRDIM, NETNR

      INTEGER N, I, J, PIDX(MAXDIM), AMIN, AMAX, NEXTNETPOINT,
     $     AAMIN, AAMAX, METHOD
      REAL P(MAXDIM), G(MAXDIM), UNVEC(MAXDIM), W(MAXDIM),
     $     V(MAXDIM,MAXDIM) 
      REAL PERCENT, EPS, TMP
      PARAMETER ( PERCENT = 0.01, EPS = 1E-6, METHOD = 2 )
      LOGICAL SOLVEDIR

      RUBBERNET = .TRUE.

C
C        Check that the NET array is long enough.
C
      IF ( NETNR.GE.(NRDIM*2**NRDIM + 2*NRDIM*NRDIM) ) THEN

         TMP = 0.0
         DO N = 1, NRDIM
            P(N) = MINP(N)
C
C              Replace the first column of the VARMAT with the errors.
C
            VARMAT(N) = SQRT( MAX( 0.0, VARMAT((N-1)*NRDIM+N) ) )
            TMP = TMP + VARMAT(N)**2
         ENDDO
         TMP = SQRT(TMP)

C        WRITE (*,*) 'err = ',  (VARMAT(N), N = 1, NRDIM)

         DO N = 2, NRDIM
            DO I = 1, NRDIM
               IF ( (N-1+I).GT.NRDIM ) THEN
                  VARMAT((N-1)*NRDIM+I) = VARMAT(N-1+I-NRDIM)
               ELSE
                  VARMAT((N-1)*NRDIM+I) = VARMAT(N-1+I)
               ENDIF
            ENDDO
            VARMAT((N-1)*NRDIM+N) = -1.0*VARMAT((N-1)*NRDIM+N)
C           WRITE (*,*)'vec',n,'=',(VARMAT((N-1)*NRDIM+I), I = 1, NRDIM)
         ENDDO
         CALL SVDCMP( VARMAT, NRDIM, NRDIM, MAXDIM, MAXDIM, W, V )
C        WRITE (*,*) 'w = ',(w(n), n = 1, nrdim)
         DO N = 1, NRDIM
C           WRITE (*,*)'vec',n,'=',(varmat((N-1)*NRDIM+I), i = 1, nrdim)
         ENDDO

C
C           Find the first 2*NRDIM points around the minimum along the
C           (+/-) directions.
C
         DO I = 1, 2
            DO N = 1, NRDIM
C
C                 Use the rows of the correlation matrix to as the
C                 initial directions to solve the contour equation.
C                 If you do not understand the meaning of this, think
C                 again!
C
               DO J = 1, NRDIM
                  UNVEC(J) = (-1.0)**I*VARMAT(NRDIM*(N-1)+J)
               ENDDO
C               DO J = N, NRDIM
C                  UNVEC(J) = (-1.0)**(I+1)*VARMAT(NRDIM*(N-1)+J)
C               ENDDO
C               CALL NORMVEC(UNVEC, 1, NRDIM)

C
C                 Displace the minimum slightly to make it possible
C                 to solve the contour equation.
C
               DO J = 1, NRDIM
C                  IF ( ABS(MINP(J)).LE.EPS ) THEN
C                     P(J) = PERCENT*UNVEC(J) + MINP(J)
C                  ELSE
C                     P(J) = PERCENT*MINP(J)*UNVEC(J) + MINP(J)
C                  ENDIF
                  P(J) = MINP(J) + TMP*UNVEC(J)
                  G(J) = MINP(J)
               ENDDO
C              write (*,*) 'p =',(p(J), j = 1, nrdim)
C              write (*,*) 'd =',(unvec(j), j = 1, nrdim)

               CALL NORMVEC(UNVEC, 1, NRDIM)
               IF ( SOLVEDIR(FUNC, CON, NRDIM, UNVEC, ACC, P, G, METHOD) 
     $              ) THEN

C
C                    Store the point.
C
                  DO J = 1, NRDIM
                     NET((2*N+I-3)*NRDIM+J) = P(J)
                  ENDDO
               ELSE
                  RUBBERNET = .FALSE.
               ENDIF

            ENDDO
         ENDDO

C
C           Check that the contour equation was solved for the
C           first 2*NRDIM points. If this was not successful then
C           it is not possible to proceed since there is no base net.
C
         IF ( RUBBERNET ) THEN
C
C              Calculate the offset in the NET array, where the next
C              point can be stored (the first 2*NRDIM points that each
C              need NRDIM indices have already been stored).
C
            AMIN = 2*NRDIM*NRDIM + 1

C
C              To increase the resolution of the net, the contour equation
C              has to be solved between the first points of the net. The
C              number of such points are 2^NRDIM. In the case of a
C              two-dimensional function the points will lie on the lines
C              that ties the first 2*NRDIM points together, and in case
C              of a three-dimensional function they will lie on
C              the surface that is defined by the three points surrounding
C              it.
C
            DO N = 1, 2**NRDIM
C
C                 Select a valid set of NRDIM points from which the next
C                 direction where the contour equation will be solved can
C                 be calculated. The index of the first coordinate of the
C                 point is used as reference.
C
               DO I = 1, NRDIM
                  PIDX(I) = 1 + (I-1)*2*NRDIM + NRDIM*INT(1 + 
     $                 0.5*SIN(2.0**(I-NRDIM)*3.141592
     $                 *(REAL(N)-0.5*(1.0 + 2.0**NRDIM)) ) )
               ENDDO

C
C                 Increase the density of the net by inserting a point
C                 in the middle of the NRDIM points just selected.
C                 NEXTNETPOINT returns the index of the last coordinate of
C                 last point stored in NET, and AMIN + 1 will be used
C                 as the start value in the net for next subnet.
C
               AAMIN = AMIN
               AMAX = INT(REAL(NETNR - AMIN + 1)/REAL(2**NRDIM - N + 1))
     $              + (AMIN - 1)
               AMIN = NEXTNETPOINT( FUNC, CON, ACC, NET, AAMIN, AMAX,
     $              NRDIM, PIDX, SFUNC ) + 1
            ENDDO
C
C              Return the last index where a point has been stored.
C
            NETNR = AMIN - 1
         ENDIF
      ELSE
         RUBBERNET = .FALSE.
      ENDIF

      RETURN
      END
C     *** End of FUNCTION RUBBERNET() ***


C     ****i* contourfinder/NEXTNETPOINT *
C
C   NAME
C     NEXTNETPOINT() -- the contour equation is solved for a point
C                       that lies between the N points given to
C                       the function.
C
C   DESCRIPTION
C     This is a recursive function that for a point that lies on the
C     line/surface/hypersurface defined by the N points that are
C     passed to the function, solves the contour equation along a line
C     that is orthogonal to the line/surface/hypersurface. The result is
C     stored in NET.
C
C     The NRDIM + 1 points are then passed to the function INCNETDEN(),
C     which calls NEXTNETPOINT() again to increase the density of the net.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-08-20
C
C   USAGE
C     NEXTNETPOINT( FUNC, CON, ACC, NET, AMIN, AMAX, N, ACC, PIDX, 
C                   SFUNC )
C
C   INPUTS
C     FUNC  - is a REAL function that takes three arguments: FUNC( P, G, N )
C             where P and G are REAL vectors that are at least N elements
C             long (the dimension of the function/number of variables).
C             FUNC should calculate the value of the function in the point P
C             and should also return the gradiant in that point in G.
C     CON   - is the REAL value that the function should have on the contour,
C             that is the the value for which the equation should be solved.
C             This must of course be greater than the value of the function
C             in the minimum point. The contour corresponding to this value
C             must be closed.
C     ACC   - is the accuracy that should be used when the contours are
C             found. In most cases 1E-2 or 1E-3 is enough.
C     NET   - is a REAL vector that is at least contains AMAX elements.
C     AMIN  - is an the index position in NET where the first point can
C             be stored. AMIN must be less than AMAX.
C     AMAX  - is the last index position in NET that can be used. AMAX must
C             be greater than AMIN.
C     N     - is dimension of the space, that is the number of variables
C             in the function FUNC.
C     PIDX  - is a vector that is not shorter than N, and contains the
C             indices for the first coordinate of the N points that spans
C             the line/surface/hypersurface that was mentioned in the
C             description.
C     SFUNC - a function with the arguments SUBFUN( P, NET, NRDIM ), where
C             P is an INTEGER vector specifing the first coordinates of
C             the NRDIM points that defines a subnet. This function will
C             be called for the smallest subnets.
C
C   SIDE EFFECTS
C     The calculated contour points are stored in NET.
C
C   RESULT
C     The last index used in the NET vector is returned.
C
C   NOTES
C     This function is used in a recursievly together with INCNETDEN().
C
C   SEE ALSO
C     INCNETDEN()
C
C   TODO
C     * Impliment a better way to calculate the position of the new point.
C
C     ***
      INTEGER FUNCTION NEXTNETPOINT( FUNC, CON , ACC, NET, AMIN, AMAX,
     $     NRDIM, PIDX, SFUNC )
      IMPLICIT NONE
      REAL NET(*), FUNC, CON , ACC
      EXTERNAL FUNC
      LOGICAL SFUNC
      EXTERNAL SFUNC
      INTEGER AMIN, AMAX, NRDIM, PIDX(*)
      
      INCLUDE 'confind.inc'
      INTEGER N, I, J, INFO, INCNETDEN, AAMIN, AAMAX, METHOD
      REAL TMP, RVAL, A(MAXDIM,MAXDIM), W(MAXDIM), V(MAXDIM,MAXDIM), 
     $     P(MAXDIM), G(MAXDIM), D(MAXDIM), RAND, NORM, TOTNORM
      REAL EPS
      PARAMETER ( EPS = 1.0E-10, METHOD = 2 )
      LOGICAL SOLVEDIR, OKGO, INCFURTHER, SLASK

      DO I = 1, NRDIM
         P(I) = 0.0
      ENDDO

C
C        Randomize a point on the line/surface/hypersurface that is
C        defined by the points in the PIDX vector.
C
      OKGO = .FALSE.
      DO WHILE ( .NOT. OKGO )
         OKGO = .TRUE.
         DO I = 1, NRDIM-1
            G(I) = RAND()
         ENDDO
            
C
C           This constraint is introduced to make sure that the
C           new point will still lie inside the vector range. To
C           avoid that this constraint will affect the unitary
C           distribution.
C
         TMP = 1.0
         IF (NRDIM.GE.2) THEN
            DO I = 2, NRDIM-1
               TMP = TMP - G(I-1)
               IF ( G(I) .GT. TMP ) OKGO = .FALSE.
            ENDDO
         ENDIF
      ENDDO

      TOTNORM = 0.0
      DO N = 2, NRDIM
         TMP = RAND() - 0.5
         NORM = 0.0
         DO I = 1, NRDIM
C
C              Calculate the vectors that span the surface/line
C              in which the NRDIM points lie. Also the sum of the
C              vector lengths are calculated to check if it is
C              worth to increase the density of the net.
C
            A(N-1,I) = NET(PIDX(N) - 1 + I) - NET(PIDX(1) - 1 + I)
            NORM = NORM + A(N-1,I)**2
            
C
C              Use the randomized points to calculate the point.
C
C           P(I) = P(I) + G(N-1)*A(N-1,I)

C
C              Calculate the point by slightly distort the centre
C              of gravity point in the surface/hypersurface.
C

            IF ( NRDIM .GT.2 ) THEN
               P(I) = P(I) + 0.3333*A(N-1,I) + TMP*0.1*A(N-1,I)
            ELSE
               P(I) = P(I) + 0.5*A(N-1,I)
            ENDIF
         ENDDO
         TOTNORM = TOTNORM + SQRT(NORM)
      ENDDO
      DO N = 1, NRDIM
         P(N) = P(N) + NET(PIDX(1) - 1 + N)
      ENDDO

C
C        Check if it is possible and neccessary to increase the net
C        density any further.
C

      IF (  (AMIN + NRDIM).LE.AMAX .AND. INCFURTHER( TOTNORM, PIDX,
     $     NET, NRDIM, ACC, FUNC )  ) THEN
C
C           Perform a Singular Value Decomposition (SVD) on the
C           A-matrix, to find a vector orthogonal to the ones that
C           span the line/surface/hypersurface.
C
C           REFERENCE: Numerical Recepies, sec 2.6.
C

         CALL SVDCMP(A,NRDIM-1,NRDIM,MAXDIM,MAXDIM,W,V)

C
C           The last column in V that corresponds to w=0
C           is the solution of the homogenous equation A*X=0.
C
         DO N = 1, NRDIM
            IF ( ABS(W(N)).LT.EPS ) THEN
               DO I = 1, NRDIM
                  D(I) = V(I,N)
               ENDDO
            ENDIF
         ENDDO
         
C
C           Solve the equation along the calculated direction.
C
         IF ( SOLVEDIR( FUNC, CON, NRDIM, D, ACC, P, G, METHOD ) )
     $        THEN
C
C              Store the point.
C
            DO N = 1, NRDIM
               NET(AMIN - 1 + N) = P(N)
            ENDDO
               
C
C              Increase the net density.
C
            AAMIN = AMIN + NRDIM
            AAMAX = AMAX
            NEXTNETPOINT = INCNETDEN( FUNC, CON, NRDIM, ACC, NET,
     $           PIDX, AMIN, AAMIN, AAMAX, SFUNC )
            
         ELSE
            write (*,*) 'oops'
            write (*,*) (D(N), N = 1, nRDIM)
            write (*,*) (P(N), N = 1, NRDIM)
C
C              If the equation could not be solved along the line
C              we do not try to increase the density of the net, and
C              return AMIN-1 as the last index used in the NET.
C
            NEXTNETPOINT = AMIN - 1
            SLASK = SFUNC( PIDX, NET, NRDIM )
         ENDIF
      ELSE
C
C           No calculations were made and AMIN-1 is returned as the
C           last index used of the net.
C
         SLASK = SFUNC( PIDX, NET, NRDIM )
         NEXTNETPOINT = AMIN - 1
      ENDIF

      RETURN
      END
C     *** End of FUNCTION NEXTNETPOINT() ***

C     ****i* contourfinder/INCFURTHER *
C
C   NAME
C     INCFURTHER() -- determines if the density of the net should
C                     be increased.
C
C   DESCRIPTION
C     Determines if the density of the net should be increased any
C     further or if the resolution of the subnet is enough.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-08-24
C
C   USAGE
C     INCFURTHER( TOTLENGTH, PIDX, NET, NRDIM, ACC )
C
C   INPUTS
C     TOTLENGTH - the total length of the vectors that span the
C                 subnet.
C     PIDX      - an INTEGER vector of length NRDIM or greater where
C                 the elements are the indices in NET of the first
C                 coordinates of the points that defines the subnet.
C     NET       - the vector that contains all the points of the net.
C     NRDIM     - the dimension of the problem.
C     ACC       - the requested accuracy of the net.
C
C   RESULT
C     If the requested accuracy has been achieved false is returned.
C
C   TODO
C     Consider the difference between the gradients in the points
C     of the subnet, but this will increase the function calls to FUNC,
C     and may not make the routine to run faster if the function calls
C     require heavy calculations.
C
C     ***
      LOGICAL FUNCTION INCFURTHER( TOTLENGTH, PIDX, NET, NRDIM, ACC,
     $     FUNC )
      IMPLICIT NONE
      REAL TOTLENGTH, ACC, NET(*), FUNC
      EXTERNAL FUNC
      INTEGER PIDX(*), NRDIM

      INCLUDE 'confind.inc'
      INTEGER N, I, J
      REAL TMP, ANGLE, P(MAXDIM), G(MAXDIM), GRAD(MAXDIM**2),
     $     NORM(MAXDIM)

C
C        Calculate the gradients in the NRDIM points.
C
      DO N = 1, NRDIM
C
C           Extract the point.
C
         DO I = 1, NRDIM
            P(I) = NET(PIDX(N)-1+I)
         ENDDO
C
C           Calculate the gradient.
C
         TMP = FUNC( P, G )
C
C           Store the gradient, and calculate the length of it.
C
         TMP = 0.0
         DO I = 1, NRDIM
            GRAD((N-1)*NRDIM+I) = G(I)
            TMP = TMP + G(I)**2
         ENDDO
         NORM(N) = SQRT(TMP)
      ENDDO

C
C        Find the largest angle between the gradients.
C
      ANGLE = 0.0
      DO N = 1, NRDIM-1
         DO I = N+1, NRDIM
C
C              Calculate the dot-product between the gradients.
C
            TMP = 0.0
            DO J = 1, NRDIM
               TMP = TMP + GRAD((N-1)*NRDIM+J)*GRAD((I-1)*NRDIM+J)
            ENDDO
C
C              Calculate the angle.
C
            TMP = ACOS( TMP/(NORM(N)*NORM(I)) )
            IF ( TMP .GT. ANGLE ) ANGLE = TMP
         ENDDO
      ENDDO

C      TMP = REAL(10*NRDIM**2)*ACC
C      IF ( TOTLENGTH .GT. (REAL(10*NRDIM**2)*ACC) ) THEN
      IF ( ANGLE .GT. NRDIM**2*5.0*3.14159265/180.0 ) THEN
         INCFURTHER = .TRUE.
      ELSE
         INCFURTHER = .FALSE.
      ENDIF

      RETURN
      END
C     *** End of FUNCTION INCFURTHER ***

C     ****i* contourfinder/INCNETDEN *
C
C   NAME
C     INCNETDEN() -- the density of the net is increased.
C
C   DESCRIPTION
C     This is a recursive function that for a point that lies on the
C     line/surface/hypersurface passed by the program, calculates N
C     new subnets and uses NEXTNETPOINT() to solve the contour
C     equation for those nets.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-08-20
C
C   USAGE
C     INCNETDEN( FUNC, CON, N, ACC, NET, PIDX, MIDP, AMIN, AMAX, SFUN )
C
C   INPUTS
C     FUNC - is a REAL function that takes three arguments: FUNC( P, G, N )
C            where P and G are REAL vectors that are at least N elements
C            long (the dimension of the function/number of variables).
C            FUNC should calculate the value of the function in the point P
C            and should also return the gradiant in that point in G.
C     CON  - is the REAL value that the function should have on the contour,
C            that is the the value for which the equation should be solved.
C            This must of course be greater than the value of the function
C            in the minimum point. The contour corresponding to this value
C            must be closed.
C     N    - is dimension of the space, that is the number of variables
C            in the function FUNC.
C     ACC  - is the accuracy that should be used when the contours are
C            found. In most cases 1E-2 or 1E-3 is enough.
C     NET  - is a REAL vector that is at least contains AMAX elements.
C     PIDX - is a vector that is not shorter than N, and contains the
C            indices for the first coordinate of the N points.
C     MIDP - is the index for the first coordinate that lies between the
C            points in PIDX and on the line/surface/hypersurface that
C            is spanned by those points.
C     AMIN - is an the index position in NET where the first point can
C            be stored. AMIN must be less than AMAX.
C     AMAX - is the last index position in NET that can be used. AMAX must
C            be greater than AMIN.
C     SFUN - a function with the arguments SUBFUN( P, NET, NRDIM ), where
C            P is an INTEGER vector specifing the first coordinates of
C            the NRDIM points that defines a subnet. This function will
C            be called for the smallest subnets.
C
C
C   SIDE EFFECTS
C     The calculated contour points are stored in NET.
C
C   RESULT
C     The last index used in the NET vector is returned.
C
C   NOTES
C     This function is used in a recursively together with NEXTNETPOINT().
C
C   SEE ALSO
C     NEXTNETPOINT()
C
C   TODO
C     * Impliment a better way to calculate the position of the new point.
C     * Impliment a routine that calculate the length/area/volume of each
C       subnets and devides the points on that bases.
C
C     ***
      INTEGER FUNCTION INCNETDEN( FUNC, CON, NRDIM, ACC, NET, PIDX,
     $           MIDP, AMIN, AMAX, SFUN )
      IMPLICIT NONE
      REAL FUNC, CON, ACC, NET(*)
      EXTERNAL FUNC
      LOGICAL SFUN
      EXTERNAL SFUN
      INTEGER MIDP, AMIN, NRDIM, AMAX, PIDX(*)

      INCLUDE 'confind.inc'
      INTEGER N, I, J, PG(MAXDIM), NEXTNETPOINT, AAMAX, AAMIN, TMP

C
C        The mid-point should be part of each new subnet created
C        by the function.
C
      PG(1) = MIDP

C
C        Create NRDIM new subnets. The subnets are created by
C        constructing a circular que.
C
      TMP = AMIN
      DO N = 1, NRDIM
         DO I = 1, NRDIM - 1
            IF ( (I + (N - 1)).LE.NRDIM ) THEN
               PG(I+1) = PIDX(I + (N - 1))
            ELSE
               PG(I+1) = PIDX(I + (N - 1) - NRDIM)
            ENDIF
         ENDDO

         AAMIN = TMP
C
C           The last index used by the previous subnet is used as the AMIN
C           for the next one.
C
         AAMAX = INT(REAL((AMAX - AAMIN + 1))/REAL(NRDIM - N + 1)) + 
     $        (AAMIN - 1)
         TMP = NEXTNETPOINT( FUNC, CON, ACC, NET, AAMIN, AAMAX, 
     $     NRDIM, PG, SFUN ) + 1
      ENDDO

C
C        Return the last index used in NET.
C
      INCNETDEN = TMP - 1

      RETURN
      END
C     *** End of FUNCTION INCNETDEN() ***

C     ****i* contourfinder/DUMMYFUNC *
C
C   NAME
C     DUMMYFUNC() -- does nothing.
C
C   DESCRIPTION
C     This function is passed to RUBBERNET() by CONFIND(), as a
C     substitute for the function passed to SCONFIND().
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-08-27
C
C   USAGE
C     DUMMYFUNC( PIDX, NET, NRDIM )
C
C   INPUTS
C     PIDX  - an INTEGER vector.
C     NET   - a REAL vector.
C     NRDIM - an INTEGER.
C
C   RESULT
C     Always returns TRUE.
C
C   USED BY
C     CONFIND()
C
C     ***
      LOGICAL FUNCTION DUMMYFUNC( PIDX, NET, NRDIM )
      IMPLICIT NONE
      INTEGER PIDX(*), NRDIM
      REAL NET(*)
      DUMMYFUNC = .TRUE.
      RETURN
      END
C     *** End of FUNCTION DUMMYFUNC() ***

C     END OF FILE
