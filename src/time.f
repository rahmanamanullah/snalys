C     ****f* date/SEC2STR *
C
C   NAME
C     SEC2STR() -- converts seconds to a string.
C
C   DESCRIPTION
C     Converts seconds to a string of the type HHHH:MM:SS.XX
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-09-13
C
C   USAGE
C     CALL SEC2STR( SECONDS , TIMESTR )
C
C   INPUT
C     SECONDS - the number of seconds to convert (of type DOUBLE PRECISION)
C     TIMESTR - string of type CHARACTER*24
C
C   SIDE EFFECTS
C     On return TIMESTR contains the string.
C
C     ***
      SUBROUTINE SEC2STR( SECONDS , TIMESTR )
      IMPLICIT NONE
      DOUBLE PRECISION SECONDS
      CHARACTER*24 TIMESTR

      INTEGER MAX, NR
      PARAMETER ( MAX = 6 )
      INTEGER*4 H, M, S, SS, HOURS(MAX), N, I

      S = SECONDS
      H = S/3600
      M = INT((SECONDS/3600.0 - DBLE(H))*60.0)
      S = SECONDS - DBLE(H)*3600.0 - DBLE(M)*60.0
      SS = INT((SECONDS - DBLE(H)*3600.0 - DBLE(M)*60.0 - DBLE(S))*100)
      
      NR = 0
      DO WHILE ( H .GT. 0 .AND. NR .LE. MAX )
         NR = NR + 1
         HOURS(NR) = INT((REAL(H)/10.0 - REAL(H/10))*10.0)
         H = H/10
         WRITE (*,*) H , NR
      ENDDO

      N = 1
      IF ( NR .GT. 0 ) THEN
         IF ( NR .EQ. 1 ) THEN
            TIMESTR(1:1) = '0'
            N = 2
         ENDIF
         DO I = 1, NR
            TIMESTR(I:I) = CHAR(ICHAR('0') + HOURS(NR-I+1))
         ENDDO
         N = NR + 1
      ELSE
         TIMESTR(1:2) = '00'
         N = 3
      ENDIF

      TIMESTR(N:N) = ':'
      N = N + 1

      CALL DATEINT2STR( M , TIMESTR(N:N+1) )
      N = N + 2
      TIMESTR(N:N) = ':'
      N = N + 1
      CALL DATEINT2STR( S , TIMESTR(N:N+1) )
      N = N + 2
      TIMESTR(N:N) = '.'
      N = N + 1
      CALL DATEINT2STR( SS , TIMESTR(N:N+1) )
      N = N + 2

      DO I = N, 24
         TIMESTR(I:I) = ' '
      ENDDO
      END
C
C End of SEC2STR
C

C     ****f* date/DATEINT2STR *
C
C   NAME
C     DATEINT2STR -- converts an date INTEGER to a CHARACTER array.
C
C   DESCRIPTION
C     Converts a date INTEGER to a CHARACTER string for printing, so
C     that for example february is printed '02' instead of '2'.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-10-25
C
C   USAGE
C     CALL DATEINT2STR( ITEM , STR )
C
C   INPUTS
C     ITEM - the INTEGER date
C     STR - A CHARACTER*(2) string.
C
C   SIDE EFFECTS
C     STR contains the date on return.
C
C     ***
      SUBROUTINE DATEINT2STR(ITEM,STR)
      IMPLICIT NONE

      INTEGER ITEM
      CHARACTER*(2) STR

      IF ( ITEM.LT.10 ) THEN
         STR(1:1) = '0'
         STR(2:2) = CHAR(ITEM+ICHAR('0'))
      ELSE
         STR(1:1) = CHAR(INT(ITEM/10)+ICHAR('0'))
         STR(2:2) = CHAR(ITEM-INT(ITEM/10)*10 + ICHAR('0'))
      ENDIF
      END
C
C End of DATEINT2STR
C


C     ****f* date/UNIX2C *
C
C   NAME
C     UNIX2C -- converts Unix system time to date/time INTEGER array.
C
C   DESCRIPTION
C     Converts Unix system time to date/time INTEGER array.
C
C   AUTHOR
C     Clive Page, Leicester University, UK.
C
C   CREATION DATE
C     1995-05-02
C
C   USAGE
C     CALL UNIX2C( UTIME , IDATE )
C
C   INPUTS
C     UTIME - INTEGER that is the unix time, i.e. seconds since 1970.0
C     IDATE - INTEGER vector with 6 elements.
C
C   SIDE EFFECTS
C     IDATE will contain: 1=year, 2=month, 3=date, 4=hour, 5=minute,
C     6=secs.
C
C   NOTE
C     There are plenty of intrinsic functions that does this in almost
C     all unices today.
C
C     ***
      SUBROUTINE UNIX2C (UTIME, IDATE)
      IMPLICIT NONE
      INTEGER UTIME, IDATE(6)

      INTEGER MJDAY, NSECS
      REAL DAY
C
C       Note the MJD algorithm only works from years 1901 to 2099.
C
      MJDAY    = INT(UTIME/86400 + 40587)
      IDATE(1) = 1858 + INT( (MJDAY + 321.51) / 365.25)
      DAY      = AINT( MOD(MJDAY + 262.25, 365.25) ) + 0.5
      IDATE(2) = 1 + INT(MOD(DAY / 30.6 + 2.0, 12.0) )
      IDATE(3) = 1 + INT(MOD(DAY,30.6))
      NSECS    = MOD(UTIME, 86400)
      IDATE(6) = MOD(NSECS, 60)
      NSECS    = NSECS / 60
      IDATE(5) = MOD(NSECS, 60)
      IDATE(4) = NSECS / 60
      END
C
C End of UNIX2C
C



C
C EOF
C
