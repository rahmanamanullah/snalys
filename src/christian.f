C     ****h* snalys/random [1.0] *
C   NAME
C     random -- a package for generating random numbers.
C
C   DESCRIPTION
C     The package generates random numbers.
C
C   AUTHOR
C     Christian Walck, (walck@physto.se)
C
C
C
C   CREATION DATE
C     1989-12-21
C
C     ***


C     ****f* random/RNUNIF *
C
C   NAME
C     RNUNIF() -- Universal pseudorandom number generator.
C
C   DESCRIPTION
C     Universal pseudorandom number generator by George Marsaglia & Arif
C     Zaman.
C
C     Initialize by a call to RNINIT and make use of, if wanted, the
C     utilities to pack/unpack seeds or use a clock-start.
C
C   AUTHOR
C     Christian Walck, (walck@physto.se)
C
C   CREATION DATE
C     1989-12-21
C
C   USAGE
C     RNUNIF( DUMMY )
C
C   INPUTS
C     DUMMY - A dummy variable.
C
C   RESULT
C     A boolean value, indicating if the contour was found, is returned.
C
C   SEE ALSO
C     George Marsaglia & Arif Zaman, FSU-SCRI-87-50, Florida State
C     University, 1987 (later on published together with a third author,
C     Wai Wan Tsang, in Stat. & Prob. Lett. 9 (1990) 35-39).
C
C     RNPACK(), RNINIT(), RNCLCK()
C
C     ***
      FUNCTION RNUNIF ( IDUM )
      COMMON /RNCMMN/ LAG1, LAG2, U(97), C, CD, CM
C
C     Floating point arithmetics version
C
    1 UNI = U(LAG1) - U(LAG2)
      IF ( UNI .LT. 0.0 ) UNI = UNI + 1.0
      U(LAG1) = UNI
      LAG1 = LAG1 - 1
      IF ( LAG1 .EQ. 0 ) LAG1 = 97
      LAG2 = LAG2 - 1
      IF ( LAG2 .EQ. 0 ) LAG2 = 97
      C = C - CD
      IF ( C .LT. 0.0 ) C = C + CM
      UNI = UNI - C
      IF ( UNI .LT. 0.0 ) UNI = UNI + 1.0
C     Avoid returning zero
      IF ( UNI .EQ. 0.0 ) GO TO 1
      RNUNIF = UNI
      RETURN
      END
C
C     End of RNUNIF()
C



C     ****f* random/RNPACK *
C
C   NAME
C     RNPACK() -- Pack the four initial seeds I, J, K and L into one
C                 integer ISEED.
C
C   DESCRIPTION
C     Pack the four initial seeds I, J, K and L into one integer ISEED.
C
C   AUTHOR
C     Christian Walck, (walck@physto.se)
C
C   CREATION DATE
C     1989-12-21
C
C   USAGE
C     RNPACK ( I, J, K, L, ISEED )
C
C   INPUTS
C     I, J, K - integers in the range 2-177
C     L       - integer in the range 0-168
C     
C
C   OUTPUT
C     ISEED is set to a value between 0 and 921350143.
C
C   SEE ALSO
C     RNINIT(), RNCLCK()
C
C   HISTORY
C     1991-01-14 minor bug fixed.
C
C   ERRORS
C     Writes an error message if one of the input integers are out
C     of range.
C
C     ***
      SUBROUTINE RNPACK ( I, J, K, L, ISEED )
C
      IF ( I.LE.1 .OR. I.GE.178 ) GO TO 999
      IF ( J.LE.1 .OR. J.GE.178 ) GO TO 999
      IF ( K.LE.1 .OR. K.GE.178 ) GO TO 999
      IF ( L.LT.0 .OR. L.GE.169 ) GO TO 999
C
C     Pack seeds densely into ISEED using bases 176, 176, 176 and 169
C     for the four seeds, respectively.
C
      ISEED = 169 * ( 176 * ( 176*(I-2) + (J-2) ) + (K-2) ) + L
      RETURN
C
C     Abort run if illegal values are used in the call
C
  999 WRITE (6,1000) I, J, K, L
 1000 FORMAT (//' RNPACK: error - illegal start seeds',4I12//)
      STOP
      END
C
C     End of RNPACK().
C




C     ****f* random/RNINIT *
C
C   NAME
C     RNINIT() -- Initialize random number generator by
C                 Marsaglia & Zaman (& Tsang).
C
C   DESCRIPTION
C     The initial values for U(97) are generated by using a 3-lag
C     Fibonacci generator with 3 initial seeds and a congruential
C     generator with one initial seed.
C
C     - The values of the three initial seeds for the 3-lag Fibonacci
C       generator may be in the range 1 to 178. This generator has a
C       period of 32044 in about half of all cases. In about a quarter
C       of all cases each the period is 16022 and 8011. In a few cases
C       with combinations like 1, 1, 178 or 178, 1, 178 the period is
C       only 1, 2 or 4. The authors proposes to avoid the combination
C       1, 1, 1 but we choose to be somewhat more restrictive and use
C       only values from 2 to 177 for the three seeds.
C     - The initial seed for the congruential generator may be in the
C       range 0 to 168. The period of this generator is 169 regardless
C       of the initial seed value.
C
C     The four initial seeds are packed into one integer which can be
C     regarded as a four digit number with a variable base for each
C     digit. This in order to get a one-to-one correspondence between
C     the four seeds and the big integer without leaving any holes.
C     With I, J, K being the seeds of the Fibonacci generator and L the
C     seed of the congruential generator we pack them as:
C        ISEED = (I-2)*176*176*169 + (J-2)*176*169 + (K-2)*169 + L
C     ISEED can thus have the range 0 to 169*176**3-1=921350143 all
C     values giving a legal and unique I, J, K, L combination.
C
C   AUTHOR
C     Christian Walck, (walck@physto.se)
C
C   CREATION DATE
C     1989-12-21
C
C   USAGE
C     RNINIT ( ISEED )
C
C   INPUTS
C     ISEED is in the range 0 to 169*176**3-1=921350143.
C
C   SIDE EFFECTS
C     Initialized to use RNUNIF().
C
C   SEE ALSO
C     RNPACK(), RNCLCK() (to create ISEED)
C
C
C   ERRORS
C     Writes an error message if ISEED is out of range.
C
C     ***
      SUBROUTINE RNINIT ( ISEED )
      COMMON /RNCMMN/ LAG1, LAG2, U(97), C, CD, CM
C
      LAG1 = 97
      LAG2 = 33
C
C     Abort run if ISEED value is out of its valid range
C
      IF ( ISEED.LT.0 .OR. ISEED.GT.921350143 ) THEN
         WRITE (6,1001) ISEED
 1001    FORMAT (//' RNINIT: error - ISEED out of range, ISEED=',I11
     +          ,'), run aborted.'//)
         STOP
      END IF
C
C     Unpack seeds for 3-lag Fibonacci and congruential generator
C
      I = MOD(ISEED/169/176/176,176) + 2
      J = MOD(ISEED/169/176,176) + 2
      K = MOD(ISEED/169,176) + 2
      L = MOD(ISEED,169)
C
C     Floating point arithmetics version
C
      DO 200 II = 1, 97
         S = 0.0
         T = 0.5
         DO 100 JJ = 1, 24
            M = MOD(MOD(I*J,179)*K,179)
            I = J
            J = K
            K = M
            L = MOD(53*L+1,169)
            IF ( MOD(L*M,64) .GE. 32 ) S = S + T
            T = 0.5 * T
  100    CONTINUE
         U(II) = S
  200 CONTINUE
      C  =   362436.0 / 16777216.0
      CD =  7654321.0 / 16777216.0
      CM = 16777213.0 / 16777216.0
      RETURN
      END
C
C     End of RNINIT().
C



C     ****f* random/RNCLCK *
C
C   NAME
C     RNCLCK() -- Clock start of random number generator.
C
C   DESCRIPTION
C     Calculate number of seconds from 1 January 1989 modulus 921350144
C     which may be used as start seed for the sequence. This implies
C     that one avoids using the same start sequence in any simulation
C     unless two jobs are ran in parallel during the same second!
C     The code previously used the CERN library routine DATIMH to
C     obtain the system date and time. This, however, is a unnecessary
C     detour in achieving this. A new subroutine RNDATE has been
C     supplied which directly uses system functions for date and time
C     on most computers.
C
C     Working in several laboratories one may, as has been done for
C     DELPHI, construct a value from a lab-code and the number of
C     seconds past during the current year.
C
C   AUTHOR
C     Christian Walck, (walck@physto.se)
C
C   CREATION DATE
C     1989-12-21
C
C   USAGE
C     RNCLCK ( ISEED )
C
C   OUTPUT
C     ISEED is set to a value between 0 and 921350143.
C
C   SEE ALSO
C     RNINIT(), RNPACK()
C
C   HISTORY
C     1990-09-22 revised by Ch.Walck
C     1994-09-05 bug fix by Ch.Walck
C     2000-10-25 removed the initial message, R. Amanullah (rahman@physto.se)
C
C     ***
      SUBROUTINE RNCLCK ( ISEED )
      SAVE IMD
      DIMENSION IMD(12)
      DATA IMD / 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 /
C
C     Get system date and time
C
      CALL RNDATE ( JY, JM, JD, IH, IM, IS )
C
C     Calculate number of days in previous years from 1 January 1989
C
      IY = 1900 + JY
      IF ( JY .LT. 89 ) IY = 2000 + JY
      IDY = 0
      IF ( IY .GT. 1989 ) THEN
         DO 100 I = 1989, IY-1
            JDY = 365
C           Code valid until 2100 which is not a leap-year (2000 is! 940905)
            IF ( MOD(I,4).EQ.0 ) JDY = 366
            IDY = IDY + JDY
  100    CONTINUE
      END IF
C
C     Number of days from 1 January until 1st current month
C     (treat leap-years as well)
C
      LY = 0
C     Code valid until 2100 which is not a leap-year (2000 is! 940905)
      IF ( MOD(JY,4).EQ.0 ) LY = 1
      ISHIFT = IMD(JM)
      IF ( JM .GT. 2 ) ISHIFT = IMD(JM) + LY
C
C     Seconds since January 1st 1989 at 00:00:00 (=0)
C
      ISEED = IS + 60 * ( IM + 60 * ( IH + 24 * ( JD-1+ISHIFT+IDY ) ) )
C
C     Bring packed seed into range (this modulus will zero the ISEED at
C     18:35:44 on March 13 2018 and at 13:11:28 on May 24 2047, ... )
C
      ISEED = MOD(ISEED,921350144)
C

C     WRITE (6,1003) JY, JM, JD, IH, IM, IS, ISEED
 1003 FORMAT (/' RNCLCK: date and time is',I3.2,'-',I2.2,'-',
     $     I2.2,I3.2,':',I2.2,':',I2.2,', clock start at seed',I10/)
      RETURN
      END
C
C     End of RNCLCK()
C


C     ****i* time/RNDATE *
C
C   NAME
C     RNDATE() -- Routine to give current date and time.
C
C   DESCRIPTION
C     Routine to give current date and time. This routine replaces the
C     old DATIMH routine from the CERN library which was an unnecessary
C     detour in achieving this information. The routine is not fully
C     tested on all computers but hopefully works anyway. If any problem
C     is encountered the self-tag DATIMH activates the old code.
C
C   AUTHOR
C     Christian Walck, (walck@physto.se)
C
C   CREATION DATE
C     1990-09-22
C
C   USAGE
C     RNDATE ( JY, JM, JD, IH, IM, IS )
C
C   OUTPUT
C     JY - Year (modulus 100)        IH - Hours
C     JM - Month                     IM - minutes
C     JD - Date                      IS - seconds
C
C   HISTORY
C     1999-02-09 Replacing working code giving Y2k-warnings/cw
C     2010-02-11 Updated this using DATE_AND_TIME F95 standard
C
C     ***
      SUBROUTINE RNDATE ( JY, JM, JD, IH, IM, IS )
      INTEGER STIME, T(8), TIME
C      EXTERNAL LTIME, TIME

      CALL DATE_AND_TIME(VALUES=T)
      JY = MOD(T(1),100)
      JM = T(2)
      JD = T(3)
      IH = T(5)
      IM = T(6)
      IS = T(7)

C      STIME = TIME()
C      CALL LTIME ( STIME, T )
C      JY  = MOD(T(6),100)
C      JM  = T(5) + 1
C      JD  = T(4)
C      IH  = T(3)
C      IM  = T(2)
C      IS  = T(1)
C     icc = ic

      RETURN
      END
C
C     End of RNDATE()
C


C
C     End of file
C
