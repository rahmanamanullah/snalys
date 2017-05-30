C     ****h* snalys/date [1.0] *
C   NAME
C     date -- a package to get the current date.
C
C   DESCRIPTION
C     This package reads the system time and returns it in a
C     human readble format.
C
C   AUTHOR
C     Christian Walck, (walck@physto.se)
C
C   CREATION DATE
C     1990-03-22
C
C     ***


C     ****f* date/GETDATE *
C
C   NAME
C     GETDATE() -- Routine to give current date and time.
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
C     GETDATE ( JY, JM, JD, IH, IM, IS )
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
      SUBROUTINE GETDATE ( JY, JM, JD, IH, IM, IS )
      INTEGER STIME, T(8), TIME

      CALL DATE_AND_TIME(VALUES=T)
      JY = MOD(T(1),100)
      JM = T(2)
      JD = T(3)
      IH = T(5)
      IM = T(6)
      IS = T(7)

      RETURN
      END
C
C     End of GETDATE()
C

C
C     End of file
C
