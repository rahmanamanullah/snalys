C     ****f* snalys/RESFITS *
C
C   NAME
C     RESFITS -- save the results in a fits file.
C
C   DESCRIPTION
C     The results are saved in a fits file, where the ML-function
C     is saved as an array extension.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2003-09-23
C
C   USAGE
C     CALL RESFITS( STIME, F_DEN, GRI, GERR )
C
C   INPUTS
C     STIME     - The start time for the execution (in unix format)
C     F_DEN     - the grid array of type REAL
C     GRI       - the best fit from the GRIDSEARCH.
C     GERR      - the estimated errors from GRIDSEARCH.
C
C   TODO
C     - The result of the fit stored in the fits header together
C       with used cpu time and so on.
C     - Add so that the input for the fit is stored as a second
C       extention (ASCII table)
C     - Add all parameters for the simulation to header of the
C       ascii table
C
C   BUGS
C     If the parameter list contains more than nine parameters, this
C     function will probably crash.
C
C   HISTORY
C     2004-01-08 Added all parameters for the fit to the header of
C                the primary extention.
C     2004-10-27 Added the saving of the values as header keywords of
C                scriptM and alpha for the best fit.
C
C     ***
      SUBROUTINE RESFITS( STIME, F_DEN, GRI, GERR )
      IMPLICIT NONE
      INTEGER*4 STIME, TIME
      EXTERNAL TIME
      REAL F_DEN(*), GRI(*), GERR(2,*)

      INCLUDE 'io.inc'
      INCLUDE 'error.inc'
      INCLUDE 'cosmo.inc'
      INCLUDE 'supernova.inc'
      INCLUDE 'operation.inc'
      INCLUDE 'probability.inc'
      
      INTEGER STATUS, BLOCKSIZE, UNIT, BITPIX, J, I, NELEMENTS, 
     $     group, fpixel
      INTEGER NAXES(NR_PAR), NAXIS
      REAL VERSION, C_NR_GET, C_START_GET, C_END_GET
      LOGICAL SIMPLE, EXTEND

      STATUS = 0

      CALL FTVERS(VERSION)
      IF ( VERSION.LT.2.47 ) WRITE (STDERR,'(A)')
     $     'resfits: warning version of cfitsio should be >=2.47'

      CALL FTSDFILE( OUTFILE , STATUS )
      CALL FTGIOU( UNIT , STATUS )
C
C        Create the new empty FITS file.  The blocksize parameter is a
C        historical artifact and the value is ignored by FITSIO.
C
      BLOCKSIZE = 1
      CALL FTINIT( UNIT, OUTFILE, BLOCKSIZE, STATUS )

C
C        Initialize parameters about the FITS image.
C        BITPIX = 16 means that the image pixels will consist of 16-bit
C        integers.  The size of the image is given by the NAXES values. 
C        The EXTEND = TRUE parameter indicates that the FITS file
C        may contain extensions following the primary array.
C
      SIMPLE = .TRUE.
      BITPIX = -32
      NAXIS = NR_VAR
      NELEMENTS = 1
      DO J = 1, NR_VAR
         NAXES(J) = INT(C_NR_GET(J))
         NELEMENTS = NELEMENTS*NAXES(J)
      ENDDO
      EXTEND = .TRUE.

C
C        Write the required header keywords to the file
C
      CALL FTPHPR(UNIT,SIMPLE,BITPIX,NAXIS,NAXES,0,1,EXTEND,STATUS)
      DO J = 1, NR_VAR
         CALL FTPKYS(UNIT,'LABEL'//CHAR(ICHAR('0') + J),P_NAMES(J),
     $        'Label of parameter '//CHAR(ICHAR('0') + J),STATUS)
         CALL FTPKYE(UNIT,'START'//CHAR(ICHAR('0') + J),
     $        C_START_GET(J),-3,
     $        'Start value for '//CHAR(ICHAR('0') + J),STATUS)
         CALL FTPKYE(UNIT,'END'//CHAR(ICHAR('0') + J),
     $        C_END_GET(J),-3,
     $        'End value for '//CHAR(ICHAR('0') + J),STATUS)
      ENDDO

      CALL FTPKYJ(UNIT,'NRSN',SN_USED,'The number of SNe used',STATUS)

      CALL FTPKYE(UNIT,'INTDISP',GAUSS_SIGMA,-3,
     $     'Assumed intrinsic dispersion',STATUS)
      IF ( WHICH_PDF .EQ. P_GAUSS) THEN
         CALL FTPKYS(UNIT,'PDF','Gaussian',
     $        'Assumed magnitude distribution',STATUS)
      ELSEIF (WHICH_PDF .EQ. P_GAUSSLIN ) THEN
         CALL FTPKYS(UNIT,'PDF','Gauss w. lin',
     $        'Assumed magnitude distribution',STATUS)
      ENDIF

      CALL FTPKYL(UNIT,'DOGRID',OPER(DO_GRID),
     $     'Grid search minimisation',STATUS)
      CALL FTPKYL(UNIT,'DOGENE',OPER(DO_GENETIC),
     $     'Minimisation using Genetic method',STATUS)
      CALL FTPKYL(UNIT,'DOPOWE',OPER(DO_POWELL),
     $     'Minimisation using Powell method',STATUS)
      CALL FTPKYL(UNIT,'DODAVI',OPER(DO_DAVIDON),
     $     'Minimisation using Davidon method',STATUS)
      CALL FTPKYL(UNIT,'ADDLDM',OPER(ADD_LDM),
     $     'Consider lensing effects',STATUS)
      CALL FTPKYL(UNIT,'ADDGDM',OPER(ADD_LDM),
     $     'Consider grey dust',STATUS)
      CALL FTPKYL(UNIT,'ADDHOSTDM',OPER(ADD_HOSTDM),
     $     'Consider the host galaxy extinction',STATUS)
      CALL FTPKYL(UNIT,'ADDMIXDM',OPER(ADD_MIXDM),
     $     'Consider photon-axion oscillations',STATUS)
      CALL FTPKYL(UNIT,'SUBINT',OPER(SUB_INTSIG),
     $     'Subtract intrinsic dispersion',STATUS)
      CALL FTPKYL(UNIT,'ZERROR',OPER(Z_ERROR),
     $     'Add redshift dependent error',STATUS)
      CALL FTPKYL(UNIT,'FLAT',OPER(FLAT),
     $     'Assume flat universe',STATUS)
      CALL FTPKYL(UNIT,'ZLIMIT',OPER(ZLIMIT),
     $     'Apply limits in the input redshifts',STATUS)
      IF ( OPER(ZLIMIT) ) THEN
         CALL FTPKYE(UNIT,'ZMIN',ZRANGE(ZLOW),-3,
     $        'Lower redshift limit',STATUS)
         CALL FTPKYE(UNIT,'ZMAX',ZRANGE(ZHIG),-3,
     $        'Upper redshift limit',STATUS)
      ENDIF
      CALL FTPKYL(UNIT,'QPOW',OPER(MATE_EQST),
     $     'Fitting alpha in the Peebles-Ratra potential',STATUS)

      CALL FTPKYL(UNIT,'DOMULTI',OPER(DO_MULTI),
     $     'Have the multiple method been used',STATUS)
      IF ( OPER(DO_MULTI) ) THEN
         CALL FTPKYE(UNIT,'JMAX',JBANDMAX,-3,
     $        'The faintness limit of discovery in J-band',STATUS)
         CALL FTPKYE(UNIT,'IMAX',IBANDMAX,-3,
     $        'The faintness limit of discovery in I-band',STATUS)
         CALL FTPKYS(UNIT,'LENSTYPE',LENSTYPE,
     $        'Lens model',STATUS)
         CALL FTPKYE(UNIT,'ETDELAY',MULTI_ETDEL,-3,
     $        'Uncertainty in days for the time delay',STATUS)
         CALL FTPKYE(UNIT,'EIMSEP',MULTI_EIMSEP,-3,
     $        'Uncertainty in arcsecs for the image sep.',STATUS)
         CALL FTPKYE(UNIT,'EIMRAT',MULTI_EIMRAT,-3,
     $        'Uncertainty between the ratio',STATUS)
         CALL FTPKYE(UNIT,'EZL',MULTI_EZL,-3,
     $        'Uncertainty in redshift of the lens',STATUS)
         CALL FTPKYE(UNIT,'EMODEL',MULTI_EMODEL,-3,
     $        'Uncertainty in the lens model',STATUS)
         CALL FTPKYE(UNIT,'EVEL',MULTI_EVEL,-3,
     $        'Uncertainty in ',STATUS)
      ENDIF

      CALL FTPKYL(UNIT,'XDIM',OPER(XDIM),
     $     'Have cosmology for X-dimension been fitted',STATUS)

      
      CALL FTPKYL(UNIT,'DOJACKE',OPER(DO_JACKE),
     $     'Have the Alam et al. parametrization been fitted',STATUS)

      CALL FTPKYL(UNIT,'DOLINDER',OPER(DO_LINDER),
     $     'Have the Linder parametrization been fitted',STATUS)

      CALL FTPKYL(UNIT,'DOHP',OPER(DO_HP),
     $     'Have the Hannestad parametrization been fitted',STATUS)

      IF ( OPER(FITMSC) ) THEN
         CALL FTPKYE(UNIT,'BESTSCM',BEST_SCRIPTM,-5,
     $        'The value of scriptM for the best fitted point',STATUS)
      ENDIF
      IF ( OPER(FITSALP) ) THEN
         CALL FTPKYE(UNIT,'BESTALPH',BEST_ALPHA,-5,
     $        'The value of alpha for the best fitted point',STATUS)
      ENDIF

C
C        Write the array to the FITS file. The last letter of the
C        subroutine name defines the datatype of the array argument; 
C        ('I' = I*2, 'J' = I*4 , 'E' = Real*4, 'D' = Real*8). The ND
C        array is treated as a single 1-D array with
C        NAXIS1 * NAXIS2 * ... * NAXISN total number of pixels.
C        GROUP is seldom used parameter that should almost always be
C        set = 1. FPIXEL is the first pixel.
C
      GROUP = 1
      FPIXEL = 1
      CALL FTPPRE(UNIT,GROUP,FPIXEL,NELEMENTS,F_DEN,STATUS)


C
C        The FITS file must always be closed before exiting the program.
C        Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
C
      CALL FTCLOS(UNIT, STATUS)
      CALL FTFIOU(UNIT, STATUS)

C
C        Check for any errors, and if so print out error
C        messages.
C
      IF ( STATUS .GT. 0 ) CALL FTSPRINTE(STATUS)


      RETURN
      END
C
C   END OF RESFITS
C


C     ****i* snalys/FTSDFILE *
C
C   NAME
C     FTSDFILE -- A simple little routine to delete a FITS file
C
C   DESCRIPTION
C     Delete the file if it already exists, so we can then recreate it.
C
C   AUTHOR
C     Stolen from cfitsio cookbook.f
C
C   CREATION DATE
C     2003-09-24
C
C   USAGES
C     CALL FTSDFILE( FILENAME , STATUS )
C
C   INPUTS
C     FILENAME - The name of the file
C     STATUS   - Integer CFITSIO status, must be zero for file to be
C                deleted
C
C     ***
      SUBROUTINE FTSDFILE( FILENAME , STATUS )
      IMPLICIT NONE

      INTEGER STATUS, UNIT, BLOCKSIZE
      CHARACTER*(*) FILENAME

C
C        Simply return if status is greater than zero
C
      IF ( STATUS .GT. 0 ) RETURN


C
C        Get an unused Logical Unit Number to use to
C        open the FITS file
C
      CALL FTGIOU( UNIT, STATUS )

C
C        Try to open the file, to see if it exists
C
      CALL FTOPEN( UNIT, FILENAME, 1, BLOCKSIZE, STATUS )

      IF ( STATUS .EQ. 0 ) THEN
C
C           File was opened;  so now delete it 
C
         CALL FTDELT(UNIT,STATUS)
      ELSE IF (STATUS .EQ. 103) THEN
C
C           File doesn't exist, so just reset status to zero
C           and clear errors
C
         STATUS = 0
         CALL FTCMSG
      ELSE
C
C           there was some other error opening the file; delete
C           the file anyway
C
         STATUS = 0
         CALL FTCMSG
         CALL FTDELT(UNIT,STATUS)
      END IF
C
C        Free the unit number for later reuse
C
      CALL FTFIOU(UNIT, STATUS)
      END
C
C END OF FTSDFILE
C



C     ****i* snalys/FTSPRINTE *
C
C   NAME
C     FTSPRINTE -- Prints error messages from CFITSIO
C
C   DESCRIPTION
C     This subroutine prints out the descriptive text corresponding to the
C     error status value and prints out the contents of the internal
C     error message stack generated by FITSIO whenever an error occurs.
C
C   AUTHOR
C     Stolen from cfitsio cookbook.f
C
C   CREATION DATE
C     2003-09-24
C
C   USAGES
C     CALL FTSPRINTE( STATUS )
C
C   INPUTS
C     STATUS   - Integer CFITSIO status
C
C     ***
      SUBROUTINE FTSPRINTE(STATUS)
      IMPLICIT NONE

      INTEGER STATUS
      CHARACTER ERRTEXT*30,ERRMESSAGE*80

      INCLUDE 'error.inc'

C
C        Check if status is OK (no error); if so, simply return
C
      IF (STATUS .LE. 0) RETURN

C
C        The FTGERR subroutine returns a descriptive 30-character text
C        string that corresponds to the integer error status number.
C        A complete list of all the error numbers can be found in the
C        back of the FITSIO User's Guide.
C
      CALL FTGERR(STATUS,ERRTEXT)
      WRITE(STDERR,'(AI2A)') 'FITSIO Error Status =',status,': ',errtext

C
C        FITSIO usually generates an internal stack of error messages
C        whenever an error occurs. These messages provide much more
C        information on the cause of the problem than can be provided by
C        the single integer error status value. The FTGMSG subroutine
C        retrieves the oldest message from the stack and shifts any
C        remaining messages on the stack down one position. FTGMSG is
C        called repeatedly until a blank message is returned, which
C        indicates that the stack is empty. Each error message may be up to
C        80 characters in length. Another subroutine, called FTCMSG, is
C        available to simply clear the whole error message stack in cases
C        where one is not interested in the contents.
C
      CALL FTGMSG(ERRMESSAGE)
      DO WHILE (ERRMESSAGE .NE. ' ')
          WRITE (STDERR,'(A)') ERRMESSAGE
          CALL FTGMSG(ERRMESSAGE)
      ENDDO
      END
C
C  END OF FTSPRINTE
C
