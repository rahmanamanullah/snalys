C
C     PRECALCS.F
C
C     Calculations that should be made for each supernova before
C     the negative log-likelihood function can be calculated.
C
C     Rahman Amanullah (rahman@physto.se), 2001-10-28
C


C     ****f* snalys/NONCOSDEP *
C
C   NAME
C     NONCOSDEP - Calculating values that are cosmology independent.
C
C   DESCRIPTION
C     Calculating values that are cosmology independent and have to be
C     calculated for each supernova before the likelihood function can
C     be obtained.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-10-28
C
C   USAGE
C     NONCOSDEP( NR,INTSIG,COSFP)
C
C   INPUTS
C     NR     - The number of supernovae that should be used for the fit.
C     INTSIG - The intrinsic magnitude error that is going to be used in
C              this fit.
C     COSFP  - REAL value giving the fraction of point masses that will
C              be used for this fit.
C             
C   RESULT
C     A LOGICAL value is returned to confirm that reading of the file
C     was successful.
C
C   SIDE EFFECTS
C     The COMMON vectors SNZS, SNMS and SNVS are calculated.
C     SNZS - the input red shifts
C     SNMS - the adjusted input magnitudes for each red shift in SNZS.
C     SNSS - the stretches for each red shift in SNZS
C     SNVS - the magnitudes variances for each red shift in SNZS.
C     SNIX - integer vector containing indices of the sorted SNZS
C            vector in ascending order.
C
C   USES
C     QSORT
C
C   TODO
C     * Make this function independent of the includes, since the includes
C       are now only used in init.f and this function.
C     * Consider if the gausslin.inc should be located where it is and
C       if the whole probability directory should still remain.
C
C   HISTORY
C     2003-03-31 Added SNSS, stretch.
C     2005-02-16 Added HOSTDM, MIXDM
C     2006-09-14 Added redshift errors (Jakob J)
C
C     ***
      SUBROUTINE NONCOSDEP(NR, INTSIG, COSFP)
      IMPLICIT NONE
      INTEGER NR
      REAL COSFP, INTSIG

      INCLUDE 'error.inc'
      INCLUDE 'operation.inc'
      INCLUDE 'supernova.inc'
      INCLUDE 'probability.inc'
      INCLUDE 'gausslin.inc'

C
C       These vectors must be passed to LOGMLFUNC in a COMMON block since
C       the minimization routines assumes that LOGMLFUNC only takes two
C       arguments.
C
      INCLUDE 'precalcs.inc'

      REAL SERROR
      EXTERNAL SERROR
      INTEGER N
C
C       The Counter must be passed in a COMMON block to GAUSSLRATIO
C
      INTEGER GLRNR
      COMMON /GAUSSLINRATIO/GLRNR

      REAL RTSAFE
      EXTERNAL GAUSSLRATIO, RTSAFE

C
C              Calculate vectors for different parameters that are not
C              cosmology dependent.
C
      DO N = 1, NR
C
C           Uncertainties needed for the multi images fit.
C
         IF ( OPER(DO_MULTI) ) THEN
            SN_DATA(ETDEL,N) = MULTI_ETDEL
            SN_DATA(EIMSEP,N) = MULTI_EIMSEP
            SN_DATA(ERATIO,N) = MULTI_EIMRAT
            SN_DATA(EZL,N) = MULTI_EZL
            SN_DATA(EVEL,N) = MULTI_EVEL
         ENDIF

         GLRNR = N
         SNZS(N) = SN_DATA(ZS,N)
         SNSS(N) = SN_DATA(STRETCH,N)
C     Redshift errors added by Jakob J 2006-09-14
         SNDZ(N)=SN_DATA(ZERR,N)


C
C                 Calculate the parameter vectors required for the GAUSSLIN
C                 distribution function.
C
         IF ( WHICH_PDF .EQ. P_GAUSSLIN ) THEN
C                  Check that the fraction of point masses is within
C                  tabulated within boundaries.
            IF ( COSFP.GE.GLFP(1).AND.COSFP.LE.GLFP(GLFPNR) ) THEN
C                       Check that the redshift is within the tabulated
C                       boundaries for the parameters.
               IF ( SNZS(N).GT.GLZ(1) .AND.
     $              SNZS(N).LT.GLZ(GLZNR) ) THEN
                  CALL SIMPBIINT(GLFP,GLZ,GLSIGMA,GLFPNR,GLZNR,
     $                 COSFP,SNZS(N),GAUSSLINPARS(GLPS,N))
                  CALL SIMPBIINT(GLFP,GLZ,GLM0,GLFPNR,GLZNR,COSFP,
     $                 SNZS(N),GAUSSLINPARS(GLPM0,N))
                  CALL SIMPBIINT(GLFP,GLZ,GLB,GLFPNR,GLZNR,COSFP,
     $                 SNZS(N),GAUSSLINPARS(GLPB,N))

C
C                       Calculate the approximation point.
C
                  GAUSSLINPARS(GLMLIM,N) = RTSAFE(GAUSSLRATIO,
     $                 -10., -1.E-5 , 1.E-5 )
               ELSE
                  WRITE (STDERR,'(3A)') 'warning: the red shift ',
     $                 'of one or more supernovae are not within ',
     $                 'tabulated boundaries.'
               ENDIF
            ELSE
               WRITE (STDERR,'(3A)') 'error: the fraction of ',
     $              'point masses are not within tabulated ',
     $              'boundaries.'
            ENDIF
         ENDIF
C
C                 Calculate the sigma for each supernova.
C
         SNVS(N) = INTSIG**2

C
C                 Add the measurement uncertainty
C
         IF ( OPER(USEMERR) ) THEN
            SNVS(N) = SNVS(N) + SERROR(N)**2
         ENDIF

C     
C                 Adjust the input magnitudes depending on the conditions
C                 for the fit.
C
C                 For the SNOC simulated supernovae the gravitational lens
C                 effects and an intrinsic spread in the magnitude can be
C                 added (see the general SNOC documentation for details).
C                 In other cases it will not mattar weather these are added
C                 or not, since they are set to zero in the init function.
C
         SNMS(N) = SN_DATA(BZMAG,N)
         IF ( OPER( ADD_LDM ) ) SNMS(N) = SNMS(N) + SN_DATA(LDM,N)
         IF ( OPER( ADD_GDM ) ) SNMS(N) = SNMS(N) + SN_DATA(GDM,N)
         IF ( OPER( ADD_HOSTDM ) ) THEN
            SNMS(N) = SNMS(N) + SN_DATA(HOSTDM,N)
            SNVS(N) = SNVS(N) + SN_DATA(DHOSTDM,N)**2
         ENDIF
         IF ( OPER( ADD_MIXDM ) ) SNMS(N) = SNMS(N) + SN_DATA(MDM,N)
         IF ( OPER( SUB_INTSIG ) ) SNMS(N) = SNMS(N) - SN_DATA(ISG,N)
      ENDDO

C
C           Sort the red shift in ascending order.
C
      CALL QSORT( NR, SNZS, SNIX )

      END
C
C     END OF NONCOSDEP
C


C     ****f* snalys/SERROR *
C
C   NAME
C     SERROR() -- The measurement error
C
C   DESCRIPTION
C     Returns the magnitude measurement error for a specified
C     supernova.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-09-19
C
C   INPUTS
C     I - the number of the supernova in the data array.
C
C   RESULT
C     The measurement error is returned.
C
C   USES
C     supernova.inc
C
C   HISTORY
C     - 2006-09-14 I removed the contribution from the redshift error.
C       This contribution to the error is cosmology dependent and 
C       will be calculated in mlfunc.f
C     ***
      REAL FUNCTION SERROR(I)
      IMPLICIT NONE
      INTEGER I

      INCLUDE 'supernova.inc'

      SERROR = SQRT( SN_DATA(MERR,I)**2)
C      SERROR = SQRT( SN_DATA(MERR,I)**2  + SN_DATA(ZERR,I)**2 )
      RETURN
      END
C
C End of SERROR()
C

      SUBROUTINE COEFBOUNDFUNC()
      IMPLICIT NONE

      INCLUDE 'cosmo.inc'

      REAL LOGMLFUNC, P(NR_PAR), G(NR_PAR)
      REAL XINT(2)
      DATA XINT(1) , XINT(2) / 0.0 , 3.0 /
C
C        The number of points to be used to fit the coefficients.
C
      INTEGER NRP, N
      PARAMETER ( NRP = 20 )

      
      

      RETURN
      END
C
C End of COEFBOUNDFUNC()
C


C     ****f* snalys/GAUSSLRATIO *
C
C   NAME
C     GAUSSLRATIO() -- The ratio for the GAUSSLIN function
C
C   DESCRIPTION
C     When the P_GAUSSLIN pdf has been choosen solving this equation
C     will give the magnitude value for which the approximation to
C     GAUSSLIN distribution function can be used. Before this function
C     can be run, the array GAUSSLINPARS has to be calculated.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-11-12
C
C   USAGE
C     CALL GAUSSLRATIO( M , F , D )
C
C   INPUTS
C     M - input REAL magnitude
C     F - REAL scalar
C     D - REAL scalar
C
C   SIDE EFFECTS
C     F - contains the function value in the point M on return.
C     D - contains the derivative of the function in the point M.
C
C   USES
C     supernova.inc
C
C   SEE ALSO
C     See Snalys note 1 for details.
C
C     ***
      SUBROUTINE GAUSSLRATIO( M , F , D )
      IMPLICIT NONE
      REAL M, F, D
      
      INCLUDE 'supernova.inc'
C
C       The percentage at which the approximation should be used.
C
      REAL XI, S, M0, B, MM, SIGMA, MDIFF
      PARAMETER (XI = 1.E-5)
      PARAMETER (S = 1.5)
      INTEGER GLRNR
      COMMON /GAUSSLINRATIO/GLRNR

      REAL ELIN, DELIN

      B = GAUSSLINPARS(GLPB,GLRNR)
      M0 = GAUSSLINPARS(GLPM0,GLRNR)
      MM = SN_DATA(BZMAG,GLRNR)
      SIGMA = GAUSSLINPARS(GLPS,GLRNR)

      MDIFF = M

      F = ALOG(XI*ELIN(MDIFF, SIGMA**2, B, M0)) +
     $     (MDIFF - M0)**2/(2*SIGMA**2)
C     F = ALOG(XI*B*ABS(MDIFF)) + S*(MDIFF)*ALOG(10.) + 
C    $     (MDIFF - M0)**2/(2*SIGMA**2)

C     D = XI*DLELIN(MDIFF, SIGMA**2, B, M0) + (MDIFF - M0)/SIGMA**2
      D = DELIN(MDIFF, SIGMA**2, B, M0)/ELIN(MDIFF, SIGMA**2, B, M0) +
     $     (MDIFF - M0)/SIGMA**2

      RETURN
      END
C
C End of GAUSSLRATIO
C
