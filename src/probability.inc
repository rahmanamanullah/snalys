C     ****v* snalys/WHICH_PDF [1.0] *
C   NAME
C     WHICH_PDF -- the probability density function
C
C   DESCRIPTION
C     This variables represents what probability density function
C     to use for the magnitudes in the analysis. Currently only
C     Gaussians are implemented.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-03-22
C
C   NOTE
C     The value of this variable should be set in init.f
C
C     ***
      INTEGER WHICH_PDF, P_GAUSS, P_GAUSSLIN, VERYTMP
      COMMON /WHICH_PDF/WHICH_PDF
      PARAMETER ( P_GAUSS = 1, P_GAUSSLIN = 2 )

C
C     The massprior storage structure. The first value is the mean
C     the second the sigma.
C
      INTEGER MPRNRPAR, MPRMEAN , MPRSTD
      PARAMETER ( MPRMEAN = 1 , MPRSTD = 2 , MPRNRPAR = 2 )
      REAL MPRPARS(MPRNRPAR)
      COMMON /MAPRIOR/MPRPARS


      INTEGER NRSIG
      COMMON /NRSIG/NRSIG

C     ****d* snalys/SIGZ [1.0] *
C   NAME
C     SIGZ -- values to keep the same prob. content for higher dimensions
C
C   DESCRIPTION
C     This array contains factors that the obtained sigmas should
C     be multiplied with to keep the same joint probability content
C     in higher dimensions for a Gaussian. In one dimension the
C     probability content within one sigma is 68.3%. This is also the
C     case for each parameter in a multi-dimensional fit, but not
C     for the joint probability. This values contain the 68.3%, 90%,
C     95% and 99% correction factors.
C
C   EXAMPLE
C     If a joint probabilty of 68.3% of the measurements is required
C     within the given error in four-dimensions the SIGZ(4,1) factor
C     will do the job. If instead 90% is required use SIGZ(4,2).
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-03-22
C
C     ***
      REAL SIG_Z(15,4)	
      DATA (SIG_Z(VERYTMP,1),VERYTMP = 1,15)
     $     /1.00,1.51,1.88,2.17,2.43,2.65,2.84,3.05,3.23,3.40,3.56,
     $     3.71,3.85,3.99,4.13/
      DATA (SIG_Z(VERYTMP,2),VERYTMP = 1,15)
     $     /1.64,2.15,2.50,2.79,3.04,3.26,3.47,3.66,3.83,4.00,4.16,
     $     4.31,4.45,4.59,4.72/
      DATA (SIG_Z(VERYTMP,3),VERYTMP = 1,15)
     $     /1.96,2.45,2.80,3.08,3.33,3.55,3.75,3.94,4.11,4.28,4.44,
     $     4.59,4.73,4.87,5.00/
      DATA (SIG_Z(VERYTMP,4),VERYTMP = 1,15)
     $     /2.58,3.03,3.37,3.64,3.88,4.10,4.30,4.48,4.66,4.82,4.97,
     $     5.12,5.26,5.40,5.53/
      
C
C End of probability.inc
C
