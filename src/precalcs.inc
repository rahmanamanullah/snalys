C     ****v* snalys/SNVECS [1.0] *
C   NAME
C     SNVECS -- pre calculated vectors that are used for calculating the
C               likelihood function.
C
C   DESCRIPTION
C     Vectors that contains values for each supernova, such as red shifts
C     magnitude variances, luminosity distances, etc.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   USES
C     MAX_SN
C
C   CREATION DATE
C     2001-10-28
C
C   HISTORY
C   - 2006-09-14 Added SNDZ containing redshift errors
C     ***
      REAL SNZS,   ! The red shifts from the data file
     $     SNMS,   ! The magnitudes from the data file (adjusted)
     $	   SNSS,   ! The stretches form the data file
     $     SNVS,   ! The magnitude variances for each supernova
     $     SNM,    ! The calculated magnitudes
     $     SNDL,   ! The calculated luminosity distances
     $     SNDZ    ! Error in redshift from data file
      INTEGER SNIX ! The sorted indices of SNZS in ascending order
      COMMON /NONCOSVECS/SNZS(MAX_SN), SNMS(MAX_SN), SNVS(MAX_SN),
     $     SNM(MAX_SN), SNDL(MAX_SN), SNSS(MAX_SN), SNDZ(MAX_SN)
      COMMON /SNIX/SNIX(MAX_SN)
C
C EOF
C
