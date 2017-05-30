C     ****d* snalys/MAX_BIN [1.0] *
C   NAME
C     MAX_BIN -- the maximum number of bins
C
C   DESCRIPTION
C     The maximum number of red shifts bins that can be used.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-08-10
C
C     ***
      INTEGER MAX_BIN
      PARAMETER ( MAX_BIN = 1000 )


C     ****v* snalys/NR_BIN [1.0] *
C   NAME
C     NR_BIN -- the number of red shift bins
C
C   DESCRIPTION
C     The number of red shift bins that should be used if the
C     events are binned instead of used as individual events.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-08-10
C
C   NOTE
C     This value should be set in init.f
C
C     ***
      INTEGER NR_BIN
      COMMON /ZNRBIN/NR_BIN



C     ****v* snalys/ZINT [1.0] *
C   NAME
C     ZINT -- the interval for the binned red shifts
C
C   DESCRIPTION
C     A REAL vector containing the range for the ZBINs.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-08-10
C
C     ***
      REAL ZINT(2)
      COMMON /ZBINZINT/ZINT




C     ****v* snalys/ZBINS [1.0] *
C   NAME
C     ZBINS -- the red shift for each bin
C
C   DESCRIPTION
C     REAL vector that contains the red shift mid point for each bin.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-08-10
C
C     ***	
      REAL ZBINS(MAX_BIN)
      COMMON /ZBINZBINS/ZBINS



C     ****v* snalys/MEASMAG [1.0] *
C   NAME
C     MEASMAG -- the mean magnitude for each bin
C
C   DESCRIPTION
C     REAL vector that contains the mean magnitudes for each bin, that
C     is the magnitude each red shift in ZBINS.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-08-10
C
C     ***	
      REAL MEASMAG(MAX_BIN)
      COMMON /ZBINMEAS/MEASMAG
      


C     ****v* snalys/SIGMA_MAG [1.0] *
C   NAME
C     SIGMA_MAG -- the magnitude uncertainty for each bin
C
C   DESCRIPTION
C     The assumed uncertainty in magnitude fore each red shift
C     bin.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-08-10
C
C     ***	
      REAL SIGMA_MAG
      COMMON /ZBINSIGMA/SIGMA_MAG


C     ****v* snalys/CORR_COSMO [1.0] *
C   NAME
C     CORR_COSMO -- the cosmology to be used for the calculations
C
C   DESCRIPTION
C     A REAL vector containing the cosmology that should be used to
C     calculate the values in MEASMAG from the red shifts in ZBINS.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2000-08-10
C
C     ***
      REAL CORR_COSMO(10)
      COMMON /ZBINCORRCOS/CORR_COSMO


C
C End of zbin.inc
C

