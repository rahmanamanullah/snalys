C
C     PARAMETERS
C


C     ****d* snalys/ABSMAG [1.0] *
C   NAME
C     ABSMAG -- the absolute magnitude of supernovae
C
C   DESCRIPTION
C     The absolute magnitude of a type Ia and type II supernovae
C     in B according to Dahlen & Fransson
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-11-04
C
C     ***
      REAL MAGIA,MAGIIN,MAGIBC,MAGIIP,MAGIIL,MAG87A
      REAL SMAGIA,SMAGIIN,SMAGIBC,SMAGIIP,SMAGIIL,SMAG87A
      PARAMETER (MAGIA  = -18.75) !H=0.85 WITH STRECTCH CORRECTION
      PARAMETER (SMAGIA  = 0.2)
      PARAMETER (MAGIIP = -17.41) !H=0.5 SIG=1.39
      PARAMETER (MAGIIN = -19.58) !H=0.5 SIG=0.6  
      PARAMETER (MAGIIL = -17.93) !H=0.5 SIG=0.51
      PARAMETER (MAGIBC = -17.99) !H=0.5 SIG=0.39
      PARAMETER (MAG87A = -14.98) !H=0.5 SIG=1.0
      PARAMETER (SMAGIIP = 1.39)
      PARAMETER (SMAGIIN = 0.6) 
      PARAMETER (SMAGIIL = 0.51)
      PARAMETER (SMAGIBC = 0.39)
      PARAMETER (SMAG87A = 1.0) 


C     ****d* snalys/MAX_SN [1.0] *
C   NAME
C     MAX_SN -- the maximum numer of events
C
C   DESCRIPTION
C     The maximum number of supernova events allowed by the program.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-03-22
C
C     ***
      INTEGER MAX_SN
      PARAMETER ( MAX_SN = 100000 )


C     ****d* snalys/ENR [1.0] *
C   NAME
C     ENR -- the number of data entries read by the program
C
C   DESCRIPTION
C     The number of entries read for each supernova, e.g. magnitude
C     red shift, contribution from lensing etc.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-03-22
C
C   HISTORY
C     2001-09-19 Changed the value to allow 6 entries.
C     2003-03-31 Added stretch.
C     2005-02-16 Added MDM,HOSTDM,DHOSTDM
C     2005-08-08 Added GDM
C
C     ***

      INTEGER ENR
      PARAMETER ( ENR = 29 )
C
C       The indices of these variables in sn_data.
C
      INTEGER ZS,BZMAG,ISG,LDM,MERR,ZERR
      PARAMETER (ZS=1,BZMAG=2,ISG=3,LDM=4,MERR=5,ZERR=6)

      INTEGER TDEL,IMSEP,RATIO,VEL,ZL
      PARAMETER (TDEL=7,IMSEP=8,RATIO=9,VEL=10,ZL=11)

      INTEGER ETDEL,EIMSEP,ERATIO,EVEL,EZL,LENSRAT,KCORBI,KCORBJ
      PARAMETER (ETDEL=12,EIMSEP=13,ERATIO=14,EVEL=15,EZL=16)
      PARAMETER (LENSRAT=17,KCORBI=18,KCORBJ=19 )

      INTEGER STRETCH
      PARAMETER (STRETCH=20)

      INTEGER MDM,HOSTDM,DHOSTDM,GDM
      PARAMETER (MDM=21,HOSTDM=22,DHOSTDM=23,GDM=24)

      INTEGER DSTRETCH,COLOR,DCOLOR,X1,DX1
      PARAMETER (DSTRETCH=25,COLOR=26,DCOLOR=27,X1=28,DX1=29)

C     ****d* snalys/ENAMES [1.0] *
C   NAME
C     ENAMES -- the names of the data entries
C
C   DESCRIPTION
C     CHARACTER array containing the names of the entries to be read
C     in the data file.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-03-22
C
C   HISTORY
C     2001-09-13 Added an extra entry name 'dzs' for the uncertainty
C                in redshift.
C     2003-03-31 Added stretch
C     2005-02-16 Added 'hostdm', 'dhostdm', 'mixdm'
C     2005-08-08 Added 'greydm'
C     2010-02-01 Added dstretch, color, dcolor
C     
C     ***
      CHARACTER*(9) ENAMES(ENR)
      INTEGER DINGEL2
      DATA (ENAMES(DINGEL2),DINGEL2=1,ENR)/'zs','BZmag','intsig',
     $     'linsdm','dBZmag','dzs',
     $     'midel','imsep','secrat','vlmi','zlmi',
     $     'emidel','eimsep','esecrat','evlmi','ezlmi',
     $     'lensrat','KcorBI','KcorBJ','stretch',
     $     'mixdm','hostdm','dhostdm','greydm',
     $     'dstretch','color','dcolor','x1','dx1'/


C
C       The indices in the ZRANGE vector for the lower and higher
C       redshift range.
C
      INTEGER ZRGPAR
      REAL ZLOW, ZHIG
      PARAMETER (ZLOW=1, ZHIG=2, ZRGPAR=2)
C
C     VARIABLES
C

C     ****v* snalys/SN_NR [1.0] *
C   NAME
C     SN_NR -- number of events in data file
C
C   DESCRIPTION
C     The number of supernovae events stored in the data
C     structure.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-03-22
C
C   NOTE
C     The value of this variable should be set in init.f.
C
C     ***

C     ****v* snalys/NR_SUPER [1.0] *
C   NAME
C     NR_SUPER -- number of events that is requested
C
C   DESCRIPTION
C     The number of supernovae events that the user wish to
C     use.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-03-22
C
C   NOTE
C     The value of this variable should be set in init.f.
C
C     ***


C     ****v* snalys/SN_USED [1.0] *
C   NAME
C     SN_USED -- number of events that is requested
C
C   DESCRIPTION
C     The number of supernova events that were actually used in the
C     analysis.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-03-22
C
C     ***


C     ****v* snalys/SN_DATA [1.0] *
C   NAME
C     SN_DATA -- the supernova data array
C
C   DESCRIPTION
C     A REAL matrix of dimensions (ENR,MAX_SN) that contains all the
C     data for all the supernovae. Each columns represent one supernova.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-03-22
C
C   NOTE
C     This array should be built in init.f
C
C     ***


C     ****v* snalys/ZRANGE [1.0] *
C   NAME
C     ZRANGE -- the requested red shift range
C
C   DESCRIPTION
C     a REAL vector that contains the min and max red shifts that
C     that sets the limited range of events if such a range is
C     requested.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-03-22
C
C   NOTE
C     These values should be set in init.f
C
C     ***

      INTEGER SN_NR, NR_SUPER, SN_USED
      COMMON /ML_INT/SN_NR, NR_SUPER, SN_USED
      REAL SN_DATA, ZRANGE, GAUSS_SIGMA, MULTI_ETDEL,
     $     MULTI_EIMSEP, MULTI_EIMRAT, MULTI_EZL, MULTI_EVEL,
     $     MULTI_EMODEL
      COMMON /SN_DATA/SN_DATA(ENR,MAX_SN),
     $     ZRANGE(ZRGPAR), GAUSS_SIGMA,MULTI_ETDEL,
     $     MULTI_EIMSEP, MULTI_EIMRAT, MULTI_EZL, MULTI_EVEL,
     $     MULTI_EMODEL



C     ****v* snalys/GAUSSLINPARS [1.0] *
C   NAME
C     GAUSSLINPARS -- values of the three pdf parameters in GAUSSLIN
C
C   DESCRIPTION
C     REAL array containing the values of the three pdf parameters in
C     GAUSSLIN for each red shift, and for a specified amount of
C     fraction of point masses in the universe.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-10-13
C
C   NOTE
C     These values should be set in init.f
C
C     ***
      REAL GAUSSLINPARS
      COMMON /GAUSSLINPARS/GAUSSLINPARS(4,MAX_SN)
      INTEGER GLPM0, GLPS, GLPB, GLMLIM
      PARAMETER ( GLPM0 = 1 , GLPS = 2 , GLPB = 3 , GLMLIM = 4 )

      CHARACTER*3 LENSTYPE
      COMMON /LENSMODEL/ LENSTYPE

      REAL IBANDMAX,JBANDMAX
      COMMON /DISCOVERY/ IBANDMAX, JBANDMAX



C     ****v* snalys/BOSPRIOR [1.0] *
C   NAME
C     BOSPRIOR -- priors from baryon oscillations
C
C   DESCRIPTION
C     Cosmology prior from baryon oscillations that are applied
C     to the luminosity distance.
C
C   AUTHOR
C     Ariel Goobar (ariel@physto.se)
C
C   CREATION DATE
C     2005-08-25
C
C   NOTE
C     These values should be set in init.f
C
C     ***
      REAL BOSPRI,DBOSPR
      COMMON /BOSPRIOR/ BOSPRI,DBOSPR

C     ****v* snalys/CMBPRIOR [1.0] *
C   NAME
C     CMBRIOR -- priors from CMB
C
C   DESCRIPTION
C     Cosmology prior from CMB that are applied
C     to the luminosity distance.
C
C   AUTHOR
C     Ariel Goobar (ariel@physto.se)
C
C   CREATION DATE
C     2006-10-18
C
C   NOTE
C     These values should be set in init.f
C
C     ***
      REAL CMBPRI,DCMBPR
      COMMON /CMBPRIOR/ CMBPRI,DCMBPR

C     ****v* snalys/OMPRIOR [1.0] *
C   NAME
C     OMRIOR -- priors on Omega_M
C
C   DESCRIPTION
C     Cosmology prior on Omega_Mthat are applied
C     to the luminosity distance
C
C   AUTHOR
C     Ariel Goobar (ariel@physto.se)
C
C   CREATION DATE
C     2006-10-24
C
C   NOTE
C     These values should be set in init.f
C
C     ***
      REAL OMPRI,DOMPR
      COMMON /OMPRIOR/ OMPRI,DOMPR


C
C End of supernova.inc
C
