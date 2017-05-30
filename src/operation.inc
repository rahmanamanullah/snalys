      INTEGER DO_GRID ,    ! Do a grid search
     $     DO_GENETIC ,    ! Do a minimization according to Genetic's method.
     $     DO_POWELL ,     ! Do a minimization according to Powell's method.
     $     DO_DAVIDON ,    ! Do a minimization according to Davidon's method.
     $     ADD_LDM ,       ! Add the contr. of grav. lensing to the magnitudes.
     $     ADD_HOSTDM ,    ! Add the contr. from the host galaxy to the magn.
     $     ADD_MIXDM ,     ! Add the contr. from photon-axion mixing to magn.
     $     NR_OP ,         ! The number of operation booleans.
     $     SUB_INTSIG,     ! If true the SNe are not intr. spread.
     $     MAG_DIST,       ! Calcualte the magnitude residuals.
     $     QUIET,          ! Nothing should be written STDOUT
     $     FLAT,           ! Flat universe assumed
     $     CON_FIND,       ! Use the fancy contour finder (not implemented)
     $     ZLIMIT ,        ! Apply limits in red shift.
     $     MASS_PRIOR ,    ! Apply a Omega_M prior function (remove this).
     $     MATE_EQST ,     ! Fit the Eriksson quintessence model
     $     FITMSC ,        ! Fit the M_script parameter.
     $	   FITSALP ,       ! Fit the stretch alpha parameter.
     $     DO_MULTI,       ! Fit cosmo para from multiple image splitting
     $     XDIM,           ! Fit cosmo para for X-dimension
     $     NRMIN,
     $	   FITS,           ! Output will be saved in a fits file
     $     DO_JACKE,       ! Fit "jacke"-cosmology
     $	   DO_LINDER,	   ! Use Linder parmetrization
     $	   DO_WILT,	   ! Use Wiltshire's cosmology
     $     DO_BOS,         ! Use prior from Baryon Oscillations
     $     DO_CMB,         ! Use prior from CMB
     $     DO_OM,          ! Use prior on OMEGA_M
     $     DO_HP,          ! Fit EOS according to Hannestad parametrization
     $     ADD_GDM,        ! Add the contr. of grey dust to the  magnitudes.
     $	   Z_ERROR,	   ! Add redshift dependent error
     $     USEMERR         ! Consider given measurement errors
      PARAMETER ( NR_OP = 30 , NRMIN = 5 )
      PARAMETER ( DO_GRID=1, DO_GENETIC=2, DO_POWELL=3, DO_DAVIDON=4,
     $     CON_FIND=5, ADD_LDM=6, SUB_INTSIG=7, MAG_DIST=8, QUIET=9,
     $     FLAT=10, ZLIMIT=11 , MASS_PRIOR=12,
     $     MATE_EQST=13, FITMSC=14, DO_MULTI=15, FITSALP=16, XDIM=17, 
     $     FITS=18, DO_JACKE=19, DO_HP=20, ADD_HOSTDM=21,
     $     ADD_MIXDM=22,ADD_GDM=23, DO_BOS=24, DO_LINDER=25, Z_ERROR=26,
     $     DO_CMB=27, DO_OM=28, USEMERR=29, DO_WILT=30)

C     ****v* snalys/OPER [1.0] *
C   NAME
C     OPER -- what operations should be performed
C
C   DESCRIPTION
C     A core LOGICAL vector where the performance of the software is
C     decided.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-03-22
C
C   NOTE
C     The values of this vector is set in init.f
C
C   HISTORY
C     2001-09-07 - added the possibility to fit Eriksson's quintessence.
C     2003-09-23 - added the possibliity to store output as fits.
C     2010-02-01 - added DO_WILT
C
C     ***
      LOGICAL OPER(NR_OP)
      COMMON /OPERATION/OPER

C
C End of operation.inc.
C
