C     ****d* snalys/NR_PAR [1.0] *
C   NAME
C     NR_PAR -- the maximum number of parameters
C
C   DESCRIPTION
C     This constant defines the maximum number of estimators.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-03-22
C
C   HISTORY
C     2003-03-31 - Added SALP, increased NR_PAR with 1.
C     2003-09-19 - Added O_RC, increased NR_PAR with 1.
C     2004-01-15 - Added B_1, B_2, increased NR_PAR with 2.
C     2004-02-16 _ Added B_inv, B_4, increased NR_PAR with 2.
C     2004-05-17 - Added W0, W1, QHP and ZX, increased NR_PAR with 4.
C     2006-09-13 - Added L1 and L2, increased NR_PAR with 2. (Jakob J)
C     2010-02-01 - Added FV0 and H0BAR, increased NR_PAR with 2 (Rahman)
C
C     ***
      INTEGER NR_PAR
      PARAMETER ( NR_PAR=22 )



C     ****d* snalys/MAX_LEN [1.0] *
C   NAME
C     MAX_LEN -- the maximum number of grid points
C
C   DESCRIPTION
C     This constant defines the maximum number of possible points
C     in the grid search.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-03-22
C
C     ***
      INTEGER MAX_LEN
      PARAMETER ( MAX_LEN=25000000 )



C     ****v* snalys/NR_VAR [1.0] *
C   NAME
C     NR_VAR -- the number of estimators
C
C   DESCRIPTION
C     The number of cosmological parameters that is to be estimated.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-03-22
C
C     ***
      INTEGER NR_VAR
      COMMON /NR_VAR/NR_VAR




C     ****v* snalys/COS_INT [1.0] *
C   NAME
C     COS_INT -- the parameters intervals
C
C   DESCRIPTION
C     In this matrix the intervals for the cosmological parameters,
C     used in the grid_search() is defined.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-03-22
C
C   NOTE
C     The dimension of this matrix is NR_PARx10.
C
C     ***
      REAL COS_INT
      COMMON /ML_COS/COS_INT(NR_PAR,10)

C     ****v* snalys/COSFP [1.0] *
C   NAME
C     COSFP -- fraction of point masses in the Universe
C
C   DESCRIPTION
C     This variable holds the fraction of point masses (in
C     percent )in the Universe.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-10-13 
C
C     ***
      REAL COSFP
      COMMON /COSFP/COSFP



C     ****v* snalys/INDEXVARIABLES [1.0] *
C   NAME
C     INDEXVARIABLES -- simple integers, one for each parameter
C
C   DESCRIPTION
C     These variable defines which row in the cosmology interval
C     that represents what parameter. They are also used to keep
C     track of the dimension in the grid array.
C
C     O_M   - Omega Mass
C     O_X   - Omega Lambda
C     A_X   - Equation of state parameter w = w_0
C     A_2   - If a linear time dependent equation of state parameter
C             is wanted w = w_0 + w_1*z
C     M_SC  - The intercept often called "M-Script".
C     ALP   - The quintessence parameter alpha.
C     SALP  - The stretch parameter alpha.
C     H_0   - The Hubble parameter used for the DO_MULTI mode.
C     O_RC  - Omega_rc = 1/(2r_cH_0)**2 for the 5DIM mode
C     B_1   - Linear coeff. in Taylor expansion of dark energy density.
C     B_2   - Quadratic coeff. in Taylor expansion of dark energy density.
C     B_inv - Inverse coeff. (x^-1) in Taylor expansion of dark energy density.
C     B_4   - Quartic coeff. in Taylor expansion of dark energydensity. 
C     L_1   - Linder parametrization w(z)=L_1+L_2*z/(1+z)
C     L_2
C     W_0   - Late time EOS state asymptote in Hannestad parametrization (HP).
C     W_1   - Early time EOS state asymptote in (HP).
C     Q_HP  - Transition speed in (HP).
C     Z_X   - Transition redshift in (HP).
C     X_DA  - parameter in higher dimensions
C     FV_0  - Wiltshare parameter
C     H_0bar- Wiltshare parameter
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-03-22
C
C   NOTE
C     The values of these should be set once and for all in init.f.
C     PLEASE KEEP THE UNDERSCORE NOTATION WHILE ADDING NEW VARIABLES.
C
C   HISTORY
C     2001-09-07 - added ALP.
C     2002-01-08 - added H_0
C     2003-03-31 - added SALP
C     2003-09-18 - added O_RC (Ariel G.)
C     2003-09-19 - removed the obsolete MSLIMS and ASLIMS
C     2004-01-15 - added the parameters B_1, B_2 (Jakob J.)
C     2004-02-16 - added the parameters B_inv, B_4 (Jakob J.)
C     2004-05-17 - added the parameters W_0, W_1, Q_HP, Z_X (Edvard M.)
C     2005-10-19 - added X_DA (Ariel)
C     2006-09-13 - added L_1, L_2 (Jakob J)
C     2010-02-01 - added FV_0, H_0BAR (Rahman)
C
C     ***
      INTEGER O_M,O_X,A_X,A_2,M_SC,ALP,H_0,SALP,O_RC,B_1,B_2,B_inv,
     $     B_4,W_0,W_1,Q_HP,Z_X,X_DA,L_1,L_2,FV_0,H_0BAR
      CHARACTER*(10) P_NAMES(NR_PAR)
      COMMON /COS_IDX/O_M,O_X,A_X,A_2,M_SC,ALP,H_0,SALP,O_RC,
     $     B_1,B_2,B_inv,B_4,W_0,W_1,Q_HP,Z_X,X_DA,L_1,L_2,
     $     FV_0,H_0BAR,P_NAMES



C     ****v* snalys/COMMON_COSMOLOGY [1.0] *
C   NAME
C    COMMON_COSMOLOGY -- values of cosmological parameters.
C
C   DESCRIPTION
C     These common variabels for the cosmological parameters are used
C     when the luminosity distance is calculated since the integration
C     function requires that only the integration variable may be passed
C     to the integrand.
C     
C     GLOB_OMM    - Omega Mass 
C     GLOB_OMX	  - Omega Lambda
C     GLOB_AX	  - Equation of state parameter w = w_0
C     GLOB_A2	  - If a linear time dependent equation of state parameter
C                   is wanted w = w_0 + w_1*z
C     GLOB_ALPHA  - The quintessence parameter alpha
C     GLOB_SALP   - The stretch alpha
C     GLOB_MSCRIPT- The intercept often called "M-Script".
C     GLOB_H0     - The Hubble parameter (only in use for the DO_MULTI
C                   option)
C     GLOB_ORC    - Omega_rc (only in use for 5DIM option)
C     GLOB_B1     - Linear coeff.
C     GLOB_B2     - Quadratic coeff.
C     GLOB_Binv   - Inverse coeff.
C     GLOB_B4     - Quartic coeff.
C     GLOB_W0     - Late time EOS state asymptote in (HP)
C     GLOB_W1     - Early time EOS state asymptote in (HP)
C     GLOB_QHP    - Transition speed in (HP)
C     GLOB_ZX     - Transistion redshift in (HP)
C     GLOB_XDA    - alpha parameter in xdim case
C     GLOB_L1	  - Linder parametrization w(z)=L_1+L_2*z/(1+z)
C     GLOB_L2
C     GLOB_FV0    - Wiltshare
C     GLOB_H0BAR
C     
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-03-22
C
C   NOTE
C     The values of this variables should only be altered by
C     UPDATE_GLOBALS, or in RVARVAL or RMPAR if the parameters is not going
C     to be fitted.
C
C   HISTORY
C     2001-09-07 - Added GLOB_ALPHA
C     2001-09-27 - Added GLOB_OMK
C     2002-01-08 - Added GLOB_H0
C     2003-03-31 - Added GLOB_SALP
C     2003-09-18 - Added GLOB_ORC (Ariel)
C     2004-01-15 - Added GLOB_B1, GLOB_B2 (Jakob)
C     2004-02-16 - Added GLOB_Binv, GLOB_B4
C     2004-05-17 - Added GLOB_W0, GLOB_W1, GLOB_Q, GLOB_ZX (Edvard)
C     2005-10-18 - Added GLOB_XDA (malc & ariel)  
C     2006-09-13 - Added GLOB_L1, GLOB_L2 (Jakob J)
C     2010-02-01 - Added GLOB_FV0, GLOB_H0BAR (Rahman)
C
C     ***
      DOUBLE PRECISION GLOB_OMM,GLOB_OMX,GLOB_AX,GLOB_A2,GLOB_MSCRIPT,
     $     GLOB_ALPHA, GLOB_OMK, GLOB_H0, GLOB_SALP, GLOB_ORC, GLOB_B1,
     $     GLOB_B2, GLOB_Binv, GLOB_B4, 
     $     GLOB_W0, GLOB_W1, GLOB_QHP, GLOB_ZX,GLOB_XDA, 
     $	   GLOB_L1, GLOB_L2,
     $     GLOB_FV0, GLOB_H0BAR
      COMMON /GLOBAL_COS/GLOB_OMM,GLOB_OMX,GLOB_AX,GLOB_A2,GLOB_MSCRIPT,
     $     GLOB_ALPHA, GLOB_OMK, GLOB_H0, GLOB_SALP, GLOB_ORC, GLOB_B1,
     $     GLOB_B2, GLOB_Binv, GLOB_B4,
     $     GLOB_W0, GLOB_W1, GLOB_QHP, GLOB_ZX,GLOB_XDA,
     $	   GLOB_L1,GLOB_L2,
     $     GLOB_FV0, GLOB_H0BAR



C     ****v* snalys/NUISANCE [1.0] *
C   NAME
C     NUISANCE
C
C   DESCRIPTION
C     The structure for saving the value of the nuisance parameters
C     for the best fit. Any nuisance parameters should be added here.
C     
C     BEST_SCRIPTM - The best scriptM value is saved here.
C     THIS_SCRIPTM - A temporary value where the best value of
C                    scriptM for the current marginalization is stored.
C     BEST_ALPHA   - See BEST_SCRIPTM
C     THIS_ALPHA   - See THIS_SCRIPTM
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2004-10-27
C
C     ***
      REAL BEST_SCRIPTM, THIS_SCRIPTM, BEST_ALPHA,
     $     THIS_ALPHA
      COMMON /NUISANCE/BEST_SCRIPTM, THIS_SCRIPTM, THIS_ALPHA,
     $     BEST_ALPHA
      
C
C End of cosmo.inc
C
