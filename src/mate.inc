C     ****** snalys/MATEINC [1.0] *
C   NAME
C     MATEINC -- common variables used in MATE()
C
C   AUTHOR
C     Martin Eriksson (mate@physto.se)
C
C   CREATION DATE
C     2001-09-07
C
C     ***
      INTEGER GLOB_N_ALPHA, GLOB_N_OMM, GLOB_N_ZP
      PARAMETER (GLOB_N_ALPHA=10, GLOB_N_OMM=11, GLOB_N_ZP=69)
      DOUBLE PRECISION GLOB_WINDEX(GLOB_N_ALPHA+GLOB_N_OMM+GLOB_N_ZP),
     $	GLOB_WCUBE(GLOB_N_ALPHA, GLOB_N_OMM, GLOB_N_ZP)
      CHARACTER*(256) WCUBEFILE
      COMMON /QUINTESSENCE/GLOB_WINDEX, GLOB_WCUBE, WCUBEFILE
C
C End of mate.inc
C
