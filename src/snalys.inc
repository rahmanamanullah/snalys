cc snalys.inc
cc
cc This is an include file for the ml_analysis program, where common
cc parameters are defined, and variables that are only set once
cc in the program are declared.
cc
cc Rahman Amanullah, rahman@physto.se, 000324

      REAL PI
      PARAMETER (PI = 3.14159265)


C     ****v* snalys/CALCDERIV [1.0] *
C   NAME
C     CALCDERIV -- calculated the derivatives in the minimization
C
C   DESCRIPTION
C     This LOGICAL variable determines whether the derivatives of the
C     likelihood function with respect to the cosmological parameters
C     should be calculated.
C
C   AUTHOR
C     Rahman Amanullah (rahman@physto.se)
C
C   CREATION DATE
C     2001-03-24
C
C   NOTE
C     This variable should be set in init.f
C
C     ***
      LOGICAL CALCDERIV
      COMMON /CALCDERIV/CALCDERIV

C
C End of snalys.inc.
C
