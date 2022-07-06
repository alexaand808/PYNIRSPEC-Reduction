C JACCOM.INC -------------------------------------------------------------------
      COMMON / JACCOM / NJAC, IGSJAC, IATJAC, ILOJAC, IUPJAC, ITNJAC
C
C VERSION
C      04-FEB-99  AD  Add ILOJAC, IUPJAC
C      01-JAN-99  AD  Original.
C
C DESCRIPTION
C     Jacobian data required for RFM calculations.
C     Loaded by JACCHK.
C
      INTEGER  NJAC           ! No. retrieval elements for Jacobian calc.
      INTEGER  IGSJAC(MAXJAC) ! Index of absorb.species for each Jac.Rtvl.El.
      INTEGER  IATJAC(MAXJAC) ! Atmos level# of each Jacobian rtvl. element
      INTEGER  ILOJAC(MAXJAC) ! Lower Atmos level# for perturbation=0
      INTEGER  IUPJAC(MAXJAC) ! Upper Atmos level# for perturbation=0
      INTEGER  ITNJAC(MAXTAN,MAXJAC) ! Index of extra Tan.pth for calc.
                                     ! dRad(ITAN)/dRtv(IJAC)
C
C ----------------------------------------------------------------------------
