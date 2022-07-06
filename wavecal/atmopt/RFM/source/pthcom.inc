C PTHCOM.INC ------------------------------------------------------------------
      COMMON / PTHCOM / NPTH, NCLC, IATPTH, IDRPTH, IGSPTH, 
     &                  ISCPTH, ITNPTH, IVJPTH, NLNPTH, 
     &          AMTPTH, PREPTH, PPAPTH, PSIPTH, RAYPTH, TEMPTH,
     &                  CLCPTH
C
C VERSION
C     03-MAY-05  AD  Add IVJPTH
C     02-JAN-00  AD  Add IDRPTH, PSIPTH. Remove redundant IDIPTH 
C     31-JAN-99  AD  Remove JACPTH
C     02-JAN-99  AD  Add JACPTH
C     10-JUL-98  AD  Remove redundant IFLPTH
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Common variables associated with RFM paths (=atm.path * absorbers)
C     Set up in TANPTH.FOR, reshuffled in PTHSRT.FOR (unless CLCFLG).
C
      INTEGER  NPTH           ! No. of paths used
      INTEGER  NCLC           ! No. of paths for explicit line calculations
      INTEGER  IATPTH(MAXPTH) ! Atmospheric profile layer# for this path
      INTEGER  IDRPTH(MAXPTH) ! 0=downward from obs, or irrelevant, 1=upward
      INTEGER  IGSPTH(MAXPTH) ! Absorbing gas# for this path
      INTEGER  ISCPTH(MAXPTH) ! Index of path to be used for scaling absorption
      INTEGER  ITNPTH(MAXPTH) ! Tangent Height# for this path
      INTEGER  IVJPTH(MAXPTH) ! Index of Jacobian for perturbed VT path (or 0)
      INTEGER  NLNPTH(MAXPTH) ! No.lines stored in cyclic buffer 
      REAL     AMTPTH(MAXPTH) ! Path gas amount [ kmol / cm^2 ] 
      REAL     PREPTH(MAXPTH) ! Pressure [atm]
      REAL     PPAPTH(MAXPTH) ! Partial pressure [atm]
      REAL     PSIPTH(MAXPTH) ! Horiz.angle [deg] of low alt end of path
      REAL     RAYPTH(MAXPTH) ! Ray length in path [km]
      REAL     TEMPTH(MAXPTH) ! Temperature [K]	
      LOGICAL  CLCPTH(MAXPTH) ! T=perform line-by-line calc.for this path
C
C -----------------------------------------------------------------------------
