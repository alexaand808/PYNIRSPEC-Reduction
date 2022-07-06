C TPSCOM.INC ------------------------------------------------------------------
      COMMON / TPSCOM / NTPS, IDGTPS, ISOTPS, IDXTPS, NPTTPS, 
     &                  PT1TPS, PTDTPS, REFTPS, QFCTPS
C
C VERSION
C     24-FEB-03  AD  Original.
C
C DESCRIPTION
C     Common variables associated with RFM tabulated TIPS data.
C     Loaded by TPSFIL. NTPS initialised in INPTPS.
C
      INTEGER      NTPS            ! No. of tabulated TIPS sets stored
      INTEGER      IDGTPS(MAXTPS)  ! Gas ID of each set of TIPS data
      INTEGER      ISOTPS(MAXTPS)  ! Isotope of each set of TIPS data
      INTEGER      NPTTPS(MAXTPS)  ! No.tabulation points for TIPS data
      INTEGER      IDXTPS(MAXISO,MAXHLN) ! Index of tabulated data,0=none
      REAL         PT1TPS(MAXTPS)  ! 1st tabulation temperature.
      REAL         PTDTPS(MAXTPS)  ! Tabulation temperature increment.
      REAL         REFTPS(MAXTPS)  ! TIPS value at reference temperature
      REAL         QFCTPS(MAXTPT,MAXTPS)  ! Tabulated TIPS values
C
C -----------------------------------------------------------------------------
