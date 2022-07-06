C ATMCOM.INC ------------------------------------------------------------------
      COMMON / ATMCOM / NATM, CHGATM, 
     &                  HGTATM, PREATM, TEMATM, VMRATM, DNSATM,
     &                  RFRATM, LNDATM, LNPATM, LNRATM, LNVATM, DSHATM
C
C VERSION
C     20-JAN-03  AD  Add CHGATM
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Atmospheric profile data for RFM.
C     Values are defined at base of each layer
C     Loaded by INPATM and ATMAUX
C
      INTEGER NATM                  !  No. atmospheric layers for profiles
      LOGICAL CHGATM(MAXATM)        !  T=atmospheric layer changed
      REAL    HGTATM(MAXATM)        !  Altitudes of each layer [km]
      REAL    PREATM(MAXATM)        !  Pressures of each layer [mb]
      REAL    TEMATM(MAXATM)        !  Temperatures of each layer [K]
      REAL    VMRATM(MAXATM,MAXGAS) !  Volume mixing ratios [ppv]
      REAL    DNSATM(MAXATM)        !  Number Density [/cm^3]
      REAL    RFRATM(MAXATM)        !  Refractivity profile (ATMRFR)
      REAL    LNDATM(MAXATM)        !  Ln ( DNSATM )
      REAL    LNPATM(MAXATM)        !  Ln ( PREATM )
      REAL    LNRATM(MAXATM)        !  Ln ( RFRATM )
      REAL    LNVATM(MAXATM,MAXGAS) !  Ln ( VMRATM )
      REAL    DSHATM(MAXATM)        !  Density Scale Height [km]
C
C -----------------------------------------------------------------------------
