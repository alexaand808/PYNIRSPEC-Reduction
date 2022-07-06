      REAL FUNCTION TPSLKP ( IDGAS, ISO, TEM )
C
C VERSION
C     07-AUG-13  AD  Remove redundant local variables ENDSEC, RECORD
C     13-FEB-03  AD  Original.
C
C DESCRIPTION
C     Look up tabulated TIPS factor
C     Called by QTFCT for each molecule/isotope with tabulated data.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER IDGAS !  I  HITRAN gas ID
      INTEGER ISO   !  I  HITRAN isotope ID
      REAL    TEM   !  I  Path temperature [K]
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES 
      INCLUDE 'tpscom.inc' ! Tabulated TIPS data
C
C LOCAL VARIABLES
      INTEGER      IPT      ! Lower intvl index in tabulation temperature axis
      INTEGER      ITPS     ! Index of molec/isotope in *TPS variables
      REAL         XPT      ! Position of Ref.Temp in temperature axis intvl.
      REAL         XTEM     ! Position of Ref.Temp in temperature axis
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      ITPS = IDXTPS(ISO,IDGAS)
      XTEM = 1.0 + ( TEM - PT1TPS(ITPS) ) / PTDTPS(ITPS)
C If Temperature outside tabulated range, limit to end values and print
C warning message to terminal
      IF ( XTEM .LT. 1.0 .OR. XTEM .GT. NPTTPS(ITPS) ) THEN
        WRITE ( *, * ) 
     &    'W-TPSLKP: TEM outside tabulated TIPS, IDGAS,ISO,TEM=',
     &     IDGAS, ISO, TEM
        XTEM = MIN ( XTEM, FLOAT ( NPTTPS(ITPS) ) )
        XTEM = MAX ( XTEM, 1.0 )
      END IF
      IPT = INT ( XTEM )
      XPT = XTEM - FLOAT ( IPT )
      TPSLKP = QFCTPS(IPT,ITPS) * ( 1.0 - XPT ) +
     &         QFCTPS(IPT+1,ITPS) * XPT
C
      END
