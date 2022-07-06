      SUBROUTINE NTECLC ( PRE, TEM, GAMMA, ANLTE, CNLTE, SQ )
C
C VERSION
C     25-OCT-02  AD  Remove checks on GAS ID before calling QTNTE.
C     14-DEC-01  AD  Set INTE=0 to avoid warning messages
C     26-OCT-00  AD  Allow for VT=0
C     09-AUG-00  AD  Change TEMNTE to VT-KT rather than absolute VT
C     09-AUG-99  AD  Make explicit DP components of RUPP, RLOW calc.
C     21-OCT-97  AD  Change ANLTE,CNLTE to DOUBLE PRECISION
C     04-APR-97  AD  Check QPRNTE flag before calc. Qtot.
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     21-SEP-96  AD  Revised to use GENLN2 v.4 non-LTE
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Calculate various non-LTE parameters for line
C     Called by ADJUST.
C
      IMPLICIT NONE
C
C ARGUMENTS
      REAL             PRE    !  I  Path pressure [atm]
      REAL             TEM    !  I  Path kinetic temperature [K]
      DOUBLE PRECISION GAMMA  !  I  exp - ( hcv/kT )
      DOUBLE PRECISION ANLTE  !  O  Non-LTE Correction factor for k absorption
      DOUBLE PRECISION CNLTE  !  O  Non-LTE Correction factor for c absorption
      REAL             SQ     !  O  Ratio tot.partition sum @296K/tps@path.tem
C
      EXTERNAL
     &  LOOKUP ! General purpose interpolation routine
     &, QTNTE  ! Calculate total internal partition sums for non-LTE.
     &, QTFCT  ! Calculate total internal partition sums for line data.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'hitcom.inc' ! HITRAN line database variables
      INCLUDE 'ntecom.inc' ! Vibrational levels for Nonlte calcs
C
C LOCAL CONSTANTS 
      REAL          PRELIM          ! Press [atm] defining high alt limit
        PARAMETER ( PRELIM = 1.0E-5 )
      REAL          TEMLIM          ! Default Max Vib.Temp for high alt.emission
        PARAMETER ( TEMLIM = 275.0 )
      REAL          FACLIM          ! Min VT/KT ratio indisting. from VT=0
        PARAMETER ( FACLIM = 0.01 )
C
C LOCAL VARIABLES
      INTEGER INTE         ! Index of Non-LTE profile
      INTEGER J            ! Save starting index for LOOKUP
      REAL    ENBOT        ! Energy [/cm] of lower level
      REAL    ENTOP        ! Energy [/cm] of upper level
      REAL    LNP          ! Log(p/mb) for interpolating Vib.Temp. profiles
      REAL    QVTEM        ! Vib.Part.Fn for path
      REAL    TEMDEF       ! Default vibrational temperature
      REAL    VTBOT        ! Vib.Temp of lower level      
      REAL    VTTOP        ! Vib.Temp of upper level      
      DOUBLE PRECISION RLOW ! Ratio of lower state population to LTE value
      DOUBLE PRECISION RUPP ! Ratio of upper state population to LTE value
      LOGICAL BOTDAT       ! T = Vib.Temp available for lower level
      LOGICAL TOPDAT       ! T = Vib.Temp available for upper level
        DATA J / 1 /
        SAVE J
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      INTE = 0             ! Unnecessary, but avoids warning messages
C
      IF ( TEM .GT. TEMLIM .AND. PRE .LT. PRELIM ) THEN
        TEMDEF = TEMLIM
      ELSE
        TEMDEF = TEM
      END IF
      LNP = LOG ( PRE * ATMB )
C
      BOTDAT = ILSNTE .NE. 0
      IF ( BOTDAT ) THEN
        CALL LOOKUP ( LNP, VTBOT, LNPATM, TEMNTE(1,ILSNTE), NATM, J, 1)
        VTBOT = VTBOT + TEM
        IF ( QPRNTE(ILSNTE) ) 
     &    CALL LOOKUP (LNP, QVTEM, LNPATM, QFNNTE(1,ILSNTE), NATM, J, 1)
        ENBOT  = ENGNTE(ILSNTE)
        INTE = ILSNTE
      ELSE
        VTBOT = TEMDEF
        ENBOT = ELS
      END IF               
      IF ( VTBOT .GT. FACLIM * TEM ) THEN
        RLOW = DEXP ( -C2 * DBLE(ENBOT) * 
     &              DBLE(TEM-VTBOT) / DBLE(TEM*VTBOT) )
      ELSE
        RLOW = 0.0D0
      END IF
C
C Note that this routine is only called if IUSNTE or ILSNTE .NE. 0, so either 
C TOPDAT or BOTDAT or both will be TRUE.
C If there is no data for the top level, there must be data for the bottom
C level so assume the top level is in LTE with the bottom (RUPP=RLOW)
C Note that QFNNTE is the same for all states of any molecule,isopte combination
C so no problem if QVTEM is determined from both top and bottom levels.
C
      TOPDAT = IUSNTE .NE. 0
      IF ( TOPDAT ) THEN
        CALL LOOKUP ( LNP, VTTOP, LNPATM, TEMNTE(1,IUSNTE), NATM, J, 1)
        VTTOP = VTTOP + TEM
        IF ( QPRNTE(IUSNTE) ) 
     &    CALL LOOKUP (LNP, QVTEM, LNPATM, QFNNTE(1,IUSNTE), NATM, J, 1)
        ENTOP  = ENGNTE(IUSNTE)
        INTE = IUSNTE
        IF ( VTTOP .GT. FACLIM * TEM ) THEN
          RUPP = DEXP ( -C2 * DBLE(ENTOP) * 
     &                  DBLE(TEM-VTTOP) / DBLE(TEM*VTTOP) )
        ELSE 
          RUPP = 0.0D0
        END IF
      ELSE
        RUPP = RLOW
      END IF
C
      ANLTE = ( RLOW - RUPP * GAMMA ) / ( 1.D0 - GAMMA )
      CNLTE = RUPP
C
C Partition functions  
C
      IF ( QPRNTE(INTE) ) THEN
        CALL QTNTE ( IDGAS, ISO, TEM, QVTEM, SQ )
      ELSE
        CALL QTFCT ( IDGAS, ISO, TEM, SQ )
      END IF
C
      END
