      SUBROUTINE INPATM ( LUNDRV, LUNATM, FAIL, ERRMSG )
C
C VERSION
C     19-SEP-12  AD  Change NAMATM from C*80 to C*200
C     26-DEC-99  AD  Add call to ATMGRA      
C     22-JAN-98  AD  Remove redundant EXTERNAL declarations: RFMLOG, INPFIL
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     28-SEP-96  AD  Remover NGAS argument to ATMAUX
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Read atm.profile files from driver table following *ATM marker
C     Called by RFMINP once.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      LUNDRV  !  I  LUN for Driver File
      INTEGER      LUNATM  !  I  LUN for Atmospheric profiles 
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  ATMAUX ! Set up auxiliary profiles of atmospheric parameters
     &, ATMFIL ! Read file containing atmospheric profiles.
     &, ATMGRA ! Interpolate to fill 2-D atmospheric profile field.
     &, ATMLAY ! Subdivide atmospheric profile layers
     &, NXTFLD ! Load next field from section of RFM driver file
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'rfmfil.inc' ! Standard filenames of RFM I/O files
C
C COMMON VARIABLES
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'flgcom.inc' ! RFM option flags
C
C LOCAL VARIABLES
      INTEGER  IGAS            ! Gas profile counter (1:NGAS)
      INTEGER  LENGTH          ! Length of filename read from driver file
      LOGICAL  GOTPRE          ! TRUE if a pressure profile has been stored
      LOGICAL  GOTTEM          ! TRUE if a temperature profile has been stored
      LOGICAL  GOTGAS(MAXGAS)  ! TRUE if a vmr profile has been stored
      LOGICAL  FIRST           ! TRUE until first profile of first file read in
      CHARACTER*200 NAMATM     ! Name of ATM file read from Driver Table
        DATA FIRST / .TRUE. /
        DATA GOTTEM / .FALSE. /
        DATA GOTGAS / MAXGAS*.FALSE. /
        DATA GOTPRE / .FALSE. /
        SAVE FIRST, GOTTEM, GOTGAS, GOTPRE
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Jump here for each new field in *ATM section. 
C
  100 CALL NXTFLD ( LUNDRV, NAMATM, LENGTH, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      IF ( LENGTH .NE. 0 ) THEN                  ! Marker for end of section
        CALL ATMFIL ( LUNATM, NAMATM(1:LENGTH), FIRST, GOTPRE, 
     &                GOTTEM, GOTGAS, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        FIRST = .FALSE.
        GOTO 100
      END IF
C
C Check all required profiles have been loaded
C
      FAIL = .TRUE.
      IF ( .NOT. GOTPRE ) THEN
        ERRMSG = 'F-INPATM: No Pressure profile specified'
        RETURN
      ELSE IF ( .NOT. GOTTEM ) THEN
        ERRMSG = 'F-INPATM: No Temperature profile specified'
        RETURN
      ELSE
        DO IGAS = 1, NGAS
          IF ( .NOT. GOTGAS(IGAS) ) THEN
            ERRMSG = 'F-INPATM: No '//CODGAS(IGAS)//
     &               ' profile specified'
            RETURN
          END IF
        END DO
      END IF
      FAIL = .FALSE.      
C
C If horizontal gradients used, interpolate to fill 2-D field
C
      IF ( GRAFLG ) CALL ATMGRA
C
C Sub-divide profile layers if required
C
      IF ( LAYFLG ) THEN 
        CALL ATMLAY ( FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
C
C Calculate auxiliary profiles
C
      CALL ATMAUX 
C
      FAIL = .FALSE.
      RETURN
C
      END
