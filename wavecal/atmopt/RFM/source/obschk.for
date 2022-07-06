      SUBROUTINE OBSCHK ( ALTTST, FAIL, ERRMSG )
C
C VERSION
C     11-DEC-07  AD  Fix Bug#68 - returning with FAIL=TRUE when actually OK
C     30-DEC-03  AD  Replace OBSATM with ATMLEV
C     18-MAR-01  AD  Remove NEWLEV argument. Remove CHKLIM
C     18-JUL-98  AD  Original.
C
C DESCRIPTION
C     Check and set observer altitude.
C     Called by INPOBS following *OBS marker in Driver table.
C
      IMPLICIT NONE
C
C ARGUMENTS
      REAL         ALTTST !  I  Observer altitude to be tested and loaded
      LOGICAL      FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  ATMLEV ! Find/insert atmospheric level for given altitude.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'obscom.inc' ! Observer Position
C
C LOCAL VARIABLES 
      INTEGER  IATM        ! Index of atmospheric level for observer
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      FAIL = .FALSE.
C
      IF ( ALTTST .LT. HGTATM(1) ) THEN
        WRITE ( ERRMSG, '(A,F10.3,A)' )  
     &    'F-OBSCHK: Observer Altitude below base of atmosphere, =',
     &    HGTATM(1), ' [km]'              
        FAIL = .TRUE.
      ELSE IF ( ALTTST .GT. HGTATM(NATM) ) THEN ! Observer above top of atm.
        IF ( ZENFLG ) THEN                      ! Fatal error if looking upward
          WRITE ( ERRMSG, '(A,F10.3,A)' )
     &    'F-OBSCHK: Observer Altitude above top of atmosphere, =',
     &    HGTATM(NATM), ' [km]'              
          FAIL = .TRUE.
        ELSE                                    ! OK for other cases
          ALTOBS = ALTTST
          IATOBS = NATM
        END IF
      ELSE 
C Note that ATMLEV itself adjusts IATOBS so cannot use this as argument
        CALL ATMLEV ( ALTTST, IATM, FAIL, ERRMSG )  
        ALTOBS = ALTTST
        IATOBS = IATM
      END IF
C
      END
