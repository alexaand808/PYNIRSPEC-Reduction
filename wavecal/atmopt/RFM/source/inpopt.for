      SUBROUTINE INPOPT ( LUNDRV, FAIL, ERRMSG )
C
C VERSION
C     07-AUG-13  AD  Change FIELD from C*80 to C*200
C     08-JUN-01  AD  Also check for '*' if LOS flag enabled
C     03-JAN-99  AD  Check for '*' if Jacobians required
C     14-JUL-98  AD  Minor change to error message for extra field
C     04-DEC-97  AD  Add GOTOPT to check for duplicate section
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Read user-defined Optical Depth filenames from RFM driver table.
C     Called by RFMINP following *OPT marker in Driver table.
C     If OPTFLG is TRUE read filename for OPT data.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      LUNDRV !  I  LUN for driver table
      LOGICAL      FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  NXTFLD ! Load next field from section of RFM driver file
     &, RFMLOG ! Write text message to RFM log file
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'outcom.inc' ! RFM output file data
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL VARIABLES
      INTEGER       IPT     !  Pointer to '*' in supplied filename
      INTEGER       LENGTH  !  No.characters in FIELD
      CHARACTER*200 FIELD   !  Field read from driver table
      CHARACTER*132 MESSGE  !  Message sent to log file
      LOGICAL       GOTOPT  !  True if section already read
        DATA GOTOPT / .FALSE. /
        SAVE GOTOPT
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      IF ( GOTOPT ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-INPOPT: Duplicate *OPT section in Driver Table'
        RETURN
      ELSE
        GOTOPT = .TRUE.
      END IF
C
      IF ( .NOT. OPTFLG ) THEN
        MESSGE='W-INPOPT: Ignoring *OPT section - no files required'
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        LENGTH = 1
        DO WHILE ( LENGTH .NE. 0 )
          CALL NXTFLD ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        END DO
      ELSE
C
        CALL NXTFLD ( LUNDRV, NAMOPT, LENGTH, FAIL, ERRMSG )
        IF ( LENGTH .GT. 0 ) THEN
          IPT = INDEX ( NAMOPT(1:LENGTH), '*' )
          MESSGE = 'I-INPOPT: User-supplied OPT filename: '
     &             //NAMOPT(1:LENGTH)
          CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
          IF ( NTAN .GT. 1 .AND. IPT .EQ. 0 ) THEN
            FAIL = .TRUE.
            ERRMSG = 'F-INPOPT: No ''*'' symbol in OPT filename '//
     &        'for inserting tangent height info'
            RETURN
          ELSE IF ( NSPC .GT. 1 .AND. IPT .EQ. 0 ) THEN
            FAIL = .TRUE.
            ERRMSG = 'F-INPOPT: No ''*'' symbol in OPT filename '//
     &        'for inserting spectral range info'
            RETURN
          ELSE IF ( ( JACFLG .OR. LOSFLG ) .AND. IPT .EQ. 0 ) THEN
            FAIL = .TRUE.
            ERRMSG = 'F-INPOPT: No ''*'' symbol in OPT filename '//
     &        'for inserting Jacobian element info'
            RETURN
          END IF
        ELSE
          FAIL = .TRUE.
          ERRMSG = 'F-INPOPT: No filename supplied in *OPT section'
          RETURN
        END IF
        CALL NXTFLD ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG )
        IF ( LENGTH .NE. 0 ) THEN
          FAIL = .TRUE.
          LENGTH = MIN ( LENGTH, 30 )
          ERRMSG = 'F-INPOPT: Unexpected extra field in *OPT '//
     &             'section: '//FIELD(1:LENGTH)
          RETURN
        END IF
      END IF
C
      FAIL = .FALSE.
C
      END
