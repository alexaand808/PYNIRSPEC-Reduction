      SUBROUTINE INPNTE ( LUNDRV, LUNNTE, FAIL, ERRMSG )
C
C VERSION
C     16OCT13 AD Simplified. Remove GETNTE argument and local GOTNTE check.
C     09-JAN-08  AD  Don't initialise NNTE here - do it in RFMDEF
C     23-JUL-03  AD  Add NTEDEF
C     14-JUL-98  AD  Replace NXTREC by NXTFLD
C     04-DEC-97  AD  Add GOTNTE to check for duplicate section
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     20-SEP-96  AD  Original.
C
C DESCRIPTION
C     Read non-LTE files from driver table.
C     Called by RFMINP once if *NTE section found in Driver Table and NTEFLG=T
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      LUNDRV  !  I  LUN for Driver File
      INTEGER      LUNNTE  !  I  Temporary LUN for non-LTE data
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  NTEDEF ! Use default .nte  filename to find any missing files.
     &, NTEFIL ! Read file containing non-LTE data
     &, NXTFLD ! Load next field from section of RFM driver file
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C LOCAL VARIABLES
      INTEGER       LENGTH ! Length of filename read from driver file
      CHARACTER*80  NAMDEF ! .nte filename template
      CHARACTER*80  NAMNTE ! Name of NTE file read from Driver Table
      LOGICAL       GOTDEF ! True if filename template loaded
C
C DATA STATEMENTS
      DATA GOTDEF / .FALSE. /
      SAVE GOTDEF
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      CALL NXTFLD ( LUNDRV, NAMNTE, LENGTH, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C Jump here for each new field in *NTE section. 
      DO WHILE ( LENGTH .NE. 0 )
        IF ( INDEX ( NAMNTE(1:LENGTH), '*' ) .NE. 0 ) THEN
          IF ( GOTDEF ) THEN
            ERRMSG = 'F-INPNTE: *NTE section contains more than '//
     &                 'one filename template'
            FAIL = .TRUE.
          ELSE
            GOTDEF = .TRUE.
            NAMDEF = NAMNTE(1:LENGTH)
          ENDIF
        ELSE
          CALL NTEFIL ( LUNNTE, NAMNTE(1:LENGTH), FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        END IF
        CALL NXTFLD ( LUNDRV, NAMNTE, LENGTH, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END DO
C
C Load any remaining files using filename template
      IF ( GOTDEF ) CALL NTEDEF ( LUNNTE, NAMDEF, FAIL, ERRMSG )
C
      END
