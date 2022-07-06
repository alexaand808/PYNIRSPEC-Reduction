      SUBROUTINE RFMLOG ( MESSGE, FAIL, ERRMSG )
C
C VERSION
C     18-AUG-11  AD  Fix bug with writing continuous text string > 80 chars.
C     14-JAN-03  AD  Limit to IEND to avoid writing <CR> (Linux compatible)
C     20-FEB-01  AD  Write out C*80 messages rather than C*79
C     06-JAN-99  AD  Ensure this copes with completely blank MESSGE
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Write text message to RFM log file.
C     Limits output to 80 characters, uses ... if continuation required
C
      IMPLICIT NONE
C
C ARGUMENTS
      CHARACTER*(*) MESSGE  !  I  Message to be written
      LOGICAL       FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG  !  O  Error message returned if FAIL is TRUE
C
C GLOBAL CONSTANTS
      INCLUDE 'logcom.inc' ! Runlog file data
C
C LOCAL VARIABLES
      INTEGER  IEND  ! Pointer to location of last character in MESSGE
      INTEGER  ILIM  ! Pointer to end of MESSGE to be written out next
      INTEGER  IOS   ! Saved value of IOSTAT for error messages
      INTEGER  ISTA  ! Pointer to start of MESSGE to be written out next
      LOGICAL  BLANK ! T=Found an blank character
      LOGICAL  MORE  ! T=more text remains to be written
      CHARACTER*80 TEXT80 ! Text records written to Log file
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Check for warning message - begins with: 'W-'
C
      IF ( MESSGE(1:2) .EQ. 'W-' ) THEN
        NWMLOG = NWMLOG + 1
      END IF
C
C First find extent of text in MESSGE string - set IEND to point to last char.
C
      BLANK = .TRUE.
      IEND = LEN ( MESSGE ) + 1
      DO WHILE ( BLANK .AND. IEND .GT. 1 )
        IEND = IEND - 1
        BLANK = MESSGE(IEND:IEND) .EQ. ' ' 
     &          .OR. ICHAR(MESSGE(IEND:IEND)) .EQ. 13      ! <CR> character
      END DO
C
      MORE = .TRUE.
      ISTA = 1
      DO WHILE ( MORE )
        TEXT80 = ' '
        IF ( IEND - ISTA + 1 .GT. 80 ) THEN
          MORE = .TRUE.
          DO ILIM = ISTA + 76, ISTA, -1
            IF ( MESSGE(ILIM:ILIM) .EQ. ' ' ) GOTO 100
          END DO
          ILIM = ISTA + 76
  100     TEXT80 = MESSGE(ISTA:ILIM)
          TEXT80(78:80) = '...'
          WRITE ( LUNLOG, '(A80)', IOSTAT=IOS ) TEXT80
          ISTA = ILIM + 1
        ELSE
          MORE = .FALSE.
          TEXT80 = MESSGE(ISTA:IEND)
          WRITE ( LUNLOG, '(A80)', IOSTAT=IOS ) TEXT80
        END IF
        IF ( IOS .NE. 0 ) THEN
          FAIL = .TRUE.
          WRITE ( ERRMSG, '(A,I11)' ) 
     &     'F-RFMLOG: Failed to write message to LOG file. IOSTAT=',IOS
          RETURN
        END IF
      END DO
C
      FAIL = .FALSE.
C
      END
