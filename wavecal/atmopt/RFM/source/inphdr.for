      SUBROUTINE INPHDR ( LUNDRV, FAIL, ERRMSG )
C
C VERSION
C     07-AUG-13  AD  Allow for records up to C*200
C     30-APR-03  AD  Remove <CR> character from header (Linux compatible) 
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Read RFM Driver Header record
C     Called once by RFMINP.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      LUNDRV  !  I  LUN for Driver File
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  NXTREC ! Load next record from RFM input file.
     &, RFMLOG ! Write text message to RFM log file
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'outcom.inc' ! RFM output file data
C
C LOCAL VARIABLES
      INTEGER       ICOM    ! Index of any comment '!' marker in record
      INTEGER       IEND    ! Last non-blank character in header
      LOGICAL       ENDSEC  ! T=end of section detected
      CHARACTER*200 RECORD  ! Record read from Driver Table
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C First non-blank record should be text for Header
C
      CALL NXTREC ( LUNDRV, RECORD, ENDSEC, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      IF ( ENDSEC ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-INPHDR: No header text found in *HDR section'
        RETURN
      END IF
C
      ICOM = INDEX ( RECORD, '!' )
      IF ( ICOM .NE. 0 ) RECORD = RECORD(1:ICOM)
C
C Find last printable character in record
      IEND = LEN ( RECORD )
      DO WHILE ( IEND .GE. 1 .AND. 
     &             ( RECORD(IEND:IEND) .EQ. ' ' .OR.          ! space
     &               ICHAR ( RECORD(IEND:IEND) ) .EQ. 13 ) )  ! <CR> character
        IEND = IEND - 1
      END DO
      IF ( IEND .EQ. 0 ) THEN
        TXTHDR = ' '
      ELSE
        TXTHDR = RECORD(1:IEND)
      END IF
      CALL RFMLOG ( '    '//TXTHDR, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
C Shouldn't be any more non-blank records in *HDR section
C
      CALL NXTREC ( LUNDRV, RECORD, ENDSEC, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
      IF ( .NOT. ENDSEC ) THEN               ! Found something else unexpected 
        FAIL = .TRUE.
        ERRMSG = 
     &    'F-INPHDR: Unexpected extra record in *HDR section: '//
     &    RECORD(1:25)//'...'        
        RETURN
      END IF
C
      END
