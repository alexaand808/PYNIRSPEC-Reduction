      SUBROUTINE NXTREC ( LUN, RECORD, ENDSEC, FAIL, ERRMSG )
C
C VERSION
C      07-AUG-13  AD  Rewritten and Allow record size up to C*200
C      27-MAR-07  AD  Test for <CR> (CHAR(13)) to find end of record
C      03-MAR-97  AD  Version 3.
C      01-OCT-96  AD  Version 2.
C      01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Load next record from RFM input file.
C     General purpose.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER       LUN     !  I  Logical Unit Number
      CHARACTER*(*) RECORD  !  O  Last record read from file
      LOGICAL       ENDSEC  !  O  T=reached end of section
      LOGICAL       FAIL    !  O  T=Fatal error detected
      CHARACTER*80  ERRMSG  !  O  Error message written if FAIL is TRUE
C
C LOCAL VARIABLES
      CHARACTER*210 REC210  ! Extended version of RECORD
      INTEGER       IPT     ! Pointer to each character in record
      INTEGER       IEND    ! Pointer to last part before end-of-line comments 
      INTEGER       IOS     ! Saved value of IOSTAT for error message
      INTEGER       LREC    ! Length of input RECORD string
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      FAIL = .FALSE.
      ENDSEC = .FALSE.
      LREC = LEN ( RECORD )
C Check that RECORD is only meant to accommodate up to C*200
      IF ( LREC .GT. 200 ) 
     &  STOP 'F-NXTREC: unexpectedly long RECORD argument'
C
  100 READ ( LUN, '(A)', IOSTAT=IOS, ERR=900, END=200 ) REC210
      IF ( REC210(1:1) .NE. '*' ) THEN         ! Not yet reached next section
        IEND = INDEX ( REC210, '!' ) - 1

        IF ( IEND .EQ. 0 ) GOTO 100          ! First character was !
        IF ( IEND .EQ. -1 .OR. IEND .GT. LREC ) THEN
          DO IPT = LREC+1, MIN(210,LREC+20)  ! Check for extra non-comment chars
            IF ( REC210(IPT:IPT) .NE. ' ' .AND. 
     &           REC210(IPT:IPT) .NE. '!'       ) THEN
              FAIL = .TRUE.
              ERRMSG = 'F-NXTREC: Extra characters beyond expected '//
     &                 'end of record'//REC210(IPT:MIN(210,IPT+20))
              RETURN
            END IF
          END DO
          IEND = LREC
        END IF
        RECORD = REC210(1:IEND)  
        DO IPT = 1, IEND
          IF ( RECORD(IPT:IPT) .NE. ' ' ) RETURN ! Found non-blank character
        END DO
        GOTO 100          ! No non-blank characters, so discard and get next
      END IF
C
C Reached end of section
      BACKSPACE ( LUN, IOSTAT=IOS, ERR=900 ) ! NB: No backspace at end-of-file
  200 ENDSEC = .TRUE.
      RETURN
C
  900 FAIL = .TRUE.
      WRITE ( ERRMSG, '(A,I11)' )
     &  'F-NXTREC: Failed to read next record of file. IOSTAT=', IOS
      END
