      SUBROUTINE ACCLOG ( TEXT, FAIL, ERRMSG )
C
C VERSION
C     18-NOV-04  AD  Return with FAIL=TRUE if size of LOGMSG exceeded
C                    And initialise FAIL=FALSE at start of code
C     01-JUL-02  AD  Original.
C
C DESCRIPTION
C     Accumulate and send text record to log file.
C     General purpose routine
C     If TEXT = ' ' the accumulate message is sent and resets.
C
      IMPLICIT NONE
      SAVE
C
C ARGUMENTS
      CHARACTER*(*) TEXT   !  I  Part of message, or ' ' to send/clear
      LOGICAL       FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  RFMLOG ! Write text message to RFM log file.
C
C LOCAL VARIABLES
      INTEGER I            ! Pointer to first non-blank character in TEXT
      INTEGER IPT          ! Pointer to start of next text section in LOGMSG
      INTEGER J            ! Pointer to last non-blank character in TEXT
      INTEGER JPT          ! Pointer to end of next text section in LOGMSG
      INTEGER L            ! Length of TEXT argument
      CHARACTER*200 LOGMSG ! Local buffer for accumulating text message
C
C DATA STATEMENTS
      DATA IPT / 1 / 
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      FAIL = .FALSE.
C
C Send/clear LOGMSG
      IF ( TEXT .EQ. ' ' ) THEN
        CALL RFMLOG ( LOGMSG(1:IPT-1), FAIL, ERRMSG )
        IPT = 1
        RETURN
      END IF
C
C Find location of first and last non-blank characters in TEXT
C Since TEXT isn't blank, there ought to be at least one character
      L = LEN ( TEXT )
      I = 1
      DO WHILE ( TEXT(I:I) .EQ. ' ' )
        I = I + 1
        IF ( I .EQ. L+1 ) STOP 'F-LOGACC: Logical error#1'
      END DO
      J = L
      DO WHILE ( TEXT(J:J) .EQ. ' ' )
        J = J - 1
        IF ( J .EQ. 0 ) STOP 'F-LOGACC: Logical error#2'
      END DO
C
C TEXT is too long to store in local buffer LOGMSG
      IF ( J - I + 1 .GT. 200 ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-LOGACC: text string too long'
        RETURN
      END IF
C
C No.characters written to LOGMSG = (J-I+1) + (1 space)  
 100  CONTINUE
      JPT = IPT + J - I + 1
C
C If LOGMSG overflows, written current contents to LOG file and reset
      IF ( JPT .GT. 200 ) THEN
        CALL RFMLOG ( LOGMSG(1:IPT-1), FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        IPT = 1
        GOTO 100
      END IF
C
      LOGMSG(IPT:JPT) = TEXT(I:J)//' '
      IPT = JPT + 1
C
      END
