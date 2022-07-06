      SUBROUTINE TXTFLD ( RECLEN, RECORD, NFIELD, ISTA, IEND )
C
C VERSION
C     14-JAN-03  AD  Add test for <CR> (Linux compatible)  
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Identify start and end points of text field in record.
C     Given a character string containing a number of alpha-numeric fields
C     separated by commas or spaces, this function returns the location of the
C     first and last characters of the specified field. 
C     Returns a value 0 if reqd field# is not found before the end of the 
C     string. End value limited to record length.
C
C     The following characters are regarded as spaces:
C       Quote Marks ['] (allows text strings to be optionally enclosed in ' ')
C       Commas [,]      (to conform to free format rules)
C       Tabs ( ICHAR( ) = 9 )
C       <CR> ( ICHAR( ) = 13 )
C     Other punctuation: treated as characters.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER             RECLEN  !  I  Length of record to be searched
      CHARACTER*(*)       RECORD  !  I  ASCII Record to be searched
      INTEGER             NFIELD  !  I  Field# required
      INTEGER             ISTA    !  O  Start location of field
      INTEGER             IEND    !  O  End location of field
C
C LOCAL VARIABLES
      INTEGER        I,J     !  Pointers within record
      INTEGER        IFIELD  !  Field# identified so far
      INTEGER        LEN     !  Effective record length up to ! mark
      LOGICAL        LSPACE  !  TRUE if pointer in space between fields
      LOGICAL        LFIELD  !  TRUE if pointer within a field
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Initialise at start of each call
C
      IFIELD = 0
      LFIELD = .FALSE.
C
C Check for any exclamation marks and truncate effective record length
C
      LEN = INDEX ( RECORD, '!' ) - 1
      IF ( LEN .EQ. -1 ) LEN = RECLEN 
C
C Now identify text fields
C
      DO I = 1, LEN
        LSPACE = ( RECORD(I:I) .EQ. ' '  .OR.            ! = space character
     &             RECORD(I:I) .EQ. ','  .OR.            ! = comma
     &             RECORD(I:I) .EQ. '''' .OR.            ! = single quote mark
     &       ICHAR(RECORD(I:I)) .EQ. 9   .OR.            ! = tab
     &       ICHAR(RECORD(I:I)) .EQ. 13  )               ! = <CR>
        IF ( .NOT. LFIELD .AND. .NOT. LSPACE ) THEN      ! found start of field
          IFIELD = IFIELD + 1
          IF ( IFIELD .EQ. NFIELD ) THEN                 ! start of reqd field#
            ISTA = I
            DO J = I+1, LEN
              LSPACE = ( RECORD(J:J) .EQ. ' '  .OR.      ! = space character
     &                   RECORD(J:J) .EQ. ','  .OR.      ! = comma
     &                   RECORD(J:J) .EQ. '''' .OR.      ! = single quote mark
     &             ICHAR(RECORD(J:J)) .EQ. 9   .OR.      ! = tab
     &             ICHAR(RECORD(J:J)) .EQ. 13 )          ! = <CR>
              IF ( LSPACE ) THEN
                IEND = J - 1
                RETURN
              END IF
            END DO
            IEND = LEN
            RETURN
          END IF
        END IF
        LFIELD = .NOT. LSPACE 
      END DO
C
C If the program gets this far, the reqd field# was not found in the record
C
      ISTA = 0
      IEND = 0
      END
