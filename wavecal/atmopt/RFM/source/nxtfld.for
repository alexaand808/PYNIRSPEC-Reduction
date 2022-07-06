      SUBROUTINE NXTFLD ( LUN, FIELD, LENGTH, FAIL, ERRMSG )
C
C VERSION
C     10OCT13 AD Add EXTERNAL declarations
C     09DEC11 AD Change RECORD from C*120 to C*200
C     17AUG11 AD Change RECORD from C*80 to C*120
C     22JAN98 AD (Standard F77) Put DATA statements *after* declarations
C     03MAR97 AD Version 3.
C     01OCT96 AD Version 2.
C     01SEP96 AD Version 1.
C
C DESCRIPTION
C     Load next field from section of RFM driver file.
C     General purpose.
C     If NEWREC = TRUE
C       Read next record
C       Read first field of record
C     Else
C       Read next field of current record
C     End if
C     If next record begins with section marker ('*'), backspaces so record
C     can be read again and sets LENGTH=0 for output, and NEWREC=TRUE for next
C     call.
C
      IMPLICIT NONE
      SAVE
C
      EXTERNAL
     &  NXTREC ! Load next record from RFM input file.
     &, TXTFLD ! Identify start and end points of text field in record
C
C ARGUMENTS
      INTEGER       LUN     !  I  Logical Unit Number
      CHARACTER*(*) FIELD   !  O  Field extracted from record
      INTEGER       LENGTH  !  O  Length of extracted field
      LOGICAL       FAIL    !  O  T=Fatal error detected
      CHARACTER*80  ERRMSG  !  O  Error message written if FAIL is TRUE
C
C LOCAL VARIABLES
      INTEGER      IEND    ! Pointer to end of field in record
      INTEGER      IFLD    ! Field counter within each record
      INTEGER      ISTA    ! Pointer to start of field in RECORD
      LOGICAL      ENDSEC  ! TRUE if end of section encountered
      LOGICAL      NEWREC  ! TRUE if new record to be read
      CHARACTER*200 RECORD  ! Record from Driver Table
C
C DATA STATEMENTS
              DATA NEWREC / .TRUE. /
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      FAIL = .FALSE.
      ENDSEC = .FALSE.
C
C Read next record from file if required
C
  100 IF ( NEWREC ) THEN 
        CALL NXTREC ( LUN, RECORD, ENDSEC, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        IF ( ENDSEC ) THEN                              ! Reached end of section
          LENGTH = 0
          RETURN
        END IF
        NEWREC = .FALSE.
        IFLD = 0
      END IF
C
C Read next field
C
      IFLD = IFLD + 1
      CALL TXTFLD ( 200, RECORD, IFLD, ISTA, IEND )
      IF ( ISTA .EQ. 0 ) THEN        ! no (more) fields in record - get next rec
        NEWREC = .TRUE.
        GOTO 100
      END IF
C
      FIELD = RECORD(ISTA:IEND)
      LENGTH = IEND - ISTA + 1

      END
