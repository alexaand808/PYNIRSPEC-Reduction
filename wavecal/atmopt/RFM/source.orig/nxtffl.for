      SUBROUTINE NXTFFL ( LUNDRV, LUNSEC, ONEREC, 
     &                    FIELD, LENGTH, FAIL, ERRMSG )
C
C VERSION
C     10OCT13 AD Original. Adapted from NXTFLD
C
C DESCRIPTION
C     Load next field from section of RFM drv, expanding filenames.
C     General purpose.
C     This version of NXTFLD will first attempt to open any field as a filename
C     and, if successful, return subsequent fields from that file
C
      IMPLICIT NONE
      SAVE
C
      EXTERNAL
     &  CHKFIL ! Check if argument represents an existing filename
     &, NXTFL2 ! Load next field from secondary input file.
     &, NXTREC ! Load next record from RFM input file.
     &, OPNFIL ! Open ASCII input file, and skip/Log any initial comments.
     &, TXTFLD ! Identify start and end points of text field in record
      LOGICAL CHKFIL
C
C ARGUMENTS
      INTEGER       LUNDRV  !  I  Logical Unit Number for Driver Table
      INTEGER       LUNSEC  !  I  LUN for secondary input file
      LOGICAL       ONEREC  !  I  T=limit input to one record from Driver Table
      CHARACTER*(*) FIELD   !  O  Field extracted from record
      INTEGER       LENGTH  !  O  Length of extracted field
      LOGICAL       FAIL    !  O  T=Fatal error detected
      CHARACTER*80  ERRMSG  !  O  Error message written if FAIL is TRUE
C
C LOCAL VARIABLES
      INTEGER       IEND    ! Pointer to end of field in record
      INTEGER       IFLD    ! Field counter within each record
      INTEGER       ISTA    ! Pointer to start of field in RECORD
      LOGICAL       ENDSEC  ! TRUE if end of section encountered
      LOGICAL       FILINP  ! TRUE if currently taking input from LUNSEC
      LOGICAL       NEWREC  ! TRUE if new record to be read
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
  100 IF ( NEWREC ) THEN 
        CALL NXTREC ( LUNDRV, RECORD, ENDSEC, FAIL, ERRMSG )
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
      IF ( FILINP ) THEN             ! Currently reading data from file
        CALL NXTFL2 ( LUNSEC, FIELD, LENGTH, FAIL, ERRMSG ) 
        IF ( FAIL ) RETURN
        IF ( LENGTH .NE. 0 ) RETURN  ! Exit with next field from file
        CLOSE ( LUNSEC )             ! Reached end-of-file, go back to lundrv
        FILINP = .FALSE.
      END IF
C
C Read next field from LUNDRV
      IFLD = IFLD + 1
      CALL TXTFLD ( 200, RECORD, IFLD, ISTA, IEND )
C
C If no more fields, either return (ONEREC) or read next record
      IF ( ISTA .EQ. 0 ) THEN        ! no (more) fields in record - get next rec
        NEWREC = .TRUE.              ! new record from LUNDRV required
        LENGTH = 0                   ! flag for no current FIELD loaded
        IF ( ONEREC ) RETURN
        GOTO 100
      END IF
C
C Identified next field in LUNDRV record
      FIELD = RECORD(ISTA:IEND)
C
C If this is a valid filename, open it, set FILINP=TRUE and back to top to
C read first field from the file
      IF ( CHKFIL ( LUNSEC, FIELD ) ) THEN
        CALL OPNFIL ( LUNSEC, FIELD, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        FILINP = .TRUE. 
        GOTO 100
      END IF    
C
C If not a filename, return FIELD and its LENGTH
      LENGTH = IEND - ISTA + 1

      END
