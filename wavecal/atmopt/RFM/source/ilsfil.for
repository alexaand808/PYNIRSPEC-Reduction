      SUBROUTINE ILSFIL ( LUNILS, NAMILS, FAIL, ERRMSG )
C
C VERSION
C     25-APR-02  AD  Allow multiple ILS functions within single file  
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Open, read ILS data file and close.
C     Called by INPILS for each ILS file listed in *ILS section of Driver Table
C
      IMPLICIT NONE 
C
C ARGUMENTS
      INTEGER       LUNILS  !  I  LUN for ILS data file
      CHARACTER*(*) NAMILS  !  I  Name of ILS data file
      LOGICAL       FAIL    !  O  T=A fatal error was detected
      CHARACTER*80  ERRMSG  !  O  Error message written if FAIL is TRUE
C      
      EXTERNAL
     &  FILILS ! Check ILS Function and load /ILSCOM/
     &, NXTREC ! Load next record from RFM input file.
     &, OPNFIL ! Open ASCII input file, and skip/Log any initial comments.
     &, RFMLOG ! Write text message to RFM log file.
     &, TXTFLD ! Identify start and end points of text field in record
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
C
C LOCAL VARIABLES      
      INTEGER          IDUMMY  ! Pointer to end of field in RECORD
      INTEGER          IOS     ! Saved value of IOSTAT for error message
      INTEGER          ISTA    ! Pointer to start of field in RECORD
      LOGICAL          ENDSEC  ! T=reached end of file, F=continue
      CHARACTER*80     MESSGE  ! Text message for Log file
      CHARACTER*80     RECORD  ! Record read from ILS file
      DOUBLE PRECISION WNOLOW  ! Lower Wno for application of ILS function
      DOUBLE PRECISION WNOUPP  ! Upper Wno for application of ILS function
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      CALL OPNFIL ( LUNILS, NAMILS, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
      ENDSEC = .FALSE.
      DO WHILE ( .NOT. ENDSEC )
        CALL NXTREC ( LUNILS, RECORD, ENDSEC, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
C      
        IF ( .NOT. ENDSEC ) THEN
          CALL TXTFLD ( 80, RECORD, 4, ISTA, IDUMMY )
          IF ( ISTA .NE. 0 ) THEN       ! Wavenumber range applicable
            READ ( RECORD(ISTA:), *, IOSTAT=IOS, ERR=900 ) 
     &        WNOLOW, WNOUPP
          ELSE
            WNOLOW = 0.D0
            WNOUPP = 0.D0
            MESSGE = 
     &        'I-ILSFIL: Loading general ILS Fn for any Wno.range'
            CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
            IF ( FAIL ) RETURN
          END IF
          BACKSPACE ( LUNILS, IOSTAT=IOS, ERR=900 )
          CALL FILILS ( LUNILS, WNOLOW, WNOUPP, FAIL, ERRMSG )
         IF ( FAIL ) RETURN
        END IF
      END DO
C
  200 CLOSE ( LUNILS, IOSTAT=IOS, ERR=900 )
      FAIL = .FALSE.
      RETURN
C
  900 FAIL = .TRUE.
      WRITE ( ERRMSG, '(A,I11)' ) 
     &  'F-ILSFIL: I/O failure on ILS file. IOSTAT=', IOS
      END
