      SUBROUTINE ATMPSI ( NAMATM, PSI, FAIL, ERRMSG )
C
C VERSION
C     03-JAN-08  AD  Original.
C
C DESCRIPTION
C     Extract Psi angle from brackets following .atm filename
C     Called by FILATM for each file in *ATM section with '(' character
C     Also strips off '(...)' from NAMATM leaving just the filename.
C
      IMPLICIT NONE
C
C ARGUMENTS
      CHARACTER*(*) NAMATM    ! I/O Name of file 
      REAL          PSI       !  O  Horizontal angle
      LOGICAL       FAIL      !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG    !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  RFMLOG ! Write text message to RFM log file.
C
C LOCAL CONSTANTS
      REAL PSIMAX                   ! Max magnitude of horizontal angle [deg]
        PARAMETER ( PSIMAX = 90.0 ) ! STRGRA in MAKNAM limits PSI to <100deg.
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
C
C LOCAL VARIABLES
      INTEGER      IEND    ! Location of end of text field in RECORD 
      INTEGER      IOS     ! Saved value of IOSTAT for error message
      INTEGER      ISTA    ! Location of start of text field in RECORD 
      CHARACTER*9  ERRSTR  ! Piece of error message
      CHARACTER*80 MESSGE  ! Text message for LOG file
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C This routine should only be called if a '(' has been found in NAMATM
      ISTA = INDEX (  NAMATM, '(' )
      IF ( ISTA .EQ. 0 ) STOP 'F-ATMPSI: Logical Error'
C
      IF ( .NOT. GRAFLG ) THEN
        MESSGE = 'W-ATMPSI: GRA flag disabled so '//
     &             'profile Horiz.angles ignored.'
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        PSI = 0.0               ! Ensure that PSI still has a defined value
      ELSE
        IEND = INDEX ( NAMATM, ')' )
        IF ( IEND .LE. ISTA ) THEN
          ERRSTR = NAMATM(ISTA:MIN(LEN(NAMATM),ISTA+8))
          ERRMSG= 'F-ATMPSI: No '')'' character '//
     &            'following ''('' for hor.angle of .atm file: '
     &            //ERRSTR
          FAIL = .TRUE.
          RETURN
        END IF
        READ ( NAMATM(ISTA+1:IEND-1), *, IOSTAT=IOS ) PSI
        IF ( IOS .NE. 0 ) THEN
          ERRSTR = NAMATM(ISTA:MIN(ISTA+8,IEND))
          WRITE ( ERRMSG, '(A,I11)' ) 
     &      'F-ATMPSI: Unable to read horiz.angle from '''//
     &      ERRSTR//'''. IOSTAT=', IOS
          FAIL = .TRUE.
          RETURN
        ELSE IF ( ABS(PSI) .GT. PSIMAX ) THEN
          WRITE ( ERRMSG, '(A,F7.3,A)' )
     &      'F-ATMPSI: |Psi angle| > allowed maximum magnitude=',
     &       PSIMAX, ' [deg]'
          FAIL = .TRUE.
          RETURN
        ELSE
          WRITE ( MESSGE, '(A,F7.3)' ) 
     &      'I-ATMPSI: Loading horizontal gradient profiles at PSI=', 
     &      PSI
          CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        END IF
      END IF
      NAMATM = NAMATM(1:ISTA-1)
C
      END
