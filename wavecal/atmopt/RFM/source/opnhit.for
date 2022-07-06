      SUBROUTINE OPNHIT ( LUNHIT, NAMHIT, FAIL, ERRMSG )
C
C VERSION
C     13-AUG-13  AD  Original. Formerly part of hitfil.for
C
C DESCRIPTION
C     Open HITRAN binary line data file 
C     Called once by HITFIL if *HIT section of driver table is present.
C     This first attempts to open the file using the RECLEN parameter from
C     reclen.inc as the value for RECL, if this fails then retries with the
C     alternative value (88 or 22).
C     'fail' is determined by the first data record of the file, which should
C     be a forward pointer block, not having the initial value of -7
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER       LUNHIT  ! I/O LUN for HIT File/next free LUN
      CHARACTER*(*) NAMHIT  !  I  Name of HIT file
      LOGICAL       FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  RFMLOG ! Write text message to RFM log file
C
C GLOBAL CONSTANTS
      INCLUDE 'reclen.inc' ! Record Length for HITRAN binary data files
C
C LOCAL VARIABLES
      INTEGER       IFORM       ! binary file format ID
      INTEGER       IOS         ! Saved value of IOSTAT for error messages
      INTEGER       IREC1       ! First data record# in file
      INTEGER       LSTAT       ! Used to identify HITRAN Fwd.Pointer block
      INTEGER       RECL2       ! Alternative value of RECL to be tried
      CHARACTER*80  MESSGE      ! Text message sent to LOG file 
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Send message to LOG file saying which file is about to be opened
      MESSGE = 'I-OPNHIT: Opening HITRAN File: '//NAMHIT
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
C and the RECLEN value being used
      WRITE ( MESSGE, '(A,I2)' ) 
     &  'I-OPNHIT: Open using RECL (from reclen.inc)=', RECLEN
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
      OPEN ( UNIT=LUNHIT, FILE=NAMHIT, STATUS='OLD', 
     &       ACCESS='DIRECT', RECL=RECLEN, IOSTAT=IOS, ERR=910 )
C
C Read header and extract first/last record#
C It seems different RECL parameters both give the same values here
C but just in case, assume any read error is due to incorrect RECL
      READ ( LUNHIT, REC=1, IOSTAT=IOS, ERR=100 )
     &  LSTAT, IFORM, IREC1
C
C Read first data record, should have LSTAT=-7 since it is a fwd pointer block
C Again it seems both RECL values give IOSTAT=0 but just in case, jump to 100
C in case any read error is due to incorrect RECL. However only the correct
C RECL value will lead to LSTAT being correctly read as -7
      READ ( LUNHIT, REC=IREC1, IOSTAT=IOS, ERR=100 ) LSTAT
C
      IF ( LSTAT .EQ. -7 ) THEN
        MESSGE = 'I-OPNHIT: File opened successfully'
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        RETURN         ! exit with file opened successfully using orig RECLEN
      END IF
C
C Jump here if problems with original RECL value
  100 CONTINUE
      CLOSE ( LUNHIT ) 
      IF ( RECLEN .EQ. 22 ) THEN
        RECL2 = 88
      ELSE
        RECL2 = 22 
      END IF    

C and the RECLEN value being used
      WRITE ( MESSGE, '(A,I2)' ) 
     &  'W-OPNHIT: Read error. Try reopening file using RECL=', RECL2
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
      OPEN ( UNIT=LUNHIT, FILE=NAMHIT, STATUS='OLD', 
     &       ACCESS='DIRECT', RECL=RECL2, IOSTAT=IOS, ERR=910 )
C
C Read header and extract first/last record#
C Assume any error this time is fatal
      READ ( LUNHIT, REC=1, IOSTAT=IOS, ERR=900 )
     &  LSTAT, IFORM, IREC1
C
C Read first data record, should have LSTAT=-7 since it is a fwd pointer block
      READ ( LUNHIT, REC=IREC1, IOSTAT=IOS, ERR=900 ) LSTAT, iform
C
      IF ( LSTAT .EQ. -7 ) THEN
        MESSGE = 'I-OPNHIT: File opened successfully'
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
      ELSE                              ! no explanation for this
        WRITE ( ERRMSG, '(A,I11)' ) 
     &    'F-OPNHIT: Expected LSTAT=-7 but found LSTAT=', LSTAT
        FAIL = .TRUE.
      END IF
      RETURN   ! exit with alternative RECLEN
C
  900 WRITE ( ERRMSG, '(A,I11)') 
     &  'F-OPNHIT: Read failure on HITRAN file. IOSTAT=', IOS
      FAIL = .TRUE.
      RETURN
C
  910 WRITE ( ERRMSG, '(A,I11)') 
     &  'F-OPNHIT: Open failure on HITRAN file. IOSTAT=', IOS
      FAIL = .TRUE.
      RETURN
C
      END
