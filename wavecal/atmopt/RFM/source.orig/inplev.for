      SUBROUTINE INPLEV ( LUNDRV, LUNLEV, FAIL, ERRMSG )
C
C VERSION
C     16OCT13 AD Simplified. Remove GETLEV argument and local GOTLEV test
C     10OCT13 AD Bug#83: Rewritten
C     01JAN04 AD Original.
C
C DESCRIPTION
C     Read altitudes for intermediate level outputs.
C     Called by RFMINP once if *LEV section in Driver Table and LEVFLG set TRUE.
C     These are temporarily stored in *JAC arrays to use the fact that these
C     arrays are automatically updated if extra atmos. levels are inserted
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      LUNDRV  !  I  LUN for Driver File
      INTEGER      LUNLEV  !  I  Spare LUN for opening LEV.file of altitudes
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  ACCLOG ! Accumulate and send text record to log file.
     &, LEVCHK ! Check range of level data
     &, NXTFFL ! Load next field from section of RFM drv, expanding filenames
C
C GLOBAL VARIABLES
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C      
C COMMON VARIABLES
      INCLUDE 'jaccom.inc' ! Jacobian data
C
C LOCAL VARIABLES
      INTEGER      IOS     ! Saved value of IOSTAT for error messages
      INTEGER      LENGTH  ! No.characters in FIELD
      REAL         ALTLEV  ! Altitude level [km] for output
      CHARACTER*80 FIELD   ! Field extracted from RECORD
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Initialise *JAC variables used to store information
      NJAC = 0
C
C Read first field in *LEV section
      CALL ACCLOG ( 'I-INPLEV: ', FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      CALL NXTFFL ( LUNDRV, LUNLEV, .FALSE., 
     &              FIELD, LENGTH, FAIL, ERRMSG )
C From here, check through each field in section
      DO WHILE ( LENGTH .NE. 0 )
        READ ( FIELD, *, IOSTAT=IOS ) ALTLEV
        IF ( IOS .NE. 0 ) THEN    
          FAIL = .TRUE.
          LENGTH = MIN ( LENGTH, 30 )     ! Restrict for error message
          WRITE ( ERRMSG, '(A,A)' ) 
     &      'F-INPLEV: Error reading *LEV altitude from field=', 
     &      FIELD(1:LENGTH)
          RETURN
        END IF
        CALL LEVCHK ( ALTLEV, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        CALL ACCLOG ( FIELD, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        CALL NXTFFL ( LUNDRV, LUNLEV, .FALSE., 
     &                FIELD, LENGTH, FAIL, ERRMSG )
      END DO
      CALL ACCLOG ( ' ', FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
      IF ( NJAC .EQ. 0 ) THEN
        FAIL = .TRUE.
        ERRMSG = 
     &    'F-INPLEV: No data found in *LEV section of driver table'
        RETURN
      END IF
C
      END
