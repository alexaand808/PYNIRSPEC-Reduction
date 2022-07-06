      SUBROUTINE INPJAC ( LUNDRV, LUNJAC, FAIL, ERRMSG )
C
C VERSION
C     10OCT13 AD Bug#83: Rewritten. Remove GETJAC argument and local GOTJAC
C     22AUG13 AD Bug#75: Reformat error message for unexpected extra field
C     25JUN07 AD Bug#64: Increase TXTFLD from C*8 to C*11 max
C     01JAN04 AD Increase size of TARGET array to allow for isotopes
C     03SEP03 AD Allow for 'SFCTEM' and 'SFCEMS' jacobians
C     18JAN00 AD Use CDUMMY rather than ' ' as argument to JACCHK
C     04FEB99 AD Modified definitions of perturbation levels
C     02JAN99 AD Original.
C
C DESCRIPTION
C     Read details of Jacobians to be calculated
C     Called by RFMINP once if JACFLG set TRUE.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      LUNDRV  !  I  LUN for Driver File
      INTEGER      LUNJAC  !  I  Spare LUN for opening Jac.file of altitudes
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  ACCLOG ! Accumulate and send text record to log file.
     &, JACCHK ! Check target & Alt.range of Jacobian data
     &, NXTFFL ! Load next field from section of RFM drv, expanding filenames
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'jaccom.inc' ! Jacobian data
C
C LOCAL VARIABLES
      INTEGER      IOS     ! Saved value of IOSTAT for error messages
      INTEGER      LENTGT  ! Length of filled TARGET string
      INTEGER      LENGTH  ! Length of field from driver table or jac file
      REAL         ALTRTV  ! Retrieval altitude [km] for Jacobians
      CHARACTER*1  CDUMMY  ! Dummy character ' '
      CHARACTER*80 FIELD   ! Field extracted from RECORD
      CHARACTER*11 TARGET  ! Retrieval species for Jacobian
      LOGICAL      SFCPTB  ! T=surface perturbation, F=atmospheric perturbation
C DATA STATEMENTS      
      DATA CDUMMY / ' ' /         ! Shouldn't get altered by JACCHK
      SAVE CDUMMY
C
C EXECUTABLE CODE ------------------------------------------------------------
C
C Each record should start with the name of a species (with profile already 
C loaded), or 'TEM', or 'PRE'
      NJAC = 0
      CALL NXTFFL ( LUNDRV, LUNJAC, .TRUE., 
     &              TARGET, LENTGT, FAIL, ERRMSG ) 
      IF ( FAIL ) RETURN
      DO WHILE ( LENTGT .GT. 0 ) 
C All recognised TARGET values have maximum 11 characters
        IF ( LENTGT .GT. 11 ) THEN
          FAIL = .TRUE.
          LENTGT = MIN ( LENTGT, 30 )               ! Restrict for error message
          ERRMSG = 'F-INPJAC: Unrecognised TARGET='//TARGET(1:LENTGT)
          RETURN
        END IF
C
C Check species (TARGET is converted to lower case within JACCHK)
        CALL JACCHK ( TARGET, .FALSE., 0.0, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        SFCPTB = TARGET .EQ. 'sfctem' .OR. TARGET .EQ. 'sfcems'
        CALL ACCLOG ( 'I-INPJAC: '//TARGET(1:LENTGT), FAIL, ERRMSG )
        IF ( FAIL ) RETURN
C Load first altitude field (if any) from same record
        CALL NXTFFL ( LUNDRV, LUNJAC, .TRUE., 
     &                FIELD, LENGTH, FAIL, ERRMSG ) 
C If no altitude fields in record, column ptb or homog.path
        IF ( LENGTH .EQ. 0 ) THEN      
          CALL JACCHK ( CDUMMY, .TRUE., 0.0, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        ELSE IF ( HOMFLG .OR. SFCPTB ) THEN        ! No altitude info required 
          LENGTH = MIN ( LENGTH, 8 )             ! Restrict for error message
          ERRMSG = 'F-INPJAC: Unexpected extra field'//
     &      ' in *JAC record for '//TARGET(1:LENTGT)//
     &      ', text='//FIELD(1:LENGTH)
          FAIL = .TRUE.
          RETURN
        END IF
C At this point first altitude is loaded into FIELD
        DO WHILE ( LENGTH .GT. 0 )          ! Repeat for subsequent alt fields
          READ ( FIELD, *, IOSTAT=IOS ) ALTRTV
          IF ( IOS .NE. 0 ) THEN
            FAIL = .TRUE.
            LENGTH = MIN ( LENGTH, 15 )     ! Restrict for error message
            WRITE ( ERRMSG, '(A,A,A,A)' ) 
     &        'F-INPJAC: Error reading *JAC altitude for ', 
     &        TARGET(1:LENTGT), ' from field=', FIELD(1:LENGTH)
            RETURN
          END IF
          CALL JACCHK ( CDUMMY, .FALSE., ALTRTV, FAIL, ERRMSG )
            IF ( FAIL ) RETURN
          CALL ACCLOG ( FIELD, FAIL, ERRMSG )
            IF ( FAIL ) RETURN
          CALL NXTFFL ( LUNDRV, LUNJAC, .TRUE., 
     &                  FIELD, LENGTH, FAIL, ERRMSG ) 
        END DO
        CALL JACCHK ( ' ', .TRUE., ALTRTV, FAIL, ERRMSG ) 
          IF ( FAIL ) RETURN
        CALL ACCLOG ( ' ', FAIL, ERRMSG )
          IF ( FAIL ) RETURN
C Read first field in next record from *JAC section
        CALL NXTFFL ( LUNDRV, LUNJAC, .TRUE., 
     &              TARGET, LENTGT, FAIL, ERRMSG ) 
        IF ( FAIL ) RETURN
      END DO
C
      IF ( NJAC .EQ. 0 ) THEN
        FAIL = .TRUE.
        ERRMSG = 
     &    'F-INPJAC: No data found in *JAC section of driver table'
        RETURN
      END IF
C
      FAIL = .FALSE.
      END
