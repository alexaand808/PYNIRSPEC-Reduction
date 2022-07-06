      SUBROUTINE INPOBS ( LUNDRV, FAIL, ERRMSG )
C
C VERSION
C     13OCT13 AD Simplified, remove GETOBS argument and local GETOBS check
C     18MAR01 AD Remove GETFOV argument. 
C                    Move FOVTAN, NADPTH, TANPTH, TANHGT to INPCHK.
C                    Remove NEWLEV argument from OBSCHK
C     02JAN00 AD Also read PSIOBS.
C     18JUL98 AD Original.
C
C DESCRIPTION
C     Read user-defined observer position
C     Called by RFMINP following *OBS marker in Driver table.
C     If OBSFLG is TRUE read observer altitude [km]
C     If GRAFLG is also TRUE, read observer horizontal angle.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      LUNDRV !  I  LUN for driver table
      LOGICAL      FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  ENDCHK ! Check end of Driver Table section has been reached.
     &, NXTFLD ! Load next field from section of RFM driver file
     &, OBSCHK ! Check and set observer altitude.
     &, RFMLOG ! Write text record from RFM input file.
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'obscom.inc' ! Observer Position
C
C LOCAL CONSTANTS
      REAL PSIMAX          ! Max abs. value of observer horiz.angle 
        PARAMETER ( PSIMAX = 90.0 ) ! 90 deg from reference profile
C
C LOCAL VARIABLES
      INTEGER       IOS     !  Saved value of IOSTAT for error message
      INTEGER       LENGTH  !  No.of characters written in FIELD 
      REAL          ALTTST  !  Observer Altitude read from Driver table
      CHARACTER*80  FIELD   !  Data field from RECORD
      CHARACTER*132 MESSGE  !  Message sent to log file
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      CALL NXTFLD ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
      IF ( LENGTH .EQ. 0 ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-INPOBS: No data in *OBS section of Driver Table'
        RETURN
      END IF
C
      READ ( FIELD, *, IOSTAT=IOS ) ALTTST
      IF ( IOS .NE. 0 ) THEN
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,I11)' ) 
     &    'F-INPOBS: Error reading Observer Altitude.'
     &    //' IOSTAT=', IOS
        RETURN
      END IF 
C
      WRITE ( MESSGE, '(A,G12.5,A)' )
     &    'I-INPOBS: Read Observer Altitude =', ALTTST, ' [km]'
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
C Check next field (should only be PSIOBS - Observer Horizontal angle).
      IF ( GRAFLG ) THEN
        CALL NXTFLD ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
C
        IF ( LENGTH .EQ. 0 ) THEN
          FAIL = .TRUE.
          ERRMSG = 
     &      'F-INPOBS: No observer horizontal angle supplied '//
     &      '(required with GRA flag)'
          RETURN
        ENDIF
C
        READ ( FIELD, *, IOSTAT=IOS ) PSIOBS
        IF ( IOS .NE. 0 ) THEN
          FAIL = .TRUE.
          WRITE ( ERRMSG, '(A,I11)' ) 
     &      'F-INPOBS: Error reading Obs.horiz.angle, IOSTAT=', IOS
        ELSE IF ( ABS ( PSIOBS ) .GT. PSIMAX ) THEN
          FAIL = .TRUE.
          WRITE ( ERRMSG, '(A,G12.5,A,G12.5,A)' )
     &      'F-INPOBS: Abs(Obs.Horiz.angle)=', ABS(PSIOBS),
     &      ' > PSIMAX=', PSIMAX, ' [deg]'
        ELSE 
          WRITE ( MESSGE, '(A,G12.5,A)' ) 
     &      'I-INPOBS: Read Observer Horiz.Angle =',PSIOBS,' [deg]'
          CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        END IF
        IF ( FAIL ) RETURN
C
        IF ( PSIOBS .LT. 0.0 .OR. PSIOBS .GT. 45.0 ) THEN
          MESSGE = 'W-INPOBS: Unusual value of Obs.Horiz.Angle'
          CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        END IF
C
      END IF
C
C Check valid altitude and insert extra profile level in atmosphere if reqd.
C
      CALL OBSCHK ( ALTTST, FAIL, ERRMSG )
      IF ( FAIL ) RETURN

C Check for any superfluous fields in section         
      CALL ENDCHK ( LUNDRV, 'OBS', FAIL, ERRMSG )
C
      END
