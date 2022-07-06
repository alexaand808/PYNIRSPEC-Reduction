      SUBROUTINE INPCRV ( LUNDRV, FAIL, ERRMSG )
C
C VERSION
C     16OCT13 AD Simplify: Remove GETCRV argument and local GOTCRV. Add ENDCHK
C     27APR00 AD Convert RADMIN, RADMAX from SP to DP
C     14JUL98 AD Replace NXTREC with NXTFLD
C     04DEC97 AD Add GOTCRV to check for duplicate section
C     03MAR97 AD Version 3.
C     01OCT96 AD Version 2.
C     18SEP96 AD Original.
C
C DESCRIPTION
C     Read user-defined local radius of curvature.
C     Called once by RFMINP following *CRV marker in Driver table if CRVFLG = T
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
     &, RFMLOG ! Write text message to RFM log file
C
C COMMON VARIABLES
      INCLUDE 'crvcom.inc' ! Local radius of curvature
C
C LOCAL CONSTANTS
      DOUBLE PRECISION RADMIN  ! Minimum sensible radius for any object
        PARAMETER ( RADMIN = 2000.0D0 )   ! cf Mars=3500km ? Titan=?
      DOUBLE PRECISION  RADMAX  ! Maximum sensible radius for any object
        PARAMETER ( RADMAX = 100000.0D0 ) ! cf Jupiter=80000km ?
C
C LOCAL VARIABLES
      INTEGER       IOS     !  Saved value of IOSTAT for error message
      INTEGER       LENGTH  !  Length of Field
      CHARACTER*80  FIELD   !  Field read from driver table
      CHARACTER*132 MESSGE  !  Message sent to log file
C
C EXECUTABLE CODE --------------------------------------------------------------
C
C Read first field in *CRV section
      CALL NXTFLD ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      READ ( FIELD, *, IOSTAT=IOS ) RADCRV
      IF ( IOS .NE. 0 ) THEN
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,I11)' ) 
     &    'F-INPCRV: Error reading Radius of Curvature value.'
     &    //' IOSTAT=', IOS
        RETURN
      ELSE 
        WRITE ( MESSGE, '(A,G10.3)' )
     &    'I-INPCRV: Using local radius of curvature=', RADCRV
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        IF ( RADCRV .LT. RADMIN .OR. RADCRV .GT. RADMAX ) THEN
          FAIL = .TRUE.
          WRITE ( ERRMSG, '(A,F8.0,A,F8.0)' )
     &      'F-INPCRV: Radius of Curvature outside allowed range: ',
     &      RADMIN, ' to ', RADMAX
          RETURN
        END IF
      END IF
C 
C Shouldn't be any more fields in *CRV section
      CALL ENDCHK ( LUNDRV, 'CRV', FAIL, ERRMSG )
C
      END
