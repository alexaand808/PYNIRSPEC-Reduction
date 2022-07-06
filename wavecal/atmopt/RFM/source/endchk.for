      SUBROUTINE ENDCHK ( LUNDRV, KEY, FAIL, ERRMSG )
C
C VERSION
C     16OCT13 AD Original.

C DESCRIPTION
C     Check end of Driver Table section has been reached.
C     General purpose routine.
C     Returns with FAIL=TRUE if any fields found before end of current section.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      LUNDRV !  I  LUN for driver table
      CHARACTER*3  KEY    !  I  Driver Table section
      LOGICAL      FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  NXTFLD ! Load next field from section of RFM driver file
C
C LOCAL VARIABLES
      INTEGER       LENGTH  !  No. of characters in FIELD
      CHARACTER*80  FIELD   !  Field read from driver table
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      CALL NXTFLD ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
      IF ( LENGTH .NE. 0 ) THEN
        FAIL = .TRUE.
        LENGTH = MIN ( LENGTH, 30 )
        ERRMSG = 'F-ENDCHK: Unexpected extra field in *'//KEY//
     &           ' section: '//FIELD(1:LENGTH)
      END IF
C
      END
