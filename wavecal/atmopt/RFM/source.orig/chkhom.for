      SUBROUTINE CHKHOM ( TANTST, FAIL, ERRMSG )
C
C VERSION
C      18-MAR-01  AD  Remove CHKWID.
C      17-JUL-98  AD  Original.
C
C DESCRIPTION
C     Check path length for HOMogeneous path calculation
C     Called by TANCHK for each path with HOM flag enabled.
C
      IMPLICIT NONE 
C
C ARGUMENTS
      REAL         TANTST  !  I  Value of homogeneous path length to be checked
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG  !  O  Error message returned if FAIL is TRUE
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      FAIL = .TRUE.
C
      IF ( TANTST .LE. 0.0 ) THEN
        WRITE ( ERRMSG, '(A,G12.3)' ) 'F-CHKHOM: Specified Homog.'//
     &    'Path Length .LE. 0, value =', TANTST
        RETURN
      END IF
C
      FAIL = .FALSE.
C
      END 
