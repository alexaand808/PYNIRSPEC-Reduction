      SUBROUTINE CHKFLX ( HGTTST, FAIL, ERRMSG )
C
C VERSION
C     27-APR-00  AD  Original. 
C
C DESCRIPTION
C     Check output levels for Flux calculations.
C     Called by TANCHK for each path with FLX flag enabled. 
C     Also inserts any additional profile levels that may be required.
C
      IMPLICIT NONE 
C
C ARGUMENTS
      REAL         HGTTST  !  I  Value of level to be checked
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG  !  O  Error message returned if FAIL is TRUE
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes      
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      FAIL = .TRUE.
C
      IF ( HGTTST .LT. HGTATM(1) ) THEN
        WRITE ( ERRMSG, '(A,G12.3,A,G12.3,A)' )
     &    'F-CHKFLX: Level = ', HGTTST, 
     &    ' km below base of atmosphere (', HGTATM(1), ' km)'
        RETURN
      ELSE IF ( HGTTST .GT. HGTATM(NATM) ) THEN
        WRITE ( ERRMSG, '(A,G12.3,A,G12.3,A)' )
     &    'F-CHKFLX: Level = ', HGTTST, 
     &    ' km above top of atmosphere (', HGTATM(NATM), ' km)'
        RETURN
      END IF
C
      FAIL = .FALSE.
C
      END 
