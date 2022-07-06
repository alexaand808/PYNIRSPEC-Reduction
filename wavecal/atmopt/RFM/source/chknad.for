      SUBROUTINE CHKNAD ( TANTST, USRELE, FAIL, ERRMSG )
C
C VERSION
C      18-JUL-98  AD  Original.
C
C DESCRIPTION
C     Check path for NADir/ZENith viewing mode
C     Called by TANCHK for each path with NAD or ZEN flags enabled. 
C
      IMPLICIT NONE 
C
C ARGUMENTS
      REAL         TANTST  !  I  Value of sec(theta) or elev. to be checked
      LOGICAL      USRELE  !  I  T=TANTST is elevation angle, F=sec(zen.angle)
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG  !  O  Error message returned if FAIL is TRUE
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
C
C LOCAL CONSTANTS
      REAL          ELEMIN ! Minimum elevation angle [deg] for NAD/ZEN views
        PARAMETER ( ELEMIN = 0.1 ) ! 0.1 equivalent to ~600 airmasses
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      FAIL = .TRUE.
C
C For zenith or nadir viewing TANTST can either be elevation angle (restricted
C to +/- 90 degrees) or airmass/sec(zen.angle) ), restricted to .GE. 1.
C
      IF ( USRELE ) THEN 
C
        IF ( ZENFLG .AND. 
     &       ( TANTST .GT. 90.0 .OR. TANTST .LT. ELEMIN ) ) THEN
          WRITE ( ERRMSG, '(A,F4.1,A,G12.3)' ) 
     &      'F-CHKNAD: Specified Elev. Ang. outside ZEN range ',
     &      ELEMIN, ':+90deg, value=', TANTST
          RETURN
C
        ELSE IF ( NADFLG .AND.
     &       ( TANTST .LT. -90.0 .OR. TANTST .GT. -ELEMIN ) ) THEN
          WRITE ( ERRMSG, '(A,F4.1,A,G12.3)' ) 
     &      'F-CHKNAD: Specified Elev. Ang. outside NAD range ',
     &      -ELEMIN, ':-90deg, value=', TANTST
          RETURN
        END IF
C
      ELSE 
        IF ( TANTST .LT. 1.0 ) THEN
          WRITE ( ERRMSG, '(A,G12.3)' ) 'F-CHKNAD: Specified NADir'//
     &      '/ZENith view sec(theta) < 1, value =', TANTST
          RETURN
        END IF
      END IF
      FAIL = .FALSE.
C
      END 
