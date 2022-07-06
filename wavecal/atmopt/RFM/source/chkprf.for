      SUBROUTINE CHKPRF ( TYPE, NLEV, PROFIL, FAIL, ERRMSG )
C
C VERSION
C     15OCT13 AD Original.
C
C DESCRIPTION
C     Check atmospheric profile on input.
C     Called by FILPRF for every profile loaded.
C     Requirements for each profile type are as follows
C     'HGT' all values must be monotonically increasing
C     'TEM' all values must be positive
C     'PRE' all values must be positive and monotonically decreasing
C     '[gas]'  all values must be zero or positive unless aerosol
C
      IMPLICIT NONE
C
C ARGUMENTS
      CHARACTER*(*) TYPE      !  I  Profile type: 'HGT', 'TEM', 'PRE' or [gas].
      INTEGER       NLEV      !  I  No. of profile levels
      REAL          PROFIL(*) !  I  Profile
      LOGICAL       FAIL      !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG    !  O  Error message returned if FAIL is TRUE
C
C LOCAL VARIABLES
      INTEGER ILEV        ! Profile level counter
      REAL    PREV        ! Previous (lower alt) profile value
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      FAIL = .TRUE.
C
C Check that altitudes increase monotonically
      IF ( TYPE .EQ. 'HGT' ) THEN
        PREV = PROFIL(1)
        DO ILEV = 2, NLEV
          IF ( PROFIL(ILEV) .LE. PREV ) THEN
            WRITE ( ERRMSG, '(A,A,A,G10.3)' ) 'F-CHKPRF: ', TYPE, 
     &        ' profile contains non-increasing altitude value=', 
     &        PROFIL(ILEV)
            RETURN
          END IF
          PREV = PROFIL(ILEV)
        END DO
C Check that temperatures are all positive
      ELSE IF ( TYPE .EQ. 'TEM' ) THEN
        DO ILEV = 1, NLEV
          IF ( PROFIL(ILEV) .LE. 0.0 ) THEN
            WRITE ( ERRMSG, '(A,A,A,G10.3)' ) 'F-CHKPRF: ', TYPE, 
     &        ' profile contains non-positive temperature value=', 
     &        PROFIL(ILEV)
            RETURN
          END IF
        END DO
C Check that pressures are all positive and decrease monotonically
      ELSE IF ( TYPE .EQ. 'PRE' ) THEN
        PREV = PROFIL(1)*1.1
        DO ILEV = 1, NLEV
          IF ( PROFIL(ILEV) .LE. 0.0 ) THEN
            WRITE ( ERRMSG, '(A,A,A,G10.3)' ) 'F-CHKPRF: ', TYPE,
     &        ' profile contains non-positive pressure value=', 
     &        PROFIL(ILEV)
            RETURN
          ELSE IF ( PROFIL(ILEV) .GE. PREV ) THEN
            WRITE ( ERRMSG, '(A,A,A,G10.3)' ) 'F-CHKPRF: ', TYPE,
     &        ' profile contains non-decreasing pressure value=', 
     &        PROFIL(ILEV)
            RETURN
          END IF
          PREV = PROFIL(ILEV)
        END DO
C Anything else, apart from aerosol, assume VMR and check all values 0 or +ve
      ELSE IF ( TYPE .NE. 'aerosol' ) THEN
        DO ILEV = 1, NLEV
          IF ( PROFIL(ILEV) .LT. 0.0 ) THEN
            WRITE ( ERRMSG, '(A,A,A,G10.3)' ) 'F-CHKPRF: ', TYPE, 
     &        ' profile contains negative VMR value=', PROFIL(ILEV)
            RETURN
          END IF
        END DO
      END IF
      FAIL = .FALSE.
C
      END
