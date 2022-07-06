      SUBROUTINE TANHGT ( FAIL, ERRMSG )
C
C VERSION
C     03-SEP-03  AD  Set SFCTAN
C     18-MAR-01  AD  Add CHKLIM
C     29-DEC-99  AD  Use GRACNV instead of TANCNV if GRAFLG enabled.
C     18-JUL-98  AD  Original.
C
C DESCRIPTION
C     Calculate full tangent height information
C     Called once by INPCHK.
C     Only calculates for nominal output tangent heights since the information
C     for additional FOV-convolution tangent paths is calculated during FOVTAN.
C
      IMPLICIT NONE
C
C ARGUMENTS
      LOGICAL      FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  CHKLIM ! Check path for Limb-viewing mode
     &, GRACNV ! Convert tangent point specifications for 2D atmosphere.
     &, TANCNV ! Convert tangent point specifications
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL VARIABLES
      INTEGER  ITAN   ! Tangent path counter
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Check that all user-specified tangent heights are valid for limb ray-tracing
      DO ITAN = 1, NTAN
        CALL CHKLIM ( USRTAN(ITAN), USRELE, USRGEO, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END DO
C
      IF ( USRELE ) THEN
        DO ITAN = 1, NTAN
          ELETAN(ITAN) = USRTAN(ITAN)
          IF ( GRAFLG ) THEN 
            CALL GRACNV ( 1, ELETAN(ITAN), GEOTAN(ITAN), HGTTAN(ITAN),
     &                    SZNTAN(ITAN), PSITAN(ITAN), SFCTAN(ITAN) )
          ELSE
            CALL TANCNV ( 1, ELETAN(ITAN), GEOTAN(ITAN), HGTTAN(ITAN),
     &                    SZNTAN(ITAN), SFCTAN(ITAN) )
          END IF
        END DO
      ELSE IF ( USRGEO ) THEN
        DO ITAN = 1, NTAN
          GEOTAN(ITAN) = USRTAN(ITAN)
          IF ( GRAFLG ) THEN
            CALL GRACNV ( 2, ELETAN(ITAN), GEOTAN(ITAN), HGTTAN(ITAN),
     &                    SZNTAN(ITAN), PSITAN(ITAN), SFCTAN(ITAN) )
          ELSE
            CALL TANCNV ( 2, ELETAN(ITAN), GEOTAN(ITAN), HGTTAN(ITAN),
     &                    SZNTAN(ITAN), SFCTAN(ITAN) )
          END IF
        END DO
      ELSE
        DO ITAN = 1, NTAN
          HGTTAN(ITAN) = USRTAN(ITAN)
          IF ( GRAFLG ) THEN
            CALL GRACNV ( 3, ELETAN(ITAN), GEOTAN(ITAN), HGTTAN(ITAN),
     &                    SZNTAN(ITAN), PSITAN(ITAN), SFCTAN(ITAN) )
          ELSE
            CALL TANCNV ( 3, ELETAN(ITAN), GEOTAN(ITAN), HGTTAN(ITAN),
     &                    SZNTAN(ITAN), SFCTAN(ITAN) )
          END IF
        END DO
      END IF
C
      END
