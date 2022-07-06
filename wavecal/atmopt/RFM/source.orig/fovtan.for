      SUBROUTINE FOVTAN ( FAIL, ERRMSG )
C
C VERSION
C     03-SEP-03  AD  Add SFCTAN
C     14-FEB-03  AD  Reduce error message to C*79
C     26-SEP-00  AD  Reduce TOLFAC from 0.1 to 0.01,
C                    Take ABS(altfov differences) for TOLFOV
C                    Set TOLFOV=0 if CLC flag selected
C     29-DEC-99  AD  Add GRACNV instead of TANCNV if GRAFLG used 
C     11-AUG-99  AD  Remove redundant atmcom.inc
C     24-JAN-99  AD  Replace OFMFLG by FVZFLG
C     18-JUL-98  AD  Allow for Angular FOV tabulation, also describe FOV
C                    relative to geom.tangent points rather than refracted
C                    (although continue to use refr. with OFM flag)
C     16-JAN-98  AD  Correction to error message if FOV Hgt too low.
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Add tangent heights to allow FOV convolution.
C     Called once by INPCHK if FOV flag is enabled. 
C
      IMPLICIT NONE 
C
C ARGUMENTS
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  CHKLIM ! Check path for Limb-viewing mode
     &, GRACNV ! Convert tangent point specifications for 2D atmosphere.
     &, RFMLOG ! Write text message to RFM log file.
     &, TANCNV ! Convert tangent point specifications
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C LOCAL CONSTANTS
      REAL          TOLFAC        ! Fractional tolerance allowed on FOV spacing
        PARAMETER ( TOLFAC = 0.01 )
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'fovcom.inc' ! FOV data
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL VARIABLES
      INTEGER      IFOV           ! Counter for FOV points
      INTEGER      ITAN           ! Counter for output tangent heights
      INTEGER      JTAN           ! Counter for all tangent heights
      INTEGER      MODE           ! Mode# for TANCNV conversion
      LOGICAL      GEOFOV         ! T=Treat FOV as Geom.Tangent Heights.
      LOGICAL      USETAN         ! T=Use existing tan.pth for FOV conv.
      REAL         TOLFOV(MAXFOV) ! Height tolerance for FOV convolution points
      REAL         TSTFOV         ! Value of FOV-convolution tan.path
      CHARACTER*80 MESSGE         ! Text message for Log file
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Set up tangent height tolerances for FOV convolution
C NB NFOV must be at least 3 to accomodate two zeros plus one non-zero value.
C
      IF ( CLCFLG ) THEN          ! If CLC flag, only use exact matches
        DO IFOV = 1, NFOV
          TOLFOV(IFOV) = 0.0
        END DO
      ELSE                        ! Otherwise use TOLFOV * spacing
        DO IFOV = 2, NFOV-1
          TOLFOV(IFOV) = TOLFAC * 
     &      MIN ( ABS ( ALTFOV(IFOV+1) - ALTFOV(IFOV) ), 
     &            ABS ( ALTFOV(IFOV) - ALTFOV(IFOV-1) )  )
        END DO
        TOLFOV(1)    = TOLFOV(2)
        TOLFOV(NFOV) = TOLFOV(NFOV-1)
      END IF
C
C Switch off calculation for set of nominal tangent paths - may be switched on
C again if path required for FOV convolution
C
      DO ITAN = 1, NTAN
        CLCTAN(ITAN) = .FALSE.
      END DO
C
C For each output tangent height (1:NTAN) check if any already listed will do
C for FOV convolution, otherwise add to list (NTAN+1:MTAN)
C
      GEOFOV = .NOT. ( ELEFOV .OR. FVZFLG )
      DO ITAN = 1, NTAN
        DO IFOV = 1, NFOV
          IF ( ELEFOV ) THEN          ! FOV expressed as elevation angles ...
            TSTFOV = ELETAN(ITAN) + ALTFOV(IFOV)
          ELSE IF ( GEOFOV ) THEN     ! FOV expressed as geom.tan.hts
            TSTFOV = GEOTAN(ITAN) + ALTFOV(IFOV)
          ELSE                        ! FOV expressed as ref.tan.hts (for OFM)
            TSTFOV = HGTTAN(ITAN) + ALTFOV(IFOV)
          END IF
          DO JTAN = 1, MTAN
            IF ( ELEFOV ) THEN
              USETAN = ABS ( ELETAN(JTAN) - TSTFOV ) .LE. TOLFOV(IFOV)  
            ELSE IF ( GEOFOV ) THEN
              USETAN = ABS ( GEOTAN(JTAN) - TSTFOV ) .LE. TOLFOV(IFOV) 
            ELSE
              USETAN = ABS ( HGTTAN(JTAN) - TSTFOV ) .LE. TOLFOV(IFOV) 
            END IF
            IF ( USETAN ) THEN
              IFVTAN(IFOV,ITAN) = JTAN
              CLCTAN(JTAN) = .TRUE.
              GOTO 100
            END IF
          END DO
C No matches in list so far, so add another and check limits.
          IF ( MTAN .GE. MAXTAN ) THEN
            FAIL = .TRUE. 
            WRITE ( ERRMSG, '(A,I11)' ) 
     &        'F-FOVTAN: No space for extra Tan.Hgt for '//
     &        'FOV Convolution, >MAXTAN=', MAXTAN
            RETURN
          ELSE 
            CALL CHKLIM ( TSTFOV, ELEFOV, GEOFOV, FAIL, ERRMSG )
            IF ( FAIL ) RETURN
          END IF
          MTAN = MTAN + 1
          IF ( ELEFOV ) THEN
            ELETAN(MTAN) = TSTFOV
            MODE = 1
          ELSE IF ( GEOFOV ) THEN
            GEOTAN(MTAN) = TSTFOV
            MODE = 2
          ELSE
            HGTTAN(MTAN) = TSTFOV
            MODE = 3
          END IF
          IF ( GRAFLG ) THEN
            CALL GRACNV ( MODE, ELETAN(MTAN), GEOTAN(MTAN), 
     &                    HGTTAN(MTAN), SZNTAN(MTAN), PSITAN(MTAN), 
     &                    SFCTAN(MTAN) )
          ELSE
            CALL TANCNV ( MODE, ELETAN(MTAN), GEOTAN(MTAN), 
     &                    HGTTAN(MTAN), SZNTAN(MTAN), SFCTAN(MTAN) )
          END IF
          IFVTAN(IFOV,ITAN) = MTAN
          CLCTAN(MTAN) = .TRUE.
  100     CONTINUE
        END DO
      END DO
C
      WRITE ( MESSGE, '(A,I11)' )
     &  'I-FOVTAN: No.extra tangent paths required for '//
     &  'FOV Convolution=', MTAN - NTAN
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
C
      END
