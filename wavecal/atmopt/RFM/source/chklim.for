      SUBROUTINE CHKLIM ( TANTST, USRELE, USRGEO, FAIL, ERRMSG )
C
C VERSION
C     03-SEP-03  AD  Add SFCTAN 
C     18-MAR-01  AD  Remove GETOBS
C     06-JAN-00  AD  Add GRACNV if GRAFLG enabled
C     18-JUL-98  AD  Original.
C
C DESCRIPTION
C     Check path for Limb-viewing mode
C     Called by TANHGT and FOVTAN for each limb path
C
      IMPLICIT NONE 
C
C ARGUMENTS
      REAL         TANTST  !  I  Value of sec(theta) or elev. to be checked
      LOGICAL      USRELE  !  I  T=TANTST is elevation angle, F=actual
      LOGICAL      USRGEO  !  I  T=TANTST is geom.tangent height, F=actual.
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  GRACNV ! Convert tangent point specifications for 2D atmosphere
     &, TANCNV ! Convert tangent point specifications
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'crvcom.inc' ! Local radius of curvature
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'obscom.inc' ! Observer Position
C
C LOCAL VARIABLES
      REAL    ELETAN    ! Elevation angle [deg]
      REAL    GEOTAN    ! Geometric (projected) tangent height [km]
      REAL    HGTTAN    ! Refracted (actual) tangent height [km]
      REAL    PSITAN    ! Horizontal angle of Refracted tangent pt [deg]
      DOUBLE PRECISION SZNTAN ! Sine(Zenith Angle) at tangent pt or surface
      LOGICAL SFCTAN    ! T=intersects surface, F=doesn't (dummy here)
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      FAIL = .TRUE.
C
      IF ( USRGEO ) THEN 
        GEOTAN = TANTST
C
        IF ( GEOTAN .GT. HGTATM(NATM) ) THEN
          WRITE ( ERRMSG, '(A,G12.3,A,F10.3,A)' )
     &      'F-CHKLIM: Tangent Height=', GEOTAN, 
     &      '[km] > Top of Atmosphere=', HGTATM(NATM), '[km]'
          RETURN
C
        ELSE IF ( SFCFLG .AND. GEOTAN .LT. -RADCRV ) THEN
          WRITE ( ERRMSG, '(A,G12.3,A,F10.3,A)' )
     &      'F-CHKLIM: Tangent Height=', GEOTAN, 
     &      '[km] < Centre of Earth=', -RADCRV, '[km]'
          RETURN
C
        ELSE IF ( .NOT. SFCFLG .AND. GEOTAN .LT. HGTATM(1) ) THEN
          WRITE ( ERRMSG, '(A,G12.3,A,F10.3,A)' )
     &      'F-CHKLIM: Tangent Height=', GEOTAN, 
     &      '[km] < Base of Atmosphere=', HGTATM(1), '[km]'
          RETURN
C
        ELSE
          IF ( GRAFLG ) THEN
            CALL GRACNV ( 2, ELETAN, GEOTAN, HGTTAN, 
     &                       SZNTAN, PSITAN, SFCTAN )
          ELSE
            CALL TANCNV ( 2, ELETAN, GEOTAN, HGTTAN, SZNTAN, SFCTAN )
          END IF
C
          IF ( ABS(SZNTAN) .LT. 1.0D0 .AND. .NOT. SFCFLG ) THEN
            WRITE ( ERRMSG, '(A,F8.3,A,F9.3,A)' )
     &        'F-CHKLIM: Tan.Ht.=', GEOTAN,
     &        ' [km] < atmos. after refraction, alt=', HGTTAN, ' [km]' 
            RETURN
C            
          END IF
        END IF
C
      ELSE IF ( USRELE ) THEN
        ELETAN = TANTST
C
        IF ( ABS ( ELETAN ) .GT. 90.0 ) THEN
          WRITE ( ERRMSG, '(A,G12.3)' ) 'F-CHKLIM: Specified Elev.'//
     &        ' angle outside range +/-90deg, value=', ELETAN
          RETURN
C
        ELSE IF ( ELETAN .GT. 0.0 ) THEN
          IF ( ALTOBS .GE. HGTATM(NATM) ) THEN
            ERRMSG = 'F-CHKLIM: Positive Elev.Angles not allowed '//
     &               'if OBServer above atmosphere'
            RETURN
          END IF
C
        ELSE
          IF ( GRAFLG ) THEN
            CALL GRACNV ( 1, ELETAN, GEOTAN, HGTTAN, 
     &                       SZNTAN, PSITAN, SFCTAN )
          ELSE
            CALL TANCNV ( 1, ELETAN, GEOTAN, HGTTAN, SZNTAN, SFCTAN )
          END IF
C
          IF ( GEOTAN .GT. HGTATM(NATM) ) THEN
            GEOTAN = MIN ( GEOTAN, 999999. )
            WRITE ( ERRMSG, '(A,F8.4,A,F10.3,A)' )
     &        'F-CHKLIM: Elev.Ang=', ELETAN,
     &        ' [deg] has geom.tan.pt. > atmos., alt=', GEOTAN, 
     &        ' [km]'
            RETURN
C
          ELSE IF ( .NOT. SFCFLG .AND. GEOTAN .LE. HGTATM(1) ) THEN
            WRITE ( ERRMSG, '(A,F8.4,A,F10.3,A)' )
     &        'F-CHKLIM: Elev.Ang=', ELETAN,
     &        ' [deg] has geom.tan.pt. < atmos., alt=', GEOTAN, 
     &        ' [km]'
            RETURN
C
          ELSE IF ( .NOT. SFCFLG .AND. HGTTAN .LE. HGTATM(1) ) THEN
            WRITE ( ERRMSG, '(A,F8.4,A,F10.3,A)' )
     &        'F-CHKLIM: Elev.Ang=', ELETAN,
     &        ' [deg] has refr.tan.pt. < atmos., alt=', HGTTAN, 
     &        ' [km]'
            RETURN
C
          END IF
        END IF
C
      ELSE        
        HGTTAN = TANTST
C
        IF ( HGTTAN .GT. HGTATM(NATM) ) THEN
          WRITE ( ERRMSG, '(A,G12.3,A,F10.3,A)' )
     &      'F-CHKLIM: Tangent Height=', HGTTAN, 
     &      '[km] > Top of Atmosphere=', HGTATM(NATM), '[km]'
          RETURN
C
        ELSE IF ( HGTTAN .LT. HGTATM(1) ) THEN
          WRITE ( ERRMSG, '(A,G12.3,A,F10.3,A)' )
     &      'F-CHKLIM: Tangent Height=', HGTTAN, 
     &      '[km] < Base of Atmosphere=', HGTATM(1), '[km]'
          RETURN
C
        END IF
C
      END IF
C
      FAIL = .FALSE.
C
      END
