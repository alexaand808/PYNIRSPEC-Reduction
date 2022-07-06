      SUBROUTINE GRACNV ( MODE, ELETAN, GEOTAN, HGTTAN, 
     &                          SZNTAN, PSITAN, SFCTAN )
C
C VERSION
C     03-SEP-03  AD  Add SFCTAN argument
C     28-NOV-00  AD  Fix floating zero divide in Mode=3
C     19-JAN-00  AD  Original.
C
C DESCRIPTION
C     Convert tangent point specifications for 2D atmosphere.
C     Called by CHKLIM, FOVTAN, and TANHGT
C
C     The ELETAN variable is only set if OBSFLG =.T. (Observer alt.specified)
C     In the case of rays intersecting the surface, HGTTAN is set = surface
C     altitude and SZNTAN is set to the sine(Zenith Angle) at that point.
C     For +ve ELETAN values, HGTTAN is set to observer altitude, and GEOTAN set
C     to projected tangent altitude in opposite direction (-ELETAN).
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER          MODE   !  I  #1=ELETAN,#2=GEOTAN,#3=HGTTAN init.defined 
      REAL             ELETAN ! I/O  Elevation Angle [deg]
      REAL             GEOTAN ! I/O  Geometric (=Projected) tangent alt. [km]
      REAL             HGTTAN ! I/O  Refracted (=Actual) tangent altitude [km]
      DOUBLE PRECISION SZNTAN !  O   Sine(Zenith Angle) at t.p. or surface  
      REAL             PSITAN !  O   LOS angle at t.p. or surface [deg]
      LOGICAL          SFCTAN !  O   T=intersects surface, F=doesn't
C
      EXTERNAL
     &  RAYGRA ! Ray-tracing in a 2-D atmosphere (z,psi).
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data  
      INCLUDE 'crvcom.inc' ! Local radius of curvature
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'obscom.inc' ! Observer Position
C
C LOCAL CONSTANTS
      INTEGER MAXTER                  ! Max. no. of iterations for tan point
        PARAMETER ( MAXTER = 10 )     ! Arbitrary value
      REAL ALTTOL                     ! Tolerance [km] for refrac.tan.ht.calc
        PARAMETER ( ALTTOL = 0.0003 ) ! 0.0003 = 30cm
      REAL DPATH                      ! Path increment for ray-tracing
        PARAMETER ( DPATH = 1.0 )     ! 1km seems to give accuracy limit
C
C LOCAL VARIABLES
      INTEGER  ITER   ! Iteration counter
      REAL     DS     ! Max Path length [km] for each ray-trace segment
      REAL     DZENDH ! Rate of change of ZENOBS with HGTTST
      REAL     HGTPRV ! Previous value of HGTTST
      REAL     HGTTOA ! Height [km] of top of atmosphere
      REAL     HGTTST ! Test value iterated towards HGTTAN
      REAL     PSITOA ! Horizontal angle [deg] of top of atmosphere (+ve)
      REAL     ZENOBS ! Zenith angle [deg] at observer
      REAL     ZENPRV ! Previous estimate of ZENOBS
      REAL     ZENTAN ! Zenith angle [deg] at actual tan.ht (nominally 270)
      REAL     ZENTOA ! Zenith angle [deg] at top of atmosphere
      DOUBLE PRECISION COSELE ! Cosine(Elevation Angle)
      DOUBLE PRECISION COSTOA ! Cosine(Elevation Angle at top of atmosphere)
      DOUBLE PRECISION RADGEO ! Radial distance [km] of geometric tangent point
      DOUBLE PRECISION RADOBS ! Radial distance [km] of observer
      DOUBLE PRECISION RADTAN ! Radial distance [km] of tangent point
      DOUBLE PRECISION RADTOA ! Radial distance [km] of top of atmosphere
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      DS = DPATH       ! 1.0 = 1km  max path length for ray tracing`
C
C Mode#1: ELETAN specified (and OBS = TRUE)
      IF ( MODE .EQ. 1 ) THEN           
        IF ( .NOT. OBSFLG ) STOP 'F-GRACNV: Logical Error#1'
        COSELE = COS ( DBLE ( ELETAN * DTORAD ) )
        RADOBS = DBLE ( ALTOBS ) + RADCRV 
        RADGEO = RADOBS * COSELE
        GEOTAN = SNGL ( RADGEO - RADCRV )   ! For both +/- values of ELETAN
        IF ( ELETAN .GT. 0.0 ) THEN         ! Upward viewing
          HGTTAN = ALTOBS
          SZNTAN = COSELE
        ELSE
          ZENOBS = 270.0 + ELETAN
          ZENTAN = 270.0
          CALL RAYGRA ( ALTOBS, PSIOBS, ZENOBS, 3, DS,
     &                  HGTTAN, PSITAN, ZENTAN )
          IF ( ZENTAN .EQ. 270.0 ) THEN ! Do this to avoid precision problems
            SZNTAN = -1.0D0      
          ELSE
            SZNTAN = SIN ( DBLE ( ZENTAN * DTORAD ) )
          END IF
        END IF
C
C Mode#2: GEOTAN specified 
      ELSE IF ( MODE .EQ. 2 ) THEN
        IF ( OBSFLG ) THEN                      ! Observer location specified
          RADGEO = DBLE ( GEOTAN ) + RADCRV 
          RADOBS = DBLE ( ALTOBS ) + RADCRV 
          COSELE = RADGEO / RADOBS
          ELETAN = -SNGL ( ACOS ( COSELE ) ) / DTORAD  ! -ve for down viewing
          ZENOBS = 270.0 + ELETAN
          ZENTAN = 270.0
          CALL RAYGRA ( ALTOBS, PSIOBS, ZENOBS, 3, DS, 
     &                  HGTTAN, PSITAN, ZENTAN )
          IF ( ZENTAN .EQ. 270.0 ) THEN ! Do this to avoid precision problems
            SZNTAN = -1.0D0      
          ELSE
            SZNTAN = SIN ( DBLE ( ZENTAN * DTORAD ) )
          END IF
        ELSE
          HGTTOA = HGTATM(NATM)
          RADTOA = DBLE ( HGTTOA ) + RADCRV
          RADGEO = DBLE ( GEOTAN ) + RADCRV
          PSITOA = SNGL ( ACOS ( RADGEO / RADTOA ) ) / DTORAD
          ZENTOA = 270.0 - PSITOA                  ! Assume PSITOA is positive
          ZENTAN = 270.0
          CALL RAYGRA ( HGTTOA, PSITOA, ZENTOA, 3, DS, 
     &                  HGTTAN, PSITAN, ZENTAN )
          IF ( ZENTAN .EQ. 270.0 ) THEN ! Do this to avoid precision problems
            SZNTAN = -1.0D0      
          ELSE
            SZNTAN = SIN ( DBLE ( ZENTAN * DTORAD ) )
          END IF
        END IF
C
C Mode#3: HGTTAN specified
      ELSE IF ( MODE .EQ. 3 ) THEN
        IF ( OBSFLG ) THEN
          RADOBS = DBLE ( ALTOBS ) + RADCRV
          RADTAN = DBLE ( HGTTAN ) + RADCRV
          COSELE = RADTAN / RADOBS         ! Geometric tan.hgt for first guess
          ZENOBS = 270.0 - SNGL ( ACOS ( COSELE ) ) / DTORAD   !-ve elevation
          ZENPRV = 270.0
          HGTPRV = ALTOBS
          DO ITER = 1, MAXTER
            ZENTAN = 270.0                  
            DS = DPATH
            CALL RAYGRA ( ALTOBS, PSIOBS, ZENOBS, 3, DS,
     &                    HGTTST, PSITAN, ZENTAN )
            IF ( ABS ( HGTTST - HGTTAN ) .LT. ALTTOL ) GOTO 100  ! Converged
            IF ( ABS ( HGTTST - HGTPRV ) .LT. ALTTOL ) GOTO 100  ! Converged
            DZENDH = ( ZENOBS - ZENPRV ) / ( HGTTST - HGTPRV )
            ZENPRV = ZENOBS
            HGTPRV = HGTTST
            ZENOBS = ZENPRV + ( HGTTAN - HGTTST ) * DZENDH 
          END DO
          WRITE ( *, * ) 'W-GRACNV: Max.iterations exceeded, residual=',
     &      HGTTAN - HGTTST
          ZENOBS = 0.5 * ( ZENOBS + ZENPRV ) ! Assuming oscillating about soln
  100     CONTINUE
          IF ( ZENOBS .GT. 270.0 ) STOP 'F-GRACNV: Logical Error#2'
          ELETAN = ZENOBS - 270.0                          ! Should be -ve
          COSELE = COS ( DBLE ( ELETAN * DTORAD ) )
          RADGEO = RADOBS * COSELE 
          GEOTAN = SNGL ( RADGEO - RADCRV )
          SZNTAN = 1.0D0
        ELSE
          PSITAN = 0.0
          ZENTAN = 90.0                   ! 90 since heading *towards* observer
          HGTTOA = HGTATM(NATM)
          CALL RAYGRA ( HGTTAN, PSITAN, ZENTAN, 1, DS,
     &                  HGTTOA, PSITOA, ZENTOA )
          RADTOA = DBLE ( HGTTOA ) + RADCRV
          COSTOA = SIN ( DBLE ( ZENTOA * DTORAD ) )
          ELETAN = - SNGL ( ACOS ( COSTOA ) ) / DTORAD
          RADGEO = RADTOA * COSTOA
          GEOTAN = SNGL ( RADGEO - RADCRV )
          SZNTAN = 1.0D0
        END IF
      ELSE
        STOP 'F-GRACNV: Logical Error#3'      ! Mode .NE. 1, 2 or 3
      END IF
C
      SFCTAN = SFCFLG .AND. HGTTAN .EQ. HGTATM(1)
C        
      END
