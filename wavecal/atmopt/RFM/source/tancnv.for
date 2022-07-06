      SUBROUTINE TANCNV ( MODE, ELETAN, GEOTAN, HGTTAN, SZNTAN, SFCTAN )
C
C VERSION
C     03-SEP-03  AD  Add SFCTAN argument
C     29-DEC-99  AD  Rewritten.
C     26-NOV-99  AD  Remove Logic Error#2 trap
C     04-NOV-99  AD  Change convergence failure from fatal to warning.
C     11-AUG-99  AD  Use D.P. for local variables 
C     29-JUL-98  AD  Original.
C
C DESCRIPTION
C     Convert tangent point specifications
C     Called by CHKLIM, FOVTAN, and TANHGT
C     This uses the relationships: n r sin(z) = n t = constant along a 
C     refracted path, where 
C       n = refractive index, n = 1 + R where R is refractivity.
C       r = radial distance from the centre of curvature 
C       z = zenith angle of ray (=angle of incidence/refraction) =90 at Tan.Pt.
C       t = radial distance of projected tangent point.
C     The function F is defined as 
C               F = n(z) z - n_obs t_obs      where F = 0 if z = z_tan
C     This is solved iteratively: 
C               z' = z - F/(dF/dz)  where dF/dz =  1 + N(z) + z (dN/dz)
C     Since refractivity is proportional to air density, it is assumed that 
C     the refractivity gradient is given by
C               dN/dz = -N/H      where H(z) is the density scale height
C
C     Note that single precision is adequate for 1 metre resolution
C     when adding Earth's radius to altitude terms, but not for adding 
C     refractivity to unity to get refractive index.
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
      DOUBLE PRECISION SZNTAN !  O  Sine(Zenith Angle) at t.p. or surface  
      LOGICAL          SFCTAN !  O  T=intersects surface, F=doesn't
C
      EXTERNAL
     &  TANITR ! Iterative determination of tangent point.
     &, VALATM ! Function to interpolate value from atmospheric profiles
      REAL VALATM
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
C LOCAL VARIABLES
      REAL     HGTTOA ! Height [km] of top of atmosphere + 1km
      DOUBLE PRECISION RADGEO ! Radial distance [km] of geometric tangent point
      DOUBLE PRECISION RADOBS ! Radial distance [km] of observer
      DOUBLE PRECISION RADTAN ! Radial distance [km] of refracted tangent point
      DOUBLE PRECISION RFROBS ! Refractivity at obs.location (or zero if space)
      DOUBLE PRECISION RFRTAN ! (Iterated) refractivity at tangent point
      DOUBLE PRECISION COSELE ! Cosine(Elevation Angle)
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Mode#1: ELETAN specified (and OBSFLG is TRUE).
      IF ( MODE .EQ. 1 ) THEN
        IF ( .NOT. OBSFLG ) STOP 'F-TANCNV: Logical Error#1'
        COSELE = COS ( DBLE ( ELETAN * DTORAD ) )
        RADOBS = DBLE ( ALTOBS ) + RADCRV
        RADGEO = RADOBS * COSELE            ! Either up or down viewing
        GEOTAN = SNGL ( RADGEO - RADCRV ) 
        IF ( ELETAN .GT. 0.0 ) THEN         ! Upward viewing
          HGTTAN = ALTOBS
          SZNTAN = COSELE
        ELSE                                ! Downward viewing
          CALL TANITR ( ALTOBS, GEOTAN, HGTTAN, SZNTAN )
        END IF
C
C Mode#2: GEOTAN specified 
      ELSE IF ( MODE .EQ. 2 ) THEN
        IF ( OBSFLG ) THEN
          RADGEO = DBLE ( GEOTAN ) + RADCRV
          RADOBS = DBLE ( ALTOBS ) + RADCRV
          COSELE = RADGEO / RADOBS
          ELETAN = - SNGL ( ACOS ( COSELE ) ) / DTORAD
          CALL TANITR ( ALTOBS, GEOTAN, HGTTAN, SZNTAN )
        ELSE
          HGTTOA = HGTATM(NATM) + 1.0        ! Ensure refractivity=0.0
          CALL TANITR ( HGTTOA, GEOTAN, HGTTAN, SZNTAN )
        END IF
C          
C Mode#3: HGTTAN specified 
      ELSE IF ( MODE .EQ. 3 ) THEN
        RADTAN = DBLE ( HGTTAN ) + RADCRV
        IF ( OBSFLG ) THEN
          RADOBS = DBLE ( ALTOBS ) + RADCRV
          IF ( GEOFLG ) THEN
            RADGEO = RADTAN
          ELSE
            RFRTAN = DBLE ( VALATM ( HGTTAN, 'RFR' ) )
            RFROBS = DBLE ( RFRATM(IATOBS) ) 
            RADGEO = ( 1.0D0 + RFRTAN ) * RADTAN / ( 1.0D0 + RFROBS )
          END IF
          GEOTAN = SNGL ( RADGEO - RADCRV )
          COSELE = RADGEO / RADOBS
          ELETAN = - SNGL ( ACOS ( COSELE ) ) / DTORAD
          SZNTAN = 1.0D0 
        ELSE
          IF ( GEOFLG ) THEN
            RADGEO = RADTAN
          ELSE
            RFRTAN = DBLE ( VALATM ( HGTTAN, 'RFR' ) )
            RADGEO = ( 1.0D0 + RFRTAN ) * RADTAN 
          END IF
          GEOTAN = SNGL ( RADGEO - RADCRV )
          SZNTAN = 1.0D0 
        END IF
      ELSE
        STOP 'F-TANCNV: Logical Error#2'
      END IF
C
      SFCTAN = SFCFLG .AND. HGTTAN .EQ. HGTATM(1)
C
      END
