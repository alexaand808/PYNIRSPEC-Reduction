      SUBROUTINE TANITR ( HGTINI, HGTGEO, HGTTAN, SZNTAN )
C
C VERSION
C     31-DEC-99  AD  Original. Previously part of TANCNV.
C
C DESCRIPTION
C     Iterative determination of tangent point.
C     Called by TANCNV
C     This uses the relationships: n r sin(z) = n t = constant along a 
C     refracted path, where 
C       n = refractive index, n = 1 + R where R is refractivity.
C       r = radial distance from the centre of curvature 
C       z = zenith angle of ray (=angle of incidence/refraction) =90 at Tan.Pt.
C       t = radial distance of projected tangent point.
C     The function F is defined as 
C               F = n(z) z - n_ini t_ini      where F = 0 if z = z_tan
C     This is solved iteratively: 
C               z' = z - F/(dF/dz)  where dF/dz =  1 + N(z) + z (dN/dz)
C     Since refractivity is proportional to air density, it is assumed that 
C     the refractivity gradient is given by
C               dN/dz = -N/H      where H(z) is the density scale height
C
      IMPLICIT NONE
C
C ARGUMENTS
      REAL HGTINI              !  I  Initial altitude [km]
      REAL HGTGEO              !  I  Initial geometric tangent height [km]
      REAL HGTTAN              !  O  Tangent height or surface height [km]
      DOUBLE PRECISION SZNTAN  !  O  Zenith angle at tan.pt or surface [deg]
C
      EXTERNAL
     &  VALATM ! Function to interpolate value from atmospheric profiles
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
C LOCAL CONSTANTS
      INTEGER       MAXTER ! Max. number of iterations for tangent point
        PARAMETER ( MAXTER = 10 ) ! Arbitrary value
C
      REAL          ALTTOL ! Tolerance [km] for refrac.tan.ht.calculations
        PARAMETER ( ALTTOL = 0.0003 ) ! 0.0003 = 30cm
C
C LOCAL VARIABLES
      INTEGER  ITER   ! Iteration counter
      REAL     HGTPRV ! Value of ALTTAN on previous iteration
      REAL     HGTSFC ! Height [km] of surface (=base of atmosphere)
      REAL     HGTTST ! (Iterated) refracted tangent point altitude
      DOUBLE PRECISION DFBYDZ ! dF/dz
      DOUBLE PRECISION DRBYDZ ! d(Refractivity)/dz, assuming Rf.has dens.sc.ht.
      DOUBLE PRECISION RADGEO ! Radial distance [km] of geometric tangent point
      DOUBLE PRECISION RADSFC ! Radial distance [km] of surface
      DOUBLE PRECISION RADTST ! (Iterated) Radial distance [km] of refr.t.pt.
      DOUBLE PRECISION RFRINI ! Refractivity at obs.location (or zero if space)
      DOUBLE PRECISION RFRSFC ! Refractivity at base of atmosphere
      DOUBLE PRECISION RFRTST ! (Iterated) refractivity at tangent point
      DOUBLE PRECISION F      ! Iterated function for tangent point (see above)
      DOUBLE PRECISION FCONST ! Constant term in calculation of F
      DOUBLE PRECISION SZNSFC ! Sin(Zen.angle) at surface
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      HGTSFC = HGTATM(1)
      RADSFC = DBLE ( HGTSFC ) + RADCRV
      RADGEO = DBLE ( HGTGEO ) + RADCRV
C
      IF ( GEOFLG ) THEN
        RFRSFC = 0.0D0
        RFRINI = 0.0D0
      ELSE
        RFRSFC = DBLE ( RFRATM(1) ) 
        RFRINI = DBLE ( VALATM ( HGTINI, 'RFR' ) )
      END IF
C
      SZNSFC = ( 1.0D0 + RFRINI ) * RADGEO / ( 1.0D0 + RFRSFC ) / RADSFC 
C
      IF ( SZNSFC .LE. 1.0D0 ) THEN      ! Ray intersects surface
        SZNTAN = SZNSFC
        HGTTAN = HGTSFC
      ELSE IF ( GEOFLG ) THEN
        SZNTAN = 1.0D0
        HGTTAN = HGTGEO
      ELSE
        SZNTAN = 1.0D0
        HGTTST = HGTGEO                        ! First guess
        FCONST  = ( 1.D0 + RFRINI ) * RADGEO 
        DO ITER = 1, MAXTER
          RFRTST = DBLE ( VALATM ( HGTTST, 'RFR' ) )
          RADTST = DBLE ( HGTTST ) + RADCRV
          F = ( 1.0D0 + RFRTST ) * RADTST - FCONST
          DRBYDZ = - RFRTST / DBLE ( VALATM ( HGTTST, 'DSH' ) )
          DFBYDZ = 1.D0 + RFRTST + ( RADTST * DRBYDZ )
          HGTPRV = HGTTST
          HGTTST = MAX ( HGTTST - SNGL ( F / DFBYDZ ), HGTSFC )
          IF ( ABS ( HGTPRV - HGTTST ) .LE. ALTTOL ) THEN   
            HGTTAN = HGTTST
            GOTO 100
          END IF
        END DO
C If routine fails to converge, assume this is due to small oscillations 
C about the true value 
        WRITE ( *, * ) 'W-TANITR: Failed to converge, residual=',
     &     HGTTST - HGTPRV
        HGTTAN = 0.5 * ( HGTTST + HGTPRV )
 100    CONTINUE
      END IF
C
      END
