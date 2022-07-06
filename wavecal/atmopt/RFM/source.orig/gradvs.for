      SUBROUTINE GRADVS ( P, R, T, DPDS, DRDS, DTDS )
C
C VERSION
C     26-JAN-01  AD  Correction: GEOFLG set RFRIDX=1 (not 0).
C     03-JAN-00  AD  Original.
C
C DESCRIPTION
C     Derivatives for 2-D ray-tracing.
C
      IMPLICIT NONE
C
C ARGUMENTS
      DOUBLE PRECISION P    !  I  D.P. horizontal angle [rad]
      DOUBLE PRECISION R    !  I  Radius coordinate of path [km]
      DOUBLE PRECISION T    !  I  D.P. zenith angle [rad]
      DOUBLE PRECISION DPDS !  O  Rate of change of P with S
      DOUBLE PRECISION DRDS !  O  Rate of change of R with S
      DOUBLE PRECISION DTDS !  O  Rate of change of T with S
C
      EXTERNAL
     &  VALGRA ! Function interpolate value from 2D atmospheric field.
      REAL VALGRA
C
C GLOBAL CONSTANTS
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
C
C COMMON VARIABLES
      INCLUDE 'crvcom.inc' ! Local radius of curvature
      INCLUDE 'flgcom.inc' ! RFM option flags
C
C LOCAL VARIABLES
      REAL             HGT  ! Height [km] of path segment
      REAL             PSI  ! Horizontal angle [deg] of path
      DOUBLE PRECISION DNDP ! Rate of change of refractive index with P
      DOUBLE PRECISION DNDR ! Rate of change of refractive index with R
      DOUBLE PRECISION REFRAC ! Refractivity (=Ref.index - 1)
      DOUBLE PRECISION RFRIDX ! Refractive index
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      HGT = SNGL ( R - RADCRV )
      PSI = SNGL ( P ) / DTORAD
      IF ( GEOFLG ) THEN
        DNDR = 0.0D0
        DNDP = 0.0D0
        RFRIDX = 1.0D0
      ELSE
        REFRAC = DBLE ( VALGRA ( HGT, PSI, 'RFR', 0 ) ) 
        RFRIDX = 1.0D0 + REFRAC
        DNDR = - REFRAC / DBLE ( VALGRA ( HGT, PSI, 'DSH', 0 ) )
        DNDP = DBLE ( VALGRA ( HGT, PSI, 'DRP', 0 ) )
      END IF
C
      DRDS = COS ( T )
      DPDS = SIN ( T ) / R
      DTDS = -SIN ( T ) * ( 1.0D0/R + DNDR/RFRIDX ) + 
     &        COS ( T ) * DNDP / ( RFRIDX * R )
C
      END
