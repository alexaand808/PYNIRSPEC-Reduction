      SUBROUTINE BBFOPT ( XTAN0, XTAN1, NFIN, OPTDEP, BBFWGT )
C
C VERSION
C     26-MAY-00  AD  Convert from SP to DP
C     11-APR-99  AD  Additional factor to give CG weak limit.
C     29-DEC-98  AD  Original
C
C DESCRIPTION
C     Calculate interpolation factor for Planck Function across layer
C     If XTAN0 = 0: interpreted as Zen/Nad viewing (=no layer curvature).
C     These weights (W) are to be applied in the equation:
C         DRad = Tau * ( B0 + W * ( BCG - B0 ) ) * ( 1 - EXP(-Opt) )
C     Where Drad = radiance contribution from layer, 
C            Tau = transmission from observer to near end of layer
C             B0 = Planck Fn at near end of layer
C            BCG = Absorber-weighted (CG) Planck Function
C            Opt = Optical depth of layer (OPTDEP)
C              W = Weight (BBFWGT) returned by routine, 0 .LE. W .LE. 1
C     
C     XTAN0 is distance from edge of layer nearest to the observer to tangent pt
C     XTAN1 is distance from edge of layer remote from the observer to tan.pt.
C     If both XTAN0 = XTAN1 = 0, interpreted as zero curvature.
C     XTAN1,XTAN0 are distances in arbitrary units,
C
      IMPLICIT NONE
C
C ARGUMENTS
      REAL     XTAN0      !  I  Distance [*] from near edge of layer to tan. pt
      REAL     XTAN1      !  I  Distance [*] from far edge of layer to tan. pt
      INTEGER  NFIN       !  I  Number of spectral points to calculate
      DOUBLE PRECISION OPTDEP(*)  !  I  Spectrum of optical depths
      DOUBLE PRECISION BBFWGT(*)  !  O  Spectrum of weights
C
C LOCAL CONSTANTS
      DOUBLE PRECISION CHIMIN          ! Minimum optical depth for full calc
        PARAMETER ( CHIMIN = 0.01D0 )  ! If CHI .LT. CHIMIN use approximation
      REAL          XLIMIT          ! Ignore curvature if layer thickness 
        PARAMETER ( XLIMIT = 0.01 ) ! < XLIMIT * distance to tangent point
C
C LOCAL VARIABLES
      INTEGER          IFIN   !  Counter for spectral points
      DOUBLE PRECISION FACTOR !  1 / ( 2 * CFACT - 1 )
      DOUBLE PRECISION CFACT  !  Curv. param (2:infty =tangent layer:no curv).
      DOUBLE PRECISION CHI    !  Optical depth of layer 
      DOUBLE PRECISION RDTAU  !  Recip. of layer transm 1 / ( 1 - exp(-CHI) )
C
C EXECUTABLE CODE --------------------------------------------------------------
C
C Both XTAN0 and XTAN1 should be either +ve or zero
      IF ( XTAN0 .LT. 0.0 ) STOP 'F-BBFOPT: Logical error#0'
      IF ( XTAN1 .LT. 0.0 ) STOP 'F-BBFOPT: Logical error#1'
C
C Curvature is ignored if both XTAN1=XTAN0=0, or if the difference between them
C is small compared to distance to tangent point (ie thin layer and/or remote 
C from tangent point).
      IF ( ABS(XTAN0-XTAN1) .LE. XLIMIT * MIN(XTAN0,XTAN1) ) THEN ! No Curve
        DO IFIN = 1, NFIN
          IF ( OPTDEP(IFIN) .LT. CHIMIN ) THEN       ! Use limiting value
            BBFWGT(IFIN) = 1.0D0
          ELSE                                       ! Explicit calculation
            CHI = OPTDEP(IFIN)
            RDTAU = 1.0D0 / ( 1.0D0 - EXP ( - CHI ) )
            BBFWGT(IFIN) = 2.0D0 + 2.0D0/CHI - 2.0D0*RDTAU
          END IF
        END DO
C
      ELSE                                           ! No curvature
C For near-side tangent layer XTAN1=0 so CFACT=1, FACTOR=1
C For far-side tangent layer, XTAN0=0 so CFACT=0, FACTOR=-1
C For other near-side layers, CFACT > 1 so +1 < FACTOR < +1/XLIMIT (approx)
C For other far-side layers,  CFACT < 0 so -1 > FACTOR > -1/XLIMIT
        CFACT = DBLE ( XTAN0 / ( XTAN0 - XTAN1 ) )
        FACTOR = 1.0D0 / ( 2.0D0 * CFACT - 1.0D0 )
        DO IFIN = 1, NFIN
          IF ( OPTDEP(IFIN) .LT. CHIMIN ) THEN       ! Use limiting value
            BBFWGT(IFIN) = 1.0D0
          ELSE                                       ! Explicit calculation
            CHI = OPTDEP(IFIN)
            RDTAU = 1.0D0 / ( 1.0D0 - EXP ( - CHI ) )
            BBFWGT(IFIN) = 6.0D0 / ( 3.0D0 + FACTOR ) 
     &        * ( 1.0D0 + 1.0D0/CHI - RDTAU + 
     &            ( FACTOR / CHI**2 ) 
     &            * ( 2.0D0 * CHI * RDTAU - CHI - 2.0D0 ) )
          END IF
        END DO
      END IF
C
      END
