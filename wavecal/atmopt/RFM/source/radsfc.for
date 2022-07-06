      SUBROUTINE RADSFC ( ITAN )
C
C VERSION
C     23-APR-12  AD  Allow for spectrally varying surface emissivity
C                    Remove REFLCT argument
C     03-MAY-05  AD  Add jdxcon.inc
C     01-JAN-04  AD  Add LEV flag
C     03-SEP-03  AD  Allow for surface parameter jacobians
C     07-NOV-00  AD  Check for tau=0 before taking log
C     26-MAY-00  AD  Convert from S.P. to D.P.
C     27-APR-00  AD  Correction: recalculate ABSTAN for reflected ray
C     09-AUG-99  AD  Use explicit DBLE( ) in calculation of RADTAN
C     30-MAR-99  AD  Correction: allow for Jacobian calculations.
C     30-DEC-98  AD  Original. Extracted from RFMRAD.
C
C DESCRIPTION
C     Surface contribution to radiative transfer calculation
C     Called by RFMRAD for ray path intersecting the surface.
C     This assumes RADTAN and TRATAN already contains integrated radiances
C     and transmission from the observer to the surface
C     If surface has unity emissivity, TRATAN and ABSTAN are unchanged and 
C     contain only the atmospheric values from the observer to the surface.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER   ITAN       !  I  Tangent path index
C
      EXTERNAL
     &  LEVUPD ! Set up calculations for intermediate output levels
     &, PLANCK ! Calculate Planck emission spectrum at specified temperature 
     &, SFCEMS ! Interpolate surface emissivity to fine grid
C
C GLOBAL CONSTANTS
      INCLUDE 'jdxcon.inc' ! Indices for non-VMR Jacobians
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
      INCLUDE 'rfmcon.inc' ! Constants for RFM calculations
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'jaccom.inc' ! Jacobian data
      INCLUDE 'sfccom.inc' ! Surface parameters
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL VARIABLES
      INTEGER   IFIN       ! Fine mesh grid point counter (1:NFIN)
      INTEGER   IJAC       !  Counter for Jacobian elements
      INTEGER   JTAN       ! Index of original plus Jacobian tan.paths
      INTEGER   KTAN       !  Index of tan.paths used for intermed.level o/p
      INTEGER   KTAN1, KTAN2 ! Range of KTAN values used for current ITAN
      REAL      TEM        ! Surface temperature plus any perturbation
      REAL      TEMPTB     ! Surface temperature perturbation
      DOUBLE PRECISION ABS    ! Surface absorptivity
      DOUBLE PRECISION EMS    ! Surface emissivity including perturbation
      DOUBLE PRECISION BBFFIN(MAXFIN) ! Planck emission
      DOUBLE PRECISION EMSFIN(MAXFIN) ! Surface emissivity
      DOUBLE PRECISION EMSPTB ! Surface emissivity perturbation
C
C EXECUTABLE CODE ------------------------------------------------------------
C
      IJAC = 0
      JTAN = ITAN
      EMSPTB = 0.0D0
      TEMPTB = 0.0
      CALL SFCEMS ( NFIN, WNOFIN, EMSFIN )
C
  100 CONTINUE
      TEM = TEMSFC + TEMPTB
      CALL PLANCK ( TEM, NFIN, WNOFIN, BBFFIN )
      DO IFIN = 1, NFIN
        EMS = EMSFIN(IFIN) + EMSPTB
        RADTAN(IFIN,JTAN) = RADTAN(IFIN,JTAN) + 
     &    TRATAN(IFIN,JTAN) * BBFFIN(IFIN) * EMS
        ABS = 1.0 - EMS
        IF ( ABS .NE. 0.0D0 ) THEN
          TRATAN(IFIN,JTAN) = TRATAN(IFIN,JTAN) * ABS
          ABSTAN(IFIN,JTAN) = ABSTAN(IFIN,JTAN) - 
     &      LOG ( MAX ( DBLE(ARGMIN), ABS ) ) 
        END IF
      END DO
C
      IF ( LEVFLG ) THEN
        CALL LEVUPD ( .FALSE., ITAN, 0, 0, KTAN1, KTAN2 )
        DO KTAN = KTAN1, KTAN2
          DO IFIN = 1, NFIN
            EMS = EMSFIN(IFIN) + EMSPTB
            RADTAN(IFIN,KTAN) = RADTAN(IFIN,KTAN) + 
     &        TRATAN(IFIN,KTAN) * BBFFIN(IFIN) * EMS 
            ABS = 1.0 - EMS
            IF ( ABS .NE. 0.0D0 ) THEN
              TRATAN(IFIN,KTAN) = TRATAN(IFIN,KTAN) * ABS
              ABSTAN(IFIN,KTAN) = ABSTAN(IFIN,KTAN) - 
     &          LOG ( MAX ( DBLE(ARGMIN), ABS ) ) 
            END IF
          END DO
        END DO
      END IF
C
      IF ( JACFLG ) THEN
        DO WHILE ( IJAC .LT. NJAC )
          IJAC = IJAC + 1
          JTAN = ITNJAC(ITAN,IJAC)
          IF ( IGSJAC(IJAC) .EQ. MAXGAS + JDXSFT ) THEN ! Perturb surface temp.
            TEMPTB = PTBSFT                        ! PTBSFT in rfmcon.inc
          ELSE
            TEMPTB = 0.0
          END IF
          IF ( IGSJAC(IJAC) .EQ. MAXGAS + JDXSFE ) THEN ! Ptb sfc emissivity
            EMSPTB = PTBSFE
          ELSE
            EMSPTB = 0.0D0
          END IF
          IF ( JTAN .NE. 0 ) GOTO 100
        END DO
      END IF
C
      END 
