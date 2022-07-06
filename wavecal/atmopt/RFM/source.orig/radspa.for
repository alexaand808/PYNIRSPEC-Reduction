      SUBROUTINE RADSPA ( ITAN )
C
C VERSION
C     14-SEP-11  AD  Original.
C
C DESCRIPTION
C     Cosmic background contribution to radiative transfer calc.
C     Called by RFMRAD for ray path intersecting reaching top of atmosphere.
C     This assumes RADTAN and TRATAN already contains integrated radiances
C     and transmission from the observer to space
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER   ITAN       !  I  Tangent path index
C
      EXTERNAL
     &  PLANCK ! Calculate Planck emission spectrum at specified temperature 
C
C GLOBAL CONSTANTS
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
C      INCLUDE 'rfmcon.inc' ! Constants for RFM calculations
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'jaccom.inc' ! Jacobian data
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL VARIABLES
      INTEGER   IFIN       ! Fine mesh grid point counter (1:NFIN)
      INTEGER   IJAC       !  Counter for Jacobian elements
      INTEGER   JTAN       ! Index of original plus Jacobian tan.paths
      DOUBLE PRECISION BBFFIN(MAXFIN) ! Planck emission
C
C EXECUTABLE CODE ------------------------------------------------------------
C
      IF ( TEMSPA .EQ. 0.0 ) RETURN
C
      CALL PLANCK ( TEMSPA, NFIN, WNOFIN, BBFFIN )
C
      DO IFIN = 1, NFIN
        RADTAN(IFIN,ITAN) = RADTAN(IFIN,ITAN) + 
     &    TRATAN(IFIN,ITAN) * BBFFIN(IFIN) 
      END DO
C
      IF ( JACFLG ) THEN
        IJAC = 0
        DO WHILE ( IJAC .LT. NJAC )
          IJAC = IJAC + 1
          JTAN = ITNJAC(ITAN,IJAC)
          IF ( JTAN .NE. 0 ) THEN
            DO IFIN = 1, NFIN
              RADTAN(IFIN,JTAN) = RADTAN(IFIN,JTAN) + 
     &        TRATAN(IFIN,JTAN) * BBFFIN(IFIN) 
            END DO
          END IF
        END DO
      END IF
C
      END 
