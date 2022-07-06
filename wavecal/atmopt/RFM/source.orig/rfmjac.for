      SUBROUTINE RFMJAC 
C
C VERSION
C     16-JAN-04  AD  Allow for TRA+MTX option
C     02-JAN-99  AD  Original.
C
C DESCRIPTION
C     Calculate Jacobian spectra.
C     Called by RFM for each widemesh interval if JAC option enabled.
C     It is assumed that all perturbed and unperturbed spectra have been
C     calculated already: this module just replaces the perturbed spectra
C     by the Jacobians.
C
      IMPLICIT NONE
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
C
C COMMON VARIABLES
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'jaccom.inc' ! Jacobian data
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL VARIABLES
      INTEGER   ITAN  ! Tangent path counter
      INTEGER   IJAC  ! Jacobian element counter
      INTEGER   IFIN  ! Fine mesh spectral point counter
      INTEGER   JTAN  ! Index of tan.pt. for ITAN,IJAC Jacobian calculation
C
C EXECUTABLE CODE ------------------------------------------------------------
C
      DO ITAN = 1, MTAN
        DO IJAC = 1, NJAC
          JTAN = ITNJAC(ITAN,IJAC)
          IF ( JTAN .NE. 0 ) THEN
            IF ( MTXFLG ) THEN  ! TRATAN(*,ITAN) stores absorp.coefficient
              DO IFIN = 1, NFIN ! and TRATAN(*,JTAN) already contains info
                RADTAN(IFIN,JTAN) = RADTAN(IFIN,JTAN)-RADTAN(IFIN,ITAN)
                ABSTAN(IFIN,JTAN) = ABSTAN(IFIN,JTAN)-ABSTAN(IFIN,ITAN)
              END DO
            ELSE
              DO IFIN = 1, NFIN
                RADTAN(IFIN,JTAN) = RADTAN(IFIN,JTAN)-RADTAN(IFIN,ITAN)
                ABSTAN(IFIN,JTAN) = ABSTAN(IFIN,JTAN)-ABSTAN(IFIN,ITAN)
                TRATAN(IFIN,JTAN) = TRATAN(IFIN,JTAN)-TRATAN(IFIN,ITAN)
              END DO
            END IF
          END IF
        END DO
      END DO
C
      END
