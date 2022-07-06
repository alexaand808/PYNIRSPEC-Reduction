      SUBROUTINE RFMLOS 
C
C VERSION
C     17-DEC-01  AD  Correction: check for NTAN=1, not ITAN=1
C     09-JUN-01  AD  Original.
C
C DESCRIPTION
C     Calculate LOS Jacobian spectra.
C     Called by RFM for each widemesh interval if LOS option enabled.
C     This is the last stage before writing output spectra, so it is assumed
C     that nominal (FOV,ILS-convolved) output spectra are already calculated 
C     and located in elements 1:NTAN, so additional NTAN spectra are generated.
C     Profile jacobians may also be calculated and stored in NTAN+1:LTAN, but 
C     only relevant for setting the offset ITAN.
C
      IMPLICIT NONE
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'rfmcon.inc' ! Constants for RFM calculations
C
C COMMON VARIABLES
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL VARIABLES
      INTEGER   ITAN  ! Tangent path counter
      INTEGER   IFIN  ! Fine mesh spectral point counter
      REAL      F1,F2,F3 ! Weight Factors for lower,medium,upper tangent paths
      REAL      Z21,Z32,Z31 ! (+ve) differences in tangent points [km or deg]
C
C EXECUTABLE CODE ------------------------------------------------------------
C
C If LOS flag selected, INPTAN checks for at least 2 tangent heights
      IF ( NTAN .EQ. 1 ) STOP 'F-RFMLOS: Logical error#1'
C
C INPTAN and RFMPTB shouldhave already checked that LTAN+NTAN fits MAXTAN.
      IF ( LTAN + NTAN .GT. MAXTAN ) STOP 'F-RFMLOS: Logical error#2'
C
C For spectra at tangent altitude limits use simple difference as gradient
      F1 = PTBLOS / ( USRTAN(2) - USRTAN(1) )                
      DO IFIN = 1, NFIN                              ! lowest tan.hgt
        RADTAN(IFIN,LTAN+1) = F1 * ( RADTAN(IFIN,2) - RADTAN(IFIN,1) )
        ABSTAN(IFIN,LTAN+1) = F1 * ( ABSTAN(IFIN,2) - ABSTAN(IFIN,1) )
        TRATAN(IFIN,LTAN+1) = F1 * ( TRATAN(IFIN,2) - TRATAN(IFIN,1) )
      END DO
      F1 = PTBLOS / ( USRTAN(NTAN) - USRTAN(NTAN-1) )
      DO IFIN = 1, NFIN                              ! highest tan.hgt
        RADTAN(IFIN,LTAN+NTAN) = 
     &    F1 * ( RADTAN(IFIN,NTAN) - RADTAN(IFIN,NTAN-1) )
        ABSTAN(IFIN,LTAN+NTAN) = 
     &    F1 * ( ABSTAN(IFIN,NTAN) - ABSTAN(IFIN,NTAN-1) )
        TRATAN(IFIN,LTAN+NTAN) = 
     &    F1 * ( TRATAN(IFIN,NTAN) - TRATAN(IFIN,NTAN-1) )
      END DO
C
C For other cases, use quadratic fit to three points
      DO ITAN = 2, NTAN-1                   
        Z21 = USRTAN(ITAN) - USRTAN(ITAN-1)
        Z32 = USRTAN(ITAN+1) - USRTAN(ITAN)
        Z31 = Z32 + Z21
        F1  = -Z32/Z21/Z31 * PTBLOS
        F2  = (Z32-Z21)/Z21/Z32 * PTBLOS
        F3  = Z21/Z31/Z32 * PTBLOS
        DO IFIN = 1, NFIN
          RADTAN(IFIN,LTAN+ITAN) = F1 * RADTAN(IFIN,ITAN-1) + 
     &                             F2 * RADTAN(IFIN,ITAN)   +
     &                             F3 * RADTAN(IFIN,ITAN+1) 
          ABSTAN(IFIN,LTAN+ITAN) = F1 * ABSTAN(IFIN,ITAN-1) + 
     &                             F2 * ABSTAN(IFIN,ITAN)   +
     &                             F3 * ABSTAN(IFIN,ITAN+1) 
          TRATAN(IFIN,LTAN+ITAN) = F1 * TRATAN(IFIN,ITAN-1) + 
     &                             F2 * TRATAN(IFIN,ITAN)   +
     &                             F3 * TRATAN(IFIN,ITAN+1) 
        END DO
      END DO
C
      END
