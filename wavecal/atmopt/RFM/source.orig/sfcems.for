      SUBROUTINE SFCEMS ( NFIN, WNOFIN, EMSFIN )
C
C VERSION
C     24-APR-12  AD  Original.
C
C DESCRIPTION
C     Interpolate surface emissivity to fine grid
C     Called by RADSFC and RFMFLX for each wide-mesh interval if SFC flag is set
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER          NFIN       !  I  Number of wno intervals
      DOUBLE PRECISION WNOFIN(*)  !  I  Wavenumber grid
      DOUBLE PRECISION EMSFIN(*)  !  O  Surface emissivity on fine grid
C                  
      EXTERNAL
     &  DLOOKP  !  General purpose interpolation routine
C
C COMMON CONSTANTS
      INCLUDE 'rfmsiz.inc' ! MIPAS RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'sfccom.inc' ! Surface parameters
C
C LOCAL VARIABLES
      INTEGER  IFIN        ! Fine mesh counter
      INTEGER  IGUESS      ! Guess index for interpolation
C
C DATA STATEMENTS
      DATA IGUESS / 1 / 
      SAVE IGUESS
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      IF ( NSFC .GT. 1 ) THEN 
        DO IFIN = 1, NFIN  
          CALL DLOOKP ( WNOFIN(IFIN), EMSFIN(IFIN), WNOSFC, EMSSFC, 
     &                  NSFC, IGUESS )
        END DO
      ELSE
        DO IFIN = 1, NFIN
          EMSFIN(IFIN) = EMSSFC(1)
        END DO
      END IF
C
      END
