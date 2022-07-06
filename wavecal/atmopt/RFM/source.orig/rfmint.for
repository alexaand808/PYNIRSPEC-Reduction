      SUBROUTINE RFMINT 
C
C VERSION
C     09-JAN-03  AD  Skip interpolation if OUTGRD is set TRUE
C     04-MAY-01  AD  Correction: Loop to LTAN, not NTAN
C     26-MAY-00  AD  Convert RADTAN, ABSTAN from S.P. to D.P.
C     09-AUG-99  AD  Set DW = WNRFIN * 1.0D0 instead of just WNRFIN * 1.0
C     17-DEC-97  AD  Original.
C
C DESCRIPTION    
C     Interpolates irregular to regular fine resolution grid.
C     Called by RFM for each widemesh interval if the GRD option is enabled.
C                  
      IMPLICIT NONE
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'grdcom.inc' ! Irregular grid
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL VARIABLES 
      INTEGER II,IJ,IK,IL  ! Indices for interpolation points
      INTEGER IFIN         ! Counter for regular grid points
      INTEGER IFIN1,IFIN2  ! Limiting indices of calc. reg.grid points 
      INTEGER IIRG         ! Counter for irregular grid points
      INTEGER ITAN         ! Tangent path counter
      INTEGER JFIN         ! Upper index for interpolated points
      INTEGER NIRG         ! No.Irregular grid points within this interval
      DOUBLE PRECISION DW  ! Limit [/cm] for Wavenumber rounding error 
      DOUBLE PRECISION WNO ! Wavenumber [/cm]
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      DW = WNRFIN * 0.1D0
C
C Check that an irregular grid was used with this spectral range
      IF ( NGRD .EQ. 0 ) RETURN   
C
C Check if irregular grid is to be directly interpolated to output grid
      IF ( OUTGRD ) RETURN
C 
C Change NFIN from number of irregular grid point to number of regular points
C 
      NIRG = NFIN
      NFIN = NINT ( ( WN2FIN - WN1FIN ) / WNRFIN )
C
C Copy all calculated (irreg) grid points to regular grid using same arrays
C NB: do this in reverse order to avoid overwriting
C
      DO IIRG = NIRG, 1, -1
        IFIN = 1 + NINT ( ( WNOFIN(IIRG) - WN1FIN ) / WNRFIN )
        DO ITAN = 1, LTAN
          RADTAN(IFIN,ITAN) = RADTAN(IIRG,ITAN)
          TRATAN(IFIN,ITAN) = TRATAN(IIRG,ITAN)
          ABSTAN(IFIN,ITAN) = ABSTAN(IIRG,ITAN)
        END DO
      END DO
      IFIN1 = 1 + NINT ( ( WNOFIN(1) - WN1FIN ) / WNRFIN )
      IFIN2 = 1 + NINT ( ( WNOFIN(NIRG) - WN1FIN ) / WNRFIN )
C
C Set values at any fine grid points outside range of irregular grid to zero
C
      DO IFIN = 1, IFIN1 - 1
        DO ITAN = 1, LTAN
          RADTAN(IFIN,ITAN) = 0.0D0
          TRATAN(IFIN,ITAN) = 0.0D0
          ABSTAN(IFIN,ITAN) = 0.0D0
        END DO
      END DO
      DO IFIN = IFIN2 + 1, NFIN
        DO ITAN = 1, LTAN
          RADTAN(IFIN,ITAN) = 0.0D0
          TRATAN(IFIN,ITAN) = 0.0D0
          ABSTAN(IFIN,ITAN) = 0.0D0
        END DO
      END DO
C
C Interpolate the missing points using the calculated points.
C NB: WNOFIN(NIRG) still contains list of wavenumbers of calculated points.
C
      IIRG = 1
      II = 0
      IJ = 0
      IK = IFIN1
      IL = IFIN1 + 1
      WNO = WNOFIN(1) - WNRFIN
      IF ( TWOGRD ) THEN
        JFIN = IFIN2 - 1
      ELSE
        JFIN = IFIN2 - 2
      END IF
C
      DO IFIN = IFIN1, JFIN
        WNO = WNO + WNRFIN
C
C If next grid point has already been calculated, increment pointers to define 
C interpolation points for next missing point(s)
C
        IF ( ABS ( WNO - WNOFIN(IIRG) ) .LT. DW ) THEN
          II = IJ
          IJ = IK
          IF ( TWOGRD ) THEN
            IK = 1 + NINT ( ( WNOFIN(IIRG+1) - WN1FIN ) / WNRFIN )
          ELSE
            IK = IL
            IL = 1 + NINT ( ( WNOFIN(IIRG+2) - WN1FIN ) / WNRFIN )
          END IF
          IIRG = IIRG + 1
        ELSE
          CALL INTFIN ( IFIN, II, IJ, IK, IL )
        END IF
      END DO
C
C Set wavenumbers of regular fine grid
C
      DO IFIN = 1, NFIN
        WNOFIN(IFIN) = WN1FIN + ( IFIN - 1 ) * WNRFIN
      END DO
C
      END 
