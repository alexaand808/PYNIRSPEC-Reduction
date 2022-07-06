      SUBROUTINE RFMGRD ( IWID )
C
C VERSION
C     07-AUG-13  AD  Remove redundant local variable R
C     17-AUG-06  AD  Adjust WN1FIN, WN2FIN by WNDFIN
C     04-APR-06  AD  Change DFLOAT to DBLE
C     23-FEB-03  AD  Avoid adjusting WN1LIM,WN2LIM if direct convolution or
C                    TAB options being used
C     03-FEB-03  AD  Redefine WN1FIN to incorporate WNOOFF 
C     09-JAN-03  AD  If OUTGRD=TRUE don't add boundary points.
C     09-AUG-99  AD  Use explicit Double Precision
C                    Change WNOWID from Real to D.P.
C     21-SEP-98  AD  Correction: calculation of penultimate irreg.grid pt
C     31-JUL-98  AD  Improvement to calculation of WNOOFF 
C     07-JUL-98  AD  Rewritten
C     26-JAN-98  AD  (Standard F77) Put DATA statements *after* declarations
C     15-DEC-97  AD  Original.
C
C DESCRIPTION    
C     Sets up fine-resolution grid for each widemesh interval
C     Called by RFM for each widemesh interval.
C     If the GRD option is enabled, sets NFIN and wavenumber grid points
C     If not enabled, just sets uniform grid.
C                  
      IMPLICIT NONE
C
C ARGUMENTS      
      INTEGER IWID    !  I  Index of current wide mesh interval [IPW]
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'grdcom.inc' ! Irregular grid
      INCLUDE 'outcom.inc' ! RFM output file data
      INCLUDE 'widcom.inc' ! Wide mesh data
C
C LOCAL VARIABLES
      INTEGER IFIN     !  Counter for fine-mesh intervals
      INTEGER IGRD     !  Pointer to stored irregular grid
      DOUBLE PRECISION DW     ! Negligible wno [/cm] to allow for rounding
      DOUBLE PRECISION WN1LIM ! Lower wno limit [/cm] for mapping irreg.grid
      DOUBLE PRECISION WN2LIM ! Upper wno limit [/cm] for mapping irreg.grid
C
C DATA STATEMENTS
      DATA IGRD / 1 /
      SAVE IGRD
C
C EXECUTABLE CODE -------------------------------------------------------------
C 
C The widemesh will be an integer multiple of the output and finemesh grids,
C but the widemesh grid is adjusted to start at an integral wavenumber, so 
C there may be some offset between the widemesh boundary and the 
C finemesh/output grid - if there is, adjust by WNDFIN to match the o/p grid.
      WN1FIN  = WNOWID(IWID-1) + WNDFIN
      WN2FIN  = WNOWID(IWID) + WNDFIN
      DW      = WNRFIN * 0.1D0
C
C If no irreg. grid used for this interval, just flag each fine mesh point for
C calculation, maintaining value of NFIN set in SPCFIN.
      IF ( .NOT. GRDFLG .OR. NGRD .EQ. 0 ) THEN
        DO IFIN = 1, NFIN
          WNOFIN(IFIN) = WN1FIN + (IFIN-1) * WNRFIN
        END DO
      ELSE               
C
C Irreg grid defined for this spectral range        
C
        IF ( IWID .EQ. 1 ) IGRD = 1
        WN1LIM = WN1FIN 
        WN2LIM = WN2FIN - WNRFIN
        IFIN = 0
C
C Unless directly interpolating to output grid, it is necessary to insert
C one or two additional grid points at each end of the widemesh interval 
C to constrain the interpolation to the full fine grid.
        IF ( .NOT. ( TABFLG .OR. OUTGRD ) ) THEN
          IFIN = 1
          WN1LIM = WN1LIM + WNRFIN
          WN2LIM = WN2LIM - WNRFIN
          WNOFIN(1) = WN1FIN
          IF ( .NOT. TWOGRD ) THEN
            WN1LIM = WN1LIM + WNRFIN
            WN2LIM = WN2LIM - WNRFIN
            IFIN = 2
            WNOFIN(2) = WN1FIN + WNRFIN
          END IF
        END IF
C
C Increment grid point counter until next grid point lies inside wno.range
C IGRD points to the next grid point to be used.
C NB: possible that this leaves no points at all, leaving IGRD=NGRD+1
        DO WHILE ( IGRD .LE. NGRD .AND. 
     &               WNOGRD(IGRD) .LT. WN1LIM - DW )
          IGRD = IGRD + 1
        END DO
C
C Fill in any grid points within this interval
C
        DO WHILE ( IGRD .LE. NGRD .AND. WNOGRD(IGRD) .LE. WN2LIM + DW )
          IFIN = IFIN + 1
          WNOFIN(IFIN) = WNOGRD(IGRD)
          IGRD = IGRD + 1
        END DO
C
C Set grid point(s) at upper edge 
C
        IF ( .NOT. ( TABFLG .OR. OUTGRD ) ) THEN
          IF ( .NOT. TWOGRD ) THEN
            IFIN = IFIN + 1
            WNOFIN(IFIN) = WN2FIN - 2 * WNRFIN
          END IF
          IFIN = IFIN + 1
          WNOFIN(IFIN) = WN2FIN - WNRFIN
        END IF
C
        NFIN = IFIN
C
      END IF
C
      END 
