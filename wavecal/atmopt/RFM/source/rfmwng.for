      SUBROUTINE RFMWNG ( IWID )
C
C VERSION
C     17-DEC-01  AD  Initialise AI, BI, YI0 to avoid warning messages.
C     08-SEP-99  AD  Check for peak of inverse quadratic outside WM interval
C                    rather than peak of original quadratic
C     09-AUG-99  AD  Make U = SNGL ( ) explicitly
C     24-JAN-99  AD  Replace IQDFLG by logical opposite QADFLG.
C     22-DEC-97  AD  Correction: check for ABS ( A ) .GT. ... instead of .GE.
C     17-DEC-97  AD  Remove calculation of regular fine mesh grid points 
C                    Express U using actual wavenumber rather than grid index
C     03-MAR-97  AD  Version 3.
C     27-FEB-97  AD  Revise tests for IQD - check position of peak
C     15-FEB-97  AD  Abandon IQD fit if quadratic less than zero
C     21-DEC-96  AD  Check for quadratic values below zero
C     01-OCT-96  AD  Version 2.
C     28-SEP-96  AD  Add option of Inverse Quadratic Fit (IQDFLG)
C     20-SEP-96  AD  Rename VIBFLG to NTEFLG
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION    
C     Interpolate tabulated wide mesh absorption across fine mesh.
C     Called by RFM for each widemesh interval.
C     The accumulated line-wing contributions for each wide-mesh interval are
C     interpolated quadratically to the fine-mesh grid.
C     Also initialises ABSFIN, CNTFIN arrays
C
C     Fitted quadratic is Y = Y0 + B*U + A*U*(U-1)
C     where B = Y2 - Y0 and A = 2 * ( Y0 + Y2 - 2 * Y1 )
C
C     The default option is to try to fit an `inverse quadratic', ie fit a 
C     quadratic to the reciprocals of the 3 interpolation points (which gives a
C     better fit for Lorentz line wings). 
C     But this inversion can only be used safely if certain criteria are met:
C     (1) all three fitted points must have values > ABSMIN
C     (2) there has to be a significant curvature component ( U**2 coeff.)
C     (3) the max/min has to lie outside the widemesh interval (this is to 
C         prevent unrealistic large absorption peaks being created based on
C         three low-absorption points)
C     If any of these tests fail, a normal quadratic fit is applied directly to
C     the interpolation points. The normal quadratic is *ALWAYS* applied if the
C     QAD option is enabled.
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
      INCLUDE 'pthcom.inc' ! RFM path data
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'widcom.inc' ! Wide mesh data
C
C LOCAL CONSTANTS
      REAL          ABSMIN            ! Min value for inverse quadratic
        PARAMETER ( ABSMIN = 1.0E-30 )
      REAL          CRVMIN            ! Min curvature for inverse quad.fit
        PARAMETER ( CRVMIN = 1.0E-3 )
C
C LOCAL VARIABLES
      INTEGER IFIN     !  Counter for fine-mesh intervals
      INTEGER IPTH     !  Path counter [IPATH]
      REAL    A, B     !  Differences between Yn values
      REAL    AI,BI    !  Differences between YIn values
      REAL    F        !  Value of quadratic function 
      REAL    U        !  Scaled Value (0:1) of Fine-Mesh points
      REAL    UPEAK    !  Position of quadratic peak in scaled grid units
      REAL    Y0,Y1,Y2 !  Values of quadratic functions at Half-Mesh pts
      REAL    YI0,YI1,YI2 !  Values of inv.quadratic functions at Half-Mesh pts
      LOGICAL USEIQD   !  T=apply inverse quadratic fit algorithm
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C These variables initialised just to avoid warning messages
      AI = 0.0
      BI = 0.0
      YI0 = 0.0
C
C Loop over all calculated paths
C 
      DO IPTH = 1, NCLC    
        Y0 = ABSWID(1,IWID,IPTH)
        Y1 = ABSWID(2,IWID,IPTH)
        Y2 = ABSWID(3,IWID,IPTH)
        B  = Y2 - Y0
        A  = 2.0 * ( Y0 + Y2 - 2.0*Y1 )
C
        USEIQD = .NOT. QADFLG
        IF ( USEIQD ) THEN                   ! test conditions for using IQD fit
          USEIQD = MIN ( Y0, Y1, Y2 ) .GT. ABSMIN 
          IF ( USEIQD ) THEN
            YI0 = 1.0 / Y0
            YI1 = 1.0 / Y1
            YI2 = 1.0 / Y2
            BI  = YI2 - YI0
            AI  = 2.0 * ( YI0 + YI2 - 2.0*YI1 )
            USEIQD = ( ABS ( AI ) .GT. CRVMIN * ABS ( BI ) )
            IF ( USEIQD ) THEN
              UPEAK = 0.5 * ( AI - BI ) / AI
              USEIQD = ( UPEAK .LE. 0.0 .OR. UPEAK .GE. 1.0 )
            END IF
          END IF
        END IF
C
        IF ( USEIQD ) THEN
          DO IFIN = 1, NFIN                     ! Loop over fine mesh frequency
            U = SNGL ( ( WNOFIN(IFIN) - WN1FIN ) / DELWID )
            F = YI0 + U * ( BI + AI * (U - 1.0) ) ! guaranteed .GE. ABSMIN
            ABSFIN(IFIN,IPTH) = 1.0 / F
          END DO
        ELSE                                    ! Normal quadratic fit
          DO IFIN = 1, NFIN                     ! Loop over fine mesh frequency
            U = SNGL ( ( WNOFIN(IFIN) - WN1FIN ) / DELWID )
            ABSFIN(IFIN,IPTH) = MAX ( Y0 + U * ( B + A*(U-1.0) ), 0.0 )
          END DO
        END IF
C
        IF ( NTEFLG ) THEN
          Y0 = CNTWID(1,IWID,IPTH)
          Y1 = CNTWID(2,IWID,IPTH)
          Y2 = CNTWID(3,IWID,IPTH)
          B  = Y2 - Y0
          A  = 2.0 * ( Y0 + Y2 - 2.0*Y1 )
          DO IFIN = 1, NFIN
            U = SNGL ( ( WNOFIN(IFIN) - WN1FIN ) / DELWID )
            CNTFIN(IFIN,IPTH) = MAX ( Y0 + U * ( B + A*(U-1.0) ), 0.0 )
          END DO
        END IF
      END DO
C   
      END
