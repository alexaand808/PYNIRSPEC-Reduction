      SUBROUTINE CTMH2O ( IPTH ) 
C
C VERSION
C     13-APR-10  AD  Change to MT_CKD_2.5
C     14-MAR-06  AD  Change to MT_CKD_1.10
C     12-JAN-04  AD  Change to MT_CKD_1.00
C     05-MAR-02  AD  Change to CKD 2.41
C     02-AUG-99  AD  Force SNGL( ) to C2 in calculation of C11
C     24-JAN-99  AD  Remove GL2FLG, modules FINCOM.INC and CTMGL2
C     03-MAR-97  AD  Version 3.
C     18-DEC-96  AD  Convert to using CKD continuum by default
C     01-OCT-96  AD  Version 2.
C     19-SEP-96  AD  Modified to incorporate CKD 2.1 algorithms
C     01-SEP-96  AD  Version 1.
C 
C DESCRIPTION    
C     Calculate H2O continuum absorption across entire widemesh grid.
C     Called by RFMCTM for each Path containing H2O.
C     This version uses the MT_CKD continuum, previous version of this 
C     subroutine using the CKD continuum is now renamed CTMCKD.
C     v2.5 algorithms extracted from AER program contnm.f 
C     
C REFERENCES
C     http://rtweb.aer.com/
C     Mlawer, M.J., D.C. Tobin, and S.A. Clough, 
C     A Revised Perspective on the Water Vapor Continuum:  The MT_CKD Model, 
C     in preparation for JQSRT, 2003.
C
C     The MT_CKD_1.1 differs from 1.0 in having an extra fudge factor applied
C     to the foreign broadening between approx 100-600 cm-1. 
C
C     The MT_CKD_2.5 differs from 2.5 in having different fudge factors applied
C     to both self and foreign broadening.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER IPTH   !  I  Path# for calculation
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
C
C COMMON VARIABLES
      INCLUDE 'widcom.inc' ! Wide mesh data
      INCLUDE 'pthcom.inc' ! RFM path data
C
C LOCAL CONSTANTS
C Parameters for self-broadening correction
      REAL    BETA                 
        PARAMETER ( BETA = 350.0 ) ! Central Wno [cm-1] for correction 
      REAL    F1                   
        PARAMETER ( F1 = 0.25 )
      REAL    N_S                   
        PARAMETER ( N_S = 6.0 ) 
C Parameters for foreign-broadening correction
      REAL    BETA1
        PARAMETER ( BETA1 = 57.83 )
      REAL    BETA2
        PARAMETER ( BETA2 = 630.0 )
      REAL    C_1
        PARAMETER ( C_1 = -0.42 )
      REAL    C_2
        PARAMETER ( C_2 = 0.3 )
      REAL    F0
        PARAMETER ( F0 = 0.06 )
      REAL    HWSQ1
        PARAMETER ( HWSQ1 = 240.0**2 )
      REAL    N_1
        PARAMETER ( N_1 = 8.0 )
      REAL    N_2
        PARAMETER ( N_2 = 8.0 )
      REAL    V0F1
        PARAMETER ( V0F1 = 255.67 )  ! Central Wno [cm-1] for correction
C Parameters associated with tabulated data in h2omtc.inc
      REAL    DELWNO               ! Waveno. steps for H2O Continuum data
        PARAMETER ( DELWNO = 10.0 )
      REAL    WNOLOW               ! Lower Waveno. for H2O Continuum data
        PARAMETER ( WNOLOW = 0.0 )       ! (Upper Waveno.= 20000.0 )
      INTEGER NWNO                 ! No.Waveno. steps for H2O Contin.data
        PARAMETER ( NWNO = 2001 )
C
C LOCAL VARIABLES      
      INTEGER I      ! Counter for loops in DATA statements
      INTEGER IW,JW  ! Indices on H2O ctm wno axis below,above Half-WM point.
      INTEGER IX     ! Index of point in XFCREV array below Half-WM point.
      INTEGER IWD2   ! Counter for half-wide-mesh grid (0:NWID*2)
      INTEGER JWID   ! Wide mesh interval#
      LOGICAL EVEN   ! T = Half-wide mesh point coincides with wide mesh point
      REAL    C11    ! 0.5 C2/T [cm], where C2 is a radiation constant 
      REAL    CTMPTH ! Continuum absorption for current path
      REAL    CTWFRN ! FASCOD-correction applied to Foreign-broadening CWFRN 
      REAL    CTWSLF ! Self-broadening interp. to Half-WM and Path Temperature.
      REAL    CWFRN  ! Continuum Foreign Coeff. interpolated to Half-WM point.
      REAL    CW260  ! Continuum 260K Self-Broad Coeff.interpolated to Half-WM 
      REAL    CW296  ! Continuum 296K Self-Broad Coeff.interpolated to Half-WM 
      REAL    DW,EW  ! Fraction of tabulation intvl below,above Half WM pt.
      REAL    DX     ! Fraction of XFAC array intvl below Half WM pt.
      REAL    FSCAL  ! Scaling factor for foreign broadening 100-600cm-1 region
      REAL    H2O296(NWNO) ! Self-broadening Coefficients at 296 K
      REAL    H2O260(NWNO) ! Self-broadening Coefficients at 260 K
      REAL    H2OFRN(NWNO) ! Foreign broadening coefficients
      REAL    SFAC   ! Interpolated value of XFAC for current wavenumber
      REAL    VDSQ1, VDMSQ1  ! functions of Wno association with f.b.correction
      REAL    VF1, VF2, VMF1 ! functions of Wno associated with f.b.correction
      REAL    WNO    ! Wavenumber of current half-wide-mesh point
      REAL    XW     ! Position of half-WM point on H2O continm tabulation axis
      REAL    XX     ! Position of half-WM point on XFCREV tabulation axis
      REAL    XFCREV(0:14) ! Self-broadening modification factor 820-960@10cm-1
      REAL    XFCRV1(120)  ! Self-broad. factor 2000-3190@10cm-1
      INCLUDE 'h2omtc.inc' ! CKD Continuum data (unchanged with MT_CKD v2.5)
C
C XFCREV is unchanged from v1.10 to v2.5
      DATA XFCREV / 1.003, 1.009, 1.015, 1.023, 1.029, 1.033, 1.037,
     &       1.039, 1.040, 1.046, 1.036, 1.027, 1.01, 1.002,  1.00  /
C
C XFCRV1 is new with v2.5
      DATA XFCRV1 / 
     &  1.000, 1.040, 1.080, 1.120, 1.160, 1.200, 1.240, 1.280, 1.318,
     &  1.357, 1.404, 1.453, 1.499, 1.553, 1.608, 1.674, 1.746, 1.818,
     &  1.899, 1.984, 2.078, 2.174, 2.276, 2.385, 2.502, 2.624, 2.747,
     &  2.883, 3.018, 3.170, 3.321, 3.473, 3.635, 3.803, 3.974, 4.144,
     &  4.327, 4.500, 4.703, 4.887, 5.102, 5.286, 5.498, 5.701, 5.935,
     &  6.155, 6.405, 6.633, 6.892, 7.115, 7.397, 7.650, 7.917, 8.177,
     &  8.437, 8.704, 8.953, 9.192, 9.428, 9.644, 9.821, 9.954, 10.11,
     &  10.17, 10.21, 10.26, 10.29, 10.28, 10.26, 10.20, 10.15, 10.16, 
     &  10.25, 10.02, 9.965, 10.01, 9.934, 9.847, 9.744, 9.566, 9.436, 
     &  9.181, 8.872, 8.547, 8.155, 7.730, 7.261, 6.777, 6.271, 5.807, 
     &  5.313, 4.845, 4.444, 4.074, 3.677, 3.362, 3.087, 2.826, 2.615, 
     &  2.385, 2.238, 2.148, 1.979, 1.939, 1.773, 1.696, 1.642, 1.569, 
     &  1.510, 1.474, 1.425, 1.375, 1.322, 1.272, 1.230, 1.180, 1.130, 
     &  1.080, 1.040, 1.000/
C                                            
C EXECUTABLE CODE -------------------------------------------------------------
C
      C11 = 0.5 * SNGL ( C2 ) / TEMPTH(IPTH)
C
      EVEN = .FALSE.
      DO IWD2 = 0, 2*NWID
        WNO = SNGL ( WNOWD2(IWD2) )
        EVEN = .NOT. EVEN
        XW = 1.0 + ( WNO - WNOLOW ) / DELWNO 
C
C Interpolate in wavenumber for each tabulated temperature
C
        IF ( XW .GE. 1.0 .AND. XW .LT. FLOAT(NWNO) ) THEN
          IW = INT ( XW )
          JW = IW + 1
          DW = XW - IW
          EW = 1.0 - DW
          CW296 = EW * H2O296(IW) + DW * H2O296(JW)
          CW260 = EW * H2O260(IW) + DW * H2O260(JW)
          CWFRN = EW * H2OFRN(IW) + DW * H2OFRN(JW)
C
C Interpolate Self-CTM modification factor to current wavenumber
C
          IF ( WNO .GT. 820.0 .AND. WNO .LT. 960.0 ) THEN  ! same as v1.10
            XX = ( WNO - 820.0 ) / 10.0 
            IX = INT ( XX )           ! IX should be in range 0:13
            DX = XX - IX
            SFAC = ( 1.0 - DX ) * XFCREV(IX) + DX * XFCREV(IX+1)
          ELSE IF ( WNO .GT. 2000.0 .AND. WNO .LT. 3190.0 ) THEN
            XX = ( WNO - 1990.0 ) / 10.0
            IX = INT ( XX )           ! IX should be in range 1:119
            DX = XX - IX
            SFAC = ( 1.0 - DX ) * XFCRV1(IX) + DX * XFCRV1(IX+1)
          ELSE
            SFAC = 1.0
          END IF
C
C Additional modification for region near 350cm-1
          SFAC = SFAC * ( 1.0 + ( F1 / ( 1.0 + ( WNO / BETA )**N_S ) ) )    
C
C Exponential extrapolation between temperatures
C 2.77777E-2 = 1/(296.0 - 260.0)
          CTWSLF = SFAC * CW296 * 
     &      (CW260/CW296)**( 2.77777E-2 * ( 296.0 - TEMPTH(IPTH) ) )
C
C Correction to foreign continuum 
C Calculating scaling factor correction for foreign broadening ~100-600cm-1
          VDSQ1  = ( WNO - V0F1 )**2
          VDMSQ1 = ( WNO + V0F1 )**2
          VF1    = ( ( WNO - V0F1 ) / BETA1 )**N_1
          VMF1   = ( ( WNO + V0F1 ) / BETA1 )**N_1
          VF2    = ( WNO / BETA2 )**N_2
          FSCAL = 1.0 + ( F0 + C_1 * 
     &            ( ( HWSQ1 / ( VDSQ1 + HWSQ1 + VF1 ) ) + 
     &              ( HWSQ1 / ( VDMSQ1 + HWSQ1 + VMF1 ) )  ) )
     &            / ( 1.0 + C_2 * VF2 )

          CTWFRN = CWFRN * FSCAL
C
          CTMPTH = AMTPTH(IPTH) * AVOG * 1.0E-20 
     &             * WNO * TANH ( C11 * WNO )  ! radiation term in contnm.f
     &             * TEMREF / ( TEMPTH(IPTH) * PREREF ) 
     &             * ( PPAPTH(IPTH) * CTWSLF + 
     &                     ( PREPTH(IPTH) - PPAPTH(IPTH) ) * CTWFRN )
        ELSE
          CTMPTH = 0.0
        END IF
C
        JWID = 1 + IWD2 / 2 
        IF ( EVEN ) THEN 
          IF ( JWID .LE. NWID ) THEN         ! true except when IWD2 = 2*NWID
            ABSWID(1,JWID,IPTH) = ABSWID(1,JWID,IPTH) + CTMPTH
            CNTWID(1,JWID,IPTH) = CNTWID(1,JWID,IPTH) + CTMPTH
          END IF
          IF ( JWID .GT. 1 ) THEN               ! true except when IWD2 = 0, 1
            ABSWID(3,JWID-1,IPTH) = ABSWID(3,JWID-1,IPTH) + CTMPTH
            CNTWID(3,JWID-1,IPTH) = CNTWID(3,JWID-1,IPTH) + CTMPTH
          END IF
        ELSE 
          IF ( CTMPTH .NE. 0.0 ) ICTWID(JWID,IGSPTH(IPTH)) = 1
          ABSWID(2,JWID,IPTH) = ABSWID(2,JWID,IPTH) + CTMPTH
          CNTWID(2,JWID,IPTH) = CNTWID(2,JWID,IPTH) + CTMPTH
        END IF
      END DO
C
      END
