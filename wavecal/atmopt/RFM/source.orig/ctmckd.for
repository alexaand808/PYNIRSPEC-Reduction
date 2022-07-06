      SUBROUTINE CTMCKD ( IPTH ) 
C
C VERSION
C     12-JAN-04  AD  Original. Old version of CTMH2O (05MAR02).
C 
C DESCRIPTION    
C     Calculate CKD H2O continuum absorption across entire widemesh grid.
C     Called by RFMCTM for each Path containing H2O if CKD switch is set TRUE.
C     NB: The RFM default is to use the new MT_CKD continuum in CTMH2O.
C     
C     This uses CKD v2.41. To use CKD v2.1 set the flag USE241 = FALSE
C     
C REFERENCES
C     CLOUGH S.A.,KNEIZYS F.X.,DAVIES R.,GAMACHE R., TIPPING R. (1980)
C       Theoretical line shape for H2O vapour: Application to the continuum
C       Atmospheric Water Vapour, Eds. A.Deepak,T.D.Wilkeson, L.H.Ruhnke, 
C       Academic Press, New York
C
C     CLOUGH S.A.,KNEIZYS F.X.,ROTHMAN L.S.,GALLERY W.O. (1981), 
C       Atmospheric spectral transmittance and radiance: FASCOD1B
C       SPIE 277 Atmospheric Transmission  152-166
C
C     CLOUGH S.A.,KNEIZYS F.X.,ROTHMAN L.S.,ANDERSON G.P.,SHETTLE E.P. (1987)
C       Current issues in infrared atmospheric transparency
C       International meeting on Atmospheric Transparency for Satellite 
C       Applications, 15-19 Sept. 1986 Capri, Italy. Ed. G.V. Silvestrini.CUEN.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER IPTH   !  I  Path# for calculation
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'rfmcon.inc' ! Constants for RFM calculations
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
C
C COMMON VARIABLES
      INCLUDE 'widcom.inc' ! Wide mesh data
      INCLUDE 'pthcom.inc' ! RFM path data
C
C LOCAL CONSTANTS
      REAL    WNOLOW               ! Lower Waveno. for H2O Continuum data
        PARAMETER ( WNOLOW = 0.0 )       ! (Upper Waveno.= 20000.0 )
      REAL    DELWNO               ! Waveno. steps for H2O Continuum data
        PARAMETER ( DELWNO = 10.0 )
      INTEGER NWNO                 ! No.Waveno. steps for H2O Contin.data
        PARAMETER ( NWNO = 2001 )
C
C LOCAL VARIABLES      
      INTEGER I      ! Counter for loops in DATA statements
      INTEGER IW,JW  ! Indices on H2O ctm wno axis below,above Half-WM point.
      INTEGER IX     ! Index of point in XFAC array below Half-WM point.
      INTEGER IWD2   ! Counter for half-wide-mesh grid (0:NWID*2)
      INTEGER JWID   ! Wide mesh interval#
      LOGICAL EVEN   ! T = Half-wide mesh point coincides with wide mesh point
      LOGICAL USE241 ! T= Use v2.41, F= use v2.1
      REAL    A1,A2,A3 ! Factors in computing absorption coefficient
      REAL    C11    ! 0.5 C2/T [cm], where C2 is a radiation constant 
      REAL    CTMPTH ! Continuum absorption for current path
      REAL    CTWFRN ! FASCOD-correction applied to Foreign-broadening CWFRN 
      REAL    CTWSLF ! Self-broadening interp. to Half-WM and Path Temperature.
      REAL    CWFRN  ! Continuum Foreign Coeff. interpolated to Half-WM point.
      REAL    CW260  ! Continuum 260K Self-Broad Coeff.interpolated to Half-WM 
      REAL    CW296  ! Continuum 296K Self-Broad Coeff.interpolated to Half-WM 
      REAL    DW,EW  ! Fraction of tabulation intvl below,above Half WM pt.
      REAL    DX     ! Fraction of XFAC array intvl below Half WM pt.
      REAL    SFAC   ! Interpolated value of XFAC for current wavenumber
      REAL    V2     ! (Wavenumber-offset)**2
      REAL    WNO    ! Wavenumber of current half-wide-mesh point
      REAL    XW     ! Position of half-WM point on H2O continm tabulation axis
      REAL    XX     ! Position of half-WM point on XFAC tabulation axis
      REAL    H2O296(NWNO) ! Self-broadening Coefficients at 296 K
      REAL    H2O260(NWNO) ! Self-broadening Coefficients at 260 K
      REAL    H2OFRN(NWNO) ! Foreign broadening coefficients
      REAL    XFAC(0:50)   ! Self-broadening modification factor 700-1200cm-1
      INCLUDE 'h2ockd.inc' ! CKD Continuum data
c
      DATA ( XFAC(I), I = 0, 50 ) /
     &    1.00000,1.01792,1.03767,1.05749,1.07730,1.09708,
     &    1.10489,1.11268,1.12047,1.12822,1.13597,1.14367,
     &    1.15135,1.15904,1.16669,1.17431,1.18786,1.20134,
     &    1.21479,1.22821,1.24158,1.26580,1.28991,1.28295,
     &    1.27600,1.26896,1.25550,1.24213,1.22879,1.21560,
     &    1.20230,1.18162,1.16112,1.14063,1.12016,1.10195,
     &    1.09207,1.08622,1.08105,1.07765,1.07398,1.06620,
     &    1.05791,1.04905,1.03976,1.02981,1.00985,1.00000,
     &    1.00000,1.00000,1.00000/
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      USE241 = .TRUE.             ! Switch on CKD v2.41
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
          IF ( WNO .LE. 700.0 .OR. WNO .GE. 1200.0 ) THEN
            SFAC = 1.0
          ELSE
            XX = ( WNO - 700.0 ) / 10.0 
            IX = INT ( XX )
            DX = XX - IX
            SFAC = ( 1.0 - DX ) * XFAC(IX) + DX * XFAC(IX+1)
          END IF
C
C Clough et al. use exponential extrapolation between temperatures
C 2.77777E-2 = 1/(296.0 - 260.0)
C
          CTWSLF = SFAC * CW296 * 
     &      (CW260/CW296)**( 2.77777E-2 * ( 296.0 - TEMPTH(IPTH) ) )
C
C FASCOD correction to self continuum (1.9.85); factor of 0.7667 at 1050 cm-1
C VOS2=1050, FACTRS2=-0.2333, HWSQ2 = 200**2
C
          V2 = ( WNO - 1050.0 )**2
          CTWSLF = CTWSLF * ( 1.0 - 0.2333 * 4.0E4 / 
     &                        ( V2 + 4.0E4 ) )
C
C CKD 2.1 Additional correction terms for self-continuum
C V0S3=1310, FACTRS3=0.15, HWSQ3=120**2, BETAS3=5E-6
C
          V2 = ( WNO - 1310.0 )**2
          CTWSLF = CTWSLF * ( 1.0 - 0.15 * 14400.0 / 
     &      ( V2 * ( 1.0 + 5.0E-6 * V2 ) + 14400.0 ) )
C
C CKD 2.41 further correction (mostly affecting < 100cm-1)
C V0S1=0, FACTRS1=0.688 (CKD 2.2 used FACTRS1=0.3), HWSQ1=100**2, BETAS1=1E-4
          IF ( USE241 ) THEN
            V2 = WNO**2
            CTWSLF = CTWSLF * ( 1.0 + 0.688 * 1.0E4 /
     &        ( V2 * ( 1.0 + 1.0E-4 * V2 ) + 1.0E4 ) )
          END IF
C
C Correction to foreign continuum 
C
          CTWFRN = CWFRN
C V0F=1130, FACTRF=-0.97, HWSQF=330**2=1.089E5, BETAF=8E-11, 
          V2 = ( WNO - 1130.0 )**2
          CTWFRN = CTWFRN * ( 1.0 - 0.97 * 1.089E5 / 
     &      ( V2 * ( 1.0 + 8.0E-11 * V2 * V2 ) + 1.089E5 ) )
C
          IF ( USE241 ) THEN
C CKD 2.41 correction VOF1=350.0, FACTRF1=-0.7, HWSQF1=200**2, BETAF1=5E-9 
            V2 = ( WNO - 350.0 )**2
            CTWFRN = CTWFRN * ( 1.0 - 0.7 * 4.0E4 /
     &        ( V2 * ( 1.0 + 5.0E-9 * V2 * V2 ) + 4.0E4 ) )
C CKD 2.41 correction V0F1A=630.0, FACTRF1a=0.75, HWSQF1a=65**2, BETAF1a=2E-8
            V2 = ( WNO - 630.0 )**2
            CTWFRN = CTWFRN * ( 1.0 + 0.75 * 4225.0 /
     &        ( V2 * ( 1.0 + 2.0E-8 * V2 * V2 ) + 4225.0 ) )
C CKD 2.41 correction V0F3=1975.0, FACTRF3=-0.65, HWSQF3=250**2, BETAF3=5E-6
            V2 = ( WNO - 1975.0 )**2
            CTWFRN = CTWFRN * ( 1.0 - 0.65 * 62500.0 /
     &        ( V2 * ( 1.0 + 5.0E-6 * V2 ) + 62500.0 ) )     ! NB V^4, not V^6
          ELSE
C CKD 2.1 correction VOF2=1900, FACTRF2=-0.6, HWSQF2=150**2, BETAF2=3E-6
            V2 = ( WNO - 1900.0 )**2
            CTWFRN = CTWFRN * ( 1.0 - 0.6 * 22500.0 / 
     &        ( V2 * ( 1.0 + 3.0E-6 * V2 ) + 22500.0 ) )     ! NB V^4, not V^6
          END IF
C
          A1 = WNO * AMTPTH(IPTH) * AVOG * TANH ( C11 * WNO )
          A2 = TEMREF / ( TEMPTH(IPTH) * PREREF )
          A3 = 1.0E-20 * ( PPAPTH(IPTH) * CTWSLF + 
     &                     ( PREPTH(IPTH) - PPAPTH(IPTH) ) * CTWFRN )
          CTMPTH = A1 * A2 * A3          
        ELSE
          CTMPTH = 0.0
        END IF
C
        JWID = 1 + IWD2 / 2 
        IF ( EVEN ) THEN 
          IF ( JWID .LE. NWID ) THEN            ! true except when IWD2 = 2*NWID
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
