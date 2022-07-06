C $Id: ctmco2.for,v 1.13 1997/03/03 17:17:03 adarby Exp $
      SUBROUTINE CTMCO2 ( IPTH ) 
C
C VERSION
C      04-AUG-99  AD  Use  SNGL(WNOWD2) instead of SNGL(WNOWD2-WNOLOW)
C                     Comment-out redundant parameter WNOUPP
C      03-MAR-97  AD  Version 3.
C      26-FEB-97  AD  Correction to accumulation index (produced noisy spectra)
C      01-OCT-96  AD  Version 2.
C      01-SEP-96  AD  Version 1.
C 
C DESCRIPTION    
C     Calculate CO2 continuum absorption on fine mesh.
C     Called by RFMCTM for each Path containing CO2.
C     
C REFERENCES
C     CLOUGH S.A.,KNEIZYS F.X.,ROTHMAN L.S.,GALLERY W.O. (1981), 
C       Atmospheric spectral transmittance and radiance: FASCOD1B
C       SPIE 277 Atmospheric Transmission  152-166
C
C     CLOUGH S.A.,KNEIZYS F.X.,ROTHMAN L.S.,ANDERSON G.P.,SHETTLE E.P. (1987)
C       Current issues in infrared atmospheric transparency
C       International meeting on Atmospheric Transparency for Satellite 
C       Applications, 15-19 Sept. 1986 Capri, Italy. Ed. G.V. Silvestrini. CUEN.
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
      REAL          WNOLOW               ! Lower Waveno. for CO2 Continuum data
        PARAMETER ( WNOLOW = 0.0 )       
      REAL          DELWNO               ! Waveno. steps for CO2 Continuum data
        PARAMETER ( DELWNO = 2.0 )
      INTEGER       NWNO                 ! No.Waveno. steps for CO2 Contin.data
        PARAMETER ( NWNO = 2001 )        ! Upper Waveno. is 4000cm-1
C
C LOCAL VARIABLES      
      INTEGER I      ! Counter for loops in DATA statements
      INTEGER IW,JW  ! Indices on CO2 ctm wno axis below,above Half-WM point.
      INTEGER IWD2   ! Counter for half-wide-mesh grid (0:NWID*2)
      INTEGER JWID   ! Wide mesh interval#
      LOGICAL EVEN   ! T = Half-wide mesh point coincides with wide mesh point
      REAL    CTMPTH ! Continuum absorption for current path
      REAL    CTW    ! Continuum absorption coefficient interp to Temp, Wno.
      REAL    CW230,CW260,CW296 ! Temp. Coeffs interpolated to Half-WM points
      REAL    DT230,DT260,DT296 ! Diff.of path temp from each tabulation temp
      REAL    DW,EW  ! Fraction of tabulation intvl below,above Half WM pt.
      REAL    XW     ! Position of half-WM point on CO2 continuum tabulation axis
      REAL    CO2296(NWNO)  ! Coefficients at 296 K
      REAL    CO2260(NWNO)  ! Coefficients at 260 K
      REAL    CO2230(NWNO)  ! Coefficients at 230 K
      INCLUDE 'co2dat.inc'  ! Coefficient data 
      SAVE    CO2296, CO2260, CO2230
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      EVEN = .FALSE.
C 
      DO IWD2 = 0, 2*NWID
        EVEN = .NOT. EVEN
        XW = 1.0 + ( SNGL ( WNOWD2(IWD2) ) - WNOLOW ) / DELWNO 
C
C Interpolate in wavenumber for each tabulated temperature
C
        IF ( XW .GE. 1.0 .AND. XW .LT. FLOAT(NWNO) ) THEN
          IW = INT ( XW )
          JW = IW + 1
          DW = XW - IW
          EW = 1.0 - DW
          CW296 = EW * CO2296(IW) + DW * CO2296(JW)
          CW260 = EW * CO2260(IW) + DW * CO2260(JW)
          CW230 = EW * CO2230(IW) + DW * CO2230(JW)
C
C Lagrangian interpolation in temperature
C
          DT230 = TEMPTH(IPTH) - 230.0
          DT260 = TEMPTH(IPTH) - 260.0
          DT296 = TEMPTH(IPTH) - 296.0
          CTW   =  5.050505E-4 * DT260 * DT296 * CW230
     &           - 9.259259E-4 * DT230 * DT296 * CW260
     &           + 4.208754E-4 * DT230 * DT260 * CW296
          CTMPTH = AMTPTH(IPTH) * ( PREPTH(IPTH) / PREREF ) * CTW
        ELSE
          CTMPTH = 0.0
        END IF
C
        JWID = 1 + IWD2 / 2 
        IF ( EVEN ) THEN 
          IF ( JWID .LE. NWID  ) THEN           ! true except when IWD2 = 2*NWID
            ABSWID(1,JWID,IPTH) = ABSWID(1,JWID,IPTH) + CTMPTH
            CNTWID(1,JWID,IPTH) = CNTWID(1,JWID,IPTH) + CTMPTH
          END IF
          IF ( JWID .GT. 1  ) THEN              ! true except when IWD2 = 0, 1
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
C-------------------------------------------------------------------------------
C                                NOTICE
C
C     This software module is part of the MIPAS Reference Forward Model supplied
C     to the European Space Agency under ESTEC contract no. 11886/96/NL/GS.
C        
C     All rights, title and interest in and to the copyright
C     and any other proprietary right, is held by
C     The University Corporation for Atmospheric Research, having a
C     correspondence address of P.O. Box 3000, Boulder, Colorado 80307, USA.
C
C     However, note that all inquiries concerning the MIPAS
C     Reference Forward Model should be submitted to the The Manager (Projects),
C     AOPP,  Clarendon Laboratory, Parks Road, Oxford, OX1 3PU, UK.
C     (Tel: +44-8165-272900,    Fax: +44-1865-272924).  
C
C-------------------------------------------------------------------------------
