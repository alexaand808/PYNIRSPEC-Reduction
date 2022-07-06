C $Id: ctmo2.for,v 1.11 1997/03/03 17:17:05 adarby Exp $
      SUBROUTINE CTMO2 ( IPTH ) 
C
C VERSION
C      22-JAN-98  AD  (Standard F77) Put DATA statements *after* declarations
C      03-MAR-97  AD  Version 3.
C      26-FEB-97  AD  Replace Timofeyev continuum with Thibault
C      01-OCT-96  AD  Version 2.
C      01-SEP-96  AD  Version 1.
C 
C DESCRIPTION    
C     Calculate O2 continuum absorption across entire widemesh grid.
C     Called by RFMCTM for each Path containing O2.
C     The continuum absorption/length s(/cm) is derived from
C
C        s(p,T) = (p/po)^2 (To/T)^2  v Bx10^-6  exp(beta(1/Tr - 1/T))
C 
C     where B and beta are tabulated functions of wavenumber, To=273K, po is
C     std.atmospheric pressure, V is the vmr of O2, Tr=296K
C     B is tabulated in (/cm /Am^2) where 1 Amagat is the ratio of the molar
C     density to the value at STP (hence the p/po To/T terms). 
C     Therefore path absorption k = s x, where x is the path length.
C     For consistency with other RFM modules, this is rearranged to be a product
C     of the path absorber amount u (kmol/cm^2), given by u = (vpx/(RT))
C     This gives a calculation 
C
C        k(p,T) = u p To^2 R / (po T) B exp(beta(1/Tr - 1/T))
C
C     where p is in atm, R is J/kmol/K, po is Pa. Note that the 10^6 scaling 
C     term for B is cancelled by the cm^3/m^3 conversion.
C
C REFERENCES
C     THIBAULT, F, V.MENOUX, R.Le DOUCEN, L. ROSENMANN, J-M. HARTMANN and 
C     Ch. BOULET
C       Infrared collision-induced absorption by O2 near 6.4um for 
C       atmospheric applications: measurements and empirical modelling.
C     Applied Optics, 36, 563-567 (1997).
C     Note that the tabulation in the paper extends from 1365-1800, but here
C     values have been added at 1360 and 1805 (0 for B, end values for BETA)
C     to give a smooth transition to zero absorption.
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
      REAL          WNOLOW            ! Lower wavenumber of continuum data [/cm]
        PARAMETER ( WNOLOW = 1360.0 ) ! Upper wavenumber = 1805.0
      REAL          DELWNO            ! Wavenumber interval of Continm.data[/cm]
        PARAMETER ( DELWNO = 5.0 )
      INTEGER       NWNO              ! No.continuum wavenumber tabulation pts.
        PARAMETER ( NWNO = 90 )
C
C LOCAL VARIABLES      
      INTEGER IW,JW  ! Indices on CO2 ctm wno axis below,above Half-WM point.
      INTEGER IWD2   ! Counter for half-wide-mesh grid (0:NWID*2)
      INTEGER JWID   ! Wide mesh interval#
      LOGICAL EVEN   ! T = Half-wide mesh point coincides with wide mesh point
      REAL    AFACT  ! Absorber amount scaling for B
      REAL    BETAW  ! BETA term interpolated to wavenumber
      REAL    BW     ! B term interpolated to wavenumber
      REAL    CTMPTH ! Continuum absorption for current path
      REAL    DW,EW  ! Fraction of tabulation intvl below,above Half WM pt.
      REAL    TFACT  ! Temperature scaling for BETA
      REAL    XW     ! Position of half-WM point on CO2 continuum tabulation axis
      REAL    B(NWNO)  ! [cm-1 Amagat-2] ( & values to be multiplied by 10^-6)
         SAVE B
      REAL   BETA(NWNO) !  [K]
         SAVE BETA
C
C DATA STATEMENTS
        DATA B / 
     &    0.000, 0.061, 0.074, 0.084, 0.096, 0.120, 0.162, 0.208, 0.246, 
     &    0.285, 0.314, 0.380, 0.444, 0.500, 0.571, 0.673, 0.768, 0.853, 
     &    0.966, 1.097, 1.214, 1.333, 1.466, 1.591, 1.693, 1.796, 1.922, 
     &    2.037, 2.154, 2.264, 2.375, 2.508, 2.671, 2.847, 3.066, 3.417, 
     &    3.828, 4.204, 4.453, 4.599, 4.528, 4.284, 3.955, 3.678, 3.477, 
     &    3.346, 3.290, 3.251, 3.231, 3.226, 3.212, 3.192, 3.108, 3.033, 
     &    2.911, 2.798, 2.646, 2.508, 2.322, 2.130, 1.928, 1.757, 1.588, 
     &    1.417, 1.253, 1.109, 0.990, 0.888, 0.791, 0.678, 0.587, 0.524,
     &    0.464, 0.403, 0.357, 0.320, 0.290, 0.267, 0.242, 0.215, 0.182,
     &    0.160, 0.146, 0.128, 0.103, 0.087, 0.081, 0.071, 0.064, 0.000/
        DATA BETA / 
     &    467.0, 467.0, 400.0, 315.0, 379.0, 368.0, 475.0, 521.0, 531.0, 
     &    512.0, 442.0, 444.0, 430.0, 381.0, 335.0, 324.0, 296.0, 248.0, 
     &    215.0, 193.0, 158.0, 127.0, 101.0,  71.0,  31.0,  -6.0, -26.0, 
     &    -47.0, -63.0, -79.0, -88.0, -88.0, -87.0, -90.0, -98.0, -99.0,
     &   -109.0,-134.0,-160.0,-167.0,-164.0,-158.0,-153.0,-151.0,-156.0,
     &   -166.0,-168.0,-173.0,-170.0,-161.0,-145.0,-126.0,-108.0, -84.0, 
     &    -59.0, -29.0,   4.0,  41.0,  73.0,  97.0, 123.0, 159.0, 198.0, 
     &    220.0, 242.0, 256.0, 281.0, 311.0, 334.0, 319.0, 313.0, 321.0, 
     &    323.0, 310.0, 315.0, 320.0, 335.0, 361.0, 378.0, 373.0, 338.0, 
     &    319.0, 346.0, 322.0, 291.0, 290.0, 350.0, 371.0, 504.0, 504.0/
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      AFACT = AMTPTH(IPTH) * PREPTH(IPTH) / TEMPTH(IPTH) * 
     &            (273.0)**2 * RGAS / 101325.0 
      TFACT = 1.0/TEMREF - 1.0/TEMPTH(IPTH) 
C
      EVEN = .FALSE.
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
          BW = EW * B(IW) + DW * B(JW)
          BETAW = EW * BETA(IW) + DW * BETA(JW)
          CTMPTH = AFACT * BW * EXP ( BETAW * TFACT ) 
        ELSE
          CTMPTH = 0.0
        END IF
C
        JWID = 1 + IWD2 / 2 
        IF ( CTMPTH .NE. 0.0 ) ICTWID(JWID,IGSPTH(IPTH)) = 1
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
