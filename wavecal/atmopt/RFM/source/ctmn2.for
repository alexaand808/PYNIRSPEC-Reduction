C $Id: ctmn2.for,v 1.12 1997/03/03 17:17:05 adarby Exp $
      SUBROUTINE CTMN2 ( IPTH ) 
C
C VERSION
C      22-JAN-98  AD  (Standard F77) Put DATA statements *after* declarations
C      11-NOV-97  AD  Correction to comment on Upper Wavenumber
C      03-MAR-97  AD  Version 3.
C      26-FEB-97  AD  Replace Menoux continuum with Lafferty
C      01-OCT-96  AD  Version 2.
C      01-SEP-96  AD  Version 1.
C 
C DESCRIPTION    
C     Calculate N2 continuum absorption across entire widemesh grid.
C     Called by RFMCTM for each Path containing N2.
C     The continuum absorption/length s(/cm) is derived from
C
C        s(p,T) = v (p/po)^2 (To/T)^2 (v+(1-v)E(T)) Bx10^-6 exp(beta(1/Tr-1/T))
C 
C     where B and beta are tabulated functions of wavenumber, To=273K, po is
C     std.atmospheric pressure, V is the vmr of N2, Tr=296K, and E(T) is the
C     (wavenumber-independent) relative efficiency of O2 broadening to N2,
C
C        E(T) = 1.294 - 0.4545 (T/Tr)
C
C     B is tabulated in (/cm /Am^2) where 1 Amagat is the ratio of the molar
C     density to the value at STP (hence the p/po To/T terms). 
C     Therefore path absorption k = s x, where x is the path length.
C     For consistency with other RFM modules, this is rearranged to be a product
C     of the path absorber amount u (kmol/cm^2), given by u = (vpx/(RT))
C     This gives a calculation 
C
C        k(p,T) = u p To^2 R / (po T) (v-(1-v)E(T)) B 10^6 exp(beta(1/Tr - 1/T))
C
C     where p is in atm, R is J/kmol/K, po is Pa. Note that the extra 10^6 
C     scaling term for the cm^3/m^3 conversion.
C
C REFERENCES
C     LAFFERTY, W.J, A.M. SOLODOV, A. WEBER, W.B. OLSON and J-M HARTMANN
C       Infrared collision-induced absorption by N2 near 4.3um for 
C       atmospheric applications: measurements and empirical modelling.
C     Applied Optics, 35, 5911-5917 (1996).
C     Note that the tabulation in the paper extends from 2125-2600, but here
C     values have been added at 2120 and 2605 (0 for B, end values for BETA)
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
        PARAMETER ( WNOLOW = 2120.0 ) ! Upper wavenumber = 2605.0
      REAL          DELWNO            ! Wavenumber interval of Continm.data[/cm]
        PARAMETER ( DELWNO = 5.0 )
      INTEGER       NWNO              ! No.continuum wavenumber tabulation pts.
        PARAMETER ( NWNO = 98 )
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
      REAL    EFACT  ! Factor to allow for Foreign (ie O2) broadening 
      REAL    TFACT  ! Temperature scaling for BETA
      REAL    VMR    ! N2 VMR for path
      REAL    XW     ! Position of half-WM point on CO2 continuum tabulation axis
      REAL    B(NWNO)  ! [cm-1 Amagat-2] 
         SAVE B
      REAL   BETA(NWNO) !  [K]
         SAVE BETA
C
C DATA STATEMENTS
        DATA B / 
     & 0.000E-07, 0.445E-07, 0.522E-07, 0.646E-07, 0.775E-07, 0.903E-07,    
     & 0.106E-06, 0.121E-06, 0.137E-06, 0.157E-06, 0.175E-06, 0.201E-06,    
     & 0.230E-06, 0.259E-06, 0.295E-06, 0.326E-06, 0.366E-06, 0.405E-06,    
     & 0.447E-06, 0.492E-06, 0.534E-06, 0.584E-06, 0.624E-06, 0.667E-06,    
     & 0.714E-06, 0.726E-06, 0.754E-06, 0.784E-06, 0.809E-06, 0.842E-06,   
     & 0.862E-06, 0.887E-06, 0.911E-06, 0.936E-06, 0.976E-06, 0.103E-05,   
     & 0.111E-05, 0.123E-05, 0.139E-05, 0.161E-05, 0.176E-05, 0.194E-05,   
     & 0.197E-05, 0.187E-05, 0.175E-05, 0.156E-05, 0.142E-05, 0.135E-05,   
     & 0.132E-05, 0.129E-05, 0.129E-05, 0.129E-05, 0.130E-05, 0.132E-05,   
     & 0.133E-05, 0.134E-05, 0.135E-05, 0.133E-05, 0.131E-05, 0.129E-05,   
     & 0.124E-05, 0.120E-05, 0.116E-05, 0.110E-05, 0.104E-05, 0.996E-06,     
     & 0.938E-06, 0.863E-06, 0.798E-06, 0.726E-06, 0.655E-06, 0.594E-06,    
     & 0.535E-06, 0.474E-06, 0.424E-06, 0.377E-06, 0.333E-06, 0.296E-06,    
     & 0.263E-06, 0.234E-06, 0.208E-06, 0.185E-06, 0.167E-06, 0.147E-06,    
     & 0.132E-06, 0.120E-06, 0.109E-06, 0.985E-07, 0.908E-07, 0.818E-07,    
     & 0.756E-07, 0.685E-07, 0.614E-07, 0.583E-07, 0.577E-07, 0.500E-07,    
     & 0.432E-07, 0.000E-07/
        DATA BETA / 
     &   802.0, 802.0, 761.0, 722.0, 679.0, 646.0, 609.0, 562.0, 511.0, 
     &   472.0, 436.0, 406.0, 377.0, 355.0, 338.0, 319.0, 299.0, 278.0,
     &   255.0, 233.0, 208.0, 184.0, 149.0, 107.0,  66.0,  25.0, -13.0,
     &   -49.0, -82.0,-104.0,-119.0,-130.0,-139.0,-144.0,-146.0,-146.0,
     &  -147.0,-148.0,-150.0,-153.0,-160.0,-169.0,-181.0,-189.0,-195.0,
     &  -200.0,-205.0,-209.0,-211.0,-210.0,-210.0,-209.0,-205.0,-199.0,
     &  -190.0,-180.0,-168.0,-157.0,-143.0,-126.0,-108.0, -89.0, -63.0,
     &   -32.0,   1.0,  35.0,  65.0,  95.0, 121.0, 141.0, 152.0, 161.0,
     &   164.0, 164.0, 161.0, 155.0, 148.0, 143.0, 137.0, 133.0, 131.0,
     &   133.0, 139.0, 150.0, 165.0, 187.0, 213.0, 248.0, 284.0, 321.0,
     &   372.0, 449.0, 514.0, 569.0, 609.0, 642.0, 673.0, 673.0 /
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      AFACT = AMTPTH(IPTH) * PREPTH(IPTH) / TEMPTH(IPTH) * 
     &            (273.0)**2 * RGAS / 101325.0 * 1.0E6
      TFACT = 1.0/TEMREF - 1.0/TEMPTH(IPTH) 
      VMR   = PPAPTH(IPTH)/PREPTH(IPTH)
      EFACT = VMR + (1.0-VMR) * (1.294 - 0.4545 * TEMPTH(IPTH) / TEMREF)
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
          CTMPTH = AFACT * BW * EFACT * EXP ( BETAW * TFACT ) 
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
