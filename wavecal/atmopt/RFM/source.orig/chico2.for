C $Id: chico2.for,v 1.8 1997/03/03 17:17:01 adarby Exp $
      SUBROUTINE CHICO2 ( IPTH, NWNO, DWNO, CHIFAC )
C
C VERSION
C      22-JAN-98  AD  (Standard F77) Put DATA statements *after* declarations
C      03-MAR-97  AD  Version 3.
C      01-OCT-96  AD  Version 2.
C      01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Calculate Chi-Factor for CO2 sub-Lorentzian lineshape.
C     Called by SUBSHP for CO2 lines.
C
C REFERENCES
C     C.Cousin,R.Le Doucen,C.Boulet and A.Henry, (1985)
C     Temperature dependence of the absorption in the region beyond the 4.3 
C     mic band head of CO2. 2: N2 and O2 broadening.
C     Appl.Opt. 24 (22) 3899-3907
C
C     Pure CO2 chi factor: Appl.Opt. 24 (6) 897-906 (1985)
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER          IPTH          !  I  Index of path
      INTEGER          NWNO          !  I  No. wavenumber points to be evaluted
      DOUBLE PRECISION DWNO(NWNO)    !  I  Array of Wavenumbers [/cm]
      REAL             CHIFAC(NWNO)  !  O  Chi-factor
C
C GLOBAL CONSTANTS
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants      
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'adjcom.inc' ! Path-Adjusted line data
      INCLUDE 'pthcom.inc' ! RFM path data
C
C LOCAL CONSTANTS
      REAL          DFRQ         ! Wavenumber increment [cm-1] for chi table
        PARAMETER ( DFRQ = 0.5 )
      INTEGER       NFRQ         ! No.wavenumber points in chi tabulation
        PARAMETER ( NFRQ = 503 )
      REAL          FRQMAX    ! Max.wavenumber [cm-1] in chi tabulation
        PARAMETER ( FRQMAX = ( NFRQ + 1 ) * DFRQ )
      INTEGER       NTEM         ! No.temperatures in chi tabulation
        PARAMETER ( NTEM = 5 )
C
C LOCAL VARIABLES
      INTEGER  IFRQ,JFRQ   ! Lower,upper indices for wavenumber interpolation
      INTEGER  ITEM,JTEM   ! Lower,upper indices for.broadening interpolation
      INTEGER  IWNO        ! Wavenumber array counter
      REAL     DFOR,DSLF   ! Weights for foreign-,self-broadening terms
      REAL     FFI,FFJ     ! Weights for lower,upper wavenumber interpolation
      REAL     FRQ         ! Distance [cm-1] from line centre
      REAL     TEM         ! Temperature [K] used for interpolation
      REAL     TFI, TFJ    ! Weights for lower,upper for.broadened temperatures
      REAL     TSI, TSJ    ! Weights for lower,upper self broadened temperatures
      REAL     CHITEM(NTEM) ! Temperatures (1:3=foreign, 4:5=self)
        SAVE   CHITEM
      REAL     CHI(NFRQ,NTEM) ! Chi tabulation for CO2
        SAVE   CHI
C
C DATA STATEMENTS
      INCLUDE 'chidat.inc' ! Data statements for CHI values
        DATA   CHITEM / 193.0, 238.0, 296.0, 218.0, 296.0 /
C
C EXECUTABLE CODE --------------------------------------------------------------
C
C Foreign and self broadening fractions 
C
      DFOR = ( PREPTH(IPTH) - PPAPTH(IPTH) ) / PREPTH(IPTH)
      DSLF = 1.0 - DFOR
C
C Temperature interpolation (foreign broadened)
C 
      TEM = MAX ( MIN ( TEMPTH(IPTH), CHITEM(3) ), CHITEM(1) ) 
      IF ( TEM .GE. CHITEM(2) ) THEN
        ITEM = 2
      ELSE
        ITEM = 1
      END IF
      JTEM = ITEM + 1
      TFJ = ( TEM - CHITEM(ITEM) ) / ( CHITEM(JTEM) - CHITEM(ITEM) )
      TFI = ( 1.0 - TFJ ) * DFOR
      TFJ = TFJ * DFOR
C
C Temperature interpolation (self broadened)
C
      TEM = MAX ( MIN ( TEMPTH(IPTH), CHITEM(5) ), CHITEM(4) ) 
      TSJ = ( TEM - CHITEM(4) ) / ( CHITEM(5) - CHITEM(4) )
      TSI = ( 1.0 - TSJ ) * DSLF
      TSJ = TSJ * DSLF
C
      DO IWNO = 1, NWNO 
C
C Frequency interpolation
C 
        FRQ = MIN ( SNGL ( ABS ( DWNO(IWNO) - WNOADJ ) ), FRQMAX )  
        IFRQ = MIN ( NFRQ - 1, INT ( FRQ / DFRQ ) + 1 )
        JFRQ = IFRQ + 1
        FFJ = FRQ - IFRQ * DFRQ 
        FFI = 1.0 - FFJ
C
        CHIFAC(IWNO) = 
     &      TFI * ( FFI * CHI(IFRQ,ITEM) + FFJ * CHI(JFRQ,ITEM) )
     &    + TFJ * ( FFI * CHI(IFRQ,JTEM) + FFJ * CHI(JFRQ,JTEM) )
     &    + TSI * ( FFI * CHI(IFRQ,4)    + FFJ * CHI(JFRQ,4) )
     &    + TSJ * ( FFI * CHI(IFRQ,5)    + FFJ * CHI(JFRQ,5) )
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
