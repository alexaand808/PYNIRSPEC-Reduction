C $Id: spccom.inc,v 1.10 1997/03/03 17:18:29 adarby Exp $
C SPCCOM.INC -------------------------------------------------------------------
      COMMON / SPCCOM / WNLSPC, WNUSPC, WNRSPC, WMNSPC, WMXSPC, WMDSPC,
     &                  NSPC, ILSSPC, LABSPC
C
C VERSION
C      18-MAR-01  AD  Remover RMXSPC
C      23-APR-99  AD  Change LABSPC from C*6 to C*8
C      03-MAR-97  AD  Version 3.
C      01-OCT-96  AD  Version 2.
C      18-SEP-96  AD  Add RMXSPC
C      01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Required RFM spectral range & resolution.
C     Loaded and checked by INPSPC
C
      DOUBLE PRECISION WNLSPC(MAXSPC) ! Specified Lower wavenumber limit [/cm]
      DOUBLE PRECISION WNUSPC(MAXSPC) ! Specified Upper wavenumber limit [/cm]
      DOUBLE PRECISION WNRSPC(MAXSPC) ! Calcuated Wavenumber resolution [/cm]
      DOUBLE PRECISION WMNSPC         ! Lower Wavenumber reqd for any range
      DOUBLE PRECISION WMXSPC         ! Upper Wavenumber reqd for any range
      DOUBLE PRECISION WMDSPC         ! Mid-Wavenumber 0.5*(WMNSPC+WMXSPC)
      INTEGER          NSPC           ! No.tabulated Spc.ranges
      INTEGER          ILSSPC(MAXSPC) ! Index of reqd ILS Fn. (set by ILSCHK)
      CHARACTER*8      LABSPC(MAXSPC) ! Codes for tabulated ranges
C
C ------------------------------------------------------------------------------
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
