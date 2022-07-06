C $Id: xsccom.inc,v 1.7 1997/03/03 17:18:31 adarby Exp $
C XSCCOM.INC -------------------------------------------------------------------
      COMMON / XSCCOM / WN1XSC, WN2XSC, 
     &                  LXSC, NXSC, IGSXSC, IFLXSC, IRCXSC
C
C VERSION
C     17-APR-01  AD  Add LXSC. Move DP variables to front so remove XSCDBL
C     31-AUG-99  AD  Convert WN1XSC, WN2XSC from REAL to D.P. and put in
C                    separate common block XSCDBL
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Common variables associated with RFM X-Secn data.
C
      DOUBLE PRECISION WN1XSC(MAXXSC)  ! Lower Wno. [/cm] of each X/S data set
      DOUBLE PRECISION WN2XSC(MAXXSC)  ! Upper Wno [/cm] of each X/S data set
      INTEGER      LXSC            ! Index of currently loaded X/S data set
      INTEGER      NXSC            ! No. of line data files to be read
      INTEGER      IGSXSC(MAXXSC)  ! Index of Gas for each X/S data set
      INTEGER      IFLXSC(MAXXSC)  ! Index of X/S File containing X/s data set
      INTEGER      IRCXSC(MAXXSC)  ! Starting Record#
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
