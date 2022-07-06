C $Id: xflcom.inc,v 1.9 1997/03/03 17:18:30 adarby Exp $
C XFLCOM.INC -------------------------------------------------------------------
      COMMON / XFLCOM / NXFL, LUNXFL, GASXFL
C
C VERSION
C      03-MAR-97  AD  Version 3.
C      16-JAN-97  AD  Add GASXFL
C      01-OCT-96  AD  Version 2.
C      17-SEP-96  AD  Remove redundant IRCXFL
C      01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Common variables associated with RFM cross-section data files
C     Loaded by INPXFL.
C
      INTEGER      NXFL            ! No. of X/S data files to be read
      INTEGER      LUNXFL(MAXXFL)  ! LUN for use with all files
      LOGICAL      GASXFL(MAXGAS)  ! T=Some X/S data for gas in at least 1 file
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
