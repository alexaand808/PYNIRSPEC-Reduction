C $Id: adjcom.inc,v 1.8 1997/03/03 17:18:07 adarby Exp $
C ADJCOM.INC -------------------------------------------------------------------
      COMMON / ADJCOM / WNOADJ, ANTADJ, CNTADJ,
     &                  WIDADJ, STRADJ, DOPADJ, YMXADJ
C
C VERSION
C      21-OCT-97  AD  Convert ANTADJ, CNTADJ to DOUBLE PRECISION 
C      03-MAR-97  AD  Version 3.
C      01-OCT-96  AD  Version 2.
C      01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Common variables associated with path-adjusted line data 
C     Set by module ADJUST. 
C
      DOUBLE PRECISION WNOADJ ! Wavenumber [cm-1]
      DOUBLE PRECISION ANTADJ ! Non-lte factors for k abs
      DOUBLE PRECISION CNTADJ ! Non-lte factors for c abs
      REAL             WIDADJ ! half width [cm-1]
      REAL             STRADJ ! line strength  [cm-1]
      REAL             DOPADJ ! Doppler half width [cm-1]
      REAL             YMXADJ ! Y mixing coeffient
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
