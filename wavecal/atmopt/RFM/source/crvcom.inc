C $Id: crvcom.inc,v 1.5 1997/03/03 17:18:11 adarby Exp $
C CRVCOM.INC -------------------------------------------------------------------
      COMMON / CRVCOM / RADCRV
C
C VERSION
C     29-DEC-99  AD  Convert RADCRV from SP to DP
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     17-SEP-96  AD  Original.
C
C DESCRIPTION
C     Local radius of curvature.
C     (for inclusion in ray-tracing calculation).
C     Initially set to REARTH in RFMDEF, modified by INPCRV if CRVFLG enabled.
C
      DOUBLE PRECISION RADCRV  ! Radius of curvature [km]
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
