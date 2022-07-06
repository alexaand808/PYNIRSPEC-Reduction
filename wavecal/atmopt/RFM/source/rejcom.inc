C $Id: rejcom.inc,v 1.7 1997/03/03 17:18:26 adarby Exp $
C REJCOM.INC -------------------------------------------------------------------
      COMMON / REJCOM / WIDREJ, FINREJ
C
C VERSION
C      03-MAR-97  AD  Version 3.
C      01-OCT-96  AD  Version 2.
C      01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Minimum line strength limits
C     (for inclusion of line in wide and fine mesh calculations)
C     Initially set to 0.0 in RFMDEF, modified by INPREJ if REJFLG enabled.
C     Typical values are 1.0e-7 (WIDREJ) and 1.0e-8 (FINREJ)
C
      REAL     WIDREJ  !  Min.Line strength for wide mesh calculation
      REAL     FINREJ  !  Min.Line strength for fine mesh calculation
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
