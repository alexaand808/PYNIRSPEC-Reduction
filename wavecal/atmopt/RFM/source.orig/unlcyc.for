C $Id: unlcyc.for,v 1.9 1997/03/03 17:18:02 adarby Exp $
      SUBROUTINE UNLCYC ( ICYC )
C
C VERSION
C      03-MAR-97  AD  Version 3.
C      01-OCT-96  AD  Version 2.
C      01-OCT-96  AD  Unload BSLQ from cyclic buffer
C      21-SEP-96  AD  Rename IUSVIB,ILSVIB to IUSNTE,ILSNTE
C      01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Unload required line from Cyclic buffer.
C     Called by RFMFIN.
C
      IMPLICIT NONE
C
C  ARGUMENTS      
      INTEGER      ICYC    !  I  Index of line in cyclic buffer
C
C GLOBAL CONSTANT 
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'cyccom.inc' ! Cyclic line data buffers
      INCLUDE 'hitcom.inc' ! HITRAN line database variables
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      WNUM   = WNOCYC(ICYC) 
      ISO    = ISOCYC(ICYC)
      IDGAS  = IDGCYC(ICYC)
      STREN  = STRCYC(ICYC)
      ABROAD = ABRCYC(ICYC) 
      SBROAD = SBRCYC(ICYC) 
      ELS    = ELSCYC(ICYC) 
      ABCOEF = ABCCYC(ICYC) 
      TSP    = TSPCYC(ICYC) 
      IUSGQ  = IUSCYC(ICYC) 
      ILSGQ  = ILSCYC(ICYC) 
      IUSNTE = IUVCYC(ICYC)
      ILSNTE = ILVCYC(ICYC)
      BSLQ   = BLQCYC(ICYC)
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
