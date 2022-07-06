C $Id: ilscom.inc,v 1.7 1997/03/03 17:18:19 adarby Exp $
C ILSCOM.INC -------------------------------------------------------------------
      COMMON / ILSCOM / WNLILS, WNUILS, FNCILS, PT1ILS, PTDILS, PT2ILS,
     &                  NILS, NPTILS 
C
C VERSION
C     27-APR-00  AD  Convert FNCILS,PT1ILS,PTDILS to DP. Add PT2ILS
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Instrument Lineshape functions.
C     Loaded and checked by FILILS.
C
      DOUBLE PRECISION WNLILS(MAXILS) ! Lower wavenumber for ILS applicability
      DOUBLE PRECISION WNUILS(MAXILS) ! Upper wavenumber for ILS applicability
      DOUBLE PRECISION FNCILS(MAXILP,MAXILS)  ! Tabulated ILS functions
      DOUBLE PRECISION PT1ILS(MAXILS) ! Initial Pts for each ILS function
      DOUBLE PRECISION PT2ILS(MAXILS) ! Upper Pts for each ILS function
      DOUBLE PRECISION PTDILS(MAXILS) ! Pt.spacing for each ILS function
      INTEGER          NILS           ! No. ILS functions stored
      INTEGER          NPTILS(MAXILS) ! No. tabulation pts for each ILS Fn.
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
