C $Id: setnte.for,v 1.5 1997/03/03 17:17:54 adarby Exp $
      SUBROUTINE SETNTE
C
C VERSION
C      03-MAR-97  AD  Version 3.
C      01-OCT-96  AD  Version 2.
C      20-SEP-96  AD  Original.
C
C DESCRIPTION
C     Set IUSNTE and ILSNTE in /HITCOM/ if line affected by non-LTE.
C     Called by REAHIT for each line if NTEFLG enabled
C    
      IMPLICIT NONE
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'hitcom.inc' ! HITRAN line database variables
      INCLUDE 'ntecom.inc' ! Non-LTE data
C
C LOCAL VARIABLES
      INTEGER INTE         ! Counter for Non-LTE datasets
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      IUSNTE = 0
      ILSNTE = 0
      DO INTE = 1, NNTE        ! Loop over all non-LTE vibrational levels
        IF ( IDGNTE(INTE) .EQ. IDGAS .AND.
     &       ISONTE(INTE) .EQ. ISO        ) THEN     ! Matches current isotope
          IF ( IGQNTE(INTE) .EQ. IUSGQ ) THEN        ! Upper Glob.Quant.match
            IUSNTE = INTE
          ELSE IF ( IGQNTE(INTE) .EQ. ILSGQ ) THEN   ! Lower Glob.Quant.match
            ILSNTE = INTE
          END IF
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
