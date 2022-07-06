C $Id: lorshp.for,v 1.9 1997/03/03 17:17:30 adarby Exp $
      SUBROUTINE LORSHP ( NWNO, DWNO, ABSORP )
C
C VERSION
C      03-MAR-97  AD  Version 3.
C      01-OCT-96  AD  Version 2.
C      01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Calculate Lorentz Line shape.
C     Called by RFMFIN and RFMWID.
C     Uses path-adjusted line data in /ADJCOM/
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER          NWNO          !  I  No. wavenumber points to be evaluted
      DOUBLE PRECISION DWNO(NWNO)    !  I  Array of Wavenumbers [/cm]
      REAL             ABSORP(NWNO)  !  O  Absorption 
C
C GLOBAL CONSTANTS
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
C
C COMMON VARIABLES
      INCLUDE 'adjcom.inc' ! Path-Adjusted line data
C
C LOCAL VARIABLES
      INTEGER  IWNO ! Wavenumber array counter
      REAL     TOP  ! Numerator of Lorentz expression
      REAL     W2   ! Width^2
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      TOP = STRADJ * WIDADJ / PI                   
      W2 = WIDADJ**2
C
      DO IWNO = 1, NWNO 
        ABSORP(IWNO) = TOP / 
     &    ( SNGL ( ( DWNO(IWNO) - WNOADJ )**2 ) + W2 )
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
