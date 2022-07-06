C $Id: sub25w.for,v 1.2 1997/03/03 17:17:59 adarby Exp $
      SUBROUTINE SUB25W ( NWNO, ABSORP )
C
C VERSION
C      03-MAR-97  AD  Version 3.
C     18-DEC-96  AD  Original
C
C DESCRIPTION
C     Subtract absorption value at 25cm-1 from line centre 
C     Called by RFMFIN, RFMWID for H2O lines if CTMFLG .AND. .NOT. GL2FLG
C     Required for use with the CKD H2O continuum coefficients.
C
C ARGUMENTS
      INTEGER          NWNO          !  I  No. wavenumber points to be evaluted
      REAL             ABSORP(NWNO)  ! I/O Absorption coefficient
C
C GLOBAL CONSTANTS
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
C
C COMMON VARIABLES
      INCLUDE 'adjcom.inc' ! Path-Adjusted line data

C LOCAL VARIABLES
      INTEGER  IWNO   ! Wavenumber array counter
      REAL     ABS25W ! Absorption evaluated at 25cm-1 from line centre
C
C EXECUTABLE CODE --------------------------------------------------------------
C
C ABS25W is just the Lorentzian line shape evaluated at 25cm-1 from centre
C
      ABS25W = STRADJ * WIDADJ / PI / ( 625.0 + WIDADJ**2 )
C
C Assume CKD continuum allows for line contributions to absorption .LE. 0 
C
      DO IWNO = 1, NWNO
        ABSORP(IWNO) = ABSORP(IWNO) - ABS25W
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
