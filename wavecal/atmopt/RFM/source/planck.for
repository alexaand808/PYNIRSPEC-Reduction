C $Id: planck.for,v 1.9 1997/03/03 17:17:37 adarby Exp $
      SUBROUTINE PLANCK ( TEM, NWNO, WNOLST, EMS )
C
C VERSION
C     26-MAY-00  AD  Convert EMS from SP to DP
C     09-AUG-99  AD  Use local D.P. DARGMN instead of S.P. ARGMIN 
C                    Make explicit EMS = SNGL( ) 
C     15-DEC-98  AD  Avoid for zero divide if Wno=0.0
C     29-JUL-98  AD  Check for overflows for low temperatures
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Calculate Planck emission spectrum at specified temperature
C     Called by RFMRAD.
C
      IMPLICIT NONE
C
C ARGUMENTS 
      REAL             TEM          ! Temperature [K]
      INTEGER          NWNO         ! No.of wavenumber points
      DOUBLE PRECISION WNOLST(NWNO) !  List of wavenumbers [/cm]
      DOUBLE PRECISION EMS(NWNO)    !  Planck emission 
C
C GLOBAL CONSTANTS
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
C
C LOCAL CONSTANTS
      DOUBLE PRECISION DARGMN   ! Minimum value for argument of exponent.
        PARAMETER ( DARGMN = 1.0D-38 )
C
C LOCAL VARIABLES
      INTEGER          IWNO  ! Wavenumber counter
      DOUBLE PRECISION ARG   ! Intermediate calculation
      DOUBLE PRECISION EFACT ! Term for exponential 
C
C EXECUTABLE CODE --------------------------------------------------------------
C
C Limit temperature to a minimum value of 1e-38/Wmax/C2 so that argument of 
C exponent is limited to magnitude of 1e-38
      EFACT = -C2 / MAX ( DBLE(TEM), DARGMN*WNOLST(NWNO)*C2 )
      DO IWNO = 1, NWNO
        ARG = EXP ( EFACT * WNOLST(IWNO) )
        EMS(IWNO) = C1 * WNOLST(IWNO)**3 * ARG 
     &              / MAX ( 1.0D0 - ARG, DARGMN )
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
