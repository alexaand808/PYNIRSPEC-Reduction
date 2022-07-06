C $Id: spcfin.for,v 1.8 1997/03/03 17:17:56 adarby Exp $
      SUBROUTINE SPCFIN ( ISPC )
C
C VERSION
C     04-JAN-99  AD  Allow for AVGFLG as well a ILSFLG
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     17-SEP-96  AD  Remove redundant RFMCON.INC
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Set Fine Mesh grid parameters from spectral range/resln.
C     Called by INPCHK for each required spectral range.
C     WNRSPC is resolution parameter, either as no.pts/wavenumber (WNRSPC>1)
C       or as actual resolution (<1).
C     If not using ILS convolution, WNRSPC refers to finemesh spacing directly.
C     If using ILS convolution, WNRSPC refers to final output spectra resolution
C       so set finemesh to closest sub-multiple of WNRSPC to NOMFIN pts/wno.
C 
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      ISPC    !  I  Index of spectral range 
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
      INCLUDE 'widcom.inc' ! Wide mesh data
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Calculation of WNRFIN should match that in ILSCHK.FOR
      IF ( ILSFLG .OR. AVGFLG ) THEN                     ! Using ILS 
        IF ( WNRSPC(ISPC) .GT. 1.0 ) THEN                ! N.Pts/wavenumber
          WNRFIN = 1.D0 / ( WNRSPC(ISPC) * NINT ( NOMFIN /WNRSPC(ISPC)))
        ELSE                                             ! Resolution
          WNRFIN = WNRSPC(ISPC) / NINT ( NOMFIN * WNRSPC(ISPC) )
        END IF
      ELSE                                               ! Not using ILS
        IF ( WNRSPC(ISPC) .GT. 1.0 ) THEN                ! N.Pts/wavenumber
          WNRFIN = 1.D0 / WNRSPC(ISPC)
        ELSE                                             ! Resolution
          WNRFIN = WNRSPC(ISPC)
        END IF
      END IF
C
      NFIN = NINT ( DELWID / WNRFIN )
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
