C $Id: cyccom.inc,v 1.10 1997/03/03 17:18:11 adarby Exp $
C CYCCOM.INC -------------------------------------------------------------------
      COMMON / CYCCOM / WNOCYC, STRCYC, ABRCYC, 
     &                  SBRCYC, ELSCYC, ABCCYC, TSPCYC, 
     &                  NLIN, ICYC1, IDGCYC, ISOCYC, IUSCYC, ILSCYC, 
     &                  IUVCYC, ILVCYC, BLQCYC
C
C VERSION
C      03-MAR-97  AD  Version 3.
C      01-OCT-96  AD  Version 2.
C      01-OCT-96  AD  Add BLQCYC to save BSLQ
C      01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Common variables associated with Cyclic line data buffers 
C
      DOUBLE PRECISION WNOCYC(MAXCYC) ! Cyclic buf of line wnos. [cm-1]
      REAL      STRCYC(MAXCYC) ! Line strength  [cm-1./(kg.moles.cm-2)]@296K
      REAL      ABRCYC(MAXCYC) ! Air-broad halfwidth  (HWHM) [cm-1/atm] @ 296K
      REAL      SBRCYC(MAXCYC) ! Self-broad halfwidth (HWHM) [cm-1/atm] @ 296K.
      REAL      ELSCYC(MAXCYC) ! Lower-state energy [cm-1]
      REAL      ABCCYC(MAXCYC) ! Coeff.of temp.dependance of air-broadened HW
      REAL      TSPCYC(MAXCYC) ! Transition shift due to pressure
      INTEGER   NLIN           ! No.lines currently stored
      INTEGER   ICYC1          ! Index for lowest wavenumber line
      INTEGER   IDGCYC(MAXCYC) ! HITRAN Gas ID
      INTEGER   ISOCYC(MAXCYC) ! Isotope Number (1=most abundant,2=2nd,3 etc)
      INTEGER   IUSCYC(MAXCYC) ! Upper state global quanta index.
      INTEGER   ILSCYC(MAXCYC) ! Lower state global quanta index
      INTEGER   IUVCYC(MAXCYC) ! Index of Vib.Temp profile affecting upper level
      INTEGER   ILVCYC(MAXCYC) ! Index of Vib.Temp profile affecting lower level
      CHARACTER*9 BLQCYC(MAXCYC) ! Lower State local quanta
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
