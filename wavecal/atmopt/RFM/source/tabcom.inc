C TABCOM.INC -------------------------------------------------------------------
      COMMON / TABCOM / NLPTAB, NTMTAB, OFFTAB,
     &                  LP1TAB, LPDTAB, TM1TAB, TMDTAB
C
C VERSION
C     28-APR-99  AD  Add OFFTAB
C     22-OCT-97  AD  Original.
C
C DESCRIPTION
C     Tabulation axes for tabulated absorption coefficients (TAB option).
C     Loaded by CHKTAB, used by TABPTH.
C
      INTEGER NLPTAB                !  No. pressure points
      INTEGER NTMTAB                !  No. temperature points
      LOGICAL OFFTAB                !  T=apply temperature offset profile
      REAL    LP1TAB                !  Ln(p/mb) of first point
      REAL    LPDTAB                !  Ln(p/mb) increment
      REAL    TM1TAB                !  Temperature [K] of first point
      REAL    TMDTAB                !  Temperature increment [K]
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
