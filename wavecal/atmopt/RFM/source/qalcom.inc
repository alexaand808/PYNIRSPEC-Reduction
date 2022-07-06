C QALCOM.INC -------------------------------------------------------------------
C
      COMMON / QALCOM / NBND, IUSQAL, ILSQAL, USEQAL
C
C VERSION
C      02-OCT-97  AD  Add MAXISO dimension. Remove ISOQAL
C      21-JUL-97  AD  Original.
C
C DESCRIPTION
C     List of band/isotope qualifiers for line data.
C     Variables only valid if USEQAL for gas is set TRUE.
C     USEQAL set by GASCHK, other variables set by GASQAL.
C     Used in module REAHIT.
C     Note: since the number of isotopes is relatively small (compared to the
C     number of possible GQ numbers) it is possible to pre-assign a logical
C     flag to each isotope/gas combination. However, for the lower and upper
C     Global Quantum numbers it is necessary to construct a list of those that
C     are required. ILS,IUS values of 0 imply use any HITRAN index value.
C
      INTEGER     NBND(MAXISO,MAXGAS)         ! No. band qualifiers
      INTEGER   ILSQAL(MAXBND,MAXISO,MAXGAS)  ! Selected band lower GQ numbers
      INTEGER   IUSQAL(MAXBND,MAXISO,MAXGAS)  ! Selected band upper GQ numbers
      LOGICAL   USEQAL(MAXGAS)                ! T=at least one qualifier set
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
