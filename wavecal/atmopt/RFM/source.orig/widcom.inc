C $Id: widcom.inc,v 1.9 1997/03/03 17:18:30 adarby Exp $
C WIDCOM.INC -------------------------------------------------------------------
      COMMON / WIDCOM / WN1WID, WN2WID, WNLWID, WNUWID, DELWID, 
     &                  WNOWID, WNOWD2, 
     &                  NWID, NILWID, NOLWID, NLNWID, ICTWID, 
     &                  ABSWID, CNTWID
C
C VERSION
C     13-APR-00  AD  Add NILWID, NOLWID
C     09-AUG-99  AD  Change WNOWID from Real to D.P.
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Common variables associated with RFM widemesh calculations.
C     Assumes NWID widemesh intervals, indexed 1:NWID,
C     NWID+1 widemesh grid points, indexed 0:NWID.
C
      DOUBLE PRECISION WN1WID      ! WNO at lower limit of widemesh grid
      DOUBLE PRECISION WN2WID      ! WNO at upper limit of widemesh grid
      DOUBLE PRECISION WNLWID      ! lowest WNO required for any lines 
      DOUBLE PRECISION WNUWID      ! highest WNO required for any lines 
      DOUBLE PRECISION DELWID      ! Width of Widemesh intervals
      DOUBLE PRECISION WNOWID(0:MAXWID) ! Wide mesh points [/cm] 
      DOUBLE PRECISION WNOWD2(0:MAXWD2) !  WNO at each widemesh parab. point
C
      INTEGER  NWID                    ! No. of wide mesh intervals 
      INTEGER  NILWID(MAXGAS)          ! No.lines within spec range for each gas
      INTEGER  NOLWID(MAXGAS)          ! No.lines outside spec range for each gas
      INTEGER  NLNWID(MAXWID,MAXGAS)   ! No.lines included for each w.m. itvl
      INTEGER  ICTWID(MAXWID,MAXGAS)   ! 1 = continuum data used for wm itvl
      REAL     ABSWID(3,MAXWID,MAXCLC) ! Widepass abs. by parabolic point.
      REAL     CNTWID(3,MAXWID,MAXCLC) ! Widepass nonlte factor absorption by
C                                 parabolic point, wide mesh, and path.
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
