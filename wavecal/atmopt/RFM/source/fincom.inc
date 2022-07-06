C FINCOM.INC ------------------------------------------------------------------
      COMMON / FINCOM / WNOFIN, WNDFIN, WNLFIN, WNUFIN, WN1FIN, 
     &                  WN2FIN, WNRFIN,
     &                  NFIN, NOMFIN, IXSFIN, ABSFIN, CNTFIN
C
C VERSION
C     17-AUG-06  AD  Add WNDFIN 
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     17-SEP-96  AD  Add NOMFIN
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Common variables associated with fine mesh calculation grid. 
C     Assumes NFIN fine intervals, indexed 1:NWID for particular IWID.
C
      DOUBLE PRECISION WNOFIN(MAXFIN)   ! Fine mesh points [/cm] 
      DOUBLE PRECISION WNDFIN           ! Offset of fine mesh rel.to widemesh
      DOUBLE PRECISION WNLFIN, WNUFIN   ! Lower, Upper Wnos for fine mesh calc.
      DOUBLE PRECISION WN1FIN, WN2FIN   ! Lower, Upper Wnos at widemesh bounds
      DOUBLE PRECISION WNRFIN           ! Fine Mesh resolution
      INTEGER          NFIN             ! Actual No. of fine mesh grid points
      INTEGER          NOMFIN           ! Nominal No.pts/cm-1 for fine mesh
      INTEGER          IXSFIN(MAXGAS)   ! 1 = X/S data used for current itvl
      REAL     ABSFIN(MAXFIN,MAXCLC) ! Fine mesh absorption 
      REAL     CNTFIN(MAXFIN,MAXCLC) ! Fine mesh path non-LTE factor 
C
C -----------------------------------------------------------------------------
