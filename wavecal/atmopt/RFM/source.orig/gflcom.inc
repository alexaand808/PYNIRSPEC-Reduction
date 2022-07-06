C GFLCOM.INC -------------------------------------------------------------------
      COMMON / GFLCOM / RESGFL, NGFL, LUNGFL, SPCGFL
C
C VERSION
C     20-MAR-01  AD  Add RESGFL 
C     09-DEC-97  AD  Original.
C
C DESCRIPTION
C     Common variables associated with RFM irregular grid (.GRD) files.
C     Updated by GRDCHK. NGFL initialised in INPGRD.
C
      DOUBLE PRECISION RESGFL      ! Resolution of GRD files (if convolving)
      INTEGER      NGFL            ! No. of GRD files to be used
      INTEGER      LUNGFL(MAXGFL)  ! LUN for each GRD file
      INTEGER      SPCGFL(MAXGFL)  ! Index of Spectral range in each GRD file
C
C ------------------------------------------------------------------------------
