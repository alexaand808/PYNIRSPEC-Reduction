C GRDCOM.INC ------------------------------------------------------------------
      COMMON / GRDCOM / WNOGRD, WGTGRD, 
     &                  NGRD, NOPGRD, IOPGRD, TWOGRD, OUTGRD, FNCGRD
C
C VERSION
C     03-JAN-03  AD  Add NOPGRD,IOPGRD,WGTGRD
C     12-DEC-97  AD  Original.
C
C DESCRIPTION
C     Common variables associated with RFM Irregular grid data.
C     These data are loaded for each spectral range by SPCGRD and GRDILS
C
      INTEGER          NGRD           ! No. wavenumber points in grid
      INTEGER          NOPGRD(MAXGRD) ! No.o/p grid points for each irr.grd.pt
      INTEGER   IOPGRD(MAXOPG,MAXGRD) ! Index of o/p grid points 
      LOGICAL          TWOGRD         ! T=grid uses two-point interpolation
      LOGICAL          OUTGRD         ! T=map directly to o/p, F=interp.first
      DOUBLE PRECISION WNOGRD(MAXGRD) ! List of wavenumbers [/cm] 
      DOUBLE PRECISION WGTGRD(MAXOPG,MAXGRD) ! Weight to each o/p grid point
      CHARACTER*3      FNCGRD         ! Interpolation function
C
C -----------------------------------------------------------------------------
