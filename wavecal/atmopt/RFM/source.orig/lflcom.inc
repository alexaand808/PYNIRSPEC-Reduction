C LFLCOM.INC -------------------------------------------------------------------
      COMMON / LFLCOM / NLFL, GASLFL, LUNLFL, NDPLFL, NDTLFL, NHDLFL, 
     &                  NSVLFL, SPCLFL, BINLFL, IRRLFL, TABLFL
C
C VERSION
C     06-JUN-00  AD  Add IRRLFL
C     04-JUN-99  AD  Add NHDLFL, BINLFL
C     28-APR-99  AD  Add NDPLFL, NDTLFL
C     20-JAN-98  AD  Add NSVLFL
C     18-JUL-97  AD  Add TABLFL
C     02-JUL-97  AD  Original.
C
C DESCRIPTION
C     Common variables associated with RFM Look-Up Table files.
C     Updated by LUTCHK. NLFL initialised in INPLUT.
C
      INTEGER      NLFL            ! No. of LUT files to be used
      INTEGER      GASLFL(MAXLFL)  ! Index of Gas in each LUT file
      INTEGER      LUNLFL(MAXLFL)  ! LUN for each LUT file
      INTEGER      NDPLFL(MAXLFL)  ! Factor for -lnp axis increments
      INTEGER      NDTLFL(MAXLFL)  ! Factor for Temp axis increments
      INTEGER      NHDLFL(MAXLFL)  ! No.header records to skip in LUT file
      INTEGER      NSVLFL(MAXLFL)  ! No.of Singular Values to use from LUT file
      INTEGER      SPCLFL(MAXLFL)  ! Index of Spectral range in each LUT file
      LOGICAL      BINLFL(MAXLFL)  ! T=Binary file, F=ASCII file
      LOGICAL      IRRLFL(MAXLFL)  ! T=Irregular Grid, F=Full Grid
      CHARACTER*3  TABLFL(MAXLFL)  ! Tabulated Function 
C
C ------------------------------------------------------------------------------
