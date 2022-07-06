C TANCOM.INC ------------------------------------------------------------------
      COMMON / TANCOM / ABSTAN, RADTAN, SZNTAN, TRATAN, 
     &                  ELETAN, GEOTAN, HGTTAN, PSITAN, USRTAN, 
     &                  LOSTAN, LTAN, MTAN, NTAN, IATTAN, IFVTAN, 
     &                  CLCTAN, SFCTAN, USRELE, USRGEO,
     &                  STRTAN
C
C VERSION
C     03-SEP-03  AD  Add SFCTAN
C     08-JUN-01  AD  Add LOSTAN
C     26-MAY-00  AD  Convert ABSTAN and RADTAN from Real to D.P.
C     28-DEC-99  AD  Add PSITAN
C     02-JAN-99  AD  Add LTAN
C     23-DEC-98  AD  Comment changes only
C     17-JUL-98  AD  Add USRTAN, ELETAN, GEOTAN, SZNTAN, STRTAN
C     10-DEC-97  AD  Comment changes only
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     List of tangent heights required for RFM calculations.
C     Loaded by INPTAN, FOVTAN.
C     SFCTAN set in CHKLIM
C     Output tangent heights are 1:NTAN (1=lowest, NTAN=highest) limb view
C     but not sorted for NAD, ZEN or HOM views.
C     Extra tangent heights for FOV calculation are in NTAN+1:MTAN
C     After FOV calculation (or without FOV calculation) extra Jacobians are
C     arranged into NTAN+1:LTAN (in RFMPTB or RFMFOV)
C
      DOUBLE PRECISION ABSTAN(MAXFIN,MAXTAN) ! Tangent path optical depth
      DOUBLE PRECISION RADTAN(MAXFIN,MAXTAN) ! Tangent path Radiance
      DOUBLE PRECISION SZNTAN(MAXTAN)        ! Sin(Zenith angle) at HGTTAN
      DOUBLE PRECISION TRATAN(MAXFIN,MAXTAN) ! Tangent path Transmission
      REAL             ELETAN(MAXTAN)        ! Elevation angles of paths
      REAL             HGTTAN(MAXTAN)        ! Lowest altitude along paths [km]
      REAL             GEOTAN(MAXTAN)        ! Projected (=Geom) tan.alt. [km]
      REAL             PSITAN(MAXTAN)        ! LOS angle of tan.pt. [deg]
      REAL             USRTAN(MAXTAN)        ! Tan.path values from Drv.Table
      INTEGER          LOSTAN                ! No. t.pths for LOS jacobians
      INTEGER          LTAN                  ! No. t.pts for o/p incl.Jacobians
      INTEGER          MTAN                  ! No. tangent pts for calculations
      INTEGER          NTAN                  ! No.nominal tan pts for output
      INTEGER          IATTAN(MAXTAN)        ! Atmos layer containing tan pt
      INTEGER          IFVTAN(MAXFOV,MAXTAN) ! Index of Tan.Hgts for FOV Conv.
      LOGICAL          CLCTAN(MAXTAN)        ! T=Radiance.calc reqd for Tan.Hgt
      LOGICAL          SFCTAN(MAXTAN)        ! T=limb ray intersects surface
      LOGICAL          USRELE                ! T=USRTAN are elevation angles
      LOGICAL          USRGEO                ! T=USRTAN are Geom.Tan.points.
      CHARACTER*6      STRTAN(MAXTAN)        ! Strings for filenames
C
C -----------------------------------------------------------------------------
