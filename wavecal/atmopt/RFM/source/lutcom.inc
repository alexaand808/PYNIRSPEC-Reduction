C LUTCOM.INC ------------------------------------------------------------------
      COMMON / LUTCOM / NLUT, IGLUT, IGOFF, IRCLUT, IVLUT, IVOFF, IXOFF, 
     &                  LUNLUT, NDPLUT, NDTLUT,
     &                  NGLUT, NLLUT, NPLUT, NTLUT, NVLUT, IDMLUT,
     &                  P1LUT, T1LUT, DPLUT, DTLUT, U, K, 
     &                  BINLUT, IRRLUT, OFFLUT, TABLUT
      COMMON / LUTDBL / V1LUT, DVLUT      
C
C VERSION
C     07-JUN-00  AD  Add IVLUT, IGLUT, IGOFF
C     04-JUN-99  AD  Add BINLUT 
C     28-APR-99  AD  Add NDPLUT, NDTLUT, OFFLUT
C     20-APR-99  AD  Put Double precision variables into separate common block
C     01-NOV-97  AD  Add NLUT, LUNLUT
C     18-JUL-97  AD  Add TABLUT
C     03-JUL-97  AD  Original.
C
C DESCRIPTION
C     Common variables associated with RFM Look-Up Table data.
C     These data are loaded for each spectral range by SPCLUT.
C
      INTEGER     NLUT          ! No.LUTs used for current spectral range
      INTEGER     IGLUT(MAXLUG) ! Wavenumber grid indices for TAB files
      INTEGER     IGOFF(MAXLUT) ! Offset for Wavenumber axis in IGLUT array
      INTEGER    IRCLUT(MAXLUT) ! Last Waveno. Record# read of Abs.Coef files
      INTEGER     IVLUT(MAXLUV) ! Wavenumber grid indices for LUTs 
      INTEGER     IVOFF(MAXLUT) ! Offset for Wavenumber axis in U matrices
      INTEGER     IXOFF(MAXLUT) ! Offset for Wavenumber axis in K matrices
      INTEGER    LUNLUT(MAXLUT) ! LUNs for LUTs using Abs.Coeff. files
      INTEGER    NDPLUT(MAXLUT) ! -lnp axis increment factor (for TAB files)
      INTEGER    NDTLUT(MAXLUT) ! Temp axis increment factor (for TAB files)
      INTEGER     NGLUT(MAXLUT) ! No. wavenumber points for Irreg.Tab files
      INTEGER     NLLUT(MAXLUT) ! No. basis vectors
      INTEGER     NPLUT(MAXLUT) ! No. pressure points 
      INTEGER     NTLUT(MAXLUT) ! No. temperature points 
      INTEGER     NVLUT(MAXLUT) ! No. wavenumber points 
      INTEGER     IDMLUT        ! 4-byte filler for alignment (dummy)
      REAL        P1LUT(MAXLUT) ! Lowest -lnp [p in mb] 
      REAL        T1LUT(MAXLUT) ! Lowest Temperature [K] 
      REAL        DPLUT(MAXLUT) ! -lnp increment
      REAL        DTLUT(MAXLUT) ! Temperature increment [K] 
      REAL     U(MAXLUV,MAXLUL) ! U-matrix (all LUTs combined)
      REAL     K(MAXLUL,MAXLUX) ! K-matrix (all LUTs combined)
      LOGICAL         BINLUT(MAXLUT) ! T=Binary file, F=ASCII file
      LOGICAL         IRRLUT(MAXLUT) ! T=irregular grid
      LOGICAL         OFFLUT(MAXLUT) ! T=assume temperature offset profile
      CHARACTER*3     TABLUT(MAXLUT) ! Tabulation function      
      DOUBLE PRECISION V1LUT(MAXLUT) ! Lowest Wavenumber [/cm] 
      DOUBLE PRECISION DVLUT(MAXLUT) ! Wavenumber increment [/cm]
C
C -----------------------------------------------------------------------------
