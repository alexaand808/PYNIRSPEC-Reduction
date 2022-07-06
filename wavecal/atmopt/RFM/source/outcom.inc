C OUTCOM.INC ------------------------------------------------------------------
      COMMON / OUTCOM / WNLOUT, WNUOUT, WNROUT, 
     &          IFLABS, IFLBBT, IFLCOO, IFLOPT, IFLRAD, IFLRJT, IFLTRA,
     &          IOFOUT, IOUT1,  IOUT2, 
     &          LUNABS, LUNBBT, LUNCOO, LUNOPT, LUNRAD, LUNRJT, LUNTAB, 
     &          LUNTRA, LUNWID, 
     &          NPWOUT, NPTOUT, 
     &          LABOUT, TXTHDR, VERSN, 
     &          NAMABS, NAMBBT, NAMCOO, NAMOPT, NAMPRF, NAMPTH, NAMRAD, 
     &          NAMRJT, NAMTAB, NAMTRA, NAMWID, 
     &          RUNID,  FILOUT
C
C VERSION
C     07-AUG-13  AD  Change NAM*, FILOUT, TXTHDR from C*80 to C*200
C     14-SEP-11  AD  Add IFLRTJ, LUNRJT, NAMRJT
C     06-JAN-04  AD  Add NAMPRF
C     09-AUG-03  AD  Add IFLBBT, LUNBBT, NAMBBT
C     24-JAN-03  AD  Add IOFOUT
C     20-JUL-02  AD  Add IFL*
C     07-APR-00  AD  Add NAMCOO, LUNCOO     
C     23-APR-99  AD  Change LABOUT from C*6 to C*8
C     02-FEB-99  AD  Add FILOUT
C     22-OCT-97  AD  Add NAMTAB, LUNTAB
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     List of RFM output file data.
C
      DOUBLE PRECISION WNLOUT    ! Lower wavenumber for output spectra
      DOUBLE PRECISION WNUOUT    ! Upper wavenumber for output spectra
      DOUBLE PRECISION WNROUT    ! Resolution for output spectra
      INTEGER          IFLABS    ! Offset index for output absorption files
      INTEGER          IFLBBT    ! Offset index for output bright.temp files
      INTEGER          IFLCOO    ! Offset index for output cooling rate files
      INTEGER          IFLOPT    ! Offset index for output optical depth files
      INTEGER          IFLRAD    ! Offset index for output radiance files
      INTEGER          IFLRJT    ! Offset index for output Rayl-J.temp files
      INTEGER          IFLTRA    ! Offset index for output transmitt. files
      INTEGER          IOFOUT    ! Offset of output grid in fine grid intvls.
      INTEGER          IOUT1     ! Lower index for output in curr.widemesh
      INTEGER          IOUT2     ! Upper index for output in curr.widemesh
      INTEGER          LUNABS    ! LUN for 1st output absorption file
      INTEGER          LUNBBT    ! LUN for 1st output bright.temp file
      INTEGER          LUNCOO    ! LUN for 1st output coolingrate file
      INTEGER          LUNOPT    ! LUN for 1st output optical depth file
      INTEGER          LUNRAD    ! LUN for 1st output radiance file
      INTEGER          LUNRJT    ! LUN for 1st output Rayl-J.temp file
      INTEGER          LUNTAB    ! LUN for 1st tabulated Abs.Coeff file
      INTEGER          LUNTRA    ! LUN for 1st transmission file
      INTEGER          LUNWID    ! LUN for 1st output Widemesh file
      INTEGER          NPWOUT    ! No.output points per widemesh interval
      INTEGER          NPTOUT    ! Total No.wavenumber points for output
      CHARACTER*8      LABOUT    ! Label for output spectral range
      CHARACTER*200    TXTHDR    ! Text header from driver table
      CHARACTER*11     VERSN     ! RFM version identifier
      CHARACTER*200    NAMABS    ! Name of output ABS file
      CHARACTER*200    NAMBBT    ! Name of output BBT file
      CHARACTER*200    NAMCOO    ! Name of output COO file
      CHARACTER*200    NAMOPT    ! Name of output OPT file
      CHARACTER*200    NAMPRF    ! Name of output PRF file
      CHARACTER*200    NAMPTH    ! Name of output PTH file
      CHARACTER*200    NAMRAD    ! Name of output RAD file
      CHARACTER*200    NAMRJT    ! Name of output RJT file
      CHARACTER*200    NAMTAB    ! Name of output TAB file
      CHARACTER*200    NAMTRA    ! Name of output TRA file
      CHARACTER*200    NAMWID    ! Name of output WID file
      CHARACTER*200    RUNID     ! Optional ID code for current run
      CHARACTER*200    FILOUT(MAXOUT) !Output filenames
C
C -----------------------------------------------------------------------------
