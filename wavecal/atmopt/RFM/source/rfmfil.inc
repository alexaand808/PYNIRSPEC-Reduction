C RFMFIL.INC ------------------------------------------------------------------
C
C VERSION
C     14-SEP-11  AD  Add DEFRJT
C     06-JAN-04  AD  Add DEFPRF
C     09-AUG-03  AD  Add DEFBBT
C     07-APR-00  AD  Add DEFCOO
C     22-OCT-97  AD  Add DEFTAB
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Standard filenames of RFM I/O files.
C
      CHARACTER*7   NAMDRV                  ! Filename of Driver File
        PARAMETER ( NAMDRV = 'rfm.drv' )
      CHARACTER*10  NAMLOG                  ! Filename of Log file
        PARAMETER ( NAMLOG = 'rfm.runlog' )
      CHARACTER*9   DEFABS                  ! Default Filename for ABS output
        PARAMETER ( DEFABS = 'abs_*.asc' )
      CHARACTER*9   DEFBBT                  ! Default Filename for BBT output
        PARAMETER ( DEFBBT = 'bbt_*.asc' )
      CHARACTER*9   DEFCOO                  ! Default Filename for COO output
        PARAMETER ( DEFCOO = 'coo_*.asc' )
      CHARACTER*9   DEFOPT                  ! Default Filename for OPT output
        PARAMETER ( DEFOPT = 'opt_*.asc' )
      CHARACTER*8   DEFPRF                  ! Default Filename for PRF output
        PARAMETER ( DEFPRF = 'prf*.asc' )
      CHARACTER*9   DEFPTH                  ! Default Filename for PTH output
        PARAMETER ( DEFPTH = 'pth_*.asc' )
      CHARACTER*9   DEFRAD                  ! Default Filename for RAD output
        PARAMETER ( DEFRAD = 'rad_*.asc' )
      CHARACTER*9   DEFRJT                  ! Default Filename for RJT output
        PARAMETER ( DEFRJT = 'rjt_*.asc' )
      CHARACTER*9   DEFTAB                  ! Default Filename for TAB output
        PARAMETER ( DEFTAB = 'tab_*.asc' )
      CHARACTER*9   DEFTRA                  ! Default Filename for TRA output
        PARAMETER ( DEFTRA = 'tra_*.asc' )
      CHARACTER*9   DEFWID                  ! Default Filename for WID output
        PARAMETER ( DEFWID = 'wid_*.asc' )
C
C -----------------------------------------------------------------------------
