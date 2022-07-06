      PROGRAM RFM
C
C VERSION
C     16-JAN-14  AD  v4.31 official
C     10-JAN-14  AD  v4.31 - Allow HOM+SFC, Set BrO=MolID#50
C     18-NOV-13  AD  v4.31 - pre-release
C     01-SEP-13  AD  v4.30 official
C     01-APR-08  AD  v4.29 test - Bug#70 (inpspc.for)
C     01-JAN-08  AD  v4.28 official
C     01-JAN-07  AD  v4.27 official
C     01-JAN-06  AD  v4.26 official
C     01-JAN-05  AD  v4.25 official
C     22-SEP-03  AD  v4.23 official
C     10-MAR-03  AD  v4.22 official
C     24-OCT-02  AD  v4.21 official - update Gamache non-lte TP fns qtnte.for
C     27-MAR-02  AD  v4.20 official (no changes)
C     20-SEP-01  AD  v4.11 official.
C     27-DEC-00  AD  v4.10
C     07-NOV-00  AD  v4.09 Pre-release of v4.10
C     14-APR-00  AD  v4.04_GRA bug fixes to WID and SFC
C     08-NOV-99  AD  v4.03 Minor bug fixes
C     05-NOV-99  AD  v4.02 General minor bug-fixes
C     12-OCT-99  AD  v4.01 Correction to handling of SFCFLG with limb viewing
C     01-OCT-99  AD  Version 4.00
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     MIPAS Reference Forward Model (RFM).
C     See http://www.atm.ox.ac.uk/RFM/ for details
C
      IMPLICIT NONE
C
      EXTERNAL
     &  RFMCTM ! Calculate continuum absorption on fine mesh
     &, RFMDEF ! Set RFM default filenames
     &, RFMFIN ! Calculate absorption on fine mesh
     &, RFMFLX ! Perform RFM Flux calculation
     &, RFMFOV ! Apply field-of-view convolution.
     &, RFMGRA ! Construct ray paths through 2D atmosphere.
     &, RFMGRD ! Sets up fine-resolution grid for each widemesh interval
     &, RFMILS ! Apply instrument line shape convolution
     &, RFMINP ! Open RFM Driver file and read input data
     &, RFMINT ! Interpolates irregular to regular fine resolution grid.
     &, RFMJAC ! Calculate Jacobian spectra
     &, RFMLOG ! Write text message to RFM log file
     &, RFMLOS ! Calculate LOS Jacobian spectra.
     &, RFMLUT ! Calculate the absorption using Look-Up Tables.
     &, RFMNAD ! Calculate nadir/zenith path parameters through atmosphere
     &, RFMOPN ! Open RFM output files and write headers
     &, RFMPRF ! Open PRF files and write out internal profile(s)
     &, RFMPTB ! Set up list of pertubations required for RFM Jacobian calc
     &, RFMPTH ! Open PTH files and write RFM ray path diagnostics 
      EXTERNAL
     &  RFMRAD ! Perform RFM radiance calculation
     &, RFMRAY ! Construct ray paths through atmosphere
     &, RFMSCA ! Determine which path absorption calculations can be scaled.
     &, RFMSPC ! Initialise RFM for each new spectral range
     &, RFMTAB ! Calculate the absorption using external Abs.Coeff. Tables.
     &, RFMWID ! Perform wide mesh absorption calculations
     &, RFMWNG ! Interpolate tabulated wide mesh absorption across fine mesh
     &, RFMWRT ! Write RFM output data for current widemesh interval.
     &, RFMXSC ! Calculate the absorption due to cross-section species
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'rfmfil.inc' ! Standard filenames of RFM I/O files
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'logcom.inc' ! Runlog file data
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
      INCLUDE 'widcom.inc' ! Wide mesh data
      INCLUDE 'outcom.inc' ! RFM output file data
C
C LOCAL CONSTANTS
      INTEGER       LUNRUN                  ! LUN for RUNLOG File
        PARAMETER ( LUNRUN = 1 )            
      INTEGER       LUNDRV                  ! LUN for Driver File
        PARAMETER ( LUNDRV = 2 )            
C
C LOCAL VARIABLES
      INTEGER       IOS        ! Saved value of IOSTAT for checks
      INTEGER       ISPC       ! Spectral range counter
      INTEGER       IWID       ! Wide-mesh counter
      INTEGER       LUN        ! LUN counter
      INTEGER       LUNNXT     ! Next available LUN
      INTEGER       LUNSAV     ! Save LUN up to start of spectral range loop
      LOGICAL       FAIL       ! Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG     ! Error message written if FAIL is TRUE
      CHARACTER*80  RFMTXT     ! Message identifying program and version
      CHARACTER*132 MESSGE     ! Other messages for Log file
        DATA  LUNNXT / 7 /
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C VERSN is a CHARACTER*11 field which appears in the header record of all 
C output files. Any local modifications to the RFM should be assigned a 
C different VERSN value, eg '4.00_OXF_13' for Oxford modification#13 to 
C "official" RFM v 4.00 
C
      VERSN = '4.31'           ! Loaded into /OUTCOM/ 
C
C Print program and version identifier to screen
C
      RFMTXT = 'R-RFM:    Program RFM v'//VERSN
      WRITE ( *, * ) RFMTXT 
C
C Prompt for any run identifier label to be appended to all output file
C
C      WRITE ( *, * ) 
C     &  'Optional ID to be appended to filenames (<CR>=none):'
C      READ ( *, '(A)' ) RUNID
C      IF ( RUNID .NE. ' ' ) THEN
C        MESSGE = 'R-RFM:    Filename append string='//RUNID
C        WRITE ( *, * ) MESSGE
C      END IF
C
C Open Runlog file (Transfer LUN to RFMLOG via /LOGCOM/LUNLOG)
C
      LUNLOG = LUNRUN
      NWMLOG = 0                  ! Initialise no.warning messages
      OPEN ( UNIT=LUNLOG, FILE=NAMLOG//RUNID, STATUS='UNKNOWN', 
     &       IOSTAT=IOS )
      IF ( IOS .NE. 0 ) THEN
        WRITE ( ERRMSG, '(A,I11)' ) 
     &  'F-RFM: Failed to open LOG file. IOSTAT=', IOS
        FAIL = .TRUE.
        GOTO 900
      END IF
      CALL RFMLOG ( RFMTXT, FAIL, ERRMSG )  ! Write RFM identifier to log file
      IF ( FAIL ) GOTO 900
C
C Set default filenames
C
      CALL RFMDEF
C
C Load input data
C
      CALL RFMINP ( LUNDRV, LUNNXT, FAIL, ERRMSG )
      IF ( FAIL ) GOTO 900
C
C Write out profile diagnostics if required
C 
      IF ( PRFFLG ) THEN 
        CALL RFMPRF ( LUNNXT, FAIL, ERRMSG )
        IF ( FAIL ) GOTO 900
      END IF
C
C Construct ray paths
C
      IF ( NADFLG .OR. ZENFLG .OR. FLXFLG ) THEN
        CALL RFMNAD 
      ELSE IF ( .NOT. HOMFLG .AND. .NOT. TABFLG ) THEN 
        IF ( GRAFLG ) THEN
          CALL RFMGRA
        ELSE
          CALL RFMRAY 
        END IF
        IF ( .NOT. CLCFLG ) THEN
          CALL RFMSCA ( FAIL, ERRMSG )
          IF ( FAIL ) GOTO 900
        END IF
      END IF
C
C Write out path diagnostics if required
C 
      IF ( PTHFLG ) THEN 
        CALL RFMPTH ( LUNNXT, FAIL, ERRMSG )
        IF ( FAIL ) GOTO 900
      END IF
C
C Calculate perturbed paths if performing Jacobian calculations
C
      IF ( JACFLG ) THEN
        CALL RFMPTB ( FAIL, ERRMSG )
        IF ( FAIL ) GOTO 900
      END IF
C
C Loop over all required spectral ranges
C
      LUNSAV = LUNNXT
      DO ISPC = 1, NSPC
C
C Initialise frequency grids for current range
C
        CALL RFMSPC ( ISPC, FAIL, ERRMSG )
        IF ( FAIL ) GOTO 900
C
C Open output files
C
        CALL RFMOPN ( LUNNXT, ISPC, FAIL, ERRMSG )
        IF ( FAIL ) GOTO 900
        MESSGE = 'R-RFM:    Commencing Wide Mesh Calculation'
        IF ( LABSPC(ISPC) .NE. ' ' ) 
     &    MESSGE = MESSGE(1:42)//' for Spectral Range='//LABSPC(ISPC)
        WRITE ( *, * ) MESSGE(1:80)
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) GOTO 900
C
C Calculate absorption at wide-mesh points
C
        CALL RFMWID ( FAIL, ERRMSG )
        IF ( FAIL ) GOTO 900
C
C Include continuum contribution to absorption
C
        IF ( CTMFLG ) CALL RFMCTM 
C
C Main loop over wide mesh. Fine mesh data is stored for current interval only
C  
        MESSGE = 'R-RFM:    Commencing Fine Mesh Calculation'
        WRITE ( *, * ) MESSGE(1:80)
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) GOTO 900
        DO IWID = 1, NWID 
          WRITE (*,*) IWID, NWID
C
C Set up fine mesh grid for this interval (regular or irregular)
C
          CALL RFMGRD ( IWID )
C         
C Interpolate line wing and continuum absorption to on fine mesh
C       
          CALL RFMWNG ( IWID )
C
C Calculate cross-section absorption on fine mesh
C
          CALL RFMXSC ( FAIL, ERRMSG )
          IF ( FAIL ) GOTO 900
C
C Calculate absorption using Look-up tables on fine mesh
C
          IF ( LUTFLG ) THEN
            CALL RFMTAB ( FAIL, ERRMSG )
            IF ( FAIL ) GOTO 900
            CALL RFMLUT 
          END IF
C
C Calculate line absorption on fine mesh
C
          CALL RFMFIN ( IWID, FAIL, ERRMSG )
          IF ( FAIL ) GOTO 900
C
C Calculate radiance on fine mesh (not required for TAB output)
C
          IF ( FLXFLG ) THEN
            CALL RFMFLX
          ELSE IF ( .NOT. TABFLG ) THEN
            CALL RFMRAD 
          END IF
C
C Calculate Jacobian spectra if required
C
          IF ( JACFLG .OR. MTXFLG ) CALL RFMJAC
C
C Convolve with FOV function if required
C
          IF ( FOVFLG ) CALL RFMFOV
C
C Interpolate to regular grid if required
C
          IF ( GRDFLG .AND. .NOT. TABFLG ) CALL RFMINT 
C
C Convolve with ILS function if required
C
          IF ( ILSFLG .OR. AVGFLG ) 
     &      CALL RFMILS ( IWID, ISPC, FAIL, ERRMSG )
          IF ( FAIL ) GOTO 900
C
C Calculate any LOS Jacobians
          IF ( LOSFLG ) CALL RFMLOS
C
C Output radiances/transmissions for each tangent point
C
          CALL RFMWRT ( IWID, FAIL, ERRMSG )
          IF ( FAIL ) GOTO 900
C
        END DO
C
C If performing ILS convolution, output lags calc. by one widemesh interval
C so output last interval
C
        IF ( ILSFLG .OR. AVGFLG ) THEN        
          CALL RFMILS ( NWID+1, ISPC, FAIL, ERRMSG )
          IF ( FAIL ) GOTO 900
          IF ( LOSFLG ) CALL RFMLOS
          CALL RFMWRT ( NWID+1, FAIL, ERRMSG )
          IF ( FAIL ) GOTO 900
        END IF
C
C Close any files specific to this spectral range
C
        DO LUN = LUNNXT-1, LUNSAV, -1
          CLOSE ( LUN, IOSTAT=IOS )
          IF ( IOS .NE. 0 ) THEN
            FAIL = .TRUE.
            WRITE ( ERRMSG, '(A,I3,A,I11)' ) 
     &        'F-RFM:   Error closing LUN', LUN, '. IOSTAT=', IOS 
            GOTO 900
          END IF
        END DO
        LUNNXT = LUNSAV
C
      END DO
C
C Program termination
C
  900 IF ( FAIL ) THEN
        WRITE (*,*) ERRMSG
        CALL RFMLOG ( ERRMSG, FAIL, ERRMSG )
      ELSE
        IF ( NWMLOG .GT. 0 ) THEN
          WRITE ( MESSGE, '(A,I11,A)' ) 
     &      'R-RFM:', NWMLOG, ' Warning messages in RUNLOG file'
          WRITE ( *, * ) MESSGE(1:80)
        END IF
        MESSGE = 'R-RFM:    Successful completion.'
        WRITE ( *, * ) MESSGE(1:80)
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) GOTO 900
      END IF
C
      END 
