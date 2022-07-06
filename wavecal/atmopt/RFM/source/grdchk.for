      SUBROUTINE GRDCHK ( LUNGRD, ISPC, FAIL, ERRMSG )
C
C VERSION
C     07-AUG-13  AD  Change FORMAT descriptors to avoid ifort warnings
C     17-AUG-06  AD  Allow irr.grid to be specified in GHz
C     20-MAR-01  AD  Add test for same wavenumber resln for irreg grids
C     05-JUN-00  AD  Suppress Altitude range check if using TABFLG as well.
C     05-AUG-99  AD  Change names of dummy variables to make use obvious:
C                    HEADER to CDUMMY, WNO1 to DDUMMY, IDUM to IDUMMY 
C     18-DEC-97  AD  Original.
C
C DESCRIPTION
C     Check contents of GRD file.
C     Called by GRDFIL for each file covering a required spectral range.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER       LUNGRD  ! I/O LUN for GRD File/next free LUN
      INTEGER       ISPC    !  I  Spectral Range#
      LOGICAL       FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  RFMLOG ! Write text message to RFM log file.
C
C GLOBAL CONSTANTS
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gflcom.inc' ! Grid file data
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL CONSTANTS
      CHARACTER*39 FNCSTR 
        PARAMETER ( FNCSTR='lin qad cub 1li 1qa 1cu 1sq lor lnl lnc' )
      DOUBLE PRECISION WNOTOL ! Tolerance for wavenumber spacing comparison
        PARAMETER ( WNOTOL = 0.01D0 )
      DOUBLE PRECISION WNOFAC ! GHz to Wno conversion factor
        PARAMETER ( WNOFAC = 1.0D7 / VLIGHT )   ! VLIGHT in phycon.inc
C
C LOCAL VARIABLES
      INTEGER       I       ! Counter
      INTEGER       IDUMMY  ! Dummy integer
      INTEGER       IGFL    ! Counter for GRD files already saved
      INTEGER       IOS     ! Saved value of IOSTAT for error messages
      INTEGER       LEFT    ! No. Hex characters left to read from GRD file
      INTEGER       NPTS    ! Number of original grid points in GRD data
      INTEGER       NUSE    ! No. grid points to be used in GRD data
      REAL          HGTMAX  ! Upper altitude limit for GRD file
      REAL          HGTMIN  ! Lower altitude limit for GRD file
      DOUBLE PRECISION WNOD ! Wavenumber increment of grid in GRD file 
      DOUBLE PRECISION DDUMMY ! Lowest Wno. of grid in GRD file (dummy)
      CHARACTER*40  CONFIG  ! Configuration information
      CHARACTER*9   ERRSTR  ! Record identifier for error messages
      CHARACTER*3   FNC     ! Interpolation method
      CHARACTER*80  CDUMMY  ! Header record read from GRD file (dummy)
      CHARACTER*80  MESSGE  ! Text message sent to LOG file
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Check that a grid hasn't already been found for this spectral range 
C
      DO IGFL = 1, NGFL
        IF ( SPCGFL(IGFL) .EQ. ISPC ) THEN
          FAIL = .TRUE.
          ERRMSG = 'F-GRDCHK: Two GRD files specified '//
     &             'for same spectral range:'//LABSPC(ISPC)
          RETURN
        END IF
      END DO
C
C Check total number of GRD files within array dimension MAXGFL 
C
      IF ( NGFL .EQ. MAXGFL ) THEN
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,I11)' )
     &    'F-GRDCHK: Too many GRDs for one RFM run. '//
     &    'Max = MAXGFL in RFMSIZ.INC=', MAXGFL
        RETURN
      ELSE
        NGFL = NGFL + 1
      END IF
C
C Read configuration info from GRD file and print to RUNLOG file
C
      ERRSTR = 'Config'
      READ ( LUNGRD, '(A)', ERR=900, IOSTAT=IOS ) CONFIG
      MESSGE = 'I-GRDCHK: Configuration info: '//CONFIG
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
C Read file header record - ignored
C
      ERRSTR = 'File Hdr'
      READ ( LUNGRD, '(A)', ERR=900, IOSTAT=IOS ) CDUMMY
C
C Read interpolation function and check against list of recognised functions
C
      ERRSTR = 'Interp.Fn'
      READ ( LUNGRD, '(A)', ERR=900, IOSTAT=IOS ) FNC
      IF ( MOD ( INDEX ( FNCSTR, FNC ), 4 ) .NE. 1 ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-GRDCHK: Unrecognised interpolation method='//FNC
        RETURN
      END IF
C
      ERRSTR = 'Wno.Dims'
      READ ( LUNGRD, *, ERR=900, IOSTAT=IOS ) NPTS, NUSE, DDUMMY, WNOD
      IF ( NPTS .LT. 0 ) THEN           ! Convert GHz to cm-1
        NPTS = -NPTS
        WNOD = WNOD * WNOFAC
      ENDIF
C
C If applying convolution, all irreg grids must be at same resolution 
C (check this here), which is also same as the fine grid 
C (check in FINCHK in case it is modified by *FIN section).
      IF ( ILSFLG .OR. AVGFLG ) THEN
        IF ( NGFL .EQ. 1 ) THEN
          RESGFL = WNOD
        ELSE
          IF ( ABS ( RESGFL - WNOD ) .GT. WNOTOL * RESGFL ) THEN
            FAIL = .TRUE.
            ERRMSG = 'F-GRDCHK: With ILS or AVG convol,'//
     &               'all irreg grids must be on same spacing' 
            RETURN
          END IF
        END IF
      END IF
C
C Check size of irregular grid fits within MAXGRD dimension
C
      IF ( NUSE .GT. MAXGRD ) THEN
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,I11,A,I11)' )
     &    'F-GRDCHK: No.pts in Irreg. grid=', NUSE,
     &    ' > MAXGRD in RFMSIZ.INC =', MAXGRD
        RETURN
      END IF
C
C Check spacing of irregular grid fits within MAXFIN pts/wavenumber
C
      IF ( WNOD .LT. 1.0D0/DBLE(MAXFIN) ) THEN
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,1PE9.2,A,I11,A)' )
     &    'F-GRDCHK: Spacing of irreg.grid=', WNOD,
     &    ' cm-1 < limit, 1/(MAXFIN=', MAXFIN, ')'
        RETURN
      END IF
C
C Read altitude range.
C For limb calculations, check tangent height range is within altitude range 
C of grid (warn if not)
C 
      ERRSTR = 'Alt.Range'
      READ ( LUNGRD, *, ERR=900, IOSTAT=IOS ) HGTMIN, HGTMAX
      IF ( .NOT. NADFLG .AND. .NOT. ZENFLG .AND. 
     &     .NOT. HOMFLG .AND. .NOT. TABFLG        ) THEN
        IF ( HGTTAN(1) .LT. HGTMIN ) THEN
          WRITE ( MESSGE, '(A,F6.1,A)' )
     &      'W-GRDCHK: Using GRD file below nominal minimum tangent '//
     &      'height,=', HGTMIN, ' km'
          CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        END IF
        IF ( HGTTAN(NTAN) .GT. HGTMAX ) THEN
          WRITE ( MESSGE, '(A,F6.1,A)' )
     &      'W-GRDCHK: Using GRD file above nominal maximum tangent '//
     &      'height,=', HGTMAX, ' km'
          CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        END IF
      END IF
C
C Read GRD to check for any errors
C
      ERRSTR = 'Grid Pts'
      LEFT = 1 + (NPTS-1)/4
      DO WHILE ( LEFT .GT. 0 )
        READ ( LUNGRD, '(50Z1)', ERR=900, IOSTAT=IOS ) 
     &    ( IDUMMY, I = 1, MIN ( 50, LEFT ) )
        LEFT = LEFT - 50
      END DO
C
      SPCGFL(NGFL) = ISPC
      LUNGFL(NGFL) = LUNGRD
      FAIL = .FALSE.
      RETURN
C
  900 WRITE ( ERRMSG, '(A,I11)') 
     &  'F-GRDCHK: Error reading '//ERRSTR//' record from GRD file.'//
     &  ' IOSTAT=', IOS
      FAIL = .TRUE.
C
      END
