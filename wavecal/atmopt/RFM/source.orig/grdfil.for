      SUBROUTINE GRDFIL ( LUNGRD, NAMGRD, FAIL, ERRMSG )
C
C VERSION
C     17-AUG-06  AD  Interpret -ve NPTS as grid info in GHz rather than cm-1
C     20-MAR-01  AD  Add check on grid resolution for each spectral range
C     05-DEC-00  AD  Comment out READONLY
C     05-JUL-99  AD  Rename HEADER to CDUMMY, NUSE to IDUMMY.
C     18-DEC-97  AD  Original.
C
C DESCRIPTION
C     Open GRD file and check contents
C     Called by INPGRD for each file specified in *GRD section of driver table.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER       LUNGRD  ! I/O LUN for GRD File/next free LUN
      CHARACTER*(*) NAMGRD  !  I  Name of GRD file
      LOGICAL       FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  GRDCHK ! Check Look-Up Table data.
     &, RFMLOG ! Write text message to RFM log file.
C
C GLOBAL CONSTANTS
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C LOCAL CONSTANTS
      DOUBLE PRECISION WNOTOL ! Wavenumber tolerance as fraction of o/p resln
        PARAMETER ( WNOTOL = 0.01D0 )
      DOUBLE PRECISION WNOFAC ! GHz to Wno conversion factor
        PARAMETER ( WNOFAC = 1.0D7 / VLIGHT )   ! VLIGHT in phycon.inc
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
C
C LOCAL VARIABLES
      INTEGER       IOS     ! Saved value of IOSTAT for error messages
      INTEGER       ISPC    ! Spectral range counter (1:NSPC)
      INTEGER       NPTS    ! No.of wavenumber grid points in GRD file 
      INTEGER       IDUMMY  ! No.of grid points used (dummy)
      LOGICAL       REJECT  ! T=ignore file - no useful data
      DOUBLE PRECISION WNOD ! Wavenumber increment [cm-1] of grid 
      DOUBLE PRECISION WNO1 ! Wavenumber [cm-1] of first grid pt (lowest Wno).
      DOUBLE PRECISION WNO2 ! Wavenumber [cm-1] of last grid pt (highest Wno).
      CHARACTER*80  MESSGE  ! Text message sent to LOG file
      CHARACTER*80  CDUMMY  ! Header records read from GRD file (dummy)
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Send message to LOG file saying which file is about to be opened
C
      MESSGE = 'I-GRDFIL: Opening GRD File: '//NAMGRD
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
C Open File 
C
      OPEN ( UNIT=LUNGRD, FILE=NAMGRD, STATUS='OLD', 
     &        IOSTAT=IOS, ERR=900 )
      REJECT = .TRUE.
C
C First just extract wavenumber range to see if this file is useful.
C (If useful, other details will be checked by GRDCHK).
C
      READ ( LUNGRD, '(A)', ERR=900, IOSTAT=IOS ) CDUMMY 
      READ ( LUNGRD, '(A)', ERR=900, IOSTAT=IOS ) CDUMMY
      READ ( LUNGRD, '(A)', ERR=900, IOSTAT=IOS ) CDUMMY
      READ ( LUNGRD, *, ERR=900, IOSTAT=IOS ) NPTS, IDUMMY, WNO1, WNOD
C
C Check if grid is in cm-1 or GHz and warn if different from GHZ Flag status
      IF ( NPTS .LT. 0 ) THEN
        NPTS = -NPTS
        WNO1 = WNO1 * WNOFAC
        WNOD = WNOD * WNOFAC
        IF ( .NOT. GHZFLG ) THEN 
          MESSGE = 'W-GRDFIL: GRD specified in GHz - converting to cm-1'
          CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        END IF
      ELSE
        IF ( GHZFLG ) THEN
          MESSGE = 'W-GRDFIL: GRD specified in cm-1, not GHz'
          CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        ENDIF
      ENDIF          
C
      WNO2 = WNO1 + (NPTS-1) * WNOD
C
C Check if grid completely spans any selected output spectral range 
C NB: same grid may be usable for more than one spectral range
C
      DO ISPC = 1, NSPC
        IF ( WNO1 .LE. WNLSPC(ISPC) .AND. 
     &       WNO2 .GE. WNUSPC(ISPC) .AND.
C if not using any spectral convolution, grid spacing must also match output
C spacing for this spectral range
     &       ( ILSFLG .OR. AVGFLG .OR.           ! either using convolution, or
     &         ABS ( WNOD - WNRSPC(ISPC) ) .LE.  ! check resln
     &         WNOTOL * WNRSPC(ISPC) ) ) THEN    ! OK, so want to use this GRD
          REWIND ( LUNGRD, ERR=900, IOSTAT=IOS )
          CALL GRDCHK ( LUNGRD, ISPC, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
          MESSGE = 'I-GRDFIL: Grid assigned to spectral range:'
     &               //LABSPC(ISPC)
          CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
          REWIND ( LUNGRD, ERR=900, IOSTAT=IOS )
          REJECT = .FALSE.
        END IF
      END DO
C
      IF ( REJECT ) THEN
        CLOSE ( LUNGRD, IOSTAT=IOS, ERR=900 )
        MESSGE = 'W-GRDFIL: File ignored - grid not applicable'
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      ELSE
        LUNGRD = LUNGRD + 1
      END IF
C
      FAIL = .FALSE.
      RETURN
C
  900 WRITE ( ERRMSG, '(A,I11)') 
     &  'F-GRDFIL: I/O failure on GRD file. IOSTAT=', IOS
      FAIL = .TRUE.
C
      END
