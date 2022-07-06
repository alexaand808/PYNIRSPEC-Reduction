      SUBROUTINE SPCGRD ( ISPC, FAIL, ERRMSG )
C
C VERSION
C     17-AUG-06  AD  Allow for -NPT (ie grd in GHz). REWIND LUNGRD 
C     13-JAN-03  AD  Add GRDILS
C     02-JAN-03  AD  Add outcom.inc since *OUT variables now set beforehand
C     27-APR-00  AD  Convert ILS fn to DP
C     11-AUG-99  AD  Change variable names HEADER,NUSE to CDUMMY, IDUMMY.
C                    Make D.P. operations explicit
C     24-JUL-99  AD  Modify CONFIG from C*40 to C*35
C     23-APR-99  AD  Adapt error messages for C*8 LABSPC 
C     22-JAN-98  AD  Replace VAX IBITS function with explicit calculation
C     18-DEC-97  AD  Original.
C
C DESCRIPTION
C     Initialise GRD data for each new spectral range.
C     Called RFMSPC once for each spectral range.
C     Note that all array dimensions have been checked previously by GRDCHK.
C     This routine finds which GRD file (if any) contains the grid points for
C     the current spectral range and loads the grid in GRDCOM.INC.
C     If no grid exists, then NGRD is set to zero.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER       ISPC    !  I  Current spectral range
      LOGICAL       FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  GRDILS ! Calc. Irreg. grid direct contribs to convolved output grid. 
     &, RFMLOG ! Write text message to RFM log file.
C
C GLOBAL CONSTANTS
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C LOCAL CONSTANTS
      DOUBLE PRECISION WNOFAC ! GHz to Wno conversion factor
        PARAMETER ( WNOFAC = 1.0D7 / VLIGHT )   ! VLIGHT in phycon.inc
C
C COMMON VARIABLES
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'ilscom.inc' ! Instrument Lineshape functions.
      INCLUDE 'gflcom.inc' ! Grid file data
      INCLUDE 'grdcom.inc' ! Irregular grid 
      INCLUDE 'outcom.inc' ! RFM output file data
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
C
C LOCAL VARIABLES
      INTEGER I              ! Counter
      INTEGER IBIT           ! Bit counter in Hex character
      INTEGER IGFL           ! GRD file counter
      INTEGER IILS           ! Index of ILS data
      INTEGER IOS            ! Saved value of IOSTAT for error messages
      INTEGER IPT            ! Index of file grid 
      INTEGER IPT1,IPT2      ! Index range of points to be read from file grid
      INTEGER LEFT           ! No.points left to read in GRD file
      INTEGER LUNGRD         ! LUN for GRD data
      INTEGER NPT            ! No.points in grid in GRD file
      INTEGER IDUMMY         ! No.points used in grid in GRD file (dummy)
      INTEGER VEC(50)        ! Record of grid data read from GRD file
      INTEGER ITWO           ! Powers of two
      CHARACTER*25     CONFIG ! GRD file configuration info
      CHARACTER*80     CDUMMY ! Header record read from GRD file (dummy)
      CHARACTER*80     MESSGE ! Message for RFM log file
      DOUBLE PRECISION DW     ! Negligible wavenumber [/cm] for rounding errors
      DOUBLE PRECISION WNO    ! Wavenumber [/cm] of each grid point in GRD file
      DOUBLE PRECISION WNOD   ! Wavenumber increment [/cm] of grid in GRD file
      DOUBLE PRECISION WNOMAX ! Highest Wno [/cm] required from GRD 
      DOUBLE PRECISION WNOMIN ! Lowest Wno [/cm] required from GRD 
      DOUBLE PRECISION WNO1   ! Lowest Wavenumber [/cm] of grid in GRD file
      DOUBLE PRECISION WNO2   ! Highest Wavenumber [/cm] of grid in GRD file
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Set min/max wavenumber points required according to output range for the 
C spectral region
C
      WNOMAX = WNUOUT
      WNOMIN = WNLOUT
      DW     = 0.1D0 * WNRFIN
C
C If using ILS convolution, min/max wavenumber points required are extended by
C width of the instrument line shape either side of output range
C
      IF ( ILSFLG ) THEN
        IILS = ILSSPC(ISPC)
        WNOMIN = WNOMIN + PT1ILS(IILS) 
        WNOMAX = WNOMAX + PT2ILS(IILS)
      ELSE IF ( AVGFLG ) THEN          ! Increment by WNROUT 
        WNOMIN = WNOMIN - WNROUT
        WNOMAX = WNOMAX + WNROUT
      END IF
C
C Loop over all stored GRD files to find which file (if any) contains grid 
C covering current spectral range
C
      OUTGRD = .FALSE.
      DO IGFL = 1, NGFL
        IF ( SPCGFL(IGFL) .EQ. ISPC ) THEN
          LUNGRD = LUNGFL(IGFL)
          REWIND ( LUNGRD, IOSTAT=IOS, ERR=900 )
          READ ( LUNGRD, '(A)', IOSTAT=IOS, ERR=900 ) CONFIG
          WRITE ( MESSGE, '(A,A,A)' )
     &      'I-SPCGRD: Spc.Rng=', LABSPC(ISPC),
     &      ' use GRD Config='//CONFIG
          CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
          READ ( LUNGRD, '(A)', IOSTAT=IOS, ERR=900 ) CDUMMY
          READ ( LUNGRD, '(A)', IOSTAT=IOS, ERR=900 ) FNCGRD
          TWOGRD = INDEX ( 'lin 1li 1sq lor lnl', FNCGRD ) .NE. 0
          READ ( LUNGRD, *, IOSTAT=IOS, ERR=900 )
     &      NPT, IDUMMY, WNO1, WNOD
          IF ( NPT .LT. 0 ) THEN     ! Convert from GHz to cm-1
            NPT = -NPT
            WNO1 = WNO1 * WNOFAC
            WNOD = WNOD * WNOFAC
          END IF
          READ ( LUNGRD, '(A)', IOSTAT=IOS, ERR=900 ) CDUMMY
C
C Although GRD data will cover output range, there is no guarantee it will also
C cover the extra width required for the ILS convolution, so check this and 
C send a warning if grid not wide enough (but continue assuming contributions
C beyond edges of grid are negligible)
C
          IF ( ILSFLG .OR. AVGFLG ) THEN
            WNO2 = WNO1 + ( NPT - 1 ) * WNOD
            IF ( WNO1 .GT. WNOMIN + DW .OR. WNO2 .LT. WNOMAX - DW ) THEN
              WRITE ( MESSGE, '(A,A,A)' )
     &          'W-SPCGRD: ILS wider than GRD, Spc.Rng=', 
     &          LABSPC(ISPC), '. Assume zero absorp.outside grid'
              CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
              IF ( FAIL ) RETURN
            END IF
          END IF
C
C Find lowest calc.pt in terms of index in file grid
          IPT = 1 + NINT ( ( WNOMIN - WNO1 ) / WNOD )
C
C If lowest point above start of file grid then ensure lowest point(s) marked
C for calculation and only read in calc. grid points above that
          IF ( IPT .GT. 1 ) THEN
            WNOGRD(1) = ( IPT - 1 ) * WNOD + WNO1
            NGRD = 1 
            IPT1 = IPT + 1
            IF ( .NOT. TWOGRD ) THEN
              WNOGRD(2) = WNOGRD(1) + WNOD
              NGRD = 2
              IPT1 = IPT1 + 1
            END IF
          ELSE
            NGRD = 0
            IPT1 = 1
          END IF
C
C Find highest calc.pt in terms of index in file grid
          IPT = 1 + NINT ( ( WNOMAX - WNO1 ) / WNOD )
C
C If before end of file grid, set last point to be read below that.
          IF ( IPT .LT. NPT ) THEN
            IPT2 = IPT - 1
            IF ( .NOT. TWOGRD ) IPT2 = IPT2 - 1
          ELSE
            IPT2 = NPT
          END IF
C
C Load selected grid points into /GRDCOM/
          WNO = WNO1 - WNOD
          IPT = 0
          LEFT = 1 + ( NPT - 1 ) / 4
          DO WHILE ( LEFT .GT. 0 )
            READ ( LUNGRD, '(50Z1)', IOSTAT=IOS, ERR=900 ) 
     &        ( VEC(I), I = 1, MIN ( 50, LEFT ) )
            DO I = 1, MIN ( 50, LEFT )
              ITWO = 16
              DO IBIT = 3, 0, -1
                ITWO = ITWO / 2                              ! 8, 4, 2, 1
                IPT = IPT + 1
                WNO = WNO + WNOD
C               IF ( IBITS ( VEC(I), IBIT, 1 ) .EQ. 1 .AND.  ! VAX function
                IF ( MOD ( VEC(I), 2*ITWO ) / ITWO .EQ. 1 .AND. 
     &               IPT .GE. IPT1 .AND. IPT .LE. IPT2 ) THEN
                  NGRD = NGRD + 1
                  WNOGRD(NGRD) = WNO
                END IF
              END DO
            END DO
            LEFT = LEFT - 50
          END DO
C
C If stopped before end of file grid, need to ensure calculated point(s) set at
C the upper limit 
          IF ( IPT2 .LT. NPT ) THEN
            IF ( TWOGRD ) THEN
              WNOGRD(NGRD+1) = IPT2 * WNOD + WNO1
              NGRD = NGRD + 1
            ELSE
              WNOGRD(NGRD+1) = IPT2 * WNOD + WNO1
              WNOGRD(NGRD+2) = WNOGRD(NGRD+1) + WNOD
              NGRD = NGRD + 2
            END IF
          END IF
C
C Attempt direct interpolation from irregular grid to output grid if possible
          IF ( ILSFLG .OR. AVGFLG ) THEN
            CALL GRDILS ( ISPC, FAIL, ERRMSG )
            IF ( FAIL ) RETURN
          END IF            
C
C Normal exit with grid data loaded
          FAIL = .FALSE.
          RETURN
        END IF
      END DO
C
      MESSGE = 'I-SPCGRD: Spc.Rng='//LABSPC(ISPC)//
     &      ': no GRD data applicable so using full grid'
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
C Normal exit with no GRD data available for this spectral range
      RETURN
C
  900 WRITE ( ERRMSG, '(A,I11)')
     &  'F-SPCGRD: I/O failure on GRD file. IOSTAT=', IOS
      FAIL = .TRUE.
C
      END
