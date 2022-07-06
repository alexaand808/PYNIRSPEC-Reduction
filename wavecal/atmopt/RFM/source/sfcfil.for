      SUBROUTINE SFCFIL ( LUNSFC, FILSFC, FAIL, ERRMSG )
C
C VERSION
C     23-APR-12  AD  Original.
C
C DESCRIPTION
C     Open, read surface emissivity data file and close.
C     Called by INPSFC if filename listed in *SFC section of Driver Table.
C
      IMPLICIT NONE 
C
C ARGUMENTS
      INTEGER       LUNSFC  !  I  LUN for SFC data file
      CHARACTER*(*) FILSFC  !  I  Name of SFC data file
      LOGICAL       FAIL    !  O  T=A fatal error was detected
      CHARACTER*80  ERRMSG  !  O  Error message written if FAIL is TRUE
C      
      EXTERNAL
     &  OPNFIL ! Open ASCII input file, and skip/Log any initial comments.
C
C GLOBAL CONSTANTS
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'sfccom.inc' ! SFC Data
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
C
C LOCAL VARIABLES      
      INTEGER ISFC         ! Counter for SFC tabulation points
      INTEGER JSFC         ! Min of 2 or MAXSFC
      INTEGER IOS          ! Saved value of IOSTAT for error message
      LOGICAL GHZSFC       ! T=.sfc emissivity tabulated in GHz, F=cm-1
      DOUBLE PRECISION FACTOR ! +ve or -ve depending on wno axis
      CHARACTER*80 WRNMSG  ! Warning message
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      CALL OPNFIL ( LUNSFC, FILSFC, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C      
      READ ( LUNSFC, *, IOSTAT=IOS, ERR=900 ) NSFC
      GHZSFC = NSFC .LT. 0     ! flag for emissivity tabulated in GHz 
      NSFC = ABS ( NSFC )
C
      IF ( NSFC .LT. 1 ) THEN
        WRITE ( ERRMSG, '(A,I11)' )
     &    'F-SFCFIL: No.tabulated SFC pts <1, value=', NSFC
        FAIL = .TRUE.
        RETURN
      ELSE IF ( NSFC .GT. MAXSFC ) THEN
        WRITE ( ERRMSG, '(A,I11,A,I11,A)' )
     &    'F-SFCFIL: No.tabulated SFC pts=', NSFC, 
     &    ' > MAXSFC=', MAXSFC, ' in rfmsiz.inc'
        FAIL = .TRUE.
        RETURN
      END IF
C
      READ ( LUNSFC, *, IOSTAT=IOS, ERR=900 ) 
     &  ( WNOSFC(ISFC), ISFC = 1, NSFC )
      READ ( LUNSFC, *, IOSTAT=IOS, ERR=900 ) 
     &  ( EMSSFC(ISFC), ISFC = 1, NSFC )
C
      CLOSE ( LUNSFC, IOSTAT=IOS, ERR=900 )
C
C If NSFC < 0, flag for values tabulated in GHz axis, so warn if RFM not in
C GHz mode and vice-versa
      FAIL = .FALSE.
      IF ( GHZSFC .AND. .NOT. GHZFLG ) THEN 
        DO ISFC = 1, NSFC
          WNOSFC(ISFC) = WNOSFC(ISFC) * 1.0D7 / VLIGHT   ! approx *1/30
        END DO
        WRNMSG = 'F-SFCFIL: SFC specified in GHz - converting to cm-1'
        CALL RFMLOG ( WRNMSG, FAIL, ERRMSG )
      ELSE IF ( .NOT. GHZSFC .AND. GHZFLG ) THEN ! conversion at output stage
        WRNMSG = 'F-SFCFIL: SFC specified in cm-1 - converting to GHz'
        CALL RFMLOG ( WRNMSG, FAIL, ERRMSG )
      END IF
      IF ( FAIL ) RETURN
C
C Check valid SFC values
      DO ISFC = 1, NSFC
        IF ( EMSSFC(ISFC) .LT. 0.0D0 .OR. EMSSFC(ISFC) .GT. 1.0D0 ) THEN
          ERRMSG = 'F-SFCFIL: Surface Emissivity outside range 0:1'
          FAIL = .TRUE.
          RETURN
        END IF
      END DO
C
C Check that spectral axis either increases or decreases monotonically
      IF ( NSFC .GT. 1 ) THEN 
        JSFC = MIN ( MAXSFC, 2 )  ! allow code to compile even if MAXSFC=1
        FACTOR = WNOSFC(JSFC) - WNOSFC(1)      ! might be zero, so fatal error
        DO ISFC = 2, NSFC       
          IF ( FACTOR * ( WNOSFC(ISFC) - WNOSFC(ISFC-1) ) .LE. 0 ) THEN
            ERRMSG = 'F-SFCFIL: Emissivity spectral grid '//
     &               'is not monotonic'
            FAIL = .TRUE.
            RETURN
          END IF
        END DO
      END IF
C
C Check emissivity values span entire spectral range required for RFM run
C Warn if not
      IF ( WMNSPC .LT. MIN ( WNOSFC(1), WNOSFC(NSFC) ) .OR.
     &     WMXSPC .GT. MAX ( WNOSFC(1), WNOSFC(NSFC) )       ) THEN
        WRNMSG = 'W-SFCFIL: Emissivity values do not cover required'//
     &    ' spectral range.'
        CALL RFMLOG ( WRNMSG, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
C
      FAIL = .FALSE.
      RETURN
C
  900 FAIL = .TRUE.
      WRITE ( ERRMSG, '(A,I11)' ) 
     &  'F-SFCFIL: I/O failure on SFC file. IOSTAT=', IOS
      END
