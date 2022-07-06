      SUBROUTINE OPNOUT ( LUNNXT, ISPC, TITLE, NAMTMP, NFLOUT, 
     &                    FAIL, ERRMSG )
C
C VERSION
C     20-SEP-12  AD  Change NAMTMP, NAMFIL, NAMLOW from C*80 to C*200
C     04-APR-06  AD  Change 'X' to '1X' in output format
C     06-JAN-04  AD  Add IGRA argument to MAKNAM
C     01-JAN-04  AD  Add LEV flags
C     09-AUG-03  AD  Allow for GHZ option
C     15-AUG-02  AD  Correct previous modifications to preserve LUN ordering
C     20-JUL-02  AD  Add NFLOUT argument, revise handling of LUNs
C     14-DEC-01  AD  Ensure LUN,LUNOUT initialised to avoid warnings
C     08-JUN-01  AD  Allow for LOS Jacobian spectra, add LOSNAM
C                    remove CARRIAGECONTROL='LIST'
C     21-JAN-01  AD  Remove redundant initial space from header records.
C     06-DEC-00  AD  Comment out CARRIAGECONTROL='LIST' 
C     22-APR-00  AD  Add FLX and MTX options
C     18-JAN-00  AD  Documentation change only
C     03-FEB-99  AD  Allow for LUNFLG
C     03-JAN-99  AD  Allow for Jacobian spectra
C     29-OCT-98  AD  Change header formt from F18.12,F19.16,F18.12 to 3F18.10
C     18-JUL-98  AD  Allow for Elev.angle and Geom.Tan.Hgt path specs.
C     10-MAR-97  AD  Add CARRIAGECONTROL='LIST'
C     03-MAR-97  AD  Version 3.
C     15-JAN-97  AD  Add Zenith path option. 
C                    Change tangent height to airmass for ZEN, NAD paths
C     01-OCT-96  AD  Version 2.
C     18-SEP-96  AD  Avoid converting filename to lower case
C                    Search for ".asc" rather than "asc" (also .bin)
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Open spectral output files for current spectral range
C     Called by RFMOPN for each spectral range
C     Sequence of LUN assignments is:
C        TAN#1, TAN#2, ... TAN#NTAN, 
C        dTAN#1/dJAC#1, dTAN#2/dJAC#1, ... dTAN#NTAN/dJAC#1
C        ...
C        dTAN#1/dJAC#NJAC, dTAN#2/dJAC#NJAC, ... dTAN#NTAN/dJAC#NJAC
C        dTAN#1/dLOS, dTAN#2/dLOS, ... dTAN#NTAN/dLOS
C     This matches output sequence in RFMWRT. 
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER       LUNNXT  ! I/O Next available LUN
      INTEGER       ISPC    !  I  Spectral range index 
      CHARACTER*(*) TITLE   !  I  Code for output type
      CHARACTER*200 NAMTMP  ! I/O Filename template (converted to lower case)
      INTEGER       NFLOUT  ! I/O No.of output files generated (LUNFLG option)
      LOGICAL       FAIL    !  O  Set TRUE if a fatal error occurs
      CHARACTER*80  ERRMSG  !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  JACNAM ! Construct filename for RFM Jacobian files.
     &, LOCASE ! Convert text string to lower case.
     &, LOSNAM ! Construct filename for RFM LOS-Jacobian files.
     &, MAKNAM ! Construct filename for RFM output files.
     &, RFMLOG ! Write text message to RFM log file.
C
C GLOBAL CONSTANTS
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'rfmcon.inc' ! Constants for RFM calculations
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'jaccom.inc' ! Jacobian data
      INCLUDE 'outcom.inc' ! RFM output file data
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL VARIABLES
      INTEGER       IFLOUT ! Output file counter
      INTEGER       IJAC   ! Counter for Jacobian elements
      INTEGER       IOS    ! Value of IOSTAT saved for error messages.
      INTEGER       IPT    ! Pointer to start of substring
      INTEGER       ITAN   ! Tangent height counter
      INTEGER       LUN    ! Actual value of Logical Unit Number.
      INTEGER       NFLLEV ! No.files for intermediate level calcs (LEVFLG)
      LOGICAL       LOSOPN ! T=open file for LOS Jacobian
      DOUBLE PRECISION GHZFAC ! Wno to GHz conversion factor
      CHARACTER*80  HEADR1 ! Header record for output files
      CHARACTER*80  HEADR2 ! 
      CHARACTER*80  HEADR3 ! titles, or Jacobian information (JACFLG)
      CHARACTER*15  LABEL  ! Label for output data 
      CHARACTER*132 MESSGE ! Text message sent to Log file.
      CHARACTER*200 NAMFIL ! Name of file actually opened (incl. RUNID)
      CHARACTER*200 NAMLOW ! Lower case version of filename 
      CHARACTER*32  PTHTXT ! Text describing viewing geometry
      CHARACTER*13  TYPTXT ! Text describing type of spectrum
      CHARACTER*10  WNOGHZ ! 'Wavenumber' or 'Frequency ' in HEADR3
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      CALL LOCASE ( NAMTMP, NAMLOW )
C
C Convert '.asc' to '.bin' in filename, or vice-versa, if appropriate 
C
      IF ( BINFLG ) THEN
        IPT = INDEX ( NAMLOW, '.asc' )
        IF ( IPT .NE. 0 ) NAMTMP(IPT:IPT+3) = '.bin'
      ELSE
        IPT = INDEX ( NAMLOW, '.bin' )
        IF ( IPT .NE. 0 ) NAMTMP(IPT:IPT+3) = '.asc'
      END IF
C
C Set according to whether output is in wavenumber or GHz
      IF ( GHZFLG ) THEN
        WNOGHZ = 'Frequency '
        GHZFAC = VLIGHT * 1.0D-7         ! VLIGHT in phycon.inc
      ELSE
        WNOGHZ = 'Wavenumber'
        GHZFAC = 1.0D0
      END IF
C
C Set header records for output files ('123456789' replaced by tangent height)
C
      TYPTXT = TITLE
      IF ( HOMFLG ) THEN
        PTHTXT = 'Homog. Path Length =123456789 km'
      ELSE IF ( FLXFLG ) THEN
        PTHTXT = 'for Altitude Level =123456789 km' 
      ELSE IF ( NADFLG ) THEN
        PTHTXT = 'Nadir Path Airmass =123456789   '
      ELSE IF ( ZENFLG ) THEN
        PTHTXT = 'Zenith Path Airmass=123456789   '
      ELSE IF ( USRELE ) THEN
        PTHTXT = 'View Elevation Ang.=123456789 dg'
      ELSE IF ( USRGEO ) THEN
        PTHTXT = 'Limb Geom.Tang.Hgt =123456789 km'
      ELSE                                    ! Normal limb-viewing geometry
        PTHTXT = 'Limb Path Tang.Hgt =123456789 km'
      END IF
      IF ( FLXFLG ) THEN
        HEADR1 = '! Flux '//TYPTXT//' calc.'//PTHTXT//
     &           ' by RFM v.'//VERSN
      ELSE
        HEADR1 = '! '//TYPTXT//' calc. for '//PTHTXT//
     &           ' by RFM v.'//VERSN
      END IF
      IPT = INDEX ( HEADR1, '123456789' )
      HEADR2 = '!'//TXTHDR
      HEADR3 = '!No.Pts  Lower_'//WNOGHZ//'   Delta_'//WNOGHZ//'  '//
     &         'Upper_'//WNOGHZ//'  Label'
C
      IF ( LABOUT .EQ. ' ' ) THEN   ! No spectral range label,so use type label
        LABEL = ''''//TITLE//''''
      ELSE 
        LABEL = ''''//LABOUT//''''
      END IF
C
C If using LUN flag, fix LUN for output files and
C check that there is enough space to store output files
      IF ( LUNFLG ) THEN
        LUN = LUNNXT
        IFLOUT = NTAN            ! Number of additional files required
        IF ( JACFLG .OR. MTXFLG .OR. LEVFLG ) 
     &       IFLOUT = IFLOUT + NTAN * NJAC
        IF ( LOSFLG ) IFLOUT = IFLOUT + NTAN
        IF ( NFLOUT + IFLOUT .GT. MAXOUT ) THEN         
          FAIL = .TRUE.
          WRITE ( ERRMSG, '(A,I11)' ) 'F-OPNOUT: No.output files '//
     &      '> MAXOUT in rfmsiz.inc, =', MAXOUT
          RETURN
        END IF
      END IF
      NFLLEV = 0
C
C Note that HEADR3 is modified by JACNAM if JAC flag enabled 
      DO ITAN = 1, NTAN       
        HEADR3 = '!No.Pts  Lower_'//WNOGHZ//'   Delta_'//WNOGHZ//'  '//
     &           'Upper_'//WNOGHZ//'  Label'
        IJAC = 0
        LOSOPN = .FALSE.
        CALL MAKNAM ( ISPC, 0, ITAN, 0, NAMTMP, NAMFIL )
C
C Establish file/LUN counter for non Jacobian spectra
        IF ( LUNFLG ) THEN
          IFLOUT = NFLOUT + ITAN 
        ELSE
          LUN = LUNNXT + ITAN - 1
        END IF
C
  100   CONTINUE      ! Repeat from here for each Jacobian spectrum
        IF ( LUNFLG ) FILOUT(IFLOUT) = NAMFIL
        MESSGE = 'I-OPNOUT: Opening output file: '//NAMFIL
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
C
C Open file with appropriate STATUS, FORM according to NEWFLG, BINFLG
C
        IF ( NEWFLG ) THEN
          IF ( BINFLG ) THEN
            OPEN ( UNIT=LUN, FILE=NAMFIL, STATUS='NEW',
     &             FORM='UNFORMATTED', IOSTAT=IOS, ERR=900 )
          ELSE 
            OPEN ( UNIT=LUN, FILE=NAMFIL, STATUS='NEW',
     &             FORM='FORMATTED', IOSTAT=IOS, ERR=900 )
          END IF
        ELSE
          IF ( BINFLG ) THEN
            OPEN ( UNIT=LUN, FILE=NAMFIL, STATUS='UNKNOWN',
     &             FORM='UNFORMATTED', IOSTAT=IOS, ERR=900 )
          ELSE
            OPEN ( UNIT=LUN, FILE=NAMFIL, STATUS='UNKNOWN',
     &             FORM='FORMATTED', IOSTAT=IOS, ERR=900 )
          END IF
        END IF
C
        WRITE ( HEADR1(IPT:IPT+8), '(F9.3)' ) USRTAN(ITAN) 
        IF ( BINFLG ) THEN
          WRITE ( LUN, IOSTAT=IOS, ERR=900 ) HEADR1
          WRITE ( LUN, IOSTAT=IOS, ERR=900 ) HEADR2
          WRITE ( LUN, IOSTAT=IOS, ERR=900 ) HEADR3
          WRITE ( LUN, IOSTAT=IOS, ERR=900 )
     &      NPTOUT, WNLOUT*GHZFAC, WNROUT*GHZFAC, WNUOUT*GHZFAC, LABEL
        ELSE
          WRITE ( LUN, '(A)', IOSTAT=IOS, ERR=900 ) HEADR1
          WRITE ( LUN, '(A)', IOSTAT=IOS, ERR=900 ) HEADR2
          WRITE ( LUN, '(A)', IOSTAT=IOS, ERR=900 ) HEADR3
          WRITE ( LUN, '(I7, 3F18.10, 1X, A)', IOSTAT=IOS, ERR=900)
     &      NPTOUT, WNLOUT*GHZFAC, WNROUT*GHZFAC, WNUOUT*GHZFAC, LABEL
        END IF
C
        IF ( LUNFLG ) CLOSE ( LUN, IOSTAT=IOS, ERR=900 )
C
C If calculating Jacobians, open file for each Jacobian element
        IF ( ( JACFLG .OR. MTXFLG ) .AND. IJAC .LT. NJAC ) THEN
          IJAC = IJAC + 1
          CALL JACNAM ( ISPC, ITAN, IJAC, NAMTMP, NAMFIL, HEADR3 )
          IF ( LUNFLG ) THEN
            IFLOUT = IFLOUT + NTAN
          ELSE
            LUN = LUN + NTAN
          END IF
          GOTO 100
        END IF
C
C If calculating intermediate output levels, variable number of outputs for
C each nominal tangent path (unlike Jacobians)
        IF ( LEVFLG ) THEN
  200     CONTINUE
          IJAC = IJAC + 1
          IF ( ITNJAC(ITAN,IJAC) .NE. 0 ) THEN
            CALL JACNAM ( ISPC, ITAN, IJAC, NAMTMP, NAMFIL, HEADR3 )
            NFLLEV = NFLLEV + 1
            IF ( LUNFLG ) THEN
              IFLOUT = NFLOUT + NTAN + NFLLEV
            ELSE
              LUN = LUNNXT + NTAN - 1 + NFLLEV
            END IF
            GOTO 100
          ELSE 
            IF ( IJAC .LT. NJAC ) GOTO 200
          END IF
        END IF
C
C If calculating pointing Jacobians, open file for each tangent altitude
        IF ( LOSFLG .AND. .NOT. LOSOPN ) THEN
          CALL LOSNAM ( ISPC, ITAN, NAMTMP, NAMFIL, HEADR3 )
          LOSOPN = .TRUE.
          IF ( LUNFLG ) THEN
            IFLOUT = IFLOUT + NTAN
          ELSE
            LUN = LUN + NTAN
          END IF
          GOTO 100
        END IF 
C
      END DO
C
C Update next available LUN or number of files opened
      IF ( LUNFLG ) THEN
        NFLOUT = IFLOUT
      ELSE
        LUNNXT = LUN + 1
      END IF
C
      FAIL = .FALSE.
      RETURN
C
  900 FAIL = .TRUE.
      WRITE ( ERRMSG, '(A,I11)' ) 
     &  'F-OPNOUT: I/O failure on output file. IOSTAT=', IOS
      END
