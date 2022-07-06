      SUBROUTINE OPNTAB ( LUNNXT, ISPC, FAIL, ERRMSG )
C
C VERSION
C     07-AUG-13  AD  Remove redundant local variable IISO
C     23-MAY-06  AD  Put format statement on to single line
C     06-APR-06  AD  Change 'X' to '1X' in format
C     10-JUN-05  AD  Change output format statement to allow more space for NV 
C                    and ensure space between each field
C     25-MAR-04  AD  Extra argument to MAKNAM. Allow for isotopes.
C     17-DEC-01  AD  Correction: ensure irreg.grd written to all tab files
C     05-DEC-00  AD  Comment out CARRIAGECONTROL='LIST'
C     07-JUN-00  AD  Allow for irregular grid output as well.
C     04-JUL-99  AD  Use I4 format for NP and NT (was I3, I5)
C     04-JUN-99  AD  Write WNL and WNORES as DOUBLE PRECISION
C                    Use C*80 for output header records
C     23-APR-99  AD  Adapt for C*8 label instead of C*6 
C     05-NOV-97  AD  Add extra header record with column headers (NL NV etc)
C                    Correct use of NEWFLG so header also written
C     04-NOV-97  AD  Check WNRSPC for no.pts/wno instead of wno resolution.
C                    Also correct spectral range for ISPC .GT. 1
C     31-OCT-97  AD  Original.
C
C DESCRIPTION
C     Open TAB files for current spectral range
C     Called by RFMOPN for each spectral range.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      LUNNXT  ! I/O Next available LUN
      INTEGER      ISPC    !  I  Spectral range index
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error occurs
      CHARACTER*80 ERRMSG  !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  LOCASE ! Convert text string to lower case.
     &, MAKNAM ! Construct filename for RFM output files.
     &, RFMLOG ! Write text message to RFM log file.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'rfmcon.inc' ! Constants for RFM calculations
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'grdcom.inc' ! Irregular grid
      INCLUDE 'outcom.inc' ! RFM output file data
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
      INCLUDE 'tabcom.inc' ! Table axes for tabulated absorb.coeffs.
C
C LOCAL CONSTANTS
      DOUBLE PRECISION WNOTOL     ! Tolerance as fraction of fine grid spacing 
        PARAMETER ( WNOTOL = 0.01D0 )
C
C LOCAL VARIABLES
      INTEGER       ICHAR  ! Character position in irregular grid code
      INTEGER       IDX    ! HITRAN index of molecule
      INTEGER       IFACT  ! Factor to negate NWNO to flag irreg.grid data
      INTEGER       IGAS   ! Absorber counter
      INTEGER       IGRD   ! Counter for irreg.grid output points
      INTEGER       IOS    ! Value of IOSTAT saved for error messages.
      INTEGER       IPT    ! Pointer to start of substring
      INTEGER       ISUM   ! Sum for irregular grid construction
      INTEGER       ITWO   ! Number establishing bit in irreg.grid code
      INTEGER       IWNO   ! Wavenumber counter on fine grid spacing
      INTEGER       LUN    ! Local value of Logical Unit Number.
      INTEGER       NWNO   ! No.of wavenumber pts to be output (-ve if irreg)
      CHARACTER*9   CODTAB ! Gas code + optional isotopic ID (eg 'i1')
      CHARACTER*80  HEADR1 ! Header record for output files
      CHARACTER*80  HEADR2 !
      CHARACTER*80  HEADR3 ! 
      CHARACTER*132 MESSGE ! Text message sent to Log file.
      CHARACTER*80  NAMFIL ! Name of file actually opened (incl. RUNID)
      CHARACTER*80  NAMLOW ! Lower-case version of NAMTAB
      CHARACTER*17  OUTSTR ! ID Label for output file
      CHARACTER*50  REC50  ! Record for coded irregular grid
      DOUBLE PRECISION WNO    ! Wavenumber [cm-1]
      DOUBLE PRECISION WNORES ! Wavenumber resolution [cm-1]
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      CALL LOCASE ( NAMTAB, NAMLOW )
C
C Convert '.asc' to '.bin' in filename, or vice-versa, if appropriate
C
      IF ( BINFLG ) THEN
        IPT = INDEX ( NAMLOW, '.asc' )
        IF ( IPT .NE. 0 ) NAMTAB(IPT:IPT+3) = '.bin'
      ELSE
        IPT = INDEX ( NAMLOW, '.bin' )
        IF ( IPT .NE. 0 ) NAMTAB(IPT:IPT+3) = '.asc'
      END IF
C
      IF ( WNRSPC(ISPC) .LT. 1 ) THEN
        WNORES = WNRSPC(ISPC)
      ELSE
        WNORES = 1.D0/WNRSPC(ISPC)
      END IF
      NWNO = 1 + NINT ( ( WNUSPC(ISPC) - WNLSPC(ISPC) ) / WNORES )
      IF ( GRDFLG ) THEN
        IFACT = -1          ! Use NWNO = -NWNO to flag irregular grid data
      ELSE
        IFACT = 1
      END IF             
C
      LUNTAB = LUNNXT 
      LUNNXT = LUNNXT + NGAS         ! NGAS files will be created
      DO IGAS = 1, NGAS
        LUN = LUNTAB + ( IGAS - 1 ) 
        CALL MAKNAM ( ISPC, IGAS, 0, 0, NAMTAB, NAMFIL )
        MESSGE = 'I-OPNTAB: Opening output file: '//NAMFIL
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( NEWFLG ) THEN
          IF ( BINFLG ) THEN
            OPEN ( UNIT=LUN, FILE=NAMFIL, STATUS='NEW',
     &             FORM='UNFORMATTED', IOSTAT=IOS, ERR=900 )
          ELSE
            OPEN ( UNIT=LUN, FILE=NAMFIL, STATUS='NEW',
     &           FORM='FORMATTED', IOSTAT=IOS, ERR=900 )
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
        CODTAB = ' '
        IPT = INDEX ( CODGAS(IGAS)//' ', ' ' ) - 1
        CODTAB = CODGAS(IGAS)(1:IPT)
        IF ( ISOMOL(IDXGAS(IGAS)) ) 
     &    WRITE ( CODTAB(IPT+1:IPT+2), '(A,I1)' ) 'i', IDIGAS(IGAS)
        HEADR1 = '! '//CODTAB//
     &           ' Tabulated Absorp.Coeff. '
     &           //'created by RFM v.'// VERSN//'                '
        HEADR2 = '!'//TXTHDR(1:79)
        HEADR3 = '!NL     NV    V1        DV      NP  P1'//
     &      '        DP        NT     T1        DT'
C
C Find out if this gas has been split into isotopes
        IDX = IDXGAS(IGAS)
        IF ( ISOMOL(IDX) ) THEN
          WRITE ( OUTSTR, '(A8,1X,I2,A,I1,1X,A3)' ) 
     &      LABOUT, IDX, '.', IDIGAS(IGAS), 'LIN'
        ELSE
          WRITE ( OUTSTR, '(A8,1X,I2,1X,A3)' ) LABOUT, IDX, 'LIN'
        END IF
        IF ( BINFLG ) THEN
          WRITE ( LUN, IOSTAT=IOS, ERR=900 ) HEADR1
          WRITE ( LUN, IOSTAT=IOS, ERR=900 ) HEADR2
          WRITE ( LUN, IOSTAT=IOS, ERR=900 ) HEADR3
          WRITE ( LUN, IOSTAT=IOS, ERR=900 ) OUTSTR
          WRITE ( LUN, IOSTAT=IOS, ERR=900 ) 
     &      0, IFACT*NWNO, WNLSPC(ISPC), WNORES,       !0=flag not LUT file
     &      NLPTAB, LP1TAB, LPDTAB, NTMTAB, TM1TAB, TMDTAB
        ELSE
          WRITE ( LUN, '(A)', IOSTAT=IOS, ERR=900 ) HEADR1
          WRITE ( LUN, '(A)', IOSTAT=IOS, ERR=900 ) HEADR2
          WRITE ( LUN, '(A)', IOSTAT=IOS, ERR=900 ) HEADR3
          WRITE ( LUN, '(A)', IOSTAT=IOS, ERR=900 ) OUTSTR
          WRITE ( LUN, 
     &    '(I2,1X,I7,2(1X,F9.4),1X,I3,2(1X,F9.6),1X,I3,2(1X,F9.3))', 
     &            IOSTAT=IOS, ERR=900 )
     &      0, IFACT*NWNO, WNLSPC(ISPC), WNORES,       !0=flag not LUT file
     &      NLPTAB, LP1TAB, LPDTAB, NTMTAB, TM1TAB, TMDTAB
        END IF
C
        IF ( GRDFLG ) THEN
          IGRD = 1
          ISUM = 0
          ITWO = 16
          ICHAR = 0
          WNO = WNLSPC(ISPC)
          DO IWNO = 1, NWNO
            ITWO = ITWO/2
            IF ( ABS ( WNOGRD(IGRD) - WNO ) .LT. WNOTOL*WNORES ) THEN
              ISUM = ISUM + ITWO
              IGRD = IGRD + 1
            END IF
            IF ( ITWO .EQ. 1 .OR. IWNO .EQ. NWNO ) THEN
              ICHAR = ICHAR + 1
              WRITE ( REC50(ICHAR:ICHAR), '(Z1)' ) ISUM
              ISUM = 0
              ITWO = 16
              IF ( ICHAR .EQ. 50 .OR. IWNO .EQ. NWNO ) THEN
                IF ( BINFLG ) THEN
                  WRITE ( LUN, IOSTAT=IOS, ERR=900 ) REC50
                ELSE
                  WRITE ( LUN, '(A50)', IOSTAT=IOS, ERR=900 ) REC50
                END IF
                REC50 = ' '
                ICHAR = 0
              END IF
            END IF
            WNO = WNO + WNORES
          END DO
          IF ( IGRD .NE. NGRD+1 ) STOP 'F-OPNTAB: Logical Error'
        END IF
C
      END DO
C
      FAIL = .FALSE.
      RETURN
C
  900 FAIL = .TRUE.
      WRITE ( ERRMSG, '(A,I11)' ) 
     &  'F-OPNTAB: I/O failure on output file. IOSTAT=', IOS
      END
