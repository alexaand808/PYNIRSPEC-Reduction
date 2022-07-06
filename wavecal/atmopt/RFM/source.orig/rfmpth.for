      SUBROUTINE RFMPTH ( LUNPTH, FAIL, ERRMSG )
C
C VERSION
C     23-APR-12  AD  Use RFLSFC rather than EMSSFC to set OBS path
C     19-OCT-11  AD  Bug#80: Initialise NSEG2
C     16-JAN-04  AD  Set IDIR argument to IDXPTH=0 unless GRAFLG=T
C     06-JAN-04  AD  Add IGRA argument to MAKNAM
C     28-APR-03  AD  Correction: if OBS flag start with IATM1 = IATOBS-1 
C     29-MAY-02  AD  Correction: if GRA flag enabled, write diagnostics for
C                    both up and down paths for second and subsequent absorbers
C     05-DEC-00  AD  Comment out CARRIAGECONTROL='LIST' 
C     20-JAN-00  AD  GRA option. Add extra argument to IDXPTH
C                    Remove 'X' as 1st char from write statement formats
C     09-AUG-99  AD  Add SNGL() to ASIN(SZNTAN) in write statements.
C     29-JUL-98  AD  Add extra header records giving viewing geometry 
C     14-JUL-98  AD  Comment change only. 
C     11-MAR-97  AD  Corrections for use with NAD/ZEN flags
C     10-MAR-97  AD  Add CARRIAGECONTROL='LIST'
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Open PTH files and write RFM ray path diagnostics.
C     Called once by RFM if PTH flag selected.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      LUNPTH  !  I  Next available LUN
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error occurs
      CHARACTER*80 ERRMSG  !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  IDXPTH ! Indexing function for Pth# from Tan.pt#, Atm.Layer# and Gas#
     &, MAKNAM ! Construct filename for RFM output files.
     &, RFMLOG ! Write text message to RFM log file.
      INTEGER IDXPTH
C
C GLOBAL CONSTANTS
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'crvcom.inc' ! Local radius of curvature
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'obscom.inc' ! Observer Position
      INCLUDE 'outcom.inc' ! RFM output file data
      INCLUDE 'pthcom.inc' ! RFM path data
      INCLUDE 'sfccom.inc' ! Surface parameters
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL VARIABLES
      INTEGER       IATM   ! Atmospheric layer counter
      INTEGER IATM1, IATM2 ! Limits for atmospheric level counter
      INTEGER       ICLC   ! Flag for whether LBL calc.performed on path
      INTEGER       IDIR   ! Direction pointer
      INTEGER       IGAS   ! Absorber counter
      INTEGER       IOS    ! Value of IOSTAT saved for error messages.
      INTEGER       IPTH   ! Counter for RFM paths
      INTEGER       ITAN   ! Tangent height counter
      INTEGER NSEG1, NSEG2 ! No. segments for each tangent path.
      REAL          AMTSUM ! Total absorber mass [kmol/cm^2]
      REAL          HGT    ! Lowest altitude for path segment
      REAL          RAYSUM ! Total pathlength [km]
      CHARACTER*132 MESSGE ! Text message sent to Log file.
      CHARACTER*80  NAMFIL ! Name of file actually opened (incl. RUNID)
      CHARACTER*80  RECORD ! Text record written out
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      DO ITAN = 1, NTAN
C
C Construct filename and open file
C
        CALL MAKNAM ( 0, 0, ITAN, 0, NAMPTH, NAMFIL )
        MESSGE = 'I-RFMPTH: Opening output file: '//NAMFIL
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        IF ( NEWFLG ) THEN
          OPEN ( UNIT=LUNPTH, FILE=NAMFIL, STATUS='NEW', 
     &           IOSTAT=IOS, ERR=900 )
        ELSE
          OPEN ( UNIT=LUNPTH, FILE=NAMFIL, STATUS='UNKNOWN', 
     &           IOSTAT=IOS, ERR=900 )
        END IF
C
C Write File Header 
C
        IF ( HOMFLG ) THEN 
          WRITE ( LUNPTH, '(A,F8.3,A)', IOSTAT=IOS, ERR=900 ) 
     &      '! Path diagnostics calculated for Length =', 
     &      USRTAN(ITAN), ' [km] by RFM v.'//VERSN
          WRITE ( LUNPTH, '(A,A)', IOSTAT=IOS, ERR=900 ) '!',TXTHDR
          IATM1 = 1
          IATM2 = 1
          IDIR  = 1
          NSEG1 = 1
          NSEG2 = 0
        ELSE IF ( ZENFLG .OR. NADFLG ) THEN
          IF ( USRELE ) THEN
            WRITE ( LUNPTH, '(A,F8.3,A)', IOSTAT=IOS, ERR=900 ) 
     &        '! Path diagnostics calculated for Ele.Ang=', 
     &        USRTAN(ITAN), ' [dg] by RFM v.'//VERSN
          ELSE
            WRITE ( LUNPTH, '(A,F8.3,A)', IOSTAT=IOS, ERR=900 ) 
     &        '! Path diagnostics calculated for AirMass=', 
     &        USRTAN(ITAN), '      by RFM v.'//VERSN
          END IF
          WRITE ( LUNPTH, '(A,A)', IOSTAT=IOS, ERR=900 ) '!',TXTHDR
          IF ( ZENFLG ) THEN
            IATM1 = 1
            IATM2 = NATM - 1
            IF ( OBSFLG ) IATM1 = MIN ( IATOBS, NATM-1 )
            IDIR = 1
            NSEG1 = IATM2 - IATM1 + 1
            NSEG2 = 0
          ELSE                              ! NADFLG
            IATM1 = NATM - 1
            IATM2 = 1
            IF ( OBSFLG .AND. .NOT. RFLSFC ) 
     &        IATM1 = MIN ( IATOBS-1, NATM-1 ) 
            IDIR = -1
            NSEG1 = IATM1 - IATM2 + 1
            NSEG2 = 0
          END IF
        ELSE
          IF ( USRELE ) THEN
            WRITE ( LUNPTH, '(A,F8.3,A)', IOSTAT=IOS, ERR=900 ) 
     &        '! Path diagnostics calculated for Ele.Ang=', 
     &        ELETAN(ITAN), ' [dg] by RFM v.'//VERSN
          ELSE IF ( USRGEO ) THEN
            WRITE ( LUNPTH, '(A,F8.3,A)', IOSTAT=IOS, ERR=900 ) 
     &        '! Path diagnostics calculated for Geo.Tan=', 
     &        GEOTAN(ITAN), ' [km] by RFM v.'//VERSN
          ELSE
            WRITE ( LUNPTH, '(A,F8.3,A)', IOSTAT=IOS, ERR=900 ) 
     &        '! Path diagnostics calculated for Tan Hgt=', 
     &        HGTTAN(ITAN), ' [km] by RFM v.'//VERSN
          END IF
          WRITE ( LUNPTH, '(A,A)', IOSTAT=IOS, ERR=900 ) '!',TXTHDR
C For limb-viewing, write out extra couple of header records giving additional
C information on viewing/tangent point geometry
          WRITE ( LUNPTH, '(A)', IOSTAT=IOS, ERR=900 )
     &      '!  Rfr.Tan   Geo.Tan   Tan.Zen   Tan.Psi   Rad.Crv   '//
     &      'Obs.Ele   Obs.Alt   Obs.Psi'
          RECORD = '!'
          WRITE ( RECORD(2:10), '(F9.3)' ) HGTTAN(ITAN)
          WRITE ( RECORD(11:20), '(F10.3)' ) GEOTAN(ITAN)
          WRITE ( RECORD(21:30), '(F10.3)' ) 
     &      SNGL(ASIN(SZNTAN(ITAN)))/DTORAD
          IF ( GRAFLG ) WRITE ( RECORD(31:40), '(F10.3)' ) PSITAN(ITAN)
          WRITE ( RECORD(41:50), '(F10.3)' ) RADCRV
          IF ( OBSFLG ) THEN
            WRITE ( RECORD(51:60), '(F10.3)' ) ELETAN(ITAN)
            WRITE ( RECORD(61:70), '(F10.3)' ) ALTOBS
            IF ( GRAFLG ) WRITE ( RECORD(71:80), '(F10.3)' ) PSIOBS
          END IF
          WRITE ( LUNPTH, '(A)' ) RECORD
          IATM1 = NATM-1
          IATM2 = IATTAN(ITAN)
          IDIR = -1
          NSEG1 = IATM1 - IATM2 + 1  ! top of atmosphere to tangent level
          NSEG2 = 0                              ! Bug#80 fix: add this line
          IF ( GRAFLG ) THEN  ! Allow diagnostics for both down and up parts 
            IF ( OBSFLG ) THEN
              IATM1 = MIN ( IATOBS-1, NATM-1 ) 
              IF ( IATTAN(ITAN) .EQ. IATOBS ) THEN ! Upview, skip downward loop
                IATM1 = IATOBS
                IATM2 = IATM1 + 1               
              ELSE
                NSEG2 = NSEG1
                NSEG1 = IATM1 - IATM2 + 1  ! Add obs to tan.level
              END IF
            ELSE
              NSEG2 = NSEG1
            END IF
          END IF
        END IF
C
        WRITE ( LUNPTH, * ) NGAS, NSEG1, NSEG2, ' = NGas, NSeg1, NSeg2'
C
C For each absorber....
C
        DO IGAS = 1, NGAS
          WRITE ( LUNPTH, '(A)', IOSTAT=IOS, ERR=900 ) CODGAS(IGAS)
          IF ( GRAFLG ) THEN
            WRITE ( LUNPTH, '(A)', IOSTAT=IOS, ERR=900 ) 
     &        'Lev  Zlow[km]  Psi[dg] Temp[K]  Press[mb]  VMR[ppv]'//
     &        '   Amt[kmol/cm2] Len.[km] Clc'
          ELSE
            WRITE ( LUNPTH, '(A)', IOSTAT=IOS, ERR=900 ) 
     &        'Lev  Zlow[km]  Zen[dg] Temp[K]  Press[mb]  VMR[ppv]'//
     &        '   Amt[kmol/cm2] Len.[km] Clc'
          END IF
C
  100     CONTINUE
          RAYSUM = 0.0
          AMTSUM = 0.0
          DO IATM = IATM1, IATM2, IDIR
            IF ( GRAFLG ) THEN
              IPTH = IDXPTH ( ITAN, IATM, IGAS, IDIR )
            ELSE
              IPTH = IDXPTH ( ITAN, IATM, IGAS, 0 )
            END IF
            IF ( CLCPTH(IPTH) ) THEN
              ICLC = 1
            ELSE
              ICLC = 0
            END IF
            IF ( IATM .EQ. IATTAN(ITAN) ) THEN
              HGT = HGTTAN(ITAN)
            ELSE
              HGT = HGTATM(IATM)
            END IF
            WRITE ( LUNPTH, '(I3, 3F9.3, 1P3E12.5, 0PF10.3, I3)',
     &                IOSTAT=IOS, ERR=900 )
     &        IATM, HGT, PSIPTH(IPTH), TEMPTH(IPTH), PREPTH(IPTH)*ATMB,
     &        PPAPTH(IPTH)/PREPTH(IPTH), AMTPTH(IPTH),
     &        RAYPTH(IPTH), ICLC
            RAYSUM = RAYSUM + RAYPTH(IPTH)
            AMTSUM = AMTSUM + AMTPTH(IPTH)
          END DO
          WRITE ( LUNPTH, '(48X,A,1PE12.5,0PF10.3)', IOSTAT=IOS, 
     &      ERR=900 ) 'Total:', AMTSUM, RAYSUM
          IF ( GRAFLG ) THEN
            IF ( IDIR .EQ. -1 ) THEN ! Repeat from tan.pt to t.o.a.
              IDIR = +1
              IATM1 = IATTAN(ITAN)
              IATM2 = NATM - 1
              GOTO 100
            ELSE                     ! IDIR=+1, ie tan.pt to t.o.a.
              IDIR = -1              ! Reset downward path for next IGAS
              IATM1 = NATM - 1
              IATM2 = IATTAN(ITAN)
            END IF            
          END IF            
        END DO
        CLOSE ( LUNPTH, IOSTAT=IOS, ERR=900 )
      END DO
      FAIL = .FALSE.
      RETURN
C
  900 FAIL = .TRUE.
      WRITE ( ERRMSG, '(A,I11)' ) 
     &  'F-RFMPTH: I/O failure on output file. IOSTAT=', IOS
      END
