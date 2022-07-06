      SUBROUTINE FILXFL ( LUNXSC, REJECT, FAIL, ERRMSG )
C
C VERSION
C     14-AUG-13  AD  Rewritten and simplified using XFLHDR
C                    No longer runs TRIANG to check triangulation
C     10-DEC-01  AD  Ensure JGAS is initialised
C     30-JAN-01  AD  Allow for HITRAN 2000 format files
C     04-JAN-00  AD  Allow for AVG option.
C     07-OCT-99  AD  Set DWNILS=1.0 instead of WIDILS
C     08-SEP-99  AD  Check for each spectral range, and if TRIANG works.      
C                    Also check NXSP for each p,T table
C     03-MAR-97  AD  Version 3.
C     16-JAN-97  AD  Set GASXFL
C     13-JAN-97  AD  Correction to check on NXFL v MAXXFL (was v.MAXHFL)
C     01-OCT-96  AD  Version 2.
C     17-SEP-96  AD  Remove redundant IRCXFL
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Check X/S data input file and store in /XFLCOM/
C     Called by XSCFIL for each new file.
C     Note that X/S data only required within widemesh grid so test for range
C     WMNSPC, WMXSPC
C
      IMPLICIT NONE 
C
C ARGUMENTS
      INTEGER      LUNXSC   !  I  LUN for current X/S file
      LOGICAL      REJECT   !  O  T=Reject this X/S file
      LOGICAL      FAIL     !  O  T=A fatal error was detected
      CHARACTER*80 ERRMSG   !  O  Error message written if FAIL is TRUE
C      
      EXTERNAL
     &  IDGNEW ! Convert old RFM index for .xsc data to new value
     &, RFMLOG ! Write text message to RFM log file.
     &, XFLHDR ! Extract header info for next spectral range from .xsc files
      INTEGER IDGNEW
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags 
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
      INCLUDE 'xflcom.inc' ! X-section data files
      INCLUDE 'xsccom.inc' ! X-section data 
C
C LOCAL VARIABLES      
      INTEGER          IDGAS    ! ID of first gas in forward pointer record
      INTEGER          IDGOLD   ! Old value of IDGAS read from file
      INTEGER          IGAS     ! Gas counter
      INTEGER          IOS      ! Saved value of IOSTAT for error message
      INTEGER         IREC1     ! Record# for start of spectral range
      INTEGER          ISPC     ! Counter for RFM spectral ranges
      INTEGER    IXSC, JXSC     ! Indices of X/S data sets
      INTEGER          NREC     ! No.records in current spectral range
      INTEGER          NXSP     ! Total No.X/S data points in each X/S data set
      INTEGER          NXST     ! No. of (p,T) tables X/S data set 
      DOUBLE PRECISION DWNILS   ! Width to allow for ILS function (if any)
      DOUBLE PRECISION WNUXSC   ! Max Wno [/cm] for X/S data
      DOUBLE PRECISION WNLXSC   ! Min Wno [/cm] for X/S data
      LOGICAL          USEXSC   ! T=Use X/S data, F= don't use
      LOGICAL          WARN     ! T=warn if old index, F=don't warn
      CHARACTER*100    HEADER   ! Header Text for X/S data
      CHARACTER*80     MESSGE   ! Text message for Log file
C
C EXECUTABLE CODE --------------------------------------------------------------

      REJECT = .TRUE.
      WARN = .TRUE.
      IF ( ILSFLG .OR. AVGFLG ) THEN
        DWNILS = 1.0D0       ! Allow extra width for ILS convolution
      ELSE
        DWNILS = 0.0D0
      END IF
C
C Find which molecules and wavenumber ranges are contained in file
      IREC1 = 1
 100  CONTINUE
        CALL XFLHDR ( LUNXSC, IDGAS, WNLXSC, WNUXSC, 
     &                  NXST, NXSP, NREC, HEADER, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        IF ( NREC .EQ. 0 ) GOTO 300         ! reached end-of-file, exit 100 loop
C
C Convert to new RFM index (>100 for x/s files) if required 
        IF ( IDGAS .LT. 100 ) THEN 
          IDGOLD = IDGAS
          IDGAS = IDGNEW ( IDGOLD )
          IF ( IDGAS .GT. 100 .AND. WARN ) THEN ! successful conversion
            WARN = .FALSE.      ! just one warning per file
            WRITE ( MESSGE, '(A,I2,A,I3)' )
     &        'W-FILXFL: File contains old molecular index=', IDGOLD,
     &        ' - assigning new index=', IDGAS
            CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
            IF ( FAIL ) RETURN
          END IF
        END IF
C
C Check that IDGAS is in the expected range for cross-section data
        IF ( IDGAS .LT. 100 .OR. IDGAS .GT. MAXMOL ) THEN
          WRITE ( ERRMSG, '(A,I11,A,I4)' )
     &      'F-FILXFL: File contains Molec#ID=', IDGAS,
     &      ' out of valid range 100:', MAXMOL
          FAIL = .TRUE.
        END IF
C
C See if x/s molecule and spectral range are required by the RFM
        USEXSC = .FALSE.
        IGAS = IGSMOL(IDGAS)
        IF ( IGAS .GT. 0 ) THEN  ! this molecule is required by RFM
          GASXFL(IGAS) = .TRUE.  ! at least one x/s file supplied for this gas
          DO ISPC = 1, NSPC
            IF ( WNLXSC .LE. WNUSPC(ISPC) + DWNILS .AND.
     &           WNUXSC .GE. WNLSPC(ISPC) - DWNILS       ) 
     &        USEXSC = .TRUE.    ! data required for at least one spect.range
          END DO
        END IF
C
C Data will be used 
        IF ( USEXSC ) THEN
          MESSGE = 'I-FILXFL: Using x/s data: '//HEADER(1:50)
          CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
C
C Using this spectral range, so check relevant array sizes adequate
          FAIL = .TRUE.
          IF ( NXSC .EQ. MAXXSC ) THEN 
            ERRMSG = 'F-FILXFL: No.X-Sec data arrays reqd > '//
     &        'MAXXSC in rfmsiz.inc'          
          ELSE IF ( NXST .GT. MAXXST ) THEN
            WRITE ( ERRMSG, '(A,I11,A,I11,A)' )
     &      'F-FILXFL: No.X/S (p,T) tables =', NXST,
     &      ' > MAXXST =', MAXXST, ' in rfmsiz.inc'
          ELSE IF ( NXSP .GT. MAXXSP ) THEN
            WRITE ( ERRMSG, '(A,I11,A)' )
     &      'F-FILXFL: No.X/S data pts =', NXSP, 
     &      ' > MAXXSP =', MAXXSP, ' in rfmsiz.inc' 
          ELSE
            FAIL = .FALSE.
          END IF
          IF ( FAIL ) RETURN
C
C Since at least one spectral range required, mark this file as useful
          REJECT = .FALSE.
C
C Insert details into list ordered by lower wavenumber
C Note that chkxsc.for will later check for any ambiguous overlaps
          NXSC = NXSC + 1
          JXSC = NXSC
          DO IXSC = 1, NXSC - 1
            IF ( WN1XSC(IXSC) .GT. WNLXSC ) THEN
              DO JXSC = NXSC, IXSC + 1, -1 
                IFLXSC(JXSC) = IFLXSC(JXSC-1) 
                IGSXSC(JXSC) = IGSXSC(JXSC-1) 
                WN1XSC(JXSC) = WN1XSC(JXSC-1) 
                WN2XSC(JXSC) = WN2XSC(JXSC-1) 
                IRCXSC(JXSC) = IRCXSC(JXSC-1) 
              END DO
              JXSC = IXSC
              GOTO 200
            END IF
          END DO                    
  200     IFLXSC(JXSC) = NXFL + 1   ! NXFL doesn't yet include this file
          IGSXSC(JXSC) = IGAS
          WN1XSC(JXSC) = WNLXSC
          WN2XSC(JXSC) = WNUXSC
          IRCXSC(JXSC) = IREC1
        END IF
        IREC1 = IREC1 + NREC
      GOTO 100
C
  300 IF ( REJECT ) THEN
        MESSGE = 
     &  'W-FILXFL: Ignore X/S file, no reqd gases in reqd waveno.range'
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      ELSE 
        IF ( NXFL .EQ. MAXXFL ) THEN
          FAIL = .TRUE.
          WRITE ( ERRMSG, '(A,I4,A)' ) 
     &      'F-FILXFL: No.different X/S files reqd > MAXXFL=', MAXXFL,
     &      ' in rfmsiz.inc'
          RETURN
        ELSE
          REWIND ( LUNXSC, IOSTAT=IOS, ERR=900 )
          NXFL = NXFL + 1
          LUNXFL(NXFL) = LUNXSC
        END IF
      END IF
      FAIL = .FALSE.
      RETURN
C
  900 FAIL = .TRUE.
      WRITE ( ERRMSG, '(A,I11)' ) 
     &  'F-FILXFL: I/O failure on X/S file. IOSTAT=', IOS

      END
