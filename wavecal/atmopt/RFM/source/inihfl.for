      SUBROUTINE INIHFL ( WNOREQ, FAIL, ERRMSG )
C
C  VERSION
C     07-AUG-13  AD  Replace AI, REF (hitcom.inc) by local variable SPARE9
C     09-JAN-08  AD  Remove special case of HDO
C     14-JAN-03  AD  Add test on SHPGAS in case using LUTs for instead
C     01-DEC-02  AD  Bug fix: prevent trying to access record#0
C     13-JUN-02  AD  Allow for case no HITRAN file being used
C     11-JUN-02  AD  Simplify: only one HITRAN file. Remove IHFL argument
C     16-FEB-02  AD  Allow for HDO
C     01-OCT-96  AD  Use FPTEMP(1:14) rather than FWDPTR(1:MAXHLN) to avoid
C                    indexing bug in Sun FORTRAN compiler 
C     21-SEP-96  AD  Rename IUSVIB, ILSVIB to IUSNTE, ILSNTE
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION    
C     Initialise the HITRAN line data file.
C     Called by RFMWID and RFMFIN once for each spectral range.
C     Set the forward pointers for each required gas
C     and set current record to required wavenumber
C
      IMPLICIT NONE
C
C ARGUMENTS      
      DOUBLE PRECISION WNOREQ   !  I  Lowest wavenumber of interest [cm-1]
      LOGICAL          FAIL     !  O  Set TRUE if a fatal error is found
      CHARACTER*80     ERRMSG   !  O  Error message written if FAIL 
C
C GLOBAL CONSTANTS
      INCLUDE 'idxcon.inc' ! HITRAN molecular indices
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'hflcom.inc' ! HITRAN line data file
      INCLUDE 'hitcom.inc' ! line database communications
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'shpcon.inc' ! Line Shape codes
C
C LOCAL VARIABLES
      INTEGER IDG1      ! HITRAN gas ID for start of Fwd Pointer block
      INTEGER IDG       ! HITRAN gas ID
      INTEGER IGAS      ! Counter for gases required for RFM
      INTEGER IDUMMY    ! Dummy read
      INTEGER IOS       ! Value of IOSTAT for error messages
      INTEGER IPTR      ! Counter for forward pointers
      INTEGER IREC      ! Record#
      INTEGER FPTEMP(14) ! Forward pointers from HITRAN file
      INTEGER K,L       ! Record pointers
      INTEGER NGSREQ    ! Number of different gases required from file
      LOGICAL USEIDG(MAXPTR)    ! T=Need this Gas ID from current file
      CHARACTER*9 SPARE9 ! Spare bytes in binary file record
      DOUBLE PRECISION DDUMMY    ! Dummy variable for read
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      IF ( LUNHFL .EQ. 0 ) THEN  ! No HITRAN line data being used
        FAIL = .FALSE.
        RETURN 
      END IF
C
      WNOHFL = 0.0D0         ! Set when actual data is read in by REACYC
C       
C Binary search for 1st record: follows code in the IBM SLUP routine:
C K is first record, L is last. 
C
      K = IR1HFL
      L = IR2HFL
      DO WHILE ( K + 1 .LT. L )
        IREC = (K + L) / 2
        READ ( LUNHFL, REC=IREC, IOSTAT=IOS, ERR=900 ) 
     &    IDUMMY, IDUMMY, IDUMMY, WNUM           ! Only need WNUM for search
C
        IF ( WNUM .LT. WNOREQ ) THEN   ! Reqd WN higher, try half way up to L
          K = IREC 
        ELSE                     ! Reqd WN lower or =, try half way down to K
          L = IREC 
        ENDIF
C
c If K & L differ by 1, have found required location where the K .LT. DWNLOW 
C and L .GE. DWNLOW, Choose record L (note K<L). If there is more than 1 
C record at exactly DWNLOW, will finish pointing to first. 
C Forward pointer block are labelled with wavenumber of next line to exploit 
C this. 
C
      END DO
      IRCHFL = L
C
C Now set the initial forward pointers for each gas
C First determine which HITRAN Gas IDs are required from this file
C (IDG < MAXPTR for all reqd line molecules checked earlier in HITFIL)
      DO IDG = 1, MAXPTR
        USEIDG(IDG) = .FALSE.
      END DO
      NGSREQ = 0
      DO IGAS = 1, NGAS
        IDG = IDXGAS(IGAS)
C If IDG .LE. MAXHLN then HITRAN line molecule
C Assume SHPLUT > SHPXSC > SHPCTM, so >SHPCTM defines using LUTs for gas
        IF ( IDG .LE. MAXHLN .AND. SHPGAS(IGAS) .LT. SHPCTM ) THEN 
          IF ( .NOT. USEIDG(IDG) ) THEN   ! Only count once for each gas
            NGSREQ = NGSREQ + 1           ! although IGAS can include different
            USEIDG(IDG) = .TRUE.          ! isotopes of the same gas
          END IF
        END IF
      END DO
C
C Step back until all pointers are set 
      IREC = IRCHFL
C
C It is possible that IREC is positioned at the start of the forward pointer
C block at the start of the file, so need to advance to the end of this block
C so that pointers for each molecule can be read from the block
 100  CONTINUE
      READ ( LUNHFL, REC=IREC, IOSTAT=IOS, ERR=900 ) 
     &  LSTAT, IDG, IDUMMY
      IF ( LSTAT .EQ. -7 ) THEN
        IREC = IREC + 1
        GOTO 100
      END IF
C
      DO WHILE ( NGSREQ .GT. 0 ) 
        IREC = IREC - 1
C There should be a forward pointer for every gas 
        IF ( IREC .EQ. 0 ) STOP 'F-INIHFL: Logical error'
        READ ( LUNHFL, REC=IREC, IOSTAT=IOS, ERR=900 ) 
     &    LSTAT, IDG, IDUMMY
        IF ( LSTAT .GE. 10 .AND. USEIDG(IDG) ) THEN  ! Found a molec. line rec.
          READ ( LUNHFL, REC=IREC, IOSTAT=IOS, ERR=900 )        ! Load /HITCOM/
     &      LSTAT, IDGAS, ISO, WNUM, STREN, TPROB, ABROAD, SBROAD,
     &      ELS, ABCOEF, TSP, IUSGQ, ILSGQ, USLQ, BSLQ, SPARE9, IFWDPT
          IGAS = IGSMOL(IDG)
          IFPHFL(IGAS) = IREC + IFWDPT
          USEIDG(IDG) = .FALSE.
          NGSREQ = NGSREQ - 1
        ELSE IF ( LSTAT .EQ. -7 ) THEN         ! Found a forward pointer record
          READ ( LUNHFL, REC=IREC, ERR=900, IOSTAT=IOS )
     &      LSTAT, IDG1, IDUMMY, DDUMMY, ( FPTEMP(IPTR), IPTR = 1, 14 )
          DO IPTR = 1, 14
            IDG = IDG1 + IPTR - 1
            IF ( USEIDG(IDG) ) THEN
              IGAS = IGSMOL(IDG)
              IFPHFL(IGAS) = IREC + FPTEMP(IPTR)
              USEIDG(IDG) = .FALSE.
              NGSREQ = NGSREQ - 1
            END IF
          END DO
        END IF   
      END DO          
C
C Also use this opportunity to initialise numbers which will then remain zero
C if LTE calculation being used
      IUSNTE = 0
      ILSNTE = 0
C
      FAIL = .FALSE.
      RETURN
C
  900 WRITE ( ERRMSG, '(A,I11,A,I11)' ) 
     &  'F-INIHFL: Failed to read HITRAN file, rec#:', IREC, 
     &  '. IOSTAT=', IOS
      FAIL = .TRUE.
C
      END
