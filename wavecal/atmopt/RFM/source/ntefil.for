       SUBROUTINE NTEFIL ( LUNNTE, NAMNTE, FAIL, ERRMSG )
C
C VERSION
C     05-JUN-09  AD  Allow free-format for record beginning '***'
C     20-FEB-01  AD  Extend VT profiles to top of atmosphere if necessary
C                    (previously regarded as fatal error condition)
C     09-AUG-99  AD  Rename dummy variable NISO to IDUMMY
C     01-JUL-98  AD  Force supplied NTE profiles to reach top of atmosphere
C                    Remove IATMHI argument to subroutines.
C     04-APR-97  AD  Correct format of NLEV>MAXLEV error message.
C                    Change MAXLEV from 100 to 121 
C                    Read in NQPR rather than calculate
C                    Add QFNCOM.INC to initialise NQFN
C                    Remove FIRST argument to NTEQFN
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     21-SEP-96  AD  Original. Based on GENLN2 module NTEPRO.
C
C DESCRIPTION
C     Read NTE Vibrational Temperature Data from file
C     Called by INPNTE for each file specified in *NTE section.
C
C     File format is:
C     *** MODEL NSET NGAS NLEV NQPR (A3,5I6)
C     %then IG=1,NGAS lines of:
C     JGAS(IG) NISO(IG) (3X,2I6)
C     %then pressure levels (mb)
C     PLEV(ILEV),ILEV=1,NLEV (5(1PE10.3))    
C     %then kinetic temperature (K)
C     IS  (=0)(I6)
C     TK(ILEV),ILEV=1,NLEV (5F10.3)
C     %then for each vib.part.func.:
C     JGAS(IG) JISO(IS,IG) (2I6)
C     QV(ILEV,IS,IG),ILEV=1,NLEV (5(1PE10.3))
C     %then ISET=1,NSET sections of:
C     ISET IGAS ISO IVSGQ ENERGY (4I6,F12.4)
C     TVIB(IL,ISET),ILEV=1,NLEV (5F10.3)    
C
      IMPLICIT NONE 
C
C ARGUMENTS
      INTEGER       LUNNTE !  I  Lun for NTE file
      CHARACTER*(*) NAMNTE !  I  Name of NTE file
      LOGICAL       FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  NTEQFN ! Store Vib.Partition Function data read from NTE file
     &, NTETEM ! Store Vib.Temp profiles read from NTE file.
     &, NXTREC ! Load next record from RFM input file.
     &, OPNFIL ! Open ASCII input file, and skip/Log any initial comments.
     &, RFMLOG ! Write text message to RFM log file.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'ntecom.inc' ! Non-LTE data
      INCLUDE 'qfncom.inc' ! Non-LTE Qvib profiles
C
C LOCAL CONSTANTS
      INTEGER       MAXLEV         ! Max.no levels in supplied profiles
        PARAMETER ( MAXLEV = 201 ) ! Should be comparable or .GT. MAXATM
C 
C LOCAL VARIABLES
      INTEGER      IATMLO   ! Lower values of IATM for interpoln.
      INTEGER      IDISO    ! Istope# 
      INTEGER      IDMOL    ! HITRAN ID code for molecule
      INTEGER      IDUMMY   ! Dummy integer 
      INTEGER      ILEV     ! Counter for File profile levels
      INTEGER      IMOL     ! Counter for different species
      INTEGER      IOS      ! IOSTAT saved for error message
      INTEGER      IQPR     ! Counter for Partition Fn profiles in current file
      INTEGER      ISET     ! Counter for vib.temp profiles in current file
      INTEGER      IVSGQ    ! Global Quantum No of Vibrational state
      INTEGER      LNTE     ! Last value of NNTE (to check for any change)
      INTEGER      NLEV     ! No.of profile levels in file
      INTEGER      NQPR     ! Total No.of different partition functions in file
      INTEGER      NMOL     ! No.of different species in file
      INTEGER      NSET     ! No.of Vib.Temp profiles in current file
      REAL         ENERGY   ! Energy of vibrational state
      REAL         LNPLEV(MAXLEV) ! Vertical scale converted to log(p/atm)
      REAL         QFNLEV(MAXLEV) ! Partition Fn profile read from file
      REAL         TKNLEV(MAXLEV) ! Kinetic temperature [K] profile from file
      REAL         TVBLEV(MAXLEV) ! Vib. temperature [K] profile from file
      LOGICAL      ENDSEC   ! Set TRUE when '*' marker or EOF reached in file
      CHARACTER*80 MESSGE   ! Text message sent to Runlog file
      CHARACTER*80 RECORD   ! Record read from file
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      FAIL = .FALSE.
      LNTE = NNTE
C
      CALL OPNFIL ( LUNNTE, NAMNTE, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
C Find model header record
C 
      ENDSEC = .FALSE.
      DO WHILE ( .NOT. ENDSEC )
        CALL NXTREC ( LUNNTE, RECORD, ENDSEC, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END DO
C
C First three characters of RECORD should be '***'
      READ ( LUNNTE, '(A)', IOSTAT=IOS, ERR=900 ) RECORD
      READ ( RECORD(4:), *, IOSTAT=IOS, ERR=900 ) 
     &  IDUMMY, NSET, NMOL, NLEV, NQPR         ! IDUMMY is actually MODEL#
C
C Read in gas, isotope pairs for which there is vib. tem. data
C
      DO IMOL = 1, NMOL 
        READ ( LUNNTE, *, IOSTAT=IOS, ERR=900 ) IDMOL, IDUMMY ! IDUMMY=NISO
      END DO
C
C Check NLEV fits local array space
      IF ( NLEV .GT. MAXLEV ) THEN
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,I11,A,I6)' )
     &    'F-NTEFIL: No.profile levels =', NLEV,
     &    ' > Local array dimension MAXLEV=', MAXLEV
        RETURN
      ELSE IF ( NLEV .LT. 1 ) THEN
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,I11)' )
     &    'F-NTEFIL: Value of NLEV must be GE 1, actual value=', NLEV
        RETURN
      END IF
C
C Read in pressure profile for this model and convert from p/mb to ln(p/mb)
      READ ( LUNNTE, *, IOSTAT=IOS, ERR=900 ) 
     &  ( LNPLEV(ILEV), ILEV = 1, NLEV )                   ! p in mb
      DO ILEV = 1, NLEV
        LNPLEV(ILEV) = LOG ( LNPLEV(ILEV) ) 
      END DO   
C
C Check there is some overlap between pressure profile and local atmos.profile
      IF ( LNPLEV(1) .LT. LNPATM(NATM) ) THEN
        MESSGE = 
     &   'W-NTEFIL: ignoring NTE file - '//
     &   'above specified ATMosphere for calculations'
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        GOTO 800
      ELSE IF ( LNPLEV(NLEV) .GT. LNPATM(NATM) ) THEN
        MESSGE = 'W-NTEFIL: Extending NTE VT profiles '//
     &           'to top of specified ATMosphere'
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF  
      IATMLO = 1
      DO WHILE ( LNPATM(IATMLO) .GT. LNPLEV(1) )
        IATMLO = IATMLO + 1
      END DO
C
      READ ( LUNNTE, *, IOSTAT=IOS, ERR=900 ) ISET     ! should be ISET=0 
      READ ( LUNNTE, *, IOSTAT=IOS, ERR=900 )          ! kinetic temp.profile
     &  ( TKNLEV(ILEV), ILEV = 1, NLEV )
C
C Read vibrational partition function profiles
C
      NQFN = 0
      IF ( NQPR .GT. 0 ) THEN
        CALL NXTREC ( LUNNTE, RECORD, ENDSEC, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        BACKSPACE ( LUNNTE, IOSTAT=IOS, ERR=900 )
        DO IQPR = 1, NQPR
          READ ( LUNNTE, *, IOSTAT=IOS, ERR=900 ) IDMOL, IDISO
          READ ( LUNNTE, *, IOSTAT=IOS, ERR=900 ) 
     &      ( QFNLEV(ILEV), ILEV = 1, NLEV )
          CALL NTEQFN ( IDMOL, IDISO, IATMLO, 
     &                NLEV, LNPLEV, QFNLEV, FAIL, ERRMSG )  
          IF ( FAIL ) RETURN
        END DO
      END IF
C
C Read vibrational temperature profiles 
C       
      CALL NXTREC ( LUNNTE, RECORD, ENDSEC, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      BACKSPACE ( LUNNTE, IOSTAT=IOS, ERR=900 )
      DO ISET = 1, NSET
        READ ( LUNNTE, *, IOSTAT=IOS, ERR=900 ) 
     &    IDUMMY, IDMOL, IDISO, IVSGQ, ENERGY               ! IDUMMY = ISET
        READ ( LUNNTE, *, IOSTAT=IOS, ERR=900 )
     &    ( TVBLEV(ILEV), ILEV = 1, NLEV )
        DO ILEV = 1, NLEV
          TVBLEV(ILEV) = TVBLEV(ILEV) - TKNLEV(ILEV)
        END DO
        CALL NTETEM ( IDMOL, IDISO, IVSGQ, IATMLO, 
     &                NLEV, LNPLEV, TVBLEV, ENERGY, FAIL, ERRMSG ) 
        IF ( FAIL ) RETURN
C
      END DO
C
C Check if file contained any useful data (NNTE will have been incremented)
      IF ( NNTE .EQ. LNTE ) THEN
        MESSGE = 
     &   'W-NTEFIL: ignoring NTE file - contains no useful data'
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
C
  800 CONTINUE
      CLOSE ( LUNNTE, IOSTAT=IOS, ERR=900 )
      FAIL = .FALSE.
      RETURN
C
  900 FAIL = .TRUE.
      WRITE ( ERRMSG, '(A,I11)' )
     &  'F-NTEFIL: I/O error on NTE data. IOSTAT=', IOS
      END
