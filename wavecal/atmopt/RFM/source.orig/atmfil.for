      SUBROUTINE ATMFIL ( LUNATM, NAMATM, FIRST, GOTPRE, GOTTEM, 
     &                    GOTGAS, FAIL, ERRMSG )
C
C VERSION
C     07-AUG-13  AD  Remove redundant ATMNTE from EXTERNAL list
C     09-JAN-08  AD  Add ATMPSI to extract PSI angle
C                    Use CHKGAS to test for molecule name
C                    Remove CFCGAS and move ATMISO
C     06-JAN-04  AD  Check that PSI value is less than max absolute size
C     20-JUL-02  AD  Insist that all profiles have PSI values with GRA flag
C     06-MAR-02  AD  Correction: only check ISOGAS if line species
C     16-FEB-02  AD  Add ATMISO.
C                    Increase LABEL from C*8 to C*20
C     06-JUN-01  AD  Add ATMCHK
C     05-DEC-00  AD  Add ERRSTR to avoid variable char.str in concatenation
C     27-APR-00  AD  Remove redundant variable NEWPSI
C     20-JAN-00  AD  Allow for specification of PSI= with GRA flag.
C     21-JUN-99  AD  Bug fix: Check FAIL status after FILPRF
C     31-DEC-98  AD  Comment change only
C     27-MAR-97  AD  Add CFCGAS to convert any CFC profile labels.
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Read file containing atmospheric profiles.
C     Called by INPATM for each file listed in *ATM section of driver table.
C     The logic is
C       If First File then
C         Read No.levels and set NATM
C         Use first profile to set heights in HGTATM 
C         Copy any other profiles direct to *ATM arrays
C       Else 
C         Read No.levels and set NLEV
C         Use first profile to set heights in HGTLEV
C         Read any other profiles into PRFLEV
C         Interpolate any other profiles from PRFLEV to *ATM with HGTLEV,HGTATM
C       End if
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER       LUNATM    !  I  LUN for Atmospheric profiles 
      CHARACTER*(*) NAMATM    !  I  Name of file 
      LOGICAL       FIRST     !  I  TRUE=use these profiles to define grid
      LOGICAL       GOTTEM    ! I/O Set TRUE if Temperature profile read in
      LOGICAL       GOTPRE    ! I/O Set TRUE if Pressure profile read in
      LOGICAL       GOTGAS(*) ! I/O Set TRUE if vmr profile read in
      LOGICAL       FAIL      !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG    !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  ATMCHK ! Check atmospheric profiles are reasonable.
     &, ATMISO ! Create separate isotope profiles if specified in .atm file
     &, ATMPSI ! Extract Psi angle from brackets following .atm filename
     &, CHKGAS ! Check if GAS is valid molecule and get Hitran/RFM index
     &, FILPRF ! Read atmospheric profile from file into /ATMCOM/.
     &, LOCASE ! Convert text string to lower case.
     &, NXTREC ! Load next record from RFM input file.
     &, OPNFIL ! Open ASCII input file, and skip/Log any initial comments.
     &, PRFGRA ! Load horizontal gradient profiles into /GRACOM/.
     &, RFMLOG ! Write text message to RFM log file.
     &, TXTFLD ! Identify start and end points of text field in record
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
C
C LOCAL VARIABLES
      INTEGER      IBRACK  ! Location of '(' following .atm filename
      INTEGER      IDXVIB  ! Encoded vibrational level information
      INTEGER      IEND    ! Location of end of text field in RECORD 
      INTEGER      IGAS    ! Gas profile counter (1:NGAS)
      INTEGER      IGQ     ! Global Quantum Index (currently dummy variable)
      INTEGER      ILEV    ! Profile point counter
      INTEGER      IMOL    ! RFM/HITRAN index of recognised molecule, else 0
      INTEGER      IOS     ! Saved value of IOSTAT for error message
      INTEGER      ISO     ! Isotope# (0=no isotope specified)
      INTEGER      ISTA    ! Location of start of text field in RECORD 
      INTEGER      NLEV    ! No. levels in supplied profiles
      REAL         PSI     ! Horizontal angle
      REAL         RDUMMY  ! Dummy real for read
      LOGICAL      ANYUSE  ! T = at least one profile used from file
      LOGICAL      ENDSEC  ! T = next profile or end-of-file encountered
      LOGICAL      GETHGT  ! Set TRUE if next profile has to be *hgt 
      CHARACTER*20 LABEL   ! Label identifying profile
      CHARACTER*7  GAS     ! Possible Molecule name extracted from LABEL(2:)
      CHARACTER*80 MESSGE  ! Text message for LOG file
      CHARACTER*80 RECORD  ! Line of Driver file or External ATM file
C
C DATA STATEMENTS
      DATA PSI / 0.0 /     ! Assume psi=0 for initial prof.unless over-ridden.
      SAVE PSI
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Check for any horizontal angle PSI following .atm filename
      IBRACK = INDEX (  NAMATM, '(' )
      IF ( IBRACK .NE. 0 ) THEN
        CALL ATMPSI ( NAMATM, PSI, FAIL, ERRMSG )  ! Returns PSI=0 if not GRAFLG
        IF ( FAIL ) RETURN
        IF ( FIRST .AND. GRAFLG .AND. PSI .NE. 0.0 ) THEN   
          CALL PRFGRA ( 0.0, -2, FAIL, ERRMSG )    ! Reserve space for psi=0.0
          IF ( FAIL ) RETURN
        END IF
      ELSE IF ( GRAFLG ) THEN
        ERRMSG = 'F-ATMFIL: With GRA flag all .atm files '//
     &           'must have a horizontal angle appended'
        FAIL = .TRUE.
        RETURN
      END IF
C
      CALL OPNFIL ( LUNATM, NAMATM, FAIL, ERRMSG ) 
      IF ( FAIL ) RETURN
      ANYUSE = FIRST
      GETHGT = .TRUE.
C
C 1st non-comment record in file must be no. levels (integer)
C 
      READ ( LUNATM, *, IOSTAT=IOS, ERR=900 ) NLEV
C
C Jump here for each new profile in file. Next non-comment record should start 
C with a '*' (ie '*HGT','*TEM', '*PRE', or *species) which results 
C in NXTREC backspacing and setting ENDSEC=TRUE. End-of-file results in same.
C
  100 CALL NXTREC ( LUNATM, RECORD, ENDSEC, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      IF ( .NOT. ENDSEC ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-ATMFIL: Unexpected record: '//RECORD(1:48)//'...'
        RETURN
      END IF
C
C End-of-file is OK as exit condition (alternative to '*END' marker)
C
      READ ( LUNATM, '(A80)', IOSTAT=IOS, END=800, ERR=900 ) RECORD
      CALL TXTFLD ( 80, RECORD, 1, ISTA, IEND )  ! '*' guarantees ISTA .GE. 1
      IEND = MIN ( ISTA+19, IEND )   ! Just in case incorrect label is > C*20
      LABEL = RECORD(ISTA:IEND)
      CALL LOCASE ( LABEL, LABEL )   ! All internal species codes lower case
C
C First profile in each file (GETHGT=TRUE) should be altitude
C
      IF ( GETHGT ) THEN
        IF ( LABEL .EQ. '*hgt' ) THEN
          CALL RFMLOG ( '    '//RECORD, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
          CALL FILPRF ( LUNATM, FIRST, NLEV, -2, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        ELSE
          FAIL = .TRUE.
          ERRMSG = 'F-ATMFIL: 1st profile wasn''t *HGT, it was:'//LABEL
          RETURN
        END IF
        GETHGT = .FALSE.
C Check for *END marker
      ELSE IF ( LABEL .EQ. '*end' ) THEN 
        GOTO 800
C Check for pressure profile
      ELSE IF ( LABEL .EQ. '*pre' ) THEN                 
        CALL RFMLOG ( '    '//RECORD, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        CALL FILPRF ( LUNATM, FIRST, NLEV, -1, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        IF ( GRAFLG ) THEN
          CALL PRFGRA ( PSI, -1, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        ELSE IF ( GOTPRE ) THEN                          ! Overwritten previous
          MESSGE = 'W-ATMFIL: Pressure profile superseded'
          CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        END IF
        GOTPRE = .TRUE.                               ! Flag as loaded
        ANYUSE = .TRUE.
C Check for temperature profile
      ELSE IF ( LABEL .EQ. '*tem' ) THEN           
        CALL RFMLOG ( '    '//RECORD, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        CALL FILPRF ( LUNATM, FIRST, NLEV, 0, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        IF ( GRAFLG ) THEN
          CALL PRFGRA ( PSI, 0, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        ELSE IF ( GOTTEM ) THEN                          ! Overwritten previous
          MESSGE = 'W-ATMFIL: Temperature profile superseded'
          CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        END IF
        GOTTEM = .TRUE.                                  ! Flag as loaded
        ANYUSE = .TRUE.
C Check for recognised absorbing molecule from those loaded in *GAS section
      ELSE                                              
        CALL CHKGAS ( LABEL(2:), GAS, IMOL, ISO, IGQ, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        IF ( IMOL .EQ. 0 ) THEN                      ! Not a recognised molecule
          READ ( LUNATM, *, IOSTAT=IOS, ERR=900 )    ! Read past profile data
     &          ( RDUMMY, ILEV = 1, NLEV)
        ELSE IF ( IGSMOL(IMOL) .EQ. 0 ) THEN         ! Not using this molecule
          READ ( LUNATM, *, IOSTAT=IOS, ERR=900 )    ! Read past profile data
     &          ( RDUMMY, ILEV = 1, NLEV)
        ELSE IF ( IGQ .NE. 0 .AND. NTEFLG ) THEN     ! VT Profile
          CALL RFMLOG ( '    '//RECORD, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
          IDXVIB = IMOL * 1000000 + ISO * 1000 + IGQ
          CALL FILPRF ( LUNATM, FIRST, NLEV, IDXVIB, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
          ANYUSE = .TRUE.
        ELSE                                         ! Use this VMR profile

          IGAS = IGSMOL(IMOL)
C Establish separate atmospheric profile for this isotopomer if required
          IF ( ISO .GT. 0 ) THEN                     ! Isotopic profile
            CALL ATMISO ( IMOL, ISO, FAIL, ERRMSG )  ! if new, increase NGAS
            IF ( FAIL ) RETURN
            IGAS = ISOGAS(ISO,IGAS)                  ! Isotopic profile
          END IF
C
          CALL RFMLOG ( '    '//RECORD, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
          CALL FILPRF ( LUNATM, FIRST, NLEV, IGAS, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
          IF ( GRAFLG ) THEN
            CALL PRFGRA ( PSI, IGAS, FAIL, ERRMSG )
            IF ( FAIL ) RETURN
          ELSE IF ( GOTGAS(IGAS) ) THEN                              
            MESSGE =
     &        'W-ATMFIL: '//CODGAS(IGAS)//' profile superseded'
            CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
            IF ( FAIL ) RETURN
          END IF
          GOTGAS(IGAS) = .TRUE.                               
          ANYUSE = .TRUE.
        END IF              
      END IF               
      GOTO 100             ! Get next profile from file
C
  800 CONTINUE
      IF ( .NOT. ANYUSE ) THEN
        MESSGE = 'W-ATMFIL: No useful profiles found in file'
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
      CLOSE ( LUNATM, IOSTAT=IOS, ERR=900 )
      FAIL = .FALSE.
C
C Check currently loaded profiles are reasonable
      CALL ATMCHK ( GOTTEM, GOTPRE, GOTGAS, FAIL, ERRMSG ) 
      RETURN
C
  900 WRITE ( ERRMSG, '(A,I11)' )
     &  'F-ATMFIL: I/O failure reading Atm.Profile file. IOSTAT=', IOS
      FAIL = .TRUE.
C
      END
