      SUBROUTINE LUTCHK ( LUN, IGAS, ISPC, NL, NPT, NWN, NHEAD,
     &                    BINFIL, IRRFIL, TAB, QAL, FAIL, ERRMSG )
C
C VERSION
C     06-JUN-00  AD  Add IRRFIL argument, allow for irreg.grid files
C     31-JUL-99  AD  Remove unused variable NYTOT
C     04-JUN-99  AD  Add BINFIL, NHEAD arguments to allow for binary files
C     29-APR-99  AD  Allow qualifier data appended to TAB files
C     23-APR-99  AD  Change error messages to accommodate C*8 LABSPC
C     16-AUG-98  AD  Correction to error message for >MAXLUL
C     19-JUN-98  AD  Modify to handle uncompressed LUTs differently.
C     22-JAN-98  AD  (Standard F77) Put DATA statements *after* declarations
C     20-JAN-98  AD  Add QAL argument and check value
C     18-JUL-97  AD  Add TAB argument and check value
C     08-JUL-97  AD  Remove redundant IOS.
C     04-JUL-97  AD  Original.
C
C DESCRIPTION
C     Check Look-Up Table data.
C     Called by LUTFIL for each usaful file specified in *LUT section.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER       LUN     !  I  LUN for LUT File (stored, not used for I/O)
      INTEGER       IGAS    !  I  Index of Gas for LUT file
      INTEGER       ISPC    !  I  Index of Spectral range using LUT 
      INTEGER       NL      !  I  No.of basis vectors in LUT
      INTEGER       NPT     !  I  No.of p,T points in LUT
      INTEGER       NWN     !  I  No.of Wavenumber points in LUT
      INTEGER       NHEAD   !  I  No.of header records that can be skipped
      LOGICAL       BINFIL  !  I  T=Binary file, F=ASCII file
      LOGICAL       IRRFIL  !  I  T=Irregular grid, F=full grid
      CHARACTER*3   TAB     !  I  LUT Tabulation function
      CHARACTER*(*) QAL     !  I  Optional qualifier appended to filename
      LOGICAL       FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  LOCASE ! Convert text string to lower case.
     &, RFMLOG ! Write text message to RFM log file.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'lflcom.inc' ! Look-Up Table file data
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
C
C LOCAL VARIABLES
      INTEGER   ILFL   ! Counter for LUT files
      INTEGER   IOS    ! Value of IOSTAT saved for error messages
      INTEGER   IPT    ! Pointer to '(' at start of qualifier string
      INTEGER   JPT    ! Pointer to ')' at end of qualifier string
      INTEGER   KPT    ! Pointer to ':' within qualifier string
      INTEGER   NIGTOT ! Total No. of Tab file Irreg.Grid pts reqd in Spc.range
      INTEGER   NSV    ! Integer read from qualifier string
      INTEGER   NLUTOT ! Total No. LUTs required for Spc.range
      INTEGER   NXTOT  ! Total No. of Cmp.LUT p,T points required in Spc.range
      INTEGER   NIVTOT ! Total No. of Wno points required in Spc.range
      INTEGER   NPTLFL(MAXLFL)  ! No. of p,T pts reqd for each file
         SAVE   NPTLFL
      INTEGER   NWNLFL(MAXLFL)  ! No. of Wno pts reqd for each file
         SAVE   NWNLFL
      CHARACTER*80 LOGMSG ! Message sent to Log file
      LOGICAL   CMPLFL(MAXLFL)  ! T=Compressed LUT, F=Uncompressed LUT
         SAVE   CMPLFL
C
C DATA STATEMENTS
           DATA NWNLFL / MAXLFL * 0 /
           DATA NPTLFL / MAXLFL * 0 /
           DATA CMPLFL / MAXLFL * .TRUE. /
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      FAIL = .TRUE.
C
C Check total number of LUT files within array dimension MAXLFL
C
      IF ( NLFL .EQ. MAXLFL ) THEN
        WRITE ( ERRMSG, '(A,I11)' )
     &    'F-LUTCHK: Too many LUTs for one RFM run. '//
     &    'Max = MAXLFL in RFMSIZ.INC=', MAXLFL
        RETURN
      ELSE
        NLFL = NLFL + 1
      END IF
C
C Check if compressed or uncompressed file 
C NB: this needs to be determined here and stored in this module since this 
C can only be determined when NL value is supplied for each file. Normally
C NL is also saved as NSVLFL but this could be set to 0 by using qualifiers
C so NSVLFL=0 is not sufficient to determine whether or not a LUT is compressed.
C
      CMPLFL(NLFL) = NL .GT. 0
C
C Check that Tabulation function is handled by RFMLUT
C
      CALL LOCASE ( TAB, TAB )
      IF ( TAB .EQ. 'lin' .OR. TAB .EQ. '4rt' .OR. TAB .EQ. 'log' ) THEN
        TABLFL(NLFL) = TAB
      ELSE
        ERRMSG = 
     &    'F-LUTCHK: Unrecognised LUT tabulation function, TAB='//TAB
        RETURN
      END IF
C
C Check any qualifier information can be successfully decoded
C
      IPT = INDEX ( QAL, '(' )           
      IF ( IPT .EQ. 0 ) THEN             ! No qualifier data, just set NSVLFL=NL
        NSVLFL(NLFL) = NL
        NDPLFL(NLFL) = 1
        NDTLFL(NLFL) = 1
      ELSE                                       ! some qualifier data appended
        JPT = INDEX ( QAL, ')' ) 
        IF ( JPT .EQ. 0 .OR. QAL(JPT+1:) .NE. ' ' )  THEN
          ERRMSG = 'F-LUTCHK: Failed to find '')'' at end of '//
     &             'qualifier attached to LUT filename'
          RETURN
        END IF
        IF ( CMPLFL(NLFL) ) THEN                 ! Compressed LUT file
          READ ( QAL(IPT+1:JPT-1), *, IOSTAT=IOS ) NSV 
          IF ( IOS .NE. 0 ) THEN
            WRITE ( ERRMSG, '(A,I11)' )
     &        'F-LUTCHK: Failed to read contents of LUT Qualifier.'//
     &        ' IOSTAT=', IOS
            RETURN
          END IF
C
C Successfully read contents of qualifier, so save and send message to logfile
          IF ( NSV .GT. NL ) THEN
            WRITE ( LOGMSG, '(A,I11,A,I11,A)' )
     &        'W-LUTCHK: Including all ', NL, 
     &        ' basis vectors in file (<No.Rqd=', NSV, ')'
            NSVLFL(NLFL) = NL
          ELSE IF ( NSV .GE. 0 ) THEN
            WRITE ( LOGMSG, '(A,I11,A,I11,A)' )
     &        'I-LUTCHK: Including ', NSV, ' out of ', NL, 
     &        ' basis vectors from LUT file'
            NSVLFL(NLFL) = NSV
          ELSE IF ( NSV .GE. -NL ) THEN
            WRITE ( LOGMSG, '(A,I11,A,I11,A)' )
     &        'I-LUTCHK: Excluding ', -NSV, ' out of ', NL, 
     &        ' basis vectors from LUT file'
            NSVLFL(NLFL) = NL + NSV
          ELSE 
            WRITE ( LOGMSG, '(A,I11,A,I11,A)' )
     &        'W-LUTCHK: Excluding all ', NL,  
     &        ' basis vectors in file (<No.Rqd=', -NSV, ')'
            NSVLFL(NLFL) = 0
          END IF
          CALL RFMLOG ( LOGMSG, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
C If the number of singular vectors is zero, it is desirable to force this to
C give zero absorption, in which case a 'log' tabulation needs to be alterered
C (otherwise k=exp(0.0) would give unity absorption)
          IF ( NSVLFL(NLFL) .EQ. 0 .AND. TABLFL(NLFL) .EQ. 'log' ) THEN
            LOGMSG = 'W-LUTCHK: No basis vectors rqd so assuming '//
     &        '''lin'' k-tabulation for 0 absorption'
            TABLFL(NLFL) = 'lin'
            CALL RFMLOG ( LOGMSG, FAIL, ERRMSG )
            IF ( FAIL ) RETURN
          END IF
C
C Read qualifier associated with a TAB input file
        ELSE
          KPT = INDEX ( QAL, ':' )
          IF ( KPT .EQ. 0 ) THEN
            ERRMSG = 'F-LUTCHK: No '':'' character found '//
     &               'in qualifier attached to uncompressed file.'
            RETURN
          END IF
          QAL(KPT:KPT) = ' '
          READ ( QAL(IPT+1:JPT-1), *, IOSTAT=IOS ) 
     &      NDPLFL(NLFL), NDTLFL(NLFL)
          IF ( IOS .NE. 0 ) THEN
            WRITE ( ERRMSG, '(A,I11)' )
     &        'F-LUTCHK: Failed to read 2 integers in qualifier.'//
     &        ' IOSTAT=', IOS
            RETURN
          END IF
C
C Successfully read contents of qualifier, so save and send message to logfile
          WRITE ( LOGMSG, '(A,I11,A,I11)' )
     &      'I-LUTCHK: Reducing -lnp,tem axes by factors ', 
     &      NDPLFL(NLFL), ',', NDTLFL(NLFL)
          CALL RFMLOG ( LOGMSG, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        END IF 
      END IF
C
C Check no.of basis vectors within array dimension MAXLUL
C
      FAIL = .TRUE.
      IF ( NSVLFL(NLFL) .GT. MAXLUL ) THEN 
        WRITE ( ERRMSG, '(A,I11,A,I11)' )
     &    'F-LUTCHK: No. basis vec.reqd. NL=', NL, 
     &    '> MAXLUL in RFMSIZ.INC,=', MAXLUL
        RETURN
      END IF
C
      NLUTOT = 0
      NXTOT = 0
      NIGTOT = 0
      NIVTOT = 0
C
C Find which previously included LUTs apply to the same spectral range
C
      DO ILFL = 1, NLFL-1
        IF ( SPCLFL(ILFL) .EQ. ISPC ) THEN
          NLUTOT = NLUTOT + 1
          NXTOT  = NXTOT  + NPTLFL(ILFL)                        ! All LUTs
          IF ( CMPLFL(ILFL) ) THEN
            NIVTOT = NIVTOT + NWNLFL(ILFL)    ! Compressed LUT
          ELSE IF ( IRRLFL(ILFL) ) THEN       ! Irr.Grid and uncompressed LUT
            NIGTOT = NIGTOT + NWNLFL(ILFL)
          END IF
C
C Check that another LUT file doesn't contain the same molecule/Spc.range
C
          IF ( GASLFL(ILFL) .EQ. IGAS ) THEN 
            WRITE ( ERRMSG, '(A,A,A,A)' )
     &        'F-LUTCHK: Two Look-Up Tables supplied for gas=', 
     &        CODGAS(IGAS), ' for Spc. Range=', LABSPC(ISPC)
            RETURN
          END IF
        END IF
      END DO
      GASLFL(NLFL) = IGAS
C
C Check total number of LUTs required for Spc.range within array dim MAXLUT.
C (NLUTOT is number excluding the current file, so add 1)
C
      IF ( NLUTOT + 1 .GT. MAXLUT ) THEN
        WRITE ( ERRMSG, '(A,A,A,I11)' )
     &    'F-LUTCHK: Too many LUTs for Spc.Rng=', LABSPC(ISPC),
     &    ' >MAXLUT in RFMSIZ.INC=', MAXLUT
        RETURN
      ELSE
        SPCLFL(NLFL) = ISPC
      END IF
C
C For compressed files, both (p,T) and Waveno.dimensions are stored internally
C but for uncompressed files (eg generated with the TAB option) only (p,T) has 
C to be stored 
C
C Check total number of p,T points for Spc.range within array dimensions
      IF ( NPT + NXTOT .GT. MAXLUX ) THEN
        WRITE ( ERRMSG, '(A,A,A,I11)' )
     &    'F-LUTCHK: Total p,T pts for Spc.Rng=', LABSPC(ISPC),
     &    ' >MAXLUX in RFMSIZ.INC=', MAXLUX
        RETURN
      ELSE IF ( CMPLFL(NLFL) .AND. NWN + NIVTOT .GT. MAXLUV ) THEN
C Check total number of Wno points for Spc.range within array dimension MAXLUV
        WRITE ( ERRMSG, '(A,A,A,I11)' )
     &    'F-LUTCHK: Total Wno pts for Spc.Rng=', LABSPC(ISPC),
     &    ' >MAXLUV in RFMSIZ.INC=', MAXLUV
        RETURN
      ELSE IF ( .NOT. CMPLFL(NLFL) .AND. NWN+NIGTOT .GT. MAXLUG ) THEN
C Check total number of Wno points for Spc.range within array dimension MAXLUG
        WRITE ( ERRMSG, '(A,A,A,I11)' )
     &    'F-LUTCHK: Total Wno pts for Spc.Rng=', LABSPC(ISPC),
     &    ' >MAXLUG in RFMSIZ.INC=', MAXLUG
        RETURN
      END IF
C
      NPTLFL(NLFL) = NPT
      NWNLFL(NLFL) = NWN
      LUNLFL(NLFL) = LUN
      NHDLFL(NLFL) = NHEAD
      BINLFL(NLFL) = BINFIL
      IRRLFL(NLFL) = IRRFIL
C
      FAIL = .FALSE.
C
      END
