      SUBROUTINE CHKGAS ( INPSTR, GAS, IMOL, ISO, IGQ, FAIL, ERRMSG )
C
C VERSION
C     14-AUG-13  AD  Simplified using MOLIDX subroutine
C     23-JUL-09  AD  Correction: only check XSCMOL for line molecules
C                    Correction: calculation of IEND
C     19-JUN-08  AD  Allow for 'x' appended to gas = use .xsc data instead
C     14-APR-08  AD  Add CH3OH as Molec#39 (=methanol)
C                        C2F6 (CFC-116) as Molec#82 (=hexafluoromethane)
C                        SF5CF3 as Molec#83 
C                        PAN as Molec#84 (=peroxyacetyl nitrate)
C                        CH3CN  as Molec#85 (=methyl cyanide or acetonitrile)
C     08-JAN-08  AD  Original.
C
C DESCRIPTION
C     Check for valid molecule name and get Hitran/RFM index 
C     RFM general purpose subroutine.
C     This routine performs the following
C     (1) Copies INPSTR as far as any '(' character to GAS and converts GAS to 
C         lower case
C     (2) GAS is first compared to standard molecule names
C         If gas = 'clono2q', return as 'clono2' but with HITRAN IMOL=35 value
C         similarly for other sf6, cf4 etc
C     (3) If (2) not found, checks gas against some common abbreviations
C         for isotopomers. If found sets GAS=standard name and ISO=N where N is
C         the isotope#
C     (4) If (2) and (3) not found, checks if GAS can be interpreted as an
C         integer value representing the HITRAN/RFM index IMOL directly
C     If the string is a recognised molecule, this routine returns IMOL>0,
C     and if not recognised, then IMOL=0.
C     (5) If there is a character '(' in the INPSTR, the contents of the first
C         pair of brackets is returned as ISO and, if there is a second pair,
C         their contents are returned as IGQ.
C
      IMPLICIT NONE
C
C ARGUMENTS
      CHARACTER*(*) INPSTR  !  I  String to be tested 
      CHARACTER*7   GAS     !  O  Standard RFM name of gas
      INTEGER       IMOL    !  O  HITRAN/RFM Index of gas, or 0 if unrecognised
      INTEGER       ISO     !  O  Isotope#n, or else 0 if not isotopomer
      INTEGER       IGQ     !  O  Global Quantum index, or else 0 
      LOGICAL       FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG  !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  LOCASE ! Convert text string to lower case
     &, MOLIDX ! Give molecule name for HITRAN/RFM index, or vice-versa
     &, RFMLOG ! Write text message to log file.
C
C LOCAL VARIABLES
      INTEGER      ILEN           ! Length of non-blank part of INPSTR 
      INTEGER      IEND           ! Position of ')' in input GAS
      INTEGER      ISTA           ! Position of '(' in input GAS
      INTEGER      JMOL           ! Molecule index read from GAS string
      INTEGER      IOS            ! Value of IOSTAT for READing string
      CHARACTER*7  MOLEC          ! Molecule returned for given JMOL
      CHARACTER*20 ERRSTR         ! Part of INPSTR used in error messages
      CHARACTER*80 MESSGE         ! Messages for Log file
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      IMOL = 0
      ISO = 0
      IGQ = 0
      GAS = ' '
      FAIL = .FALSE.
C
      ILEN = INDEX ( INPSTR, ' ' ) - 1       ! Position of last non-blank char.
      IF ( ILEN .EQ. -1 ) ILEN = LEN ( INPSTR )
      ERRSTR = INPSTR(1:MIN(ILEN,20))        ! Limit to C*20 for error messages
      IEND = ILEN
C
C (1) If input string contains a bracket, just extract text before bracket
      ISTA = INDEX ( INPSTR, '(' )
      IF ( ISTA .GT. 0 ) IEND = ISTA - 1
C
C Maximum length of a recognised molecule name is C*7
      IF ( IEND .EQ. 0 .OR. IEND .GT. 7 ) RETURN
      GAS = INPSTR(1:IEND)

C Convert GAS to lower case
      CALL LOCASE ( GAS, GAS )
C
C (2) Assume that GAS is a molecule name:
      IMOL = 0
      CALL MOLIDX ( IMOL, GAS, FAIL, ERRMSG )
      IF ( FAIL ) RETURN  ! only FAIL if ambiguity over line or x/s data
      IF ( IMOL .GT. 0 ) GOTO 100
C
C (3) Assume that GAS is isotopologue abbreviation
      IF ( GAS .EQ. 'hdo' ) THEN
        IMOL = 1
        ISO = 4
        MESSGE = 'I-CHKGAS: Converting ''hdo'' to h2o(4)'
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        GOTO 100 
      ELSE IF ( GAS .EQ. 'ch3d' ) THEN
        IMOL = 6
        ISO = 3
        MESSGE = 'I-CHKGAS: Converting ''ch3d'' to ch4(3)'
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        GOTO 100
      END IF
C 
C (4) See if GAS can be read as an integer index within the valid range
      JMOL = 0
      IOS = 1
      READ ( GAS, *, IOSTAT=IOS ) JMOL
      IF ( IOS .EQ. 0 ) THEN
        CALL MOLIDX ( JMOL, MOLEC, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        IF ( MOLEC .NE. ' ' ) THEN
          GAS = MOLEC
          WRITE ( MESSGE, '(A,I3,A,A)' )
     &        'I-CHKGAS: Translating HITRAN/RFM Index=', JMOL,
     &        ' as gas=', GAS
          CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
          IMOL = JMOL
          GOTO 100
        END IF
      END IF
C
C If program gets to this point, unable to identify molecule
      FAIL = .FALSE.
      IMOL = 0
      RETURN
C
 100  CONTINUE
C
      IF ( IEND .EQ. ILEN ) RETURN        ! Normal exit for just molecule ID
C
C (5) This second part deals with any appended '(..)(..)' strings to the 
C molecule name, interpreting the contents of the first bracket as an isotope#
C and the second bracket as a Global Quantum level index. These either come
C from the *ATM or *JAC sections of the driver table, in the *GAS section the
C qualifiers are removed before this is called.
C
C If isotopomeric abbreviation used, isotope is already defined so that anything
C in brackets is interpreted as the Global Quantum index
      IF ( ISO .NE. 0 ) GOTO 200
C
      IEND = INDEX ( INPSTR, ')' )
      IF ( IEND .LT. ISTA+2 ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-CHKGAS: Unable to read isotopic number: '//ERRSTR
        RETURN
      END IF
      READ ( INPSTR(ISTA+1:IEND-1), *, IOSTAT=IOS ) ISO ! Read isotope number
      IF ( IOS .NE. 0 ) THEN
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,I11)' ) 
     &    'F-CHKGAS: Bad isotope# in string: '//ERRSTR//'. IOSTAT=', 
     &    IOS
        RETURN
      END IF
      IF ( IEND .EQ. ILEN ) RETURN       ! Normal exit for molecule + isotope
C
C GQI should follow on immediately from isotope# as '(..)'
 200  CONTINUE
      ISTA = IEND + 1
      IF ( INPSTR(ISTA:ISTA) .NE. '(' ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-CHKGAS: Unable to parse absorber name: '//ERRSTR
        RETURN
      END IF
      IEND = IEND + INDEX ( INPSTR(ISTA:), ')' ) 
      IF ( IEND .LT. ISTA+2 ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-CHKGAS: Unable to read GQ index from: '//ERRSTR
        RETURN
      END IF
      READ ( INPSTR(ISTA+1:IEND-1), *, IOSTAT=IOS ) IGQ ! Read GQ Index
      IF ( IOS .NE. 0 ) THEN
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,I11)' ) 
     &    'F-CHKGAS: Bad Global Quantum# in string: '//ERRSTR//
     &    '. IOSTAT=', IOS
        RETURN
      END IF
C
      IF ( IEND .EQ. ILEN ) RETURN   ! Normal exit for molecule + isotope + GQI
C
      FAIL = .TRUE.
      ERRMSG = 'F-CHKGAS: Unable to parse string: '//ERRSTR
C
      END
