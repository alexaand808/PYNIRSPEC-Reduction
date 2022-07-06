      SUBROUTINE LUTFIL ( LUNLUT, NAMLUT, FAIL, ERRMSG )
C
C VERSION
C     04-APR-06  AD  Change 'X' to '1X' in format 
C     25-OCT-04  AD  Bug Fix: Check isotopic LUTs actually required
C     26-MAR-04  AD  Allow for isotopic LUTs
C     21-JUN-01  AD  Assume binary if ".bin" or ".BIN" in filename 
C     05-DEC-00  AD  Comment out READONLY
C     08-JUN-00  AD  Allow for Irregular grids
C     09-AUG-99  AD  Make checks on V1,V2 v. WNLSPC,WNUSPC explicitly D.P. 
C                    Cross check MWCODE v. LABSPC and warn if different
C                    Change names of variables P1,T1 to PDUMMY,TDUMMY
C     28-JUN-99  AD  Additional binary file check on MWCODE_ID_TAB record
C     05-JUN-99  AD  Allow for ESA format LUT files as well. 
C                    Allow for binary LUT file
C                    Write first record(s) to LOG file.
C                    Count number of header records in file
C     24-MAY-99  AD  Check for 0 values of DP,DV,DT
C     23-APR-99  AD  Change MWCODE from C*6 to C*8
C     20-JAN-98  AD  Check for Qualifiers after LUT filename
C     18-JUL-97  AD  Read TAB from LUT file and pass to LUTCHK
C     03-JUL-97  AD  Original.
C
C DESCRIPTION
C     Open Look-Up Table file and check contents
C     Called by INPLUT for each file specified in *LUT section of driver table.
C     Keeps LUN open for all usable LUT files.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER       LUNLUT  ! I/O LUN for LUT File/next free LUN
      CHARACTER*(*) NAMLUT  ! I/O Name of LUT file (+optional qualifier)
      LOGICAL       FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  LUTCHK ! Check Look-Up Table data.
     &, LUTGRD ! Read coded irregular grid data and turn into integer indices.
     &, RFMLOG ! Write text message to RFM log file.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
C
C LOCAL CONSTANTS
      DOUBLE PRECISION WNOTOL    ! Waveno. tolerance as fraction of LUT spacing
        PARAMETER ( WNOTOL = 0.01D0 ) 
C
C LOCAL VARIABLES
      INTEGER       ID     ! HITRAN Molecule ID contained in LUT file
      INTEGER       IDUMMY(1) ! Dummy array
      INTEGER       IGAS   ! Index for gas
      INTEGER       IREC   ! Record counter
      INTEGER       IOS    ! Saved value of IOSTAT for error messages
      INTEGER       IPT    ! Pointer to '(' marking start of qualifier data
      INTEGER       ISO    ! Isotope# (or -1 if generic LUT)
      INTEGER       ISPC   ! Spectral range counter (1:NSPC)
      INTEGER       LPT    ! Length of original NAMLUT string
      INTEGER       NWN    ! No.of wavenumber points in LUT/TAB file
      INTEGER       NHEAD  ! No.of header records that can be skipped
      INTEGER       NL     ! No.of basis vectors file
      INTEGER       NP     ! No.of pressure points 
      INTEGER       NPT    ! NP * NT
      INTEGER       NREC   ! Counter for number of records read
      INTEGER       NT     ! No.of temperature points
      INTEGER       NV     ! No.of wavenumber points          
      LOGICAL       ANYQAL ! T=qualifiers appended to end of file name
      LOGICAL       BINFIL ! T=Binary file, F=ASCII file
      LOGICAL       IRRFIL ! T=Irregular grid, F=Full grid
      LOGICAL       REJECT ! T=ignore file - no useful data
      REAL          DP     ! Increment in -lnp
      REAL          DT     ! Temperature increment [K]
      REAL          PDUMMY ! Lower values of -lnp [p in mb] (not used)
      REAL          TDUMMY ! Lower values of T [K] (not used)
      CHARACTER*132 MESSGE ! Text message sent to LOG file
      CHARACTER*8   MWCODE ! MW code contained in LUT file
      CHARACTER*14  QAL    ! Qualifier string (if any) appended to filename
      CHARACTER*80  RECORD ! Record read from LUT file
      CHARACTER*3   TAB    ! Code for tabulated function
      DOUBLE PRECISION V1,V2  ! Lower,Upper wavenumbers of LUT [/cm]
      DOUBLE PRECISION DV     ! Wavenumber increment [/cm]
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Separate any qualifier string from the filename
C
      ANYQAL = INDEX ( NAMLUT, '(' ) .GT. 1
      IF ( ANYQAL ) THEN
        IPT = INDEX ( NAMLUT, '(' )
        LPT = LEN ( NAMLUT )
        QAL = NAMLUT(IPT:LPT)
        NAMLUT = NAMLUT(1:IPT-1)
      ELSE
        QAL = ' '
      END IF
C
C Send message to LOG file saying which file is about to be opened
C
      MESSGE = 'I-LUTFIL: Opening LUT File: '//NAMLUT
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
C Open File. Even if binary file, open as ASCII so that headers can be read
C and written out to LOG file without worrying about length of record
C
      OPEN ( UNIT=LUNLUT, FILE=NAMLUT, STATUS='OLD', 
     &       IOSTAT=IOS, ERR=900 )
      REJECT = .TRUE.
C
C If filename contains ".bin" or ".BIN", assume binary file
      IF ( INDEX ( NAMLUT, '.bin' ) .NE. 0 .OR. 
     &     INDEX ( NAMLUT, '.BIN' ) .NE. 0      ) THEN
        MESSGE = 
     &  'I-LUTFIL: Filename contains .bin or .BIN so assume binary'
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        IOS = 1         ! Simulates an I/O error, required for switch to binary
        GOTO 100
      END IF
C
C Test file by reading down to end of dimensions record: assume that if any 
C errors occur this could be a binary file.
C
      READ ( LUNLUT, '(A)', ERR=100, IOSTAT=IOS ) RECORD
      READ ( LUNLUT, '(A)', ERR=100, IOSTAT=IOS ) RECORD
      NREC = 2
      DO WHILE ( INDEX ( RECORD(1:3), '#' ) .NE. 0 .OR.
     &           INDEX ( RECORD(1:3), '!' ) .NE. 0      )
        READ ( LUNLUT, '(A15)', ERR=900, IOSTAT=IOS ) RECORD(1:15)
        NREC = NREC + 1
      END DO
      NHEAD = NREC - 1 
      READ ( RECORD, '(A8,1X,I2,1X,A3)', ERR=100, IOSTAT=IOS ) 
     &  MWCODE, ID, TAB       ! should be OK if isotope ID added as well
      IF ( RECORD(9:9) .NE. ' ' .OR. 
     &     ( RECORD(12:12) .NE. ' ' .AND. 
     &       RECORD(12:12) .NE. '.'       ) ) THEN
        IOS = 1
        GOTO 100
      END IF
      READ ( LUNLUT, *, ERR=100, IOSTAT=IOS ) 
     &  NL, NV, V1, DV, NP, PDUMMY, DP, NT, TDUMMY, DT
C
 100  CONTINUE
      IF ( IOS .NE. 0 ) THEN    ! Error reading, assume due to unformatted data
        NHEAD = 3               ! Assume all binary files have 3 header records
        BINFIL = .TRUE.
        MESSGE = 'I-LUTFIL: Reading as a binary file'
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        CLOSE ( LUNLUT )
        OPEN ( UNIT=LUNLUT, FILE=NAMLUT, STATUS='OLD', 
     &         FORM='UNFORMATTED', IOSTAT=IOS, ERR=900 )
        DO IREC = 1, NHEAD
          READ ( LUNLUT ) RECORD          ! Skip header records
          IF ( IREC .EQ. 1 ) THEN
            CALL RFMLOG ( RECORD, FAIL, ERRMSG )
            IF ( FAIL ) RETURN
          END IF
        END DO
        READ ( LUNLUT ) RECORD(1:15)      ! Re-read MWCODE, ID, TAB record
        IF ( RECORD(12:12) .EQ. '.' ) THEN     ! Isotope code included
          BACKSPACE ( LUNLUT ) 
          READ ( LUNLUT ) RECORD(1:17)    ! Re-read MWCODE, ID+ISO, TAB record
        END IF          
      ELSE
        BINFIL = .FALSE.
        REWIND ( LUNLUT )   ! Reposition to start of file
C
C Write first record to log file 
        READ ( LUNLUT, '(A)', ERR=900, IOSTAT=IOS ) RECORD
        NREC = 1
        CALL RFMLOG ( RECORD, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        IF ( INDEX ( RECORD(1:3), '!' ) .EQ. 0 ) THEN              ! ESA format
          READ ( LUNLUT, '(A)', ERR=900, IOSTAT=IOS ) RECORD
          NREC = NREC + 1
          CALL RFMLOG ( RECORD, FAIL, ERRMSG )  ! Write 2nd rec to log file too
          IF ( FAIL ) RETURN
        END IF
        DO IREC = NREC, NHEAD       ! Load MW code etc into RECORD
          READ ( LUNLUT, '(A)', ERR=900, IOSTAT=IOS ) RECORD
        END DO
      END IF
C
C Read MW code, HITRAN ID and TABulation function from first record
C
      IF ( RECORD(12:12) .EQ. '.' ) THEN
        READ ( RECORD, '(A8,1X,I2,1X,I1,1X,A3)', ERR=900, IOSTAT=IOS ) 
     &    MWCODE, ID, ISO, TAB
      ELSE
        READ ( RECORD, '(A8,1X,I2,1X,A3)', ERR=900, IOSTAT=IOS ) 
     &    MWCODE, ID, TAB
        ISO = -1
      END IF
C
C Check HITRAN ID is within valid range
C
      IF ( ID .LT. 1 .OR. ID .GT. MAXMOL ) THEN
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,I11)' )
     &    'F-LUTFIL: File contains unrecognised ID value, =', ID
        RETURN
C
C Check if this species is required
C
      ELSE IF ( IGSMOL(ID) .NE. 0 ) THEN       ! this is a required species
        IGAS = IGSMOL(ID)
C
C Check for discrepancy between isotopic/generic profiles and LUTs
        IF ( .NOT. ISOMOL(ID) ) THEN           ! generic profile
          IF ( ISO .NE. -1 ) THEN              ! isotopic LUT
            WRITE ( MESSGE, '(A)' ) 'W-LUTFIL: generic profile '//
     &        'specified in *ATM sect but isotopic LUT supplied'
            CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
            IF ( FAIL ) RETURN
            REJECT = .TRUE.
            GOTO 200
          END IF
        ELSE                                   ! isotopic profile
          IF ( ISO .EQ. -1 ) THEN              ! generic LUT
            WRITE ( MESSGE, '(A)' ) 'W-LUTFIL: isotopic profile '//
     &        'specified in *ATM sect but generic LUT supplied'
            CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
            IF ( FAIL ) RETURN
            REJECT = .TRUE.
            GOTO 200
          ELSE IF ( ISO .GT. NISGAS(IGAS) ) THEN  ! NB: ISO was read '(I1)'
            WRITE ( ERRMSG, '(A,I1,A,A)' ) 
     &        'F-LUTFIL: isotope#', ISO, 
     &        ' is invalid for gas=', CODGAS(IGAS)
            FAIL = .TRUE.
            RETURN
          ELSE IF ( ISO .GT. 0 .AND.        ! Isotopic LUT but no isotopic prfl
     &              ISOGAS(ISO,IGAS) .EQ. ISOGAS(0,IGAS) ) THEN
            WRITE ( MESSGE, '(A,A,A,I1)' ) 
     &        'W-LUTFIL: No isotopic profile specified for gas=', 
     &        CODGAS(IGAS), ' isotope#', ISO
            CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
            IF ( FAIL ) RETURN
            REJECT = .TRUE.
            GOTO 200
          ELSE
            IGAS = ISOGAS(ISO,IGAS)       ! Reset IGAS for this isotope
          END IF
        END IF
C
        IF ( BINFIL ) THEN
          READ ( LUNLUT, ERR=900, IOSTAT=IOS ) 
     &      NL, NV, V1, DV, NP, PDUMMY, DP, NT, TDUMMY, DT
        ELSE
          READ ( LUNLUT, *, ERR=900, IOSTAT=IOS ) 
     &      NL, NV, V1, DV, NP, PDUMMY, DP, NT, TDUMMY, DT
        END IF
C
C Check axis parameters are legal
        FAIL = .TRUE.
        IF ( NV .EQ. 0 ) THEN
          ERRMSG = 'F-LUTFIL: NV must be .NE. 0'
        ELSE IF ( NP .LE. 0 ) THEN  
          ERRMSG = 'F-LUTFIL: NP must be .GE. 1'
        ELSE IF ( NT .LE. 0 ) THEN  
          ERRMSG = 'F-LUTFIL: NT must be .GE. 1'
        ELSE IF ( DV .EQ. 0.0 .AND. NV .NE. 1 ) THEN
          ERRMSG = 'F-LUTFIL: DV must be non-zero unless NV=1'
        ELSE IF ( DP .EQ. 0.0 .AND. NP .NE. 1 ) THEN
          ERRMSG = 'F-LUTFIL: DP must be non-zero unless NP=1'
        ELSE IF ( DT .EQ. 0.0 .AND. NT .NE. 1 ) THEN
          ERRMSG = 'F-LUTFIL: DT must be non-zero unless NT=1'
        ELSE 
          FAIL = .FALSE.
        END IF
        IF ( FAIL ) RETURN
C
        IRRFIL = NV .LT. 0
        NV = ABS ( NV ) 
        V2 = V1 + ( NV - 1 ) * DV
        NPT = NP * NT
C
C Check if LUT completely spans a selected output spectral range 
C NB: same LUT may be usable for more than one spectral range
C
        DO ISPC = 1, NSPC
          NWN = NV
          IF ( V1 .LE. WNLSPC(ISPC)+WNOTOL*DV .AND. 
     &         V2 .GE. WNUSPC(ISPC)-WNOTOL*DV   ) THEN   ! OK, so use this LUT
            IF ( IRRFIL ) THEN
              CALL LUTGRD ( LUNLUT, BINFIL, 1, NWN, 
     &                      IDUMMY, FAIL, ERRMSG )
              IF ( FAIL ) RETURN
            END IF
            CALL LUTCHK ( LUNLUT, IGAS, ISPC, NL, NPT, NWN, NHEAD, 
     &                    BINFIL, IRRFIL, TAB, QAL, FAIL, ERRMSG )
            IF ( FAIL ) RETURN
            REJECT = .FALSE.
            IF ( MWCODE .NE. LABSPC(ISPC) .AND. 
     &           LABSPC(ISPC) .NE. ' '          ) THEN
              MESSGE = 'W-LUTFIL: LUT label='//MWCODE//
     &          ' .NE. spc.range label='//LABSPC(ISPC)
              CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
              IF ( FAIL ) RETURN
            END IF
          END IF
        END DO
      END IF
C
 200  CONTINUE
      IF ( REJECT ) THEN
        CLOSE ( LUNLUT, IOSTAT=IOS, ERR=900 )
        MESSGE = 'W-LUTFIL: File ignored - Look-Up Table not applicable'
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      ELSE
        LUNLUT = LUNLUT + 1
      END IF
C
      FAIL = .FALSE.
      RETURN
C
  900 WRITE ( ERRMSG, '(A,I11)') 
     &  'F-LUTFIL: I/O failure on LUT file. IOSTAT=', IOS
      FAIL = .TRUE.
C
      END
