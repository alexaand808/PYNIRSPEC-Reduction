      PROGRAM HITBIN
C
C VERSION (update VERSID)
C     30JUN09 AD 1.20 HITRAN2008: allow for SBROAD in F5.4 or F5.3 (2008)
C     09MAR02 AD 1.10 DUPCHK: if two lines have equal status change from 
C                       termination condition to warning.
C                       Remove redundant MAXNEW parameter.
C     09SEP01 AD 1.00 Original. Based on HITLIN.
C
C DESCRIPTION
C     Convert HITRAN sequential ASCII file to direct access binary file
C     For single file this is functionally the same as HITLIN.
C     Differs from HITLIN in that HITLIN allows merging of several ASCII files
C     simultaneously into single binary file, whereas this program allows 
C     sequential merging of ASCII files into existing binary.
C     This program also removes duplicate records rather than setting LSTAT -ve
C     Note: the record length of the unformatted output file depends
C           on the no of bytes per word of the machine being used. 
C           System specific RECL is used in the output binary file 
C           open statement
C
C     Variables for each transition record (100 characters in input records)
C              Input Output Description
C      LSTAT     -     I    Priority of transition information.
C      IGAS     I2     I    Molecule number.
C      ISO      I1     I    Isotope number (=1 most abundant, 2=second etc)
C      WNUM    F12.6   D    Line frequency [cm-1].
C      STREN   E10.3   R    Line strength  I=[cm-1./(molec.cm-2)] @ 296K.
C                                          O=[cm-1./(kg.moles.cm-2)] @296K
C                           (Scale input by Avogadro No.to avoid underflows.)
C      TPROB   E10.3   R    Transition probability [Debyes2].
C      ABROAD   F5.4   R    Air-broadened halfwidth  (HWHM) [cm-1/atm] @ 296K.
C      SBROAD   F5.4   R    Self-broadened halfwidth (HWHM) [cm-1/atm] @ 296K.
C               F5.3   R    HITRAN v2008 onwards: changed format
C      ELS     F10.4   R    Lower-state energy [cm-1].
C      ABCOEF   F4.2   R    Coefficient of temperature dependance of ABROAD
C      TSP      F8.6   R    Transition shift due to pressure 
C                           (F8.6 not F8.5 as stated in HITRAN ref., although
C                            some files have been written using F8.5 following
C                            HITRAN ref - could lead to confusion!)
C      IUSGQ    I3     I     Upper state global (=Vib) quanta index 
C      ILSGQ    I3     I     Lower state global quanta index.
C      USLQ     A9     A9    Upper state local (=Rot) quanta. 
C      BSLQ     A9     A9    Lower state local quanta. 
C      AI      3I1     A3    Accuracy indices for freq, strength, halfwidth 
C      REF     3I2     A6    Indices for lookup of references for freq,strength
C                            halfwidth, not yet used (AI+REF combined into A9)
C      IFWDPT    -     I     Forward pointer on data line.
C
C     LSTAT values:
C        -7 forward pointer block
C        -2 file termination record
C         0 file header record
C         4 other header records
C        >9 line transition
C
      IMPLICIT NONE
C
C LOCAL CONSTANTS
      INTEGER LUNASC              ! LUN for ASCII (input) file
        PARAMETER ( LUNASC = 1 )
      INTEGER LUNNEW              ! LUN for new binary (output) file
        PARAMETER ( LUNNEW = 2 )
      INTEGER LUNOLD              ! LUN for old binary (input) file
        PARAMETER ( LUNOLD = 3 )
      INTEGER MAXNFP              ! Max no.of forward pointers per record
        PARAMETER ( MAXNFP = 14 ) ! Don't change this
      INTEGER NFPREC              ! No. forward pointer records in a f.p.block
        PARAMETER ( NFPREC = 4 )  
      INTEGER MAXMOL              ! Max HITRAN molecule index allowed
        PARAMETER ( MAXMOL = NFPREC*MAXNFP ) ! Change NFPREC to increase MAXMOL
      INTEGER NRECFP              ! No of records between f.p. blocks
        PARAMETER ( NRECFP = 200 )
      REAL AVOG                       ! Avogadro's number (*1e-26)
        PARAMETER ( AVOG = 6.0221367 )
      DOUBLE PRECISION WNOCHK         ! Range to check for duplicate lines
        PARAMETER ( WNOCHK = 0.1D0 )  
      DOUBLE PRECISION WNOEND         ! Large Wno Flag for end of file
        PARAMETER ( WNOEND = 2.0D6 )  ! NB: larger than WNOBIG in USERIN
C
C LOCAL VARIABLES
      INTEGER IDUMMY          ! Dummy integer
      INTEGER IEXASC          ! STREN exponent read from ASCII file
      INTEGER IFPBIN          ! Forward pointer read from binary file
      INTEGER IFPREC          ! Record# within forward pointer block
      INTEGER ILSASC, ILSBIN  ! ILSGQ from ASCII, binary files
      INTEGER IMOL            ! Molecule counter
      INTEGER IOFF            ! Offset for molecule# in f.p. records
      INTEGER IRCMOL(MAXMOL)  ! Next record# containing transition for each mol
      INTEGER IREC,IRECO      ! Record# in new,old binary file
      INTEGER IREC1,IREC1O    ! First record# in new,old file after headers
      INTEGER IREC2,IREC2O    ! Last record# in new,old file (termination rec)
      INTEGER IRECFP          ! First record# in new file after last f.p.block
      INTEGER ISOASC,ISOBIN   ! ISO from ASCII, binary files
      INTEGER IUSASC,IUSBIN   ! IUSGQ from ASCII, binary files
      INTEGER JREC            ! Secondary record counter
      INTEGER LSTASC,LSTBIN   ! LSTAT associated with ASCII, binary files
      INTEGER MOLASC,MOLBIN   ! IGAS from ASCII, binary files
      INTEGER NEWSTT          ! Status# of ASCII data
      INTEGER RECLEN          ! Record length of binary files
      LOGICAL INSERT          ! T=insert current line, F=ignore current line
      LOGICAL MERGE           ! T=merge ASCII and binary files, F=ASCII only
      LOGICAL USEMOL(MAXMOL)  ! T=get molecule from ASCII file, F=ignore
      REAL    ABCASC, ABCBIN  ! ABCOEF from ASCII, binary files
      REAL    ABRASC, ABRBIN  ! ABROAD from ASCII, binary files
      REAL    ELSASC, ELSBIN  ! ELS from ASCII, binary files
      REAL    SBRASC, SBRBIN  ! SBROAD from ASCII, binary files
      REAL    STRASC, STRBIN  ! STREN from ASCII, binary files
      REAL    TPRASC, TPRBIN  ! TPROB from ASCII, binary files
      REAL    TSPASC, TSPBIN  ! TSP from ASCII, binary files
      DOUBLE PRECISION WNOASC,WNOBIN ! WNUM from ASCII, binary file
      DOUBLE PRECISION WNOUPP ! Last wavenumber written/read
      DOUBLE PRECISION WNO1   ! Lowest wavenumber to select from ASCII file
      DOUBLE PRECISION WNO2   ! Highest wavenumber to select from ASCII file
      CHARACTER*9 ARFASC,ARFBIN ! AI + REF from ASCII, binary files
      CHARACTER*48 HEAD48       ! User header for binary file
      CHARACTER*84 HEADR2       ! HITBIN header record for binary file
      CHARACTER*9 LSLASC,LSLBIN ! BSLQ from ASCII, binary files
      CHARACTER*5 SBRSTR        ! SBROAD read as character string from ASCII
      CHARACTER*9 USLASC,USLBIN ! USLQ from ASCII, binary files
      CHARACTER*4 VERSID        ! Version identifier
C
C DATA STATEMENTS
      DATA USEMOL / MAXMOL * .TRUE. /       ! Default: use all molecules
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      VERSID = '1.10'
      WRITE ( *, '(A)' ) 'R-HITBIN: Running HITBIN v'//VERSID
C Get user inputs and open files
      CALL USERIN ( LUNASC, LUNOLD, LUNNEW, MAXMOL, WNO1, WNO2, 
     &              USEMOL, MERGE, NEWSTT, RECLEN, HEAD48 )
C
      WRITE ( HEADR2, '(A,I3,A,A)' ) 'RECL=', RECLEN, 
     &  ' Converted to binary format by HITBIN v.', VERSID
C
C Copy header records to new file
      CALL NEWHDR ( LUNASC, LUNOLD, LUNNEW, HEADR2, MERGE, IREC )
C
      WRITE ( *, '(A)' ) 
     &  'I-HITBIN: Finding first record of ASCII file...'
C Load first record within required range from ASCII file
      READ ( LUNASC, '(I2,X,F12.6)' ) MOLASC, WNOASC
      DO WHILE ( WNOASC .LT. WNO1 .OR. .NOT. USEMOL(MOLASC) )
        READ ( LUNASC, '(I2,X,F12.6)', END=100 ) MOLASC, WNOASC
      END DO
 100  CONTINUE
      IF ( WNOASC .GE. WNO1 .AND. WNOASC .LE. WNO2 ) THEN      ! use this file
        BACKSPACE ( LUNASC )       
        READ ( LUNASC, 1000 ) 
     &    MOLASC, ISOASC, WNOASC, STRASC, IEXASC, TPRASC,
     &    ABRASC, SBRSTR, ELSASC, ABCASC, TSPASC, 
     &    IUSASC, ILSASC, USLASC, LSLASC, ARFASC
 1000   FORMAT( I2, I1, F12.6, F6.3, 1X, I3, 1PE10.3, 
     &          0PF5.4, A5, F10.4, F4.2, F8.6, 2I3, 2A9, A9 )
C To avoid overflow problems, STRASC and IEXASC read separately from ASCII
C file, and scaled by Avogadro's number
          STRASC = AVOG * STRASC * ( 10.0**( IEXASC + 26 ) )
          LSTASC = NEWSTT
C Self-Broadening could be F5.4 or F5.3 - this should cope with either
          READ ( SBRSTR, * ) SBRASC
      ELSE                                                   ! ignore file
        STOP 'F-HITBIN: ASCII file outside reqd wavenumber range'
      END IF
C
      IREC1 = IREC                        ! Save first data record
C
      WRITE ( *, '(A)' ) 'I-HITBIN: Writing new binary file...'
C Begin output data with an empty forward pointer block
      DO JREC = 1, NFPREC 
        WRITE ( LUNNEW, REC=IREC ) -7
        IREC = IREC + 1
      END DO
      IRECFP = IREC 
C
C If merging with old binary file, find first line record in old file
      IF ( MERGE ) THEN
        READ ( LUNOLD, REC=1 ) IDUMMY, IDUMMY, IREC1O, IREC2O 
        IRECO = IREC1O
        READ ( LUNOLD, REC=IRECO ) LSTBIN
        DO WHILE ( LSTBIN .LT. 10 )
          IRECO = IRECO + 1
          READ ( LUNOLD, REC=IRECO ) LSTBIN
        END DO
        READ ( LUNOLD, REC=IRECO ) 
     &    LSTBIN, MOLBIN, ISOBIN, WNOBIN, STRBIN, TPRBIN, 
     &    ABRBIN, SBRBIN, ELSBIN, ABCBIN, TSPBIN, IUSBIN, 
     &    ILSBIN, USLBIN, LSLBIN, ARFBIN, IFPBIN
      ELSE
        WNOBIN = WNOEND
      END IF
C
C Repeat for each record, taking lowest wavenumber 
      DO WHILE ( WNOBIN .LT. WNOEND .OR. WNOASC .LT. WNOEND ) 
C Check if forward pointer block required
         IF ( IREC - IRECFP .EQ. NRECFP .OR.           ! insert f.p. block
     &        IREC .EQ. IRECFP-NFPREC        ) THEN    ! reinsert f.p. block
          DO JREC = 1, NFPREC  
            WRITE ( LUNNEW, REC=IREC ) -7
            IREC = IREC + 1
          END DO
          IRECFP = IREC
        END IF
        INSERT = .TRUE.
C Copy record from ASCII file to new file
        IF ( WNOASC .LT. WNOBIN ) THEN
          IF ( WNOBIN - WNOASC .LE. WNOCHK ) 
     &      CALL DUPCHK ( LUNNEW, IREC, WNOCHK, LSTASC, MOLASC, 
     &        ISOASC, WNOASC, IUSASC, ILSASC, USLASC, LSLASC, INSERT )
          IF ( INSERT ) THEN
            WRITE ( LUNNEW, REC=IREC ) 
     &        LSTASC, MOLASC, ISOASC, WNOASC, STRASC, TPRASC,
     &        ABRASC, SBRASC, ELSASC, ABCASC, TSPASC, IUSASC, 
     &        ILSASC, USLASC, LSLASC, ARFASC, 0          ! 0 = forward pointer
            IREC = IREC + 1
            WNOUPP = WNOASC
          END IF
 200      CONTINUE
          WNOASC = WNOEND                        ! Flag for EOF
          READ ( LUNASC, 1000, END=300 ) 
     &      MOLASC, ISOASC, WNOASC, STRASC, IEXASC, TPRASC,
     &      ABRASC, SBRSTR, ELSASC, ABCASC, TSPASC, IUSASC,
     &      ILSASC, USLASC, LSLASC, ARFASC
          IF ( WNOASC .GT. WNO2 ) THEN
            WNOASC = WNOEND
            GOTO 300
          END IF
          IF ( .NOT. USEMOL(MOLASC) ) GOTO 200
C To avoid overflow problems, STRASC and IEXASC read separately from ASCII
C file, and scaled by Avogadro's number
          STRASC = AVOG * STRASC * ( 10.0**( IEXASC + 26 ) )
          READ ( SBRSTR, * ) SBRASC
          LSTASC = NEWSTT
 300      CONTINUE
C Copy record from old binary file to new file
        ELSE                                               ! WNOBIN .LE. WNOASC
          IF ( WNOASC - WNOBIN .LE. WNOCHK )
     &      CALL DUPCHK ( LUNNEW, IREC, WNOCHK, LSTBIN, MOLBIN, 
     &        ISOBIN, WNOBIN, IUSBIN, ILSBIN, USLBIN, LSLBIN, INSERT )
          IF ( INSERT ) THEN
            WRITE ( LUNNEW, REC=IREC ) 
     &        LSTBIN, MOLBIN, ISOBIN, WNOBIN, STRBIN, TPRBIN,
     &        ABRBIN, SBRBIN, ELSBIN, ABCBIN, TSPBIN, IUSBIN, 
     &        ILSBIN, USLBIN, LSLBIN, ARFBIN, IFPBIN
            IREC = IREC + 1
            WNOUPP = WNOBIN
          END IF
          LSTBIN = 0
          DO WHILE ( LSTBIN .LT. 10 ) 
            IRECO = IRECO + 1
            IF ( IRECO .GT. IREC2O ) THEN
              WNOBIN = WNOEND
              GOTO 400
            END IF
            READ ( LUNOLD, REC=IRECO ) LSTBIN
          END DO
          READ ( LUNOLD, REC=IRECO ) 
     &      LSTBIN, MOLBIN, ISOBIN, WNOBIN, STRBIN, TPRBIN, 
     &      ABRBIN, SBRBIN, ELSBIN, ABCBIN, TSPBIN, IUSBIN, 
     &      ILSBIN, USLBIN, LSLBIN, ARFBIN, IFPBIN
 400      CONTINUE
        END IF
        IF ( MOD ( IREC, 100000 ) .EQ. 0 ) 
     &    WRITE ( *, '(A,I8,A,F12.6)' ) 
     &    'I-HITBIN: Record#', IREC, ' Wavenumber=', WNOUPP
      END DO
C
C Save record after last data record
      IREC2 = IREC
C
C Write header record
      WRITE ( LUNNEW, REC=1 ) 0, MAXMOL, IREC1, IREC2, HEAD48
C
C Write termination record
      WRITE ( LUNNEW, REC=IREC2 ) -2, 0, 0, WNOUPP, (0,IMOL=1,14)
      WRITE ( *, '(A,I8,A,F12.6)' )
     &  'I-HITBIN: Last Record#', IREC2, ' Wavenumber=', WNOUPP
C
C Calculate forward pointers
C
      WRITE ( *, '(A)' ) 'I-HITBIN: Writing forward pointers...'
      DO IMOL = 1, MAXMOL
        IRCMOL(IMOL) = IREC2
      END DO
      IFPREC = NFPREC
      DO IREC = IREC2-1, IREC1, -1 
        READ ( LUNNEW, REC=IREC ) LSTBIN
        IF ( LSTBIN .EQ. -7 ) THEN
          IOFF = ( IFPREC - 1 ) * MAXNFP
          WRITE ( LUNNEW, REC=IREC ) LSTBIN, IOFF+1, 0, WNOUPP, 
     &      ( IRCMOL(IMOL)-IREC, IMOL = IOFF+1, IOFF+MAXNFP ) 
          IF ( IFPREC .EQ. 1 ) THEN
            IFPREC = NFPREC
          ELSE
            IFPREC = IFPREC - 1
          END IF
        ELSE
          READ ( LUNNEW, REC=IREC ) 
     &      LSTBIN, MOLBIN, ISOBIN, WNOBIN, STRBIN, TPRBIN, 
     &      ABRBIN, SBRBIN, ELSBIN, ABCBIN, TSPBIN, IUSBIN, 
     &      ILSBIN, USLBIN, LSLBIN, ARFBIN, IFPBIN
          IFPBIN = IRCMOL(MOLBIN) - IREC 
          IRCMOL(MOLBIN) = IREC
          WNOUPP = WNOBIN
          WRITE ( LUNNEW, REC=IREC ) 
     &      LSTBIN, MOLBIN, ISOBIN, WNOBIN, STRBIN, TPRBIN, 
     &      ABRBIN, SBRBIN, ELSBIN, ABCBIN, TSPBIN, IUSBIN, 
     &      ILSBIN, USLBIN, LSLBIN, ARFBIN, IFPBIN
        END IF
        IF ( MOD ( IREC, 100000 ) .EQ. 0 ) 
     &    WRITE ( *, '(A,I8,A,F12.6)' ) 
     &    'I-HITBIN: Record#', IREC, ' Wavenumber=', WNOUPP
      END DO
C
      STOP  'R-HITBIN: Successful completion'
      END
C
      SUBROUTINE USERIN ( LUNASC, LUNOLD, LUNNEW, MAXMOL, WNO1, WNO2, 
     &                    USEMOL, MERGE, NEWSTT, RECLEN, HEAD48 )
C
C DESCRIPTION
C     Get user-inputs for HITBIN and open files
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER LUNASC         !  I  LUN for ASCII file
      INTEGER LUNOLD         !  I  LUN for old binary file
      INTEGER LUNNEW         !  I  LUN for new binary file
      INTEGER MAXMOL         !  I  Max.number of different molecules
      DOUBLE PRECISION WNO1  !  O  Lower wavenumber for selection from ASCII 
      DOUBLE PRECISION WNO2  !  O  Upper wavenumber for selection from ASCII
      LOGICAL USEMOL(MAXMOL) ! I/O T=use molecule, F=don't, assume=T on input
      LOGICAL MERGE          !  O  T=merge with old binary, F=ASCII input only
      INTEGER NEWSTT         !  O  Status for new line data
      INTEGER RECLEN         !  O  Record length (bytes or 4-byte words)
      CHARACTER*48 HEAD48    !  O  Header for new binary file
C
C LOCAL CONSTANTS
      INTEGER MAXIDX                 ! Max no of user-supplied indices
        PARAMETER ( MAXIDX = 100 )
      DOUBLE PRECISION WNOBIG        ! Larger than any reasonable wavenumber
        PARAMETER ( WNOBIG = 1.0D6 ) ! NB smaller than WNOEND in MAIN
C
C LOCAL VARIABLES
      INTEGER      IDUMMY  ! Dummy integer
      INTEGER      IIDX    ! Counter for list of molecule IDs
      INTEGER      IMOL    ! HITRAN Molecule ID
      INTEGER      IRECO   ! Record# in old file
      INTEGER      IREC1O  ! First data record in old file
      INTEGER      IREC2O  ! Last data record in old file
      INTEGER      LSTATO  ! Status of old data record
      INTEGER      LSTMAX  ! Maximum value of LSTATO
      INTEGER      MOLIDX(MAXIDX) ! List of molecule IDs from user
      LOGICAL      USE     ! T=use, F=exclude molecule
      CHARACTER*80 FILNAM  ! Name of file
      CHARACTER*80 RECORD  ! Record read from terminal
C
C DATA STATEMENTS
      DATA MOLIDX / MAXIDX * 0 /
      DATA IDUMMY   / 0 /
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Get name of ASCII file and open it
      WRITE ( *, '($,A)' ) 'Input HITRAN ASCII file: '
      READ ( *, '(A)' ) FILNAM
      OPEN ( UNIT=LUNASC, FILE=FILNAM, STATUS='OLD' )
C
C Get wavenumber range for ASCII file (default = use all)
      WRITE ( *, '($,A)' ) 'Wavenumber range (cm-1) [<CR>=all]: '
      READ ( *, '(A)' ) RECORD
      IF ( RECORD .EQ. ' ' ) THEN
        WNO1 = -1.0D0
        WNO2 = WNOBIG
      ELSE
        READ ( RECORD, * ) WNO1, WNO2
      END IF
C
C Get list of molecules to include/exclude (default = include all)
      WRITE ( *, '($,A)' ) 
     &  'HITRAN ID#s of gases to use (-ve=exclude) [<CR>=use all]: '
      READ ( *, '(A)' ) RECORD
      READ ( RECORD, *, END=100 ) 
     &  ( MOLIDX(IIDX), IIDX = 1, MAXIDX ), IDUMMY
      IF (IDUMMY.NE.0) STOP 'F-USERIN: No.indices > MAXIDX dimension'
 100  CONTINUE
      IF ( MOLIDX(1) .GT. 0 ) THEN  ! Include only user-supplied molecules
        DO IMOL = 1, MAXMOL         ! so change default to FALSE for all mols.
          USEMOL(IMOL) = .FALSE.
        END DO
      END IF
      IIDX = 1
      DO WHILE ( MOLIDX(IIDX) .NE. 0 )
        IF ( MOLIDX(IIDX) * MOLIDX(1) .LT. 0 ) 
     &    STOP 'F-USERIN: All ID#s must be positive or all negative'
        IMOL = ABS ( MOLIDX(IIDX) )
        USE  = MOLIDX(IIDX) .GT. 0
        IF ( IMOL .GT. MAXMOL ) 
     &    STOP 'F-USERIN: ID# > highest allowed index MAXMOL' 
        IF ( USE .EQV. USEMOL(IMOL) ) 
     &    STOP 'F-USERIN: Repeated ID# in list'
        USEMOL(IMOL) = USE
        IIDX = IIDX + 1
      END DO
C
C Get record length for new binary file
      WRITE ( *, '($,A)' ) 
     &  'Record Length for binary file (22=DEC, 88=others) [<CR>=88]: '
      READ ( *, '(A)' ) RECORD
      RECLEN = 88
      IF ( RECORD .NE. ' ' ) READ ( RECORD, * ) RECLEN
      IF ( RECLEN .NE. 22 .AND. RECLEN .NE. 88 ) 
     &  STOP 'F-USERIN: Value must be either 22 or 88'
C
C Get name of the new binary file to be created
      WRITE ( *, '($,A)' ) 'New binary file: '
      READ ( *, '(A)' ) FILNAM
      OPEN ( UNIT=LUNNEW, FILE=FILNAM, STATUS='NEW', ACCESS='DIRECT',
     &       RECL=RECLEN )
C
C Get header for new file 
      WRITE ( *, '($,A)' ) 'Header for new file (up to 48 chars): '
      READ ( *, '(A)' ) HEAD48
C      
C Get name of any existing binary file for merge
      WRITE ( *, '($,A)' ) 
     &  'Merge with existing binary file, filename [<CR>=none]: '
      READ ( *, '(A)' ) FILNAM
      MERGE = FILNAM .NE. ' '
      IF ( MERGE ) 
     &  OPEN ( UNIT=LUNOLD, FILE=FILNAM, STATUS='OLD', ACCESS='DIRECT',
     &         RECL=RECLEN )
C
C Get status to be associated with lines from ASCII data file
      IF ( MERGE ) THEN
        NEWSTT = 1          ! Use 1 as flag for increment old status by 1
        WRITE ( *, '($,A)' ) 
     &    'Status of new ASCII data (.GE.10) [<CR>=increment]: '
      ELSE
        NEWSTT = 10         ! Default lowest status for new file
        WRITE ( *, '($,A)' ) 
     &    'Status of new ASCII data (.GE.10) [<CR>=10]: '
      END IF
      READ ( *, '(A)' ) RECORD
      IF ( RECORD .NE. ' ' ) THEN
        READ ( RECORD, * ) NEWSTT
        IF ( NEWSTT .LE. 9 ) STOP 'F-USERIN: Status must be .GE. 10'
      END IF
C
      IF ( NEWSTT .EQ. 1 ) THEN
        WRITE ( *, '(A)' ) 
     &    'I-USERIN: Searching old file for highest Status...'
        READ ( LUNOLD, REC=1 ) IDUMMY, IDUMMY, IREC1O, IREC2O
        LSTMAX = 0
        DO IRECO = IREC1O, IREC2O
          READ ( LUNOLD, REC=IRECO ) LSTATO
          IF ( LSTATO .GT. LSTMAX ) THEN
            WRITE ( *, '(A,I8,A,I3)' ) 
     &        'I-USERIN: Old record#', IRECO, ' Status=', LSTATO
            LSTMAX = LSTATO
          END IF
        END DO
        NEWSTT = LSTMAX + 1
        WRITE ( *, '(A,I3)' ) 'I-USERIN: Setting new Status=', NEWSTT
      END IF
C
      END
C
      SUBROUTINE NEWHDR ( LUNASC, LUNOLD, LUNNEW, HEADR2, MERGE, IREC )
C
C DESCRIPTION
C     Copy header records into new binary file
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER LUNASC         !  I  LUN for ASCII file
      INTEGER LUNOLD         !  I  LUN for old binary file
      INTEGER LUNNEW         !  I  LUN for new binary file
      CHARACTER*84 HEADR2    !  I  HITBIN header for this file
      LOGICAL MERGE          !  I  T=merge with old binary, F=ASCII to binary
      INTEGER IREC           !  O  Record# in new file after headers
C
C LOCAL VARIABLES
      INTEGER LSTAT       ! Record ID (4=header record)
      INTEGER IDUMMY      ! Dummy integer
      CHARACTER*84 HEADER ! Header record
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Write HITBIN header (record length and version ID) to rec#2
      WRITE ( LUNNEW, REC=2 ) 4, HEADR2
C
C Copy headers from old binary file to rec#2 onwards in new binary file
      IREC = 3
      IF ( MERGE ) THEN
        WRITE ( *, '(A)' ) 
     &    'I-NEWHDR: Copying headers from binary file...'
        HEADER = ' '
        READ (LUNOLD,REC=1) LSTAT, IDUMMY, IDUMMY, IDUMMY, HEADER(1:48)
        LSTAT = 4
        DO WHILE ( LSTAT .EQ. 4 )
          WRITE ( LUNNEW, REC=IREC ) 4, HEADER      
          WRITE ( *, '(A)' ) HEADER
          IREC = IREC + 1
          READ ( LUNOLD, REC=IREC ) LSTAT, HEADER
        END DO
      END IF
C
C Then copy headers from ASCII file into new binary file
      HEADER = '!'
      WRITE ( *, '(A)' ) 'I-NEWHDR: Copying headers from ASCII file...'
      DO WHILE ( HEADER(1:1) .EQ. '!' )
        READ ( LUNASC, '(A)' ) HEADER 
        IF ( HEADER(1:1) .EQ. '!' ) THEN
          WRITE ( LUNNEW, REC=IREC ) 4, HEADER
          WRITE ( *, '(A)' ) HEADER
          IREC = IREC + 1       
        END IF
      END DO
      BACKSPACE ( LUNASC )
C
      END
C
      SUBROUTINE DUPCHK ( LUN, IREC, WNOCHK, LST, MOL, ISO, WNO, 
     &                    IUS, ILS, USL, LSL, INSERT )
C
C DESCRIPTION
C     Check back in new binary file for duplicates of current record
C     If duplicate found:
C       with lower status than current line, remove & decrease IREC by 1
C       with higher status than current line, set INSERT = FALSE
C       with same status, keep old line but issue warning.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER          LUN    !  I  LUN of new binary file
      INTEGER          IREC   ! I/O Next record of new binary file
      DOUBLE PRECISION WNOCHK !  I  Wavenumber range to check
      INTEGER          LST    !  I  LSTAT value of current line
      INTEGER          MOL    !  I  Molecule# of current line
      INTEGER          ISO    !  I  Isotope# of current line
      DOUBLE PRECISION WNO    !  I  Wavenumber of current line
      INTEGER          IUS    !  I  Upper State Global Quantum Number 
      INTEGER          ILS    !  I  Lower State Global Quantum Number 
      CHARACTER*9      USL    !  I  Upper State Local Quantum Number 
      CHARACTER*9      LSL    !  I  Lower State Local Quantum Number 
      LOGICAL          INSERT !  O  T=insert current line, F=ignore 
C
C LOCAL CONSTANTS
      INTEGER    MAXWRN       !  Max number of warning messages issued
        PARAMETER ( MAXWRN = 10 )
C
C LOCAL VARIABLES
      INTEGER  IWRN    !  Warning message counter
      INTEGER  JREC    !  Record# being checked
      INTEGER  LSTPRV, MOLPRV, ISOPRV, IUSPRV, ILSPRV, IFPPRV
      INTEGER  KREC    !  Record#
      INTEGER  LREC    !  Last line record# overwritten
      REAL     STRPRV, TPRPRV, ABRPRV, SBRPRV, ELSPRV, ABCPRV, TSPPRV
      DOUBLE PRECISION WNOPRV 
      CHARACTER*9 LSLPRV, USLPRV, ARFPRV
C
C DATA STATEMENTS
      DATA IWRN / 0 /
      SAVE IWRN
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      INSERT = .TRUE.
      JREC = IREC
      WNOPRV = WNO
      DO WHILE ( WNO - WNOPRV .LT. WNOCHK .AND. JREC .GT. 1 )
        JREC = JREC - 1
        READ ( LUN, REC=JREC ) LSTPRV
        IF ( LSTPRV .GE. 10 ) THEN
          READ ( LUN, REC=JREC ) LSTPRV, MOLPRV, ISOPRV, WNOPRV,
     &        STRPRV, TPRPRV, ABRPRV, SBRPRV, ELSPRV, ABCPRV, TSPPRV,
     &        IUSPRV, ILSPRV, USLPRV, LSLPRV, ARFPRV, IFPPRV
          IF ( MOL .EQ. MOLPRV .AND. ISO .EQ. ISOPRV .AND.
     &         IUS .EQ. IUSPRV .AND. ILS .EQ. ILSPRV .AND.
     &         USL .EQ. USLPRV .AND. LSL .EQ. LSLPRV       ) THEN ! same line
            IF ( LST .GT. LSTPRV ) THEN            ! Remove earlier line
              LREC = JREC
              DO KREC = JREC+1, IREC-1
                READ ( LUN, REC=KREC ) LSTPRV
                IF ( LSTPRV .GE. 10 ) THEN
                  READ ( LUN, REC=KREC ) 
     &              LSTPRV, MOLPRV, ISOPRV, WNOPRV, STRPRV, TPRPRV, 
     &              ABRPRV, SBRPRV, ELSPRV, ABCPRV, TSPPRV, IUSPRV, 
     &              ILSPRV, USLPRV, LSLPRV, ARFPRV, IFPPRV
                  WRITE ( LUN, REC=LREC ) 
     &              LSTPRV, MOLPRV, ISOPRV, WNOPRV, STRPRV, TPRPRV, 
     &              ABRPRV, SBRPRV, ELSPRV, ABCPRV, TSPPRV, IUSPRV, 
     &              ILSPRV, USLPRV, LSLPRV, ARFPRV, IFPPRV
                  LREC = KREC
                END IF
              END DO
              IREC = LREC 
            ELSE IF ( LST .LT. LSTPRV ) THEN            ! Ignore current line
              INSERT = .FALSE.                    
            ELSE IF ( IWRN .LT. MAXWRN .AND.  ! Warn, but keep old & new lines
     &                ( WNO .GT. 160.0D0 .OR. .NOT.  ! not O3(4) dupl.<160cm-1
     &                  ( MOL .EQ. 3 .AND. ISO .EQ. 4 ) ) ) THEN 
              IWRN = IWRN + 1                
              WRITE ( *, '(A)' ) 'W-DUPCHK: Duplicate lines ?'//
     &          ' of equal status - keeping both'
              WRITE ( *, '(A,2I2,2F12.6)' ) 
     &          'Mol,Iso,Wno1,Wno2=', MOLPRV, ISOPRV, WNOPRV, WNO
              IF ( IWRN .EQ. MAXWRN ) WRITE ( *, '(A)' ) 
     &          'W-DUPCHK: Further warnings suppressed'
            END IF                       ! else keep old & new lines
          END IF
        END IF
      END DO
C
      END
