      SUBROUTINE HITFIL ( LUNHIT, NAMHIT, FAIL, ERRMSG )
C
C VERSION
C     07-AUG-13  AD  Add OPNHIT to do actual opening of file
C                    Add IFORM and check file format identifier
C     09-JAN-08  AD  Remove special treatment of HDO
C     20-OCT-03  AD  Separate error messages for Open and Read failures
C     09-OCT-02  AD  Correct limits for checking HITRAN overlap
C     13-JUN-02  AD  Remove test on whether species present in file
C     12-JUN-02  AD  Simplify to use only single HITRAN file. Remove FILHFL
C     16-FEB-02  AD  Allow for HDO
C     11-DEC-01  AD  Remove READONLY. 
C                    Rearrange logic of reading forward pointers to allow 
C                    arbitrary number of f.p. records and check .le. MAXPTR
C     05-DEC-00  AD  Comment out READONLY
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-OCT-96  AD  Rearrange to avoid compiler bug when reading in FWDPTR
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Open HITRAN binary line data file read from Driver Table
C     Called once by INPHIT if *HIT section of driver table is present.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER       LUNHIT  ! I/O LUN for HIT File/next free LUN
      CHARACTER*(*) NAMHIT  !  I  Name of HIT file
      LOGICAL       FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  OPNHIT !  Open HITRAN binary line data file 
     &, RFMLOG ! Write text message to RFM log file
C
C GLOBAL CONSTANTS
      INCLUDE 'idxcon.inc' ! HITRAN molecular indices
      INCLUDE 'rfmcon.inc' ! Constants for RFM calculations
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'hflcom.inc' ! HITRAN line data file
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
C
C LOCAL VARIABLES
      INTEGER       FWDPTR(MAXPTR) ! Forward pointers for each HITRAN line ID
      INTEGER       FPTEMP(14)  ! Temporary buffer for forward pointers
      INTEGER       IDUMMY      ! Dummy integer for read
      INTEGER       IFORM       ! Binary file format identifier
      INTEGER       IGAS        ! Absorber counter
      INTEGER       IHLN        ! HITRAN line molecule index 
      INTEGER       IHLN1       ! HITRAN line molecule index of first pointer
      INTEGER       IOS         ! Saved value of IOSTAT for error messages
      INTEGER       IPTR        ! Pointer counter (1:14)
      INTEGER       IREC        ! Record counter
      INTEGER       LSTAT       ! Used to identify HITRAN Fwd.Pointer block
      CHARACTER*48  LABEL       ! Label read from HITRAN file
      CHARACTER*111 MESSGE      ! Text message sent to LOG file (C*31 + C*80)
      DOUBLE PRECISION DDUMMY   ! Dummy R*8
C
C DATA STATEMENTS
      DATA FWDPTR / MAXPTR * 0 /
      SAVE FWDPTR
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Open file
      CALL OPNHIT ( LUNHIT, NAMHIT, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
C Read header and extract first/last record#
      READ ( LUNHIT, REC=1, IOSTAT=IOS, ERR=900 )
     &  LSTAT, IFORM, IR1HFL, IR2HFL, LABEL
C 
C Print File label to LOG file
      MESSGE = 'I-HITFIL: HITRAN File Label='//LABEL
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
C Check format - data generated with old version of HITBIN had IFORM=56,
C newer data has IFORM=1 but format is the same, and the only one currently
C handled by the RFM
      IF ( IFORM .NE. 1 .AND. IFORM .NE. 56 ) THEN
        WRITE ( ERRMSG, '(A,I11)' )
     &    'F-HITFIL: Unexpected format ID in rec#1 of binary file,'//
     &    ' value=', IFORM 
        FAIL = .TRUE.
        RETURN
      END IF
C
C Find wavenumber range by reading first,last records
      READ ( LUNHIT, REC=IR1HFL, IOSTAT=IOS, ERR=900 )
     &  LSTAT, IDUMMY, IDUMMY, WN1HFL
      READ ( LUNHIT, REC=IR2HFL, IOSTAT=IOS, ERR=900 )
     &  LSTAT, IDUMMY, IDUMMY, WN2HFL
C
C Check HITRAN wavenumber range overlaps required spectral range(s). 
      IF ( WN1HFL .GT. WMXSPC .OR. WN2HFL .LT. WMNSPC ) THEN
        WRITE ( ERRMSG, '(A,2F10.4)' )
     &    'F-HITFIL: HITRAN file outside reqd waveno range, =',
     &    WN1HFL, WN2HFL
        FAIL = .TRUE.
        RETURN
      ELSE IF ( WN1HFL .GT. MAX ( WMNSPC-FWIND, 0.0D0 ) .OR. 
     &          WN2HFL .LT. WMXSPC+FWIND ) THEN
        WRITE ( MESSGE, '(A,F7.1,A)' )
     &    'W-HITFIL: HITRAN file may not include all lines within',
     &    FWIND, ' cm-1 of output'
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
C
C Initialise pointers for all different line species to point to last
C record in file and set IFPHFL=0 for any IGAS values used to indicate isotopes
      DO IGAS = 1, NGAS
        IHLN = IDXGAS(IGAS)
        IF ( IHLN .LE. MAXHLN .AND. IGAS .EQ. IGSMOL(IHLN) ) THEN
          IFPHFL(IGAS) = IR2HFL
        ELSE
          IFPHFL(IGAS) = 0
        END IF
      END DO
C
C Read Forward Pointers from 1st n records from start of file
      IREC = IR1HFL
      READ ( LUNHIT, REC=IREC, ERR=900, IOSTAT=IOS ) LSTAT
      IF ( LSTAT .NE. -7 ) THEN                 ! Check this is FP record
        WRITE ( ERRMSG, '(A,I11,A,I11)' )
     &  'F-HITFIL: Expected Fwd Ptr (LSTAT=-7) for Rec#', IREC,
     &  ', got LSTAT=', LSTAT
        FAIL = .TRUE.
        RETURN
      END IF
C NB: Can't read into FWDPTR directly since read & use of IHLN1 index in same
C READ statement causes problems to some (eg Sun) compilers, hence use FPTEMP.
      DO WHILE ( LSTAT .EQ. -7 )
        READ ( LUNHIT, REC=IREC, ERR=900, IOSTAT=IOS )
     &    LSTAT, IHLN1, IDUMMY, DDUMMY, ( FPTEMP(IPTR), IPTR = 1, 14 )
        IF ( IHLN1+13 .GT. MAXPTR ) THEN
          FAIL = .TRUE.
          WRITE ( ERRMSG, '(A,I5,A)' ) 
     &      'F-HITFIL: No.forward ptrs in HITRAN file > '//
     &      'MAXPTR=', MAXPTR, ' in rfmsiz.inc' 
          RETURN
        END IF
        DO IPTR = 1, 14
          IHLN = IHLN1 - 1 + IPTR
          FWDPTR(IHLN) = FPTEMP(IPTR)
          IGAS = IGSMOL(IHLN)
C IFPHFL is checked later in ABSCHK to see if lines of each species are
C present in the file, ie that forward pointer is before the end of file
          IF ( IGAS .NE. 0 ) IFPHFL(IGAS) = FPTEMP(IPTR)
        END DO
        IREC = IREC + 1
        READ ( LUNHIT, REC=IREC, ERR=900, IOSTAT=IOS ) LSTAT
      END DO
C    
      LUNHFL = LUNHIT
C
      FAIL = .FALSE.
      RETURN
C
  900 WRITE ( ERRMSG, '(A,I11)') 
     &  'F-HITFIL: Read failure on HITRAN file. IOSTAT=', IOS
      FAIL = .TRUE.
      RETURN
C
      END
