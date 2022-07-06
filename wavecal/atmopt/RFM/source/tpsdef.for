      SUBROUTINE TPSDEF ( LUNTPS, NAMDEF, FAIL, ERRMSG )
C
C VERSION
C     23-JUL-03  AD  Original.
C
C DESCRIPTION
C     Use default .tps filename to find any missing files.
C     Called once by INPTPS if filename template found in *TPS section of 
C     driver table.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER       LUNTPS ! I/O LUN for TPS File/next free LUN
      CHARACTER*(*) NAMDEF !  I  Default name of TPS file
      LOGICAL       FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  RFMLOG ! Write text message to RFM log file.
     &, TPSFIL ! Open TIPS data file, check and load contents.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'tpscom.inc' ! Tabulated TIPS data
C
C LOCAL VARIABLES
      INTEGER       IGAS   ! Counter for absorbing species
      INTEGER       ILEN   ! Length of NAMDEF
      INTEGER       IOS    ! Saved value of IOSTAT for error messages
      INTEGER       IPT    ! Pointer to '*' character in template
      INTEGER       ITPS   ! Counter for current TPS tabulations
      INTEGER       JPT    ! Pointer to last character in absorber name
      CHARACTER*132 MESSGE ! Text message sent to LOG file
      CHARACTER*80  NAMTPS ! Name of TPS file
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Send message to LOG file saying which file is about to be opened
C
      MESSGE = 'I-TPSDEF: Checking default .tps files'
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
      IPT = INDEX ( NAMDEF, '*' )
      ILEN = LEN ( NAMDEF )
      IF ( ILEN .GT. 80 ) STOP 'F-TPSDEF: Logical error'
      
      DO IGAS = 1, NGAS                            ! Loop over absorbers
        DO ITPS = 1, NTPS                          ! Loop over TPS data
          IF ( IDGTPS(ITPS) .EQ. IGSMOL(IGAS) ) GOTO 100 ! data already loaded
        END DO
C At this point, absorber found requiring .tps data but no data yet loaded
        JPT = INDEX ( CODGAS(IGAS)//' ', ' ' ) - 1
        IF ( IPT .EQ. 1 .AND. ILEN .EQ. 1 ) THEN   ! NAMDEF = '*'  only
          NAMTPS = CODGAS(IGAS)(1:JPT)
        ELSE IF ( IPT .EQ. 1 ) THEN                ! NAMDEF = '*...'  
          NAMTPS = CODGAS(IGAS)(1:JPT)//NAMDEF(2:)       
        ELSE IF ( IPT .EQ. ILEN ) THEN             ! NAMDEF = '...*'  
          NAMTPS = NAMDEF(1:IPT-1)//CODGAS(IGAS)(1:JPT)
        ELSE                                     ! NAMDEF = '...*...' (usual)
          NAMTPS = NAMDEF(1:IPT-1)//CODGAS(IGAS)(1:JPT)//
     &             NAMDEF(IPT+1:)
        ENDIF
        MESSGE = 'I-TPSDEF: looking for file: '//NAMTPS
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        OPEN ( UNIT=LUNTPS, FILE=NAMTPS, STATUS='OLD', 
     &           IOSTAT=IOS, ERR=100 )      ! Skip if any problems opening
        CLOSE ( LUNTPS )
C Reopen and read in data
        CALL TPSFIL ( LUNTPS, NAMTPS, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
 100    CONTINUE
      END DO
C
      FAIL = .FALSE.
      RETURN
C
      END
