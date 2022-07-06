      SUBROUTINE NTEDEF ( LUNNTE, NAMDEF, FAIL, ERRMSG )
C
C VERSION
C     23-JUL-03  AD  Original.
C
C DESCRIPTION
C     Use default .nte filename to find any missing files.
C     Called once by INPNTE if filename template found in *NTE section of 
C     driver table.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER       LUNNTE ! I/O LUN for NTE File/next free LUN
      CHARACTER*(*) NAMDEF !  I  Default name of NTE file
      LOGICAL       FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  NTEFIL ! Read NTE Vibrational Temperature Data from file
     &, RFMLOG ! Write text message to RFM log file.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'ntecom.inc' ! Non-LTE data
C
C LOCAL VARIABLES
      INTEGER       IGAS   ! Counter for absorbing species
      INTEGER       ILEN   ! Length of NAMDEF
      INTEGER       IOS    ! Saved value of IOSTAT for error messages
      INTEGER       IPT    ! Pointer to '*' character in template
      INTEGER       JPT    ! Pointer to last character in absorber name
      INTEGER       INTE   ! Counter for current NTE tabulations
      CHARACTER*132 MESSGE ! Text message sent to LOG file
      CHARACTER*80  NAMNTE ! Name of NTE file
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Send message to LOG file saying which file is about to be opened
C
      MESSGE = 'I-NTEDEF: Checking default .nte files'
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
      IPT = INDEX ( NAMDEF, '*' )
      ILEN = LEN ( NAMDEF )
      IF ( ILEN .GT. 80 ) STOP 'F-NTEDEF: Logical error'
      
      DO IGAS = 1, NGAS                            ! Loop over absorbers
        DO INTE = 1, NNTE                          ! Loop over NTE data
          IF ( IDGNTE(INTE) .EQ. IGSMOL(IGAS) ) GOTO 100 ! data already loaded
        END DO
C At this point, absorber found requiring .nte data but no data yet loaded
        JPT = INDEX ( CODGAS(IGAS)//' ', ' ' ) - 1
        IF ( IPT .EQ. 1 .AND. ILEN .EQ. 1 ) THEN   ! NAMDEF = '*'  only
          NAMNTE = CODGAS(IGAS)(1:JPT)
        ELSE IF ( IPT .EQ. 1 ) THEN                ! NAMDEF = '*...'  
          NAMNTE = CODGAS(IGAS)(1:JPT)//NAMDEF(2:)       
        ELSE IF ( IPT .EQ. ILEN ) THEN             ! NAMDEF = '...*'  
          NAMNTE = NAMDEF(1:IPT-1)//CODGAS(IGAS)(1:JPT)
        ELSE                                     ! NAMDEF = '...*...' (usual)
          NAMNTE = NAMDEF(1:IPT-1)//CODGAS(IGAS)(1:JPT)//
     &             NAMDEF(IPT+1:)
        ENDIF
        MESSGE = 'I-NTEDEF: looking for file: '//NAMNTE
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        OPEN ( UNIT=LUNNTE, FILE=NAMNTE, STATUS='OLD', 
     &           IOSTAT=IOS, ERR=100 )      ! Skip if any problems opening
        CLOSE ( LUNNTE )
C Reopen and read in data
        CALL NTEFIL ( LUNNTE, NAMNTE, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
 100    CONTINUE
      END DO
C
      FAIL = .FALSE.
      RETURN
C
      END
