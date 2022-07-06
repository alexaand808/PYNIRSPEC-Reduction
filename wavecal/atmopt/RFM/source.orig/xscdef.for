      SUBROUTINE XSCDEF ( LUNXSC, NAMDEF, FAIL, ERRMSG )
C
C VERSION
C     20-SEP-12  AD  Change NAMXSC from C*80 to C*200
C     21-JUL-03  AD  Original.
C
C DESCRIPTION
C     Use default .xsc filename to find any missing files.
C     Called once by INPXSC if filename template found in *XSC section of 
C     driver table.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER       LUNXSC ! I/O LUN for XSC File/next free LUN
      CHARACTER*(*) NAMDEF !  I  Default name of XSC file
      LOGICAL       FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  FILXFL ! Check X/S data input file and store in /XFLCOM/
     &, RFMLOG ! Write text message to RFM log file.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'shpcon.inc' ! Line-shape codes
C
C COMMON VARIABLES
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'xsccom.inc' ! X-section data
C
C LOCAL VARIABLES
      INTEGER       IGAS   ! Counter for absorbing species
      INTEGER       ILEN   ! Length of NAMDEF
      INTEGER       IOS    ! Saved value of IOSTAT for error messages
      INTEGER       IPT    ! Pointer to '*' character in template
      INTEGER       JPT    ! Pointer to last character in absorber name
      INTEGER       IXSC   ! Counter for current xsc tabulations
      LOGICAL       REJECT ! T=ignore file - no useful data
      CHARACTER*132 MESSGE ! Text message sent to LOG file
      CHARACTER*200 NAMXSC ! Name of XSC file
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Send message to LOG file saying which file is about to be opened
C
      MESSGE = 'I-XSCDEF: Checking default .xsc files'
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
      IPT = INDEX ( NAMDEF, '*' )
      ILEN = LEN ( NAMDEF )
      IF ( ILEN .GT. 200 ) STOP 'F-XSCDEF: Logical error'
      
      DO IGAS = 1, NGAS                            ! Loop over absorbers
        IF ( SHPGAS(IGAS) .EQ. SHPXSC ) THEN
          DO IXSC = 1, NXSC                        ! Loop over xsc data
            IF ( IGSXSC(IXSC) .EQ. IGAS ) GOTO 100 ! data already loaded
          END DO
C At this point, absorber found requiring .xsc data but no data yet loaded
          JPT = INDEX ( CODGAS(IGAS)//' ', ' ' ) - 1
          IF ( IPT .EQ. 1 .AND. ILEN .EQ. 1 ) THEN   ! NAMDEF = '*'  only
            NAMXSC = CODGAS(IGAS)(1:JPT)
          ELSE IF ( IPT .EQ. 1 ) THEN                ! NAMDEF = '*...'  
            NAMXSC = CODGAS(IGAS)(1:JPT)//NAMDEF(2:)       
          ELSE IF ( IPT .EQ. ILEN ) THEN             ! NAMDEF = '...*'  
            NAMXSC = NAMDEF(1:IPT-1)//CODGAS(IGAS)(1:JPT)
          ELSE                                     ! NAMDEF = '...*...' (usual)
            NAMXSC = NAMDEF(1:IPT-1)//CODGAS(IGAS)(1:JPT)//
     &               NAMDEF(IPT+1:)
          ENDIF
          MESSGE = 'I-XSCDEF: looking for file: '//NAMXSC
          CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
          OPEN ( UNIT=LUNXSC, FILE=NAMXSC, STATUS='OLD', 
     &           IOSTAT=IOS, ERR=100 )      ! Skip if any problems opening
C Check file is useful, and store files in *XFL 
          CALL FILXFL ( LUNXSC, REJECT, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
          IF ( REJECT ) THEN
            CLOSE ( LUNXSC )
          ELSE 
            LUNXSC = LUNXSC + 1
          END IF
        ENDIF
 100    CONTINUE
      END DO
C
      FAIL = .FALSE.
      RETURN
C
      END
