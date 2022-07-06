      SUBROUTINE GRDDEF ( LUNGRD, NAMDEF, FAIL, ERRMSG )
C
C VERSION
C     23-JUL-03  AD  Original.
C
C DESCRIPTION
C     Use default .GRD filename to find any missing files.
C     Called once by INPGRD if filename template found in *GRD section of 
C     driver table.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER       LUNGRD ! I/O LUN for GRD File/next free LUN
      CHARACTER*(*) NAMDEF !  I  Default name of GRD file
      LOGICAL       FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  GRDFIL ! Open GRD file and check contents
     &, RFMLOG ! Write text message to RFM log file.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'gflcom.inc' ! Grid file data
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
C
C LOCAL VARIABLES
      INTEGER       IGFL   ! Counter for current GRD files
      INTEGER       ILEN   ! Length of NAMDEF
      INTEGER       IOS    ! Saved value of IOSTAT for error messages
      INTEGER       IPT    ! Pointer to '*' character in template
      INTEGER       ISPC   ! Spectral range counter
      CHARACTER*132 MESSGE ! Text message sent to LOG file
      CHARACTER*80  NAMGRD ! Name of GRD file
      INTEGER       JPT    ! Pointer to last character in spectral range label
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Send message to LOG file saying which file is about to be opened
C
      MESSGE = 'I-GRDDEF: Checking default .grd files'
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
C
      IPT = INDEX ( NAMDEF, '*' )
      ILEN = LEN ( NAMDEF )
      IF ( ILEN .GT. 80 ) STOP 'F-GRDDEF: Logical error'
      
      DO ISPC = 1, NSPC                            ! Loop over spectral ranges
        DO IGFL = 1, NGFL
          IF ( SPCGFL(IGFL) .EQ. ISPC ) GOTO 100   ! .grd file already found
        END DO
C At this point, spectral range found requiring .grd data but no data loaded
        JPT = INDEX ( LABSPC(ISPC)//' ', ' ' ) - 1
        IF ( IPT .EQ. 1 .AND. ILEN .EQ. 1 ) THEN   ! NAMDEF = '*'  only
          NAMGRD = LABSPC(ISPC)(1:JPT)
        ELSE IF ( IPT .EQ. 1 ) THEN                ! NAMDEF = '*...'  
          NAMGRD = LABSPC(ISPC)(1:JPT)//NAMDEF(2:)       
        ELSE IF ( IPT .EQ. ILEN ) THEN             ! NAMDEF = '...*'  
          NAMGRD = NAMDEF(1:IPT-1)//LABSPC(ISPC)(1:JPT)
        ELSE                                     ! NAMDEF = '...*...' (usual)
          NAMGRD = NAMDEF(1:IPT-1)//LABSPC(ISPC)(1:JPT)//NAMDEF(IPT+1:)
        ENDIF
        MESSGE = 'I-GRDDEF: looking for file: '//NAMGRD
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        OPEN ( UNIT=LUNGRD, FILE=NAMGRD, STATUS='OLD', 
     &           IOSTAT=IOS, ERR=100 )      ! Skip if any problems opening
        CLOSE ( LUNGRD )                    ! Reset for GRDFIL
C Read file contents and load if useful (incrementing LUNGRD)
        CALL GRDFIL ( LUNGRD, NAMGRD, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
 100    CONTINUE
      END DO
C
      FAIL = .FALSE.
C
      END
