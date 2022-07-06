      SUBROUTINE LUTDEF ( LUNLUT, NAMDEF, FAIL, ERRMSG )
C
C VERSION
C     25-MAR-04  AD  Allow for isotopes in filenames
C     23-JUL-03  AD  Original.
C
C DESCRIPTION
C     Use default LUT filename to find any missing files.
C     Called once by INPLUT if filename template found in *LUT section of 
C     driver table.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER       LUNLUT ! I/O LUN for LUT File/next free LUN
      CHARACTER*(*) NAMDEF !  I  Default name of LUT file
      LOGICAL       FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  LUTFIL ! Open LUT file and check contents
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
      INTEGER       IGAS   ! Absorber counter
      INTEGER       ILEN   ! Length of NAMDEF
      INTEGER       ILFL   ! Counter for stored LUTs
      INTEGER       IOS    ! Saved value of IOSTAT for error messages
      INTEGER       IPT    ! Pointer to '*' character in template
      INTEGER       ISPC   ! Spectral range counter
      INTEGER       JPT    ! Pointer to last character in spectral range label
      INTEGER       KPT    ! Pointer to last character in molecule name
      INTEGER       LPT    ! Pointer to last character in SUBSTR
      LOGICAL       CSFILE ! T=SVD-compressed MIPAS LUT, F=RFM direct tab.
      CHARACTER*132 MESSGE ! Text message sent to LOG file
      CHARACTER*80  NAMLUT ! Name of LUT file
      CHARACTER*20  SUBSTR ! Combination of spectral range label and molecule
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Send message to LOG file saying which file is about to be opened
C
      MESSGE = 'I-LUTDEF: Checking default LUT files'
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
C
      IPT = INDEX ( NAMDEF, '*' )
      ILEN = LEN ( NAMDEF )
      IF ( ILEN .GT. 80 ) STOP 'F-LUTDEF: Logical error'
      CSFILE = INDEX ( NAMDEF, 'CS_' ) .GT. 0
C      
      DO ISPC = 1, NSPC                            ! Loop over spectral ranges
        DO IGAS = 1, NGAS                          ! Loop over gases
          DO ILFL = 1, NLFL                        ! Loop over existing LUTs
            IF ( SPCLFL(ILFL) .EQ. ISPC .AND.
     &           GASLFL(ILFL) .EQ. IGAS       ) GOTO 100 ! LUT already found
          END DO
C At this point, spectral range/gas found requiring LUT data but no data loaded
C Construct spec/gas part of LUT filename, which depends on whether this is
C a SVD compressed LUT (CSFILE=TRUE) or a direct tabulation (FALSE)
          JPT = INDEX ( LABSPC(ISPC)//' ', ' ' ) - 1
          IF ( CSFILE ) THEN                        ! always has spectral label
            SUBSTR = LABSPC(ISPC)(1:JPT)//'_'
            WRITE ( SUBSTR(JPT+2:JPT+3), '(I2.2)' ) IDXGAS(IGAS)
          ELSE
            KPT = INDEX ( CODGAS(IGAS)//' ', ' ' ) - 1
            IF ( JPT .GE. 1 ) THEN
              SUBSTR = LABSPC(ISPC)(1:JPT)//CODGAS(IGAS)(1:KPT)
            ELSE                                     ! no spectral label
              SUBSTR = CODGAS(IGAS)(1:KPT)
            END IF
          END IF
          IF ( ISOMOL(IDXGAS(IGAS)) ) THEN
            LPT = INDEX ( SUBSTR//' ', ' ' ) - 1
            WRITE ( SUBSTR(LPT+1:LPT+2), '(A,I1)' ) 'i', IDIGAS(IGAS)
          END IF          
          LPT = INDEX ( SUBSTR//' ', ' ' ) - 1
          IF ( IPT .EQ. 1 .AND. ILEN .EQ. 1 ) THEN   ! NAMDEF = '*'  only
            NAMLUT = SUBSTR(1:LPT)
          ELSE IF ( IPT .EQ. 1 ) THEN                ! NAMDEF = '*...'  
            NAMLUT = SUBSTR(1:LPT)//NAMDEF(2:)       
          ELSE IF ( IPT .EQ. ILEN ) THEN             ! NAMDEF = '...*'  
            NAMLUT = NAMDEF(1:IPT-1)//SUBSTR(1:LPT)
          ELSE                                     ! NAMDEF = '...*...' (usual)
            NAMLUT = NAMDEF(1:IPT-1)//SUBSTR(1:LPT)//NAMDEF(IPT+1:)
          ENDIF
          MESSGE = 'I-LUTDEF: looking for file: '//NAMLUT
          CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
            IF ( FAIL ) RETURN
          OPEN ( UNIT=LUNLUT, FILE=NAMLUT, STATUS='OLD', 
     &             IOSTAT=IOS, ERR=100 )      ! Skip if any problems opening
          CLOSE ( LUNLUT )                    ! Reset for LUTFIL
C Read file contents and load if useful (incrementing LUNLUT)
          CALL LUTFIL ( LUNLUT, NAMLUT, FAIL, ERRMSG )
            IF ( FAIL ) RETURN
 100      CONTINUE
        END DO
      END DO
C
      FAIL = .FALSE.
C
      END
