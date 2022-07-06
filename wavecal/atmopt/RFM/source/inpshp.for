      SUBROUTINE INPSHP ( LUNDRV, FAIL, ERRMSG )
C
C VERSION
C     17OCT13 AD Simplified: remove GETSHP argument and local GOTSHP checks
C     11AUG03 AD Add VVW shape
C     16FEB02 AD Allow for different isotopes 
C     31DEC98 AD Comment change only
C     01JUN98 AD Correct warning message for CO2: warn if CTM *not* 
C                  combined with chi-factor 
C     04DEC97 AD Add GOTSHP to check for duplicate section
C     01DEC97 AD Check for SHPCTM, also legitimate continuum shapes.
C                Send info message to LOG file about shapes being used.
C     03MAR97 AD Version 3.
C     01OCT96 AD Version 2.
C     01SEP96 AD Version 1.
C
C DESCRIPTION
C     Read RFM line shapes from driver table *SHP section.
C     Called by RFMINP once if SHPFLG set TRUE.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      LUNDRV  !  I  LUN for Driver File
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  LOCASE ! Convert text string to lower case.
     &, NXTREC ! Load next record from RFM input file.
     &, RFMLOG ! Write text message to RFM log file.
     &, TXTFLD ! Identify start and end points of text field in record
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'idxcon.inc' ! HITRAN molecular indices
      INCLUDE 'shpcon.inc' ! Line Shape codes
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
C
C LOCAL VARIABLES
      INTEGER      IEND        ! Pointer to end of filename in record
      INTEGER      IFLD        ! Text Field# within record
      INTEGER      IGAS        ! Absorber counter
      INTEGER      ISHP        ! Line shape index#
      INTEGER      ISO         ! Isotope counter
      INTEGER      ISTA        ! Pointer to start of filename in record
      INTEGER      JGAS        ! Index of isotopes in *GAS arrays
      INTEGER      SHPDEF      ! Default Line shape to be used
      LOGICAL      ENDSEC      ! T=End of driver table section reached.
      CHARACTER*7  CODE        ! Absorber name (C*7 is longest name in use)
      CHARACTER*80 LOGMSG      ! Text message for Log file
      CHARACTER*80 RECORD      ! Line of Driver file or External TAN file
      CHARACTER*3  SHAPE       ! Line shape name
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Each record should start with a recognised lineshape (VOIGT,LORENTZ etc)
C
      SHPDEF = 0
      CALL NXTREC ( LUNDRV, RECORD, ENDSEC, FAIL, ERRMSG ) 
      IF ( FAIL ) RETURN
      DO WHILE ( .NOT. ENDSEC ) 
        CALL TXTFLD ( 80, RECORD, 1, ISTA, IEND )    ! Locate Field#1 (shape) 
        IEND = MIN ( ISTA+2, IEND )                  ! Limit to C*3
        SHAPE = RECORD(ISTA:IEND)
        CALL LOCASE ( SHAPE, SHAPE )      
        IF ( SHAPE .EQ. 'voi' ) THEN
          ISHP = SHPVOI
        ELSE IF ( SHAPE .EQ. 'lor' ) THEN
          ISHP = SHPLOR
        ELSE IF ( SHAPE .EQ. 'dop' ) THEN
          ISHP = SHPDOP
        ELSE IF ( SHAPE .EQ. 'chi' ) THEN
          ISHP = SHPCHI
        ELSE IF ( SHAPE .EQ. 'vvw' .OR. SHAPE .EQ. 'van') THEN
          ISHP = SHPVVW
        ELSE
          FAIL = .TRUE.                                  
          ERRMSG = 'F-INPSHP: Unrecognised Line Shape: '//SHAPE
          RETURN
        END IF
C
C Subsequent fields in same record should contain list of gases with spec.shape
        IFLD = 2
        CALL TXTFLD ( 80, RECORD, IFLD, ISTA, IEND )
        DO WHILE ( ISTA .GT. 0 ) 
          CODE = RECORD(ISTA:IEND)
          CALL LOCASE ( CODE, CODE )
C First check if '*' is found (apply shape for all remaining gases)
          IF ( CODE .EQ. '*' ) THEN
            IF ( SHPDEF .EQ. 0 ) THEN
              SHPDEF = ISHP
              LOGMSG = 'I-INPSHP: Setting default '//SHAPE//
     &          ' lineshape for any remaining gases'
              CALL RFMLOG ( LOGMSG, FAIL, ERRMSG )
              IF ( FAIL ) RETURN	   
              GOTO 100
            ELSE
              ERRMSG = 'F-INPSHP: Two different line shapes '//
     &                 'specified default gases: *'
              FAIL = .TRUE.
              RETURN
            END IF
          END IF
          DO IGAS = 1, NGAS
            IF ( CODE .EQ. CODGAS(IGAS) ) THEN
              IF ( SHPGAS(IGAS) .EQ. SHPCTM ) THEN
                ERRMSG = 'F-INPSHP: line shape specified for '//
     &            CODE//' but already set to continuum-only'
              ELSE IF ( SHPGAS(IGAS) .EQ. SHPXSC ) THEN
                ERRMSG = 'F-INPSHP: line shape specified for '//
     &            CODE//' but already set as cross-section'
              ELSE IF ( SHPGAS(IGAS) .NE. 0 ) THEN
                ERRMSG = 'F-INPSHP: Two different line shapes '//
     &                   'specified for gas:'//CODE
              ELSE                          ! SHPGAS = 0, not yet defined
                SHPGAS(IGAS) = ISHP
                DO ISO = 1, NISGAS(IGAS)        ! Also set any other isotopes
                  JGAS = ISOGAS(ISO,IGAS)
                  SHPGAS(JGAS) = SHPGAS(IGAS)
                END DO
                LOGMSG = 'I-INPSHP: Using '//SHAPE
     &            //' lineshape for gas='//CODE
                CALL RFMLOG ( LOGMSG, FAIL, ERRMSG )
                IF ( FAIL ) RETURN
                GOTO 100
              END IF
              FAIL = .TRUE.
              RETURN                     ! Exit with fatal error
            END IF
          END DO
          ERRMSG = 
     &      'F-INPSHP: Line Shape specified for unlisted gas='//CODE
          FAIL = .TRUE.
          RETURN
 100      CONTINUE
          IFLD = IFLD + 1
          CALL TXTFLD ( 80, RECORD, IFLD, ISTA, IEND )
        END DO
        CALL NXTREC ( LUNDRV, RECORD, ENDSEC, FAIL, ERRMSG ) 
        IF ( FAIL ) RETURN
      END DO
C
C Reached end of section - check all gases have codes specified. 
      DO IGAS = 1, NGAS
        IF ( SHPGAS(IGAS) .EQ. 0 ) THEN         ! No shape specified
          IF ( SHPDEF .NE. 0 ) THEN             ! A default shape has been spec
            SHPGAS(IGAS) = SHPDEF               ! ..so use default shape
            LOGMSG = 'I-INPSHP: Using default lineshape for gas='
     &               //CODGAS(IGAS)
            CALL RFMLOG ( LOGMSG, FAIL, ERRMSG )
            IF ( FAIL ) RETURN
          ELSE                                   ! No default spec.either
            ERRMSG = 
     &        'F-INPSHP: No line shape specified for gas='//CODGAS(IGAS)
            RETURN
          END IF
        END IF
C
C Issue warning if inconsistent shape + continuum used for CO2 or H2O
        IF ( CTMGAS(IGAS) .AND. SHPGAS(IGAS) .NE. SHPCTM ) THEN
          IF ( IDXGAS(IGAS) .EQ. IDXCO2 ) THEN
            IF ( SHPGAS(IGAS) .NE. SHPCHI ) THEN
              LOGMSG = 'W-INPSHP: CHI-FACTOR recommended for'//
     &                 ' use with CO2 continuum'  
              CALL RFMLOG ( LOGMSG, FAIL, ERRMSG )
              IF ( FAIL ) RETURN
            END IF
          ELSE IF ( IDXGAS(IGAS) .EQ. IDXH2O ) THEN
            IF ( SHPGAS(IGAS) .NE. SHPLOR .AND. 
     &           SHPGAS(IGAS) .NE. SHPVOI       ) THEN
              LOGMSG = 'W-INPSHP: LORENTZ or VOIGT line-shape'//
     &                 ' recommended for use with H2O continuum'
              CALL RFMLOG ( LOGMSG, FAIL, ERRMSG )
              IF ( FAIL ) RETURN
            END IF
          END IF
        END IF
C
C Issue warning if VVW correction flag also being used with VVW lineshape
        IF ( VVWFLG .AND. SHPGAS(IGAS) .EQ. SHPVVW ) THEN
          LOGMSG = 'W-INPSHP: VVW correction will be ignored '//
     &             'for lines calculated with VVW shape'
          CALL RFMLOG ( LOGMSG, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        END IF
C
C Currently, Chi-factor only implemented for CO2
        IF ( SHPGAS(IGAS) .EQ. SHPCHI .AND. 
     &       IDXGAS(IGAS) .NE. IDXCO2 ) THEN
          FAIL = .TRUE. 
          ERRMSG = 'F-INPSHP: Chi-factor only appicable to CO2'
          RETURN
        END IF
      END DO
C
      FAIL = .FALSE.
      END
