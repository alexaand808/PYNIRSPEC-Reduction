      SUBROUTINE INPTAN ( LUNDRV, LUNTAN, KEYCHK, FAIL, ERRMSG )
C
C VERSION
C     16OCT13 AD Rewritten, use NXTFFL and remove TANFIL
C                Remove NOTTAN argument for TANCHK. 
C                Remove MAXSTR parameter, add ACCLOG
C     08JUN01 AD Check at least 1 tan.ht if using LOS flag
C     18MAR01 AD Move HOMPTH, FLXPTH, NADPTH, TANPTH, TANHGT to INPCHK
C     02MAY00 AD Add FLXPTH and tests for flux calculations
C     31JUL99 AD Split info message listing tan.hgts. 
C     02JAN99 AD Initialise LTAN=NTAN
C     18JUL98 AD Add KEYCHK argument and test alternative inputs
C     10JUL98 AD Check for at least one entry in *TAN section
C     22OCT97 AD Add TABPTH and handling of TAB option
C     03MAR97 AD Version 3.
C     14JAN97 AD Add extra routine NADPTH
C     01OCT96 AD Version 2.
C     01SEP96 AD Version 1.
C
C DESCRIPTION
C     Read RFM tangent points from driver table following *TAN marker
C     Called by RFMINP once.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      LUNDRV  !  I  LUN for Driver File
      INTEGER      LUNTAN  !  I  LUN for tangent point file
      CHARACTER*4  KEYCHK  !  I  Key actually used for *TAN section
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  ACCLOG ! Accumulate and send text record to log file.
     &, NXTFFL ! Load next field from section of RFM drv, expanding filenames
     &, RFMLOG ! Write text message to RFM log file.
     &, TABPTH ! Set up table-entry paths for RFM calculations.
     &, TANCHK ! Read tangent height from string and add to /TANCOM/
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL VARIABLES
      INTEGER      LENGTH  ! Length of field read from driver table
      CHARACTER*80 FIELD   ! Tangent height or Tan.Hgt filename
      CHARACTER*80 WRNMSG  ! Text message sent to LOG file       
C
C EXECUTABLE CODE --------------------------------------------------------------
C
C Check if valid Key has been used to describe *TAN section (assuming that
C it has already been tested to be either *dim, *ele, *geo, *len, *sec or *tan)
C
      FAIL = .TRUE.
      IF ( KEYCHK .EQ. '*dim' .AND. .NOT. TABFLG ) THEN
        ERRMSG = 'F-INPTAN: *DIM section only valid with TAB flag'
      ELSE IF ( KEYCHK .EQ. '*ele' .AND. .NOT. OBSFLG ) THEN
        ERRMSG = 'F-INPTAN: *ELE section only valid with OBS flag'
      ELSE IF ( KEYCHK .EQ. '*len' .AND. .NOT. HOMFLG ) THEN
        ERRMSG = 'F-INPTAN: *LEN section only valid with HOM flag'
      ELSE IF ( KEYCHK .EQ. '*lev' .AND. .NOT. FLXFLG ) THEN
        ERRMSG = 'F-INPTAN: *LEV section only valid with FLX flag'
      ELSE IF ( KEYCHK .EQ. '*sec' .AND. FLXFLG ) THEN
        ERRMSG = 'F-INPTAN: *SEC section not valid FLX flag'
      ELSE IF ( KEYCHK .EQ. '*sec' .AND. .NOT. 
     &                  ( NADFLG .OR. ZENFLG ) ) THEN
        ERRMSG = 
     &    'F-INPTAN: *SEC section only valid with NAD or ZEN flags'
      ELSE IF ( KEYCHK .EQ. '*geo' .AND. ( ZENFLG .OR. NADFLG .OR.
     &          TABFLG .OR. HOMFLG ) ) THEN
        ERRMSG = 
     &    'F-INPTAN: *GEO not valid with ZEN, NAD, HOM or TAB flags' 
      ELSE
        FAIL = .FALSE.
      END IF
      IF ( FAIL ) RETURN
      USRELE = KEYCHK .EQ. '*ele'
      USRGEO = KEYCHK .EQ. '*geo'
C
      NTAN = 0
C TAB option has its own listing of parameters to log file in CHKTAB
      IF ( .NOT. TABFLG ) THEN
        CALL ACCLOG ( 'I-INPTAN: ', FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
      CALL NXTFFL ( LUNDRV, LUNTAN, .FALSE., 
     &              FIELD, LENGTH, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      DO WHILE ( LENGTH .NE. 0 ) 
        CALL TANCHK ( FIELD, FAIL, ERRMSG ) 
        IF ( FAIL ) RETURN
        IF ( .NOT. TABFLG ) THEN
          CALL ACCLOG ( FIELD, FAIL, ERRMSG ) 
          IF ( FAIL ) RETURN
        END IF
        CALL NXTFFL ( LUNDRV, LUNTAN, .FALSE., 
     &                FIELD, LENGTH, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END DO
C
C Special case for TAB option - check all 6 parameters supplied, set up paths 
C and exit
      IF ( TABFLG ) THEN
        IF ( NTAN .NE. 6 ) THEN
          WRITE ( ERRMSG, '(A,I1)' ) 
     &      'F-INPTAN: Not enough tabulation parameters. '//
     &      '6 required, found ', NTAN
          FAIL = .TRUE.
          RETURN
        ELSE
          NTAN = 0
          CALL TABPTH ( FAIL, ERRMSG )
          RETURN
        END IF
      ELSE
        CALL ACCLOG ( ' ', FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
C
C Check at least one tangent height supplied
      IF ( NTAN .EQ. 0 ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-INPTAN: No entries in '//KEYCHK//' section'
        RETURN
      END IF
C
C For LOS flag (line-of-sight pointing Jacobians at least 2 and preferably 3
C tangent paths are required (3 for quadratic fit)
      IF ( LOSFLG ) THEN
        IF ( NTAN .EQ. 1 ) THEN
          FAIL = .TRUE.
          ERRMSG = 'F-INPTAN: LOS flag requires at least 2 tangent hts'
          RETURN
        ELSE IF ( 2 * NTAN .GT. MAXTAN ) THEN
          FAIL = .TRUE.
          WRITE ( ERRMSG, '(A,I9,A,I9,A)' )
     &      'F-INPTAN: LOS flag requires', 2 * NTAN, 
     &      'tan. paths, >MAXTAN=', MAXTAN,' in rfmsiz.inc'
          RETURN
        ELSE IF ( NTAN .EQ. 2 ) THEN
          WRNMSG = 'W-INPTAN: Only 2 tan.hts, '//
     &             'so LOS Jacobians from linear interpolation'
          CALL RFMLOG ( WRNMSG, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        END IF
        LOSTAN = NTAN
      ELSE
        LOSTAN = 0
      END IF          
C
      MTAN = NTAN
      LTAN = NTAN
C
      END
