      SUBROUTINE QTWARN ( IDGAS, ISO, TEM )
C
C VERSION
C     01OCT13 AD Original.
C
C DESCRIPTION
C     Handle warning messages from subroutine QTFCT
C     Called multiply by QTFCT
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER IDGAS !  I  HITRAN gas ID
      INTEGER ISO   !  I  HITRAN isotope ID, or -1 = no TIPS for molecule
      REAL    TEM   !  I  Path temperature [K] (if IDGAS = -1 )
C
      EXTERNAL
     &  RFMLOG ! Write text message to RFM log file. 
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C LOCAL CONSTANTS
      INTEGER MAXWRN
        PARAMETER ( MAXWRN = 10 ) ! Max no of different warning messages
C
C COMMON VARIABLES
      INCLUDE 'gascom.inc' ! Molecular and isotope data
C
C LOCAL VARIABLES
      INTEGER      IGAS  ! Index of molecule in *GAS arrays
      INTEGER      IWRN  ! Warning counter
      INTEGER      NWRN  ! No. of different warning messages issued so far
      INTEGER      IDGWRN(MAXWRN) ! List of IDGAS already issued as warnings
      INTEGER      ISOWRN(MAXWRN) ! List of ISO already issued as warnings
      LOGICAL      FAIL   ! T=fatal error detected
      LOGICAL      WARNHI ! T=issued warning if above TIPS high temp limit 
      LOGICAL      WARNLO ! T=issued warning if below TIPS low temp limit
      CHARACTER*80 ERRMSG ! Error message written if FAIL is TRUE
      CHARACTER*80 WRNMSG ! Text for warning message
C
C DATA STATEMENTS
      DATA NWRN / 0 /
      DATA WARNLO / .FALSE. /
      DATA WARNHI / .FALSE. /
      SAVE NWRN, WARNLO, WARNHI
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C First check for warning messages due to low or high temperature limits
      IF ( TEM .LT. 60.0 ) THEN
        IF ( .NOT. WARNLO ) THEN
          WRITE ( WRNMSG, '(A,1PG10.3)' ) 
     &      'W-QTFCT: Extrapolating TIPS below 60K, T=', TEM
          WARNLO = .TRUE.
          CALL RFMLOG ( WRNMSG, FAIL, ERRMSG )
          IF ( FAIL ) GOTO 900
        END IF
        RETURN
      ELSE IF ( TEM .GT. 3010.0 ) THEN
        IF ( .NOT. WARNHI ) THEN
          WRITE ( WRNMSG, '(A,1PG10.3)' ) 
     &      'W-QTFCT: Extrapolating TIPS above 3010K, T=', TEM
          WARNHI = .TRUE.
          CALL RFMLOG ( WRNMSG, FAIL, ERRMSG )
          IF ( FAIL ) GOTO 900
        END IF
        RETURN
      END IF
C
C For other warnings, check if maximum number of warnings already issues
      IF ( NWRN .GE. MAXWRN ) RETURN                 ! reached max warnings
C
C Check if warning has been issues for this particular molec,iso combination
      DO IWRN = 1, NWRN
        IF ( IDGWRN(IWRN) .EQ. IDGAS .AND.
     &       ISOWRN(IWRN) .EQ. ISO        ) RETURN ! already warned of this
      END DO
C
C Construct warning message
      IGAS = IGSMOL(IDGAS)
      IF ( ISO .EQ. -1 ) THEN         
        WRITE ( WRNMSG, '(A,A,A)' ) 
     &    'W-QTFCT: No TIPS data for ', CODGAS(IGAS),
     &    ' - assuming T^(3/2) dependence'
      ELSE 
        WRITE ( WRNMSG, '(A,A,A,I2,A)' )    
     &    'W-QTFCT: No TIPS data for ', CODGAS(IGAS), ' Isotope#', ISO, 
     &    ' - assuming TIPS ~ Isotope#1'
      END IF
      CALL RFMLOG ( WRNMSG, FAIL, ERRMSG )
      IF ( FAIL ) GOTO 900
C
C Save details to avoid repetition
      NWRN = NWRN + 1
      IDGWRN(NWRN) = IDGAS
      ISOWRN(NWRN) = ISO
      IF ( NWRN .EQ. MAXWRN ) 
     &  WRNMSG = 'W-QTWARN: Further such warnings suppressed'
      CALL RFMLOG ( WRNMSG, FAIL, ERRMSG )
      IF ( FAIL ) GOTO 900
C
      RETURN
C
C Fatal error from RFMLOG
 900  CONTINUE
      WRITE ( *, * ) 'F-QTWARN: Fatal error from RFMLOG'
      WRITE ( *, * ) ERRMSG
      STOP
C
      END

