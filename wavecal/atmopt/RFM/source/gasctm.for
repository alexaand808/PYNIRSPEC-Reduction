      SUBROUTINE GASCTM ( CTMFLG, ANYQAL, LQS, QALSTR, FAIL, ERRMSG )
C
C VERSION
C     01-DEC-97  AD  Original
C
C DESCRIPTION
C     Initialise for Continuum data
C     Called by GASCHK for any gas with a qualifier string.
C     Also checks for (ctm) or (noctm) qualifiers appended to gas name
C     If (ctm) qualifier is present (=continuum only)
C       Check CTM flag is set (else fatal error)
C       Check CTMGAS flag is set for this molecule (else fatal error)
C       Set SHPGAS=SHPCTM (flag for continuum only shape)
C       Check no other qualifiers present (else fatal error)
C     Else if (noctm) qualifier present (=exclude continuum)
C       Check CTM flag is set (else warning)
C       Check CTMGAS is set for this molecule (else warning)
C       Set CTMGAS=.FALSE.
C       Remove this qualifier from QALSTR
C       Check if other qualifiers present, else set ANYQAL=.FALSE.
C     End if
C
      IMPLICIT NONE
C
C ARGUMENTS
      LOGICAL       CTMFLG  !  I  TRUE means enable continuum data
      LOGICAL       ANYQAL  ! I/O (Further) qualifiers present 
      INTEGER       LQS     ! I/O Length of QALSTR (1:LQS)
      CHARACTER*(*) QALSTR  ! I/O String containing list of qualifiers
      LOGICAL       FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  LOCASE ! Convert text string to lower case.
     &, RFMLOG ! Write text message to RFM log file.
C
C GLOBAL CONSTANTS
      INCLUDE 'idxcon.inc' ! HITRAN molecular indices
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'shpcon.inc' ! Line-shape codes
C
C COMMON VARIABLES
      INCLUDE 'gascom.inc' ! Molecular and isotope data
C
C LOCAL VARIABLES
      INTEGER   IPT  ! Pointer to '(' element of current qualifier 
      CHARACTER*80 LOGMSG ! Message sent to RFM log file
C
C EXECUTABLE CODE --------------------------------------------------------------
C
C Initialise continuum flags for those species with continuum data
C
      FAIL = .FALSE.
      IF ( IDXGAS(NGAS) .EQ. IDXH2O .OR.
     &     IDXGAS(NGAS) .EQ. IDXCO2 .OR.
     &     IDXGAS(NGAS) .EQ. IDXO2  .OR.
     &     IDXGAS(NGAS) .EQ. IDXN2       )  THEN
        CTMGAS(NGAS) = CTMFLG
      ELSE
        CTMGAS(NGAS) = .FALSE.
      END IF
      IF ( .NOT. ANYQAL ) RETURN
C
      CALL LOCASE ( QALSTR, QALSTR )     
C
C First check for presence of (ctm) qualifier anywhere in QALSTR
C
      IPT = INDEX ( QALSTR, '(ctm)' )
      IF ( IPT .NE. 0 ) THEN
        IF ( .NOT. CTMFLG ) THEN
          FAIL = .TRUE.
          ERRMSG = 'F-GASCTM: (ctm) qualifier set for gas='//
     &      CODGAS(NGAS)//' but CTM flag not enabled'
        ELSE IF ( .NOT. CTMGAS(NGAS) ) THEN
          FAIL = .TRUE. 
          ERRMSG = 'F-GASCTM: No continuum data available for gas='
     &      //CODGAS(NGAS)
        ELSE IF ( LQS .GT. 5 ) THEN
          FAIL = .TRUE.
          ERRMSG = 'F-GASCTM: Additional qualifiers present '//
     &        'as well as (ctm) for gas='//CODGAS(NGAS)
        ELSE
          ANYQAL = .FALSE.
          SHPGAS(NGAS) = SHPCTM
          LOGMSG = 'I-GASCTM: setting gas='//CODGAS(NGAS)//
     &        ' to continuum-only calculation'
          CALL RFMLOG ( LOGMSG, FAIL, ERRMSG )
        END IF
        RETURN                     ! exit for (ctm) qualifier
      END IF
C
C If (ctm) not present, check for (noctm) qualifier
C
      IPT = INDEX ( QALSTR, '(noctm)' ) 
      IF ( IPT .NE. 0 ) THEN
        IF ( .NOT. CTMFLG ) THEN
          LOGMSG = 'W-GASCTM: CTM Flag not set so (noctm) qualifier'//
     &      ' for gas='//CODGAS(NGAS)//' is redundant'
          CALL RFMLOG ( LOGMSG, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        END IF
        IF ( .NOT. CTMGAS(NGAS) ) THEN        
          LOGMSG = 'W-GASCTM: No continuum data for gas='//
     &      CODGAS(NGAS)//' so (noctm) qualifier is redundant'
          CALL RFMLOG ( LOGMSG, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        ELSE
          CTMGAS(NGAS) = .FALSE.
          LOGMSG = 'I-GASCTM: Excluding continuum for gas='
     &      //CODGAS(NGAS)
          CALL RFMLOG ( LOGMSG, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        END IF
C
C Check for any remaining qualifiers (isotopes and/or bands)
C 
        IF ( IPT .EQ. 1 ) THEN        ! no qualifiers before '(noctm)'...
          IF ( LQS .EQ. 7 ) THEN      ! ... and no qualifiers after '(noctm)'
            ANYQAL = .FALSE.
          ELSE                        ! ... and other qualifiers after '(noctm)'
            QALSTR = QALSTR(IPT+7:)   
          END IF
        ELSE                          ! other qualifiers before '(noctm)' ...
          IF ( LQS .EQ. IPT+6 ) THEN  ! ... and no qualifiers after '(noctm)'
            QALSTR = QALSTR(1:IPT-1)
          ELSE                        ! ... and other qualifiers after '(noctm)'
            QALSTR = QALSTR(1:IPT-1)//QALSTR(IPT+7:)
          END IF
        END IF
        LQS = LQS - 7
      END IF
C
C Exit for successful handling of (noctm) qualifier, or for neither (ctm) nor
C (noctm) qualifiers present
C
      FAIL = .FALSE.     
      END 
