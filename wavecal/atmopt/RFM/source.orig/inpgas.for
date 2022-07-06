      SUBROUTINE INPGAS ( LUNDRV, LUNGAS, GETHIT, GETXSC, FAIL, ERRMSG )
C
C VERSION
C     13-JUN-02  AD  Use SHPGAS rather than IDXGAS to test line/XSC species
C     28-MAR-01  AD  Add IAXGAS
C     21-MAR-01  AD  Set CHGGAS = TRUE. 
C     20-DEC-98  AD  Modified to include GASALL
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Read RFM absorbing gases from driver table following *GAS marker
C     Called by RFMINP once.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      LUNDRV  !  I  LUN for Driver File
      INTEGER      LUNGAS  !  I  LUN for Gas file
      LOGICAL      GETHIT  !  O  Set TRUE if HITRAN line data required
      LOGICAL      GETXSC  !  O  Set TRUE if X-Section data required
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  GASALL ! Add wildcard gases to list of absorbers
     &, GASCHK ! Check Gas name and set indices.
     &, GASFIL ! Open file and load pretabulated gases.
     &, NXTFLD ! Load next field from section of RFM driver file
     &, RFMLOG ! Write text message to RFM log file.
     &, TXTFLD ! Identify start and end points of text field in record
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'idxcon.inc' ! HITRAN molecular indices
      INCLUDE 'shpcon.inc' ! Line Shape codes
C
C LOCAL CONSTANTS
      INTEGER       MAXSTR        ! Max size of text string listing all gases
        PARAMETER ( MAXSTR = 22 + 8 * MAXGAS ) 
C
C COMMON VARIABLES
      INCLUDE 'gascom.inc' ! Molecular and isotope data
C
C LOCAL VARIABLES
      INTEGER      I1,I2          ! Pointers to start,end of name in CODGAS
      INTEGER      IEND           ! End location for text identifying gas
      INTEGER      IGAS           ! Gas counter
      INTEGER      IMOL           ! HITRAN index of molecule
      INTEGER      ISTA           ! Starting location for text identifying gas
      INTEGER      LENGTH         ! Length of field read from Driver file
      LOGICAL      NOTGAS         ! Set TRUE if gas not recognised
      CHARACTER*(MAXSTR) MESSGE   ! Text message sent to LOG file
      CHARACTER*80 NAMGAS         ! Either name of gas or file containing gases
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      NGAS = 0 
      IAXGAS = 0
      DO IMOL = 1, MAXMOL
        IGSMOL(IMOL) = 0
      END DO
C
C Loop here for field in *GAS section of driver table 
C
  100 CALL NXTFLD ( LUNDRV, NAMGAS, LENGTH, FAIL, ERRMSG ) 
      IF ( FAIL ) RETURN
      IF ( LENGTH .NE. 0 ) THEN 
        CALL GASCHK ( NAMGAS(1:LENGTH), NOTGAS, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
C
C If field was not recognised as a molecule name, it may be a filename
C
        IF ( NOTGAS ) THEN
          IF ( NAMGAS(1:1) .EQ. '*' ) THEN
            CALL GASALL ( NAMGAS(1:LENGTH), .FALSE., FAIL, ERRMSG )
          ELSE
            CALL GASFIL ( LUNGAS, NAMGAS(1:LENGTH), FAIL, ERRMSG )
          END IF
          IF ( FAIL ) RETURN
        END IF
        GOTO 100
      END IF
C
C Add any gases specified by wildcard (no effect if none specified)
C
      CALL GASALL ( ' ', .TRUE., FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
C Construct message for log file listing all identified species
C
      GETHIT = .FALSE.
      GETXSC = .FALSE.
      MESSGE = 'I-INPGAS: Using gases:'
      ISTA = 23 
      IF ( NGAS .EQ. 0 ) THEN
        MESSGE = 'W-INPGAS: No gases specified'
      ELSE
        DO IGAS = 1, NGAS
          CALL TXTFLD ( 7, CODGAS(IGAS), 1, I1, I2 )
          IEND = ISTA + I2 - I1 + 1 
          MESSGE(ISTA:IEND) = ' '//CODGAS(IGAS)(I1:I2)
          ISTA = IEND + 1
          IF ( SHPGAS(IGAS) .EQ. SHPXSC ) THEN
            GETXSC = .TRUE.
          ELSE
            GETHIT = .TRUE.
          END IF
          CHGGAS(IGAS) = .TRUE.                   ! Flag for new profile loaded
          IF ( IDXGAS(IGAS) .EQ. IDXAER ) IAXGAS = IGAS
        END DO
      END IF
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )        ! Output all gases so far
C
      END
