      SUBROUTINE REACYC ( FAIL, ERRMSG )
C
C VERSION
C     14-JAN-03  AD  Set EOF=FALSE each time
C     12-JUN-02  AD  Simplify: assume only one HITRAN file
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-OCT-96  AD  Save BSLQ in cyclic buffer
C     21-SEP-96  AD  Rename ILSVIB,IUSVIB to ILSNTE,IUSNTE
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Read HITRAN data into cyclic buffers.
C     Called by RFMFIN at start of each fine mesh calculation.
C
      IMPLICIT NONE
C
C  ARGUMENTS      
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG  !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  REAHIT ! Read HITRAN line data file for RFM
C
C GLOBAL CONSTANTS 
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'cyccom.inc' ! Cyclic line data buffers
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'hflcom.inc' ! HITRAN line data file
      INCLUDE 'hitcom.inc' ! HITRAN line database variables
C
C LOCAL VARIABLES
      INTEGER        ICYC  !  Index of line cyclic buffer
      INTEGER        ILIN  !  Counter for lines in Cyclic buffer
      LOGICAL        EOF   !  Set TRUE if end of HITRAN file reached
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Count the no. of useful lines already in buffer and label the first 
C ILIN ends up storing 1 + (no.lines which can be dropped)
C
      ILIN = 1
      ICYC = MOD ( ILIN + ICYC1 - 2, MAXCYC ) + 1  ! Lowest WNO line in buffer
      DO WHILE ( WNOCYC(ICYC) .LT. WNLFIN .AND. NLIN .GT. 0 )
        NLIN = NLIN - 1
        ILIN = ILIN + 1
        ICYC = MOD ( ILIN + ICYC1 - 2, MAXCYC ) + 1
      END DO
      ICYC1 = ICYC
C
C Note that to avoid back-stepping and rereading a line it is necessary to load
C one line beyond the required range into the cyclic buffer from each file
C
      EOF = .FALSE.
      DO WHILE ( WNOHFL .LT. WNUFIN .AND. .NOT. EOF ) 
        CALL REAHIT ( EOF, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        WNOHFL = WNUM             ! WNUM valid even for EOF record
        IF ( .NOT. EOF ) THEN
          IF ( NLIN .GE. MAXCYC ) THEN
            ERRMSG = 
     &        'F-REACYC: Cyclic line buffer full (MAXCYC too small)'
            FAIL = .TRUE.
            RETURN    
          ELSE                          ! Save new line data from /HITCOM/
            NLIN = NLIN + 1
            ICYC = MOD ( ICYC1 + NLIN - 2, MAXCYC ) + 1
            IDGCYC(ICYC) = IDGAS
            ISOCYC(ICYC) = ISO
            WNOCYC(ICYC) = WNUM
            STRCYC(ICYC) = STREN
            ABRCYC(ICYC) = ABROAD
            SBRCYC(ICYC) = SBROAD
            ELSCYC(ICYC) = ELS
            ABCCYC(ICYC) = ABCOEF
            TSPCYC(ICYC) = TSP
            IUSCYC(ICYC) = IUSGQ
            ILSCYC(ICYC) = ILSGQ
            IUVCYC(ICYC) = IUSNTE
            ILVCYC(ICYC) = ILSNTE
            BLQCYC(ICYC) = BSLQ
          END IF
        END IF
      END DO
C
      END
