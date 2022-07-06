      SUBROUTINE INPTPS ( LUNDRV, LUNTPS, FAIL, ERRMSG )
C
C VERSION
C     17OCT13 AD Simplified: Remove GETTPS and local GOTTPS check
C     23JUL03 AD Add TPSDEF
C     12FEB03 AD Original.
C
C DESCRIPTION
C     Read list of tabulated TIPS data files
C     Called by RFMINP once if *TPS key listed in RFM Driver table.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      LUNDRV  !  I  LUN for Driver File
      INTEGER      LUNTPS  !  I  Temporary LUN for TPS file
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  NXTFLD ! Load next field from section of RFM driver file
     &, TPSDEF ! Use default .tps filename to find any missing files.
     &, TPSFIL ! Open TIPS data file, check and load contents.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'tpscom.inc' ! Tabulated TIPS data
C
C LOCAL VARIABLES
      INTEGER       IHLN   ! Counter for HITRAN line molecules
      INTEGER       IISO   ! Counter for isotopes of each molecule
      INTEGER       LENGTH ! Length of field read from driver file
      CHARACTER*200 NAMDEF ! .tps filename template
      CHARACTER*200 NAMTPS ! Name of TPS file
      LOGICAL       GOTDEF ! True if filename template loaded
C
C DATA STATEMENTS
      DATA GOTDEF / .FALSE. /
      SAVE GOTDEF
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Initialise variables in TPSCOM.INC
      NTPS = 0
      DO IHLN = 1, MAXHLN
        DO IISO = 1, MAXISO
          IDXTPS(IISO,IHLN) = 0
        END DO
      END DO
C 
      CALL NXTFLD ( LUNDRV, NAMTPS, LENGTH, FAIL, ERRMSG ) 
      IF ( FAIL ) RETURN
      IF ( LENGTH .EQ. 0 ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-INPTPS: No filenames in *TPS section'
        RETURN
      END IF
C Loop here for each field of Driver Table (= name of TPS file)
      DO WHILE ( LENGTH .NE. 0 ) 
        IF ( INDEX ( NAMTPS(1:LENGTH), '*' ) .NE. 0 ) THEN
          IF ( GOTDEF ) THEN
            ERRMSG = 'F-INPTPS: *TPS section contains more than '//
     &                 'one filename template'
            FAIL = .TRUE.
            RETURN
          ELSE
            GOTDEF = .TRUE.
            NAMDEF = NAMTPS(1:LENGTH)
          ENDIF
        ELSE
          CALL TPSFIL ( LUNTPS, NAMTPS(1:LENGTH), FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        ENDIF
        CALL NXTFLD ( LUNDRV, NAMTPS, LENGTH, FAIL, ERRMSG ) 
        IF ( FAIL ) RETURN
      END DO
C
C Load any remaining files using filename template
      IF ( GOTDEF ) CALL TPSDEF ( LUNTPS, NAMDEF, FAIL, ERRMSG )
C
      END
