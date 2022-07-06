      SUBROUTINE INPXSC ( LUNDRV, LUNNXT, FAIL, ERRMSG )
C
C VERSION
C     17OCT13 AD Simplified: Remove GETXSC arguments and local GOTXSC check
C     19SEP12 AD Change NAMXSC, NAMDEF from C*80 to C*200
C     23JUL03 AD Add XSCDEF
C     13JUN02 AD Allow for empty section 
C     17APR01 AD Also set LXSC=0
C     04DEC97 AD Add GOTXSC to check for duplicate section
C     03MAR97 AD Version 3.
C     26FEB97 AD Set GETXSC FALSE if any file found
C     16JAN97 AD Initialise GASXFL
C     01OCT96 AD Version 2.
C     01SEP96 AD Version 1.
C
C DESCRIPTION
C     Read list of Molecular Cross-section data files from RFM driver table.
C     Called by RFMINP once.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      LUNDRV  !  I  LUN for Driver File
      INTEGER      LUNNXT  ! I/O Next available LUN
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  NXTFLD ! Load next field from section of RFM driver file
     &, RFMLOG ! Write text message to RFM log file.
     &, XSCCHK ! Check reqd X/S species uniquely represented in X/Sec files
     &, XSCDEF ! Use default .xsc filename to find any missing files.
     &, XSCFIL ! Open X-Section data file read from RFM driver table.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'xflcom.inc' ! X-section data files
      INCLUDE 'xsccom.inc' ! X-section data
C
C LOCAL VARIABLES
      INTEGER       IGAS        ! Gas counter
      INTEGER       LENGTH      ! Length of field read from driver file
      CHARACTER*132 MESSGE      ! Text message sent to LOG file
      CHARACTER*200 NAMDEF      ! XSC filename template
      CHARACTER*200 NAMXSC      ! Name of XSC file
      LOGICAL       GOTDEF      ! True if filename template loaded
C
C DATA STATEMENTS
      DATA GOTDEF / .FALSE. /
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      NXFL = 0
      NXSC = 0
      LXSC = 0
      DO IGAS = 1, NGAS
        GASXFL(IGAS) = .FALSE.
      END DO
C
C Read first field in driver table
      CALL NXTFLD ( LUNDRV, NAMXSC, LENGTH, FAIL, ERRMSG ) 
      IF ( FAIL ) RETURN
      IF ( LENGTH .EQ. 0 ) THEN
        IF ( LUTFLG ) THEN
          MESSGE = 'W-INPXSC: No X/S files in *XSC section '//
     &             '- assume data supplied via LUTs'
          CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        ELSE
          FAIL = .TRUE.
          ERRMSG = 'F-INPXSC: No data in *XSC section'
        END IF
        IF ( FAIL ) RETURN
      END IF
C Jump here for each subsequent field in Driver Table
      DO WHILE ( LENGTH .NE. 0 ) 
        IF ( INDEX ( NAMXSC(1:LENGTH), '*' ) .NE. 0 ) THEN
          IF ( GOTDEF ) THEN
            ERRMSG = 'F-INPXSC: *XSC section contains more than '//
     &               'one filename template'
            FAIL = .TRUE.
            RETURN
          ELSE
            GOTDEF = .TRUE.
            NAMDEF = NAMXSC(1:LENGTH)
          ENDIF
        ELSE
          CALL XSCFIL ( LUNNXT, NAMXSC(1:LENGTH), FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        END IF
        CALL NXTFLD ( LUNDRV, NAMXSC, LENGTH, FAIL, ERRMSG ) 
        IF ( FAIL ) RETURN
      END DO
C
      IF ( GOTDEF ) THEN
        CALL XSCDEF ( LUNNXT, NAMDEF, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      ENDIF
C
      CALL XSCCHK ( FAIL, ERRMSG )               ! Check no overlaps
C
      END
