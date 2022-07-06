      SUBROUTINE RFMOPN ( LUNNXT, ISPC, FAIL, ERRMSG )
C
C VERSION
C     20-SEP-12  AD  Change NAMTMP from C*80 to C*200
C     14-SEP-11  AD  Add RJT flag
C     09-AUG-03  AD  Add BBT flag
C     20-JUL-02  AD  Set file counters for LUN flag.
C     20-APR-00  AD  Add COOling rate files
C     01-OCT-99  AD  Use local variable NAMTMP as argument to subroutines 
C     02-FEB-99  AD  Add RFMSIZ.INC
C     22-OCT-97  AD  Add OPNTAB
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Open RFM output files and write headers
C     Called by RFM for each spectral range.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      LUNNXT  ! I/O Next available LUN
      INTEGER      ISPC    !  I  Spectral range number
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error occurs
      CHARACTER*80 ERRMSG  !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  OPNOUT ! Open spectral output files for current spectral range
     &, OPNTAB ! Open TAB files for current spectral range
     &, OPNWID ! Open WID files for current spectral range
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'outcom.inc' ! RFM output file data
C
C LOCAL VARIABLES
      INTEGER      NFLOUT  ! No.of output files (LUN flag option only)
      CHARACTER*200 NAMTMP ! Temporary file name
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      NFLOUT = 0
C
C If Absorption data required, open output ABS files and write headers
C
      IF ( ABSFLG ) THEN
        LUNABS = LUNNXT
        NAMTMP = NAMABS
        IFLABS = NFLOUT
        CALL OPNOUT ( LUNNXT, ISPC, 'Absorption', NAMTMP, NFLOUT,
     &                FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        NAMABS = NAMTMP
      END IF
C
C If Bright.Temp data required, open output BBT files and write headers
C
      IF ( BBTFLG ) THEN
        LUNBBT = LUNNXT 
        NAMTMP = NAMBBT
        IFLBBT = NFLOUT
        CALL OPNOUT ( LUNNXT, ISPC, 'Bright. Temp.', NAMTMP, NFLOUT,
     &                FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        NAMBBT = NAMTMP
      END IF
C
C If Cooling Rate data required, open output COO files and write headers
C
      IF ( COOFLG ) THEN
        LUNCOO = LUNNXT
        NAMTMP = NAMCOO
        IFLCOO = NFLOUT
        CALL OPNOUT ( LUNNXT, ISPC, 'Cooling Rate', NAMTMP, NFLOUT,
     &                FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        NAMCOO = NAMTMP
      END IF
C
C If Optical Depth data required, open output OPT files and write headers
C
      IF ( OPTFLG ) THEN
        LUNOPT = LUNNXT
        NAMTMP = NAMOPT
        IFLOPT = NFLOUT
        CALL OPNOUT ( LUNNXT, ISPC, 'Optical Depth', NAMTMP, NFLOUT,
     &                FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        NAMOPT = NAMTMP
      END IF
C
C If Radiance data required, open output RAD files and write headers
C
      IF ( RADFLG ) THEN
        LUNRAD = LUNNXT 
        NAMTMP = NAMRAD
        IFLRAD = NFLOUT
        CALL OPNOUT ( LUNNXT, ISPC, 'Radiance', NAMTMP, NFLOUT,
     &                FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        NAMRAD = NAMTMP
      END IF
C
C If Rayl-J.Temp data required, open output RJT files and write headers
C 
      IF ( RJTFLG ) THEN
        LUNRJT = LUNNXT 
        NAMTMP = NAMRJT
        IFLRJT = NFLOUT
        CALL OPNOUT ( LUNNXT, ISPC, 'Rayl-J.Temp.', NAMTMP, NFLOUT,
     &                FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        NAMRJT = NAMTMP
      END IF
C
C If Transmission data required, open output TRA files and write headers
C
      IF ( TRAFLG ) THEN
        LUNTRA = LUNNXT
        NAMTMP = NAMTRA
        IFLTRA = NFLOUT
        CALL OPNOUT ( LUNNXT, ISPC, 'Transmission', NAMTMP, NFLOUT,
     &                FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        NAMTRA = NAMTMP
      END IF
C
C If Tabulated absorption coeffs reqd, open output TAB files and write headers
C
      IF ( TABFLG ) THEN
        CALL OPNTAB ( LUNNXT, ISPC, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
C
C If Widemesh diagnostics required, open output WID files and write headers
C
      IF ( WIDFLG ) THEN
        CALL OPNWID ( LUNNXT, ISPC, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
C
      END
