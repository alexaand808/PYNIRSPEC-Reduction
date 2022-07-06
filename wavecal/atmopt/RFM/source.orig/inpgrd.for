      SUBROUTINE INPGRD ( LUNDRV, LUNNXT, FAIL, ERRMSG )
C
C VERSION
C     16OCT13 AD Simplified. Remove GETGRD argument & local GOTGRD
C     23JUL03 AD Add GRDDEF
C     20MAR01 AD Initialise RESGFL
C     14JUL98 AD Replace NXTREC by NXTFLD
C     18DEC97 AD Original.
C
C DESCRIPTION
C     Read list of GRD files of user-defined irregular grids
C     Called by RFMINP once if *GRD key listed in RFM Driver table and GRDFLG=T
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
     &  GRDDEF ! Use default .grd filename to find any missing files.
     &, NXTFLD ! Load next field from section of RFM driver file
     &, GRDFIL ! Open Look-Up Table file and check contents.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'gflcom.inc' ! Grid file data
C
C LOCAL VARIABLES
      INTEGER       LENGTH      ! Length of field read from driver file
      CHARACTER*80  NAMDEF      ! GRD filename template
      CHARACTER*80  NAMGRD      ! Name of GRD file
      LOGICAL       GOTDEF      ! True if filename template loaded
C
C DATA STATEMENTS
      DATA GOTDEF / .FALSE. /
      SAVE GOTDEF
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Initialise variables in GFLCOM.INC
      NGFL = 0
      RESGFL = 0.0D0
C
C Loop here for each field of Driver Table (= name of GRD file)
      CALL NXTFLD ( LUNDRV, NAMGRD, LENGTH, FAIL, ERRMSG ) 
      IF ( FAIL ) RETURN
      DO WHILE ( LENGTH .NE. 0 ) 
        IF ( INDEX ( NAMGRD(1:LENGTH), '*' ) .NE. 0 ) THEN
          IF ( GOTDEF ) THEN
            ERRMSG = 'F-INPGRD: *GRD section contains more than '//
     &                 'one filename template'
            FAIL = .TRUE.
          ELSE
            GOTDEF = .TRUE.
            NAMDEF = NAMGRD(1:LENGTH)
          ENDIF
        ELSE
          CALL GRDFIL ( LUNNXT, NAMGRD(1:LENGTH), FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        END IF
        CALL NXTFLD ( LUNDRV, NAMGRD, LENGTH, FAIL, ERRMSG ) 
        IF ( FAIL ) RETURN
      END DO
C
C Look for any spectral ranges without assigned GRD files
      IF ( GOTDEF ) CALL GRDDEF ( LUNNXT, NAMDEF, FAIL, ERRMSG )
C
      END
