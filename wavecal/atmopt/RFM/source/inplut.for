      SUBROUTINE INPLUT ( LUNDRV, LUNNXT, FAIL, ERRMSG )
C
C VERSION
C     16OCT13 AD Simplified. Remove GETLUT argument and local GOTLUT test
C     07AUG-13 AD Remove redundant local variable ILFL
C                    Remove redundant GRDDEF from EXTERNAL list
C     23JUL03 AD Add LUTDEF
C     12FEB03 AD Removed redundant variables and .inc files
C     13JUN02 AD Remove GETHIT,GETXSC arguments
C     14JUL98 AD Replace NXTREC by NXTFLD
C     10JUL98 AD Add GETHIT, GETXSC arguments and modify if required
C     04DEC97 AD Add GOTLUT to check for duplicate section
C     02JUL97 AD Original.
C
C DESCRIPTION
C     Read list of Look-Up tables of absorption coefficients.
C     Called by RFMINP once if *LUT key listed in RFM Driver table and LUTFLG=T
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
     &  LUTDEF ! Use default LUT filename to find any missing files.
     &, LUTFIL ! Open Look-Up Table file and check contents.
     &, NXTFLD ! Load next field from section of RFM driver file
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'lflcom.inc' ! Look-Up Table file data
C
C LOCAL VARIABLES
      INTEGER       LENGTH      ! Length of field read from driver file
      CHARACTER*80  NAMDEF      ! LUT filename template
      CHARACTER*80  NAMLUT      ! Name of LUT file
      LOGICAL       GOTDEF      ! True if filename template loaded
C
C DATA STATEMENTS
      DATA GOTDEF / .FALSE. /
      SAVE GOTDEF
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Initialise variables in LFLCOM.INC
      NLFL = 0
C
      CALL NXTFLD ( LUNDRV, NAMLUT, LENGTH, FAIL, ERRMSG ) 
      IF ( FAIL ) RETURN
C Loop here for each field of Driver Table (= name of LUT file)
      DO WHILE ( LENGTH .NE. 0 )
        IF ( INDEX ( NAMLUT(1:LENGTH), '*' ) .NE. 0 ) THEN
          IF ( GOTDEF ) THEN
            ERRMSG = 'F-INPLUT: *LUT section contains more than '//
     &                 'one filename template'
            FAIL = .TRUE.
          ELSE
            GOTDEF = .TRUE.
            NAMDEF = NAMLUT(1:LENGTH)
          ENDIF
        ELSE
          CALL LUTFIL ( LUNNXT, NAMLUT(1:LENGTH), FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        END IF
        CALL NXTFLD ( LUNDRV, NAMLUT, LENGTH, FAIL, ERRMSG ) 
        IF ( FAIL ) RETURN
      END DO
C
      IF ( GOTDEF ) CALL LUTDEF ( LUNNXT, NAMDEF, FAIL, ERRMSG )
C
      END
