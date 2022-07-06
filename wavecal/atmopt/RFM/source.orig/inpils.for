      SUBROUTINE INPILS ( LUNDRV, LUNILS, FAIL, ERRMSG )
C
C VERSION
C      16OCT13 AD Simplified. Remove GETILS argument and local GOTILS check.
C      18MAR01 AD Remove GETFIN argument, move ILSCHK to INPCHK
C      15MAR01 AD Add GETFIN argument to test if FIN data already loaded
C      04DEC97 AD Add GOTILS to check for duplicate section
C      03MAR97 AD Version 3.
C      01OCT96 AD Version 2.
C      01SEP96 AD Version 1.
C
C DESCRIPTION
C     Read list of ILS data files from RFM driver table.
C     Called by RFMINP once if *ILS section header found in driver table and
C     ILSFLG=TRUE
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      LUNDRV  !  I  LUN for Driver File
      INTEGER      LUNILS  !  I  Temporary LUN for reading ILS files
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  ILSFIL ! Open, read ILS data file and close
     &, NXTFLD ! Load next field from section of RFM driver file
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'ilscom.inc' ! ILS data
C
C LOCAL VARIABLES
      INTEGER       LENGTH ! Length of name of ILS file
      CHARACTER*120 NAMILS ! Name of ILS file
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      NILS = 0
C
      CALL NXTFLD ( LUNDRV, NAMILS, LENGTH, FAIL, ERRMSG ) 
      IF ( FAIL ) RETURN
      DO WHILE ( LENGTH .NE. 0 )
        CALL ILSFIL ( LUNILS, NAMILS(1:LENGTH), FAIL, ERRMSG )   
        IF ( FAIL ) RETURN
        CALL NXTFLD ( LUNDRV, NAMILS, LENGTH, FAIL, ERRMSG ) 
        IF ( FAIL ) RETURN
      END DO
C
      IF ( NILS .EQ. 0 ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-INPILS: No (suitable) ILS data provided'
      END IF
C
      END
