      SUBROUTINE INPHIT ( LUNDRV, LUNNXT, FAIL, ERRMSG )
C
C VERSION
C     16OCT13 AD Simplified. Remove GETHIT argument and local GOTHIT test
C     19SEP12 AD Change NAMHIT from C*80 to C*200
C     13JUN02 AD Initialise LUNHFL
C     12JUN02 AD Use only single HITRAN file
C     04DEC97 AD Add GOTHIT to check for duplicate section
C     03MAR97 AD Version 3.
C     26FEB97 AD Set GETHIT=FALSE if any HIT file found
C     01OCT96 AD Version 2.
C     01SEP96 AD Version 1.
C
C DESCRIPTION
C     Read HITRAN line data filename from RFM driver table.
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
     &  ENDCHK ! Check end of Driver Table section has been reached.
     &, HITFIL ! Open HITRAN line data file read from Driver Table
     &, NXTFLD ! Load next field from section of RFM driver file
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'hflcom.inc' ! HITRAN line data file
C
C LOCAL VARIABLES
      INTEGER       LENGTH ! Length of field in Driver Table
      CHARACTER*200 NAMHIT ! Name of HITRAN file 
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      LUNHFL = 0           ! Initialise = no actual data file (yet)
C
      CALL NXTFLD ( LUNDRV, NAMHIT, LENGTH, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
      IF ( LENGTH .EQ. 0 ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-INPHIT: No filename supplied in *HIT section'
        RETURN 
      END IF
C
C Read HITRAN filename
      CALL HITFIL ( LUNNXT, NAMHIT(1:LENGTH), FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      LUNNXT = LUNNXT + 1
C
C Check no more fields in this section (only one allowed, ie filename)
      CALL ENDCHK ( LUNDRV, 'HIT', FAIL, ERRMSG )
C
      END
