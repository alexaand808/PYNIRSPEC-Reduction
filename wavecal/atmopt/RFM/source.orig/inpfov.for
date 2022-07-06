      SUBROUTINE INPFOV ( LUNDRV, LUNFOV, FAIL, ERRMSG )
C
C VERSION
C     16OCT13 AD Simplified. Remove GETFOV argument.
C     18MAR01 AD Remove GETOBS argument. Move FOVTAN, TANPTH to INPCHK
C     11AUG99 AD Remove redundant flgcom.inc
C     18JUL98 AD Warn if using Altitude-based FOV with OBS flag
C                    Add GETOBS argument
C                    Add check for any contents of *FOV section
C     14JUL98 AD Replace NXTREC with NXTFLD
C     04DEC97 AD Add GOTFOV to check for duplicate section
C     03MAR97 AD Version 3.
C     01OCT96 AD Version 2.
C     01SEP96 AD Version 1.
C
C DESCRIPTION
C     Read FOV data file from RFM driver table and load.
C     Called by RFMINP once if *FOV section header found in driver table
C     and FOVFLG=TRUE
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      LUNDRV  !  I  LUN for Driver File
      INTEGER      LUNFOV  ! I/O Next available LUN
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  ENDCHK ! Check end of Driver Table section has been reached.
     &, FOVFIL ! Open, read FOV data file and close
     &, NXTFLD ! Load next field from section of RFM driver file.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'fovcom.inc' ! FOV data
C
C LOCAL VARIABLES
      INTEGER       LENGTH      ! Length of file name
      CHARACTER*80  NAMFOV      ! Name of FOV file
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      NFOV = 0
C
C Get first field in this Driver Table section
C
      CALL NXTFLD ( LUNDRV, NAMFOV, LENGTH, FAIL, ERRMSG ) 
      IF ( FAIL ) RETURN
      IF ( LENGTH .EQ. 0 ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-INPFOV: No filename supplied in *FOV section'
        RETURN
      END IF

      CALL FOVFIL ( LUNFOV, NAMFOV(1:LENGTH), FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
C Check no more fields in this section (only one allowed, ie filename)
      CALL ENDCHK ( LUNDRV, 'FOV', FAIL, ERRMSG )
C
      END
