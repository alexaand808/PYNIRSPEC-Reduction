      SUBROUTINE SKPSEC ( LUNDRV, LOGMSG, FAIL, ERRMSG )
C
C VERSION
C     10OCT13 AD Original.
C
C DESCRIPTION
C     Skip driver table section.
C     Called multiply by RFMINP.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER       LUNDRV !  I  LUN for Driver File
      CHARACTER*(*) LOGMSG !  I  Message to be printed to RFM log file 
      LOGICAL       FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     & NXTREC ! Load next record from RFM input file.
     &, RFMLOG ! Write text message to RFM log file.
C
C LOCAL VARIABLES
      CHARACTER*80 RECORD  ! Dummy record from Driver file 
      LOGICAL      ENDSEC  ! T=End of driver table section reached.
C
C EXECUTABLE CODE ------------------------------------------------------------
C
      IF ( LOGMSG .NE. ' ' ) THEN
        CALL RFMLOG ( LOGMSG, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF

      ENDSEC = .FALSE.
      DO WHILE ( .NOT. ENDSEC )
        CALL NXTREC ( LUNDRV, RECORD, ENDSEC, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END DO
C
      END
