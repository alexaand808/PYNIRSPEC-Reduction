      SUBROUTINE INPREJ ( LUNDRV, FAIL, ERRMSG )
C
C VERSION
C     17OCT13 AD Simplified. Remove GETREJ argument and local GOTREJ checks
C     14JUL98 AD Replace NXTREC with NXTFLD
C     04DEC97 AD Add GOTREJ to check for duplicate section
C     03MAR97 AD Version 3.
C     01OCT96 AD Version 2.
C     17SEP96 AD Remove redundant include file: rfmsiz.inc
C     01SEP96 AD Version 1.
C
C DESCRIPTION
C     Read user-defined line-rejection parameters from RFM driver table.
C     Called by RFMINP following *REJ marker in Driver table.
C     If REJFLG is TRUE read modified line rejection parameters
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      LUNDRV !  I  LUN for driver table
      LOGICAL      FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  ENDCHK ! Check end of Driver Table section has been reached.
     &, NXTFLD ! Load next field from section of RFM driver file
     &, RFMLOG ! Write text message to RFM log file
C
C COMMON VARIABLES
      INCLUDE 'rejcom.inc' ! Minimum Line strength limits
C
C LOCAL VARIABLES
      INTEGER       IOS     !  Saved value of IOSTAT for error message
      INTEGER       LENGTH  !  No.characters in FIELD
      CHARACTER*80  FIELD   !  Field read from driver table
      CHARACTER*132 MESSGE  !  Message sent to log file
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      CALL NXTFLD ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      READ ( FIELD, *, IOSTAT=IOS ) WIDREJ
      IF ( IOS .NE. 0 ) THEN
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,I11)' ) 
     &    'F-INPREJ: Error reading Min.line str. for Widemesh calcs.'
     &    //' IOSTAT=', IOS
        RETURN
      END IF
C
      CALL NXTFLD ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      READ ( FIELD, *, IOSTAT=IOS ) FINREJ
       IF ( IOS .NE. 0 ) THEN
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,I11)' ) 
     &    'F-INPREJ: Error reading Min.line str. for Finemesh calcs.'
     &    //' IOSTAT=', IOS
        RETURN
      END IF
C
      WRITE ( MESSGE, '(A,G9.2,A,G9.2)' )
     &  'I-INPREJ: Min.line strengths for wide,fine calcs=',
     &   WIDREJ, ',', FINREJ
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
      CALL ENDCHK ( LUNDRV, 'REJ', FAIL, ERRMSG )
C
      END
