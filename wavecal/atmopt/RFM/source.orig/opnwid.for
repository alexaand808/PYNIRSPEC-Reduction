      SUBROUTINE OPNWID ( LUNNXT, ISPC, FAIL, ERRMSG )
C
C VERSION
C     06-APR-06  AD  Remove leading 'X' in format
C     06-JAN-04  AD  Add IGRA argument to MAKNAM
C     05-DEC-00  AD  Comment out CARRIAGECONTROL='LIST' 
C     10-MAR-97  AD  Add CARRIAGECONTROL='LIST'
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Open WID files for current spectral range
C     Called by RFMOPN for each spectral range.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      LUNNXT  ! I/O Next available LUN
      INTEGER      ISPC    !  I  Spectral range index
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error occurs
      CHARACTER*80 ERRMSG  !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  MAKNAM ! Construct filename for RFM output files.
     &, RFMLOG ! Write text message to RFM log file.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'rfmcon.inc' ! Constants for RFM calculations
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'outcom.inc' ! RFM output file data
      INCLUDE 'widcom.inc' ! Wide mesh data
C
C LOCAL VARIABLES
      INTEGER       IGAS   ! Absorber counter
      INTEGER       IOS    ! Value of IOSTAT saved for error messages.
      INTEGER       LUN    ! Local value of Logical Unit Number.
      INTEGER       MGAS   ! Initial counter for WID output files (0:1)
      CHARACTER*132 MESSGE ! Text message sent to Log file.
      CHARACTER*80  NAMFIL ! Name of file actually opened (incl. RUNID)
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      LUNWID = LUNNXT
      IF ( NGAS .EQ. 1 ) THEN      ! Only write out single gas file
        MGAS = 1
      ELSE                         ! Write file for each gas, plus total file
        MGAS = 0
      END IF
      LUNNXT = LUNNXT + NGAS - MGAS + 1 ! NGAS files, +1 for total of all gases
      LUN = LUNWID - 1
      DO IGAS = MGAS, NGAS
        LUN = LUN + 1
        CALL MAKNAM ( ISPC, IGAS, 0, 0, NAMWID, NAMFIL )
        MESSGE = 'I-OPNWID: Opening output file: '//NAMFIL
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        IF ( NEWFLG ) THEN
          OPEN ( UNIT=LUN, FILE=NAMFIL, STATUS='NEW',
     &           IOSTAT=IOS, ERR=900 )
        ELSE
          OPEN ( UNIT=LUN, FILE=NAMFIL, STATUS='UNKNOWN',
     &           IOSTAT=IOS, ERR=900 )
        END IF
        IF ( IGAS .EQ. 0 ) THEN
          WRITE ( LUN, '(A)', IOSTAT=IOS, ERR=900 )
     &      '! Total Widemesh Line count created by RFM v.'// VERSN
        ELSE 
          WRITE ( LUN, '(A)', IOSTAT=IOS, ERR=900 )
     &     '! '//CODGAS(IGAS)//' Widemesh Line count created by RFM v.'
     &     // VERSN
        END IF
        WRITE ( LUN, '(A,A)', IOSTAT=IOS, ERR=900 ) '!',TXTHDR
        WRITE ( LUN, * ) NWID, WN1WID, DELWID,
     &    ' =NWid,Wno1,DWno'
        WRITE ( LUN, '(A)' ) 'Itvl#  Wavenumber Tot.Lines '//
     &    '= Wide mesh + Fine mesh  Local  X/Section Continuum'
      END DO
C
      FAIL = .FALSE.
      RETURN
C
  900 FAIL = .TRUE.
      WRITE ( ERRMSG, '(A,I11)' ) 
     &  'F-OPNWID: I/O failure on output file. IOSTAT=', IOS
      END
