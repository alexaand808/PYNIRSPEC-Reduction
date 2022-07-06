      SUBROUTINE INPFIN ( LUNDRV, FAIL, ERRMSG )
C
C VERSION
C     16OCT13 AD Simplified. Remove GETFIN argument. Add ENDCHK
C     07AUG13 AD Change FORMAT descriptors to avoid ifort warnings
C     17AUG06 AD If GHZ flag and value < 1, interpret as GHz spacing
C     18MAR01 AD Remove GETILS argument, move ILSCHK, SPCFIN to INPCHK
C     15MAR01 AD Add GETILS argument,
C     05AUG99 AD Change USRFIN from Real to Double Precisions 
C     23APR99 AD Change format to allow for C*8 LABSPC
C     27NOV98 AD Check user-resln adequate for all output grids
C     14JUL98 AD Change to use NXTFLD instead of NXTREC
C     17DEC97 AD Comment changes only
C     04DEC97 AD Add GOTFIN to check for duplicate section
C     03MAR97 AD Version 3.
C     01OCT96 AD Version 2.
C     18SEP96 AD Original.
C
C DESCRIPTION
C     Read user-defined fine-mesh resolution 
C     Called by RFMINP following *FIN marker in Driver table and FINFLG=TRUE
C     Note that this is only required if coupled with an ILS convolution -
C     otherwise the fine mesh resolution is the same as the output resolution. 
C     If a finemesh resolution is not specified, a default value is used.
C
C     Note that any user-specified resolution >1 is interpreted as pts/cm-1,
C     while values <1 are interpreted as spacing in cm-1.
C
C     If GHZ flag is enabled, then any value <1 is interpreted as GHz spacing
C     but values >1 are still pts/cm-1.
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
C GLOBAL CONSTANTS
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C LOCAL CONSTANTS
      DOUBLE PRECISION WNOFAC ! GHz to Wno conversion factor
        PARAMETER ( WNOFAC = 1.0D7 / VLIGHT )   ! VLIGHT in phycon.inc
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL VARIABLES
      INTEGER       IOS     !  Saved value of IOSTAT for error message
      INTEGER       LENGTH  !  No. of characters in FIELD
      DOUBLE PRECISION USRFIN  ! User-specified fine mesh resolution
      DOUBLE PRECISION GHZFIN  ! Original resln in GHz before conversion
      CHARACTER*80  FIELD   !  Field read from driver table
      CHARACTER*132 MESSGE  !  Message sent to log file
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      CALL NXTFLD ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      IF ( LENGTH .EQ. 0 ) THEN
        FAIL = .TRUE. 
        ERRMSG = 'F-INPFIN: No data in *FIN section of Driver Table'
        RETURN
      END IF
C
      READ ( FIELD, *, IOSTAT=IOS ) USRFIN
      IF ( IOS .NE. 0 ) THEN
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,I11)' ) 
     &    'F-INPFIN: Error reading fine-mesh resolution.'
     &    //' IOSTAT=', IOS
        RETURN
      END IF      
C
C If value < 1 and GHZ flag, read as GHz spacing and convert to cm-1 spacing
      IF ( USRFIN .LT. 1.0 .AND. USRFIN .GT. 0.0 .AND. GHZFLG ) THEN
        GHZFIN = USRFIN
        USRFIN = GHZFIN * WNOFAC
        WRITE ( MESSGE, '(A,1PG10.3,A,G14.7,A)' ) 
     &    'I-INPFIN: Converting *FIN spacing from ', GHZFIN,
     &    ' GHz to ', USRFIN, ' cm-1'
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
C
      IF ( NINT ( USRFIN ) .GT. MAXFIN .OR. 
     &     SNGL ( USRFIN ) .LT. 1.0/FLOAT(MAXFIN) ) THEN  
        FAIL = .TRUE.                              ! Also prevents USRFIN=0
        WRITE ( ERRMSG, '(A,I6,A)' )
     &    'F-INPFIN: user-specified resolution > MAXFIN=', MAXFIN,
     &    ' pts/cm-1'
        RETURN
      ELSE IF ( USRFIN .GT. 1.0D0 ) THEN
        NOMFIN = NINT ( USRFIN )
      ELSE
        NOMFIN = NINT ( 1.D0 / USRFIN ) 
      END IF
C
      WRITE ( MESSGE, '(A,I11,A)' )
     &  'I-INPFIN: Using fine mesh calculation @', NOMFIN,
     &  ' pts/cm-1'
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
C Check no more data in *FIN section
      CALL ENDCHK ( LUNDRV, 'FIN', FAIL, ERRMSG )
C
      END
