      SUBROUTINE INPSFC ( LUNDRV, LUNSFC, FAIL, ERRMSG ) 
C 
C VERSION 
C     13OCT13 AD Simplified: remove GETSFC argument and local GOTSFC checks
C     07AUG13 AD Change FORMAT descriptors to avoid ifort warnings
C     23APR12 AD Allow for spectrally varying emissivity
C                Add LUNSFC argument. 
C     23JAN02 AD Increase precision in info message G10.3 to G10.5
C     14JUL98 AD Original.  
C 
C DESCRIPTION 
C     Read RFM surface parameters from driver table 
C     Called by RFMINP following *SFC marker in Driver table.  
C     If SFCFLG is TRUE read surface temperature and emissivity.  
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      LUNDRV !  I  LUN for driver table
      INTEGER      LUNSFC !  I  Temporary LUN for surface emissivity file
      LOGICAL      FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  CHKFIL ! Check if argument represents an existing filename
     &, ENDCHK ! Check end of Driver Table section has been reached.
     &, NXTFLD ! Load next field from section of RFM driver file
     &, RFMLOG ! Write text message to RFM log file
     &, SFCFIL ! Open, read surface emissivity data file and close
      LOGICAL CHKFIL
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'sfccom.inc' ! Surface parameters
C
C LOCAL CONSTANTS
      REAL          TEMDEL  ! Max. |SFC-ATM| Temp.difference before warning
        PARAMETER ( TEMDEL = 10.0 )
      REAL          TEMMIN  ! Minimum allowed surface temperature
        PARAMETER ( TEMMIN = 0.0 )
      REAL          TEMMAX  ! Maximum allowe surface temperature
        PARAMETER ( TEMMAX = 1000.0 )
C
C LOCAL VARIABLES
      INTEGER       IOS     !  Saved value of IOSTAT for error message
      INTEGER       LENGTH  !  No.of characters written in FIELD
      CHARACTER*80  FIELD   !  Field read from driver table
      CHARACTER*132 MESSGE  !  Message sent to log file
C
C EXECUTABLE CODE --------------------------------------------------------------
C
C First field is Surface Temperature (Mandatory)
      CALL NXTFLD ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      IF ( LENGTH .EQ. 0 ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-INPSFC: No data in *SFC section'
        RETURN
      END IF
C
      READ ( FIELD, *, IOSTAT=IOS ) TEMSFC
      IF ( IOS .NE. 0 ) THEN
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,I11)' )
     &    'F-INPSFC: Error reading Surface Temperature.'
     &    //' IOSTAT=', IOS
        RETURN
      END IF
      WRITE ( MESSGE, '(A,G10.3,A)' )
     &  'I-INPSFC: Reading Surface Temperature = ', TEMSFC, ' K'
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
      IF ( TEMSFC .LT. TEMMIN .OR. TEMSFC .GT. TEMMAX ) THEN       
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,F8.0,A,F8.0)' )
     &     'F-INPSFC: Surface Temperature outside allowed range: ',
     &     TEMMIN, ' to ', TEMMAX
        RETURN
      END IF
C
      IF ( ABS ( TEMSFC-TEMATM(1) ) .GT. TEMDEL ) THEN
        WRITE ( MESSGE, '(A,F8.2,A)' )
     &   'W-INPSFC: Large SFC-ATM Temperature discontinuity, =',
     &   (TEMSFC-TEMATM(1)), ' K'
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
C
C Second field is Surface Emissivity (Optional - default value = 1, ie BB Sfc.)
      CALL NXTFLD ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG )
      IF ( LENGTH .EQ. 0 ) THEN             ! assume default BB surface
        NSFC = 1
        EMSSFC(1) = 1.0D0
        WNOSFC(1) = 0.0D0
        MESSGE = 'I-INPSFC: Assuming default Surface Emissivity=1'
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        RETURN
      END IF
C
      IF ( CHKFIL ( LUNSFC, FIELD ) ) THEN 
        CALL SFCFIL ( LUNSFC, FIELD, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      ELSE
        READ ( FIELD, *, IOSTAT=IOS ) EMSSFC(1)
        IF ( IOS .NE. 0 ) THEN
          FAIL = .TRUE.
          WRITE ( ERRMSG, '(A,I11)' )
     &      'F-INPSFC: Error reading Surface Emissivity value.'
     &    //' IOSTAT=', IOS
          RETURN
        END IF
        NSFC = 1
        WNOSFC(1) = 0.0D0
        WRITE ( MESSGE, '(A,G12.5)' )
     &    'I-INPSFC: Setting Surface Emissivity = ', EMSSFC(1)
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        IF ( EMSSFC(1) .LT. 0.D0 .OR. EMSSFC(1) .GT. 1.D0 ) THEN       
          FAIL = .TRUE.
          ERRMSG = 
     &      'F-INPSFC: Surface Emissivity outside range 0:1'
          RETURN
        END IF
      END IF
C Set flag for reflecting surface
      RFLSFC = ( NSFC .GT. 1 .OR. EMSSFC(1) .LT. 1.0D0 )
C
C Check that there are no more fields in this section
      CALL ENDCHK ( LUNDRV, 'SFC', FAIL, ERRMSG )
C
      END
