      SUBROUTINE INPNAM ( LUNDRV, TYP, FAIL, ERRMSG )
C
C VERSION
C     17OCT13 AD Original. Based on INPABS etc
C
C DESCRIPTION
C     Read user-defined output filenames from RFM driver table.
C     Called by RFMINP following *ABS, *BBT etc markers in Driver table.
C     This replaces the original different routines INPABS, INPBBT etc
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      LUNDRV !  I  LUN for driver table
      CHARACTER*3  TYP    !  I  File type (eg 'ABS', 'BBT' ) - upper case reqd
      LOGICAL      FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  ENDCHK ! Check end of Driver Table section has been reached.
     &, NXTFLD ! Load next field from section of RFM driver file
     &, RFMLOG ! Write text message to RFM log file
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'gracom.inc' ! Atmospheric 2-D field data
      INCLUDE 'outcom.inc' ! RFM output file data
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL VARIABLES
      INTEGER       IPT    !  Pointer to '*' in supplied filename
      INTEGER       LENGTH !  No.characters in FIELD
      LOGICAL       USEGAS ! T=output filename includes gas identification
      LOGICAL       USEGRA ! T=output filename includes psi-angle info
      LOGICAL       USEJAC ! T=output filename includes Jac/LOS identification
      LOGICAL       USETAN ! T=output filename includes tangent height info
      LOGICAL       USESPC ! T=output filename includes spectral range info
      CHARACTER*200 NAMFIL !  Field read from driver table
      CHARACTER*238 MESSGE !  Message sent to log file
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      CALL NXTFLD ( LUNDRV, NAMFIL, LENGTH, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      IF ( LENGTH .EQ. 0 ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-INPNAM: No filename supplied in *'//TYP//' section'
        RETURN
      END IF
C               123456789012345678901234   567   89012345678
      MESSGE = 'I-INPNAM: User-supplied '//TYP//' filename: '
     &             //NAMFIL
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
C Check if the filename has a '*' symbol and, if not, whether it is required
      IPT = INDEX ( NAMFIL, '*' )
      IF ( IPT .EQ. 0 ) THEN
        USEGAS = .FALSE.
        USEGRA = .FALSE.
        USEJAC = .FALSE.
        USETAN = .FALSE. 
        USESPC = .FALSE.
        IF ( TYP .EQ. 'PRF' ) THEN
          USEGRA = .TRUE.
        ELSE IF ( TYP .EQ. 'PTH' ) THEN
          USETAN = .TRUE.
        ELSE IF ( TYP .EQ. 'TAB' ) THEN
          USEGAS = .TRUE.
          USESPC = .TRUE.
        ELSE IF ( TYP .EQ. 'WID' ) THEN
          USEGAS = .TRUE.
          USESPC = .TRUE.
        ELSE                     ! standard spectral output files
          USETAN = .TRUE.
          USESPC = .TRUE.
          USEJAC = .TRUE.
        END IF          
C
        FAIL = .TRUE.
        IF ( USEGAS .AND. NGAS .GT. 1 ) THEN
          ERRMSG = 'F-INPNAM: No ''*'' symbol in '//TYP//
     &             ' filename for inserting gas name'
        ELSE IF ( USEGRA .AND. GRAFLG .AND. NGRA .GT. 1 ) THEN
            ERRMSG = 'F-INPNAM: No ''*'' symbol in '//TYP//
     &               ' filename for inserting psi-angle info'
        ELSE IF ( USETAN .AND. NTAN .GT. 1 ) THEN
            ERRMSG = 'F-INPNAM: No ''*'' symbol in '//TYP//
     &               ' filename for inserting tangent height info'
        ELSE IF ( USESPC .AND. NSPC .GT. 1 ) THEN
            ERRMSG = 'F-INPNAM: No ''*'' symbol in '//TYP//
     &               ' filename for inserting spectral range info'
        ELSE IF ( USEJAC .AND. ( JACFLG .OR. LOSFLG ) ) THEN
            ERRMSG = 'F-INPNAM: No ''*'' symbol in '//TYP//
     &               ' filename for inserting Jacobian element info'
        ELSE
          FAIL = .FALSE.
        END IF
        IF ( FAIL ) RETURN
      END IF
C
C If valid filename, save as appropriate common variable in outcom.inc
      IF ( TYP .EQ. 'ABS' ) THEN
        NAMABS = NAMFIL
      ELSE IF ( TYP .EQ. 'BBT' ) THEN
        NAMBBT = NAMFIL
      ELSE IF ( TYP .EQ. 'OPT' ) THEN
        NAMOPT = NAMFIL
      ELSE IF ( TYP .EQ. 'PRF' ) THEN
        NAMPRF = NAMFIL
      ELSE IF ( TYP .EQ. 'PTH' ) THEN
        NAMPTH = NAMFIL
      ELSE IF ( TYP .EQ. 'RAD' ) THEN
        NAMRAD = NAMFIL
      ELSE IF ( TYP .EQ. 'RJT' ) THEN
        NAMRJT = NAMFIL
      ELSE IF ( TYP .EQ. 'TAB' ) THEN
        NAMTAB = NAMFIL
      ELSE IF ( TYP .EQ. 'TRA' ) THEN
        NAMTRA = NAMFIL
      ELSE IF ( TYP .EQ. 'WID' ) THEN
        NAMWID = NAMFIL
      ELSE
        STOP 'F-INPNAM: Logical Error - unrecognised TYP'
      END IF
C
C Check that there are no more fields in this Driver file section
      CALL ENDCHK ( LUNDRV, TYP, FAIL, ERRMSG )
C
      END
