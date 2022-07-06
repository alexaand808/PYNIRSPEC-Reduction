      SUBROUTINE TANCHK ( TANSTR, FAIL, ERRMSG )
C
C VERSION
C     13OCT13  AD  Remove NOTTAN argument, fatal error if not valid number
C     07AUG13  AD  Remove redundant local variable GETOBS
C     12DEC07  AD  Bug#69: force minus sign instead of empty space in TANSTR
C                    when small negative tangent heights are rounded to 0.
C     18MAR01  AD  Remove CHKLIM
C     03MAY00  AD  Insert TANSTR into sorted list for all flags
C     16APR00  AD  Add CHKFLX  
C     18JUL98  AD  Allow for Elevation Angles and Geom.Tangent Heights
C                    Add CHKHOM, CHKNAD, CHKLIM. Remove CHKWID.
C     22JAN98  AD  Remove redundant EXTERNAL declaration: RFMLOG
C     22OCT97  AD  Add CHKTAB if TAB option selected.
C     03MAR97  AD  Version 3.
C     15JAN97  AD  Check for airmass for ZEN or NAD viewing
C                    Refine DHTMIN to be a fractional difference
C     01OCT96  AD  Version 2.
C     30SEP96  AD  Add '/' to list of first characters interpreted as filenam
C     29SEP96  AD  Add DHTMIN to check minimum separation of 1m
C     18SEP96  AD  Add CHKWID
C     01SEP96  AD  Version 1.
C
C DESCRIPTION
C     Check if string is legal Tangent Height.
C     Called by INPTAN for each field in *TAN section, or by TANFIL.
C     Note that this just loads reqd output tangent heights, more may be added
C     by FOVTAN (up to MTAN) for the FOV convolution.
C
      IMPLICIT NONE
C
      EXTERNAL
     &  CHKFLX ! Check output levels for Flux calculations
     &, CHKHOM ! Check path length for Homogeneous path calculation
     &, CHKNAD ! Check path for Nadir/Zenith viewing mode
     &, CHKTAB ! Check Tabulation parameters.
C
C ARGUMENTS
      CHARACTER*(*) TANSTR  !  I  Tangent string to be tested 
      LOGICAL       FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG  !  O  Error message returned if FAIL is TRUE
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL VARIABLES
      INTEGER     ILEN   ! Length of TANSTR (for error message)
      INTEGER     IOS    ! Saved value of IOSTAT
      INTEGER     ITAN   ! Counter for listed tangent heights
      INTEGER     JTAN   ! Pointer for insertion of new tangent height
      REAL        TANTST ! Number read from TANSTR as possible tangent height
      CHARACTER*6 STRTST ! Tan.value written as a string for output filename
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      FAIL = .TRUE.
C
C To be identified as a real number, the string TANSTR must be readable without
C an error
      READ ( TANSTR, *, IOSTAT=IOS ) TANTST 
      IF ( IOS .NE. 0 ) THEN
        ILEN = MIN ( LEN ( TANSTR ), 30 )
        WRITE ( ERRMSG, '(A,A)' ) 
     &    'F-TANCHK: Unreadable value in *TAN section:', TANSTR(1:ILEN)
        RETURN
      END IF
C
C If producing tabulated absorption coefficients, TANTST is a tabulation param.
C so use separate routine to perform checks and return for next parameter
      IF ( TABFLG ) THEN
        CALL CHKTAB ( TANTST, FAIL, ERRMSG )
        RETURN
C
      ELSE IF ( FLXFLG ) THEN     ! NB must be called before ZEN/NAD test
        CALL CHKFLX ( TANTST, FAIL, ERRMSG ) ! since FLX+ZEN/NAD is possible
        IF ( FAIL ) RETURN
C
      ELSE IF ( ZENFLG .OR. NADFLG ) THEN
        CALL CHKNAD ( TANTST, USRELE, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
C
      ELSE IF ( HOMFLG ) THEN
        CALL CHKHOM ( TANTST, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
C
C Limb checks depend on inputs from *OBS, *CRV and *FOV sections so are 
C postponed until INPCHK
C
      END IF           ! end case of viewing geometry
C
C Check enough array space to add extra tangent height
C
      IF ( NTAN .EQ. MAXTAN ) THEN
        WRITE ( ERRMSG, '(A,I4,A)' )
     &    'F-TANCHK: No.Tan Hgts required > MAXTAN=', MAXTAN,   
     &    ' in RFMSIZ.INC'
        RETURN
      END IF
C
C Check that this is distinguishable from other USRTAN values, and, for limb
C views, sort into ascending order of USRTAN.
C
      JTAN = 1
      IF ( ABS ( TANTST ) .GE. 99999.5 ) THEN
        ERRMSG = 'F-TANCHK: Cannot handle values .GE. 99999.5'
        RETURN
      END IF 
      IF ( TANTST .GE. 0.0 ) THEN
        IF ( NINT ( TANTST * 1000.0 ) .LT. 99999 ) THEN
          WRITE ( STRTST, '(I5.5)' ) NINT ( TANTST * 1000.0 ) 
        ELSE 
          WRITE ( STRTST, '(I5.5)' ) NINT ( TANTST )
        END IF 
      ELSE
        IF ( NINT ( TANTST * 1000.0 ) .GT. -99999 ) THEN
          WRITE ( STRTST, '(I6.5)' ) NINT ( TANTST * 1000.0 ) 
C For small negative values, there's a chance that rounding may lead to 0
C in which case ensure that there is a minus sign otherwise there will be a 
C space in the output filename
          IF ( STRTST .EQ. ' 00000' ) STRTST = '-00000'
        ELSE 
          WRITE ( STRTST, '(I6.5)' ) NINT ( TANTST )
        END IF 
      END IF
C
      DO ITAN = 1, NTAN
        IF ( STRTST .EQ. STRTAN(ITAN) ) THEN 
          WRITE ( ERRMSG, '(A,G12.3)' ) 
     &      'F-TANCHK: Repeated Tan.Hgt=', TANTST
          RETURN
        ELSE IF ( USRTAN(ITAN) .LT. TANTST ) THEN
          JTAN = ITAN + 1
        END IF                                ! Just continue if USRTAN > TANTST
      END DO
      DO ITAN = NTAN, JTAN, -1                ! Insert into sorted list
        USRTAN(ITAN+1) = USRTAN(ITAN)
        STRTAN(ITAN+1) = STRTAN(ITAN)
      END DO
      USRTAN(JTAN) = TANTST
      STRTAN(JTAN) = STRTST
      NTAN = NTAN + 1
      CLCTAN(NTAN) = .TRUE.
      FAIL = .FALSE. 
C
      END

