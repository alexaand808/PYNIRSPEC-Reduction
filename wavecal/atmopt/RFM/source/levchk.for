      SUBROUTINE LEVCHK ( LEVTST, FAIL, ERRMSG )
C
C VERSION
C     01-JAN-04  AD  Original.
C
C DESCRIPTION
C     Check if string is legal altitude level.
C     Called by INPLEV for each field in *LEV section, or by LEVFIL.
C     This loads altitude levels into IATJAC arrays, using NJAC as a counter.
C
      IMPLICIT NONE
C
C ARGUMENTS
      REAL          LEVTST  !  I  Altitude level to checked
      LOGICAL       FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  ATMLEV ! Find/insert atmospheric level for given altitude.
C
C GLOBAL CONSLEVTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'jaccom.inc' ! Jacobian data
C
C LOCAL VARIABLES
      INTEGER IATM ! Index of atmospheric level corresponding to LEVTST
      INTEGER IJAC ! Index of output levels stored by increasing altitude
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Check that value lies within atmosphere
      FAIL = .TRUE.
      IF ( ABS ( LEVTST ) .GE. 99999.5 ) THEN
        ERRMSG = 'F-LEVCHK: Cannot handle values .GE. 99999.5'
      ELSE IF ( LEVTST .GT. HGTATM(NATM) ) THEN
        ERRMSG = 'F-LEVCHK: Output level is above top of atmosphere'
      ELSE IF ( LEVTST .LT. HGTATM(1) ) THEN
        ERRMSG = 'F-LEVCHK: Output level is below base of atmosphere'
      ELSE
        FAIL = .FALSE.
      END IF
      IF ( FAIL ) RETURN
C
C Check enough array space to add extra output level.
      IF ( NJAC .EQ. MAXJAC ) THEN
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,I4,A)' )
     &    'F-LEVCHK: No.o/p Levs required > MAXJAC=', MAXJAC,   
     &    ' in rfmsiz.inc'
        RETURN
      END IF
C
C Find level in atmosphere, inserting extra level if required
C ATMLEV should adjust values stored in IATJAC if extra level is added
      CALL ATMLEV ( LEVTST, IATM, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      NJAC = NJAC + 1
C
C Insert into list ordered by increasing altitude, checking for repeated value
      IJAC = NJAC - 1
      DO WHILE ( IATJAC(IJAC) .GT. IATM ) 
        IATJAC(IJAC+1) = IATJAC(IJAC)
        IJAC = IJAC - 1
      END DO
      IF ( IATJAC(IJAC) .EQ. IATM ) THEN
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,G12.3)' ) 
     &    'F-LEVCHK: Repeated Lev.Alt=', LEVTST
        RETURN
      END IF
      IATJAC(IJAC+1) = IATM
C
      FAIL = .FALSE.
C
      END
