      SUBROUTINE OBSATM ( NEWLEV, FAIL, ERRMSG )
C
C VERSION
C     25-OCT-03  AD  Also change Jacobian levels if already loaded
C     20-JAN-03  AD  Set CHGATM for extra level
C     05-OCT-00  AD  Bug fix: only set HGTATM once, not for every IGRA 
C     24-JAN-99  AD  Replace OFMFLG by LINFLG
C     17-DEC-98  AD  Bug fix: allow for VMR=0 when interpolating
C     18-JUL-98  AD  Original.
C
C DESCRIPTION
C     Set observer parameters and add level to profile
C     Called once by OBSCHK (OBS flag) following *OBS marker in Driver table.
C     Also called multiply by FLXPTH (FLX flag) to insert extra profile levels.
C
      IMPLICIT NONE
C
C ARGUMENTS
      LOGICAL      NEWLEV !  O  Set TRUE if an extra level has been added
      LOGICAL      FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  ATMAUX ! Set up auxiliary profiles of atmospheric parameters
     &, RFMLOG ! Write text record from RFM input file.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'gracom.inc' ! Atmospheric 2-D field data
      INCLUDE 'jaccom.inc' ! Jacobian data
      INCLUDE 'obscom.inc' ! Observer Position
C
C LOCAL CONSTANTS
      REAL          TOLLEV  !  Fraction of layer thickness which can be ignored
        PARAMETER ( TOLLEV = 0.001 ) ! when position observer at layer boundary
C
C LOCAL VARIABLES
      INTEGER       IATM    !  Atmospheric profile level index
      INTEGER       IGAS    !  Absorbing species counter
      INTEGER       IGRA    !  Horizontal gradient profile counter
      INTEGER       IJAC    !  Jacobian counter
      REAL          ALPHA   !  Profile interpolation fraction (0:1)
      CHARACTER*80  MESSGE  !  Message sent to log file
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      NEWLEV = .FALSE.
      IF ( ALTOBS .GT. HGTATM(NATM) ) THEN
        IATOBS = NATM
        FAIL = .FALSE.
        RETURN
      END IF
C
      IATM = NATM - 1
      DO WHILE ( HGTATM(IATM) .GT. ALTOBS )
        IATM = IATM - 1
      END DO
      ALPHA = ( ALTOBS - HGTATM(IATM) ) / 
     &        ( HGTATM(IATM+1) - HGTATM(IATM) )
      IF ( ALPHA .LE. TOLLEV ) THEN
        IATOBS = IATM
      ELSE IF ( 1.0 - ALPHA .LE. TOLLEV ) THEN
        IATOBS = IATM + 1
      ELSE                 ! need to insert an extra layer
        IATOBS = IATM + 1
        WRITE ( MESSGE, '(A,I5)' )
     &    'I-OBSATM: Inserting extra profile level#', IATOBS
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        IF ( NATM .EQ. MAXATM ) THEN
          FAIL = .TRUE.
          ERRMSG = 'F-OBSATM: no more space, array limit '//
     &             'MAXATM in RFMSIZ.INC already used'
          RETURN
        END IF
        NATM = NATM + 1
        CHGATM(NATM) = .TRUE.
        DO IATM = NATM, IATOBS + 1, -1
          HGTATM(IATM) = HGTATM(IATM-1)
        END DO
        HGTATM(IATOBS) = ALTOBS
C
        IGRA = 1
  100   CONTINUE
        IF ( GRAFLG ) CALL MOVGRA ( IGRA, 0 )
        DO IATM = NATM, IATOBS + 1, -1
          PREATM(IATM) = PREATM(IATM-1)
          TEMATM(IATM) = TEMATM(IATM-1)
          DO IGAS = 1, NGAS
            VMRATM(IATM,IGAS) = VMRATM(IATM-1,IGAS)
          END DO
        END DO
        IATM = IATOBS - 1
        TEMATM(IATOBS) = ALPHA * TEMATM(IATOBS) + 
     &           ( 1.0 - ALPHA ) * TEMATM(IATM)
        PREATM(IATOBS) = EXP ( ALPHA * LOG ( PREATM(IATOBS) ) + 
     &                       ( 1.0 - ALPHA ) * LOG ( PREATM(IATM) ) )
        DO IGAS = 1, NGAS
          IF ( LINFLG .OR. VMRATM(IATOBS,IGAS) .LT. ARGMIN 
     &                .OR. VMRATM(IATM,IGAS) .LT. ARGMIN ) THEN ! Interp.VMR
            VMRATM(IATOBS,IGAS) =  ALPHA * VMRATM(IATOBS,IGAS) + 
     &                     ( 1.0 - ALPHA ) * VMRATM(IATM,IGAS) 
          ELSE                     ! Linear interpolation in ln(VMR) wrt Z
            VMRATM(IATOBS,IGAS) = 
     &                   EXP ( ALPHA * LOG ( VMRATM(IATOBS,IGAS) ) + 
     &                 ( 1.0 - ALPHA ) * LOG ( VMRATM(IATM,IGAS) ) )
          END IF
        END DO
        IF ( GRAFLG ) THEN
          CALL MOVGRA ( 0, IGRA )
          IF ( IGRA .LT. NGRA ) THEN
            IGRA = IGRA + 1
            GOTO 100
          ELSE
            CALL MOVGRA ( NOMGRA, 0 ) 
          END IF
        END IF
C
C If Jacobians loaded, readjust perturbation limits for any level GE IATOBS
        IF ( JACFLG ) THEN
          DO IJAC = 1, NJAC         ! NJAC = 0 if not yet loaded
            IF ( IATJAC(IJAC) .GE. IATOBS ) 
     &         IATJAC(IJAC) = IATJAC(IJAC) + 1
            IF ( ILOJAC(IJAC) .GE. IATOBS ) 
     &         ILOJAC(IJAC) = ILOJAC(IJAC) + 1
            IF ( IUPJAC(IJAC) .GE. IATOBS ) 
     &         IUPJAC(IJAC) = IUPJAC(IJAC) + 1
          END DO
        END IF
C
        CALL ATMAUX        ! Recalculate auxiliary profiles
        NEWLEV = .TRUE.
      END IF
C
      FAIL = .FALSE.
C
      END
