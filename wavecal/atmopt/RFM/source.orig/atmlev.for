      SUBROUTINE ATMLEV ( ALT, IDXATM, FAIL, ERRMSG )
C
C VERSION
C     01-JAN-04  AD  Original. Based on OBSATM, which it replaces.
C
C DESCRIPTION
C     Find/insert atmospheric level for given altitude.
C     Called by FLXPTH, JACCHK, LEVCHK and OBSCHK.
C     If the ALTitude is above the top of the atmosphere or below the base of
C       the atmosphere a fatal error message results
C     If the ALTitude exactly matches an existing atmospheric level, the index
C       IATM of that level is returned. 
C     If the ALTitude matches an existing atmospheric level within a given 
C       tolerance (set by TOLLEV) the ALTitude is adjusted to that level and 
C       the index returned.
C     If the ALTitude does not match, a new level is inserted and its index is
C       returned
C
      IMPLICIT NONE
C
C ARGUMENTS
      REAL         ALT    ! I/O Altitude
      INTEGER      IDXATM !  O  Index of atmospheric level
      LOGICAL      FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  ATMAUX ! Set up auxiliary profiles of atmospheric parameters
     &, MOVGRA ! Exchange profiles between /ATMCOM/ and /GRACOM/ arrays.
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
        PARAMETER ( TOLLEV = 0.001 ) ! if ALT close to layer boundary
C
C LOCAL VARIABLES
      INTEGER       IATM    !  Atmospheric profile level index
      INTEGER       IGAS    !  Absorbing species counter
      INTEGER       IGRA    !  Horizontal gradient profile counter
      INTEGER       IJAC    !  Jacobian counter
      REAL          ALPHA   !  Profile interpolation fraction (0:1)
      CHARACTER*80  WRNMSG  !  Message sent to log file
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Argument ALT should be checked before this subroutine is called, so this 
C error message should never arise
      IF ( ALT .GT. HGTATM(NATM) .OR. ALT .LT. HGTATM(1) ) THEN
        ERRMSG = 'F-ATMLEV: Alt is outside range of atmosphere'
        FAIL = .TRUE.
        RETURN
      END IF
C
      IATM = NATM - 1
      DO WHILE ( HGTATM(IATM) .GT. ALT )
        IATM = IATM - 1
      END DO
      ALPHA = ( ALT - HGTATM(IATM) ) / 
     &        ( HGTATM(IATM+1) - HGTATM(IATM) )
      IF ( ALPHA .LE. TOLLEV ) THEN
        IDXATM = IATM
      ELSE IF ( 1.0 - ALPHA .LE. TOLLEV ) THEN
        IDXATM = IATM + 1
      ELSE                 ! need to insert an extra layer
        IDXATM = IATM + 1
        WRITE ( WRNMSG, '(A,I5)' )
     &    'W-ATMLEV: Inserting extra profile level#', IDXATM
        CALL RFMLOG ( WRNMSG, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        IF ( NATM .EQ. MAXATM ) THEN
          FAIL = .TRUE.
          ERRMSG = 'F-ATMLEV: no more space, array limit '//
     &             'MAXATM in rfmsiz.inc already used'
          RETURN
        END IF
        NATM = NATM + 1
        CHGATM(NATM) = .TRUE.
        DO IATM = NATM, IDXATM + 1, -1
          HGTATM(IATM) = HGTATM(IATM-1)
        END DO
        HGTATM(IDXATM) = ALT
C
        IGRA = 1
  100   CONTINUE
        IF ( GRAFLG ) CALL MOVGRA ( IGRA, 0 )
        DO IATM = NATM, IDXATM + 1, -1
          PREATM(IATM) = PREATM(IATM-1)
          TEMATM(IATM) = TEMATM(IATM-1)
          DO IGAS = 1, NGAS
            VMRATM(IATM,IGAS) = VMRATM(IATM-1,IGAS)
          END DO
        END DO
        IATM = IDXATM - 1
        TEMATM(IDXATM) = ALPHA * TEMATM(IDXATM) + 
     &           ( 1.0 - ALPHA ) * TEMATM(IATM)
        PREATM(IDXATM) = EXP ( ALPHA * LOG ( PREATM(IDXATM) ) + 
     &                       ( 1.0 - ALPHA ) * LOG ( PREATM(IATM) ) )
        DO IGAS = 1, NGAS
          IF ( LINFLG .OR. VMRATM(IDXATM,IGAS) .LT. ARGMIN 
     &                .OR. VMRATM(IATM,IGAS) .LT. ARGMIN ) THEN ! Interp.VMR
            VMRATM(IDXATM,IGAS) =  ALPHA * VMRATM(IDXATM,IGAS) + 
     &                     ( 1.0 - ALPHA ) * VMRATM(IATM,IGAS) 
          ELSE                     ! Linear interpolation in ln(VMR) wrt Z
            VMRATM(IDXATM,IGAS) = 
     &                   EXP ( ALPHA * LOG ( VMRATM(IDXATM,IGAS) ) + 
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
C If Observer altitude loaded and above extra level, add one to index
        IF ( OBSFLG .AND. IATOBS .GE. IDXATM ) IATOBS = IATOBS + 1
C
C If Jacobians loaded, readjust perturbation limits for any level GE IDXATM
C NJAC should be 0 if not yet loaded, but to be safe use MIN in case NJAC is 
C uninitialised and set to a large integer
        IF ( JACFLG .OR. LEVFLG ) THEN
          DO IJAC = 1, MIN ( NJAC, MAXJAC )      ! NJAC = 0 if not yet loaded
            IF ( IATJAC(IJAC) .GE. IDXATM ) 
     &         IATJAC(IJAC) = IATJAC(IJAC) + 1
            IF ( ILOJAC(IJAC) .GE. IDXATM ) 
     &         ILOJAC(IJAC) = ILOJAC(IJAC) + 1
            IF ( IUPJAC(IJAC) .GE. IDXATM ) 
     &         IUPJAC(IJAC) = IUPJAC(IJAC) + 1
          END DO
        END IF
C
        CALL ATMAUX        ! Recalculate auxiliary profiles
      END IF
C
C Adjust altitude to match profile levels
      IF ( ALT .NE. HGTATM(IDXATM) ) THEN
        WRITE ( WRNMSG, '(A,G10.3,A,G10.3,A)' ) 
     &    'W-ATMLEV: Change altitude from', ALT, 'km to', 
     &    HGTATM(IDXATM), 'km to match profile.'
        CALL RFMLOG ( WRNMSG, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        ALT = HGTATM(IDXATM)
      END IF
C
      FAIL = .FALSE.
C
      END
