      SUBROUTINE LEVTAN ( FAIL, ERRMSG )
C
C VERSION
C     23-APR-12  AD  Use RFLSFC rather than EMSSFC to determine reflection
C     03-MAY-05  AD  Add jdxcon.inc 
C     01-JAN-04  AD  Original.
C
C DESCRIPTION
C     Set up tangent paths for intermediate output levels
C     Called once by INPCHK if LEV flag is enabled. 
C
      IMPLICIT NONE 
C
C ARGUMENTS
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  RFMLOG ! Write text message to RFM log file.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'jdxcon.inc' ! Indices for non-VMR Jacobians
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'jaccom.inc' ! Jacobian data
      INCLUDE 'obscom.inc' ! Observer Position
      INCLUDE 'sfccom.inc' ! Surface parameters
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL VARIABLES
      INTEGER      IATM    ! Atmospheric index for current output level
      INTEGER      IATLEV(MAXJAC) ! Atmospheric indices of output levels
      INTEGER      IDIR    ! Ray integration direction, -1=down, +1=up
      INTEGER      ILEV    ! Intermediate output level counter
      INTEGER      ILEV1,ILEV2 ! Start/end level indices for up/down paths
      INTEGER      ITAN    ! Counter for output tangent heights
      INTEGER      NLEV    ! Number of intermediate output levels selected
      LOGICAL      USELEV  ! T=o/p level applies to tan.path, F=doesn't.
      REAL         HGTLEV  ! Altitude [km] of output level
      CHARACTER*80 MESSGE  ! Text message for Log file
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Copy data stored temporarily in IATJAC into local IATLEV
      NLEV = NJAC
      DO ILEV = 1, NLEV
        IATLEV(ILEV) = IATJAC(ILEV)
      END DO
C
      NJAC = 0             ! Re-initialise
C LTAN Should already be initialised to NTAN in INPTAN 
      IF ( LTAN .NE. NTAN ) STOP 'F-LEVTAN: Logical Error#1'
C
C For each ray, starting at observer, make list of output levels which apply
      DO ITAN = 1, NTAN
        DO IDIR = -1, 1, 2           ! Downward integration then upward
C Need to set up IJAC, LTAN ordering to match sequence in RFMRAD so that 
C a continuous sequence of indices are accessed in each rad.transf calc.
          IF ( IDIR .EQ. -1 ) THEN
            ILEV1 = NLEV
            ILEV2 = 1
          ELSE
            ILEV1 = 1
            ILEV2 = NLEV
          END IF
C Loop over levels establishing which levels apply in down/up paths
          DO ILEV = ILEV1, ILEV2, IDIR
            IATM = IATLEV(ILEV)
            HGTLEV = HGTATM(IATM)
            IF ( ZENFLG ) THEN
              IF ( IDIR .EQ. -1 ) THEN 
                USELEV = .FALSE.
              ELSE
                USELEV = .NOT. OBSFLG .OR.         ! if not OBSFLG use all levs
     &                   ( HGTLEV .GT. ALTOBS )    ! else levs above observer
              END IF
            ELSE IF ( NADFLG ) THEN
              IF ( IDIR .EQ. -1 ) THEN
                USELEV = .NOT. OBSFLG .OR.
     &                   ( HGTLEV .LT. ALTOBS )
              ELSE
                USELEV = RFLSFC
              END IF
            ELSE
              IF ( IDIR .EQ. -1 ) THEN
                USELEV = HGTLEV .GT. HGTTAN(ITAN) .AND.
     &                   ( .NOT. OBSFLG .OR. HGTLEV .LT. ALTOBS )
              ELSE
                USELEV = HGTLEV .GE. HGTTAN(ITAN) .AND.
     &                   ( .NOT. SFCTAN(ITAN) .OR. RFLSFC )
              END IF
            END IF
C
            IF ( USELEV ) THEN         ! Add this level as new LTAN,NJAC
              IF ( NJAC .EQ. MAXJAC ) THEN 
                FAIL = .TRUE.
                ERRMSG = 'F-LEVTAN: MAXJAC in rfmsiz.inc '//
     &                   'too small to store intermediate level data'
                RETURN
              END IF
              NJAC = NJAC + 1
              IF ( LTAN .EQ. MAXTAN ) THEN
                FAIL = .TRUE.
                ERRMSG = 'F-LEVTAN: MAXTAN in rfmsiz.inc '//
     &                   'too small to store intermediate level data'
                RETURN
              END IF
              LTAN = LTAN + 1
C
              ITNJAC(ITAN,NJAC) = LTAN
              ITNJAC(LTAN,NJAC) = ITAN
C IATJAC is used to get altitude to form output filename in JACNAM
              IATJAC(NJAC) = IATM 
C ILOJAC is used in LEVINI to set layer to start radiative transfer integration
C On downward path, start integration for Lev#i from layer i-1.
              IF ( IDIR .EQ. -1 ) THEN
                ILOJAC(NJAC) = -IATJAC(NJAC) + 1
              ELSE
                ILOJAC(NJAC) = IATJAC(NJAC)
              END IF
C Setting IGSJAC results in JACNAM inserting 'down' and 'up' into o/p filenames
              IF ( IDIR .EQ. -1 ) THEN
                IGSJAC(NJAC) = MAXGAS + JDXLVD
              ELSE
                IGSJAC(NJAC) = MAXGAS + JDXLVU
              END IF
            END IF
          END DO
        END DO
      END DO
C
      WRITE ( MESSGE, '(A,I11)' )
     &  'I-LEVTAN: No.extra tangent paths reqd for '//
     &  'intermediate Lev output=', LTAN - NTAN
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
C
      END
