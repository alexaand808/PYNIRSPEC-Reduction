      SUBROUTINE ILSCHK ( FAIL, ERRMSG )
C
C VERSION
C     11-JUN-02  AD  Add tolerance for matching wavenumber boundaries
C     20-APR-01  AD  Return with FAIL=TRUE if no ILS fn found for spec.range
C     27-APR-00  AD  Get PT2ILS from ilscom.inc rather than local variable.
C     07-JUL-99  AD  Remove redundant variable ILSDEF
C     23-APR-99  AD  Compress error text to allow for C*8 LABSPC
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Assign ILS Functions to each spectral range 
C     Called once by INPCHK.
C
      IMPLICIT NONE 
C
C ARGUMENTS 
      LOGICAL          FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80     ERRMSG !  O  Error message written if FAIL is TRUE
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'ilscom.inc' ! Instrument Lineshape functions.
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
C 
C LOCAL VARIABLES
      INTEGER          IILS   ! Counter for ILS Functions
      INTEGER          ISPC   ! Counter for spectral ranges for RFM calcs.
      LOGICAL          DEFILS ! T=This is the default ILS function
      DOUBLE PRECISION DELWID ! Wide Mesh spacing for each spc.range
      DOUBLE PRECISION DELFIN ! Fine Mesh spacing for each spc.range
      DOUBLE PRECISION WNOTOL ! Tolerance [cm-1] for matching wavenumbers
C
C EXECUTABLE CODE ------------------------------------------------------------
C
      DO ISPC = 1, NSPC
C Calculation of DELWID and DELFIN need to match that in SPCWID, SPCFIN
        IF ( WNRSPC(ISPC) .GT. 1.0D0 ) THEN
          DELWID = 1.0D0
          DELFIN = 1.0D0 / (WNRSPC(ISPC)*NINT(NOMFIN/WNRSPC(ISPC)))
        ELSE
          DELWID = WNRSPC(ISPC) * NINT ( 1.0D0 / WNRSPC(ISPC) )
          DELFIN = WNRSPC(ISPC) / NINT ( NOMFIN * WNRSPC(ISPC) )
        END IF
        ILSSPC(ISPC) = 0
        WNOTOL = 0.1D0 * DELFIN
        DO IILS = 1, NILS 
          DEFILS = WNLILS(IILS) .EQ. 0.0D0 .AND. 
     &             WNUILS(IILS) .EQ. 0.0D0 
          IF ( ( WNLSPC(ISPC) .GE. WNLILS(IILS)-WNOTOL .AND.
     &           WNUSPC(ISPC) .LE. WNUILS(IILS)+WNOTOL ) 
     &                                           .OR. DEFILS ) THEN ! matches
            IF ( ABS ( PT1ILS(IILS) ) .GT. DELWID .OR.
     &           ABS ( PT2ILS(IILS) ) .GT. DELWID      ) THEN
              FAIL = .TRUE.
              ERRMSG = 'F-ILSCHK: ILS function extends wider '//
     &                 'than 1 widemesh interval from centre'
              RETURN
            END IF
            IF ( NINT ( ABS(PT1ILS(IILS))/DELFIN ) .GT. MAXILF .OR.
     &           NINT ( ABS(PT2ILS(IILS))/DELFIN ) .GT. MAXILF ) THEN
              FAIL = .TRUE.
              ERRMSG = 'F-ILSCHK: ILS function width > MAXILF '//
     &                 'in RFMSIZ.INC when mapped to fine grid'
              RETURN
            END IF
            IF ( ILSSPC(ISPC) .EQ. 0 ) THEN
              ILSSPC(ISPC) = IILS
            ELSE IF ( .NOT. DEFILS ) THEN
              FAIL = .TRUE.
              ERRMSG = 'F-ILSCHK: Ambiguous: at least 2 ILS Fns.'
     &          //' can be applied to Spc.range: '//LABSPC(ISPC)
              RETURN
            END IF
          END IF
        END DO
        IF ( ILSSPC(ISPC) .EQ. 0 ) THEN
          FAIL = .TRUE.
          ERRMSG = 'F-ILSCHK: No ILS Fn supplied to cover spectral'
     &             //' range: '//LABSPC(ISPC)
          RETURN
        END IF
      END DO
C
      FAIL = .FALSE.
      END
