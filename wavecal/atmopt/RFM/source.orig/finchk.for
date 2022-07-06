      SUBROUTINE FINCHK ( FAIL, ERRMSG )
C
C VERSION
C     13-JAN-03  AD  Bug fix: check NGFL GE 1 if testing RESGFL
C     20-MAR-01  AD  Original.
C
C DESCRIPTION
C     Check fine mesh resolution.
C     Called once by INPCHK.
C
      IMPLICIT NONE
C
C ARGUMENTS
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  CHKRES ! Check tangent pt line-widths compatible with spectral resln.
     &, SPCFIN ! Set Fine Mesh grid parameters from spectral range/resln.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'rfmcon.inc' ! Constants for RFM calculations
C
C LOCAL CONSTANTS
      DOUBLE PRECISION WNOTOL   ! Wavenumber tolerance for comparing resoln.
        PARAMETER ( WNOTOL = 0.01D0 )
C
C COMMON VARIABLES
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'flgcom.inc' ! RFM option flags      
      INCLUDE 'gflcom.inc' ! Grid file data
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
C
C LOCAL VARIABLES
      INTEGER      ISPC    ! Counter for spectral ranges
      REAL         WNRMAX  ! Maximum fine mesh spacing
      CHARACTER*31 ERRSTR  ! Common text component of error message
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      WNRMAX = 0.0         ! Ensures updated by any WNRFIN value
C
      DO ISPC = 1, NSPC
        ERRSTR = 'F-FINCHK: Range Label='//LABSPC(ISPC)//' '  ! C*31
C
C If convolving, check output spacing larger than NOMFIN spacing
        IF ( ILSFLG .OR. AVGFLG ) THEN
          IF ( ( WNRSPC(ISPC) .GT. 1.0 .AND.
     &           NINT ( NOMFIN / WNRSPC(ISPC) ) .LT. 1 ) .OR.
     &         ( WNRSPC(ISPC) .LE. 1.0 .AND.
     &           NINT ( NOMFIN * WNRSPC(ISPC) ) .LT. 1 )      ) THEN
            FAIL = .TRUE.
            WRITE ( ERRMSG, '(A)' )
     &        ERRSTR//' output grid finer than fine mesh grid'
            RETURN
          END IF
        END IF
C
C Calculate number of fine grid points per widemesh interval
        CALL SPCFIN ( ISPC )
        IF ( NFIN .GT. MAXFIN ) THEN
          FAIL = .TRUE.
          WRITE ( ERRMSG, '(A,I6,A,I6)' )
     &      ERRSTR//' No.Fine Mesh intvls rqd=', NFIN, 
     &      ' >MAXFIN=', MAXFIN
          RETURN
        END IF
        WNRMAX = MAX ( WNRMAX, SNGL ( WNRFIN ) )
C If using irregular grids, spacing must be the same as WNRFIN 
C (which may have been modified by *FIN section)
        IF ( GRDFLG .AND. ( ILSFLG .OR. AVGFLG ) .AND. NGFL .GT. 1 
     &       .AND. ABS ( WNRFIN - RESGFL ) .GT. WNOTOL * WNRFIN ) THEN
          FAIL = .TRUE.
          ERRMSG = 'F-FINCHK: fine mesh resolution inconsistent '//
     &             'with irregular grid resolution'
          RETURN
        END IF
      END DO
C
      IF ( .NOT. TABFLG .AND. .NOT. NADFLG .AND. .NOT. ZENFLG .AND. 
     &     .NOT. FLXFLG ) THEN 
        CALL CHKRES ( WNRMAX, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
C
      END
