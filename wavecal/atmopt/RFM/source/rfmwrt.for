      SUBROUTINE RFMWRT ( IWID, FAIL, ERRMSG )
C
C VERSION
C     13-SEP-11  AD  Add RJT option
C     09-AUG-03  AD  Add BBT option
C     28-JAN-03  AD  Use IOFOUT, WNOOFF instead of calculating IDXOFF
C     05-JUN-00  AD  Allow for GRDFLG+TABFLG
C     22-APR-00  AD  Allow for COOFLG
C     04-JAN-00  AD  Allow for AVGFLG
C     11-AUG-99  AD  Replace general SAVE by specific SAVE of local variables
C     01-DEC-98  AD  Add FIRSTP check
C     27-NOV-98  AD  Simplify calculation of IOUT1,IOUT2 by using NPLEFT
C     09-NOV-98  AD  Correct to use INT rather than NINT for low res o/p.
C     08-AUG-98  AD  Correct for last Widemesh interval with ILS convolution
C     07-JUL-98  AD  Set output points relative to WNOFIN(1) not WNOWID(IWID-1)
C     23-OCT-97  AD  Add WRTTAB
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Write RFM output data for current widemesh interval.
C     Called by RFM for each widemesh interval.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      IWID    !  I  Current wide-mesh interval
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error occurs
      CHARACTER*80 ERRMSG  !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  WRTOUT ! Write RFM output spectra for current widemesh interval.
     &, WRTTAB ! Write RFM tabulated absorption coefficients.
     &, WRTWID ! Write RFM widemesh line count diagnostics.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'outcom.inc' ! RFM output file data
      INCLUDE 'widcom.inc' ! Wide mesh data
C
C LOCAL CONSTANTS
      DOUBLE PRECISION TOLFAC       ! Fraction of o/p Resln for comparing Wno.
        PARAMETER ( TOLFAC = 0.001D0 )
C
C LOCAL VARIABLES
      INTEGER NPLEFT ! No.points remaining to be written out
      LOGICAL FIRSTP ! T=look for first o/p point
      DOUBLE PRECISION WNOOFF ! Wno of first grid pt within intvl matching o/p 
C
      SAVE FIRSTP, NPLEFT
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      IF ( WIDFLG .AND. IWID .LE. NWID ) THEN
        CALL WRTWID ( IWID, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
C
      IF ( TABFLG .AND. GRDFLG ) THEN   ! 1,NFIN already set for output grid
        IOUT1 = 1
        IOUT2 = NFIN
        IF ( IWID .EQ. 1 ) THEN
          DO WHILE ( WNOFIN(IOUT1) .LT. WNLOUT - TOLFAC * WNROUT )
            IOUT1 = IOUT1 + 1
            IF ( IOUT1 .GT. NFIN ) STOP 'F-RFMWRT: Logical Error#1'
          END DO
        ELSE IF ( IWID .EQ. NWID ) THEN 
	  DO WHILE ( WNOFIN(IOUT2) .GT. WNUOUT + TOLFAC * WNROUT )
            IOUT2 = IOUT2 - 1
            IF ( IOUT2 .LT. 1 ) STOP 'F-RFMWRT: Logical Error#2'
          END DO
        END IF
        CALL WRTTAB ( FAIL, ERRMSG )
        RETURN
      ELSE IF ( IWID .EQ. 1 ) THEN
        NPLEFT = NPTOUT
        FIRSTP = .TRUE.
        IF ( ILSFLG .OR. AVGFLG ) RETURN ! FIRSTP wont be cancelled with IWID=1
        WNOOFF = WN1FIN + WNRFIN * IOFOUT
        IOUT1 = 1 + NINT ( ( WNLOUT - WNOOFF ) / WNROUT ) 
        IOUT2 = MIN ( NPWOUT, IOUT1 + NPLEFT - 1 )
      ELSE IF ( ( ILSFLG .OR. AVGFLG ) .AND. FIRSTP ) THEN         ! Set IOUT1
        WNOOFF = WN1FIN + WNRFIN * IOFOUT
        IF ( IWID .LE. NWID ) THEN
          IOUT1 = 1 + 
     &       NINT ( ( WNLOUT - WNOOFF + DELWID ) / WNROUT ) 
        ELSE                          ! If IWID=NWID+1 then FINGRD not updated 
          IOUT1 = 1 + NINT ( ( WNLOUT - WNOOFF ) / WNROUT ) 
        END IF
        IOUT2 = MIN ( NPWOUT, IOUT1 + NPLEFT - 1 )
      ELSE
        IOUT1 = 1
        IOUT2 = MIN ( NPWOUT, NPLEFT )
      END IF
      IF ( IOUT1 .GT. IOUT2 ) RETURN
      FIRSTP = .FALSE.                             ! Valid IOUT1 has been set 
      NPLEFT = NPLEFT - ( IOUT2 - IOUT1 + 1 )
C
      IF ( ABSFLG ) THEN
        CALL WRTOUT ( 'ABS', FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
      IF ( BBTFLG ) THEN
        CALL WRTOUT ( 'BBT', FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
      IF ( COOFLG ) THEN
        CALL WRTOUT ( 'COO', FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
      IF ( OPTFLG ) THEN
        CALL WRTOUT ( 'OPT', FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
      IF ( RADFLG ) THEN
        CALL WRTOUT ( 'RAD', FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
      IF ( RJTFLG ) THEN
        CALL WRTOUT ( 'RJT', FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
      IF ( TRAFLG ) THEN
        CALL WRTOUT ( 'TRA', FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
C
      IF ( TABFLG .AND. IWID .LE. NWID ) THEN
        CALL WRTTAB ( FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
C
      FAIL = .FALSE.
      END
