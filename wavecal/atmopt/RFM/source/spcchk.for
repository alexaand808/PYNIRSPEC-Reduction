      SUBROUTINE SPCCHK ( FAIL, ERRMSG )
C
C VERSION
C     24OCT13 AD Increase WNLMAX from 250000 to 50000
C     14SEP11 AD Add RJTFLG
C     04APR06 AD Change DFLOAT to DBLE
C     22MAR04 AD Bug#51: Correct indices in I-SPCCHK message
C     08AUG03 AD Allow for GHZ flag
C     05SEP01 AD Increase WNLMAX from 20000 to 25000.
C     18MAR01 AD Move checks on fine resln & SPCFIN to INPSPC
C                    Remove RMXSPC calculation
C     15MAR01 AD Don't call SPCFIN if FIN flag is set TRUE.
C                    If convolving with FIN flag FALSE, check ouput resolution
C                    is not finer than NOMFIN.
C     11AUG99 AD Make D.P. operations more explicit.
C     23APR99 AD Adapt error messages to allow for C*8 LABSPC
C     15DEC98 AD Set WNLMIN to 0.001cm-1 (instead of 0cm-1)
C     02JUL98 AD Remove MIPAS-specific limitations on range/resolution.
C                    Fix ERRMSG for non-integer value of No.pts/Wno.
C                    Adjust WNUSPC to ensure matches WNLSPC and WNRSPC
C     22JAN98 AD Remove redundant EXTERNAL declaration: TXTFLD
C     17DEC97 AD Set ILSSPC = 0 to define value for SPCWID
C     03NOV97 AD Set WNRMIN = 1.0/FLOAT(NPPMAX) rather than 1.D0/DFLOAT()
C     03MAR97 AD Version 3.
C     01OCT96 AD Version 2.
C     01OCT96 AD Correction: Use NINT(WNRSPC) in some error messages 
C     24SEP96 AD Shorten error message statements to limit to C*80
C     18SEP96 AD Add calculation of RMXSPC.
C                    Add OPTFLG to list of output type checks
C     01SEP96 AD Version 1.
C
C DESCRIPTION
C     Check tabulated spectral range and resolution data.
C     Called once by INPSPC following *SPC marker.
C     Requested parameters checked against values specified in RFM requiremnts.
C     Only checks ranges for which calculations have been specified. 
C     Also removes all ranges for which no calc. specified (ie resln=0).
C
      IMPLICIT NONE
C
C ARGUMENTS
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  RFMLOG ! Write text message to RFM log file.
     &, SPCWID ! Set Wide Mesh grid parameters from spectral range/resln.
C
C GLOBAL CONSTANTS
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'rfmcon.inc' ! Constants for RFM calculations
C
C LOCAL CONSTANTS
      DOUBLE PRECISION WNUMAX        ! Maximum Waveno.required for RFM
        PARAMETER (    WNUMAX = 50000.0D0 )    ! 50 000 allows for UV
      DOUBLE PRECISION WNLMIN        ! Minimum Waveno.required for RFM
        PARAMETER (    WNLMIN = 0.001D0 )
C
C COMMON VARIABLES
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
      INCLUDE 'widcom.inc' ! Wide mesh data
C
C LOCAL VARIABLES
      INTEGER ISPC         ! Counter for all tabulated spectral ranges
      INTEGER JSPC         ! Counter for spectral ranges to be calculated
      LOGICAL NOMORE       ! T=use only one unlabelled spectral range
      CHARACTER*31 ERRSTR  ! Common text component of error message
      CHARACTER*80 MESSGE  ! Info message for Log file
      CHARACTER*3  WNOGHZ  ! Either 'Wno' (GHZFLG=F) or 'GHz' (GHZFLG=T)
      DOUBLE PRECISION WNOFAC ! GHz to Wno conversion factor
      DOUBLE PRECISION WNORES  ! Wavenumber resolution
      DOUBLE PRECISION RPT ! (No.spectral points - 1) as a real number
C
C DATA STATEMENTS
      DATA NOMORE  / .FALSE. /
      SAVE NOMORE
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      JSPC = 0
      WMXSPC = WNLMIN      ! Ensures updated by any WNUSPC+FWIND value
      WMNSPC = WNUMAX      !   "                    WNLSPC-FWIND
      WMDSPC = 0.5D0 * ( WNUMAX + WNLMIN )
      IF ( GHZFLG ) THEN 
        WNOGHZ = 'GHz'
        WNOFAC = 1.0D7 / VLIGHT          ! VLIGHT in phycon.inc
      ELSE
        WNOGHZ = 'Wno'
        WNOFAC = 1.0D0
      END IF
C
      DO ISPC = 1, NSPC
        IF ( WNRSPC(ISPC) .NE. 0.0 ) THEN    ! This range selected for calc.
          FAIL = .TRUE.                      ! Set temporarily
          ERRSTR = 'F-SPCCHK: Range Label='//LABSPC(ISPC)//' '  ! C*31
          IF ( WNLSPC(ISPC) .GT. WNUSPC(ISPC) ) THEN
            WRITE ( ERRMSG, '(A,G12.5,A,G12.5)' ) ERRSTR //
     &        ' Lower ' // WNOGHZ // '=', WNLSPC(ISPC), 
     ^        ', > Upper ' // WNOGHZ // '=', WNUSPC(ISPC)
          ELSE IF ( WNLSPC(ISPC) .LT. WNLMIN ) THEN
            WRITE ( ERRMSG, '(A,G12.5,A,G12.5)' ) ERRSTR//' Lower ' //
     &        WNOGHZ // '=', WNLSPC(ISPC), ', < Minimum=', WNLMIN
          ELSE IF ( WNUSPC(ISPC) .GT. WNUMAX ) THEN
            WRITE ( ERRMSG, '(A,G12.5,A,G12.5)' ) ERRSTR//' Upper ' //
     &        WNOGHZ // '=', WNUSPC(ISPC), ', > Maximum=', WNUMAX
          ELSE IF ( WNRSPC(ISPC) .LT. 0.0 ) THEN
            WRITE ( ERRMSG, '(A,G12.5,A)' ) ERRSTR//
     &        ' Reqd Resln =', WNRSPC(ISPC), ' is negative'
          ELSE IF ( .NOT. GHZFLG .AND. WNRSPC(ISPC) .GT. 1.0 .AND. 
     &              MOD ( WNRSPC(ISPC), 1.D0 ) .NE. 0.0 ) THEN 
            WRITE ( ERRMSG, '(A,G12.5)' ) ERRSTR//
     &        ' non-integer for No.pts/Wno,=', WNRSPC(ISPC)
          ELSE 
            FAIL = .FALSE.
          END IF
          IF ( FAIL ) RETURN
C
C Check that upper wavenumber complies with lower wavenumber and resolution:
C adjust upper wavenumber if outside 0.01 of resolution and send warning.
C
          WNORES = WNRSPC(ISPC)
          IF ( .NOT. GHZFLG .AND. WNORES .GT. 1.D0 ) 
     &      WNORES = 1.D0 / WNORES
          RPT = ( WNUSPC(ISPC) - WNLSPC(ISPC) ) / WNORES 
c          IF ( ABS ( RPT - DFLOAT(NINT(RPT)) ) 
          IF ( ABS ( RPT - DBLE(NINT(RPT)) ) 
     &         .GT. 0.01D0 * WNORES            ) THEN
            WNUSPC(ISPC) = WNLSPC(ISPC) + NINT(RPT) * WNORES
            WRITE ( MESSGE, '(A,F12.6)' ) 'W-SPCCHK: Label='//
     &        LABSPC(ISPC)//' Adjusting upper ' // WNOGHZ // 
     &        ' boundary to ', WNUSPC(ISPC)
            CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
            IF ( FAIL ) RETURN
          END IF
C
C If using an unlabelled spectral range, it must be the only one used so check
C it is the first one and set flag to stop any more
C
          IF ( NOMORE .OR. 
     &         ( LABSPC(ISPC) .EQ. ' ' .AND. JSPC .NE. 0 ) )THEN
            FAIL = .TRUE.
            ERRMSG = 'F-CHKSPC: Attempted to use unlabelled spec.range '
     &        //'within other spectral ranges'
            RETURN
          END IF
          IF ( LABSPC(ISPC) .EQ. ' ' ) NOMORE = .TRUE.
C
C Basic range, resolution parameters OK
C
C Send log message indicating current range identified
C
          WRITE ( MESSGE, '(A,3F12.5)' ) 'I-SPCCHK: Label='
     &      //LABSPC(ISPC)//' Range, Resln=', WNLSPC(ISPC), 
     &      WNUSPC(ISPC), WNRSPC(ISPC)
          CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
C
C Reindex as JSPC, convert from GHz to Wno if necessary (WNOFAC), 
C derive other params and check against wide and fine mesh array sizes
C
          JSPC = JSPC + 1
          LABSPC(JSPC) = LABSPC(ISPC)            
          WNLSPC(JSPC) = WNLSPC(ISPC) * WNOFAC
          WNUSPC(JSPC) = WNUSPC(ISPC) * WNOFAC
          WNRSPC(JSPC) = WNRSPC(ISPC) * WNOFAC 
          ILSSPC(JSPC) = 0  ! reqd.for SPCWID, overwritten when ILS data loaded
C
C Update Min/Max wavenumbers required for any spectral range
          WMNSPC = MIN ( WMNSPC, WNLSPC(JSPC) )
          WMXSPC = MAX ( WMXSPC, WNUSPC(JSPC) )
C
C Determine no. wide mesh intervals required and check array sizes
C
          CALL SPCWID ( JSPC )
          IF ( NWID .GT. MAXWID ) THEN
            FAIL = .TRUE.
            WRITE ( ERRMSG, '(A,I6,A,I6)' )
     &        ERRSTR//' No.Wide Mesh intvls rqd=', NWID, 
     &        ' >MAXWID=', MAXWID
            RETURN
          END IF
        END IF
      END DO
      NSPC = JSPC
      WMDSPC = 0.5D0 * ( WMXSPC + WMNSPC )
C
C If no spectral range specified, note this in Log file (Warning status if
C spectral output also requested, ie ABS, BBT, OPT, RAD, RJT, TRA or WID flags).
C
      IF ( NSPC .EQ. 0 ) THEN
        IF ( ABSFLG .OR. BBTFLG .OR. OPTFLG .OR. 
     &       RADFLG .OR. RJTFLG .OR. TRAFLG .OR. WIDFLG ) THEN
          MESSGE = 'W-SPCCHK: No spectral calculations specified'
        ELSE
          MESSGE = 'I-SPCCHK: No spectral calculations specified'
        END IF
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
C
      END
