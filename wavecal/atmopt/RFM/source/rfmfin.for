      SUBROUTINE RFMFIN ( IWID, FAIL, ERRMSG  )
C
C VERSION
C     20-JUL-05  AD  Only call MIXSHP for CO2 lines
C     13-MAR-04  AD  Set SUBWNG according to default isotope
C     25-FEB-04  AD  MORSE: Add check for CLCPTH
C     17-FEB-04  AD  Accumulate NLNPTH
C     11-AUG-03  AD  Add VVWSHP, VVWCOR
C     01-MAY-03  AD  Set check .LT. WNUFIN instead of .LE. WNUFIN
C     12-JUN-02  AD  Simplify: assume only one HITRAN file
C     16-FEB-02  AD  Allow for isotopes  
C     04-NOV-99  AD  Check for CTMGAS rather than CTMFLG before SUBWNG
C     24-JAN-99  AD  Remove GL2FLG
C     01-DEC-97  AD  Add test for SHPCTM
C     18-OCT-97  AD  Change VAX "TYPE" to standard FORTRAN "WRITE" 
C     03-JUL-97  AD  Add test for SHPLUT
C     03-MAR-97  AD  Version 3.
C     29-JAN-97  AD  Debug - check for IBUF out of range
C     18-DEC-96  AD  Subtract offsets if CKD H2O Continuum applied
C     01-OCT-96  AD  Version 2.
C     20-SEP-96  AD  Rename VIBFLG to NTEFLG
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Perform a fine pass over the wavenumber grid.
C     Called by RFM for each wide-mesh interval.
C
      IMPLICIT NONE
C                  
C ARGUMENTS
      INTEGER      IWID               !  I  Wide Mesh#
      LOGICAL      FAIL               !  O  True if a fatal error is detected
      CHARACTER*80 ERRMSG             !  O  Error message written if FAIL
C
      EXTERNAL
     &  ADJUST ! Adjust line parameters for path conditions
     &, AWIDTH ! Calculate average half-width of line.
     &, CHISHP ! Calculate Voigt Line shape allowing for chi-factor.
     &, DOPSHP ! Calculate Doppler Line shape.
     &, INIHFL ! Initialise the line data files
     &, LORSHP ! Calculate Lorentz Line shape
     &, MIXSHP ! Calculate Voigt Line shape allowing for line-mixing
     &, REACYC ! Read HITRAN data into cyclic buffers.
     &, SUB25W ! Subtract value at 25cm-1 from abs.coeff.
     &, UNLCYC ! Unload required line from cyclic buffer.
     &, VOISHP ! Calculate Voigt Line shape
     &, VVWCOR ! Apply Van Vleck-Weisskopf correction to line shape
     &, VVWSHP ! Calculate Van Vleck-Weisskopf Line shape.
      REAL AWIDTH
C
C COMMON CONSTANTS
      INCLUDE 'idxcon.inc' ! HITRAN molecular indices
      INCLUDE 'rfmcon.inc' ! Constants for RFM calculations
      INCLUDE 'rfmsiz.inc' ! MIPAS RFM Array sizes
      INCLUDE 'shpcon.inc' ! Line Shape codes
C
C COMMON VARIABLES
      INCLUDE 'adjcom.inc' ! Path-Adjusted line data
      INCLUDE 'cyccom.inc' ! Cyclic line data buffers
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'pthcom.inc' ! RFM path data
      INCLUDE 'rejcom.inc' ! Minimum Line strength limits
C
C LOCAL VARIABLES
      INTEGER  ICYC        ! Cyclic buffer index
      INTEGER  IDGPTH      ! HITRAN Gas ID for path
      INTEGER  IFIN        ! Fine mesh counter
      INTEGER  IFIN1,IFIN2 ! Limits of fine mesh grid for lineshape calc.
      INTEGER  IFINC       ! Fine mesh grid pt closest to line centre
      INTEGER  IGAS        ! Absorber counter
      INTEGER  ILIN        ! Line counter
      INTEGER  IPTH        ! Path counter
      INTEGER  ISHP        ! Lineshape code for gas in path
      INTEGER  NFINC       ! No. fine mesh grid pts covering line wing
      INTEGER  NREQ        ! No.fine mesh grid points reqd for lineshape calc.
      REAL     AVWDTH      ! Average line width [/cm] for gas,path
      REAL     ABS(0:MAXFIN) ! Accumulated line absorption for current path
      REAL     ABSLIN(0:MAXFIN) ! Single line absorption for current path
      REAL     CNT(0:MAXFIN) ! Accumulated non-lte factors for current path
      LOGICAL  SUBWNG      ! T = Subtract abs.coeff at 25cm-1
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Set DP frequency parameters for this interval
C    
      WNLFIN  = WN1FIN - FEXC 
      WNUFIN  = WN2FIN + FEXC 
      IFIN1 = 1
      IFIN2 = NFIN
      NREQ = NFIN 
C       
C Read lines between WNLFIN and WNUFIN into cyclic buffer
C
      IF ( IWID .EQ. 1 ) THEN
        NLIN = 0                ! Clear cyclic buffer
        ICYC1 = 1
        CALL INIHFL ( WNLFIN, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
      CALL REACYC ( FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
C Loop over paths
C
      DO IPTH = 1, NCLC
        IF ( .NOT. CLCPTH(IPTH) ) GOTO 200
        IGAS = IGSPTH(IPTH)
        IDGPTH = IDXGAS(IGAS)
        ISHP = SHPGAS(IGAS)
        IF ( ISHP .EQ. SHPXSC .OR. ISHP .EQ. SHPCTM .OR.
     &       ISHP .GE. SHPLUT ) GOTO 200
        SUBWNG = CTMGAS(ISOGAS(0,IGAS)) .AND. ( IDGPTH .EQ. IDXH2O ) 
        DO IFIN = 1, NFIN
          ABS(IFIN) = 0.0
          CNT(IFIN) = 0.0
        END DO
C
        AVWDTH = AWIDTH ( IPTH, WN1FIN )  
C
C Loop over lines stored in buffer 
C
        DO ILIN = 1, NLIN
          ICYC = MOD ( ILIN + ICYC1 - 2, MAXCYC ) + 1
C
C This is to check for out-of-range values of ICYC which may result from the
C DEC Alpha optimisation. If this happens the solution appears to be to compile
C this particular module with /NOOPT.
          if ( ICYC .le. 0 .or. ICYC .gt. MAXCYC ) 
     &      write (*,*) 'W-RFMFIN: ICYC out of range, =', ICYC
C
          IF ( IDGCYC(ICYC) .EQ. IDGPTH .AND. 
     &         ISOGAS(ISOCYC(ICYC),IGAS) .EQ. IGAS .AND. 
     &         WNOCYC(ICYC) .LT. WNUFIN ) THEN
C NB < WNUFIN allows for any line sitting exactly on a widemesh boundary
C at WNUFIN will already be included in widemesh calc so .exclude from finemesh
            NLNPTH(IPTH) = NLNPTH(IPTH) + 1
            CALL UNLCYC ( ICYC )
            CALL ADJUST ( IPTH )
C
            IF ( STRADJ .GE. FINREJ ) THEN        ! Line not too weak
              IF ( WNGFLG ) THEN
                IFINC  = NINT ( ( WNOCYC(ICYC) - WN1FIN ) / WNRFIN )
                NFINC  = NINT ( NWDCUT * AVWDTH / SNGL ( WNRFIN ) )
                IFIN1  = MAX ( 1, IFINC - NFINC )
                IFIN2  = MIN ( NFIN, IFINC + NFINC )
                NREQ   = IFIN2 - IFIN1 + 1
                IF ( NREQ .LT. 1 ) GOTO 100
              ENDIF
C
C Loop over finemesh wavenumbers (just need IFIN1,IFIN2 really)
C
              IF ( ISHP .EQ. SHPVOI ) THEN
                IF ( MIXFLG .AND. IDGPTH .EQ. IDXCO2 ) THEN 
                  CALL MIXSHP ( NREQ, WNOFIN(IFIN1), ABSLIN(IFIN1) )
                ELSE
                  CALL VOISHP ( NREQ, WNOFIN(IFIN1), ABSLIN(IFIN1) )
                END IF
              ELSE IF ( ISHP .EQ. SHPLOR ) THEN
                CALL LORSHP ( NREQ, WNOFIN(IFIN1), ABSLIN(IFIN1) )
              ELSE IF ( ISHP .EQ. SHPDOP ) THEN
                CALL DOPSHP ( NREQ, WNOFIN(IFIN1), ABSLIN(IFIN1) )
              ELSE IF ( ISHP .EQ. SHPCHI ) THEN
                CALL CHISHP ( IPTH, NREQ, WNOFIN(IFIN1), ABSLIN(IFIN1) )
              ELSE IF ( ISHP .EQ. SHPVVW ) THEN
                CALL VVWSHP ( NREQ, WNOFIN(IFIN1), ABSLIN(IFIN1) )
              ELSE 
                STOP 'F-RFMFIN: Unrecognised lineshape value'
              ENDIF
              IF ( VVWFLG .AND. ISHP .NE. SHPVVW ) 
     &          CALL VVWCOR ( NREQ, WNOFIN(IFIN1), ABSLIN(IFIN1) )
C
C If applying CKD H2O continuum need to subtract abs.coeff at 25cm
              IF ( SUBWNG ) 
     &          CALL SUB25W ( NREQ, ABSLIN(IFIN1) )
C
C NB without non-LTE, ANTADJ=CNTADJ=1
C
              IF ( NTEFLG ) THEN
                DO IFIN = IFIN1, IFIN2
                  ABS(IFIN) = ABS(IFIN) + 
     &              SNGL ( ANTADJ * DBLE ( ABSLIN(IFIN) ) )
                  CNT(IFIN) = CNT(IFIN) + 
     &              SNGL ( CNTADJ * DBLE ( ABSLIN(IFIN) ) )
                END DO
              ELSE
                DO IFIN = IFIN1, IFIN2
                  ABS(IFIN) = ABS(IFIN) + ABSLIN(IFIN)
                END DO
              END IF
C
            END IF
          END IF
  100     CONTINUE
        END DO
        DO IFIN = 1, NFIN
          ABSFIN(IFIN,IPTH) = ABSFIN(IFIN,IPTH) + ABS(IFIN)
          CNTFIN(IFIN,IPTH) = CNTFIN(IFIN,IPTH) + CNT(IFIN)
        END DO
C
  200   CONTINUE
      END DO
C
      END
