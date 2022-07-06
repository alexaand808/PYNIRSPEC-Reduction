      SUBROUTINE CHISHP ( IPTH, NWNO, DWNO, ABSORP )
C
C VERSION
C     18-NOV-04  AD  Replace global constant MAXPNT with local MAXWNO
C     03-APR-99  AD  Add SNGL ( ) to calc. of X
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-OCT-96  AD  Correction to inclusion of line-mixing
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Calculate Voigt Line shape allowing for chi-factor.
C     Called by RFMFIN and RFMWID if CHI shape selected.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER          IPTH          !  I  Path index
      INTEGER          NWNO          !  I  No. wavenumber points to be evaluted
      DOUBLE PRECISION DWNO(NWNO)    !  I  Array of Wavenumbers [/cm]
      REAL             ABSORP(NWNO)  !  O  Absorption 
C
      EXTERNAL
     &  CHICO2 ! Calculate Chi-Factor for CO2 sub-Lorentzian lineshape.
     &, HUMLCK ! Calculate Humlicek complex prob.function for Voigt Line shape.
C
C GLOBAL CONSTANTS
      INCLUDE 'idxcon.inc' ! HITRAN molecular indices
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C LOCAL CONSTANTS
      INTEGER MAXWNO       ! Max expected size of NWNO parameter
        PARAMETER ( MAXWNO = MAXFIN + ( MAXWD2 - MAXFIN ) *
     &    ( 2 * MAXWD2 / ( MAXWD2 + MAXFIN ) ) )        ! =MAX(MAXWD2,MAXFIN )
      REAL          CHILIM ! Inner limit for chi-factor only calculation
        PARAMETER ( CHILIM = 10.0 )
      REAL          MIXLIM ! Outer limit for line-mixing factor only calc.
        PARAMETER ( MIXLIM = 8.0 )
C
C COMMON VARIABLES
      INCLUDE 'adjcom.inc' ! Path-Adjusted line data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'hitcom.inc' ! HITRAN line database variables
C
C LOCAL VARIABLES
      INTEGER  IWNO        ! Wavenumber array counter
      REAL     CHIFAC(MAXWNO) ! Chi factor
      REAL     CHIMUL      ! Weighting for Chi-factor term
      REAL     DELWNO      ! Wno difference from line centre
      REAL     H0          ! REPWID * Line strength
      REAL     REPWID      ! sqrt(ln2) / Doppler width
      REAL     X(MAXWNO)   ! scaled wno from line centre (x-coord for Humlicek)
      REAL     Y           ! REPWID * Lorentz width (y-coord for Humlicek)
      COMPLEX  V(MAXWNO)   ! Complex Voigt vector
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      IF ( NWNO .GT. MAXWNO ) STOP 'F-CHISHP: Logical Error'
C
C Currently, Chi-factor only tabulated for CO2
C
      IF ( IDGAS .EQ. IDXCO2 ) THEN
        CALL CHICO2 ( IPTH, NWNO, DWNO, CHIFAC )
      ELSE
        DO IWNO = 1, NWNO
          CHIFAC(IWNO) = 1.0
        END DO
      END IF
C
C Calculate Voigt line shape
C
      REPWID = 0.83255461 / DOPADJ
      H0 = REPWID * 0.5641895 * STRADJ
      Y = WIDADJ * REPWID
      DO IWNO = 1, NWNO 
        X(IWNO) =  SNGL ( DWNO(IWNO) - WNOADJ ) * REPWID
      END DO
      CALL HUMLCK ( NWNO, X, Y, V )
C
      IF ( MIXFLG .AND. YMXADJ .NE. 0.0 ) THEN
C
C Compute absorption due to line mixing combined with chi-factor
C
        DO IWNO  = 1, NWNO 
          DELWNO = SNGL ( ABS ( DWNO(IWNO) - WNOADJ ) )
          IF ( DELWNO .LE. MIXLIM ) THEN    ! Line-mixing only
            ABSORP(IWNO) = H0 * ( REAL ( V(IWNO) ) + 
     &                     YMXADJ * AIMAG ( V(IWNO) ) )
          ELSE IF ( DELWNO .GE. CHILIM ) THEN  ! Chi-factor only
            ABSORP(IWNO) = H0 * REAL ( V(IWNO) ) * CHIFAC(IWNO) 
          ELSE
            CHIMUL =  ( DELWNO - MIXLIM ) / ( CHILIM - MIXLIM )           
            ABSORP(IWNO) = CHIMUL * H0 * CHIFAC(IWNO) * REAL ( V(IWNO))
     &           + ( 1.0 - CHIMUL ) * H0 * ( REAL ( V(IWNO) ) + 
     &                          YMXADJ * AIMAG ( V(IWNO) ) )
          END IF
        END DO
      ELSE                            
C
C Compute absorption due chi-factor alone
C
        DO IWNO  = 1, NWNO 
          ABSORP(IWNO) = H0 * REAL ( V(IWNO) ) * CHIFAC(IWNO) 
        END DO
      END IF
C
      END                                                 
