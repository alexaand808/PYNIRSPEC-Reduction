      SUBROUTINE MIXSHP ( NWNO, DWNO, ABSORP )
C
C VERSION
C     18-NOV-04  AD  Replace global constant MAXPNT with local MAXWNO
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Calculate Voigt Line shape allowing for line mixing.
C     Called by RFMFIN and RFMWID if GL2 or MIX Flag enabled.
C     (Uses the same algorithm as GENLN2 for mixed and non-mixed case).
C     The Voigt lineshape formulation:
C
C                 g(X,Y) = S * g0 * K(X,Y)
C                 g0 = 1/Ad * SQRT(ln2/pi)
C                 X = (nu - nu0)/Ad *SQRT(ln2)
C                 Y = Al/Ad *SQRT(ln2)
C                 K(X,Y) = Y/pi * 
C                 INT^(+infty)_(-infty){exp(-t**2)/[Y**2 + (X-t)**2]}dt
C
C     This routine calculates the complex probability function using a 
C     vectorized version of the Humlicek JQSRT V27 437 1982 paper. 
C     The calculation is performed for the array of x,y pairs for a given line 
C     over the fine mesh points of the current wide mesh. 
C
C     Uses path-adjusted line data in /ADJCOM/
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER          NWNO          !  I  No. wavenumber points to be evaluted
      DOUBLE PRECISION DWNO(NWNO)    !  I  Array of Wavenumbers [/cm]
      REAL             ABSORP(NWNO)  !  O  Absorption 
C
      EXTERNAL
     &  HUMLCK ! Calculate Humlicek complex prob.function for Voigt Line shape.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C LOCAL CONSTANTS
      INTEGER MAXWNO       ! Max expected size of NWNO parameter
        PARAMETER ( MAXWNO = MAXFIN + ( MAXWD2 - MAXFIN ) *
     &    ( 2 * MAXWD2 / ( MAXWD2 + MAXFIN ) ) )        ! =MAX(MAXWD2,MAXFIN)
C
C COMMON VARIABLES
      INCLUDE 'adjcom.inc' ! Path-Adjusted line data
      INCLUDE 'flgcom.inc' ! RFM option flags
C
C LOCAL VARIABLES
      INTEGER  IWNO        ! Wavenumber array counter
      REAL     H0
      REAL     REPWID   
      REAL     X(MAXWNO)
      REAL     Y
      COMPLEX  V(MAXWNO) 
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      IF ( NWNO .GT. MAXWNO ) STOP 'F-MIXSHP: Logical error'
C
C  SQRT(ln 2) = 0.83255461, 1/SQRT(PI) = 0.5641895
C
      REPWID = 0.83255461 / DOPADJ
      H0 = REPWID * 0.5641895 * STRADJ
      Y = WIDADJ * REPWID
      DO IWNO = 1, NWNO 
        X(IWNO) =  SNGL ( DWNO(IWNO) - WNOADJ ) * REPWID
      END DO
      CALL HUMLCK ( NWNO, X, Y, V )
C
C Compute absorption due to Voigt line shape
C
      IF ( MIXFLG ) THEN
        DO IWNO  = 1, NWNO 
          ABSORP(IWNO) = H0 * ( REAL ( V(IWNO) ) + 
     &                          YMXADJ * AIMAG ( V(IWNO) ) )
        END DO
      ELSE
        DO IWNO  = 1, NWNO 
          ABSORP(IWNO) = H0 * REAL ( V(IWNO) )
        END DO
      END IF
C
      END                                                 
