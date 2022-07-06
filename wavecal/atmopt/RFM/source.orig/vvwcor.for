      SUBROUTINE VVWCOR ( NWNO, DWNO, ABSORP )
C
C VERSION
C     11-SEP-03  AD  Change factor (v/v0) to (v/v0)^2
C     11-AUG-03  AD  New. Based on LORSHP  
C
C DESCRIPTION
C     Apply Van Vleck-Weisskopf correction to line shape.
C     Called by RFMFIN and RFMWID if VVW flag set true.
C     Not called if VVW lineshape already used for calculation.
C     Uses path-adjusted line data in /ADJCOM/
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER          NWNO          !  I  No. wavenumber points to be evaluted
      DOUBLE PRECISION DWNO(NWNO)    !  I  Array of Wavenumbers [/cm]
      REAL             ABSORP(NWNO)  ! I/O Absorption 
C
C GLOBAL CONSTANTS
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
C
C COMMON VARIABLES
      INCLUDE 'adjcom.inc' ! Path-Adjusted line data
C
C LOCAL VARIABLES
      INTEGER  IWNO ! Wavenumber array counter
      REAL     TOP  ! Numerator of Lorentz expression
      REAL     W2   ! Width^2
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      TOP = STRADJ * WIDADJ / PI                   
      W2 = WIDADJ**2
C
      DO IWNO = 1, NWNO 
        ABSORP(IWNO) = ( DWNO(IWNO) / WNOADJ )**2 * ( ABSORP(IWNO) + 
     &    TOP / ( SNGL ( ( DWNO(IWNO) + WNOADJ )**2 ) + W2 )   )
      END DO
C
      END                                                 
