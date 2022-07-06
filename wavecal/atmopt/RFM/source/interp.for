      SUBROUTINE INTERP ( NXD, NXG, XD, XG, YD, LOGINT, EXTRAP, YG )
C
C VERSION
C     12-JUN-04  AD  Ensure JGUESS is properly initialised and saved
C     20-DEC-00  AD  Don't extrapolate if only one point    
C     18-AUG-98  AD  Limit argument of LOG to a minimum value ARGMIN
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Interpolate data values YD(XD) to points YG(XG).
C     Called by ATMFIL.
C     If EXTRAP is set FALSE, then where XD<XG or XD>XG the last data values 
C       within XG are duplicated.
C     If EXTRAP is set TRUE, then outside points are extrapolated (if NXD>1)
C     If LOGINT is set TRUE, interpolation is linear in Log(YD). Note that this
C       also results in YD being temporarily converted to Log values
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER    NXD       !  I  No. of data points
      INTEGER    NXG       !  I  No. of grid points
      REAL       XD(NXD)   !  I  List of data coordinates
      REAL       XG(NXG)   !  I  List of grid coordinates
      REAL       YD(NXD)   ! I/O List of data points (see note on LOGINT)
      LOGICAL    LOGINT    !  I  TRUE=interpolate linearly in Log(YD)
      LOGICAL    EXTRAP    !  I  TRUE=extrapolate beyond ends of XD
      REAL       YG(NXG)   !  O  Interpolated/Extrapolated data points on grid
C
      EXTERNAL
     &  LOOKUP  !  General purpose interpolation routine 
C
C GLOBAL CONSTANTS
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
C
C LOCAL VARIABLES
      INTEGER    IXD          !  Counter for original data array
      INTEGER    JGUESS       !  Guess value to speed up LOOKUP
      INTEGER    IXG          !  Index for grid arrays
      REAL       XDMAX, XDMIN !  Max, Min values of data coordinates
      REAL       YDMAX, YDMIN !  Data values at end points XDMAX, XDMIN
C
C DATA STATEMENTS
      DATA JGUESS / 1 /
      SAVE JGUESS
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      IF ( LOGINT ) THEN
        DO IXD = 1, NXD
          YD(IXD) = LOG ( MAX ( ARGMIN, YD(IXD) ) )
        END DO
      END IF
C
      IF ( EXTRAP .AND. NXD .GT. 1 ) THEN       ! Can't extrapolate from 1 pt
        JGUESS = 1
        DO IXG = 1, NXG
          CALL LOOKUP ( XG(IXG), YG(IXG), XD, YD, NXD, JGUESS, 1 )
          IF ( LOGINT ) YG(IXG) = EXP ( YG(IXG) )
        END DO
      ELSE
        IF ( XD(NXD) .GT. XD(1) ) THEN         ! XD values ascending with index
          XDMAX = XD(NXD) 
          XDMIN = XD(1) 
          YDMAX = YD(NXD)
          YDMIN = YD(1)
        ELSE
          XDMIN = XD(NXD)
          XDMAX = XD(1)
          YDMIN = YD(NXD)
          YDMAX = YD(1)
        END IF
        DO IXG = 1, NXG
          IF ( XG(IXG) .LE. XDMIN ) THEN 
            YG(IXG) = YDMIN
            IF ( LOGINT ) YG(IXG) = EXP ( YG(IXG) )
          ELSE IF ( XG(IXG) .GE. XDMAX ) THEN
            YG(IXG) = YDMAX
            IF ( LOGINT ) YG(IXG) = EXP ( YG(IXG) )
          ELSE
            CALL LOOKUP ( XG(IXG), YG(IXG), XD, YD, NXD, JGUESS, 1 )
            IF ( LOGINT ) YG(IXG) = EXP ( YG(IXG) )
          END IF
        END DO
      END IF
C
      IF ( LOGINT ) THEN
        DO IXD = 1, NXD
          YD(IXD) = EXP ( YD(IXD) )
        END DO
      END IF

      END
