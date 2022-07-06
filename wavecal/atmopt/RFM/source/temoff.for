      REAL FUNCTION TEMOFF ( LNP )
C
C VERSION
C     29-APR-FEB-99  AD  Original.
C
C DESCRIPTION
C     Temperature offset for pretabulated (p,T) data.
C     Called by RFMLUT, RFMTAB and TABPTH. 
C
      IMPLICIT NONE
C
C ARGUMENTS
      REAL      LNP            !  I  -ln(p), p is pressure [mb]
C
C LOCAL CONSTANTS
      INTEGER       NLEV       ! No.levels in temperature offset profile
        PARAMETER ( NLEV = 5 )
C
C LOCAL VARIABLES
      INTEGER   ILEV           ! Profile level counter
      REAL      ALPHA          ! Interpolation fraction
      REAL      TEMLEV(NLEV)   ! Profile temperature levels [K]
      REAL      LNPLEV(NLEV)   ! Profile -ln(p) levels [p in mb]
C
C DATA STATEMENTS
C Corresponding p[mb] = 1013.0,    250.0,     16.0, 1.0,    0.01 
        DATA LNPLEV / -6.92067, -5.52146, -2.77259, 0.0, 4.60517 /
        DATA TEMLEV / 280.0, 220.0, 210.0, 270.0, 200.0 /
      SAVE TEMLEV, LNPLEV
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      IF ( LNP .LE. LNPLEV(1) ) THEN
        TEMOFF = TEMLEV(1)
      ELSE IF ( LNP .GE. LNPLEV(NLEV) ) THEN
        TEMOFF = TEMLEV(NLEV)
      ELSE 
        DO ILEV = 2, NLEV
          IF ( LNP .GE. LNPLEV(ILEV) ) THEN 
            ALPHA = (LNP-LNPLEV(ILEV)) / (LNPLEV(ILEV-1)-LNPLEV(ILEV))
            TEMOFF = ALPHA*TEMLEV(ILEV-1) + (1.0-ALPHA)*TEMLEV(ILEV)
            RETURN
          END IF
        END DO
        STOP 'F-TEMOFF: Logical error'
      END IF
C
      END
