      SUBROUTINE ABSCHK ( FAIL, ERRMSG )
C
C VERSION
C     13-JUN-02  Original.
C
C DESCRIPTION
C     Check information available for all absorbers.
C     Called once by INPCHK.
C
      IMPLICIT NONE
C
C ARGUMENTS
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG  !  O  Error message written if FAIL is TRUE
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'shpcon.inc' ! Line Shape codes
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'hflcom.inc' ! HITRAN line data file
      INCLUDE 'lflcom.inc' ! Look-Up Table file data
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
      INCLUDE 'xflcom.inc' ! X-section data files
C
C LOCAL VARIABLES
      INTEGER IGAS         ! Absorber counter
      INTEGER ILFL         ! Look-Up Table counter
      INTEGER ISPC         ! Spectral range counter
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      DO IGAS = 1, NGAS
        ISPC = 0
C If LUTs are being used, check if information for any gas for all spectral
C ranges is supplied by LUTs, in which case no HITRAN or XSC data is reqd
C for that gas
        IF ( LUTFLG ) THEN 
          DO ILFL = 1, NLFL
            IF ( GASLFL(ILFL) .EQ. IGAS ) ISPC = ISPC + 1
          END DO
        END IF
        IF ( ISPC .NE. NSPC ) THEN        ! Need X/S or HITRAN data
          FAIL = .TRUE.
          IF ( SHPGAS(IGAS) .EQ. SHPXSC ) THEN        ! X/S species
            IF ( .NOT. GASXFL(IGAS) ) THEN 
              ERRMSG = 'F-ABSCHK: No X/S data for Gas:'//CODGAS(IGAS)
              RETURN
            END IF
          ELSE                                        ! HITRAN line species
            IF ( LUNHFL .EQ. 0 .OR.                   ! No HITRAN file
     &           IFPHFL(IGAS) .GE. IR2HFL ) THEN      ! No lines in file
              ERRMSG = 'F-ABSCHK: No HITRAN line data for Gas:'//
     &          CODGAS(IGAS)
              RETURN
            END IF
          END IF
        END IF
      END DO
C
C Normal exit
      FAIL = .FALSE.
      END
