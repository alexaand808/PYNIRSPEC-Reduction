      SUBROUTINE MOVGRA ( IGRA, JGRA )
C
C VERSION
C     31-DEC-99  AD  Original.
C
C DESCRIPTION
C     Exchange profiles between /ATMCOM/ and /GRACOM/ arrays.
C     Called by PRFGRA and ATMLAY.
C     An index IGRA or JGRA =0 means use profiles in /ATMCOM/
C     Other indices refer to position in /GRACOM/
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER IGRA  !  I  Old index of profile data
      INTEGER JGRA  !  I  New index of profile data
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'gracom.inc' ! Atmospheric 2-D field data
C
C LOCAL VARIABLES
      INTEGER IATM  ! Atmospheric profile level counter
      INTEGER IGAS  ! Absorber counter
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      IF ( JGRA .EQ. 0 ) THEN                      ! Move IGRA to ATM
        DO IATM = 1, NATM
          TEMATM(IATM) = TEMGRA(IATM,IGRA)
          PREATM(IATM) = PREGRA(IATM,IGRA)
          DO IGAS = 1, NGAS
            VMRATM(IATM,IGAS) = VMRGRA(IATM,IGAS,IGRA)
          END DO
        END DO
      ELSE IF ( IGRA .EQ. 0 ) THEN                 ! Move ATM to JGRA
        DO IATM = 1, NATM
          TEMGRA(IATM,JGRA) = TEMATM(IATM)
          PREGRA(IATM,JGRA) = PREATM(IATM)
          DO IGAS = 1, NGAS
            VMRGRA(IATM,IGAS,JGRA) = VMRATM(IATM,IGAS)
          END DO
        END DO
      ELSE                                         ! Move IGRA to JGRA
        PSIGRA(JGRA) = PSIGRA(IGRA) 
        DO IATM = 1, NATM
          TEMGRA(IATM,JGRA) = TEMGRA(IATM,IGRA)
          FLTGRA(JGRA) = FLTGRA(IGRA)
          PREGRA(IATM,JGRA) = PREGRA(IATM,IGRA)
          FLPGRA(JGRA) = FLPGRA(IGRA)
          DO IGAS = 1, NGAS
            VMRGRA(IATM,IGAS,JGRA) = VMRGRA(IATM,IGAS,IGRA)
            FLVGRA(IGAS,JGRA) = FLVGRA(IGAS,IGRA)
          END DO
        END DO
      END IF
C
      END
