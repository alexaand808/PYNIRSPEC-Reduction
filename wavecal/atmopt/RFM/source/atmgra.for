      SUBROUTINE ATMGRA
C
C VERSION
C     26-DEC-99  AD  Original
C
C DESCRIPTION
C     Interpolate to fill 2-D atmospheric profile field.
C     Called once by INPATM after all profiles are loaded
C     For each species (including pressure and temperature) the procedure is
C     to establish the psi locations at which profile values have been supplied
C     and interpolate values for any missing internal locations, and duplicate
C     end values for any missing external locations.
C     Horizontal interpolation is linear for T, log for p (linear if p small),
C     and depends upon LIN flag for VMR as in ATMLAY. 
C     The end result should be that full profiles for p, T and VMR are set at
C     all specified psi values.
C     Also loads reference profile (psi=0) into atmcom.inc
C
      IMPLICIT NONE
C
      EXTERNAL
     &  MOVGRA ! Exchange profiles between /ATMCOM/ and /GRACOM/ arrays.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'gracom.inc' ! Atmospheric 2-D field data
C
C LOCAL VARIABLES
      INTEGER IATM ! Profile level counter
      INTEGER IGRA ! Profile location counter
      INTEGER IPRF ! Species profile counter (-1=pre,0=tem)
      INTEGER JGRA ! Index for interpolated profiles
      INTEGER LGRA ! Previous value of IGRA (or 0 if no previous value)
      REAL    AI   ! Interpolation weight for IGRA
      REAL    AL   ! Interpolation weight for LGRA
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      DO IPRF = -1, NGAS
        LGRA = 0
        DO IGRA = 1, NGRA
          IF ( ( IPRF .EQ. -1 .AND. FLPGRA(IGRA) ) .OR.        ! Profile exists
     &         ( IPRF .EQ.  0 .AND. FLTGRA(IGRA) ) .OR.
     &         ( IPRF .GT.  0 .AND. FLVGRA(MAX(IPRF,1),IGRA) ) ) THEN
            IF ( LGRA .EQ. 0 ) THEN      ! No previous (=lower psi) profile 
              DO JGRA = 1, IGRA-1        ! So just copy this to all prev.
                DO IATM = 1, NATM
                  IF ( IPRF .EQ. -1 ) THEN
                    PREGRA(IATM,JGRA) = PREGRA(IATM,IGRA)
                  ELSE IF ( IPRF .EQ. 0 ) THEN
                    TEMGRA(IATM,JGRA) = TEMGRA(IATM,IGRA)
                  ELSE
                    VMRGRA(IATM,IPRF,JGRA) = VMRGRA(IATM,IPRF,IGRA)
                  END IF
                END DO
              END DO
            ELSE                         ! Previous value also, so interpolate
              DO JGRA = LGRA+1, IGRA-1
                AL = ( PSIGRA(IGRA) - PSIGRA(JGRA) ) / 
     &               ( PSIGRA(IGRA) - PSIGRA(LGRA) )
                AI = 1.0 - AL 
                DO IATM = 1, NATM
                  IF ( IPRF .EQ. -1 ) THEN
                    IF ( PREGRA(IATM,IGRA) .LT. ARGMIN .OR.
     &                   PREGRA(IATM,LGRA) .LT. ARGMIN      ) THEN
                      PREGRA(IATM,JGRA) = AI * PREGRA(IATM,IGRA) +
     &                                    AL * PREGRA(IATM,LGRA)
                    ELSE
                      PREGRA(IATM,JGRA) = 
     &                EXP ( AI * LOG ( PREGRA(IATM,IGRA) ) +
     &                      AL * LOG ( PREGRA(IATM,LGRA) )   )
                    END IF
                  ELSE IF ( IPRF .EQ. 0 ) THEN
                    TEMGRA(IATM,JGRA) = AI * TEMGRA(IATM,IGRA) + 
     &                                  AL * TEMGRA(IATM,LGRA)
                  ELSE
                    IF ( LINFLG .OR. 
     &                   VMRGRA(IATM,IPRF,IGRA) .LT. ARGMIN .OR.
     &                   VMRGRA(IATM,IPRF,LGRA) .LT. ARGMIN      ) THEN
                      VMRGRA(IATM,IPRF,JGRA) = 
     &                  AI * VMRGRA(IATM,IPRF,IGRA) +
     &                  AL * VMRGRA(IATM,IPRF,LGRA) 
                    ELSE
                      VMRGRA(IATM,IPRF,JGRA) = 
     &                  EXP ( AI * LOG ( VMRGRA(IATM,IPRF,IGRA) ) +
     &                        AL * LOG ( VMRGRA(IATM,IPRF,LGRA) )   )
                    END IF
                  END IF
                END DO
              END DO
            END IF
            LGRA = IGRA                           ! Save this IGRA as next LGRA
          ELSE IF ( IGRA .EQ. NGRA ) THEN         ! No profile loaded at NGRA
            DO JGRA = LGRA+1, NGRA                ! So copy last LGRA onwards
              DO IATM = 1, NATM
                IF ( IPRF .EQ. -1 ) THEN
                  PREGRA(IATM,JGRA) = PREGRA(IATM,LGRA)
                ELSE IF ( IPRF .EQ. 0 ) THEN
                  TEMGRA(IATM,JGRA) = TEMGRA(IATM,LGRA)
                ELSE
                  VMRGRA(IATM,IPRF,JGRA) = VMRGRA(IATM,IPRF,LGRA)
                END IF
              END DO
            END DO
          END IF
        END DO           
      END DO
C
      CALL MOVGRA ( NOMGRA, 0 )
C
      END                     
