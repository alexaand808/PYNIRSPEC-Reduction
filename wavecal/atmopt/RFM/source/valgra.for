      REAL FUNCTION VALGRA ( HGT, PSI, FNC, IGAS )
C
C VERSION
C     23-JUL-03  AD  Use Refractivity instead of Density for 'DSH'
C     17-DEC-01  AD  Save DATM,EATM etc
C     24-MAY-01  AD  Correction: avoid extrapolation low psi values.
C     31-DEC-99  AD  Original.
C
C DESCRIPTION
C     Function interpolate value from 2D atmospheric field.
C
      IMPLICIT NONE
C
C ARGUMENTS
      REAL        HGT  !  I  Altitude [km]
      REAL        PSI  !  I  LOS angle [deg]
      CHARACTER*3 FNC  !  I  Function to be interpolated
      INTEGER     IGAS !  I  Gas# (defined if FNC = 'VMR' )
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gracom.inc' ! Atmospheric 2-D field data
C
C LOCAL VARIABLES
      INTEGER  IATM   ! Lower vertical profile index
      INTEGER  IGRA   ! Lower horizontal profile index
      INTEGER  JATM   ! Upper vertical profile index
      INTEGER  JGRA   ! Upper horizontal profile index
      REAL     DATM   ! Interpolation weight for IATM
      REAL     EATM   ! Interpolation weight for JATM
      REAL     DGRA   ! Interpolation weight for IGRA
      REAL     EGRA   ! Interpolation weight for JGRA
      REAL     HGTPRV ! Previous value of HGT
      REAL     PSIPRV ! Previous value of PSI       
C
C DATA STATEMENTS
      DATA IATM / 1 / 
      DATA IGRA / 1 /
      DATA JATM / 1 / 
      DATA JGRA / 1 /
      DATA DATM / 0.0 /
      DATA EATM / 0.0 /
      DATA DGRA / 0.0 /
      DATA EGRA / 0.0 /
      DATA HGTPRV / -1.0E20 /
      DATA PSIPRV / -1.0E20 / 
      SAVE IATM, IGRA, JATM, JGRA, DATM, EATM, DGRA, EGRA, 
     &     HGTPRV, PSIPRV
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Update IGRA index if required
      IF ( PSI .EQ. PSIPRV ) GOTO 200
 100  CONTINUE
      IF ( PSI .GT. PSIGRA(MIN(IGRA+1,NGRA)) .AND. IGRA .LT. NGRA ) THEN
        IGRA = IGRA + 1
        GOTO 100
      ELSE IF ( PSI .LT. PSIGRA(IGRA) .AND. IGRA .GT. 1 ) THEN
        IGRA = IGRA - 1
        GOTO 100
      END IF
      IF ( PSI .LT. PSIGRA(IGRA) ) THEN        ! Duplicate low psi value
        JGRA = 1
        EGRA = 0.0
      ELSE IF ( IGRA .EQ. NGRA ) THEN          ! Duplicate high psi value
        JGRA = NGRA
        EGRA = 0.0
      ELSE                                     ! Interpolate psi value
        JGRA = IGRA + 1
        EGRA = ( PSI - PSIGRA(IGRA) ) / ( PSIGRA(JGRA) - PSIGRA(IGRA) ) 
      END IF
      DGRA = 1.0 - EGRA
C
C Update IATM index if required
      IF ( HGT .EQ. HGTPRV ) GOTO 300
 200  CONTINUE
      IF ( HGT .GT. HGTATM(MIN(IATM+1,NATM)) .AND. IATM .LT. NATM ) THEN
        IATM = IATM + 1
        GOTO 200
      ELSE IF ( HGT .LT. HGTATM(IATM) .AND. IATM .GT. 1 ) THEN
        IATM = IATM - 1
        GOTO 200
      END IF
      IF ( IATM .LT. NATM ) THEN
        JATM = IATM + 1
        EATM = ( HGT - HGTATM(IATM) ) / ( HGTATM(JATM) - HGTATM(IATM) ) 
      ELSE
        EATM = 0.0
        JATM = IATM
      END IF
      DATM = 1.0 - EATM
C
  300 CONTINUE
      IF ( FNC .EQ. 'TEM' ) THEN
        VALGRA = DATM * DGRA * TEMGRA(IATM,IGRA) +
     &           DATM * EGRA * TEMGRA(IATM,JGRA) +
     &           EATM * DGRA * TEMGRA(JATM,IGRA) +
     &           EATM * EGRA * TEMGRA(JATM,JGRA)
C
      ELSE IF ( FNC .EQ. 'RFR' ) THEN
        IF ( MIN ( RFRGRA(JATM,IGRA), RFRGRA(JATM,JGRA ) ) 
     &       .GT. ARGMIN ) THEN         ! Assume RFR(IATM) > RFR(JATM)
          VALGRA = EXP ( DATM * DGRA * LOG ( RFRGRA(IATM,IGRA) ) +
     &                   DATM * EGRA * LOG ( RFRGRA(IATM,JGRA) ) +
     &                   EATM * DGRA * LOG ( RFRGRA(JATM,IGRA) ) +
     &                   EATM * EGRA * LOG ( RFRGRA(JATM,JGRA) )   )
        ELSE
          VALGRA =  DATM * DGRA * RFRGRA(IATM,IGRA) +
     &              DATM * EGRA * RFRGRA(IATM,JGRA) +
     &              EATM * DGRA * RFRGRA(JATM,IGRA) +
     &              EATM * EGRA * RFRGRA(JATM,JGRA) 
        END IF
C
      ELSE IF ( FNC .EQ. 'DSH' ) THEN              ! Density Scale Height (+ve)
        IF ( JATM .NE. IATM ) THEN
C 23JUL03: Use refractivity gradient if possible (see also ATMAUX calculation)
          IF ( MIN ( RFRGRA(JATM,IGRA), RFRGRA(JATM,JGRA) ) 
     &         .GT. ARGMIN*1.0E6 ) THEN
            VALGRA = ( HGTATM(JATM) - HGTATM(IATM) ) /
     &            ( DGRA * LOG(RFRGRA(IATM,IGRA)/RFRGRA(JATM,IGRA)) +
     &              EGRA * LOG(RFRGRA(IATM,JGRA)/RFRGRA(JATM,JGRA))   )
          ELSE IF ( MIN ( DNSGRA(JATM,IGRA), DNSGRA(JATM,JGRA) ) 
     &         .GT. ARGMIN ) THEN
            VALGRA = ( HGTATM(JATM) - HGTATM(IATM) ) /
     &            ( DGRA * LOG(DNSGRA(IATM,IGRA)/DNSGRA(JATM,IGRA)) +
     &              EGRA * LOG(DNSGRA(IATM,JGRA)/DNSGRA(JATM,JGRA))   )
          ELSE
            VALGRA = 1.0E20
          END IF
        ELSE
          VALGRA = 1.0E20
        END IF
C
      ELSE IF ( FNC .EQ. 'DRP' ) THEN                ! d(Refractivity)/d(Psi)
        IF ( IGRA .EQ. JGRA ) THEN                   ! No horizontal gradient
          VALGRA = 0.0
        ELSE                  ! Since horiz.grad is +/-, don't do log interp.
          VALGRA = ( DATM * ( RFRGRA(IATM,JGRA) - RFRGRA(IATM,IGRA) )+
     &               EATM * ( RFRGRA(JATM,JGRA) - RFRGRA(JATM,IGRA) ) )
     &             / ( PSIGRA(JGRA) - PSIGRA(IGRA) ) / DTORAD
        END IF
C
      ELSE IF ( FNC .EQ. 'VMR' ) THEN
        IF ( LINFLG .OR. 
     &       MIN ( VMRGRA(IATM,IGAS,IGRA), VMRGRA(IATM,IGAS,JGRA), 
     &             VMRGRA(JATM,IGAS,IGRA), VMRGRA(JATM,IGAS,JGRA) )
     &       .LT. ARGMIN ) THEN
          VALGRA =  DATM * DGRA * VMRGRA(IATM,IGAS,IGRA) +
     &              DATM * EGRA * VMRGRA(IATM,IGAS,JGRA) +
     &              EATM * DGRA * VMRGRA(JATM,IGAS,IGRA) +
     &              EATM * EGRA * VMRGRA(JATM,IGAS,JGRA) 
        ELSE
          VALGRA = EXP ( DATM * DGRA * LOG ( VMRGRA(IATM,IGAS,IGRA) ) +
     &                   DATM * EGRA * LOG ( VMRGRA(IATM,IGAS,JGRA) ) +
     &                   EATM * DGRA * LOG ( VMRGRA(JATM,IGAS,IGRA) ) +
     &                   EATM * EGRA * LOG ( VMRGRA(JATM,IGAS,JGRA) )  )
        END IF
C
      ELSE IF ( FNC .EQ. 'PRE' ) THEN
        IF ( MIN ( PREGRA(JATM,IGRA), PREGRA(JATM,JGRA) ) 
     &                                        .LT. ARGMIN ) THEN
          VALGRA =  DATM * DGRA * PREGRA(IATM,IGRA) +
     &              DATM * EGRA * PREGRA(IATM,JGRA) +
     &              EATM * DGRA * PREGRA(JATM,IGRA) +
     &              EATM * EGRA * PREGRA(JATM,JGRA) 
        ELSE
          VALGRA = EXP ( DATM * DGRA * LOG ( PREGRA(IATM,IGRA) ) +
     &                   DATM * EGRA * LOG ( PREGRA(IATM,JGRA) ) +
     &                   EATM * DGRA * LOG ( PREGRA(JATM,IGRA) ) +
     &                   EATM * EGRA * LOG ( PREGRA(JATM,JGRA) )   )
        END IF
C
      ELSE IF ( FNC .EQ. 'DNS' ) THEN
        IF ( MIN ( DNSGRA(JATM,IGRA), DNSGRA(JATM,JGRA) ) 
     &                                        .LT. ARGMIN ) THEN
          VALGRA =  DATM * DGRA * DNSGRA(IATM,IGRA) +
     &              DATM * EGRA * DNSGRA(IATM,JGRA) +
     &              EATM * DGRA * DNSGRA(JATM,IGRA) +
     &              EATM * EGRA * DNSGRA(JATM,JGRA) 
        ELSE
          VALGRA = EXP ( DATM * DGRA * LOG ( DNSGRA(IATM,IGRA) ) +
     &                   DATM * EGRA * LOG ( DNSGRA(IATM,JGRA) ) +
     &                   EATM * DGRA * LOG ( DNSGRA(JATM,IGRA) ) +
     &                   EATM * EGRA * LOG ( DNSGRA(JATM,JGRA) )   )
        END IF

      ELSE
        STOP 'F-VALGRA: Logical Error'
      END IF
C
      END

