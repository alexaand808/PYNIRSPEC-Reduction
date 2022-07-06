      REAL FUNCTION VALATM ( HGT, FNC )
C
C VERSION
C     29-DEC-99  AD  Original.
C
C DESCRIPTION
C     Function to interpolate value from atmospheric profiles.
C
      IMPLICIT NONE
C
C ARGUMENTS
      REAL        HGT  !  I  Altitude [km]
      CHARACTER*3 FNC  !  I  Function to be interpolated
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
C
C LOCAL VARIABLES
      INTEGER          IATM ! Lower vertical profile index
      INTEGER          JATM ! Upper vertical profile index
      REAL             DATM ! Interpolation weight for IATM
      REAL             EATM ! Interpolation weight for JATM
C
C DATA STATEMENTS
      DATA IATM / 1 / 
      SAVE IATM
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Update IATM index if required
 100  CONTINUE
      IF ( HGT .GT. HGTATM(MIN(IATM+1,NATM)) .AND. IATM .LT. NATM ) THEN
        IATM = IATM + 1
        GOTO 100
      ELSE IF ( HGT .LT. HGTATM(IATM) .AND. IATM .GT. 1 ) THEN
        IATM = IATM - 1
        GOTO 100
      END IF
      IF ( IATM .LT. NATM ) THEN
        JATM = IATM + 1
        EATM = ( HGT - HGTATM(IATM) ) / ( HGTATM(JATM) - HGTATM(IATM) ) 
      ELSE
        JATM = IATM
        EATM = 0.0
      END IF
      DATM = 1.0 - EATM
C
      IF ( FNC .EQ. 'RFR' ) THEN
        IF ( HGT .GT. HGTATM(NATM) ) THEN  ! Assume value=0 above top of atm.
          VALATM = 0.0
        ELSE IF ( RFRATM(JATM) .GT. ARGMIN ) THEN  ! Assume RFR(IATM)>RFR(JATM)
          VALATM = EXP ( DATM * LOG ( RFRATM(IATM) ) +
     &                   EATM * LOG ( RFRATM(JATM) )   )
        ELSE
          VALATM =  DATM * RFRATM(IATM) + EATM * RFRATM(JATM) 
        END IF
C
C Density Scale Height applies to each layer rather than level, so no interp.
      ELSE IF ( FNC .EQ. 'DSH' ) THEN
        VALATM = DSHATM(IATM)
C
      ELSE
        STOP 'F-VALATM: Logical Error'
      END IF
C
      END

