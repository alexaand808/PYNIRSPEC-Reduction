      SUBROUTINE PTBVIB ( IJAC )
C
C VERSION
C     03-MAY-05  AD  Original.
C
C DESCRIPTION
C     Perturb/unperturb VT profiles for Jacobian calc.
C     Called by ADJUST for each perturbed VT path element
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      IJAC    !  I  Index of VT Jacobian 
C
C GLOBAL CONSTANTS
      INCLUDE 'jdxcon.inc' ! Indices for non-VMR Jacobians
      INCLUDE 'rfmcon.inc' ! Constants for RFM calculations
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'jaccom.inc' ! Jacobian data
      INCLUDE 'ntecom.inc' ! Non-LTE data
C
C LOCAL VARIABLES
      INTEGER  IATM   ! Atmospheric level counter
      INTEGER  ILOW   ! Atmospheric level# just below perturbed levels
      INTEGER  INTE   ! Index of vibrational temperature profile
      INTEGER  IPTB   ! Atmospheric level# of Jacobian retrieval element
      INTEGER  IUPP   ! Atmospheric level# just above perturbed levels
      LOGICAL  LOWCOL ! True=all lower levels treated as column perturbation
      LOGICAL  UPPCOL ! True=all upper levels treated as column perturbation
      REAL     FACTOR ! Scaling factor for perturbation at atmos.level
      REAL     SAVNTE(MAXATM) ! Saved value of unperturbed profile
C
C Data statements to satisfy code-checking programs - values irrelevant
      DATA INTE / 0 /      
      DATA ILOW / 0 /
      DATA IUPP / 0 /
      DATA SAVNTE / MAXATM * 0.0 /
      SAVE INTE, ILOW, IUPP, SAVNTE
C
C EXECUTABLE CODE ------------------------------------------------------------
C
C Calls should alternate with saved INTE = 0 or > 0
      IF ( IJAC .GT. 0 ) THEN
        IF ( INTE .NE. 0 ) STOP 'F-PTBVIB: Logical Error#1'
        INTE = IGSJAC(IJAC) - ( MAXGAS + JDXNTE )
        IPTB = IATJAC(IJAC)
        ILOW = ILOJAC(IJAC)
        IUPP = IUPJAC(IJAC)
        LOWCOL = ILOW .EQ. 0
        UPPCOL = IUPP .EQ. NATM+1 
        DO IATM = ILOW+1, IUPP-1
          FACTOR = 1.0
          IF ( IATM .LT. IPTB .AND. .NOT. LOWCOL ) THEN
            FACTOR = ( HGTATM(IATM) - HGTATM(ILOW) ) /
     &               ( HGTATM(IPTB) - HGTATM(ILOW) )
          ELSE IF ( IATM .GT. IPTB .AND. .NOT. UPPCOL ) THEN
            FACTOR = ( HGTATM(IUPP) - HGTATM(IATM) ) /
     &               ( HGTATM(IUPP) - HGTATM(IPTB) )
          END IF
          SAVNTE(IATM) = TEMNTE(IATM,INTE)
          TEMNTE(IATM,INTE) = TEMNTE(IATM,INTE) + FACTOR*PTBNTE
        END DO
C
      ELSE IF ( IJAC .EQ. 0 ) THEN
        IF ( INTE .LE. 0 ) STOP 'F-PTBVIB: Logical Error#2'
        DO IATM = ILOW+1, IUPP-1
          TEMNTE(IATM,INTE) = SAVNTE(IATM)
        END DO
        INTE = 0
C
      ELSE         
        STOP 'F-PTBVIB: Logical Error#3'  ! INTE LE 0 should never happen
C
      END IF
C
      END        
