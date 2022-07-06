      SUBROUTINE PTBATM ( IJAC )
C
C VERSION
C     03-MAY-05  AD  Add jdxcon.inc
C     20-JAN-03  AD  Add CHGATM if perturbing profile
C     14-DEC-01  AD  Correct aerosol for p,T jacobians. Call ATMAUX.
C     22-MAR-01  AD  Set CHGGAS if perturbing profile
C     30-DEC-99  AD  Add perturbation to horiz.gradient profiles as well
C     09-AUG-99  AD  Initialse and SAVE ILOW, IUPP and SAVATM.
C     31-JUL-99  AD  Remove unused variable MGAS
C     16-JUN-99  AD  Special case for Aerosol perturbations
C     04-FEB-99  AD  Changed definitions of perturbations.
C     01-JAN-99  AD  Original.
C
C DESCRIPTION
C     Perturb/unperturb atmospheric profiles for Jacobian calc.
C     Called by RFMPTB for each Jacobian element in turn if JAC option enabled.
C     Can also be called by RADJAC for each Jacobian temperature element in
C     turn if both JAC and BFX flags are enabled.
C     Routine starts by undoing previous perturbation (except for IJAC=1,
C     which is the first call)>
C     Last call should be with IJAC .GT. NJAC, in which case the only action
C     is to undo the previous perturbation (IJAC=NJAC).
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      IJAC   !  I  Index of Jacobian element (1:NJAC+1).
C
      EXTERNAL
     &  ATMAUX ! Set up auxiliary profiles of atmospheric parameters    
C
C GLOBAL CONSTANTS
      INCLUDE 'idxcon.inc' ! HITRAN molecular indices
      INCLUDE 'jdxcon.inc' ! Indices for non-VMR Jacobians
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
      INCLUDE 'rfmcon.inc' ! Constants for RFM calculations
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'gracom.inc' ! Atmospheric 2-D field data
      INCLUDE 'jaccom.inc' ! Jacobian data
C
C LOCAL CONSTANTS
      INTEGER  MAX2DP      ! Max. number of points in 2-D atmos.field
        PARAMETER ( MAX2DP = MAXATM * MAXGRA )
C
C LOCAL VARIABLES
      INTEGER  IATM   ! Atmospheric level counter
      INTEGER  IGAS   ! Index of species to be perturbed
      INTEGER  IGRA   ! Index of horizontal gradient profiles
      INTEGER  ILOW   ! Atmospheric level# just below perturbed levels
      INTEGER  IPTB   ! Atmospheric level# of Jacobian retrieval element
      INTEGER  IUPP   ! Atmospheric level# just above perturbed levels
      INTEGER  JGAS   ! Counter for all gases
      LOGICAL  LOWCOL ! True=all lower levels treated as column perturbation
      LOGICAL  PREJAC ! True=Jacobian for pressure element
      LOGICAL  TEMJAC ! True=Jacobian for temperature element
      LOGICAL  UPPCOL ! True=all upper levels treated as column perturbation
      REAL     FACTOR ! Scaling factor for perturbation at atmos.level
      REAL   SAVATM(MAXATM) ! Saved value of unperturbed profile
      REAL   SAVAER(MAXATM) ! Saved value of unperturbed aerosol
      REAL   SAVGRA(MAXATM,MAXGRA) ! Saved value of unptb 2-D field
      REAL   SAVAXG(MAXATM,MAXGRA) ! Saved value of unptd 2D aerosol field
C
C Data statements to satisfy code-checking programs - values irrelevant
      DATA IGAS / 0 /
      DATA ILOW / 0 /
      DATA IUPP / 0 /
      DATA TEMJAC / .FALSE. /
      DATA PREJAC / .FALSE. /
      DATA SAVATM / MAXATM * 0.0 /
      DATA SAVAER / MAXATM * 0.0 /
      DATA SAVGRA / MAX2DP * 0.0 /
      SAVE IGAS, ILOW, IUPP, SAVATM, SAVAER, TEMJAC, PREJAC
C
C EXECUTABLE CODE ------------------------------------------------------------
C
C Unless first call, copy back previous atmospheric perturbations
C using previously stored values of IUPP, ILOW, IGAS
      IF ( IJAC .GT. 1 ) THEN
        IF ( IGAS .LE. NGAS ) THEN
          DO IATM = ILOW+1, IUPP-1
            VMRATM(IATM,IGAS) = SAVATM(IATM)
            IF ( VMRATM(IATM,IGAS) .GT. ARGMIN ) THEN
              LNVATM(IATM,IGAS) = LOG ( VMRATM(IATM,IGAS) )
            ELSE
              LNVATM(IATM,IGAS) = LOGMIN
            END IF
            IF ( GRAFLG ) THEN
              DO IGRA = 1, NGRA
                VMRGRA(IATM,IGAS,IGRA) = SAVGRA(IATM,IGRA)
              END DO
            END IF
          END DO
        ELSE IF ( TEMJAC ) THEN          ! Temperature
          DO IATM = ILOW+1, IUPP-1
            TEMATM(IATM) = SAVATM(IATM)
            DNSATM(IATM) = PREATM(IATM) * AVOG * 1.0E-4 
     &                     / ( RGAS * TEMATM(IATM) )
            IF ( DNSATM(IATM) .GT. ARGMIN ) THEN
              LNDATM(IATM) = LOG ( DNSATM(IATM) )
            ELSE
              LNDATM(IATM) = LOGMIN
            END IF
            IF ( GRAFLG ) THEN
              DO IGRA = 1, NGRA
                TEMGRA(IATM,IGRA) = SAVGRA(IATM,IGRA)
                DNSGRA(IATM,IGRA) = PREGRA(IATM,IGRA) * AVOG * 1.0E-4
     &                              / ( RGAS * TEMGRA(IATM,IGRA) )
              END DO
            END IF
            IF ( IGSMOL(IDXAER) .NE. 0 ) THEN       ! Aerosol extinction
              JGAS = IGSMOL(IDXAER)
              VMRATM(IATM,JGAS) = SAVAER(IATM)
              IF ( VMRATM(IATM,JGAS) .GT. ARGMIN ) THEN
                LNVATM(IATM,JGAS) = LOG ( VMRATM(IATM,JGAS) )
              ELSE
                LNVATM(IATM,JGAS) = LOGMIN
              END IF
              IF ( GRAFLG ) THEN
                DO IGRA = 1, NGRA
                  VMRGRA(IATM,JGAS,IGRA) = SAVAXG(IATM,IGRA) 
                END DO
              END IF
            END IF
          END DO  
          CALL ATMAUX                    ! Restore refractive index profile
        ELSE IF ( PREJAC ) THEN          ! Pressure
          DO IATM = ILOW+1, IUPP-1
            PREATM(IATM) = SAVATM(IATM)
            IF ( PREATM(IATM) .GT. ARGMIN ) THEN
              LNPATM(IATM) = LOG ( PREATM(IATM) )
            ELSE
              LNPATM(IATM) = LOGMIN
            END IF
            DNSATM(IATM) = PREATM(IATM) * AVOG * 1.0E-4 
     &                     / ( RGAS * TEMATM(IATM) )
            IF ( DNSATM(IATM) .GT. ARGMIN ) THEN
              LNDATM(IATM) = LOG ( DNSATM(IATM) )
            ELSE
              LNDATM(IATM) = LOGMIN
            END IF
            IF ( GRAFLG ) THEN
              DO IGRA = 1, NGRA
                PREGRA(IATM,IGRA) = SAVGRA(IATM,IGRA)
                DNSGRA(IATM,IGRA) = PREGRA(IATM,IGRA) * AVOG * 1.0E-4
     &                              / ( RGAS * TEMGRA(IATM,IGRA) )
              END DO
            END IF
            IF ( IGSMOL(IDXAER) .NE. 0 ) THEN       ! Aerosol extinction
              JGAS = IGSMOL(IDXAER)
              VMRATM(IATM,JGAS) = SAVAER(IATM)
              IF ( VMRATM(IATM,JGAS) .GT. ARGMIN ) THEN
                LNVATM(IATM,JGAS) = LOG ( VMRATM(IATM,JGAS) )
              ELSE
                LNVATM(IATM,JGAS) = LOGMIN
              END IF
              IF ( GRAFLG ) THEN
                DO IGRA = 1, NGRA
                  VMRGRA(IATM,JGAS,IGRA) = SAVAXG(IATM,IGRA) 
                END DO
              END IF
            END IF
          END DO
          CALL ATMAUX                     ! Restore refractive index profile
        END IF
      END IF
      IF ( IJAC .GT. NJAC ) THEN
        IGAS = 0           ! Re-initialise for IJAC=1 if called again
        DO JGAS = 1, NGAS
          CHGGAS(JGAS) = .FALSE.
        END DO
        RETURN
      END IF
C Assume all default path calculations are stored prior to next perturbation
      DO IATM = 1, NATM
        CHGATM(IATM) = .FALSE.
      END DO
C
C Carry out next perturbation
C
      IPTB = IATJAC(IJAC)
      IGAS = IGSJAC(IJAC)
      TEMJAC = IGAS .EQ. MAXGAS + JDXTEM
      PREJAC = IGAS .EQ. MAXGAS + JDXPRE
      ILOW = ILOJAC(IJAC)
      IUPP = IUPJAC(IJAC)
      LOWCOL = ILOW .EQ. 0
      UPPCOL = IUPP .EQ. NATM+1 
      DO IATM = ILOW+1, IUPP-1
        FACTOR = 1.0
        IF ( IATM .LT. IPTB .AND. .NOT. LOWCOL ) THEN
          FACTOR = ( HGTATM(IATM) - HGTATM(ILOW) ) /
     &             ( HGTATM(IPTB) - HGTATM(ILOW) )
        ELSE IF ( IATM .GT. IPTB .AND. .NOT. UPPCOL ) THEN
          FACTOR = ( HGTATM(IUPP) - HGTATM(IATM) ) /
     &             ( HGTATM(IUPP) - HGTATM(IPTB) )
        END IF
        CHGATM(IATM) = .TRUE.           ! Flag layer IATM as perturbed
        IF ( IGAS .LE. NGAS ) THEN
          SAVATM(IATM) = VMRATM(IATM,IGAS)
          IF ( IDXGAS(IGAS) .EQ. IDXAER ) THEN
C factor 1.0e-5 converts /km to /cm 
C NB: this different to ATMAUX since VMR is now in ppv
            VMRATM(IATM,IGAS) = VMRATM(IATM,IGAS) + 
     &        FACTOR * PTBAER * 1.0E-5 / DNSATM(IATM)
          ELSE
            VMRATM(IATM,IGAS) = VMRATM(IATM,IGAS)*(1.0 + FACTOR*PTBVMR)
          END IF
          IF ( VMRATM(IATM,IGAS) .GT. ARGMIN ) 
     &      LNVATM(IATM,IGAS) = LOG ( VMRATM(IATM,IGAS) )
          IF ( GRAFLG ) THEN 
            DO IGRA = 1, NGRA 
              SAVGRA(IATM,IGRA) = VMRGRA(IATM,IGAS,IGRA) 
              IF ( IDXGAS(IGAS) .EQ. IDXAER ) THEN
                VMRGRA(IATM,IGAS,IGRA) = VMRGRA(IATM,IGAS,IGRA) +
     &            FACTOR * PTBAER * 1.0E-5 / DNSGRA(IATM,IGRA)
              ELSE
                VMRGRA(IATM,IGAS,IGRA) = VMRGRA(IATM,IGAS,IGRA) *
     &                                   ( 1.0 + FACTOR * PTBVMR )
              END IF
            END DO
          END IF
          DO JGAS = 1, NGAS
            CHGGAS(JGAS) = .FALSE.
          END DO
          CHGGAS(IGAS) = .TRUE.
        ELSE IF ( TEMJAC ) THEN
          SAVATM(IATM) = TEMATM(IATM)
          TEMATM(IATM) = TEMATM(IATM) + FACTOR*PTBTEM
          DNSATM(IATM) = PREATM(IATM) * AVOG * 1.0E-4 
     &                  / ( RGAS * TEMATM(IATM) )
          IF ( DNSATM(IATM) .GT. ARGMIN ) THEN
            LNDATM(IATM) = LOG ( DNSATM(IATM) )
          ELSE
            LNDATM(IATM) = LOGMIN
          END IF
          IF ( GRAFLG ) THEN
            DO IGRA = 1, NGRA
              SAVGRA(IATM,IGRA) = TEMGRA(IATM,IGRA)
              TEMGRA(IATM,IGRA) = TEMGRA(IATM,IGRA) + FACTOR*PTBTEM
              DNSGRA(IATM,IGRA) = PREGRA(IATM,IGRA) * AVOG * 1.0E-4
     &                            / ( RGAS * TEMGRA(IATM,IGRA) )
            END DO
          END IF
          IF ( IGSMOL(IDXAER) .NE. 0 ) THEN       ! Aerosol extinction
            JGAS = IGSMOL(IDXAER)
            SAVAER(IATM) = VMRATM(IATM,JGAS)
            VMRATM(IATM,JGAS) = VMRATM(IATM,JGAS) * 
     &                          TEMATM(IATM) / SAVATM(IATM)
            IF ( VMRATM(IATM,JGAS) .GT. ARGMIN ) 
     &        LNVATM(IATM,JGAS) = LOG ( VMRATM(IATM,JGAS) )
            IF ( GRAFLG ) THEN
              DO IGRA = 1, NGRA
                SAVAXG(IATM,IGRA) = VMRGRA(IATM,JGAS,IGRA) 
                VMRGRA(IATM,JGAS,IGRA) = VMRGRA(IATM,JGAS,IGRA) * 
     &                                   TEMATM(IATM) / SAVATM(IATM)
              END DO
            END IF
          END IF
          DO JGAS = 1, NGAS
            CHGGAS(JGAS) = .TRUE.
          END DO
          CALL ATMAUX                      ! Perturbed refractive index 
        ELSE IF ( PREJAC ) THEN
          SAVATM(IATM) = PREATM(IATM)
          PREATM(IATM) = PREATM(IATM) * ( 1.0 + FACTOR*PTBPRE )
          IF ( PREATM(IATM) .GT. ARGMIN ) 
     &      LNPATM(IATM) = LOG ( PREATM(IATM) ) 
          DNSATM(IATM) = PREATM(IATM) * AVOG * 1.0E-4 
     &                  / ( RGAS * TEMATM(IATM) )
          IF ( DNSATM(IATM) .GT. ARGMIN ) THEN
            LNDATM(IATM) = LOG ( DNSATM(IATM) )
          ELSE
            LNDATM(IATM) = LOGMIN
          END IF
          IF ( GRAFLG ) THEN
            DO IGRA = 1, NGRA
              SAVGRA(IATM,IGRA) = PREGRA(IATM,IGRA)
              PREGRA(IATM,IGRA) = PREGRA(IATM,IGRA)*(1.0+FACTOR*PTBPRE)
              DNSGRA(IATM,IGRA) = PREGRA(IATM,IGRA) * AVOG * 1.0E-4
     &                            / ( RGAS * TEMGRA(IATM,IGRA) )
            END DO
          END IF
          IF ( IGSMOL(IDXAER) .NE. 0 ) THEN       ! Aerosol extinction
            JGAS = IGSMOL(IDXAER)
            SAVAER(IATM) = VMRATM(IATM,JGAS)
            VMRATM(IATM,JGAS) = VMRATM(IATM,JGAS) * 
     &                          SAVATM(IATM) / PREATM(IATM)
            IF ( VMRATM(IATM,JGAS) .GT. ARGMIN ) 
     &        LNVATM(IATM,JGAS) = LOG ( VMRATM(IATM,JGAS) )
            IF ( GRAFLG ) THEN
              DO IGRA = 1, NGRA
                SAVAXG(IATM,IGRA) = VMRGRA(IATM,JGAS,IGRA) 
                VMRGRA(IATM,JGAS,IGRA) = VMRGRA(IATM,JGAS,IGRA) * 
     &                                   SAVATM(IATM) / PREATM(IATM)
              END DO
            END IF
          END IF
          DO JGAS = 1, NGAS
            CHGGAS(JGAS) = .TRUE.
          END DO
          CALL ATMAUX                       ! Perturbed refractive index 
        END IF
      END DO
C Also need to allow for the fact that if profile is changed at level ILOW+1
C then the calculation in *layer* ILOW will also be affected
      IF ( ILOW .GE. 1 ) CHGATM(ILOW) = .TRUE.
C
      END
