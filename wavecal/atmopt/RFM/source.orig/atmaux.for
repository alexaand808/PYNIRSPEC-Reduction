      SUBROUTINE ATMAUX 
C
C VERSION
C     08-FEB-05  AD  Approximate Density ScHt if layers too narrow for diff.
C                    Set refractivity before calculating DSH
C     23-JUL-03  AD  Calculate DSH from refractivity rather than density
C     28-MAR-01  AD  Check IAXGAS before converting aerosol
C     05-OCT-00  AD  Bug fix: correct handling of Aerosol with GRA flag
C     31-DEC-99  AD  Also set RFR and DNS for 2-D fields
C     11-AUG-99  AD  Remove redundant flgcom.inc
C     06-APR-99  AD  Add flag to prevent repeated conversion of aer.units
C     24-JAN-99  AD  Replace RFRGL2, RFROFM with REFRAC
C                    Remove GL2FLG.
C     18-JUL-98  AD  Comment change only
C     18-OCT-97  AD  Change VAX "TYPE" to standard FORTRAN "WRITE"
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     28-SEP-96  AD  Add treatment of Aerosol "VMR". 
C                    Remove NGAS argument (enter via gascom.inc)
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Set up auxiliary profiles of atmospheric parameters
C     Called once by INPATM and possibly OBSATM 
C     These include Refractivity, Density, Log(Density), Log(Pressure) and
C     Log(VMR).
C
      IMPLICIT NONE
C
      EXTERNAL
     &  MOVGRA ! Exchange profiles between /ATMCOM/ and /GRACOM/ arrays.     
     &, REFRAC ! Calculate refractivity using Edlen formula
      REAL REFRAC
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
      INCLUDE 'idxcon.inc' ! HITRAN molecular indices
C
C LOCAL VARIABLES
      INTEGER  IATM   ! Surface counter for interpolated profiles
      INTEGER  IGAS   ! Gas counter
      INTEGER  IGRA   ! Gradient profile location counter
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      IF ( GRAFLG ) THEN
C Set Density and Refractivity for gradient profiles
C NB Since REFRAC requires profile to be loaded into /ATMCOM/ it is necessary
C to transfer profiles p,T,VMR first
        DO IGRA = 1, NGRA
          CALL MOVGRA ( IGRA, 0 )
          DO IATM = 1, NATM
            DNSGRA(IATM,IGRA) = PREATM(IATM) * AVOG * 1.0E-4 / 
     &                          RGAS / TEMATM(IATM)
            RFRGRA(IATM,IGRA) = REFRAC ( IATM )
            IF ( IAXGAS .GT. 0 ) VMRGRA(IATM,IAXGAS,IGRA) = 
     &        VMRATM(IATM,IAXGAS) * 10.0 / DNSGRA(IATM,IGRA)
          END DO
        END DO
        CALL MOVGRA ( NOMGRA, 0 )                ! Restore ref.prof to /ATMCOM/
      END IF
C
      DO IATM = 1, NATM
        IF ( PREATM(IATM) .GT. ARGMIN ) THEN
          LNPATM(IATM) = LOG ( PREATM(IATM) )
        ELSE
          LNPATM(IATM) = LOGMIN
        END IF
C
        DNSATM(IATM) = PREATM(IATM) * AVOG * 1.0E-4 / RGAS/TEMATM(IATM)
        IF ( DNSATM(IATM) .GT. ARGMIN ) THEN
          LNDATM(IATM) = LOG ( DNSATM(IATM) )
        ELSE
          LNDATM(IATM) = LOGMIN
        END IF
C
        RFRATM(IATM) = REFRAC ( IATM )
C
        IF ( RFRATM(IATM) .GT. ARGMIN ) THEN
          LNRATM(IATM) = LOG ( RFRATM(IATM) )
        ELSE
          LNRATM(IATM) = LOGMIN          
        END IF
C
        IF ( IATM .GT. 1 ) THEN
C 08FEB05: If layers are too narrow cannot calculate scale height accurately
C from difference, so use approximation RT/gM (R/gM=8.3/(9.8*29)=0.029km/K)
C Also prevent any negative density gradients
          IF ( LNDATM(IATM-1) - LNDATM(IATM) .LT. 0.01 ) THEN
            DSHATM(IATM-1) = 0.029 * TEMATM(IATM-1) 
C
C 23JUL03: Use refractivity difference to calculate scale height if possible 
C (allows for refractivity depending on H2O as well as density)
C NB similar calculation in VALGRA for 'DSH' option
          ELSE IF (  RFRATM(IATM) .GT. ARGMIN*1.0E6 ) THEN
             DSHATM(IATM-1) = ( HGTATM(IATM) - HGTATM(IATM-1) ) /
     &                        LOG ( RFRATM(IATM-1)/RFRATM(IATM) )
C If refractivity too small to calculate accurately, use density scale height
          ELSE IF ( LNDATM(IATM) .GT. LOGMIN ) THEN
            DSHATM(IATM-1) = ( HGTATM(IATM) - HGTATM(IATM-1) ) /
     &                       ( LNDATM(IATM-1) - LNDATM(IATM) )
          ELSE
            DSHATM(IATM-1) = 1.0E20
            WRITE (*,*) 'W-ATMAUX: Setting DSHATM large'
          END IF
        END IF
C
C Convert aerosol from extinction [km-1]/10^6 to ext/air.dens [cm-1*cm^3]
C Factor 10 comes from *10^6 (undo ppmv-ppv conversion in VMR) / 10^5 to 
C convert from ext/km to ext/cm. Only do this if IAXGAS > 0, which is the flag
C indicating that a new aerosol extinction profile has been set
C (NB this routine also called by OBSATM if extra layer inserted, so no extra
C conversion required)
        IF ( IAXGAS .GT. 0 ) 
     &    VMRATM(IATM,IAXGAS) = VMRATM(IATM,IAXGAS) * 10.0/DNSATM(IATM)
C
        DO IGAS = 1, NGAS
          IF ( VMRATM(IATM,IGAS) .GT. ARGMIN ) THEN
            LNVATM(IATM,IGAS) = LOG ( VMRATM(IATM,IGAS) )
          ELSE
            LNVATM(IATM,IGAS) = LOGMIN
          END IF
        END DO
C
      END DO
C
C Flag aerosol conversion as completed
      IAXGAS = 0
C
      END
