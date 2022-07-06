      SUBROUTINE ATMISO ( IMOL, ISO, FAIL, ERRMSG )
C
C VERSION
C     10JAN14 AD Remove TCOGAS, WIDGAS
C     07AUG13 AD Remove redundant variable JISO
C     07JAN08 AD Simplify assuming IMOL and ISO supplied as input arguments
C     22MAR04 AD Set ISOMOL = TRUE
C     13MAR04 AD Set CTMGAS = FALSE for isotopes
C     10FEB04 AD Check ISOGAS for IGAS. Set IDIGAS
C     30APR02 AD Correction: set CHGGAS = TRUE and set WGTGAS(1,NGAS)
C     16FEB02 AD Original.
C
C DESCRIPTION
C     Create separate isotope profiles if specified in .atm file
C     Called by ATMFIL and JACISO for each species where isotope is specified.
C
      IMPLICIT NONE
C
C ARGUMENTS 
      INTEGER       IMOL   !  I  HITRAN/RFM Molecule ID
      INTEGER       ISO    !  I  Isotope# of molecule
      LOGICAL       FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG !  O  Error message written if FAIL is TRUE 
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'gascom.inc' ! Molecular and isotope data
C
C LOCAL VARIABLES
      INTEGER       IGAS   ! Index of species in *GAS arrays
      INTEGER       IISO   ! Isotope counter
      INTEGER       JGAS   ! Counter over gases
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Find index of default isotopic concentration for species (if required)
C This routine should only be called if molecule already assigned index IGAS
      IGAS = IGSMOL(IMOL)
      IF ( IGAS .EQ. 0 ) STOP 'F-ATMISO: Logical Error#1'
      IF ( IDXGAS(IGAS) .NE. IMOL ) STOP 'F-ATMISO: Logical Error#2'
C
C Check isotope number is valid for this gas
      IF ( ISO .GT. NISGAS(IGAS) ) THEN
        WRITE ( ERRMSG, '(A,I3,A,I3,A,A)' ) 
     &    'F-ATMISO: Isotope#', ISO, ' > Max#', NISGAS(IGAS),
     &    ' for species ', CODGAS(IGAS)
        RETURN
      END IF
C
C ISOGAS initially assigned to default IGAS for all isotopes, so if still 
C IGAS then it needs to be reassigned to new "species" 
      IF ( ISOGAS(ISO,IGAS) .EQ. IGAS ) THEN
        IF ( NGAS .GE. MAXGAS ) THEN
          WRITE ( ERRMSG, '(A,I11)' ) 
     &   'F-ATMISO: MAXGAS in RFMSIZ.INC too small, value=', MAXGAS
          RETURN
        END IF
        NGAS = NGAS + 1
        IDIGAS(NGAS) = ISO
        IDXGAS(NGAS) = IDXGAS(IGAS)
        ISOGAS(ISO,IGAS) = NGAS
        ISOMOL(IDXGAS(IGAS)) = .TRUE.  
        NISGAS(NGAS) = NISGAS(IGAS)
C Copy isotope location information from default gas index IGAS to new NGAS
        DO IISO = 0, NISGAS(IGAS)
          ISOGAS(IISO,NGAS) = ISOGAS(IISO,IGAS)
          IF ( ISOGAS(IISO,IGAS) .NE. IGAS ) THEN ! Update other isotopes also
            JGAS = ISOGAS(IISO,IGAS)
            ISOGAS(ISO,JGAS) = NGAS
          END IF
        END DO
        SHPGAS(NGAS) = SHPGAS(IGAS)
        CTMGAS(NGAS) = .FALSE.
        WGTGAS(ISO,NGAS) = WGTGAS(ISO,IGAS)
        WGTGAS(1,NGAS) = WGTGAS(1,IGAS)    ! Main isotope wgt used in AWIDTH
        CODGAS(NGAS) = CODGAS(IGAS)
        CHGGAS(NGAS) = .TRUE.
        FAIL = .FALSE.
        RETURN
      END IF
      FAIL = .FALSE.
C
      END 
