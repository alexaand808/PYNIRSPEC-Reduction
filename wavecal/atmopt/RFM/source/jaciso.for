      SUBROUTINE JACISO ( IMOL, ISO, JGAS, FAIL, ERRMSG )
C
C VERSION
C     07-JAN-08  AD  Rewritten and simplified. Remove JACIGQ.
C     03-MAY-05  AD  Add JACIGQ to allow for VT Jacobians
C     01-JAN-04  AD  Original.
C
C DESCRIPTION
C     Check valid parameters for isotopic Jacobian.
C     Called by JACCHK for each target where isotope is specified.
C
      IMPLICIT NONE
C
C ARGUMENTS 
      INTEGER       IMOL   !  I  HITRAN/RFM ID for molecule
      INTEGER       ISO    !  I  Isotope# (1=most abundant, etc)
      INTEGER       JGAS   !  O  Index of isotope,gas in *GAS arrays
      LOGICAL       FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG !  O  Error message written if FAIL is TRUE 
C
      EXTERNAL
     &  ATMISO ! Create separate isotope profiles if specified in .atm file
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'gascom.inc' ! Molecular and isotope data
C
C LOCAL VARIABLES
      INTEGER       IATM   ! Counter for atmospheric profile levels
      INTEGER       IGAS   ! Main index of species in *GAS arrays
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      FAIL = .FALSE.
      IGAS = IGSMOL(IMOL)     ! default profile index for all isotopomers
      IF ( IGAS .EQ. 0 ) STOP 'F-JACISO: Logical Error#1'
      IF ( IDXGAS(IGAS) .NE. IMOL ) STOP 'F-JACISO: Logical Error#2'
C
C Check Isotope# is in valid range for target species
      IF ( ISO .LT. 1 ) THEN
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,I11)' ) 
     &    'F-JACISO: Invalid Isotope# Jacobian, value=', ISO
        RETURN
      ELSE IF ( ISO .GT. NISGAS(IGAS) ) THEN
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,I3,A,I3,A,A)' ) 
     &    'F-JACISO: Isotope#', ISO, ' > Max#', NISGAS(IGAS),
     &    ' for species=', CODGAS(IGAS)
        RETURN
      END IF
C
C Check if isotopic profile is loaded, otherwise increment list of gases as in
C ATMISO and copy "standard" profile concentration
      IF ( ISOGAS(ISO,IGAS) .NE. IGAS ) THEN       ! already loaded
        JGAS = ISOGAS(ISO,IGAS) 
      ELSE
        CALL ATMISO ( IMOL, ISO, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
C Adding new isotopic concentration will increment list of gases (last is NGAS)
        JGAS = NGAS
        DO IATM = 1, NATM
          VMRATM(IATM,JGAS) = VMRATM(IATM,IGAS)
          LNVATM(IATM,JGAS) = LNVATM(IATM,IGAS)
        END DO
      END IF
C
      END 
