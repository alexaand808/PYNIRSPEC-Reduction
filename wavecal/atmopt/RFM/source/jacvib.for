      SUBROUTINE JACVIB ( IMOL, ISO, IGQ, JGAS, FAIL, ERRMSG )
C
C VERSION
C     07-AUG-13  AD  Remove redundant local variable IOS
C     07-JAN-08  AD  Rewritten to evaluate JGAS
C     03-MAY-05  AD  Original.
C
C DESCRIPTION
C     Check valid parameters for VT Jacobian
C     Called by JACCHK if Global Quantum number assigned to molecule name in
C     *JAC section of the driver table.
C
      IMPLICIT NONE
C
C ARGUMENTS 
      INTEGER       IMOL   !  I  HITRAN/RFM ID for molecule
      INTEGER       ISO    !  I  Isotope# (1=most abundant, etc), 0=all iso.s
      INTEGER       IGQ    !  I  Global Quantum Index of vibrational level
      INTEGER       JGAS   !  O  Jacobian index encoded from input arguments
      LOGICAL       FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG !  O  Error message written if FAIL is TRUE 
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C GLOBAL VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
C
C LOCAL VARIABLES
      INTEGER       IGAS   ! Index of absorber in GASCOM arrays
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      IGAS = IGSMOL(IMOL)  ! all isotopic variants use same default index 
C
      IF ( .NOT. NTEFLG ) THEN 
        FAIL = .TRUE.
        ERRMSG = 'F-JACVIB: VT Jacobians require NTE flag to be set'
        RETURN
      END IF
C
C Read global quantum number
      IF ( IGQ .LE. 0 .OR. IGQ .GE. 100 ) THEN   ! 99 is arbitrary maximum
        WRITE ( ERRMSG, '(A,I11)' ) 
     &    'F-JACVIB: Invalid GQ Index# in IGQSTR: ', IGQ
        RETURN
      END IF
C
C Check Isotope# is in valid range for target species (0 = use all isotopes)
      IF ( ISO .LT. 0 ) THEN
        WRITE ( ERRMSG, '(A,I11)' ) 
     &    'F-JACVIB: Invalid Isotope# in Vib.Tem Jacobian, value=', ISO
      ELSE IF ( ISO .GT. NISGAS(IGAS) ) THEN
        WRITE ( ERRMSG, '(A,I3,A,I3,A,A)' ) 
     &    'F-JACVIB: Isotope#', ISO, ' > Max#', NISGAS(IGAS),
     &    ' for species=', CODGAS(IGAS)
        RETURN
      END IF
C
C This coding needs to match decoding in subroutine CHKVTJ
      JGAS = IMOL * 1000000 + ISO * 1000 + IGQ
C
      FAIL = .FALSE.
C
      END 
