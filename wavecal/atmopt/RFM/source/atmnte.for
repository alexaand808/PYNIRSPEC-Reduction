      SUBROUTINE ATMNTE ( IDXVIB, VIBATM, FAIL, ERRMSG )
C
C VERSION
C     13-AUG-13  AD  Remove redundant variable IPRF
C     09-JAN-07  AD  Original. Based on MORSE routine VIBNTE
C
C DESCRIPTION
C     Load vibrational temperature profile from .atm file
C     Called by FILPRF if NTE flag enabled and VT profile found in .atm file.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      IDXVIB  !  I  Encoded details of Vibrational level
      REAL         VIBATM  !  I  Vib.Tem profile interpolated to *ATM levels
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  RFMLOG ! Write text message to RFM log file.

C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! Dimensions for common arrays
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'ntecom.inc' ! Non-LTE data
C
C LOCAL CONSTANTS
      INTEGER      NVIB    ! No. of Vib.Temps with locally listed Energies
        PARAMETER  ( NVIB = 2 ) 
C
C LOCAL VARIABLES
      INTEGER      IATM    ! Counter for atmospheric profile levels
      INTEGER      IMOL    ! HITRAN ID for VT molecule 
      INTEGER      INTE    ! Counter for stored NTE datasets
      INTEGER      IGAS    ! Index in *GAS arrays of molecule
      INTEGER      IGQ     ! Global Quantum Index of vibrational level
      INTEGER      ISO     ! Isotope# (or 0=apply to all isotopes)
      INTEGER      IVIB    ! Counter for locally stored Energy levels info
      INTEGER      LSTVIB(NVIB) ! List of VT indices for Energy levels
      REAL         ENGVIB(NVIB) ! Energy for each Vib.level
      CHARACTER*80 MESSGE  ! Text message for log file
C
C DATA STATEMENTS
C Copy data from .nte files if other levels are required
      DATA LSTVIB / 5000002,   5001002 /
      DATA ENGVIB / 2143.2711, 2143.2711 /
      SAVE LSTVIB, ENGVIB
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      IMOL = IDXVIB/1000000     ! Extract Molecule ID, Isotope#, GQ#
      ISO  = MOD ( IDXVIB, 1000000 ) / 1000
      IGQ  = MOD ( IDXVIB, 1000 )
C
      IGAS = IGSMOL(IMOL)

C Check if a VT profile has already been loaded for this mol/iso/igq
      DO INTE = 1, NNTE
        IF ( IDXVIB .EQ. IDXNTE(INTE) ) THEN
          DO IATM = 1, NATM
            TEMNTE(IATM,INTE) = VIBATM(IATM)
          END DO
          MESSGE = 'W-ATMNTE: VT profile superseded'
          CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
          RETURN  ! Normal exit with VT profile replaced
        END IF
      END DO
C
C If not already loaded this VT, add to list
      IF ( NNTE .GT. MAXNTE ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-ATMNTE: MAXNTE in rfmsiz.inc too small'
        RETURN
      END IF
C
      NNTE = NNTE + 1
      MNTE = NNTE
      IDGNTE(NNTE) = IMOL
      ISONTE(NNTE) = ISO
      IGQNTE(NNTE) = IGQ
      IDXNTE(NNTE) = IDXVIB
      QPRNTE(NNTE) = .FALSE.
C
C Check that Gas/IGQ VTs are either marked for all isotopes or individual iso.
      DO INTE = 1, NNTE-1
        IF ( IDGNTE(INTE) .EQ. IDGNTE(NNTE) .AND.
     &       IGQNTE(INTE) .EQ. IGQNTE(NNTE) ) THEN
          IF ( ISONTE(INTE) .EQ. 0 .OR. ISONTE(NNTE) .EQ. 0 ) THEN
            FAIL = .TRUE.
            WRITE ( ERRMSG, '(A,A,A,I3,A)' ) 'F-ATMNTE: Gas ', 
     &        CODGAS(IGAS), ' GQIdx#', IGQNTE(NNTE),
     &        ' has VT flagged for both * and specific isotps.'
            RETURN
          ENDIF
        ENDIF 
      ENDDO
C
C Check that the energy level for this VT is stored here 
      DO IVIB = 1, NVIB
        IF ( IDXVIB .EQ. LSTVIB(IVIB) ) THEN
          ENGNTE(NNTE) = ENGVIB(IVIB)
          FAIL = .FALSE.
          RETURN                       ! Normal exit for new VT profile
        ENDIF
      END DO
C
      FAIL = .TRUE.
      WRITE ( ERRMSG, '(A,A,A,I1,A,I3)' ) 
     &  'F-ATMNTE: No Energy level listed for Gas=', CODGAS(IGAS), 
     &  ' Iso=', ISONTE(INTE), ' GQ#=', IGQNTE(INTE)
      RETURN
C
      END
