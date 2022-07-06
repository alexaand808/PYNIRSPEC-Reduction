      SUBROUTINE NTETEM ( IDG, ISO, IGQ, IATMLO, 
     &                    NLEV, LNPLEV, TEMLEV, ENG, FAIL, ERRMSG )
C
C VERSION
C     09-JAN-08  AD  Set IDXVIB and allow for some VT preloaded in ATM section
C     20-FEB-01  AD  Duplicte VTs to top of atmospheric profile if necessary
C     09-AUG-00  AD  Store VTs as differences from KTs rather than absolute
C     01-JUL-98  AD  Remove IATMHI argument. Assume NTE profile always extends
C                    to top of ATM profile
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     21-SEP-96  AD  Original.
C
C DESCRIPTION
C     Store Vib.Temp profiles read from NTE file.
C     Called by NTEFIL for each Vib.Temp profile in NTE file.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER       IDG        !  I  HITRAN ID for molecule
      INTEGER       ISO        !  I  Isotope#
      INTEGER       IGQ        !  I  Global Quantum#
      INTEGER       IATMLO     !  I  Lower index of ATM profile for interpolatn.
      INTEGER       NLEV       !  I  No.of profile levels
      REAL          LNPLEV(*)  !  I  ln(p/mb) profile
      REAL          TEMLEV(*)  !  I  Vib.Temp difference profile [K]
      REAL          ENG        !  I  Energy [cm-1 of Vib.level
      LOGICAL       FAIL       !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG     !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  LOOKUP ! General purpose interpolation routine
     &, RFMLOG ! Write text message to RFM log file.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'ntecom.inc' ! Non-LTE data 
      INCLUDE 'qfncom.inc' ! non-LTE Partition Functions
C
C LOCAL VARIABLES
      INTEGER      IATM    ! Counter for atmospheric profile levels
      INTEGER      IGAS    ! Index of GAS in RFM list /GASCOM/
      INTEGER      INTE    ! Counter for NTE Vib.Temp profiles
      INTEGER      IQFN    ! Counter for Partition Function profiles
      INTEGER      J       ! Saved index for interpolation LOOKUP
      CHARACTER*80 MESSGE  ! Warning message sent to log file
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      FAIL   = .FALSE.
C
C Check if molecule actually required
C
      IF ( IGSMOL(IDG) .EQ. 0 ) RETURN
      IGAS = IGSMOL(IDG)
C
C Check not already loaded. Files INTE=1:MNTE have been preloaded in the *ATM
C section - assume that these have precedence so just issue warning if 
C duplicated. However there should be no duplicates within the *NTE section
C
      DO INTE = 1, NNTE
        IF ( IDGNTE(INTE) .EQ. IDG .AND. ISONTE(INTE) .EQ. ISO .AND.
     &       IGQNTE(INTE) .EQ. IGQ ) THEN
          IF ( INTE .LE. MNTE ) THEN      ! don't overwrite VT from ATM section
            WRITE ( MESSGE, '(A,A7,A,I1,A,I2,A)' ) 
     &        'W-NTETEM: ignore VT for gas=', CODGAS(IGAS), 
     &        ', iso=', ISO, ' GQIdx=', IGQ, 
     &        ' - preloaded in *ATM section'
            CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
          ELSE             
            FAIL = .TRUE.
            WRITE ( ERRMSG, '(A,A,A,I2,A,I11)' )
     &       'F-NTETEM: Repeated Vib.Tem profile for gas=', 
     &       CODGAS(IGAS), ', ISO#=', ISO, ', GQIdx=', IGQ
          END IF
          RETURN
        END IF
      END DO
C
C Check enough space in list to add another
C
      IF ( NNTE .EQ. MAXNTE ) THEN
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,I11)' )
     &    'F-NTETEM: No.Different Vib.Temp profiles > '//
     &    'MAXNTE in RFMSIZ.INC,=', MAXNTE
        RETURN
      END IF
C
C Add to list
C
      NNTE = NNTE + 1
      IDGNTE(NNTE) = IDG
      ISONTE(NNTE) = ISO
      IGQNTE(NNTE) = IGQ
      IDXNTE(NNTE) = IDG * 1000000 + ISO * 1000 + IGQ
      ENGNTE(NNTE) = ENG
C
C Interpolate to atmospheric profile
C Assume VT=KT below specified VT profile
      DO IATM = 1, IATMLO-1
        TEMNTE(IATM,NNTE) = 0.0
      END DO
      J = 1
      DO IATM = IATMLO, NATM
        IF ( LNPATM(IATM) .LT. LNPLEV(NLEV) ) THEN ! Duplicate above VT prof.
          TEMNTE(IATM,NNTE) = TEMLEV(NLEV)
        ELSE
          CALL LOOKUP ( LNPATM(IATM), TEMNTE(IATM,NNTE), LNPLEV,
     &                  TEMLEV, NLEV, J, 1 )     ! TEMLEV is VT-KT already
        END IF
      END DO
C
C Find if there is a Vibrational Partition Function for this profile
C
      QPRNTE(NNTE) = .FALSE.
      DO IQFN = 1, NQFN
        IF ( IDGQFN(IQFN) .EQ. IDG .AND. ISOQFN(IQFN) .EQ. ISO ) THEN
          QPRNTE(NNTE) = .TRUE.
          DO IATM = 1, NATM
            QFNNTE(IATM,NNTE) = PRFQFN(IATM,IQFN)
          END DO
          RETURN            
        END IF
      END DO
C
      END
