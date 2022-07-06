      SUBROUTINE CHKSRC ( MOL, IDX, WARN, FAIL, ERRMSG )
C
C VERSION
C     02-AUG-13  AD  Original.
C
C DESCRIPTION
C     Set/Get HITRAN/RFM Index for molec. which have two options.
C     On first call for molecule this sets the index, on subsequent calls
C     this index is returned
C
      IMPLICIT NONE
      SAVE
C
C ARGUMENTS
      CHARACTER*(*) MOL    !  I  Molecule name
      INTEGER       IDX    ! I/O Initial/Assigned value of molecule index
      LOGICAL       WARN   !  I  TRUE = warn if this is initial assignment
      LOGICAL       FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  RFMLOG ! Write text message to RFM log file.
C
C LOCAL CONSTANTS
      INTEGER       MAXMOL ! Max number of molecules listed
        PARAMETER ( MAXMOL = 100 ) ! only 10 or so ambiguous molecules required
C
C LOCAL VARIABLES
      INTEGER      IMOL           ! Counter for listed molecules
      INTEGER      IDXMOL(MAXMOL) ! Stored HITRAN/RFM indices of molecules
      INTEGER      NMOL           ! Number of listed molecules
      CHARACTER*7  MOLLST(MAXMOL) ! List of molecules stored so far
      CHARACTER*80 WRNMSG         ! Warning message if default is changed
C
C DATA STATEMENTS
      DATA NMOL / 0 /
      SAVE NMOL, IDXMOL, MOLLST
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      FAIL = .FALSE.
C
C See if source already assigned for this molecule
      DO IMOL = 1, NMOL 
        IF ( MOL .EQ. MOLLST(IMOL) ) THEN
          IDX = IDXMOL(IMOL)
          RETURN
        END IF
      END DO
C
C This shouldn't happen, but just in case
      IF ( NMOL .GT. MAXMOL ) THEN
        ERRMSG = 'F-CHKSRC: Local array size MAXMOL is too small'
        FAIL = .TRUE.
        RETURN
      END IF
C 
C New molecule so assign index
      NMOL = NMOL + 1
      MOLLST(NMOL) = MOL
      IDXMOL(NMOL) = IDX
      IF ( WARN ) THEN
        IF ( IDX .LT. 100 ) THEN
          WRNMSG = 'W-CHKSRC: changing to use line data for gas='//MOL
        ELSE
          WRNMSG = 'W-CHKSRC: changing to use x/s data for gas='//MOL
        END IF
        CALL RFMLOG ( WRNMSG, FAIL, ERRMSG )
      END IF
C
      END

