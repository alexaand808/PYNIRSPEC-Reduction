      SUBROUTINE JACCHK ( TARGET, FINISH, ALT, FAIL, ERRMSG )
C
C VERSION
C     29-AUG-13  AD  Remove redundant local variables IGAS, IATM
C     25-SEP-08  AD  Bug#71: Correct check for exceeding MAXJAC
C     07-JAN-08  AD  Remove CFCGAS and use GASCHK to test for valid molecule
C     15-NOV-07  AD  Bug#66: failed to reset FAIL=TRUE after external calls
C     03-MAY-05  AD  Add jdxcon.inc
C     25-MAR-05  AD  Comment change: TARGET increased from C*8 max to C*10 max
C     01-JAN-04  AD  Allow isotopes
C     30-DEC-03  AD  Add ATMLEV to allow extra profile levels to be inserted
C     03-SEP-03  AD  Allow 'SFCTEM' and 'SFCEMS' values for TARGET
C                    Set LSTALT(2) = 0 for total column perturbations
C     05-AUG-99  AD  Slight modification of logic. Remove redundant MAXALT.
C                    Correction: change LSTALT(MAXJAC+1) to LSTALT(MAXJAC+2)
C     31-JUL-99  AD  Bug Fix: Add SAVE LSTALT. Remove unused JJAC.
C     04-FEB-99  AD  Modified definitions of perturbation levels
C     03-JAN-99  AD  Original.
C
C DESCRIPTION
C     Check details of Jacobians and set up related data.
C     Called by INPJAC for each species & altitude in *JAC sec. of drv.table.
C     (corresponding to once for each species to be retrieved)
C     Also called by JACFIL for each altitude contained in file of altitudes.
C     The TARGET field contains either
C     - the (new) species (eg 'CH4', 'TEM'); or
C     - ' ', in which case it is assumed that either
C        - FINISH = TRUE, in which case altitude list for current species has
C          been completed; or
C        - FINISH = FALSE, in which case new altitude is read from ALT.
C     If called with FINISH=TRUE before any altitudes have been loaded for the
C     current species it is assumed this is a column perturbation (for the 
C     complete atmosphere). Otherwise at least three altitudes are required to
C     define the perturbation levels required for the Jacobian calculations.
C
      IMPLICIT NONE
C
C ARGUMENTS
      CHARACTER*(*) TARGET ! I/O Target species for retrieval (max C*11)
      LOGICAL       FINISH !  I  T=No more altitudes for current species
      REAL          ALT    !  I  Retrieval altitude [km] 
      LOGICAL       FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  ATMLEV ! Find/insert atmospheric level for given altitude.
     &, CHKGAS ! Check for valid molecule name and get Hitran/RFM index
     &, JACVIB ! Check valid parameters for VT Jacobian
     &, JACISO ! Check valid parameters for isotopic Jacobian.
     &, LOCASE ! Convert text string to lower case.
C
C GLOBAL CONSTANTS
      INCLUDE 'jdxcon.inc' ! Indices for non-VMR Jacobians
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'jaccom.inc' ! Jacobian data
C
C LOCAL CONSTANTS
      REAL          ALTTOL ! Minimum altitude [km] distinguishable from zero
        PARAMETER ( ALTTOL = 0.001 ) ! 1m allows distinguishing filenames
C
C LOCAL VARIABLES
      INTEGER      IALT   ! Counter for altitude list
      INTEGER      IEND   ! Pointer to last character in TARGET
      INTEGER      IGQ    ! Global Quantum index of Vibrational level
      INTEGER      IJAC   ! Counter for jacobian rtvl. elements 
      INTEGER      IMOL   ! HITRAN/RFM index of molecule
      INTEGER      ISO    ! Isotope#1 (1=most abundant, etc, 0=all isotopes)
      INTEGER      JALT   ! Counter for altitude list
      INTEGER      JATM   ! Index of atmospheric profile level
      INTEGER      JGAS   ! Index of retrieval species in /GASCOM/ arrays
      INTEGER      LATM   ! Saved value of NATM
      INTEGER      LSTALT(MAXJAC+2) ! List of Atmos.level indices for ptb.alts.
      INTEGER      NALT   ! Number of stored altitudes for current species
      CHARACTER*7  GAS    ! Name of molecule identified for Jacobians (dummy)
      CHARACTER*80 SPECES ! Current retrieval species
C
C Data statements are only to satisfy various code-checking programs: none
C of these values is actually used in the correct sequence of calls
      DATA SPECES / ' ' /
      DATA IEND   / 0 / 
      DATA JGAS   / 0 / 
      DATA NALT   / 0 / 
      SAVE SPECES, IEND, JGAS, NALT, LSTALT
C
C EXECUTABLE CODE ------------------------------------------------------------
C
      FAIL = .FALSE.
C
C If called with TARGET .NE. ' ' this is the first call for a new retrieval
C species with the altitude list to follow on subsequent calls
C
      IF ( TARGET .NE. ' ' ) THEN      
        CALL LOCASE ( TARGET, TARGET )      ! Convert to lower case
        IF ( TARGET .EQ. 'tem' ) THEN
          JGAS = MAXGAS + JDXTEM
        ELSE IF ( TARGET .EQ. 'pre' ) THEN
          JGAS = MAXGAS + JDXPRE
        ELSE IF ( TARGET .EQ. 'sfctem' ) THEN
          IF ( .NOT. SFCFLG ) THEN
            FAIL = .TRUE.
            ERRMSG = 'F-JACCHK: SFCTEM Jacobians require '//
     &               'surface flag SFC to be enabled'
            RETURN
          END IF
          JGAS = MAXGAS + JDXSFT
        ELSE IF ( TARGET .EQ. 'sfcems' ) THEN
          IF ( .NOT. SFCFLG ) THEN
            FAIL = .TRUE.
            ERRMSG = 'F-JACCHK: SFCEMS Jacobians require '//
     &               'surface flag SFC to be enabled'
            RETURN
          END IF
          JGAS = MAXGAS + JDXSFE
        ELSE 
          CALL CHKGAS ( TARGET, GAS, IMOL, ISO, IGQ, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
          IF ( IMOL .EQ. 0 ) THEN
            FAIL = .TRUE.
            ERRMSG = 'F-JACCHK: unrecognised target parameter for'//
     &               ' Jacobian calc,='//TARGET
            RETURN
          ELSE IF ( IGSMOL(IMOL) .EQ. 0 ) THEN
            FAIL = .TRUE.
            ERRMSG = 'F-JACCHK: No profile loaded for Jacobian calc.'
     &                 //' for species='//TARGET
            RETURN
          ELSE                                  
            IF ( IGQ .GT. 0 ) THEN         ! Vib Tem Jacobian
              CALL JACVIB (  IMOL, ISO, IGQ, JGAS, FAIL, ERRMSG )
              IF ( FAIL ) RETURN
            ELSE IF ( ISO .GT. 0 ) THEN    ! Isotopic Jacobian
              CALL JACISO ( IMOL, ISO, JGAS, FAIL, ERRMSG )
              IF ( FAIL ) RETURN
            ELSE
              JGAS = IGSMOL(IMOL)
            END IF
          END IF
        END IF
        DO IJAC = 1, NJAC
          IF ( IGSJAC(IJAC) .EQ. JGAS ) THEN
            FAIL = .TRUE.
            ERRMSG = 'F-JACCHK: Repeated retrieval species '//
     &               'for Jacobian calc.='//TARGET
            RETURN
          END IF
        END DO
        IEND = INDEX ( TARGET, ' ' ) - 1
        SPECES = TARGET            
        NALT = 0                  ! Normal exit for first call with new species
C
C If called with FINISH=TRUE, transfer list of sorted alts for this species
C to JACCOM.INC
C Note: NJAC is number of previously stored Jacobians, NALT-2 is number to
C       be added to this list.
      ELSE IF ( FINISH ) THEN
        IF ( NALT .EQ. 0 ) THEN    ! Indicates column or surface perturbation
          NALT = 3
          IF ( NJAC + NALT - 2 .GT. MAXJAC ) THEN
            FAIL = .TRUE.
            WRITE ( ERRMSG, '(A,I11)' ) 'F-JACCHK: No.Jacobian '//
     &        'elements > dimension MAXJAC in RFMSIZ.INC,=', MAXJAC
            RETURN
          END IF
          IF ( SPECES(1:IEND) .EQ. 'sfctem' .OR. 
     &         SPECES(1:IEND) .EQ. 'sfcems'      ) THEN  ! surface perturbation
            LSTALT(1) = 0
            LSTALT(2) = 0
            LSTALT(3) = 0
          ELSE                                           ! column perturbation
            LSTALT(1) = 0
            LSTALT(2) = 0
            LSTALT(3) = NATM + 1 
          END IF
        ELSE IF ( NALT .LE. 2 ) THEN
          FAIL = .TRUE.
          ERRMSG = 'F-JACCHK: at least 3 altitudes required '//
     &             'for Jacobians, species='//SPECES(1:IEND)
          RETURN
        END IF
        DO IALT = 2, NALT-1
          NJAC = NJAC + 1
          ILOJAC(NJAC) = LSTALT(IALT-1)
          IATJAC(NJAC) = LSTALT(IALT)
          IUPJAC(NJAC) = LSTALT(IALT+1)
          IGSJAC(NJAC) = JGAS
        END DO          
C
C Anything else will be adding to the list of Jacobian altitudes 
C to be calculated
      ELSE
        IF ( NJAC + NALT - 2 .GE. MAXJAC ) THEN  !  NALT not yet incremented
          FAIL = .TRUE.
          WRITE ( ERRMSG, '(A,I11)' ) 'F-JACCHK: No.Jacobian elements'//
     &    ' > dimension MAXJAC in RFMSIZ.INC,=', MAXJAC
          RETURN
        END IF
        IF ( ALT .LE. HGTATM(1) - ALTTOL ) THEN
          JATM = 0
        ELSE IF ( ALT .GE. HGTATM(NATM) + ALTTOL ) THEN
          JATM = NATM + 1
        ELSE 
          LATM = NATM               ! Save current value of NATM
          CALL ATMLEV ( ALT, JATM, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
C If an extra profile level was inserted, ATMLEV should adjust previously
C stored Jacobians but not the currently stored LSTALT values
          IF ( NATM .GT. LATM ) THEN  ! T=Extra profile level was inserted
            DO IALT = 1, NALT
              IF ( LSTALT(IALT) .GE. JATM ) 
     &          LSTALT(IALT) = LSTALT(IALT) + 1
            END DO
          END IF
        END IF
        NALT = NALT + 1
C
C Insert altitude into ordered list, checking for repeated values (just for 
C the current species)
        DO IALT = 1, NALT-1
          IF ( LSTALT(IALT) .EQ. JATM ) THEN
            FAIL = .TRUE.
            WRITE ( ERRMSG, '(A,G12.5,A)' ) 'F-JACCHK: Repeated '//
     &        'altitude=', ALT, ' [km] for species='//SPECES(1:IEND)
            RETURN
          ELSE IF ( LSTALT(IALT) .GT. JATM ) THEN  ! insert into list here
            DO JALT = NALT-1, IALT, -1
              LSTALT(JALT+1) = LSTALT(JALT)
            END DO
            LSTALT(IALT) = JATM
            GOTO 300
          END IF
        END DO
        LSTALT(NALT) = JATM                        ! add to end of list
  300   CONTINUE
      END IF
C 
      FAIL = .FALSE.           ! Normal exit for all cases
      END
