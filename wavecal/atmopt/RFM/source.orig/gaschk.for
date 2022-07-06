      SUBROUTINE GASCHK ( GASSTR, NOTGAS, FAIL, ERRMSG )
C
C VERSION
C     14-AUG-13  AD  Avoid modifying GASSTR
C     07-AUG-13  AD  Remove redundant local variable IGAS
C     04-JAN-08  AD  Add CHKGAS to identify molecules
C                    Set NISGAS, GASISO locally for all molecules
C     10-FEB-03  AD  Set IDIGAS(NGAS) = 0
C     23-FEB-03  AD  Add BrO as ID=40, move GeH4 (was ID=40) to new ID=46
C     06-MAR-02  AD  Add SF6 as a xsc molecule (#64) - replaces X64
C                    Change line molecule to "sf6q"
C                    Change "clonoq" back to "clono2q" 
C     16-FEB-02  AD  Add HDO#39
C     07-FEB-01  AD  Add CFCs/HFCs #70-81
C     16-OCT-00  AD  Add GeH4(#40), C3H8(41), C2N2(42), C4H2(43), HC3N(44) and
C                    C3H4(45)
C     13-SEP-99  AD  Change "clono2q" to "clonoq" and add "c2h4" as Gas#38.
C     20-APR-99  AD  Replace CFC chemical names by F-numbers
C     21-DEC-98  AD  Add user-defined cross-sections X64 ... X69
C                    Add test for duplication of ClONO2
C     01-DEC-97  AD  Add GASCTM
C                    Allow use of HITRAN index instead of gas name
C     02-OCT-97  AD  Make repeated gas a fatal error rather than a warning
C     22-JUL-97  AD  Add check for qualifier data, new module GASQAL
C     27-MAR-97  AD  Use CFCGAS to convert any CFC codes.
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Check Gas name and set indices.
C     Called by INPGAS for each field in *GAS section.
C     Note that all strings are converted to lower case for internal use.
C
      IMPLICIT NONE
C
C ARGUMENTS
      CHARACTER*(*) GASSTR  !  I  Gas name to be tested (assumed non-blank)
      LOGICAL       NOTGAS  !  O  Set TRUE if gas not recognised
      LOGICAL       FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  CHKGAS ! Check if GAS is valid molecule and get Hitran/RFM index
     &, GASCTM ! Initialise for Continuum data
     &, GASISO ! Load isotopic data for required gases.
     &, GASQAL ! Check & load isotope/band qualifiers for each molecule.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'shpcon.inc' ! Line-shape codes
C
C COMMON VARIABLES
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'qalcom.inc' ! List of band/isotope qualifiers for line data.
C
C LOCAL VARIABLES
      INTEGER      IDUMMY  ! Dummy argument for CHKGAS
      INTEGER      IMOL    ! HITRAN index of molecule
      INTEGER      IPT     ! Pointer to start of qual.info in GASSTR
      INTEGER      ISO     ! Isotope#
      INTEGER      LPT     ! Length of non-blank part of GASSTR
      INTEGER      LQS     ! Length of QALSTR written
      LOGICAL      ANYQAL  ! Set TRUE if any qualifiers added
      CHARACTER*7  GAS     ! Gas name read from GASSTR
      CHARACTER*80 QALSTR  ! Qualifier data associated with molecule
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      FAIL   = .FALSE.
      NOTGAS = .TRUE.
C
C Set LPT to last non-blank character in GASSTR
      LPT = INDEX ( GASSTR, ' ' ) - 1
      IF ( LPT .EQ. -1 ) LPT = LEN ( GASSTR )
C
C Check for any qualifier data associated with molecule -look for '(' in string
      LQS = 1
      ANYQAL = INDEX ( GASSTR, '(' ) .GT. 1   ! Cannot start with '('
      IF ( ANYQAL ) THEN
        IPT = INDEX ( GASSTR, '(' )
        LQS = LPT + 1 - IPT
        QALSTR(1:LQS) = GASSTR(IPT:LPT)
       END IF
C
C CHKGAS will ignore any brackets
      CALL CHKGAS ( GASSTR, GAS, IMOL, ISO, IDUMMY, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      NOTGAS = ( IMOL .EQ. 0 ) 
      IF ( NOTGAS ) RETURN     ! Exit with no molecule identified
C
C Normally CHKGAS will also remove any '(...)' and return Isotope#, but in the
C *GAS section of the driver table the qualifiers can be more complicated so
C these are removed earlier. However if gas is a recognised isotopomeric 
C abbreviation ('hdo' or 'ch3d') CHKGAS could still return isotope numbers, 
C in which case pretend these were part of the original qualifier string
      IF ( ISO .GT. 0 ) THEN   ! Identified isotopomeric abbreviation
        ANYQAL = .TRUE.
        IF ( LQS .EQ. 1 ) THEN
          WRITE ( QALSTR, '(A,I1,A)' ) '(', ISO, ')'
          LQS = 3
        ELSE                   ! Qualifiers already attached
          WRITE ( QALSTR, '(A,I1,A,A)' ) '(', ISO, ')', QALSTR(1:LQS)
          LQS = LQS + 3
        END IF
      END IF
C
C Check molecule isn't already loaded
      IF ( IGSMOL(IMOL) .NE. 0 ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-GASCHK: Repeated Gas='//GAS
        RETURN
      END IF
C Check enough space to add new molecule
      IF ( NGAS .EQ. MAXGAS ) THEN  ! Not enough space for another gas
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,I3)' )
     &  'F-GASCHK: No.gases reqd > MAXGAS in RFMSIZ.INC,=', MAXGAS
        RETURN
      END IF
C Everything OK so add to list in gascom.inc
      NGAS = NGAS + 1
      IDIGAS(NGAS) = 0              ! default isotope profile
      IDXGAS(NGAS) = IMOL
      IGSMOL(IMOL) = NGAS
      ISOGAS(0,NGAS) = NGAS
      NISGAS(NGAS) = 0              ! could be changed by GASISO
      ISOMOL(IMOL) = .FALSE.
      NTEGAS(NGAS) = .FALSE.  
      CODGAS(NGAS) = GAS
      IF ( IMOL .LE. MAXHLN ) THEN           ! HITRAN Line molecule
        CALL GASISO ( FAIL, ERRMSG )         ! Load isotope,line params
        IF ( FAIL ) RETURN
        IF ( SHPFLG ) THEN                   ! Line shapes to be supplied
          SHPGAS(NGAS) = 0
        ELSE                                 ! Set default line shape
          SHPGAS(NGAS) = SHPVOI          
        END IF
        CALL GASCTM ( CTMFLG, ANYQAL, LQS,         ! May reset ANYQAL=F
     &                QALSTR(1:LQS), FAIL, ERRMSG ) 
        IF ( FAIL ) RETURN
        IF ( ANYQAL ) THEN
          USEQAL(NGAS) = .TRUE.
          CALL GASQAL ( QALSTR(1:LQS), FAIL, ERRMSG ) ! Check qual. info
          IF ( FAIL ) RETURN
        ELSE
          USEQAL(NGAS) = .FALSE.
        END IF
      ELSE                                    ! Cross-section molecule
        SHPGAS(NGAS) = SHPXSC         
        IF ( ANYQAL ) THEN
          FAIL = .TRUE.
          ERRMSG = 'F-GASCHK: qualifier data not applicable '//
     &      'to X/S molecules: '//GAS
          RETURN
        END IF
      END IF
      FAIL = .FALSE.                      ! Normal exit with new molecule added 
C
      END
