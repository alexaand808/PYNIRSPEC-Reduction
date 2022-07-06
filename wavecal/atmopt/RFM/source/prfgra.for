      SUBROUTINE PRFGRA ( PSI, IPRF, FAIL, ERRMSG )
C
C VERSION
C     31-DEC-99  AD  Original.
C
C DESCRIPTION
C     Load horizontal gradient profiles into /GRACOM/.
C     Called by ATMFIL for each profile.
C     Identifies if new location profile and copies data from current ATM
C     file into new location, profile determined by IPRF
C       IPRF=<-1 No profiles to load, just reserve space for PSI 
C       IPRF=-1  Pressure profile
C       IPRF=0   Temperature profile
C       IPRF>0   VMR profile for Absorber#IPRF
C
      IMPLICIT NONE
C
C ARGUMENTS
      REAL         PSI     !  I  Gradient angle
      INTEGER      IPRF    !  I  Profile identifier
      LOGICAL      FAIL    !  O  T=Fatal error detected
      CHARACTER*80 ERRMSG  !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  MOVGRA ! Exchange profiles between /ATMCOM/ and /GRACOM/ arrays.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'gracom.inc' ! Atmospheric 2-D field data
C
C LOCAL CONSTANTS
      REAL PSITOL          ! Tolerance for distinguishing psi angles [deg]
        PARAMETER ( PSITOL = 1.0E-6 ) ! 1e-6 deg = 900m on earth's surface 
C
C LOCAL VARIABLES
      INTEGER      IATM    ! Counter for atmospheric profile levels
      INTEGER      IGRA    ! Counter for no.different profile locations
      INTEGER      IGAS    ! Counter for absorbers
      INTEGER      JGRA    ! Secondary counter for diff.profile locations
      LOGICAL      FIRST   ! T=first call of this routine
      LOGICAL      NEWPSI  ! T=new profile location
      REAL         PSIPRV  ! Previous value of PSI
C
C DATA STATEMENTS
      DATA FIRST / .TRUE. /
      DATA IGRA / 0 /
      DATA PSIPRV / 0.0 /
      SAVE FIRST, PSIPRV, IGRA
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      IF ( FIRST ) THEN           ! First time this routine 
        NGRA = 0                  ! Initialise NGRA in /GRACOM/
        FIRST = .FALSE.
        NEWPSI = .TRUE.
      ELSE IF ( PSI .NE. PSIPRV ) THEN  ! Value of PSI has changed
        DO JGRA = 1, NGRA
          IF ( ABS(PSI-PSIGRA(JGRA)) .LT. PSITOL ) THEN ! Value used earlier
            IGRA = JGRA                                 ! save earlier index
            NEWPSI = .FALSE.                            ! not a new value
            GOTO 100
          END IF
        END DO
        NEWPSI = .TRUE.                 ! New value of PSI supplied
      ELSE
        NEWPSI = .FALSE.                ! PSI unchanged, so not new value
      END IF
 100  CONTINUE
      PSIPRV = PSI                      ! Copy this value for testing next time
C
C If a new value of PSI is supplied, check sufficient array space and insert
C into ordered list at position IGRA
      IF ( NEWPSI ) THEN
        IF ( NGRA .EQ. MAXGRA ) THEN                     ! No more array space
          FAIL = .TRUE.
          ERRMSG = 'F-PRFGRA: No. of horizontal gradient '//
     &             'profiles > MAXGRA in RFMSIZ.INC'
          RETURN
        ELSE
          NGRA = NGRA + 1                                ! Increment counter
        END IF
        IGRA = NGRA       ! assume new PSI is last (=largest) in list initially
        DO WHILE ( PSIGRA(MAX(1,IGRA-1)) .GT. PSI .AND. IGRA .GT. 1 )
          CALL MOVGRA ( IGRA-1, IGRA )                   ! Shuffe list 
          IF ( NOMGRA .EQ. IGRA-1 ) NOMGRA = IGRA
          IGRA = IGRA - 1 
        END DO
C
        PSIGRA(IGRA) = PSI
C
C If this location is near PSI=0.0, mark as the "nominal" (reference) profile
        IF ( ABS ( PSIGRA(IGRA) ) .LE. PSITOL ) NOMGRA = IGRA   
C
C Set flags for new profile location FALSE until actual profiles loaded
        FLTGRA(IGRA) = .FALSE.
        FLPGRA(IGRA) = .FALSE.
        DO IGAS = 1, NGAS
         FLVGRA(IGAS,IGRA) = .FALSE.
        END DO
      END IF
C
C From this point, IGRA points to profile location index, whether new or old.
      IF ( IPRF .EQ. -1 ) THEN                 ! Flag for pressure profile
        DO IATM = 1, NATM
          PREGRA(IATM,IGRA) = PREATM(IATM)
        END DO
        FLPGRA(IGRA) = .TRUE.
      ELSE IF ( IPRF .EQ. 0 ) THEN             ! Flag for temperature profile
        DO IATM = 1, NATM
          TEMGRA(IATM,IGRA) = TEMATM(IATM)
        END DO
        FLTGRA(IGRA) = .TRUE.
      ELSE IF ( IPRF .GE. 1 ) THEN             ! Flag for VMR profile gas#iprf
        DO IATM = 1, NATM
          VMRGRA(IATM,IPRF,IGRA) = VMRATM(IATM,IPRF)
        END DO
        FLVGRA(IPRF,IGRA) = .TRUE.
      END IF
C
      END
