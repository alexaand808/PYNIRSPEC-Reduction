      SUBROUTINE LEVUPD ( LEVINI, ITAN, IATM, IDIR, KTAN1, KTAN2 )
C
C VERSION
C     04-APR-06  AD  Correction: save IJAC1, IJAC2
C     01-JAN-04  AD  Original.
C
C DESCRIPTION
C     Set up calculations for intermediate output levels
C     Called by RADTRA for each nominal tan.path and atmospheric level
C     For a given nominal tangent path ITAN, the set of tangent paths KTAN1-2
C     assigned to intermediate outputs is gradually increased as the path 
C     integration proceeds down then up through the atmospheric layers.
C
      IMPLICIT NONE 
C
C ARGUMENTS
      LOGICAL  LEVINI      ! I/O T=initialise for new path, reset=F on exit
      INTEGER  ITAN        !  I  Index of nominal tangent path
      INTEGER  IATM        !  I  Index of atmospheric layer
      INTEGER  IDIR        !  I  Direction of ray path integ.(-1=down,+1=up)
      INTEGER  KTAN1       !  O Index of 1st intermediate o/p assoc. with ITAN
      INTEGER  KTAN2       !  O Index of latest intermed. o/p assoc. with ITAN
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'jaccom.inc' ! Jacobian data
C
C LOCAL VARIABLES
      INTEGER      IJAC    ! Last JAC index initialised for current ITAN
      INTEGER      IJAC1   ! 1st JAC index for current ITAN
      INTEGER      IJAC2   ! last JAC index for current ITAN
      INTEGER      KTAN    ! TAN index for new o/p path
      INTEGER      KTAN1L  ! Locally saved KTAN1
      INTEGER      KTAN2L  ! Locally saved KTAN2
      INTEGER      LASTAN  ! ITAN saved from previous call
C
C DATA STATEMENTS
      DATA IJAC   / 0 /    ! Value not used
      DATA IJAC1   / 0 /   ! Value not used
      DATA IJAC2   / 0 /   ! Value not used
      DATA LASTAN / 0 /    ! ITAN=0 should never occur 
      DATA KTAN1L / 0 /
      DATA KTAN2L / -1 /
      SAVE IJAC, IJAC1, IJAC2, LASTAN, KTAN1L, KTAN2L
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C For new tangent path find IJAC index of first and last output levels which
C apply for this path
      IF ( ITAN .NE. LASTAN ) THEN
        IJAC1 = 0
        DO IJAC = 1, NJAC
          IF ( ITNJAC(ITAN,IJAC) .NE. 0 ) THEN
            IF ( IJAC1 .EQ. 0 ) IJAC1 = IJAC
            IJAC2 = IJAC
          END IF
        END DO
        LASTAN = ITAN
      END IF

C Initialise for subsequent calls 
      IF ( LEVINI ) THEN
        IJAC = IJAC1 - 1
        KTAN1L = 0
        KTAN2L = -1      ! Prevents KTAN1,KTAN2 loop from executing
        LEVINI = .FALSE.
      END IF
C
C If next level matches current atmospheric layer then increment KTAN2 to 
C include KTAN index for this level
      IF ( IJAC .LT. IJAC2 ) THEN
        IF ( IDIR * IATM .EQ. ILOJAC(IJAC+1) ) THEN
          IJAC = IJAC + 1
          KTAN = ITNJAC(ITAN,IJAC)
C KTAN1,KTAN2 should span range of currently active TAN indices for o/p
          IF ( KTAN1L .EQ. 0 ) KTAN1L = KTAN
          KTAN2L = KTAN
        END IF
      END IF        
C
C Save KTAN1L and KTAN2L locally since these need to be remembered between
C successive calls to RADTRA, eg between downward and upward limb paths.
      KTAN1 = KTAN1L
      KTAN2 = KTAN2L
C
      END
