      SUBROUTINE NADPTH ( FAIL, ERRMSG )
C
C VERSION
C     17-FEB-04  AD  Initialise NLNPTH
C     04-NOV-01  AD  If ZEN+OBS, limit number of path segments calculated
C     18-MAR-01  AD  Explicitly set FAIL = FALSE before exit
C     11-AUG-99  AD  Remove redundant flgcom.inc
C     01-FEB-99  AD  Set ISCPTH = path index if not scaled (was =0).
C     14-JUL-98  AD  Comment change only.
C     22-JAN-98  AD  Remove redundant EXTERNAL declarations: PTHSRT
C     03-MAR-97  AD  Version 3.
C     14-JAN-97  AD  Original.
C
C DESCRIPTION
C     Construct list of RFM nadir paths from atmos.profiles.
C     Called by INPCHK if NAD or ZEN flags selected.
C     A path is defined for each combination of absorber and atmos layer.
C
      IMPLICIT NONE
C
C ARGUMENTS
      LOGICAL      FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG !  O  Contains error message if FAIL is TRUE.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data  
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'obscom.inc' ! Observer Position
      INCLUDE 'pthcom.inc' ! RFM path data
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL VARIABLES
      INTEGER  IGAS   !  Gas counter
      INTEGER  IPTH   !  Local path counter
      INTEGER  ITAN   !  Tangent path counter
      INTEGER  IATM   !  Atmospheric layer counter
      INTEGER  IATM1  !  Lowest layer# required for path calculation
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Special case of upward viewing with observer within atmosphere (ie above
C the surface so IATOBS in range 1:NATM-1): only need layers# IATOBS:NATM-1
      IF ( ZENFLG .AND. OBSFLG ) THEN
        NCLC = NGAS * ( NATM - IATOBS )
        IATM1 = IATOBS
      ELSE                          ! otherwise whole atmosphere required
        NCLC = NGAS * ( NATM - 1 )
        IATM1 = 1
      END IF
      IF ( NCLC .GT. MAXCLC ) THEN
        WRITE ( ERRMSG, '(A,I11,A,I7)' )
     &    'F-NADPTH: No.paths for line-by-line calc=', NCLC,
     &    ' >MAXCLC in RFMSIZ=', MAXCLC
        FAIL = .TRUE.
        RETURN
      END IF
C 
      NPTH = NCLC * NTAN
      IF ( NPTH .GT. MAXPTH ) THEN
        WRITE ( ERRMSG, '(A,I11,A,I7)' )
     &    'F-NADPTH: Total no.paths required calc=', NPTH,
     &    ' >MAXPTH in RFMSIZ=', MAXPTH
        FAIL = .TRUE.
        RETURN
      END IF
C
      IPTH = 0
      DO ITAN = 1, NTAN
        IATTAN(ITAN) = 0
        DO IGAS = 1, NGAS
          DO IATM = IATM1, NATM-1
            IPTH = IPTH + 1
            IGSPTH(IPTH) = IGAS
            ITNPTH(IPTH) = ITAN
            IATPTH(IPTH) = IATM
            NLNPTH(IPTH) = 0
            IF ( ITAN .EQ. 1 ) THEN
              CLCPTH(IPTH) = .TRUE.
              ISCPTH(IPTH) = IPTH
            ELSE
              CLCPTH(IPTH) = .FALSE.
              ISCPTH(IPTH) = IPTH - (ITAN-1)*NGAS*(NATM-IATM1)
            END IF
          END DO
        END DO
      END DO
C
      FAIL = .FALSE.
C
      END
