      SUBROUTINE HOMPTH ( FAIL, ERRMSG )
C
C VERSION
C     17-FEB-04  AD  Initialise NLNPTH
C     19-MAR-01  AD  Explicitly set FAIL=FALSE before exit
C     06-SEP-00  AD  Correction: only set NPTH,NCLC on first call 
C     01-FEB-99  AD  Set ISCPTH = same path index if not scaled (was 0). 
C     01-JAN-99  AD  Scale all paths except for first 'TAN' (path length)
C     04-SEP-98  AD  Correction: set path length from USRTAN not HGTTAN
C     14-JUL-98  AD  Comment change only.
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Set up homogeneous paths for RFM calculations.
C     Called by INPCHK if HOM flag enabled, RFMPTB if JAC flag also enabled.
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
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
C
C COMMON VARIABLES
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'atmcom.inc' ! Atmospheric profile data  
      INCLUDE 'tancom.inc' ! Tangent heights
      INCLUDE 'pthcom.inc' ! RFM path data
C
C LOCAL VARIABLES
      INTEGER  IGAS   !  Gas counter
      INTEGER  IPTH   !  Path segment counter
      INTEGER  ITAN   !  Tangent path counter
      LOGICAL  FIRST  !  T=first call, F=subsequent calls
C
C DATA STATEMENTS
      DATA FIRST / .TRUE. /
      SAVE FIRST
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Only set/check NPTH, NCLC on first call since any subsequent calls will be
C from RFMPTB which will also increment these values
      IF ( FIRST ) THEN
        NPTH = NGAS * NTAN 
        IF ( NPTH .GT. MAXPTH ) THEN
          FAIL = .TRUE.
          WRITE ( ERRMSG, '(A,I11,A,I7)' )
     &      'F-HOMPTH: No.paths reqd=', NPTH,
     &      ' >MAXPTH in RFMSIZ=', MAXPTH
          RETURN
        END IF
C
        NCLC = NGAS 
        IF ( NCLC .GT. MAXCLC ) THEN
          FAIL = .TRUE.
          WRITE ( ERRMSG, '(A,I11,A,I7)' )
     &      'F-HOMPTH: No.paths for line-by-line calc=', NCLC,
     &      ' >MAXCLC in RFMSIZ=', MAXCLC
          RETURN
        END IF
        FIRST = .FALSE.
      END IF
C
      IPTH = 0
      DO ITAN = 1, NTAN           
        DO IGAS = 1, NGAS
          IPTH = IPTH + 1
          IGSPTH(IPTH) = IGAS
          ITNPTH(IPTH) = ITAN
          IATPTH(IPTH) = 1
          NLNPTH(IPTH) = 0
          TEMPTH(IPTH) = TEMATM(1)
          PREPTH(IPTH) = PREATM(1) / ATMB
          PPAPTH(IPTH) = VMRATM(1,IGAS) * PREPTH(IPTH)
          RAYPTH(IPTH) = USRTAN(ITAN)
          AMTPTH(IPTH) = 1.0E5 / AVOG * VMRATM(1,IGAS) * DNSATM(1) 
     &                   * USRTAN(ITAN)
          IF ( ITAN .EQ. 1 ) THEN
            CLCPTH(IPTH) = .TRUE.
            ISCPTH(IPTH) = IPTH
          ELSE
            CLCPTH(IPTH) = .FALSE.
            ISCPTH(IPTH) = IPTH - ( ITAN - 1 ) * NGAS 
          END IF
        END DO
        IATTAN(ITAN) = 1
      END DO
C
      FAIL = .FALSE.         
      END
