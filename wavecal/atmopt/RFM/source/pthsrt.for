      SUBROUTINE PTHSRT
C
C VERSION
C     03-MAY-05  AD  Add IVJPTH
C     02-JAN-00  AD  Add IDRPTH, PSIPTH, and remove unused IDIPTH 
C     11-AUG-99  AD  Remove redundant flgcom.inc 
C     02-FEB-99  AD  Comment change: also called by RFMPTH
C                    Remove requirement that at least one scaled path exists
C                    Replace all ISCPTH elements anyway
C     10-JUL-98  AD  Remove redundant IFLPTH
C     03-MAR-97  AD  Version 3.
C     14-JAN-97  AD  Change description text only.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Sort Calculated paths to start of /PTHCOM/ list.
C     Called once by RFMSCA, TANPTH and (if JAC flag) RFMPTB.
C     Shuffles paths until Paths#1:NCLC are the calculated paths,
C     NCLC+1:NPTH are the uncalculated (or scaled) paths.
C
C     This is routine is also necessary if the FOV flag is employed but the
C     nominal tangent path does not constitute one of the FOV convolution paths.
C     In this case it is necessary to move the uncalculated paths (from the 
C     nominal tangent path) to the end of the list NCLC+1:NPTH 
C
C     If Jacobians are calculated, these also add new paths, some of which will
C     require separate calculation and others just scaled by absorber amount.
C
      IMPLICIT NONE
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'pthcom.inc' ! RFM path data
C
C LOCAL VARIABLES
      INTEGER  IPTH   ! Path counter
      INTEGER  JPTH   ! Secondary path counter for shuffle
      INTEGER  KPTH   ! Tertiary path counter
      INTEGER  MPTH   ! Start point for search for next Calc.path beyond NCLC
      INTEGER  IATSAV ! Atmospheric profile layer# for this path
      INTEGER  IDRSAV ! Direction for path 
      INTEGER  IGSSAV ! Absorbing gas# for this path
      INTEGER  ISCSAV ! Index of path to be used for scaling absorption
      INTEGER  ITNSAV ! Tangent Height# for this path
      INTEGER  IVJSAV ! Vibrational Temperature Jacobian
      INTEGER  NLNSAV ! No.lines stored in cyclic buffer 
      REAL     PSISAV ! Horiz.Angle [deg]
      REAL     PRESAV ! Pressure [atm]
      REAL     TEMSAV ! Temperature 
      REAL     PPASAV ! Partial pressure 
      REAL     AMTSAV ! Path gas amount [ kmol / cm^2 ] 
      REAL     RAYSAV ! Ray length in path [km]
      LOGICAL  CLCSAV ! T=perform line-by-line calc.for this path
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      MPTH = NCLC + 1
      DO IPTH = 1, NCLC
        IF ( .NOT. CLCPTH(IPTH) ) THEN             ! IPTH not a calc.path
          DO JPTH = MPTH, NPTH
            IF ( CLCPTH(JPTH) ) THEN               ! swap with IPTH
C
C Replace any path scaling indices for JPTH with IPTH
              DO KPTH = IPTH, NPTH
                IF ( ISCPTH(KPTH) .EQ. JPTH ) ISCPTH(KPTH) = IPTH
              END DO
C
              IATSAV = IATPTH(IPTH) 
              IDRSAV = IDRPTH(IPTH) 
              IGSSAV = IGSPTH(IPTH) 
              ISCSAV = ISCPTH(IPTH) 
              ITNSAV = ITNPTH(IPTH) 
              IVJSAV = IVJPTH(IPTH)
              NLNSAV = NLNPTH(IPTH) 
              PSISAV = PSIPTH(IPTH) 
              PRESAV = PREPTH(IPTH) 
              TEMSAV = TEMPTH(IPTH) 
              PPASAV = PPAPTH(IPTH) 
              AMTSAV = AMTPTH(IPTH) 
              RAYSAV = RAYPTH(IPTH) 
              CLCSAV = CLCPTH(IPTH) 
C
              IATPTH(IPTH) = IATPTH(JPTH)
              IDRPTH(IPTH) = IDRPTH(JPTH) 
              IGSPTH(IPTH) = IGSPTH(JPTH) 
              ISCPTH(IPTH) = ISCPTH(JPTH) 
              ITNPTH(IPTH) = ITNPTH(JPTH) 
              IVJPTH(IPTH) = IVJPTH(JPTH) 
              NLNPTH(IPTH) = NLNPTH(JPTH) 
              PSIPTH(IPTH) = PSIPTH(JPTH) 
              PREPTH(IPTH) = PREPTH(JPTH) 
              TEMPTH(IPTH) = TEMPTH(JPTH) 
              PPAPTH(IPTH) = PPAPTH(JPTH) 
              AMTPTH(IPTH) = AMTPTH(JPTH) 
              RAYPTH(IPTH) = RAYPTH(JPTH) 
              CLCPTH(IPTH) = CLCPTH(JPTH) 
C
              IATPTH(JPTH) = IATSAV
              IDRPTH(JPTH) = IDRSAV 
              IGSPTH(JPTH) = IGSSAV 
              ISCPTH(JPTH) = ISCSAV 
              ITNPTH(JPTH) = ITNSAV 
              IVJPTH(JPTH) = IVJSAV 
              NLNPTH(JPTH) = NLNSAV 
              PSIPTH(JPTH) = PSISAV 
              PREPTH(JPTH) = PRESAV 
              TEMPTH(JPTH) = TEMSAV 
              PPAPTH(JPTH) = PPASAV 
              AMTPTH(JPTH) = AMTSAV 
              RAYPTH(JPTH) = RAYSAV 
              CLCPTH(JPTH) = CLCSAV 
C
C Ensure next search for CLCPTH=T starts at next point in PTH list
              MPTH = JPTH + 1
              GOTO 100
            END IF
          END DO
        END IF
  100   CONTINUE
      END DO
C
      END
