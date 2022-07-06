      SUBROUTINE NTEPTH ( IJAC, LPTH )
C
C VERSION
C     03-MAY-05  AD  Original.
C
C DESCRIPTION
C     Determine paths affected by Non-LTE Jacobian.
C     Called by RFMPTB for each non-LTE Jacobian.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      IJAC   !  I  Index of Jacobian element (1:NJAC+1).
      INTEGER      LPTH   !  I  Number of unperturbed paths to consider
C
C GLOBAL CONSTANTS
      INCLUDE 'jdxcon.inc' ! Indices for non-VMR Jacobians
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'jaccom.inc' ! Jacobian data
      INCLUDE 'ntecom.inc' ! Non-LTE data
      INCLUDE 'pthcom.inc' ! RFM path data
C
C LOCAL VARIABLES
      INTEGER  IDG    ! HITRAN index of gas for non-LTE data
      INTEGER  IGAS   ! Index of species to be perturbed
      INTEGER  INTE   ! Index of non-LTE data
      INTEGER  IPTH   ! Counter over paths
      INTEGER  ISO    ! HITRAN isotope# for non-LTE data
C
C EXECUTABLE CODE ------------------------------------------------------------
C
      INTE = IGSJAC(IJAC) - ( MAXGAS + JDXNTE )
      IDG  = IDGNTE(INTE)
      ISO  = ISONTE(INTE)
      IGAS = IGSMOL(IDG)
      IGAS = ISOGAS(ISO,IGAS)   ! May or may not be separate index for isotope

      DO IPTH = 1, LPTH
        IF ( IGSPTH(IPTH) .EQ. IGAS .AND.
     &       IATPTH(IPTH) .GE. ILOJAC(IJAC) .AND.
     &       IATPTH(IPTH) .LT. IUPJAC(IJAC)        ) THEN
          IVJPTH(IPTH) = IJAC
        ELSE
          IVJPTH(IPTH) = 0
        END IF
      END DO
C
      END
