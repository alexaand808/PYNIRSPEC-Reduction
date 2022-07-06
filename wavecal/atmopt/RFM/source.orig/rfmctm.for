      SUBROUTINE RFMCTM 
C
C VERSION
C     12-JAN-04  AD  Change CTMH2O to use MT_CKD, and add CTMCKD for old vers.
C     11-DEC-97  AD  Correction to comment: called once for entire widemesh
C     01-DEC-97  AD  Add check for CTMGAS=.TRUE.
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C 
C DESCRIPTION    
C     Calculate continuum absorption on wide mesh.
C     Called by RFM once if CTM flag enabled.
C     Continuum data available for absorbers H2O, O2, N2, CO2.
C     
      IMPLICIT NONE
C
      EXTERNAL
     &  CTMCKD ! Calculate CKD H2O continuum absorption
     &, CTMCO2 ! Calculate CO2 continuum absorption
     &, CTMH2O ! Calculate H2O continuum absorption
     &, CTMO2  ! Calculate O2  continuum absorption
     &, CTMN2  ! Calculate N2  continuum absorption
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'idxcon.inc' ! HITRAN molecular indices
C
C COMMON VARIABLES
      INCLUDE 'pthcom.inc' ! RFM path data
      INCLUDE 'gascom.inc' ! Molecular and isotope data
C
C LOCAL VARIABLES
      INTEGER  IDGPTH ! HITRAN Gas ID for path absorber
      INTEGER  IGAS   ! Gas counter
      INTEGER  IPTH   ! Gas/segment Path counter
      LOGICAL  USECKD ! F=use new MT_CKD H2O continuum, T=use old CKD continuum
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      USECKD = .FALSE.  
C
      DO IPTH = 1, NCLC
        IGAS = IGSPTH(IPTH)
        IF ( CTMGAS(IGAS) ) THEN
          IDGPTH = IDXGAS(IGAS)
          IF ( IDGPTH .EQ. IDXH2O ) THEN 
            IF ( USECKD ) THEN 
              CALL CTMCKD ( IPTH )
            ELSE
              CALL CTMH2O ( IPTH )
            END IF
          ELSE IF ( IDGPTH .EQ. IDXCO2 ) THEN
            CALL CTMCO2 ( IPTH )
          ELSE IF ( IDGPTH .EQ. IDXO2 ) THEN 
            CALL CTMO2 ( IPTH )
          ELSE IF ( IDGPTH .EQ. IDXN2 ) THEN
            CALL CTMN2 ( IPTH )
          END IF
        END IF
      END DO
C
      END
