      SUBROUTINE FLXEFN ( ABSLEV, IATM )
C
C VERSION
C     16-JAN-04  AD  Original.
C
C DESCRIPTION
C     Calculate absorption coefficient at atmospheric level.
C     Called by RFMFLX for each output level with EFN flag.
C
      IMPLICIT NONE
C
C ARGUMENTS
      DOUBLE PRECISION ABSLEV(*) !  O  Absorption coefficient [/cm2]
      INTEGER IATM               !  I  Index of atmospheric layer 
C
      EXTERNAL
     &  IDXPTH ! Indexing function for Pth# from Tan.pt#, Atm.Layer# and Gas#
      INTEGER IDXPTH
C
C GLOBAL CONSTANTS
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'gascom.inc' ! Molecular and isotope data
C
C LOCAL VARIABLES
      INTEGER IFIN   ! Fine mesh grid point counter (1:NFIN)
      INTEGER IGAS   ! Absorber counter
      INTEGER IPTH   ! Path counter
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Initialise 
      DO IFIN = 1, NFIN
        ABSLEV(IFIN) = 0.0D0
      END DO
C
C Sum over absorbing species
      DO IGAS = 1, NGAS
        IPTH = IDXPTH ( 1, IATM, IGAS, 2 )
        IF ( IPTH .EQ. 0 ) STOP 'F-FLXEFN: Logical error'
C
        DO IFIN = 1, NFIN
          ABSLEV(IFIN) = ABSLEV(IFIN) + 
     &                   DBLE ( ABSFIN(IFIN,IPTH) / AVOG )
        END DO
      END DO
C
      END 
