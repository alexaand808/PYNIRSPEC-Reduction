      INTEGER FUNCTION JDXPTH ( ITAN, IATM, IGAS, IDIR )
C
C VERSION
C     15-JAN-04  AD  Remove test of GRAFLG - always check IDIR
C     10-DEC-01  AD  Save LTAN, LDIR
C     05-NOV-01  AD  Original.
C
C DESCRIPTION
C     Indexing function for Pth# from Tan.pt#, Atm.Layer# and Gas#
C     General purpose RFM function.
C     Same as IDXPTH but with separate memory of ITAN.
C
      IMPLICIT NONE 
C
C ARGUMENTS
      INTEGER   ITAN  !  I  Tangent Ray#
      INTEGER   IATM  !  I  Atmospheric Layer#
      INTEGER   IGAS  !  I  Absorber#
      INTEGER   IDIR  !  I  Direction (ignored unless GRAFLG is TRUE)
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data      
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'pthcom.inc' ! RFM path data
C
C LOCAL VARIABLES
      INTEGER   IDXTAN(MAXATM,MAXGAS) ! Store IPTH indices for ITAN,IDIR
      INTEGER   IPTH  !  Counter
      INTEGER   JATM  !  Atmosphere layer counter
      INTEGER   JGAS  !  Absorber counter
      INTEGER   LDIR  !  Previous IDIR
      INTEGER   LTAN  !  Previous ITAN
C
C DATA STATEMENTS
      DATA LTAN / 0 / 
      DATA LDIR / 0 /
      SAVE LTAN, LDIR, IDXTAN
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      IF ( ITAN .NE. LTAN .OR. IDIR .NE. LDIR ) THEN
        DO JGAS = 1, NGAS
          DO JATM = 1, NATM
            IDXTAN(JATM,JGAS) = 0
          END DO
        END DO
        DO IPTH = 1, NPTH
          IF ( ITNPTH(IPTH) .EQ. ITAN .AND. IDRPTH(IPTH) .EQ. IDIR ) 
     &      IDXTAN(IATPTH(IPTH),IGSPTH(IPTH)) = IPTH
        END DO
        LTAN = ITAN
        LDIR = IDIR
      END IF
      JDXPTH = IDXTAN(IATM,IGAS)
      END
