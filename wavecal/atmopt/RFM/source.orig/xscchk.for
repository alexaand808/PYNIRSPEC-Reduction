      SUBROUTINE XSCCHK ( FAIL, ERRMSG )
C
C VERSION
C     13-JUN-02  AD  Only check for overlaps, not on whether files exist.
C     26-APR-00  AD  Require some XSC data to be supplied for each gas
C     06-SEP-99  AD  Check for index .GE. 50 rather than 51
C     10-JUL-98  AD  Modify text of error message
C     22-JAN-98  AD  (Standard F77) Put DATA statements *after* declarations
C     03-MAR-97  AD  Version 3.
C     26-FEB-97  AD  Change no X/S data from fatal to warning
C     16-JAN-97  AD  Distinguish between fatal errors due to no X/S data for 
C                    gas, and no X/S data within required range
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Check reqd X/S species uniquely represented in X/Sec files
C     Called once by INPXSC.
C
      IMPLICIT NONE
C
C ARGUMENTS
      LOGICAL      FAIL   !  O  T = a fatal error has been detected
      CHARACTER*80 ERRMSG !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  RFMLOG ! Write text message to RFM log file.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'shpcon.inc' ! Line Shape codes
C
C COMMON VARIABLES
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'xflcom.inc' ! X-section data files
      INCLUDE 'xsccom.inc' ! X/S data sets
C
C LOCAL VARIABLES
      INTEGER      IGAS   !  Counter for required gases
      INTEGER      IXSC   !  Counter for X/S data sets
      LOGICAL      ANYXSC !  T = gas found in at least one X/C data set 
      DOUBLE PRECISION WNOLIM ! Highest Wavenumber [cm] for gas so far
                  SAVE WNOLIM
      CHARACTER*80     WRNMSG  !  Text for warning messages
C DATA STATEMENTS
        DATA WNOLIM / 0.0 /   ! Only to satisfy FLINT - value never used
C
C EXECTUABLE CODE -------------------------------------------------------------
C
C For each gas, check through sorted list for any overlaps
C            
      DO IGAS = 1, NGAS
        IF ( SHPGAS(IGAS) .EQ. SHPXSC ) THEN
          ANYXSC = .FALSE.
          DO IXSC = 1, NXSC
            IF ( IGSXSC(IXSC) .EQ. IGAS ) THEN
              IF ( ANYXSC .AND. WN1XSC(IXSC) .LT. WNOLIM ) THEN
                FAIL = .TRUE.
                ERRMSG ='F-XSCCHK: Overlapping X/S files for Gas:'
     &            //CODGAS(IGAS)
                RETURN
              ELSE
                WNOLIM = DBLE ( WN2XSC(IXSC) )
                ANYXSC = .TRUE.
              END IF
            END IF
          END DO
          IF ( .NOT. ANYXSC ) THEN
            WRNMSG ='W-XSCCHK: No X/S data within reqd range for Gas:'
     &        //CODGAS(IGAS)
            CALL RFMLOG ( WRNMSG, FAIL, ERRMSG )
            IF ( FAIL ) RETURN
          END IF
        END IF
      END DO
C
      FAIL = .FALSE.
      END
