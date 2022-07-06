      SUBROUTINE RFMRAD 
C
C VERSION
C     07-AUG-13  AD  Remove redundant local variable REFLCT
C     23-APR-12  AD  Remove REFLCT argument from RADSFC, instead use 
C                    RFLSFC in sfccom.inc
C     14-SEP-11  AD  Add RADSPA for cosmic background radiance
C     03-SEP-03  AD  Change SFCTAN from local to common variable
C     12-OCT-99  AD  Check each limb tan.path for surface intersection
C     31-JUL-99  AD  Remove unused variable ISGN 
C     30-DEC-98  AD  Split into new subroutines RADTRA and RADSFC
C                     Implement BFX flag.
C                     Correct to allow surface effects in limb-viewing mode
C     14-JUL-98  AD  Allow Observer (OBS) option
C     13-JUL-98  AD  Allow surface (SFC) option for NAD and HOM paths
C     04-APR-97  AD  Correction to allow for laser-emission - use SIGN to keep
C                     sign of ABSLAY to allow negative values
C     03-MAR-97  AD  Version 3.
C     14-JAN-97  AD  Add Zenith-viewing option
C     24-DEC-96  AD  Calculate FACTOR checking for AMTPTH(JPTH)=0 
C     12-DEC-96  AD  Revised for more accurate calc.with multiple absorbers
C     01-OCT-96  AD  Version 2.
C     20-SEP-96  AD  Rename VIBFLG to NTEFLG
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Perform RFM radiance calculation
C     Called by RFM for each widemesh interval.
C
      IMPLICIT NONE
C
      EXTERNAL
     &  RADSFC ! Surface contribution to radiative transfer calculation
     &, RADSPA ! Cosmic background contribution to radiative transfer calc.
     &, RADTRA ! Radiative transfer calculation for single ray
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data 
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'obscom.inc' ! Observer Position
      INCLUDE 'sfccom.inc' ! Surface parameters
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL VARIABLES
      INTEGER   IATMIN   !  Lowest layer# (either tang.layer or bottom layer)
      INTEGER   IATTOP   !  Top layer# of atmosphere
      INTEGER   ITAN     !  Tangent point counter
      LOGICAL   INIT     !  T=initialise for start of new ray path
      LOGICAL   UPVIEW   !  T=upward limb-view
       
C EXECUTABLE CODE -------------------------------------------------------------
C
      IATTOP = NATM-1
      DO ITAN = 1, MTAN
        IF ( CLCTAN(ITAN) ) THEN
C
C All calculations start at the observer and proceed away into the atmosphere
C
          INIT = .TRUE.      ! Set .FALSE. after calling RADTRA
C
          IF ( HOMFLG ) THEN
            CALL RADTRA ( ITAN, 1, 1, 1, INIT )
            IF ( SFCFLG ) THEN
              CALL RADSFC ( ITAN )
              IF ( RFLSFC ) CALL RADTRA ( ITAN, 1, 1, 1, INIT )
            END IF
C
          ELSE IF ( ZENFLG ) THEN 
            IF ( OBSFLG ) THEN   
              CALL RADTRA ( ITAN, IATOBS, IATTOP, 1, INIT )
            ELSE
              CALL RADTRA ( ITAN, 1, IATTOP, 1, INIT )
            END IF
            CALL RADSPA ( ITAN )
C
          ELSE IF ( NADFLG ) THEN
            IF ( OBSFLG ) THEN
              CALL RADTRA ( ITAN, (IATOBS-1), 1, -1, INIT )
            ELSE
              CALL RADTRA ( ITAN, IATTOP, 1, -1, INIT )
            END IF
            CALL RADSFC ( ITAN ) ! SFC flag=T for all NAD calculations
            IF ( RFLSFC ) THEN
              CALL RADTRA ( ITAN, 1, IATTOP, 1, INIT )
              CALL RADSPA ( ITAN )
            END IF
          ELSE                           ! limb viewing
            IATMIN = IATTAN(ITAN)
            IF ( OBSFLG ) THEN
              UPVIEW = IATMIN .EQ. IATOBS
              IF ( UPVIEW ) THEN
                CALL RADTRA ( ITAN, IATOBS, IATTOP, 1, INIT )
                CALL RADSPA ( ITAN ) 
              ELSE
                CALL RADTRA ( ITAN, (IATOBS-1), IATMIN, -1, INIT )
                IF ( SFCTAN(ITAN) ) THEN
                  CALL RADSFC ( ITAN )
                  IF ( RFLSFC ) THEN
                    CALL RADTRA ( ITAN, 1, IATTOP, 1, INIT )
                    CALL RADSPA ( ITAN ) 
                  END IF
                ELSE
                  CALL RADTRA ( ITAN, IATMIN, IATTOP, 1, INIT )
                  CALL RADSPA ( ITAN ) 
                END IF
              END IF
            ELSE
              CALL RADTRA ( ITAN, IATTOP, IATMIN, -1, INIT )
              IF ( SFCTAN(ITAN) ) THEN
                CALL RADSFC ( ITAN )
                IF ( RFLSFC ) THEN
                  CALL RADTRA ( ITAN, 1, IATTOP, 1, INIT )
                  CALL RADSPA ( ITAN ) 
                END IF
              ELSE
                CALL RADTRA ( ITAN, IATMIN, IATTOP, 1, INIT )
                CALL RADSPA ( ITAN ) 
              END IF
            END IF
C
          END IF
        END IF
      END DO
C
      END 
