      SUBROUTINE INPCHK ( FAIL, ERRMSG )
C
C VERSION
C     07-AUG-13  AD  Remove redundant SPCFIN and RFMLOG from EXTERNAL list
C     03-APR-05  AD  Add CHKVTJ
C     01-JAN-04  AD  Add LEVTAN
C     13-JUN-02  AD  Add ABSCHK
C     16-JUN-01  AD  Set FAIL=FALSE if TAB flag
C     18-MAR-01  AD  Original.
C
C DECSRIPTION
C     Cross-check inputs from more than one drv table section
C     Called once by RFMINP after all driver table info is loaded.
C     This is used to perform checks where inputs come from sections which 
C     can be in arbitrary order, so only meaningful after all driver table
C     data has been loaded.
C
      IMPLICIT NONE
C
C ARGUMENTS
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  ABSCHK ! Check information available for all absorbers.
     &, CHKVTJ ! Check vibrational temperature Jacobians
     &, FINCHK ! Check fine mesh resolution
     &, FLXPTH ! Set up paths for flux calculations.
     &, FOVTAN ! Add tangent heights to allow FOV convolution
     &, HOMPTH ! Set up homogeneous paths for RFM calculations.
     &, ILSCHK ! Assign ILS Functions to each spectral range 
     &, LEVTAN ! Set up tangent paths for intermediate output levels
     &, NADPTH ! Set up nadir paths for RFM calculations.
     &, TANHGT ! Calculate actual (refracted) tangent heights
     &, TANPTH ! Construct list of RFM paths from tang.hts and atmos.profiles.
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      CALL ABSCHK ( FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
      IF ( HOMFLG ) THEN
        CALL HOMPTH ( FAIL, ERRMSG )
      ELSE IF ( FLXFLG ) THEN       ! must be tested before NAD/ZEN flags
        CALL FLXPTH ( FAIL, ERRMSG )
      ELSE IF ( NADFLG .OR. ZENFLG ) THEN
        CALL NADPTH ( FAIL, ERRMSG )
      ELSE IF ( .NOT. TABFLG ) THEN      ! limb-viewing calculation
        CALL TANHGT ( FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        IF ( FOVFLG ) THEN
          CALL FOVTAN ( FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        END IF
        CALL TANPTH ( FAIL, ERRMSG )
      ELSE                               ! TAB calculation - no checks here
        FAIL = .FALSE. 
      END IF
      IF ( FAIL ) RETURN
C
      IF ( ILSFLG ) THEN
        CALL ILSCHK ( FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
C
      IF ( LEVFLG ) THEN
        CALL LEVTAN ( FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
C
C Check that fine resolution is not coarser than any requested output grid
      CALL FINCHK ( FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
      IF ( JACFLG .AND. NTEFLG ) THEN
        CALL CHKVTJ ( FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
C
      END
