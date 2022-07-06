C IDXCON.INC ------------------------------------------------------------------
C
C VERSION
C     23-AUG-13  AD  Change IDXAER from 50 to 100
C     09-JAN-08  AD  Remove IDXHDO
C     16-FEB-02  AD  Add IDXHDO=39
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     28-SEP-96  AD  Add IDXAER=50
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Constants associated with HITRAN molecules  required for specific tests.
C
      INTEGER       IDXH2O
        PARAMETER ( IDXH2O = 1 )
      INTEGER       IDXCO2
        PARAMETER ( IDXCO2 = 2 )
      INTEGER       IDXO2
        PARAMETER ( IDXO2 = 7 )
      INTEGER       IDXN2
        PARAMETER ( IDXN2 = 22 )
      INTEGER       IDXAER           
        PARAMETER ( IDXAER = 100 )    ! RFM-specific assignment, not HITRAN
C
C -----------------------------------------------------------------------------
