C HFLCOM.INC ------------------------------------------------------------------
      COMMON / HFLCOM / WN1HFL, WN2HFL, WNOHFL,
     &                  IR1HFL, IR2HFL, IRCHFL, IFPHFL, LUNHFL
C
C VERSION
C     11-JUN-02  AD  Reduce to single file rather than array of files   
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Common variables associated with RFM HITRAN input file
C
      DOUBLE PRECISION WN1HFL  ! Lowest WNO in each file
      DOUBLE PRECISION WN2HFL  ! Highest WNO in each file
      DOUBLE PRECISION WNOHFL  ! WNO of last line loaded into cyc.buffer
      INTEGER          IR1HFL  ! First data record of File
      INTEGER          IR2HFL  ! Last data record of File
      INTEGER          IRCHFL  ! Current record no 
      INTEGER          IFPHFL(MAXGAS)  ! Forward pointer 
      INTEGER          LUNHFL  ! LUN
C
C -----------------------------------------------------------------------------
