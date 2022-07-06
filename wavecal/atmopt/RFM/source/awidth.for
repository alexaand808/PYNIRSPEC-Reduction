      REAL FUNCTION AWIDTH ( IPTH, DWNO )
C
C VERSION
C     10JAN13 AD Use local values of TCOGAS and WIDSTP rather than from /GASCOM/
C     10-AUG-03  AD  Change VLIGHT from SNGL to DBLE  
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Calculate average half-width of line.
C     Called by RFMFIN and RFMWID.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER          IPTH   !  I  Index of Path in /PTHCOM/
      DOUBLE PRECISION DWNO   !  I  Wavenumber [cm]
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'phycon.inc' ! Math and Physical constants
C
C COMMON VARIABLES
      INCLUDE 'pthcom.inc' ! RFM path data
      INCLUDE 'gascom.inc' ! Molecular and isotope data
C
C LOCAL CONSTANTS
      INTEGER       NIDX   ! 1 + No. of defined values for TCOEFF,WIDSTP
        PARAMETER ( NIDX = 29 )      ! Defined for Mol#1-28, use#29 for others
C
C LOCAL VARIABLES
      INTEGER IDX          ! Index of path gas in local arrays
      INTEGER IGAS         ! Index of path gas in /GASCOM/ arrays
      REAL    PRESAV       ! Pressure broadened width
      REAL    DOPPAV       ! Doppler broadened width
      REAL    TCOEFF(NIDX) ! Typical Temp coefficient 
      REAL    WIDSTP(NIDX) ! Typical Half width at STP
C
C DATA STATEMENTS
      DATA TCOEFF / 
     &    6.40023E-01, 7.52456E-01, 7.59870E-01, 7.79879E-01,         ! 1-4
     &    6.90001E-01, 8.62794E-01, 5.00000E-01, 5.00000E-01,         ! 5-8
     &    5.00000E-01, 5.00000E-01, 5.00000E-01, 5.00000E-01,         ! 9-12
     &    5.00000E-01, 5.00000E-01, 6.11556E-01, 5.00000E-01,         ! 13-16
     &    5.00000E-01, 5.00000E-01, 5.00000E-01, 5.00000E-01,         ! 17-20
     &    5.00000E-01, 5.00000E-01, 5.00000E-01, 5.00000E-01,         ! 21-24
     &    5.00000E-01, 7.50000E-01, 5.00000E-01, 5.00000E-01,         ! 25-28
     &    0.5 /                       ! default for all other molecules
C
      DATA WIDSTP / 
     &    6.25677E-02, 7.11645E-02, 6.41407E-02, 7.44450E-02,         ! 1-4
     &    6.01508E-02, 5.84041E-02, 5.20206E-02, 5.27166E-02,         ! 5-8
     &    1.32993E-01, 6.53849E-02, 7.70492E-02, 1.29943E-01,         ! 9-12
     &    8.29986E-02, 2.00000E-02, 4.74378E-02, 4.59445E-02,         ! 13-16
     &    3.97692E-02, 8.49993E-02, 6.99998E-02, 1.08001E-01,         ! 17-20
     &    5.99988E-02, 6.00000E-02, 1.03198E-01, 8.00014E-02,         ! 21-24
     &    1.00002E-01, 7.08757E-02, 1.00000E-01, 7.49980E-02,         ! 25-28
     &    0.1 /                       ! default for all other molecules
C
      SAVE TCOEFF, WIDSTP
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      IGAS   = IGSPTH(IPTH)
      IDX = MIN ( IDXGAS(IGAS), NIDX )
      PRESAV = WIDSTP(IDX) * PREPTH(IPTH) / PREREF * 
     &         ( TEMREF /TEMPTH(IPTH) )**TCOEFF(IDX)
      DOPPAV = SNGL ( DWNO / VLIGHT ) * 
     &         SQRT ( R2 * TEMPTH(IPTH) / WGTGAS(1,IGAS))
      AWIDTH = 0.534310785 * PRESAV + 
     &         SQRT ( 0.216866444 * PRESAV**2 + DOPPAV**2 )
C
      END                                                              
