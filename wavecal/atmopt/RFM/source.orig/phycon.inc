C PHYCON.INC ------------------------------------------------------------------
C
C VERSION
C     14-SEP-11  AD  Add TEMSPA
C     12-NOV-04  AD  Change VLIGHT from REAL to DP.
C     27-JAN-02  AD  Correct CPMOL (was 34.968), revise some constants 
C     19-APR-00  AD  Add CPMOL
C     29-DEC-99  AD  Convert REARTH from SP to DP  
C     02-AUG-99  AD  Add "D0" exponent to C2 to force Double Precision
C     07-JAN-98  AD  Change ARGMIN from 1.0E-38 to 9.99999E-37
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Mathematical and physical constants used within MIPAS RFM.
C
C     29-JUN-02 following constants were revised to use the 1998 CODATA
C     recommended values (email from T.J.Nightingale, RAL). 
C          AVOG was  6.0221367E26 increased by 0.0001%
C          RGAS was  8314.5112    decreased by 0.0005%
C            R2 was  11526.360    decreased by 0.0005%
C            C1 was  1.1910439D-3 decreased by 0.0001%
C            C2 was  1.4387686D0  increased by 0.0005%
C REFERENCE
C     NIST: http://physics.nist.gov/cuu/Units/index.html
C
      REAL   ARGMIN                       ! Min.argument for LOG(ARGMIN)
        PARAMETER ( ARGMIN = 9.99999E-37 )! 1.0E-38 causes problems on DEC OSF
      REAL   LOGMIN                       ! Min.value for LOG(ARGMIN)
        PARAMETER ( LOGMIN = -87.49823 )
      REAL   ATMB                         ! Pressure conversion mb --> atm
        PARAMETER ( ATMB = 1013.25 )      ! [mb/atm]
      REAL   AVOG                         ! Avogadro constant 
        PARAMETER ( AVOG = 6.02214199E26) ! [/kmole] - CODATA
      REAL   CPMOL                        ! Molec.heat cap of air
        PARAMETER ( CPMOL = 29012. )      ! [J/K/kmole] N2/O2/Ar/CO2 @25C
      REAL   LN2                          ! Ln(2)
        PARAMETER ( LN2 = 0.69314718 )
      REAL   PI                           ! pi
        PARAMETER ( PI = 3.141592654 )
      REAL   DTORAD                       ! deg_to_radians conversion factor
        PARAMETER ( DTORAD = PI/180.0 )  
      REAL   PREREF                       ! Line data reference pressure 
        PARAMETER ( PREREF = 1.0 )        ! [atm]
      REAL   RGAS                         ! Gas constant
        PARAMETER ( RGAS = 8314.472 )     ! [J/K/kmole] - CODATA 
      REAL   RSQRPI                       ! 1.0/Sqrt(pi)
        PARAMETER ( RSQRPI = 0.56418958 )
      REAL   SQRLN2                       ! SQRT ( LN ( 2 ) ) 
        PARAMETER ( SQRLN2 = 0.832554611 )  
      REAL   R2                           ! 2*ln 2*k*N 
        PARAMETER ( R2 = 11526.306)       ! [k.mole] - CODATA
      REAL   TEMREF                       ! Line data reference temperature 
        PARAMETER ( TEMREF = 296.0 )      ! [K]
      REAL   TEMSPA                       ! Cosmic background temperature
        PARAMETER ( TEMSPA = 2.7 )        ! [K]
      DOUBLE PRECISION VLIGHT             ! Velocity of light 
        PARAMETER ( VLIGHT = 2.99792458D8 )! [m/s] - CODATA 
      DOUBLE PRECISION C1                 ! 1st radiation const. 
        PARAMETER ( C1 = 1.191042723D-3 ) ! [nW/(cm2.sr.(cm-1)4)] - CODATA
      DOUBLE PRECISION C2                 ! 2nd radiation const. = hc/k 
        PARAMETER ( C2 = 1.4387752D0 )    ! [K/cm-1] - CODATA
      DOUBLE PRECISION REARTH             ! Mean earth radius [km]
        PARAMETER ( REARTH = 6367.421D0 ) ! Matches GENLN2 calc for USStdAtm.
C
C -----------------------------------------------------------------------------
