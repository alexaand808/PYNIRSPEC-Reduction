C RFMCON.INC ------------------------------------------------------------------
C
C VERSION
C     03-MAY-05  AD  Add PTBNTE
C     03-SEP-03  AD  Add PTBSFT, PTBSFE
C     08-JUN-01  AD  Add PTBLOS
C     20-APR-00  AD  Add PTBSRC
C     07-OCT-99  AD  Remove WIDILS
C     02-AUG-99  AD  Add "D0" to Double Precision constants
C     16-JUN-99  AD  Add PTBAER
C     23-JAN-99  AD  Change LENLAB from 6 to 8
C     01-JAN-99  AD  Add PTBVMR, PTBTEM, PTBPRE
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     17-SEP-96  AD  Rename NOMFIN to DEFFIN
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Constants associated with RFM spectral calculations
C
      DOUBLE PRECISION FEXC           ! Window for fine mesh calcs
        PARAMETER (    FEXC = 1.0D0 )
      DOUBLE PRECISION FWIND          ! Window for wide mesh calcs
        PARAMETER (    FWIND = 25.0D0 )
      INTEGER       LENLAB            ! Length of Spc.range label
        PARAMETER ( LENLAB = 8 )
      INTEGER       DEFFIN            ! Default Nominal No.pts/cm-1 fine mesh
        PARAMETER ( DEFFIN = 2000 )   
      REAL          NWDCUT            ! No. Av.Voigt halfwidths for cut-off
        PARAMETER ( NWDCUT = 200.0 )  ! 
      REAL          PTBAER            ! Amount [/km] for aerosol perturbation
        PARAMETER ( PTBAER = 1.0E-4 ) ! 
      REAL          PTBLOS            ! Fraction for LOS perturbation
        PARAMETER ( PTBLOS = 0.001 )  ! 0.001 = 1/1000 x units in *TAN section
      REAL          PTBNTE            ! Amount [K] for vib.temp perturbation
        PARAMETER ( PTBNTE = 1.0 )    ! 1.0 = 1 K
      REAL          PTBPRE            ! Fraction for pressure perturbation
        PARAMETER ( PTBPRE = 0.01 )   ! 0.01 = 1%
      REAL          PTBSFE            ! Amount [0-1] for surface emissivity ptb
        PARAMETER ( PTBSFE = 0.01 )
      REAL          PTBSFT            ! Amount [K] for surface temperature ptb
        PARAMETER ( PTBSFT = 1.0 )
      REAL          PTBSRC            ! Fraction for Source Fn perturbation
        PARAMETER ( PTBSRC = 0.01 )   ! 0.01 = 1%
      REAL          PTBTEM            ! Amount [K] for temperature perturbation
        PARAMETER ( PTBTEM = 1.0 )    ! 1.0 = 1 K
      REAL          PTBVMR            ! Fraction for VMR perturbation
        PARAMETER ( PTBVMR = 0.01 )   ! 0.01 = 1%
C
C -----------------------------------------------------------------------------
