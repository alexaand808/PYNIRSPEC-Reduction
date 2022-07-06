      REAL FUNCTION WIDGL2 ( PRE, TEM )
C
C VERSION
C     10-AUG-03  AD  Change VLIGHT from SNGL to DBLE
C     11-AUG-99  AD  Add SNGL ( ) around WMDSPC
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Calculate typical average half-width of line.
C     Called by ATMLAY.
C     NB: this is for a "typical" gas whereas AWIDTH is for specific gases
C
      IMPLICIT NONE
C
C ARGUMENTS
      REAL     PRE  !  I  Pressure [mb]
      REAL     TEM  !  I  Temperature [K]
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
C
C LOCAL CONSTANTS
      REAL          EPS           ! 0.0990*ln2
        PARAMETER ( EPS = 0.0686 )
      REAL          TCOTYP        ! Typical temperature coeff for Press.Broad
        PARAMETER ( TCOTYP = 0.5 )
      REAL          WGTTYP        ! Typical molecular weight [g]
        PARAMETER ( WGTTYP = 36.0 )
      REAL          WIDTYP        ! Typical standard line width [/cm]
        PARAMETER ( WIDTYP = 0.1 )
C
C COMMON VARIABLES
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
C
C LOCAL VARIABLES
      REAL    PREWID !  Pressure broadened width
      REAL    DOPWID !  Doppler broadened width
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      PREWID = WIDTYP * PRE / ATMB * ( TEMREF / TEM )**TCOTYP
      DOPWID = SNGL ( WMDSPC / VLIGHT ) * SQRT ( R2 * TEM / WGTTYP )
      WIDGL2 = 0.5 * PREWID * ( 1.0 + EPS ) + 
     &         SQRT ( 0.25 * PREWID**2 * (1.0-EPS)**2 + DOPWID**2 )
C
      END                                                              
