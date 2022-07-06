      SUBROUTINE RFMDEF 
C
C VERSION
C     14-SEP-11  AD  Add NAMRJT
C     08-JAN-08  AD  Initialise NNTE
C     06-JAN-04  AD  Add NAMPRF
C     09-AUG-03  AD  Add NAMBBT
C     16-APR-00  AD  Add NAMCOO
C     22-OCT-97  AD  Add NAMTAB
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     17-SEP-96  AD  Set default radius of curvature (RADCRV) and NOMFIN
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Set RFM default filenames and Min.line strengths
C     Callled by RFM once
C
      IMPLICIT NONE
C
C GLOBAL CONSTANTS
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
      INCLUDE 'rfmcon.inc' ! Constants for RFM calculations
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES 
      INCLUDE 'crvcom.inc' ! Local radius of curvature
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'ntecom.inc' ! Non-LTE data
      INCLUDE 'outcom.inc' ! RFM output file data
      INCLUDE 'rejcom.inc' ! Minimum Line strength limits
      INCLUDE 'rfmfil.inc' ! Standard filenames of RFM I/O files
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      NAMABS = DEFABS
      NAMBBT = DEFBBT
      NAMCOO = DEFCOO
      NAMOPT = DEFOPT
      NAMPRF = DEFPRF
      NAMPTH = DEFPTH
      NAMRAD = DEFRAD
      NAMRJT = DEFRJT
      NAMTAB = DEFTAB
      NAMTRA = DEFTRA 
      NAMWID = DEFWID
C
      WIDREJ = 0.0
      FINREJ = 0.0      
C
      NOMFIN = DEFFIN
C
      NNTE = 0     ! Initialise here since NTE flag can be used without *NTE sec
C
      RADCRV = REARTH
C
      END     
