C HITCOM.INC -------------------------------------------------------------------
      COMMON / HITCOM / WNUM,
     &                  STREN, TPROB, ABROAD, SBROAD, ELS, ABCOEF, TSP, 
     &                  LSTAT, IDGAS, ISO, IUSGQ, ILSGQ, IFWDPT, 
     &                  IUSNTE, ILSNTE, 
     &                  USLQ, BSLQ
C
C VERSION
C      07-AUG-13  AD  Remove AI, REF and BLANK from common data
C      03-MAR-97  AD  Version 3.
C      01-OCT-96  AD  Version 2.
C      20-SEP-96  AD  Rename IUSVIB,ILSVIB to IUSNTE,ILSNTE
C      01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Quantites read from the HITRAN line data via the HITLIN output file.
C     IUSNTE, ILSNTE set by SETNTE.
C
      DOUBLE PRECISION WNUM  ! Line frequency [cm-1] in double precision.
C
      REAL    STREN  ! Line strength  [cm-1./(kg.moles.cm-2)]@296K. 
C                      (HITRAN value * 6.022E26 to avoid underflows.)
      REAL    TPROB  ! Transition probability [Debyes2].
      REAL    ABROAD ! Air-broad halfwidth  (HWHM) [cm-1/atm] @ 296K.
      REAL    SBROAD ! Self-broad halfwidth (HWHM) [cm-1/atm] @ 296K.
      REAL    ELS    ! Lower-state energy [cm-1].
      REAL    ABCOEF ! Coefficient of temp.dependance of air-broadened HW
      REAL    TSP    ! Transition shift due to pressure (presently empty, some 
C                      coupling coeffs. inserted).
      INTEGER LSTAT  ! Status of transition information.
      INTEGER IDGAS  ! Molecule number.(1-32)
      INTEGER ISO    ! Isotope number (1=most abundant, 2=second,).
      INTEGER IUSGQ  ! Upper state global quanta index.
      INTEGER ILSGQ  ! Lower state global quanta index.
      INTEGER IFWDPT ! Forward pointer on data line.
      INTEGER IUSNTE ! Index of non-LTE dataset affecting upper level
      INTEGER ILSNTE ! Index of non-LTE dataset affecting lower level
C
      CHARACTER*9 USLQ   ! Upper state local quanta.
      CHARACTER*9 BSLQ   ! Lower state local quanta.
C
C ------------------------------------------------------------------------------
