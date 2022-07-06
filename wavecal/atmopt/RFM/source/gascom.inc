C GASCOM.INC ------------------------------------------------------------------
C
      COMMON / GASCOM / NGAS, IAXGAS, IDIGAS, IDXGAS, IGSMOL, ISOGAS, 
     &                  NISGAS, SHPGAS, CHGGAS, CTMGAS, ISOMOL, NTEGAS,
     &                  WGTGAS, 
     &                  CODGAS
C
C VERSION
C     10JAN14 AD Remove TCOGAS, WIDGAS
C     01MAY07 AD (MORSE) Add NTEGAS
C     21MAR04 AD Add ISOMOL
C     06FEB04 AD Add IDIGAS 
C     16FEB02 AD Add ISOGAS
C     28MAR01 AD Add IAXGAS
C     22MAR01 AD Add CHGGAS
C     01DEC97 AD Add CTMGAS 
C     03MAR97 AD Version 3.
C     01OCT96 AD Version 2.
C     01SEP96 AD Version 1.
C
C DESCRIPTION
C     List of gases and associated isotope masses and abundances.
C     Loaded by INPGAS and GASISO.
C     IAXGAS is used to flag whether the aerosol profile requires conversion
C     from extinction to extinction/density (in routine ATMAUX)
C     IAXGAS=0 no aerosol, or aerosol already converted 
C     IAXGAS=1:MAXGAS index of aerosol profile requiring conversion
C
      INTEGER     NGAS                   ! No.species to be used
      INTEGER     IAXGAS                 ! Index of aerosol for conversion
      INTEGER     IDIGAS(MAXGAS)         ! HITRAN isotope index (or 0)
      INTEGER     IDXGAS(MAXGAS)         ! HITRAN index for each species 
      INTEGER     IGSMOL(MAXMOL)         ! IGAS for each HITRAN index (or 0).
      INTEGER     ISOGAS(0:MAXISO,MAXGAS)  ! IDXGAS value for each isotope
      INTEGER     NISGAS(MAXGAS)         ! No.isotopes for species
      INTEGER     SHPGAS(MAXGAS)         ! Lineshape to be used for gas
      LOGICAL     CHGGAS(MAXGAS)         ! T=vmr profile has changed
      LOGICAL     CTMGAS(MAXGAS)         ! T=use continuum for this gas
      LOGICAL     ISOMOL(MAXMOL)         ! T=molecule split by isotopes
      LOGICAL     NTEGAS(MAXGAS)         ! T=use non-LTE for molecule
      REAL        WGTGAS(MAXISO,MAXGAS)  ! Molecular wts [Atomic units] 
      CHARACTER*7 CODGAS(MAXGAS)         ! Character codes for species
C
C -----------------------------------------------------------------------------
