C SFCCOM.INC -------------------------------------------------------------------
      COMMON / SFCCOM / EMSSFC, WNOSFC, NSFC, RFLSFC, TEMSFC
C
C VERSION
C     23-APR-12  AD  Change EMSSFC to array and add WNOSFC, RFLSFC
C     13-JUL-98  AD  Original.
C
C DESCRIPTION
C     Surface parameters.
C     Read from Driver Table if SFC flag is enabled.
C     Loaded/checked by inpsfc.for
C     Accessed by sfcems.for, rfmrad.for and rfmflx.for
C
      DOUBLE PRECISION EMSSFC(MAXSFC) ! Surface Emissivity (0:1)
      DOUBLE PRECISION WNOSFC(MAXSFC) ! Wavenumbers [cm-1] for sfc.emissivity
      INTEGER          NSFC   ! Number of Emissivity values
      LOGICAL          RFLSFC ! T=reflecting surface, F=no reflection
      REAL             TEMSFC ! Surface Temperature [K]
C
C ------------------------------------------------------------------------------
