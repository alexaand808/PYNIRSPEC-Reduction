C OBSCOM.INC -------------------------------------------------------------------
      COMMON / OBSCOM / IATOBS, ALTOBS, PSIOBS
C
C VERSION
C     29-DEC-99  AD  Add PSIOBS
C     16-JUL-98  AD  Original.
C
C DESCRIPTION
C     Observer Position.
C     Set by module OBSATM
C     PSIOBS is only defined if GRA flag enabled.
C
      INTEGER  IATOBS ! Atmospheric Layer# with observer at base (1:NATM)
      REAL     ALTOBS ! Observer Altitude (above surface at 0km) [km]
      REAL     PSIOBS ! Observer horizontal angle [deg]
C
C ------------------------------------------------------------------------------
