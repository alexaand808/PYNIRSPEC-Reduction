C JDXCON.INC ------------------------------------------------------------------
C
C VERSION
C     03-MAY-05  AD  Original
C
C DESCRIPTION
C     Indices for non-VMR Jacobians
C     IGSJAC = 1:MAXGAS indicates the VMR Jacobian for gas IGAS=IGSJAC 
C     IGSJAC = MAXGAS + Index below for additional types of Jacobian
C     (upward and downward path intermediate outputs aren't actually jacobians
C      but use the same arrays)
C
      INTEGER       JDXTEM        ! Temperature
        PARAMETER ( JDXTEM = 1 )
      INTEGER       JDXPRE        ! Pressure
        PARAMETER ( JDXPRE = 2 )
      INTEGER       JDXSFT        ! Surface temperature
        PARAMETER ( JDXSFT = 3 )
      INTEGER       JDXSFE        ! Surface emissivity
        PARAMETER ( JDXSFE = 4 )
      INTEGER       JDXLVU        ! Upward path intermediate output
        PARAMETER ( JDXLVU = 5 )    
      INTEGER       JDXLVD        ! Downward path intermediate output
        PARAMETER ( JDXLVD = 6 )    
      INTEGER       JDXNTE           
        PARAMETER ( JDXNTE = 10 ) ! Offset for vib. temperature Jacobians
C                                 ! IJAC = INTE + MAXGAS + JDXNTE
C
C -----------------------------------------------------------------------------
