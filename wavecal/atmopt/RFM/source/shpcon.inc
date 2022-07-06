C SHPCON.INC ------------------------------------------------------------------
C
C VERSION
C     11-AUG-03  AD  Add SHPVVW
C     01-DEC-97  AD  Add SHPCTM, remove redundant LNOCON,LWTHCN
C     02-JUL-97  AD  Add SHPLUT
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Line shape integer codes.
C     Used in RFMFIN and RFMWID.
C
      INTEGER SHPVOI                ! Voigt line shape 
        PARAMETER ( SHPVOI = 1 )
      INTEGER SHPLOR                ! Lorentz line shape 
        PARAMETER ( SHPLOR = 2 )    
      INTEGER SHPDOP                ! Doppler line shape 
        PARAMETER ( SHPDOP = 3 )   
      INTEGER SHPCHI                ! Chi-factor scaled Voigt lineshape 
        PARAMETER ( SHPCHI = 4 ) 
      INTEGER SHPVVW                ! Van Vleck-Weisskopf lineshape 
        PARAMETER ( SHPVVW = 5 ) 
      INTEGER SHPCTM                ! Continuum data only
        PARAMETER ( SHPCTM = 8 )
      INTEGER SHPXSC                ! Cross-section data 
        PARAMETER ( SHPXSC = 9 )
      INTEGER SHPLUT                ! Offset for LUT data
        PARAMETER ( SHPLUT = 10 )   ! Needs to be highest SHP* value
C
C -----------------------------------------------------------------------------