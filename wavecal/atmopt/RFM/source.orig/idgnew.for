      INTEGER FUNCTION IDGNEW ( IDGOLD )
C
C VERSION
C     07-AUG-13  AD  Original
C
C DESCRIPTION
C     Convert old RFM index for .xsc data to new value
C     General purpose RFM function.
C     Returns IDGNEW = 0 if IDGOLD not listed.
C     Returns IDGNEW = IDGOLD if IDGOLD GE 100 (assumes already new index)
C     See molidx.for for definitive list of new index assignments
C
      IMPLICIT NONE 
C
C ARGUMENTS
      INTEGER IDGOLD !  I  Old index
C 
C GLOBAL CONSTANTS
      INCLUDE 'idxcon.inc' ! HITRAN molecular indices
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      IDGNEW = 0
      IF ( IDGOLD .EQ. 50 ) IDGNEW = IDXAER     ! aerosol
      IF ( IDGOLD .EQ. 51 ) IDGNEW = 111        ! f11
      IF ( IDGOLD .EQ. 52 ) IDGNEW = 112        ! f12
      IF ( IDGOLD .EQ. 53 ) IDGNEW = 117        ! f13
      IF ( IDGOLD .EQ. 54 ) IDGNEW = 118        ! f14
      IF ( IDGOLD .EQ. 55 ) IDGNEW = 111        ! f21
      IF ( IDGOLD .EQ. 56 ) IDGNEW = 122        ! f22
      IF ( IDGOLD .EQ. 57 ) IDGNEW = 113        ! f113
      IF ( IDGOLD .EQ. 58 ) IDGNEW = 114        ! f114
      IF ( IDGOLD .EQ. 59 ) IDGNEW = 115        ! f115
      IF ( IDGOLD .EQ. 60 ) IDGNEW = 104        ! ccl4
      IF ( IDGOLD .EQ. 61 ) IDGNEW = 101        ! clono2
      IF ( IDGOLD .EQ. 62 ) IDGNEW = 102        ! n2o5
      IF ( IDGOLD .EQ. 63 ) IDGNEW = 105        ! hno4
      IF ( IDGOLD .EQ. 64 ) IDGNEW = 103        ! sf6
C
      IF ( IDGOLD .EQ. 70 ) IDGNEW = 123        ! f123
      IF ( IDGOLD .EQ. 71 ) IDGNEW = 124        ! f124
      IF ( IDGOLD .EQ. 72 ) IDGNEW = 125        ! f141b
      IF ( IDGOLD .EQ. 73 ) IDGNEW = 126        ! f142b
      IF ( IDGOLD .EQ. 74 ) IDGNEW = 127        ! f225ca
      IF ( IDGOLD .EQ. 75 ) IDGNEW = 128        ! f225cb
      IF ( IDGOLD .EQ. 76 ) IDGNEW = 132        ! f32
      IF ( IDGOLD .EQ. 77 ) IDGNEW = 131        ! f125
      IF ( IDGOLD .EQ. 78 ) IDGNEW = 133        ! f134
      IF ( IDGOLD .EQ. 79 ) IDGNEW = 134        ! f134a
      IF ( IDGOLD .EQ. 80 ) IDGNEW = 135        ! f143a
      IF ( IDGOLD .EQ. 81 ) IDGNEW = 136        ! f152a
      IF ( IDGOLD .EQ. 82 ) IDGNEW = 116        ! f116
      IF ( IDGOLD .EQ. 83 ) IDGNEW = 106        ! sf5cf3
      IF ( IDGOLD .EQ. 84 ) IDGNEW = 145        ! pan
      IF ( IDGOLD .EQ. 85 ) IDGNEW = 142        ! ch3cn
      IF ( IDGOLD .EQ. 86 ) IDGNEW = 153        ! c6h6
      IF ( IDGOLD .EQ. 87 ) IDGNEW = 151        ! c2h6
      IF ( IDGOLD .EQ. 88 ) IDGNEW = 152        ! c3h8
      IF ( IDGOLD .EQ. 89 ) IDGNEW = 144        ! acetone
C
      IF ( IDGOLD .GE. 100 ) IDGNEW = IDGOLD    ! no change, already new index
      END
