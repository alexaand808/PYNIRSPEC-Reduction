      SUBROUTINE MOLIDX ( IDX, MOL, FAIL, ERRMSG )
C
C VERSION
C     10JAN14 AD Add #61 BrOq
C     24OCT13 AD Add new GEISA and UV cross-section molecules
C     01OCT13 AD Add new GEISA molecules. Change number 135-136
C     15AUG13 AD Add old alternative names for some CFCs etc
C     07AUG13 AD Original.
C
C DESCRIPTION
C     Give molecule name for HITRAN/RFM index, or vice-versa
C     If IDX <= 0 returns IDX for given MOL, or unchanged if MOL unrecognised 
C     If IDX > 0 returns MOL for given IDX, or ' ' if IDX is unrecognised
C
      IMPLICIT NONE
C
      EXTERNAL
     &  CHKSRC ! Set/Get HITRAN/RFM Index for molec. which have two options.
C
C ARGUMENTS
      INTEGER       IDX    ! I/O Molecule index
      CHARACTER*(*) MOL    ! I/O Molecule name
      LOGICAL       FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG !  O  Error message written if FAIL is TRUE
C
C GLOBAL CONSTANTS
      INCLUDE 'idxcon.inc' ! HITRAN molecular indices
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C HITRAN line molecules
      IF ( IDX .GT. 0 ) THEN
        MOL = ' '
        IF ( IDX .EQ. IDXH2O )   MOL = 'h2o' 
        IF ( IDX .EQ. IDXCO2 )   MOL = 'co2' 
        IF ( IDX .EQ. 3 ) THEN
          MOL = 'o3'
          CALL CHKSRC ( MOL, IDX, .FALSE., FAIL, ERRMSG )
        END IF
        IF ( IDX .EQ. 4 ) THEN
          MOL = 'n2o'
          CALL CHKSRC ( MOL, IDX, .FALSE., FAIL, ERRMSG )
        END IF
        IF ( IDX .EQ. 5 )   MOL = 'co' 
        IF ( IDX .EQ. 6 )   MOL = 'ch4'
        IF ( IDX .EQ. IDXO2 )   MOL = 'o2' 
        IF ( IDX .EQ. 8 )   MOL = 'no' 
        IF ( IDX .EQ. 9 ) THEN
          MOL = 'so2'
          CALL CHKSRC ( MOL, IDX, .FALSE., FAIL, ERRMSG )
        END IF
        IF ( IDX .EQ. 10 ) THEN
          MOL = 'no2'
          CALL CHKSRC ( MOL, IDX, .FALSE., FAIL, ERRMSG )
        END IF
        IF ( IDX .EQ. 11 )  MOL = 'nh3'
        IF ( IDX .EQ. 12 )  MOL = 'hno3'
        IF ( IDX .EQ. 13 )  MOL = 'oh' 
        IF ( IDX .EQ. 14 )  MOL = 'hf' 
        IF ( IDX .EQ. 15 )  MOL = 'hcl'
        IF ( IDX .EQ. 16 )  MOL = 'hbr'
        IF ( IDX .EQ. 17 )  MOL = 'hi' 
        IF ( IDX .EQ. 18 )  MOL = 'clo'
        IF ( IDX .EQ. 19 )  MOL = 'ocs'
        IF ( IDX .EQ. 20 ) THEN
          MOL = 'h2co'
          CALL CHKSRC ( MOL, IDX, .FALSE., FAIL, ERRMSG )
        END IF
        IF ( IDX .EQ. 21 )  MOL = 'hocl'
        IF ( IDX .EQ. IDXN2 )  MOL = 'n2' 
        IF ( IDX .EQ. 23 )  MOL = 'hcn'
        IF ( IDX .EQ. 24 )  MOL = 'ch3cl'
        IF ( IDX .EQ. 25 )  MOL = 'h2o2'
        IF ( IDX .EQ. 26 )  MOL = 'c2h2'
        IF ( IDX .EQ. 27 ) THEN
           MOL = 'c2h6'
           CALL CHKSRC ( MOL, IDX, .FALSE., FAIL, ERRMSG )
        END IF
        IF ( IDX .EQ. 28 )  MOL = 'ph3' 
        IF ( IDX .EQ. 29 )  MOL = 'cof2'
        IF ( IDX .EQ. 30 ) THEN
          MOL = 'sf6'
          CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
        IF ( IDX .EQ. 31 )  MOL = 'h2s' 
        IF ( IDX .EQ. 32 )  MOL = 'hcooh'
        IF ( IDX .EQ. 33 )  MOL = 'ho2'
        IF ( IDX .EQ. 34 )  MOL = 'o' 
        IF ( IDX .EQ. 35 ) THEN
          MOL = 'clono2'
          CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
        IF ( IDX .EQ. 36 )  MOL = 'no+'
        IF ( IDX .EQ. 37 )  MOL = 'hobr'
        IF ( IDX .EQ. 38 ) THEN
          MOL = 'c2h4'
          CALL CHKSRC ( MOL, IDX, .FALSE., FAIL, ERRMSG )
        END IF
        IF ( IDX .EQ. 39 ) THEN
          MOL = 'ch3oh'
          CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
        IF ( IDX .EQ. 40 )  MOL = 'ch3br' 
        IF ( IDX .EQ. 41 ) THEN
          MOL = 'ch3cn'
          CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
        IF ( IDX .EQ. 42 ) THEN
           MOL = 'cf4'
          CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
        IF ( IDX .EQ. 43 )  MOL = 'c4h2'
        IF ( IDX .EQ. 44 )  MOL = 'hc3n'
        IF ( IDX .EQ. 45 )  MOL = 'h2' 
        IF ( IDX .EQ. 46 )  MOL = 'cs' 
        IF ( IDX .EQ. 47 )  MOL = 'so3'
C
        IF ( IDX .EQ. 50 ) THEN
           MOL = 'bro'
           CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
C
C Additional GEISA molecules
        IF ( IDX .EQ. 51 )  MOL = 'geh4'
        IF ( IDX .EQ. 52 ) THEN
           MOL = 'c3h8'
           CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
        IF ( IDX .EQ. 53 )  MOL = 'c2n2'
        IF ( IDX .EQ. 54 )  MOL = 'c3h4'
        IF ( IDX .EQ. 55 )  MOL = 'hnc'
        IF ( IDX .EQ. 56 ) THEN
           MOL = 'c6h6'
           CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
C
C Heavy  molecules
        IF ( IDX .EQ. IDXAER ) MOL = 'aerosol'
        IF ( IDX .EQ. 101 ) THEN
          MOL = 'clono2'                    ! Chlorine Nitrate
          CALL CHKSRC ( MOL, IDX, .FALSE., FAIL, ERRMSG )
        END IF
        IF ( IDX .EQ. 102 ) MOL = 'n2o5'    ! DiNitrogen Pentoxide
        IF ( IDX .EQ. 103 ) THEN
          MOL = 'sf6'     ! Sulphur Hexafluoride
          CALL CHKSRC ( MOL, IDX, .FALSE., FAIL, ERRMSG )
        END IF
        IF ( IDX .EQ. 104 ) MOL = 'ccl4'    ! Carbon Tetrachloride
        IF ( IDX .EQ. 105 ) MOL = 'hno4'    ! Peroxynitric Acid
        IF ( IDX .EQ. 106 ) MOL = 'sf5cf3'  !Trifluoromethyl Sulphur Pentaflouride 
        IF ( IDX .EQ. 107 ) MOL = 'brono2'  ! Bromine Nitrate
        IF ( IDX .EQ. 108 ) MOL = 'cloocl'  ! Chlorine Peroxide
        IF ( IDX .EQ. 109 ) MOL = 'x109'    ! spare
        IF ( IDX .EQ. 110 ) MOL = 'x110'    ! spare
C
C CFCs
        IF ( IDX .EQ. 111 ) MOL = 'f11'     ! ccl3f
        IF ( IDX .EQ. 112 ) MOL = 'f12'     ! ccl2f2
        IF ( IDX .EQ. 113 ) MOL = 'f113'    ! c2cl3f3
        IF ( IDX .EQ. 114 ) MOL = 'f114'    ! c2cl2f4
        IF ( IDX .EQ. 115 ) MOL = 'f115'    ! c2clf5
        IF ( IDX .EQ. 116 ) MOL = 'f116'    ! c2f6
        IF ( IDX .EQ. 117 ) MOL = 'f13'     ! cclf3
        IF ( IDX .EQ. 118 ) THEN
          MOL = 'f14'                       ! cf4
          CALL CHKSRC ( MOL, IDX, .FALSE., FAIL, ERRMSG )
        END IF
C
C HCFCs
        IF ( IDX .EQ. 121 ) MOL = 'f21'     ! chcl2f
        IF ( IDX .EQ. 122 ) MOL = 'f22'     ! chclf2
        IF ( IDX .EQ. 123 ) MOL = 'f123'    ! chcl2cf3
        IF ( IDX .EQ. 124 ) MOL = 'f124'    ! chclfcf3
        IF ( IDX .EQ. 125 ) MOL = 'f141b'   ! ch3ccl2f
        IF ( IDX .EQ. 126 ) MOL = 'f142b '  ! ch3cclf2
        IF ( IDX .EQ. 127 ) MOL = 'f225ca'  ! chcl2cf2cf3
        IF ( IDX .EQ. 128 ) MOL = 'f225cb'  ! cclf2cf2chclf
C
C HFCs
        IF ( IDX .EQ. 131 ) MOL = 'f125'    ! chf2cf3
        IF ( IDX .EQ. 132 ) MOL = 'f32'     ! ch2f2
        IF ( IDX .EQ. 133 ) MOL = 'f134'    ! chf2chf2
        IF ( IDX .EQ. 134 ) MOL = 'f134a'   ! cfh2cf3
C 01OCT13: Insert 135, f143, from GEISA and renumber f143a,f152a
        IF ( IDX .EQ. 135 ) MOL = 'f143'    ! ch2fchf2 1,1,2-Trifluoroethane
        IF ( IDX .EQ. 136 ) MOL = 'f143a'   ! cf3ch3   1,1,1-Trifluoroethane
        IF ( IDX .EQ. 137 ) MOL = 'f152a'   ! ch3chf2
C New GEISA molecule
        IF ( IDX .EQ. 138 ) MOL = 'f365mfc' ! CF3CH2CF2CH3, pentaflourobutane
C
C MeHCs
        IF ( IDX .EQ. 141 ) THEN
          MOL = 'ch3oh'                     ! Methanol
          CALL CHKSRC ( MOL, IDX, .FALSE., FAIL, ERRMSG )
        END IF
        IF ( IDX .EQ. 142 ) THEN
          MOL = 'ch3cn'                     ! Acetonitrile
          CALL CHKSRC ( MOL, IDX, .FALSE., FAIL, ERRMSG )
        END IF
        IF ( IDX .EQ. 143 ) MOL = 'ch3cho'  ! Acetaldehyde
        IF ( IDX .EQ. 144 ) MOL = 'acetone' ! ch3coch3
        IF ( IDX .EQ. 145 ) MOL = 'pan'     ! ch3c(o)oono2
C
C NMCs
        IF ( IDX .EQ. 151 ) THEN
          MOL = 'c2h6'                      ! Ethane
          CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
        IF ( IDX .EQ. 152 ) MOL = 'c3h8'    ! Propane
        IF ( IDX .EQ. 153 ) THEN
          MOL = 'c6h6'                      ! Benzene
          CALL CHKSRC ( MOL, IDX, .FALSE., FAIL, ERRMSG )
        END IF
C New GEISA molecules
        IF ( IDX .EQ. 154 ) THEN
          MOL = 'c2h2'                      ! Acetylene
          CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
        IF ( IDX .EQ. 155 ) THEN
          MOL = 'c2h4'                      ! Ethylene
          CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
C
C Halocarbons
        IF ( IDX .EQ. 161 ) MOL = 'c4f8'    ! Octafluorocyclobutane
C
C UV Cross-sections
        IF ( IDX .EQ. 171 ) THEN
          MOL = 'o3'
          CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
        IF ( IDX .EQ. 172 ) THEN
          MOL = 'n2o'
          CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
        IF ( IDX .EQ. 173 ) THEN
          MOL = 'so2'
          CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
        IF ( IDX .EQ. 174 ) THEN
          MOL = 'no2'
          CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
        IF ( IDX .EQ. 175 ) THEN
          MOL = 'h2co'
          CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
        IF ( IDX .EQ. 176 ) THEN
          MOL = 'bro'
          CALL CHKSRC ( MOL, IDX, .FALSE., FAIL, ERRMSG )
        END IF
        IF ( IDX .EQ. 177 ) MOL = 'no3'
        IF ( IDX .EQ. 178 ) MOL = 'oclo'
        IF ( IDX .EQ. 181 ) MOL = 'c7h8'
        IF ( IDX .EQ. 182 ) MOL = 'oxylene'               ! o-C8H10
        IF ( IDX .EQ. 183 ) MOL = 'mxylene'               ! m-C8H10
        IF ( IDX .EQ. 184 ) MOL = 'pxylene'               ! p-C8H10
C
C Return index for given molecule name
      ELSE
C
        IF ( MOL .EQ. 'h2o'     ) IDX = IDXH2O
        IF ( MOL .EQ. 'co2'     ) IDX = IDXCO2
        IF ( MOL .EQ. 'o3'      ) THEN
          IDX = 3
          CALL CHKSRC ( MOL, IDX, .FALSE., FAIL, ERRMSG )
        END IF
        IF ( MOL .EQ. 'n2o'     ) THEN
          IDX = 4
          CALL CHKSRC ( MOL, IDX, .FALSE., FAIL, ERRMSG )
        END IF
        IF ( MOL .EQ. 'co'      ) IDX = 5
        IF ( MOL .EQ. 'ch4'     ) IDX = 6
        IF ( MOL .EQ. 'o2'      ) IDX = IDXO2
        IF ( MOL .EQ. 'no'      ) IDX = 8
        IF ( MOL .EQ. 'so2'     ) THEN
          IDX = 9
          CALL CHKSRC ( MOL, IDX, .FALSE., FAIL, ERRMSG )
        END IF
        IF ( MOL .EQ. 'no2'     ) THEN
          IDX = 10
          CALL CHKSRC ( MOL, IDX, .FALSE., FAIL, ERRMSG )
        END IF
        IF ( MOL .EQ. 'nh3'     ) IDX = 11
        IF ( MOL .EQ. 'hno3'    ) IDX = 12
        IF ( MOL .EQ. 'oh'      ) IDX = 13
        IF ( MOL .EQ. 'hf'      ) IDX = 14
        IF ( MOL .EQ. 'hcl'     ) IDX = 15
        IF ( MOL .EQ. 'hbr'     ) IDX = 16
        IF ( MOL .EQ. 'hi'      ) IDX = 17
        IF ( MOL .EQ. 'clo'     ) IDX = 18
        IF ( MOL .EQ. 'ocs'     ) IDX = 19
        IF ( MOL .EQ. 'h2co'    ) THEN
          IDX = 20
          CALL CHKSRC ( MOL, IDX, .FALSE., FAIL, ERRMSG )
        END IF
        IF ( MOL .EQ. 'hocl'    ) IDX = 21
        IF ( MOL .EQ. 'n2'      ) IDX = IDXN2
        IF ( MOL .EQ. 'hcn'     ) IDX = 23
        IF ( MOL .EQ. 'ch3cl'   ) IDX = 24
        IF ( MOL .EQ. 'h2o2'    ) IDX = 25
        IF ( MOL .EQ. 'c2h2'    ) IDX = 26
        IF ( MOL .EQ. 'c2h6'    ) THEN
          IDX = 27
          CALL CHKSRC ( MOL, IDX, .FALSE., FAIL, ERRMSG )
        END IF
        IF ( MOL .EQ. 'ph3'     ) IDX = 28
        IF ( MOL .EQ. 'cof2'    ) IDX = 29
        IF ( MOL .EQ. 'sf6q'    ) THEN
          IDX = 30
          MOL = 'sf6'
          CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
        IF ( MOL .EQ. 'h2s'     ) IDX = 31
        IF ( MOL .EQ. 'hcooh'   ) IDX = 32
        IF ( MOL .EQ. 'ho2'     ) IDX = 33
        IF ( MOL .EQ. 'o'       ) IDX = 34
        IF ( MOL .EQ. 'clono2q' ) THEN
          IDX = 35
          MOL = 'clono2'
          CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
        IF ( MOL .EQ. 'no+'     ) IDX = 36
        IF ( MOL .EQ. 'hobr'    ) IDX = 37
        IF ( MOL .EQ. 'c2h4'    ) THEN
          IDX = 38
          CALL CHKSRC ( MOL, IDX, .FALSE., FAIL, ERRMSG )
        END IF
        IF ( MOL .EQ. 'ch3ohq'  ) THEN
          IDX = 39
          MOL = 'ch3oh'
          CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
        IF ( MOL .EQ. 'ch3br'   ) IDX = 40
        IF ( MOL .EQ. 'ch3cnq'  ) THEN
          IDX = 41
          MOL = 'ch3cn'
          CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF          
        IF ( MOL .EQ. 'cf4q' .OR. MOL .EQ. 'f14q' ) THEN
          IDX = 42
          MOL = 'f14'
          CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF          
        IF ( MOL .EQ. 'c4h2'    ) IDX = 43
        IF ( MOL .EQ. 'hc3n'    ) IDX = 44
        IF ( MOL .EQ. 'h2'      ) IDX = 45
        IF ( MOL .EQ. 'cs'      ) IDX = 46
        IF ( MOL .EQ. 'so3'     ) IDX = 47
C
        IF ( MOL .EQ. 'broq' ) THEN
          IDX = 50
          MOL = 'bro'
          CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
C
C Additional GEISA line molecules
        IF ( MOL .EQ. 'geh4'    ) IDX = 51
        IF ( MOL .EQ. 'c3h8q' ) THEN
          IDX = 52
          MOL ='c3h8'
          CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
        IF ( MOL .EQ. 'c2n2'    ) IDX = 53
        IF ( MOL .EQ. 'c3h4'    ) IDX = 54
        IF ( MOL .EQ. 'hnc'     ) IDX = 55
        IF ( MOL .EQ. 'c6h6q' ) THEN
          IDX = 56
          MOL ='c6h6'
          CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
C
C Heavy molecules
        IF ( MOL .EQ. 'aerosol' ) IDX = IDXAER
        IF ( MOL .EQ. 'clono2' ) THEN          ! Chlorine nitrate
          IDX = 101
          CALL CHKSRC ( MOL, IDX, .FALSE., FAIL, ERRMSG )
        END IF
        IF ( MOL .EQ. 'n2o5'   ) IDX = 102     ! DiNitrogen Pentoxide
        IF ( MOL .EQ. 'sf6'    ) THEN          ! Sulphur Hexafluoride
          IDX = 103    
          CALL CHKSRC ( MOL, IDX, .FALSE., FAIL, ERRMSG )
        END IF
        IF ( MOL .EQ. 'ccl4'   ) IDX = 104     ! Carbon Tetrachloride
        IF ( MOL .EQ. 'hno4'   ) IDX = 105     ! Peroxynitric Acid
        IF ( MOL .EQ. 'sf5cf3' ) IDX = 106     ! Trifluoromethyl Sulphur Pentaflouride 
        IF ( MOL .EQ. 'brono2' ) IDX = 107     ! Bromine Nitrate
        IF ( MOL .EQ. 'cloocl' ) IDX = 108     ! Chlorine Peroxide
        IF ( MOL .EQ. 'x109'   ) IDX = 109     ! spare
        IF ( MOL .EQ. 'x110'   ) IDX = 110     ! spare
C
C CFCs
        IF ( MOL .EQ. 'ccl3f'   .OR. 
     &       MOL .EQ. 'cfcl3'   .OR. MOL .EQ. 'f11'  ) IDX = 111
        IF ( MOL .EQ. 'ccl2f2'  .OR. 
     &       MOL .EQ. 'cf2cl2'  .OR. MOL .EQ. 'f12'  ) IDX = 112
        IF ( MOL .EQ. 'c2cl3f3' .OR. MOL .EQ. 'f113' ) IDX = 113
        IF ( MOL .EQ. 'c2cl2f4' .OR. MOL .EQ. 'f114' ) IDX = 114
        IF ( MOL .EQ. 'c2clf5'  .OR. MOL .EQ. 'f115' ) IDX = 115
        IF ( MOL .EQ. 'c2f6'    .OR. MOL .EQ. 'f116' ) IDX = 116
        IF ( MOL .EQ. 'cclf3'   .OR. MOL .EQ. 'f13'  ) IDX = 117
        IF ( MOL .EQ. 'cf4'     .OR. MOL .EQ. 'f14'  ) THEN
          IDX = 118
          MOL = 'f14'
          CALL CHKSRC ( MOL, IDX, .FALSE., FAIL, ERRMSG )
        END IF
C
C HCFCs
        IF ( MOL .EQ. 'chcl2f'   .OR. MOL .EQ. 'f21'   ) IDX = 121
        IF ( MOL .EQ. 'chclf2'   .OR. MOL .EQ. 'f22'   ) IDX = 122
        IF ( MOL .EQ. 'chcl2cf3' .OR. MOL .EQ. 'f123'  ) IDX = 123
        IF ( MOL .EQ. 'chclfcf3' .OR. MOL .EQ. 'f124'  ) IDX = 124
        IF ( MOL .EQ. 'ch3ccl2f' .OR. MOL .EQ. 'f141b' ) IDX = 125
        IF ( MOL .EQ. 'ch3cclf2' .OR. MOL .EQ. 'f142b' ) IDX = 126

        IF ( MOL .EQ. 'chcl2cf2cf3' .OR. MOL .EQ. 'f225ca' ) IDX = 127
        IF ( MOL .EQ. 'cclf2cf2chclf' .OR. MOL .EQ. 'f225cb' ) IDX=128
C
C HFCs
        IF ( MOL .EQ. 'chf2cf3'  .OR. MOL .EQ. 'f125'  ) IDX = 131
        IF ( MOL .EQ. 'ch2f2'    .OR. MOL .EQ. 'f32'   ) IDX = 132
        IF ( MOL .EQ. 'chf2chf2' .OR. MOL .EQ. 'f134'  ) IDX = 133
        IF ( MOL .EQ. 'cfh2cf3'  .OR. MOL .EQ. 'f134a' ) IDX = 134
        IF ( MOL .EQ.                          'f143'  ) IDX = 135
        IF ( MOL .EQ. 'cf3ch3'   .OR. MOL .EQ. 'f143a' ) IDX = 136
        IF ( MOL .EQ. 'ch3chf2'  .OR. MOL .EQ. 'f152a' ) IDX = 137
        IF ( MOL .EQ.                        'f365mfc' ) IDX = 138
C
C MeHCs
        IF ( MOL .EQ. 'ch3oh'    ) THEN          ! Methanol
          IDX = 141
          CALL CHKSRC ( MOL, IDX, .FALSE., FAIL, ERRMSG )
        END IF
        IF ( MOL .EQ. 'ch3cn'    ) THEN          ! Acetonitrile
          IDX = 142
          CALL CHKSRC ( MOL, IDX, .FALSE., FAIL, ERRMSG )
        END IF
        IF ( MOL .EQ. 'ch3cho'   ) IDX = 143     ! Acetaldehyde
        IF ( MOL .EQ. 'ch3coch3' .OR. MOL .EQ. 'acetone' ) IDX = 144  ! Acetone
        IF ( MOL .EQ. 'ch3c(o)oono2' .OR. MOL .EQ. 'pan' ) IDX = 145  ! PAN
C
C NMCs
        IF ( MOL .EQ. 'c2h6x'  ) THEN            ! Ethane
          IDX = 151
          MOL = 'c2h6'
          CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
        IF ( MOL .EQ. 'c3h8'  ) THEN             ! Propane
          IDX = 152
          CALL CHKSRC ( MOL, IDX, .FALSE., FAIL, ERRMSG )
        END IF
        IF ( MOL .EQ. 'c2h2x'  ) THEN             ! Acetylene
          IDX = 154
          MOL = 'c2h2'
          CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
        IF ( MOL .EQ. 'c2h4x'  ) THEN             ! Ethylene
          IDX = 155                      
          MOL = 'c2h4'
          CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
C
C Halocarbons
        IF ( MOL .EQ. 'c4f8' ) IDX = 161         ! Octafluorocyclobutane
C
C UV Cross-sections
        IF ( MOL .EQ. 'o3x'  ) THEN
          IDX = 171                      
          MOL = 'o3'
          CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
        IF ( MOL .EQ. 'n2ox'  ) THEN
          IDX = 172                      
          MOL = 'n2o'
          CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
        IF ( MOL .EQ. 'so2x'  ) THEN
          IDX = 173                      
          MOL = 'so2'
          CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
        IF ( MOL .EQ. 'no2x'  ) THEN
          IDX = 174                      
          MOL = 'no2'
          CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
        IF ( MOL .EQ. 'h2cox'  ) THEN
          IDX = 175                      
          MOL = 'h2co'
          CALL CHKSRC ( MOL, IDX, .TRUE., FAIL, ERRMSG )
        END IF
        IF ( MOL .EQ. 'bro'  ) THEN
          IDX = 176
          CALL CHKSRC ( MOL, IDX, .FALSE., FAIL, ERRMSG )
        END IF
        IF ( MOL .EQ. 'no3'  ) IDX = 177
        IF ( MOL .EQ. 'oclo' ) IDX = 178
        IF ( MOL .EQ. 'c7h8' ) IDX = 181           ! Toluene
        IF ( MOL .EQ. 'oxylene' ) IDX = 182        ! o-C8H10
        IF ( MOL .EQ. 'mxylene' ) IDX = 183        ! m-C8H10
        IF ( MOL .EQ. 'pxylene' ) IDX = 184        ! p-C8H10

      END IF
C
      END

