      SUBROUTINE GASISO ( FAIL, ERRMSG )
C
C VERSION
C     14JAN14 AD Minor Bug: Increase NISHIT for CO2 from 10 to 11.
C                    Reduce NISHIT from 9 to 4 for SO2
C     10-JAN-14  AD  New GEISA molecules & BrO
C                    Remove TCOEFF and WIDSTP - now local in AWIDTH
C                    Check that non-zero weight is assigned
C     04-AUG-13  AD  New HITRAN 2012 molecules and isotopes
C                    Replace TCODEF, WIDDEF with explicit values in DATA
C     09-JAN-09  AD  Replace HDO with CH3OH (ID=39)
C     03-JAN-08  AD  Remove ISOMOL, IDIGAS,ISOGAS(0,*) - now set in GASCHK
C     21-MAR-04  AD  Initialise ISOMOL
C     06-FEB-04  AD  Add IDIGAS
C     30-OCT-03  AD  Add SO2 isotopes #3 (628) and #4 (636)
C     23-FEB-03  AD  Add BrO as ID=40, move GeH4 (was ID=40) to new ID=46
C     12-JUN-02  AD  Add HDO as a separate gas
C     15-FEB-02  AD  Set ISOGAS
C     27-JUL-01  AD  Add Isotope#6 for H2O (172) 
C     26-OCT-00  AD  Remove local dimension MAXHIS - use MAXISO instead.
C     16-OCT-00  AD  Add GEISA isotopes
C     03-MAR-97  AD  Version 3.
C     23-NOV-96  AD  Add extra 5 HITRAN'96 species.
C     01-OCT-96  AD  Version 2.
C     01-OCT-96  AD  Add extra CO isotope 37, wgt=30
C     21-SEP-96  AD  Assume wgts of O3 isotopes 4,5 = 49 a.u (not 52).
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Load isotopic data for required gases.
C     Called by GASCHK for each gas required with line data (ie not X/S data)
C     Assumes NGAS currently set to point to required gas in *GAS variables
C
      IMPLICIT NONE
C
C ARGUMENTS
      LOGICAL      FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG !  O  Error message written if FAIL is TRUE
C
C GLOBAL CONSTANTS
      INCLUDE 'idxcon.inc' ! HITRAN molecular indices
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'gascom.inc' ! Molecular and isotope data
C
C LOCAL VARIABLES
      INTEGER IDX                    ! HITRAN index of species
      INTEGER IISO                   ! Counter for HITRAN isotope lists
      INTEGER NISHIT(MAXHLN)         ! No.isotopes for each HITRAN line species
      REAL    WGTISO(MAXISO,MAXHLN)  ! Isotopic weights [Atomic units]
C
C DATA STATEMENTS
      DATA NISHIT / 6,11, 5, 5, 6, 4, 3, 3, 4, 1, 2, 2, 3, 2, 4, 4,  !  1-16
     &              2, 2, 6, 3, 2, 2, 3, 2, 1, 3, 2, 1, 2, 1, 3, 1,  ! 17-32
     &              1, 1, 2, 1, 2, 2, 1, 2, 1, 1, 1, 1, 2, 4, 1, 0,  ! 33-48
     &              0, 1, 1, 1, 1, 1, 1, 1 /                         ! 49-56
C H2O                                    161  181  171  162  182  172
      DATA (WGTISO(IISO,1),IISO=1, 6) / 18., 20., 19., 19., 21., 20./
      DATA (WGTISO(IISO,2),IISO=1, 11) 
C CO2       626  636  628  627  638  637  828  827  727  838  837
     &     / 44., 45., 46., 45., 47., 46., 48., 47., 46., 49., 48. / 
C O3                                   666  668  686  667  676
      DATA (WGTISO(IISO,3),IISO=1, 5) / 48., 50., 50., 49., 49./
C N2O                                  446  456  546  448  447
      DATA (WGTISO(IISO,4),IISO=1, 5) / 44., 45., 45., 46., 45./
C CO                                    26   36   28   27   38   37
      DATA (WGTISO(IISO,5),IISO=1, 6) / 28., 29., 30., 29., 31., 30./
C CH4                                  211  311  212  312
      DATA (WGTISO(IISO,6),IISO=1, 4) / 16., 17., 17., 18. /
C O2                                    66   68   67
      DATA (WGTISO(IISO,7),IISO=1, 3) / 32., 34., 33./
C NO                                    46   56   48
      DATA (WGTISO(IISO,8),IISO=1, 3) / 30., 31., 32./
C SO2                                  626  646  628* 636*    *Geisa
      DATA (WGTISO(IISO,9),IISO=1, 4) / 64., 66., 66., 65./
C NO2                                  646
      DATA (WGTISO(IISO,10),IISO=1,1) / 46./
C NH3                                 4111 5111
      DATA (WGTISO(IISO,11),IISO=1,2) / 17., 18./
C HNO3                                 146  156
      DATA (WGTISO(IISO,12),IISO=1,2) / 63., 64./
C OH                                    61   81   62
      DATA (WGTISO(IISO,13),IISO=1,3) / 17., 19., 18./
C HF                                    19   29
      DATA (WGTISO(IISO,14),IISO=1,2) / 20., 21./
C HCl                                   15   17   25   27
      DATA (WGTISO(IISO,15),IISO=1,4) / 36., 38., 37., 39./
C HBr                                   19   11   29   21
      DATA (WGTISO(IISO,16),IISO=1,4) / 80., 82., 81., 83./
C HI                                    17    27
      DATA (WGTISO(IISO,17),IISO=1,2) /128., 129./
C ClO                                   56   76
      DATA (WGTISO(IISO,18),IISO=1,2) / 51., 53./
C OCS                                  622  624  632  623  822  634* *Geisa
      DATA (WGTISO(IISO,19),IISO=1,6) / 60., 62., 61., 61., 62., 63. /
C H2CO                                 126  136  128
      DATA (WGTISO(IISO,20),IISO=1,3) / 30., 31., 32./
C HOCl                                 165  167
      DATA (WGTISO(IISO,21),IISO=1,2) / 52., 54./
C N2                                    44   45
      DATA (WGTISO(IISO,22),IISO=1,2) / 28., 29./
C HCN                                  124  134  125
      DATA (WGTISO(IISO,23),IISO=1,3) / 27., 28., 28./
C CH3Cl                                215  217
      DATA (WGTISO(IISO,24),IISO=1,2) / 50., 52./
C H2O2                                1661
      DATA (WGTISO(IISO,25),IISO=1,1) / 34./
C C2H2                                1221 1231 1222
      DATA (WGTISO(IISO,26),IISO=1,3) / 26., 27., 27./
C C2H6                                1221 1231*
      DATA (WGTISO(IISO,27),IISO=1,2) / 30., 31./
C PH3                                 1111
      DATA (WGTISO(IISO,28),IISO=1,1) / 34./
C COF2                                 269  369
      DATA (WGTISO(IISO,29),IISO=1,2) / 66., 67./
C SF6                                    29
      DATA (WGTISO(IISO,30),IISO=1,1) / 146./
C H2S                                  121  141  131
      DATA (WGTISO(IISO,31),IISO=1,3) / 34., 36., 35. /    
C HCOOH                                126
      DATA (WGTISO(IISO,32),IISO=1,1) / 46./               
C HO2                                  166
      DATA (WGTISO(IISO,33),IISO=1,1) / 33./               
C O                                      6
      DATA (WGTISO(IISO,34),IISO=1,1) / 16./               
C ClONO2                              5646 7646
      DATA (WGTISO(IISO,35),IISO=1,2) / 97., 99./          
C NO+                                   46
      DATA (WGTISO(IISO,36),IISO=1,1) / 30./               
C HOBr                                 169  161
      DATA (WGTISO(IISO,37),IISO=1,2) / 96., 98. /               
C C2H4                                 221  231  
      DATA (WGTISO(IISO,38),IISO=1,2) / 28., 29./               
C CH3OH                                2161
      DATA (WGTISO(IISO,39),IISO=1,1) / 32. /
C CH3Br                                219  211 
      DATA (WGTISO(IISO,40),IISO=1,2) / 94., 96. /
C CH3CN                                2124
      DATA (WGTISO(IISO,41),IISO=1,1) / 41. /               
C CF4                                   29
      DATA (WGTISO(IISO,42),IISO=1,1) / 88. /               
C C4H2                                 2211
      DATA (WGTISO(IISO,43),IISO=1,1) / 50. /               
C HC3N                                 1224
      DATA (WGTISO(IISO,44),IISO=1,1) / 51. /               
C H2                                   11   12 
      DATA (WGTISO(IISO,45),IISO=1,2) / 2., 3. /               
C CS                                    22   24   32   23 
      DATA (WGTISO(IISO,46),IISO=1,4) / 44., 46., 45., 45. /               
C SO3                                   26
      DATA (WGTISO(IISO,47),IISO=1,1) / 80. /               
C
C BrO                                   56
      DATA (WGTISO(IISO,50),IISO=1,1) / 96. /               

C GEISA molecules
C GeH4                                  411*
      DATA (WGTISO(IISO,51),IISO=1,1) / 77. /           
C C3C8                                  221*
      DATA (WGTISO(IISO,52),IISO=1,1) / 44. /               
C C2N2                                  224*
      DATA (WGTISO(IISO,53),IISO=1,1) / 52. /               
C C3H4                                  341*
      DATA (WGTISO(IISO,54),IISO=1,1) / 40. /               
C HNC                                   142*
      DATA (WGTISO(IISO,55),IISO=1,1) / 27. /               
C C2H6                                  2211
      DATA (WGTISO(IISO,56),IISO=1,1) / 30. /               
C
        SAVE NISHIT, WGTISO
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      IDX = IDXGAS(NGAS)
C
C Load isotope data
      IF ( IDX .LE. MAXHLN ) THEN         ! only exists for line data
        DO IISO = 1, NISHIT(IDX) 
          IF ( WGTISO(IISO,IDX) .EQ. 0.0 ) THEN
            FAIL = .TRUE.
            WRITE ( ERRMSG, '(A,A,A,I2)' )
     &        'F-GASISO: No weight defined for gas,iso ', CODGAS(NGAS),
     &        ',', IISO 
            RETURN
          END IF
          WGTGAS(IISO,NGAS) = WGTISO(IISO,IDX)
          ISOGAS(IISO,NGAS) = NGAS
        END DO
        NISGAS(NGAS) = NISHIT(IDX)
        FAIL = .FALSE.
      ELSE 
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,I11)' ) 
     &    'F-GASISO: Called with unexpected Gas ID,=', IDX
      END IF
C
      END
