      SUBROUTINE QTFCT_V42 ( IDGAS, ISO, TEM, SQ )
C
C VERSION
C     23-FEB-03  AD  Add BrO as ID=40, move GeH4 (was ID=40) to new ID=46
C     13-FEB-03  AD  Updated to use TIPS97.FOR coefficients plus 
C                    extra values supplied by Bianca Dinelli (BMD).
C                    Correction to handling of minor isotopes.
C     23-JAN-01  AD  Local version, new HNO3 coefficients
C     16-OCT-00  AD  Check for IDGAS >38 as well
C     13-SEP-99  AD  Add extra isotope for C2H4 (ID=38)
C     09-AUG-99  AD  Remove surplus brackets around QCOEFH(ISPC,4)
C     08-JUL-98  AD  Add new HNO3 coefficients (commented out at the moment)
C     03-NOV-97  AD  Remove duplicate NO high-temp data (QCOEFH)
C     25-AUG-97  AD  Add extra data QCOEFH to allow for T=500-1500K 
C     03-MAR-97  AD  Version 3.
C     12-DEC-96  AD  Also check for C2H6 (no TIPS coefficients).
C     23-NOV-96  AD  Revised for HITRAN'96 (except for HNO3 - keep HITRAN'92).
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Calculate total internal partition sums.
C     Called by ADJUST and NTECLC.
C     The ratio of total internal partition sums, Q(296K)/Q(Path Temp) is 
C     calculated for a given molecule and isotpic species. The line strength 
C     can then be adjusted by multiplying by this factor.
C
C     Based on programs TIPS and QSUMS by R.R.GAMACHE which calculate the 
C     total internal partition functions in the HITRAN select routine.
C     Reduced QCOEF to two temperature ranges (70-500K, 500-1500K).
C
C     Additional values supplied taken from "MARSCHALS Partition Sums" 
C     B.M.Dinelli, January 2003.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER IDGAS !  I  HITRAN gas ID
      INTEGER ISO   !  I  HITRAN isotope ID 
      REAL    TEM   !  I  Path temperature [K]
      REAL    SQ    !  O  Ratio of tot. part. sum at 296K to tps at pth temp
C
      EXTERNAL
     &  TPSLKP ! Look up tabulated TIPS factor
      REAL TPSLKP
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM array sizes
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
C
C LOCAL CONSTANTS
      INTEGER MAXH42              ! v42 version of MAXHLN in rfmsiz.inc
        PARAMETER ( MAXH42 = 46 ) 
      INTEGER MAXSPE              ! Total number of isotope species (all gases)
        PARAMETER ( MAXSPE = 104 )
      INTEGER MAXQCF              ! No. of TIPS polynomial coefficients
        PARAMETER ( MAXQCF = 4 )
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'tpscom.inc' ! Tabulated TIPS data
C
C LOCAL VARIABLES
      INTEGER          I          ! Counter for DATA statements
      INTEGER          ISPC       ! Counter for isotope species (all gases)
      INTEGER ISPIDG(MAXH42)      ! Starting pos-1 of each gas in the isotope array
      DOUBLE PRECISION DBLTEM     ! Double Precision TEM
      DOUBLE PRECISION QSTD       ! Q296 for selected isotope species
      DOUBLE PRECISION QT         
      DOUBLE PRECISION QCOEF(MAXSPE,MAXQCF) ! TIPS 70-500K
      DOUBLE PRECISION QCOEFH(MAXSPE,MAXQCF) ! TIPS 500-1500K 
      DOUBLE PRECISION Q296(MAXSPE) ! TIPS at reference Temperature, 296 K.

        DATA ( ISPIDG(I), I = 1, MAXH42 ) / 
     &     0, 6,15,20,25,31,34,37,40,42,43,45,46,49,
     &    50,52,54,55,57,63,66,68,69,72,74,75,77,79,
     &    80,81,82,85,86,87,88,90,91,93,95,96,98,99,
     &    100,101,102,103 /
        SAVE ISPIDG
c
c...Total internal partition sums for T>=70 to <=500 K range:
c...   H2O  --   161
      DATA (QCOEF( 1,I),I=1,4)/-.44405D+01, .27678D+00,
     &                .12536D-02,-.48938D-06/
c...   H2O  --   181
      DATA (QCOEF( 2,I),I=1,4)/-.43624D+01, .27647D+00,
     &                .12802D-02,-.52046D-06/
c...   H2O  --   171
      DATA (QCOEF( 3,I),I=1,4)/-.25767D+02, .16458D+01,
     &                .76905D-02,-.31668D-05/
c...   H2O  --   162
      DATA (QCOEF( 4,I),I=1,4)/-.23916D+02, .13793D+01,
     &                .61246D-02,-.21530D-05/
c...   H2O  --   182                       Fitted for 150-300K by B. Dinelli 
      DATA (QCOEF( 5,I),I=1,4)/-.35241D+02, .15025D+01,
     &                .60235D-02,-.23945D-05/
c...   H2O  --   172                       Fitted for 150-300K by B. Dinelli 
      DATA (QCOEF( 6,I),I=1,4)/-.20852D+03, .89400D+01,
     &                .36033D-01,-.14443D-04/
c...   CO2  --   626
      DATA (QCOEF( 7,I),I=1,4)/-.13617D+01, .94899D+00,
     &               -.69259D-03, .25974D-05/
c...   CO2  --   636
      DATA (QCOEF( 8,I),I=1,4)/-.20631D+01, .18873D+01,
     &               -.13669D-02, .54032D-05/
c...   CO2  --   628
      DATA (QCOEF( 9,I),I=1,4)/-.29175D+01, .20114D+01,
     &               -.14786D-02, .55941D-05/
c...   CO2  --   627
      DATA (QCOEF(10,I),I=1,4)/-.16558D+02, .11733D+02,
     &               -.85844D-02, .32379D-04/
c...   CO2  --   638
      DATA (QCOEF(11,I),I=1,4)/-.44685D+01, .40330D+01,
     &               -.29590D-02, .11770D-04/
c...   CO2  --   637
      DATA (QCOEF(12,I),I=1,4)/-.26263D+02, .23350D+02,
     &               -.17032D-01, .67532D-04/
c...   CO2  --   828
      DATA (QCOEF(13,I),I=1,4)/-.14811D+01, .10667D+01,
     &               -.78758D-03, .30133D-05/
c...   CO2  --   728
      DATA (QCOEF(14,I),I=1,4)/-.17600D+02, .12445D+02,
     &               -.91837D-02, .34915D-04/
c...   CO2  --   838                       Geisa isotope, no data avaialable  
      DATA (QCOEF(15,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...    O3  --   666
      DATA (QCOEF(16,I),I=1,4)/-.16443D+03, .69047D+01,
     &                .10396D-01, .26669D-04/
c...    O3  --   668
      DATA (QCOEF(17,I),I=1,4)/-.35222D+03, .14796D+02,
     &                .21475D-01, .59891D-04/
c...    O3  --   686
      DATA (QCOEF(18,I),I=1,4)/-.17466D+03, .72912D+01,
     &                .10093D-01, .29991D-04/
c...    O3  --   667
      DATA (QCOEF(19,I),I=1,4)/-.20540D+04, .85998D+02,
     &                .12667D+00, .33026D-03/
c...    O3  --   676
      DATA (QCOEF(20,I),I=1,4)/-.10148D+04, .42494D+02,
     &                .62586D-01, .16319D-03/
c...   N2O  --   446
      DATA (QCOEF(21,I),I=1,4)/ .24892D+02, .14979D+02,
     &               -.76213D-02, .46310D-04/
c...   N2O  --   456
      DATA (QCOEF(22,I),I=1,4)/ .36318D+02, .95497D+01,
     &               -.23943D-02, .26842D-04/
c...   N2O  --   546
      DATA (QCOEF(23,I),I=1,4)/ .24241D+02, .10179D+02,
     &               -.43002D-02, .30425D-04/
c...   N2O  --   448
      DATA (QCOEF(24,I),I=1,4)/ .67708D+02, .14878D+02,
     &               -.10730D-02, .34254D-04/
c...   N2O  --   447
      DATA (QCOEF(25,I),I=1,4)/ .50069D+03, .84526D+02,
     &                .83494D-02, .17154D-03/
c...    CO  --    26
      DATA (QCOEF(26,I),I=1,4)/ .27758D+00, .36290D+00,
     &               -.74669D-05, .14896D-07/
c...    CO  --    36
      DATA (QCOEF(27,I),I=1,4)/ .53142D+00, .75953D+00,
     &               -.17810D-04, .35160D-07/
c...    CO  --    28
      DATA (QCOEF(28,I),I=1,4)/ .26593D+00, .38126D+00,
     &               -.92083D-05, .18086D-07/
c...    CO  --    27                              
      DATA (QCOEF(29,I),I=1,4)/ .16376D+01, .22343D+01,
     &               -.49025D-04, .97389D-07/
c...    CO  --    38                              
      DATA (QCOEF(30,I),I=1,4)/ .51216D+00, .79978D+00,
     &               -.21784D-04, .42749D-07/
c...    CO  --    37                              
      DATA (QCOEF(31,I),I=1,4)/ .32731D+01, .46577D+01,
     &               -.69833D-04, .18853D-06/
c...   CH4  --   211                              
      DATA (QCOEF(32,I),I=1,4)/-.26479D+02, .11557D+01,
     &                .26831D-02, .15117D-05/
c...   CH4  --   311                              
      DATA (QCOEF(33,I),I=1,4)/-.52956D+02, .23113D+01,
     &                .53659D-02, .30232D-05/
c...   CH4  --   212                              
      DATA (QCOEF(34,I),I=1,4)/-.21577D+03, .93318D+01,
     &                .21779D-01, .12183D-04/
c...    O2  --    66                              
      DATA (QCOEF(35,I),I=1,4)/ .35923D+00, .73534D+00,
     &               -.64870D-04, .13073D-06/
c...    O2  --    68                              
      DATA (QCOEF(36,I),I=1,4)/-.40039D+01, .15595D+01,
     &               -.15357D-03, .30969D-06/
c...    O2  --    67                              
      DATA (QCOEF(37,I),I=1,4)/-.23325D+02, .90981D+01,
     &               -.84435D-03, .17062D-05/
c...    NO  --    46                              TIPS97.FOR modification
      DATA (QCOEF(38,I),I=1,4)/-.25296D+02, .26349D+01,
     &                .58517D-02,-.52020D-05/
c...    NO  --    56                              TIPS97.FOR modification
      DATA (QCOEF(39,I),I=1,4)/-.14990D+02, .18240D+01,
     &                .40261D-02,-.35648D-05/
c...    NO  --    48                              TIPS97.FOR modification
      DATA (QCOEF(40,I),I=1,4)/-.26853D+02, .27816D+01,
     &                .61493D-02,-.54410D-05/
c...   SO2  --   626                              
      DATA (QCOEF(41,I),I=1,4)/-.24056D+03, .11101D+02,
     &                .22164D-01, .52334D-04/
c...   SO2  --   646                              
      DATA (QCOEF(42,I),I=1,4)/-.24167D+03, .11151D+02,
     &                .22270D-01, .52550D-04/
c...   NO2  --   646                              
      DATA (QCOEF(43,I),I=1,4)/-.53042D+03, .24216D+02,
     &                .66856D-01, .43823D-04/
c...   NH3  --  4111                              TIPS97.FOR modification
      DATA (QCOEF(44,I),I=1,4)/-.62293D+02, .30915D+01,
     &                .94575D-02, .18416D-05/
c...   NH3  --  5111                              TIPS97.FOR modification
      DATA (QCOEF(45,I),I=1,4)/-.42130D+02, .20569D+01,
     &                .63387D-02, .12127D-05/
c...  HNO3  --   146      *** NB these are the HITRAN'92 HNO3 values ***
c     DATA (QCOEF(46,I),I=1,4)/-.74208D+04, .34984D+03,
c    &                .89051D-01, .39356D-02/
C These are the new (limited Temp.range) coefficients suggested by Gamache,
C via D.Edwards
c       data (qcoef(46,i),i=1,4)/-.20026D+05, .59842D+03,
c     &                -.15008D+01, .72683D-02 /
C These are coefficients for 150-300K calculated by Castelli/Dinelli
       data (qcoef(46,i),i=1,4)/-4.0963D+04, 8.8476D+02,
     &                 -2.8347D+00,  9.2646D-03 /
c...    OH  --    61                              TIPS97.FOR modification
      DATA (QCOEF(47,I),I=1,4)/  .87390D+01, .15977D+00,
     &                .38291D-03,-.35669D-06/
c...    OH  --    81                              TIPS97.FOR modification
      DATA (QCOEF(48,I),I=1,4)/ .86770D+01, .16175D+00,
     &                .38223D-03,-.35466D-06/
c...    OH  --    62                              TIPS97.FOR modification
      DATA (QCOEF(49,I),I=1,4)/ .10239D+02, .43783D+00,
     &                .10477D-02,-.94570D-06/
c...    HF  --    19                              
      DATA (QCOEF(50,I),I=1,4)/ .15486D+01, .13350D+00,
     &                .59154D-05,-.46889D-08/
c...   HCl  --    15                              
      DATA (QCOEF(51,I),I=1,4)/ .28627D+01, .53122D+00,
     &                .67464D-05,-.16730D-08/
c...   HCl  --    17                              
      DATA (QCOEF(52,I),I=1,4)/ .28617D+01, .53203D+00,
     &                .66553D-05,-.15168D-08/
c...   HBr  --    19                              
      DATA (QCOEF(53,I),I=1,4)/ .27963D+01, .66532D+00,
     &                .34255D-05, .52274D-08/
c...   HBr  --    11                              
      DATA (QCOEF(54,I),I=1,4)/ .27953D+01, .66554D+00,
     &                .32931D-05, .54823D-08/
c...    HI  --    17                              
      DATA (QCOEF(55,I),I=1,4)/ .40170D+01, .13003D+01,
     &               -.11409D-04, .40026D-07/
c...   ClO  --    56                             TIPS97.FOR modification 
      DATA (QCOEF(56,I),I=1,4)/  .90968D+02, .70918D+01,
     &                .11639D-01, .30145D-05/
c...   ClO  --    76                             TIPS97.FOR modification  
      DATA (QCOEF(57,I),I=1,4)/ .92598D+02, .72085D+01,
     &                .11848D-01, .31305D-05/ 
c...   OCS  --   622                              
      DATA (QCOEF(58,I),I=1,4)/-.93697D+00, .36090D+01,
     &               -.34552D-02, .17462D-04/
c...   OCS  --   624                              
      DATA (QCOEF(59,I),I=1,4)/-.11536D+01, .37028D+01,
     &               -.35582D-02, .17922D-04/
c...   OCS  --   632                              
      DATA (QCOEF(60,I),I=1,4)/-.61015D+00, .72200D+01,
     &               -.70044D-02, .36708D-04/
c...   OCS  --   623                       Fitted for 150-300K by B. Dinelli
      DATA (QCOEF(61,I),I=1,4)/ .48948D+02, .14114D+02,
     &               -.01290D-00, .71479D-04/
c...   OCS  --   822                              
      DATA (QCOEF(62,I),I=1,4)/-.21569D+00, .38332D+01,
     &               -.36783D-02, .19177D-04/
c...   OCS  --   634                        Geisa isotope, no data avaialable 
      DATA (QCOEF(63,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...  H2CO  --   126                              
      DATA (QCOEF(64,I),I=1,4)/-.11760D+03, .46885D+01,
     &                .15088D-01, .35367D-05/
c...  H2CO  --   136                              
      DATA (QCOEF(65,I),I=1,4)/-.24126D+03, .96134D+01,
     &                .30938D-01, .72579D-05/
c...  H2CO  --   128                              
      DATA (QCOEF(66,I),I=1,4)/-.11999D+03, .52912D+01,
     &                .14686D-01, .43505D-05/
c...  HOCl  --   165                              
      DATA (QCOEF(67,I),I=1,4)/-.73640D+03, .34149D+02,
     &                .93554D-01, .67409D-04/
c...  HOCl  --   167                              
      DATA (QCOEF(68,I),I=1,4)/-.74923D+03, .34747D+02,
     &                .95251D-01, .68523D-04/
c...    N2  --    44                              
      DATA (QCOEF(69,I),I=1,4)/ .13684D+01, .15756D+01,
     &               -.18511D-04, .38960D-07/
c...   HCN  --   124                              
      DATA (QCOEF(70,I),I=1,4)/-.13992D+01, .29619D+01,
     &               -.17464D-02, .65937D-05/
c...   HCN  --   134                              
      DATA (QCOEF(71,I),I=1,4)/-.25869D+01, .60744D+01,
     &               -.35719D-02, .13654D-04/
c...   HCN  --   125                              
      DATA (QCOEF(72,I),I=1,4)/-.11408D+01, .20353D+01,
     &               -.12159D-02, .46375D-05/
c... CH3Cl  --   215                              
      DATA (QCOEF(73,I),I=1,4)/-.91416D+03, .34081D+02,
     &                .75461D-02, .17933D-03/
c... CH3Cl  --   217                              
      DATA (QCOEF(74,I),I=1,4)/-.92868D+03, .34621D+02,
     &                .76674D-02, .18217D-03/
c...  H2O2  --  1661                              
      DATA (QCOEF(75,I),I=1,4)/-.36499D+03, .13712D+02,
     &                .38658D-01, .23052D-04/
c...  C2H2  --  1221                              
      DATA (QCOEF(76,I),I=1,4)/-.83088D+01, .14484D+01,
     &               -.25946D-02, .84612D-05/
c...  C2H2  --  1231                              
      DATA (QCOEF(77,I),I=1,4)/-.66736D+02, .11592D+02,
     &               -.20779D-01, .67719D-04/
c...  C2H6  --  1221                      Fitted for 150-300K by B. Dinelli 
      DATA (QCOEF(78,I),I=1,4)/-.90815D+04, .20823D+03,
     &                -.43870D+00, 0.00220D+00/
c...  C2H6  --  1231                       Geisa isotope, no data avaialable
      DATA (QCOEF(79,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...   PH3  --  1111                              
      DATA (QCOEF(80,I),I=1,4)/-.15068D+03, .64718D+01,
     &                .12588D-01, .14759D-04/
c...  COF2  --   269                              
      DATA (QCOEF(81,I),I=1,4)/-.54180D+04, .18868D+03,
     &               -.33139D+00, .18650D-02/
c...   SF6  --    29                      Fitted for 150-300K by B. Dinelli
      DATA (QCOEF(82,I),I=1,4)/-.17727D+07, .29936D+05,
     &               -.16257D+03, .33901D+00/
c...   H2S  --   121                              
      DATA (QCOEF(83,I),I=1,4)/-.15521D+02, .83130D+00,
     &                .33656D-02,-.85691D-06/
c...   H2S  --   141                              
      DATA (QCOEF(84,I),I=1,4)/-.15561D+02, .83337D+00,
     &                .33744D-02,-.85937D-06/
c...   H2S  --   131                              
      DATA (QCOEF(85,I),I=1,4)/-.62170D+02, .33295D+01,
     &                .13480D-01,-.34323D-05/
c... HCOOH  --   126                              
      DATA (QCOEF(86,I),I=1,4)/-.29550D+04, .10349D+03,
     &               -.13146D+00, .87787D-03/
c...   HO2  --   166                              
      DATA (QCOEF(87,I),I=1,4)/-.15684D+03, .74450D+01,
     &                .26011D-01,-.92704D-06/
c...     O  --     6                         Set to +1 so that SQ=1 is returned
      DATA (QCOEF(88,I),I=1,4)/+.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...ClONO2  --  5646                              
      DATA (QCOEF(89,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...ClONO2  --  7646                              
      DATA (QCOEF(90,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...   NO+  --    46                               
      DATA (QCOEF(91,I),I=1,4)/ .91798D+00, .10416D+01,
     &               -.11614D-04, .24499D-07/
c...  HOBr  --   169                        Fitted for 150-300K by B. Dinelli
      DATA (QCOEF(92,I),I=1,4)/-.16598D+04, .56784D+02,
     &                .10030D+00, .16975D-03/
c...  HOBr  --   161                        Fitted for 150-300K by B. Dinelli
      DATA (QCOEF(93,I),I=1,4)/-.16531D+04, .56566D+02,
     &                .09980D+00, .16976D-03/
c...  C2H4  --  221                        Geisa isotope, no data avaialable
      DATA (QCOEF(94,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...  C2H4  --  231                        Geisa isotope, no data avaialable
      DATA (QCOEF(95,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...  HDO  --  162                          Repeat of H2O 162 isotope
      DATA (QCOEF(96,I),I=1,4)/-.23916D+02, .13793D+01,
     &                .61246D-02,-.21530D-05/
c...  BrO  --  96                        MARSCHALS isotope, no data avaialable
      DATA (QCOEF(97,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...  BrO  --  16                        MARSCHALS isotope, no data avaialable
      DATA (QCOEF(98,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...  C3H8  --  221                        Geisa isotope, no data avaialable
      DATA (QCOEF(99,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...  C2N2  --  224                        Geisa isotope, no data avaialable
      DATA (QCOEF(100,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...  C4H2  --  211                        Geisa isotope, no data avaialable
      DATA (QCOEF(101,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...  HC3N  --  124                        Geisa isotope, no data avaialable
      DATA (QCOEF(102,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...  C3H4  --  341?                       Geisa isotope, no data avaialable
      DATA (QCOEF(103,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...  GeH4  --  411                        
      DATA (QCOEF(104,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
      SAVE QCOEF
C
C...Total internal partition sums for T>500 to <=1500 K range:
c...   H2O  --   161                              
      DATA (QCOEFH( 1,I),I=1,4)/-.94327D+02, .81903D+00,
     &                .74005D-04, .42437D-06/
c...   H2O  --   181                              
      DATA (QCOEFH( 2,I),I=1,4)/-.95686D+02, .82839D+00,
     &                .68311D-04, .42985D-06/
c...   H2O  --   171                              
      DATA (QCOEFH( 3,I),I=1,4)/-.57133D+03, .49480D+01,
     &                .41517D-03, .25599D-05/
c...   H2O  --   162                              
      DATA (QCOEFH( 4,I),I=1,4)/-.53366D+03, .44246D+01,
     &               -.46935D-03, .29548D-05/
c...   H2O  --   182                       Geisa isotope, no data avaialable
      DATA (QCOEFH( 5,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...   H2O  --   172                       Geisa isotope, no data avaialable
      DATA (QCOEFH( 6,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...   CO2  --   626                              
      DATA (QCOEFH( 7,I),I=1,4)/-.50925D+03, .32766D+01,
     &               -.40601D-02, .40907D-05/
c...   CO2  --   636                              
      DATA (QCOEFH( 8,I),I=1,4)/-.11171D+04, .70346D+01,
     &               -.89063D-02, .88249D-05/
c...   CO2  --   628                              
      DATA (QCOEFH( 9,I),I=1,4)/-.11169D+04, .71299D+01,
     &               -.89194D-02, .89268D-05/
c...   CO2  --   627                              
      DATA (QCOEFH(10,I),I=1,4)/-.66816D+04, .42402D+02,
     &               -.53269D-01, .52774D-04/
c...   CO2  --   638                              
      DATA (QCOEFH(11,I),I=1,4)/-.25597D+04, .15855D+02,
     &               -.20440D-01, .19855D-04/
c...   CO2  --   637                              
      DATA (QCOEFH(12,I),I=1,4)/-.14671D+05, .91204D+02,
     &               -.11703D+00, .11406D-03/
c...   CO2  --   828                              
      DATA (QCOEFH(13,I),I=1,4)/-.73235D+04, .46140D+02,
     &               -.58473D-01, .57573D-04/
c...   CO2  --   728                              
      DATA (QCOEFH(14,I),I=1,4)/-.73235D+04, .46140D+02,
     &               -.58473D-01, .57573D-04/
c...   CO2  --   838                        Geisa isotope, no data avaialable
      DATA (QCOEFH(15,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...    O3  --   666                              
      DATA (QCOEFH(16,I),I=1,4)/-.11725D+05, .66515D+02,
     &               -.96010D-01, .94001D-04/
c...    O3  --   668                              
      DATA (QCOEFH(17,I),I=1,4)/-.25409D+05, .14393D+03,
     &               -.20850D+00, .20357D-03/
c...    O3  --   686                              
      DATA (QCOEFH(18,I),I=1,4)/-.12624D+05, .71391D+02,
     &               -.10383D+00, .10106D-03/
c...    O3  --   667                              
      DATA (QCOEFH(19,I),I=1,4)/-.14000D+06, .79825D+03,
     &               -.11465D+01, .11372D-02/
c...    O3  --   676                              
      DATA (QCOEFH(20,I),I=1,4)/-.69175D+05, .39442D+03,
     &               -.56650D+00, .56189D-03/
c...   N2O  --   446                              
      DATA (QCOEFH(21,I),I=1,4)/-.12673D+05, .75128D+02,
     &               -.10092D+00, .95557D-04/
c...   N2O  --   456                              
      DATA (QCOEFH(22,I),I=1,4)/-.90045D+04, .52833D+02,
     &               -.71771D-01, .67297D-04/
c...   N2O  --   546                              
      DATA (QCOEFH(23,I),I=1,4)/-.89960D+04, .53096D+02,
     &               -.71784D-01, .67592D-04/
c...   N2O  --   448                              
      DATA (QCOEFH(24,I),I=1,4)/-.13978D+05, .82338D+02,
     &               -.11167D+00, .10507D-03/
c...   N2O  --   447                              
      DATA (QCOEFH(25,I),I=1,4)/-.79993D+05, .47265D+03,
     &               -.63804D+00, .60218D-03/
c...    CO  --    26                              
      DATA (QCOEFH(26,I),I=1,4)/ .90723D+01, .33263D+00,
     &                .11806D-04, .27035D-07/
c...    CO  --    36                              
      DATA (QCOEFH(27,I),I=1,4)/ .20651D+02, .68810D+00,
     &                .34217D-04, .55823D-07/
c...    CO  --    28                              
      DATA (QCOEFH(28,I),I=1,4)/ .98497D+01, .34713D+00,
     &                .15290D-04, .28766D-07/
c...    CO  --    27                              
      DATA (QCOEFH(29,I),I=1,4)/ .58498D+02, .20351D+01,
     &                .87684D-04, .16554D-06/
c...    CO  --    38                              
      DATA (QCOEFH(30,I),I=1,4)/ .23511D+02, .71565D+00,
     &                .46681D-04, .58223D-07/
c...    CO  --    37                              
      DATA (QCOEFH(31,I),I=1,4)/ .11506D+03, .42727D+01,
     &                .17494D-03, .34413D-06/
c...   CH4  --   211                              
      DATA (QCOEFH(32,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...   CH4  --   311                              
      DATA (QCOEFH(33,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...   CH4  --   212                              
      DATA (QCOEFH(34,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...    O2  --    66                              
      DATA (QCOEFH(35,I),I=1,4)/ .36539D+02, .57015D+00,
     &                .16332D-03, .45568D-07/
c...    O2  --    68                              
      DATA (QCOEFH(36,I),I=1,4)/ .77306D+02, .11818D+01,
     &                .38661D-03, .89415D-07/
c...    O2  --    67                              
      DATA (QCOEFH(37,I),I=1,4)/ .44281D+03, .69531D+01,
     &                .21669D-02, .53053D-06/
c...    NO  --    46                              TIPS97.FOR modification 
      DATA (QCOEFH(38,I),I=1,4)/-.78837D+02, .39173D+01,
     &                .80657D-03, .22042D-06/
c...    NO  --    56                              TIPS97.FOR modification 
      DATA (QCOEFH(39,I),I=1,4)/-.67000D+02, .27874D+01,
     &                .45181D-03, .21161D-06/
c...    NO  --    48                              TIPS97.FOR modification 
      DATA (QCOEFH(40,I),I=1,4)/-.98460D+02, .42347D+01,
     &                .71550D-03, .32213D-06/
c...   SO2  --   626                              
      DATA (QCOEFH(41,I),I=1,4)/-.21162D+05, .11846D+03,
     &               -.16648D+00, .16825D-03/
c...   SO2  --   646                              
      DATA (QCOEFH(42,I),I=1,4)/-.21251D+05, .11896D+03,
     &               -.16717D+00, .16895D-03/
c...   NO2  --   646                              
      DATA (QCOEFH(43,I),I=1,4)/-.27185D+05, .16489D+03,
     &               -.19540D+00, .22024D-03/
c...   NH3  --  4111                              TIPS97.FOR modification
      DATA (QCOEFH(44,I),I=1,4)/-.59594D+04, .32387D+02,
     &               -.40459D-01, .31843D-04/
c...   NH3  --  5111                              TIPS97.FOR modification
      DATA (QCOEFH(45,I),I=1,4)/-.39908D+04, .21674D+02,
     &               -.27090D-01, .21308D-04/
c...  HNO3  --   146  *** NB these are the HITRAN'92 HNO3 values ***
      DATA (QCOEFH(46,I),I=1,4)/ .40718D+06,-.24214D+04,
     &                .64287D+01,-.10773D-02/
c...    OH  --    61                              TIPS97.FOR modification
      DATA (QCOEFH(47,I),I=1,4)/-.88840D+01, .30202D+00,
     &               -.15565D-04, .14330D-07/
c...    OH  --    81                              TIPS97.FOR modification
      DATA (QCOEFH(48,I),I=1,4)/-.51535D+01, .29076D+00,
     &               -.72340D-05, .20702D-07/
c...    OH  --    62                              TIPS97.FOR modification
      DATA (QCOEFH(49,I),I=1,4)/-.41683D+02, .83890D+00,
     &               -.36063D-04, .38083D-07/
c...    HF  --    19                              
      DATA (QCOEFH(50,I),I=1,4)/-.36045D-01, .14220D+00,
     &               -.10755D-04, .65523D-08/
c...   HCl  --    15                              
      DATA (QCOEFH(51,I),I=1,4)/ .25039D+01, .54430D+00,
     &               -.38656D-04, .39793D-07/
c...   HCl  --    17                              
      DATA (QCOEFH(52,I),I=1,4)/ .14998D+01, .54847D+00,
     &               -.42209D-04, .41029D-07/
c...   HBr  --    19                              
      DATA (QCOEFH(53,I),I=1,4)/ .67229D+01, .66356D+00,
     &               -.33749D-04, .54818D-07/
c...   HBr  --    11                              
      DATA (QCOEFH(54,I),I=1,4)/ .67752D+01, .66363D+00,
     &               -.33655D-04, .54823D-07/
c...    HI  --    17                              
      DATA (QCOEFH(55,I),I=1,4)/ .29353D+02, .12220D+01,
     &                .10209D-04, .10719D-06/
c...   ClO  --    56                              TIPS97.FOR modification
      DATA (QCOEFH(56,I),I=1,4)/ .22662D+03, .61093D+01,
     &                .14454D-01, .16928D-06/
c...   ClO  --    76                              TIPS97.FOR modification
      DATA (QCOEFH(57,I),I=1,4)/ .23304D+03, .61805D+01,
     &                .14797D-01, .16629D-06/ 
c...   OCS  --   622                              
      DATA (QCOEFH(58,I),I=1,4)/-.54125D+04, .29749D+02,
     &               -.44698D-01, .38878D-04/
c...   OCS  --   624                              
      DATA (QCOEFH(59,I),I=1,4)/-.55472D+04, .30489D+02,
     &               -.45809D-01, .39847D-04/
c...   OCS  --   632                              
      DATA (QCOEFH(60,I),I=1,4)/-.11863D+05, .64745D+02,
     &               -.98318D-01, .84563D-04/
c...   OCS  --   623                        Geisa isotope, no data avaialable
      DATA (QCOEFH(61,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...   OCS  --   822                              
      DATA (QCOEFH(62,I),I=1,4)/-.61288D+04, .33520D+02,
     &               -.50734D-01, .43792D-04/
c...   OCS  --   634                        Geisa isotope, no data avaialable
      DATA (QCOEFH(63,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...  H2CO  --   126                              
      DATA (QCOEFH(64,I),I=1,4)/-.17628D+05, .91794D+02,
     &               -.13055D+00, .89336D-04/
c...  H2CO  --   136                              
      DATA (QCOEFH(65,I),I=1,4)/-.36151D+05, .18825D+03,
     &               -.26772D+00, .18321D-03/
c...  H2CO  --   128                              
      DATA (QCOEFH(66,I),I=1,4)/-.17628D+05, .91794D+02,
     &               -.13055D+00, .89336D-04/
c...  HOCl  --   165                              
      DATA (QCOEFH(67,I),I=1,4)/-.24164D+05, .15618D+03,
     &               -.13206D+00, .21900D-03/
c...  HOCl  --   167                              
      DATA (QCOEFH(68,I),I=1,4)/-.24592D+05, .15895D+03,
     &               -.13440D+00, .22289D-03/
c...    N2  --    44                              
      DATA (QCOEFH(69,I),I=1,4)/ .27907D+02, .14972D+01,
     &               -.70424D-05, .11734D-06/
c...   HCN  --   124                              
      DATA (QCOEFH(70,I),I=1,4)/-.78078D+03, .61725D+01,
     &               -.53816D-02, .73379D-05/
c...   HCN  --   134                              
      DATA (QCOEFH(71,I),I=1,4)/-.16309D+04, .12801D+02,
     &               -.11242D-01, .15268D-04/
c...   HCN  --   125                              
      DATA (QCOEFH(72,I),I=1,4)/-.56301D+03, .43794D+01,
     &               -.38928D-02, .52467D-05/
c... CH3Cl  --   215                              
      DATA (QCOEFH(73,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c... CH3Cl  --   217                              
      DATA (QCOEFH(74,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...  H2O2  --  1661                              
      DATA (QCOEFH(75,I),I=1,4)/-.27583D+05, .15064D+03,
     &               -.19917D+00, .16977D-03/
c...  C2H2  --  1221                              
      DATA (QCOEFH(76,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...  C2H2  --  1231                              
      DATA (QCOEFH(77,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...  C2H6  --  1221                              
      DATA (QCOEFH(78,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...  C2H6  --  1231                       Geisa isotope, no data avaialable 
      DATA (QCOEFH(79,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...   PH3  --  1111                              
      DATA (QCOEFH(80,I),I=1,4)/-.28390D+05, .14463D+03,
     &               -.21473D+00, .14346D-03/
c...  COF2  --   269                              
      DATA (QCOEFH(81,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...   SF6  --    29                              
      DATA (QCOEFH(82,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...   H2S  --   121                              
      DATA (QCOEFH(83,I),I=1,4)/-.37572D+03, .29157D+01,
     &               -.98642D-03, .24113D-05/
c...   H2S  --   141                              
      DATA (QCOEFH(84,I),I=1,4)/-.37668D+03, .29231D+01,
     &               -.98894D-03, .24174D-05/
c...   H2S  --   131                              
      DATA (QCOEFH(85,I),I=1,4)/-.15049D+04, .11678D+02,
     &               -.39510D-02, .96579D-05/
c... HCOOH  --   126                              
      DATA (QCOEFH(86,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...   HO2  --   166                              
      DATA (QCOEFH(87,I),I=1,4)/-.32576D+04, .25539D+02,
     &               -.12803D-01, .29358D-04/
c...     O  --     6                         Set to +1 so that SQ=1 is returned
      DATA (QCOEFH(88,I),I=1,4)/+.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...ClONO2  --  5646                              
      DATA (QCOEFH(89,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...ClONO2  --  7646                              
      DATA (QCOEFH(90,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...   NO+  --    46                              
      DATA (QCOEFH(91,I),I=1,4)/ .17755D+02, .99262D+00,
     &               -.70814D-05, .76699D-07/
C... HOBr   --   169                     
      DATA (QCOEFH(92,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
C... HOBr   --   161                     
      DATA (QCOEFH(93,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...  C2H4  --  221                        Geisa isotope, no data avaialable
      DATA (QCOEFH(94,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...  C2H4  --  231                        Geisa isotope, no data avaialable
      DATA (QCOEFH(95,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...  HDO  --  162                          Repeat of H2O 162 isotope
      DATA (QCOEFH(96,I),I=1,4)/-.23916D+02, .13793D+01,
     &                .61246D-02,-.21530D-05/
c...  BrO  --  96                        MARSCHALS isotope, no data available
      DATA (QCOEFH(97,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...  BrO  --  16                        MARSCHALS isotope, no data available
      DATA (QCOEFH(98,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...  C3H8  --  221                        Geisa isotope, no data avaialable
      DATA (QCOEFH(99,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...  C2N2  --  224                        Geisa isotope, no data avaialable
      DATA (QCOEFH(100,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...  C4H2  --  211                        Geisa isotope, no data avaialable
      DATA (QCOEFH(101,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...  HC3N  --  124                        Geisa isotope, no data avaialable
      DATA (QCOEFH(102,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...  C3H4  --  341?                       Geisa isotope, no data avaialable
      DATA (QCOEFH(103,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
c...  GeH4  --  411                        Geisa isotope, no data avaialable
      DATA (QCOEFH(104,I),I=1,4)/-.10000D+01, .00000D+00,
     &                .00000D+00, .00000D+00/
      SAVE QCOEFH
C
      DATA Q296/ 
c           H2O 161,         181,         171,         162,         182,
     &  .174626D+03, .176141D+03, .105306D+04, .865122D+03, .874664D+03,  
c               172,     CO2 626,         636,         628,         627,
     &  .522018D+04, .286219D+03, .576928D+03, .607978D+03, .354389D+04,
c               638,         637,         828,         728;         838;
     &  .123528D+04, .714432D+04, .323407D+03, .376700D+04,        -1.0,
c            O3 666,
     &  .348186D+04,
c               668,         686,         667,         676;     N2O 446,
     &  .746207D+04, .364563D+04, .430647D+05, .212791D+05, .499183D+04,
c               456,         546,         448,         447;       CO 26,
     &  .334938D+04, .344940D+04, .526595D+04, .307008D+05, .107428D+03,
c                36,          28,          27,          38,          37;
     &  .224704D+03, .112781D+03, .661209D+03, .236447D+03, .138071D+04,
c           CH4 211,         311,         212;       O2 66,          68,
     &  .589908D+03, .117974D+04, .477061D+04, .215726D+03, .452188D+03,
c                67;      NO  46,          56,          48;     SO2 626,
     &  .263998D+04, .113243D+04, .785200D+03, .119417D+04, .634449D+04,
c               646;     NO2 646;    NH3 4111,        5111;    HNO3 146;
C new Gamache HNO3 coeffs give 0.214111 instead of 0.206001 
C Castelli/Dinelli coeffs give 0.212832
C TIPS97.FOR gives incorrect NH3 values for new coefficients, corrected here
     &  .637321D+04, .136318D+05, .172918D+04, .115353D+04, .212832D+06,
c             OH 61,         81,           62;       HF 19;      HCl 15,
     &  .803295D+02, .808460D+02, .207108D+03, .414625D+02, .160650D+03,
c                17;      HBr 19,          11;       HI 17;      ClO 56,
     &  .160887D+03, .200165D+03, .200227D+03, .388948D+03, .328810D+04,
c                76;     OCS 622,         624,         632,         623,
     &  .334560D+04, .121746D+04, .124793D+04, .247482D+04, .494961D+04, 
c               822,         634;
     &  .130948D+04,        -1.0,
c          H2CO 126,         136,         128;    HOCl 165,         167;
     &  .268388D+04, .550322D+04, .284573D+04, .193166D+05, .196584D+05,
c             N2 44;     HCN 124,         134,         125;   CH3Cl 215,
     &  .467136D+03, .893323D+03, .183657D+04, .615046D+03, .144858D+05,
c               217;   H2O2 1661;   C2H2 1221,        1231,   C2H6 1221,
     &  .147153D+05, .767871D+04, .412519D+03, .330014D+04, .711730D+05,
c              1231;
     &         -1.0,
c          PH3 1111;    COF2 269;      SF6 29;     H2S 121,         141,
     &  .325067D+04, .697632D+05, .163662D+07, .503204D+03, .504486D+03,
c               131;   HCOOH 126;     HO2 166;         O 6; ClONO2 5646,
     &  .201546D+04, .389257D+05, .430184D+04,+.100000D+01, .212829D+07,
c              7646;      NO+ 46;    HOBr 169,         161;
     &  .218246D+07, .308855D+03, .283385D+05, .282379D+05, 
c          C2H4 221,         231;     HDO 162;      BrO 96;          16;
     &         -1.0,        -1.0, .865122D+03,        -1.0,        -1.0,
c          C3H8 221;    C2N2 224;    C4H2 211;    HC3N 124;    C3H4 341;
     &         -1.0,        -1.0,        -1.0,        -1.0,        -1.0,
c          GeH4 411
     &         -1.0 /
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      DBLTEM = DBLE ( TEM )
      IF ( IDGAS .GT. MAXH42 ) THEN  ! Assume cross-section molecule
        SQ = ( TEMREF / TEM )**1.5
        RETURN
      ELSE IF ( TPSFLG .AND. IDXTPS(ISO,IDGAS) .GT. 0 ) THEN
        SQ = TPSLKP ( IDGAS, ISO, TEM )
      ELSE
        ISPC = ISPIDG(IDGAS) + ISO     ! Position in array of molecular isotope
        QSTD = Q296(ISPC)              ! Total internal partition sum at 296K
        IF ( TEM .LT. 70.0 .OR. TEM .GT. 1500.0 ) 
     &    WRITE (*,*) 'W-QTFCT: Extrap.beyond 70-1500K limits, T=', TEM
        IF ( TEM .LE. 500.0 ) THEN
          QT   =   QCOEF(ISPC,1)  
     &           + DBLTEM * ( QCOEF(ISPC,2) 
     &           + DBLTEM * ( QCOEF(ISPC,3) 
     &           + DBLTEM *   QCOEF(ISPC,4) ) ) 
        ELSE 
          QT   =   QCOEFH(ISPC,1)  
     &           + DBLTEM * ( QCOEFH(ISPC,2) 
     &           + DBLTEM * ( QCOEFH(ISPC,3) 
     &           + DBLTEM *   QCOEFH(ISPC,4) ) ) 
        END IF
C
C For species for which no data is available, QT=-1 is returned.
C For these just use T^(-3/2)
C For monatomic O, QT=Q296=+1 is returned, so SQ=1
        IF ( QT .EQ. -1.0D0 ) THEN
          SQ = ( TEMREF / TEM )**1.5
        ELSE
          SQ = SNGL ( QSTD / QT )
        END IF
      END IF
C
      END
