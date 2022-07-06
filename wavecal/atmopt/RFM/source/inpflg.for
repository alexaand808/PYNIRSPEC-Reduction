      SUBROUTINE INPFLG ( LUNDRV, FAIL, ERRMSG )
C
C VERSION
C     08JAN13 AD Allow combination of HOM + SFC flags
C     08-AUG-13  AD  Add V42 flag
C                    Remove specific checks for obsolete flags: GL2,IQD,OFM,SCA
C     11-SEP-11  AD  Add RJT flag
C     10-AUG-05  AD  Make HOM and SFC flags incompatible
C     20-JAN-05  AD  Make LUT and NTE flags incompatible
C     06-JAN-04  AD  Add PRF flag
C     29-DEC-03  AD  Add LEV flag
C     11-AUG-03  AD  Add VVW flag
C     08-AUG-03  AD  Add BBT and GHZ flags
C     12-FEB-03  AD  Add TPS flag
C     04-NOV-01  AD  Allow FLX+VRT+OPT, forbid FLX+VRT+COO
C     31-OCT-01  AD  Add VRT flag
C     07-JUN-01  AD  Add LOS flag
C     30-MAY-00  AD  Make GRD and TAB flags compatible.
C     26-MAY-00  AD  Add DBL flag
C     24-APR-00  AD  Add COO, FLX and MTX flags
C     05-JAN-00  AD  Add GRA, AVG flags.
C     13-SEP-99  AD  Allow JAC and NTE flags to be used together
C     02-FEB-99  AD  Remove SCAFLG. Add CLCFLG, LUNFLG.
C     24-JAN-99  AD  Remove OFMFLG, GL2FLG, IQDFLG.
C                    Add LINFLG, FVZFLG, QADFLG.
C     30-DEC-98  AD  Add JAC flag
C     29-DEC-98  AD  Add BFX flag
C     18-JUL-98  AD  Require SFC flag if NAD flag used
C     13-JUL-98  AD  Add SFC flag.
C     17-DEC-97  AD  Add GRD flag.
C     21-OCT-97  AD  Add TAB flag.
C     03-JUL-97  AD  Add LUT flag.
C     03-MAR-97  AD  Version 3.
C     16-JAN-97  AD  Add new flags OBS, ZEN
C     01-OCT-96  AD  Version 2.
C     30-SEP-96  AD  Reorder checks for incompatible flags.
C     28-SEP-96  AD  Add new flag IQD
C     17-SEP-96  AD  Add new flags CRV and FIN. Rename VIB flag to NTE flag
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Read RFM option flags data following *FLG marker
C     Called once by RFMINP.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      LUNDRV  !  I  LUN for Driver File
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  LOCASE ! Convert text string to lower case
     &, NXTFLD ! Load next field from section of RFM driver file
     &, RFMLOG ! Write text message to RFM log file
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
C
C LOCAL CONSTANTS
      INTEGER       MAXFLG       ! No. Flags listed in /FLGCOM/
        PARAMETER ( MAXFLG = 51 )
C
C LOCAL VARIABLES
      INTEGER       IFLG    ! Flag counter (1:MAXFLG)
      INTEGER       IPT     ! Pointer to locations in MESSGE ENATXT strings
      INTEGER       LENGTH  ! Length of character field read from driver table
      LOGICAL       SETFLG(MAXFLG)
      CHARACTER*3   FLGLST(MAXFLG)
      CHARACTER*80  FIELD   ! Field read from driver table
      CHARACTER*3   FLG     ! Flag read from driver table
      CHARACTER*132 MESSGE  ! Message listing enabled flags for LOG file
        DATA SETFLG / MAXFLG * .FALSE. /
        DATA FLGLST / 'abs', 'avg', 'bbt', 'bfx', 'bin', 'clc', 
     &                'coo', 'crv', 'ctm', 'dbl', 'fin', 'flx', 
     &                'fov', 'fvz', 'geo', 'ghz', 'gra', 'grd', 
     &                'hom', 'ils', 'jac', 'lay', 'lev', 'lin', 
     &                'los', 'lun', 'lut', 'mix', 'mtx', 'nad', 
     &                'new', 'nte', 'obs', 'opt', 'prf', 'pth', 
     &                'qad', 'rad', 'rej', 'rjt', 'sfc', 'shp', 
     &                'tab', 'tps', 'tra', 'v42', 'vrt', 'vvw', 
     &                'wid', 'wng', 'zen' /
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Read next field in *FLG section of driver table
C
  100 CALL NXTFLD ( LUNDRV, FIELD, LENGTH, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
      IF ( LENGTH .NE. 0 ) THEN                ! Not yet reached end of section
        IF ( LENGTH .NE. 3 ) THEN
          FAIL = .TRUE.
          LENGTH = MIN ( LENGTH, 40 )
          ERRMSG = 'F-INPFLG: Not a C*3 string, field='//FIELD(1:LENGTH)
          RETURN
        END IF
        CALL LOCASE ( FIELD(1:3), FLG )
        DO IFLG = 1, MAXFLG
          IF ( FLG .EQ. FLGLST(IFLG) ) THEN
            SETFLG(IFLG) = .TRUE.
            GOTO 100
          END IF
        END DO
        FAIL = .TRUE.
        ERRMSG = 'F-INPFLG: Unrecognised flag='//FLG
        RETURN
      END IF
C
      ABSFLG = SETFLG(1)
      AVGFLG = SETFLG(2)
      BBTFLG = SETFLG(3)
      BFXFLG = SETFLG(4)
      BINFLG = SETFLG(5)
      CLCFLG = SETFLG(6)
      COOFLG = SETFLG(7)
      CRVFLG = SETFLG(8)
      CTMFLG = SETFLG(9)
      DBLFLG = SETFLG(10)
      FINFLG = SETFLG(11)
      FLXFLG = SETFLG(12)
      FOVFLG = SETFLG(13)
      FVZFLG = SETFLG(14)
      GEOFLG = SETFLG(15)
      GHZFLG = SETFLG(16)
      GRAFLG = SETFLG(17)
      GRDFLG = SETFLG(18)
      HOMFLG = SETFLG(19)
      ILSFLG = SETFLG(20)
      JACFLG = SETFLG(21)
      LAYFLG = SETFLG(22)
      LEVFLG = SETFLG(23)
      LINFLG = SETFLG(24)
      LOSFLG = SETFLG(25)
      LUNFLG = SETFLG(26)
      LUTFLG = SETFLG(27)
      MIXFLG = SETFLG(28)
      MTXFLG = SETFLG(29)
      NADFLG = SETFLG(30)
      NEWFLG = SETFLG(31)
      NTEFLG = SETFLG(32)
      OBSFLG = SETFLG(33)
      OPTFLG = SETFLG(34)
      PRFFLG = SETFLG(35)
      PTHFLG = SETFLG(36)
      QADFLG = SETFLG(37)
      RADFLG = SETFLG(38)
      REJFLG = SETFLG(39)
      RJTFLG = SETFLG(40)
      SFCFLG = SETFLG(41)
      SHPFLG = SETFLG(42)
      TABFLG = SETFLG(43)
      TPSFLG = SETFLG(44)
      TRAFLG = SETFLG(45)
      V42FLG = SETFLG(46)
      VRTFLG = SETFLG(47)
      VVWFLG = SETFLG(48)
      WIDFLG = SETFLG(49)
      WNGFLG = SETFLG(50)
      ZENFLG = SETFLG(51)
C
C Check for inconsistent combinations of flags
C
      FAIL = .TRUE.
      IF ( ABSFLG .AND. TABFLG ) THEN
        ERRMSG = 'F-INPFLG: ABS and TAB flags are incompatible'
C
      ELSE IF ( AVGFLG .AND. ILSFLG ) THEN
        ERRMSG = 'F-INPFLG: AVG and ILS flags are incompatible'
      ELSE IF ( AVGFLG .AND. OPTFLG ) THEN
        ERRMSG = 'F-INPFLG: AVG and OPT flags are incompatible'
      ELSE IF ( AVGFLG .AND. TABFLG ) THEN
        ERRMSG = 'F-INPFLG: AVG and TAB flags are incompatible'
C
      ELSE IF ( BBTFLG .AND. TABFLG ) THEN
        ERRMSG = 'F-INPFLG: BBT and TAB flags are incompatible'
C
      ELSE IF ( BFXFLG .AND. HOMFLG ) THEN
        ERRMSG = 'F-INPFLG: BFX and HOM flags are incompatible'
      ELSE IF ( BFXFLG .AND. NTEFLG ) THEN
        ERRMSG = 'F-INPFLG: BFX and NTE flags are incompatible'
      ELSE IF ( BFXFLG .AND. TABFLG ) THEN
        ERRMSG = 'F-INPFLG: BFX and NTE flags are incompatible'
C
      ELSE IF ( CLCFLG .AND. FLXFLG ) THEN
        ERRMSG = 'F-INPFLG: CLC and FLX flags are incompatible'
      ELSE IF ( CLCFLG .AND. HOMFLG ) THEN
        ERRMSG = 'F-INPFLG: CLC and HOM flags are incompatible'
      ELSE IF ( CLCFLG .AND. NADFLG ) THEN
        ERRMSG = 'F-INPFLG: CLC and NAD flags are incompatible'
      ELSE IF ( CLCFLG .AND. TABFLG ) THEN
        ERRMSG = 'F-INPFLG: CLC and TAB flags are incompatible'
      ELSE IF ( CLCFLG .AND. ZENFLG ) THEN
        ERRMSG = 'F-INPFLG: CLC and HOM flags are incompatible'
C
      ELSE IF ( COOFLG .AND. VRTFLG ) THEN
        ERRMSG = 'F-INPFLG: COO and VRT flags are incompatible'
C
      ELSE IF ( CRVFLG .AND. FLXFLG ) THEN
        ERRMSG = 'F-INPFLG: CRV and FLX flags are incompatible'
      ELSE IF ( CRVFLG .AND. HOMFLG ) THEN
        ERRMSG = 'F-INPFLG: CRV and HOM flags are incompatible'
      ELSE IF ( CRVFLG .AND. NADFLG ) THEN
        ERRMSG = 'F-INPFLG: CRV and NAD flags are incompatible'
      ELSE IF ( CRVFLG .AND. TABFLG ) THEN
        ERRMSG = 'F-INPFLG: CRV and TAB flags are incompatible'
      ELSE IF ( CRVFLG .AND. ZENFLG ) THEN
        ERRMSG = 'F-INPFLG: CRV and ZEN flags are incompatible'
C
      ELSE IF ( DBLFLG .AND. TABFLG ) THEN
        ERRMSG = 'F-INPFLG: DBL and TAB flags are incompatible'
C
C NB: COO, VRT and MTX all require FLX so no need to test these options also
      ELSE IF ( FLXFLG .AND. FOVFLG ) THEN
        ERRMSG = 'F-INPFLG: FLX and FOV flags are incompatible'
      ELSE IF ( FLXFLG .AND. GEOFLG ) THEN
        ERRMSG = 'F-INPFLG: FLX and GEO flags are incompatible'
      ELSE IF ( FLXFLG .AND. GRAFLG ) THEN
        ERRMSG = 'F-INPFLG: FLX and GRA flags are incompatible'
      ELSE IF ( FLXFLG .AND. HOMFLG ) THEN
        ERRMSG = 'F-INPFLG: FLX and HOM flags are incompatible'
      ELSE IF ( FLXFLG .AND. LEVFLG ) THEN
        ERRMSG = 'F-INPFLG: FLX and LEV flags are incompatible'
      ELSE IF ( FLXFLG .AND. LOSFLG ) THEN
        ERRMSG = 'F-INPFLG: FLX and LOS flags are incompatible'
      ELSE IF ( FLXFLG .AND. OBSFLG ) THEN
        ERRMSG = 'F-INPFLG: FLX and OBS flags are incompatible'
      ELSE IF ( FLXFLG .AND. OPTFLG .AND. .NOT. VRTFLG ) THEN
        ERRMSG = 'F-INPFLG: FLX and OPT flags also require VRT flag'
      ELSE IF ( FLXFLG .AND. PTHFLG ) THEN
        ERRMSG = 'F-INPFLG: FLX and PTH flags are incompatible'
      ELSE IF ( FLXFLG .AND. TABFLG ) THEN
        ERRMSG = 'F-INPFLG: FLX and TAB flags are incompatible'
C
      ELSE IF ( FOVFLG .AND. HOMFLG ) THEN
        ERRMSG = 'F-INPFLG: FOV and HOM flags are incompatible'
      ELSE IF ( FOVFLG .AND. NADFLG ) THEN
        ERRMSG = 'F-INPFLG: FOV and NAD flags are incompatible'
      ELSE IF ( FOVFLG .AND. TABFLG ) THEN
        ERRMSG = 'F-INPFLG: FOV and TAB flags are incompatible'
      ELSE IF ( FOVFLG .AND. ZENFLG ) THEN
        ERRMSG = 'F-INPFLG: FOV and ZEN flags are incompatible'
      ELSE IF ( FOVFLG .AND. OPTFLG ) THEN
        ERRMSG = 'F-INPFLG: FOV and OPT flags are incompatible'
C
      ELSE IF ( GEOFLG .AND. HOMFLG ) THEN
        ERRMSG = 'F-INPFLG: GEO and HOM flags are incompatible'
      ELSE IF ( GEOFLG .AND. NADFLG ) THEN
        ERRMSG = 'F-INPFLG: GEO and NAD flags are incompatible'
      ELSE IF ( GEOFLG .AND. TABFLG ) THEN
        ERRMSG = 'F-INPFLG: GEO and TAB flags are incompatible'
      ELSE IF ( GEOFLG .AND. ZENFLG ) THEN
        ERRMSG = 'F-INPFLG: GEO and ZEN flags are incompatible'
C
      ELSE IF ( GRAFLG .AND. HOMFLG ) THEN
        ERRMSG = 'F-INPFLG: GRA and HOM flags are incompatible'
      ELSE IF ( GRAFLG .AND. NADFLG ) THEN
        ERRMSG = 'F-INPFLG: GRA and NAD flags are incompatible'
      ELSE IF ( GRAFLG .AND. TABFLG ) THEN
        ERRMSG = 'F-INPFLG: GRA and TAB flags are incompatible'
      ELSE IF ( GRAFLG .AND. ZENFLG ) THEN
        ERRMSG = 'F-INPFLG: GRA and ZEN flags are incompatible'
C
      ELSE IF ( GRDFLG .AND. WNGFLG ) THEN
        ERRMSG = 'F-INPFLG: GRD and WNG flags are incompatible'
C
      ELSE IF ( HOMFLG .AND. LAYFLG ) THEN
        ERRMSG = 'F-INPFLG: HOM and LAY flags are incompatible'
      ELSE IF ( HOMFLG .AND. LEVFLG ) THEN
        ERRMSG = 'F-INPFLG: HOM and LEV flags are incompatible'
      ELSE IF ( HOMFLG .AND. LOSFLG ) THEN
        ERRMSG = 'F-INPFLG: HOM and LOS flags are incompatible'
      ELSE IF ( HOMFLG .AND. NADFLG ) THEN
        ERRMSG = 'F-INPFLG: HOM and NAD flags are incompatible'
      ELSE IF ( HOMFLG .AND. OBSFLG ) THEN
        ERRMSG = 'F-INPFLG: HOM and OBS flags are incompatible'
c      ELSE IF ( HOMFLG .AND. SFCFLG ) THEN
c        ERRMSG = 'F-INPFLG: HOM and SFC flags are incompatible'
      ELSE IF ( HOMFLG .AND. ZENFLG ) THEN
        ERRMSG = 'F-INPFLG: HOM and ZEN flags are incompatible'
C
      ELSE IF ( ILSFLG .AND. OPTFLG ) THEN
        ERRMSG = 'F-INPFLG: ILS and OPT flags are incompatible'
      ELSE IF ( ILSFLG .AND. TABFLG ) THEN
        ERRMSG = 'F-INPFLG: ILS and TAB flags are incompatible'
C
      ELSE IF ( JACFLG .AND. LEVFLG ) THEN
        ERRMSG = 'F-INPFLG: JAC and LEV flags are incompatible'
      ELSE IF ( JACFLG .AND. MTXFLG ) THEN
        ERRMSG = 'F-INPFLG: JAC and MTX flags are incompatible'
      ELSE IF ( JACFLG .AND. TABFLG ) THEN
        ERRMSG = 'F-INPFLG: JAC and TAB flags are incompatible'
C
      ELSE IF ( LEVFLG .AND. LOSFLG ) THEN
        ERRMSG = 'F-INPFLG: LEV and LOS flags are incompatible'
      ELSE IF ( LEVFLG .AND. TABFLG ) THEN
        ERRMSG = 'F-INPFLG: LEV and TAB flags are incompatible'
C 
      ELSE IF ( LOSFLG .AND. MTXFLG ) THEN
        ERRMSG = 'F-INPFLG: LOS and MTX flags are incompatible'
      ELSE IF ( LOSFLG .AND. NADFLG ) THEN
        ERRMSG = 'F-INPFLG: LOS and NAD flags are incompatible'
      ELSE IF ( LOSFLG .AND. TABFLG ) THEN
        ERRMSG = 'F-INPFLG: LOS and TAB flags are incompatible'
      ELSE IF ( LOSFLG .AND. ZENFLG ) THEN
        ERRMSG = 'F-INPFLG: LOS and ZEN flags are incompatible'
C
      ELSE IF ( LUTFLG .AND. NTEFLG ) THEN
        ERRMSG = 'F-INPFLG: LUT and NTE flags are incompatible'
C
      ELSE IF ( MTXFLG .AND. NADFLG ) THEN
        ERRMSG = 'F-INPFLG: MTX and NAD flags are incompatible'
      ELSE IF ( MTXFLG .AND. ZENFLG ) THEN
        ERRMSG = 'F-INPFLG: MTX and ZEN flags are incompatible'
C
      ELSE IF ( NADFLG .AND. TABFLG ) THEN
        ERRMSG = 'F-INPFLG: NAD and TAB flags are incompatible'
      ELSE IF ( NADFLG .AND. ZENFLG ) THEN
        ERRMSG = 'F-INPFLG: NAD and ZEN flags are incompatible'
C
      ELSE IF ( OBSFLG .AND. TABFLG ) THEN
        ERRMSG = 'F-INPFLG: OBS and TAB flags are incompatible'
C
      ELSE IF ( OPTFLG .AND. TABFLG ) THEN
        ERRMSG = 'F-INPFLG: OPT and TAB flags are incompatible'
C
      ELSE IF ( PTHFLG .AND. TABFLG ) THEN
        ERRMSG = 'F-INPFLG: PTH and TAB flags are incompatible'
C
      ELSE IF ( RADFLG .AND. TABFLG ) THEN
        ERRMSG = 'F-INPFLG: RAD and TAB flags are incompatible'
C
      ELSE IF ( RJTFLG .AND. TABFLG ) THEN
        ERRMSG = 'F-INPFLG: RJT and TAB flags are incompatible'
C
      ELSE IF ( SFCFLG .AND. TABFLG ) THEN
        ERRMSG = 'F-INPFLG: SFC and TAB flags are incompatible'
C
      ELSE IF ( SFCFLG .AND. ZENFLG ) THEN
        ERRMSG = 'F-INPFLG: SFC and ZEN flags are incompatible'
C
      ELSE IF ( TABFLG .AND. TRAFLG ) THEN
        ERRMSG = 'F-INPFLG: TAB and TRA flags are incompatible'
      ELSE IF ( TABFLG .AND. ZENFLG ) THEN
        ERRMSG = 'F-INPFLG: TAB and ZEN flags are incompatible'
C
      ELSE IF ( BFXFLG .AND. .NOT. RADFLG ) THEN
        ERRMSG = 'F-INPFLG: BFX flag also requires RAD flag enabled'
      ELSE IF ( COOFLG .AND. .NOT. FLXFLG ) THEN
        ERRMSG = 'F-INPFLG: COO flag also requires FLX flag enabled'
      ELSE IF ( FINFLG .AND. .NOT. ( ILSFLG .OR. AVGFLG ) ) THEN
        ERRMSG = 'F-INPFLG: FIN flag also requires either '//
     &           'ILS or AVG flag enabled'
      ELSE IF ( FLXFLG .AND. TRAFLG .AND. .NOT. MTXFLG 
     &          .AND. .NOT. ( NADFLG .OR. ZENFLG ) ) THEN
        ERRMSG = 'F-INPFLG: TRA+FLX-MTX flags also requires ZEN or NAD'
      ELSE IF ( FLXFLG .AND. ABSFLG .AND. .NOT. MTXFLG 
     &          .AND. .NOT. ( NADFLG .OR. ZENFLG ) ) THEN
        ERRMSG = 'F-INPFLG: ABS+FLX-MTX flags also requires ZEN or NAD'
      ELSE IF ( FLXFLG .AND. .NOT. ZENFLG .AND. .NOT. SFCFLG ) THEN
        ERRMSG = 'F-INPFLG: FLX-ZEN flags also requires SFC flag'
      ELSE IF ( FVZFLG .AND. .NOT. FOVFLG ) THEN
        ERRMSG = 'F-INPFLG: FVZ flag also requires FOV flag enabled'
      ELSE IF ( MTXFLG .AND. .NOT. FLXFLG ) THEN
        ERRMSG = 'F-INPFLG: MTX flag also requires FLX flag enabled'
      ELSE IF ( MTXFLG .AND. COOFLG .AND. .NOT. SFCFLG ) THEN
        ERRMSG = 'F-INPFLG: MTX+COO flags also requires SFC flag'
      ELSE IF ( MTXFLG .AND. RADFLG .AND. .NOT. SFCFLG ) THEN
        ERRMSG = 'F-INPFLG: MTX+RAD flags also requires SFC flag'
      ELSE IF ( NADFLG .AND. .NOT. SFCFLG ) THEN
        ERRMSG = 'F-INPFLG: NAD flag also requires SFC flag enabled'
      ELSE IF ( VRTFLG .AND. .NOT. FLXFLG ) THEN
        ERRMSG = 'F-INPFLG: VRT flag also requires FLX flag enabled'
      ELSE
        FAIL = .FALSE.
      END IF
      IF ( FAIL ) RETURN     
C
C Construct message for log file listing status of all flags
C
      MESSGE = 'I-INPFLG: Enabled flags: '      
      IPT = 22
      DO IFLG = 1, MAXFLG
        IF ( SETFLG(IFLG) ) THEN
          IPT = IPT + 4
          MESSGE(IPT:IPT+3) = FLGLST(IFLG)//' '
        END IF
      END DO
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
C
      END
