      SUBROUTINE RFMINP ( LUNDRV, LUNNXT, FAIL, ERRMSG )
C
C VERSION
C     17OCT13 AD Move tests for duplicate sections here rather than in
C                INP*, also add SKPSEC  
C     14AUG13 AD Correction: add LUNNXT argument to INPSFC
C     19SEP12 AD Allow for records up to C*200 instead of C*80
C     13SEP11 AD Add INPRJT
C     05JUL07 AD Bug#65 Allow for improperly terminated driver tables
C     06JAN04 AD Add INPPRF
C     01JAN04 AD Add INPLEV
C     09AUG03 AD Add INPBBT
C     12FEB03 AD Add GETTPS
C     13JUN02 AD Remove GETHIT, GETXSC arguments to INPLUT
C     18MAR01 AD Change arguments to INPFIN, INPFOV, INPILS. Add INPCHK
C     15MAR01 AD Add GETFIN argument to INPILS, GETILS arg to INPFIN
C     15APR00 AD Add INPCOO and allow for *LEV instead of *TAN
C     31DEC98 AD Add INPJAC
C     31JUL98 AD Add GETHIT, GETXSC arguments to INPLUT
C     18JUL98 AD Add GETOBS,GETFOV arguments to INPFOV,INPOBS respectively
C     15JUL98 AD Allow for *DIM,*ELE,*GEO,*LEN or *SEC instead of *TAN
C     14JUL98 AD Add INPOBS
C     13JUL98 AD Add INPSFC
C     18DEC97 AD Add INPGRD
C     04DEC97 AD Add check for duplicate primary sections
C     22OCT97 AD Add INPTAB
C     03JUL97 AD Add INPLUT
C     03MAR97 AD Version 3.
C     27FEB97 AD Allow lower or upper case keywords
C     01OCT96 AD Version 2.
C     20SEP96 AD Add INPCRV to handle *CRV section and INPFIN for *FIN.
C                    Rename *VIB section to *NTE section
C     01SEP96 AD Version 1.
C
C DESCRIPTION
C     Open RFM Driver file and read input data.
C     Called once by RFM.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      LUNDRV  !  I  LUN for Driver File
      INTEGER      LUNNXT  ! I/O Next free LUN 
      LOGICAL      FAIL    !  O  Set TRUE if an error is detected
      CHARACTER*80 ERRMSG  !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  INPATM ! Read RFM atmospheric profiles from driver table 
     &, INPCHK ! Cross-check inputs from more than one drv table section
     &, INPCRV ! Read user-defined local radius of curvature.
     &, INPFIN ! Read user-defined fine-mesh resolution
     &, INPFLG ! Read RFM option flags from driver table
     &, INPFOV ! Read FOV data file from RFM driver table and load.
     &, INPGAS ! Read RFM absorbing gases from driver table 
     &, INPGRD ! Read list of GRD files of user-defined irregular grids
     &, INPHDR ! Read RFM Driver Header record from driver table
     &, INPHIT ! Read HITRAN line data filename from RFM driver table.
     &, INPILS ! Read list of ILS data files from driver table
     &, INPJAC ! Read details of Jacobians to be calculated
     &, INPLEV ! Read altitude levels for intermediate outputs.
     &, INPLUT ! Read list of Look-Up tables of absorption coefficients.
     &, INPNAM ! Read user-defined output filenames from RFM driver table.
     &, INPNTE ! Read non-LTE files from driver table 
     &, INPOBS ! Read user-defined observer position
C
        EXTERNAL
     &  INPREJ ! Read minimum line strength limits from RFM driver table.
     &, INPSFC ! Read RFM surface parameters from driver table 
     &, INPSHP ! Read RFM line shape model from driver table 
     &, INPSPC ! Read RFM wavenumber range and resolution from driver table 
     &, INPTAN ! Read RFM tangent points from driver table 
     &, INPTPS ! Read tabulated Total Internal Partition Sum data
     &, INPXSC ! Read RFM input Cross-section data files from driver table
     &, LOCASE ! Convert text string to lower case.
     &, RFMLOG ! Write text message to RFM log file
     &, SKPSEC ! Skip driver table section
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmfil.inc' ! Standard filenames of RFM I/O files
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C LOCAL CONSTANTS
      INTEGER       MINKEY        ! Minimum number of keyword labels required
        PARAMETER ( MINKEY = 6 )
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'ntecom.inc' ! Non-LTE data 
C
C LOCAL VARIABLES
      INTEGER      IKEY   ! Counter for mandatory keywords
      INTEGER      IOS    ! Saved value of IOSTAT for error messages
      LOGICAL      LDUMMY ! T = record pointer just before start of section
      LOGICAL      GETHIT ! T = HITRAN data required (set in INPGAS)
      LOGICAL      GETNTE ! T = NTE data required
      LOGICAL      GETXSC ! T = X-Secn data required (set in INPGAS)  
      LOGICAL      GOTABS !1 T = Absorption spectrum filename found             
      LOGICAL      GOTBBT !2 T = Black Body Temperature spectrum filename found 
      LOGICAL      GOTCOO !3 T = Cooling Rate spectrum filename found 
      LOGICAL      GOTCRV !4 T = Radius of Curvature section found
      LOGICAL      GOTFIN !5 T = Fine-mesh section found
      LOGICAL      GOTFOV !6 T = FOV data section found
      LOGICAL      GOTGRD !7 T = GRD data section found
      LOGICAL      GOTHIT !8 T = HITRAN filename section found
      LOGICAL      GOTILS !9 T = ILS data section found
      LOGICAL      GOTJAC !10 T = Jacobian section found
      LOGICAL      GOTLEV !11 T = Level data section found
      LOGICAL      GOTLUT !12 T = Look-Up table section found
      LOGICAL      GOTNTE !13 T = Non-LTE data section found
      LOGICAL      GOTOBS !14 T = Observer Position data section found
      LOGICAL      GOTOPT !15 T = Optical thickness spectrum filename found
      LOGICAL      GOTPRF !16 T = Profile output filename found
      LOGICAL      GOTPTH !17 T = Path diagnostic filename found
      LOGICAL      GOTRAD !18 T = Radiance spectrum filename found
      LOGICAL      GOTREJ !19 T = Min.line strength limits found
      LOGICAL      GOTRJT !20 T = Rayleigh-Jeans Temp spectrum filename found
      LOGICAL      GOTSFC !21 T = Surface data section found
      LOGICAL      GOTSHP !22 T = Lineshape data section found 
      LOGICAL      GOTTAB !23 T = Abs.Coeff table filename found
      LOGICAL      GOTTPS !24 T = Tabulated TIPS data section found
      LOGICAL      GOTTRA !25 T = Transmittance spectrum filename found
      LOGICAL      GOTWID !26 T = Widemesh diagnostic filename found
      LOGICAL      GOTXSC !27 T = X-Secn data section found
      CHARACTER*200 RECORD ! Record read from Driver File
      CHARACTER*80 NAMVAR ! Driver Table name converted from param.to variable
      CHARACTER*4  KEYCHK ! Part of driver table record tested for key
      CHARACTER*4  KEYLST(MINKEY)  ! List of mandatory keywords
      CHARACTER*80 WRNMSG ! Warning message for log file
C DATA STATEMENTS
       DATA KEYLST / '*hdr','*flg','*spc','*gas','*atm','*tan' /
       DATA GOTABS / .FALSE. /  ! 1  
       DATA GOTBBT / .FALSE. /  ! 2      
       DATA GOTCOO / .FALSE. /  ! 3      
       DATA GOTCRV / .FALSE. /  ! 4      
       DATA GOTFIN / .FALSE. /  ! 5      
       DATA GOTFOV / .FALSE. /  ! 6      
       DATA GOTGRD / .FALSE. /  ! 7      
       DATA GOTHIT / .FALSE. /  ! 8      
       DATA GOTILS / .FALSE. /  ! 9  
       DATA GOTJAC / .FALSE. /  ! 10  
       DATA GOTLEV / .FALSE. /  ! 11 
       DATA GOTLUT / .FALSE. /  ! 12  
       DATA GOTNTE / .FALSE. /  ! 13  
       DATA GOTOBS / .FALSE. /  ! 14  
       DATA GOTOPT / .FALSE. /  ! 15  
       DATA GOTPRF / .FALSE. /  ! 16  
       DATA GOTPTH / .FALSE. /  ! 17  
       DATA GOTRAD / .FALSE. /  ! 18  
       DATA GOTREJ / .FALSE. /  ! 19  
       DATA GOTRJT / .FALSE. /  ! 20  
       DATA GOTSFC / .FALSE. /  ! 21  
       DATA GOTSHP / .FALSE. /  ! 22  
       DATA GOTTAB / .FALSE. /  ! 23  
       DATA GOTTPS / .FALSE. /  ! 24  
       DATA GOTTRA / .FALSE. /  ! 25 
       DATA GOTWID / .FALSE. /  ! 26 
       DATA GOTXSC / .FALSE. /  ! 27 
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Open Driver file and print first record title to log file
C
      NAMVAR = NAMDRV
      CALL OPNFIL ( LUNDRV, NAMVAR, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
C Skip any subsequent comment records. First non-comment record should start 
C with *HDR in which case LDUMMY should be set TRUE.
C
      CALL NXTREC ( LUNDRV, RECORD, LDUMMY, FAIL, ERRMSG ) 
      IF ( FAIL ) RETURN
C
      DO IKEY = 1, MINKEY
        READ ( LUNDRV, '(A)', IOSTAT=IOS, ERR=900 ) RECORD
        CALL LOCASE ( RECORD(1:4), KEYCHK )
        IF ( KEYCHK .EQ. KEYLST(IKEY) ) THEN
          CALL RFMLOG ( RECORD, FAIL, ERRMSG )         ! Copy to LOG file
          IF ( FAIL ) RETURN
          IF ( IKEY .EQ. 1 ) THEN              
            CALL INPHDR ( LUNDRV, FAIL, ERRMSG )            
            IF ( FAIL ) RETURN
          ELSE IF ( IKEY .EQ. 2 ) THEN              
            CALL INPFLG ( LUNDRV, FAIL, ERRMSG )
            IF ( FAIL ) RETURN
          ELSE IF ( IKEY .EQ. 3 ) THEN 
            CALL INPSPC ( LUNDRV, LUNNXT, FAIL, ERRMSG )
            IF ( FAIL ) RETURN
          ELSE IF ( IKEY .EQ. 4 ) THEN
            CALL INPGAS ( LUNDRV, LUNNXT, GETHIT, GETXSC, FAIL, ERRMSG )
            IF ( FAIL ) RETURN
          ELSE IF ( IKEY .EQ. 5 ) THEN
            CALL INPATM ( LUNDRV, LUNNXT, FAIL, ERRMSG )
            IF ( FAIL ) RETURN
          ELSE IF ( IKEY .EQ. 6 ) THEN
            CALL INPTAN ( LUNDRV, LUNNXT, KEYCHK, FAIL, ERRMSG )
            IF ( FAIL ) RETURN
          END IF
C For *TAN section, valid alternatives are *DIM, *ELE, *GEO, *LEN, *LEV or *SEC
        ELSE IF ( IKEY .EQ. 6 .AND. 
     &            ( KEYCHK .EQ. '*dim' .OR. KEYCHK .EQ. '*ele' .OR.
     &              KEYCHK .EQ. '*geo' .OR. KEYCHK .EQ. '*len' .OR. 
     &              KEYCHK .EQ. '*lev' .OR. KEYCHK .EQ. '*sec' ) ) THEN
          CALL RFMLOG ( RECORD, FAIL, ERRMSG )         ! Copy to LOG file
          IF ( FAIL ) RETURN
          CALL INPTAN ( LUNDRV, LUNNXT, KEYCHK, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        ELSE 
          FAIL = .TRUE.
          ERRMSG = 'F-RFMINP: Unexpected Record in Drv.file:'//
     &              RECORD(1:37)//'...'
          RETURN
        END IF
      END DO
C
C Check if NTE data section is required
      GETNTE = NTEFLG .AND. NNTE .EQ. 0 ! might have been loaded in .atm files
C
C Check for optional records, which may be mandatory if GET* variable is set
C
  100 READ ( LUNDRV, '(A)', IOSTAT=IOS, ERR=900 ) RECORD
      IF ( IOS .EQ. -1 ) THEN                       ! End of file
        RECORD = 'W-RFMINP: rfm.drv not properly terminated '//
     &    'with *END <CR> record - end assumed' 
        CALL RFMLOG ( RECORD, FAIL, ERRMSG )
        RECORD = '*END'
      ELSE
        CALL RFMLOG ( RECORD, FAIL, ERRMSG )
      END IF
      IF ( FAIL ) RETURN
      CALL LOCASE ( RECORD(1:4), KEYCHK )
C Set up error and warning messages specific for this section, in case required
      ERRMSG = 'F-RFMINP: Duplicate '//RECORD(1:4)//
     &         ' section in driver file'
      WRNMSG = 'F-RFMINP: Ignoring '//RECORD(1:4)//
     &         ' section - option not enabled'
      FAIL = .TRUE.
C
      IF ( KEYCHK .EQ. '*abs' ) THEN
        IF ( GOTABS ) RETURN    ! with FAIL=TRUE and preset ERRMSG
        GOTABS = .TRUE.
        IF ( ABSFLG ) THEN
          CALL INPNAM ( LUNDRV, 'ABS', FAIL, ERRMSG )
        ELSE
          CALL SKPSEC ( LUNDRV, WRNMSG, FAIL, ERRMSG )
        END IF  
      ELSE IF ( KEYCHK .EQ. '*bbt' ) THEN
        IF ( GOTBBT ) RETURN    ! with FAIL=TRUE and preset ERRMSG
        GOTBBT = .TRUE.
        IF ( BBTFLG ) THEN
          CALL INPNAM ( LUNDRV, 'BBT', FAIL, ERRMSG )
        ELSE
          CALL SKPSEC ( LUNDRV, WRNMSG, FAIL, ERRMSG )
        END IF  
      ELSE IF ( KEYCHK .EQ. '*coo' ) THEN
        IF ( GOTCOO ) RETURN    ! with FAIL=TRUE and preset ERRMSG
        GOTCOO = .TRUE.
        IF ( COOFLG ) THEN
          CALL INPNAM ( LUNDRV, 'COO', FAIL, ERRMSG )
        ELSE
          CALL SKPSEC ( LUNDRV, WRNMSG, FAIL, ERRMSG )
        END IF  
      ELSE IF ( KEYCHK .EQ. '*crv' ) THEN
        IF ( GOTCRV ) RETURN    ! with FAIL=TRUE and preset ERRMSG
        GOTCRV = .TRUE.
        IF ( CRVFLG ) THEN
          CALL INPCRV ( LUNDRV, FAIL, ERRMSG )
        ELSE
          CALL SKPSEC ( LUNDRV, WRNMSG, FAIL, ERRMSG ) ! preset WRNMSG
        ENDIF
      ELSE IF ( KEYCHK .EQ. '*fin' ) THEN
        IF ( GOTFIN ) RETURN    
        GOTFIN = .TRUE.
        IF ( FINFLG ) THEN
          CALL INPFIN ( LUNDRV, FAIL, ERRMSG )
        ELSE
          CALL SKPSEC ( LUNDRV, WRNMSG, FAIL, ERRMSG ) 
        ENDIF
      ELSE IF ( KEYCHK .EQ. '*fov' ) THEN
        IF ( GOTFOV ) RETURN    
        GOTFOV = .TRUE.
        IF ( FOVFLG ) THEN
          CALL INPFOV ( LUNDRV, LUNNXT, FAIL, ERRMSG )
        ELSE
          CALL SKPSEC ( LUNDRV, WRNMSG, FAIL, ERRMSG ) 
        ENDIF
      ELSE IF ( KEYCHK .EQ. '*grd' ) THEN
        IF ( GOTGRD ) RETURN    
        GOTGRD = .TRUE.
        IF ( GRDFLG ) THEN
          CALL INPGRD ( LUNDRV, LUNNXT, FAIL, ERRMSG )
        ELSE
          CALL SKPSEC ( LUNDRV, WRNMSG, FAIL, ERRMSG ) 
        ENDIF
      ELSE IF ( KEYCHK .EQ. '*hit' ) THEN
        IF ( GOTHIT ) RETURN
        GOTHIT = .TRUE.
        IF ( GETHIT ) THEN
          CALL INPHIT ( LUNDRV, LUNNXT, FAIL, ERRMSG )
        ELSE
          CALL SKPSEC ( LUNDRV, WRNMSG, FAIL, ERRMSG )
        END IF
      ELSE IF ( KEYCHK .EQ. '*ils' ) THEN
        IF ( GOTILS ) RETURN
        GOTILS = .TRUE.
        IF ( ILSFLG ) THEN
          CALL INPILS ( LUNDRV, LUNNXT, FAIL, ERRMSG )
        ELSE
          CALL SKPSEC ( LUNDRV, WRNMSG, FAIL, ERRMSG )
        END IF
      ELSE IF ( KEYCHK .EQ. '*jac' ) THEN
        IF ( GOTJAC ) RETURN    
        GOTJAC = .TRUE.
        IF ( JACFLG ) THEN
          CALL INPJAC ( LUNDRV, LUNNXT, FAIL, ERRMSG )
        ELSE
          CALL SKPSEC ( LUNDRV, WRNMSG, FAIL, ERRMSG )
        END IF
      ELSE IF ( KEYCHK .EQ. '*lev' ) THEN
        IF ( GOTLEV ) RETURN    
        GOTLEV = .TRUE.
        IF ( LEVFLG ) THEN
          CALL INPLEV ( LUNDRV, LUNNXT, FAIL, ERRMSG )
        ELSE
          CALL SKPSEC ( LUNDRV, WRNMSG, FAIL, ERRMSG )
        END IF
      ELSE IF ( KEYCHK .EQ. '*lut' ) THEN
        IF ( GOTLUT ) RETURN    
        GOTLUT = .TRUE.
        IF ( LUTFLG ) THEN
          CALL INPLUT ( LUNDRV, LUNNXT, FAIL, ERRMSG )
        ELSE
          CALL SKPSEC ( LUNDRV, WRNMSG, FAIL, ERRMSG )
        END IF
      ELSE IF ( KEYCHK .EQ. '*nte' ) THEN
        IF ( GOTNTE ) RETURN    
        GOTNTE = .TRUE.
        IF ( NTEFLG ) THEN
          CALL INPNTE ( LUNDRV, LUNNXT, FAIL, ERRMSG )
        ELSE
          CALL SKPSEC ( LUNDRV, WRNMSG, FAIL, ERRMSG )
        END IF
      ELSE IF ( KEYCHK .EQ. '*obs' ) THEN
        IF ( GOTOBS ) RETURN    
        GOTOBS = .TRUE.
        IF ( OBSFLG ) THEN
          CALL INPOBS ( LUNDRV, FAIL, ERRMSG )
        ELSE
          CALL SKPSEC ( LUNDRV, WRNMSG, FAIL, ERRMSG )
        END IF
      ELSE IF ( KEYCHK .EQ. '*opt' ) THEN
        IF ( GOTOPT ) RETURN 
        GOTOPT = .TRUE.
        IF ( OPTFLG ) THEN
          CALL INPNAM ( LUNDRV, 'OPT', FAIL, ERRMSG )
        ELSE
          CALL SKPSEC ( LUNDRV, WRNMSG, FAIL, ERRMSG )
        END IF 
      ELSE IF ( KEYCHK .EQ. '*prf' ) THEN
        IF ( GOTPRF ) RETURN 
        GOTPRF = .TRUE.
        IF ( PRFFLG ) THEN
          CALL INPNAM ( LUNDRV, 'PRF', FAIL, ERRMSG )
        ELSE
          CALL SKPSEC ( LUNDRV, WRNMSG, FAIL, ERRMSG )
        END IF  
      ELSE IF ( KEYCHK .EQ. '*pth' ) THEN
        IF ( GOTPTH ) RETURN 
        GOTPTH = .TRUE.
        IF ( PTHFLG ) THEN
          CALL INPNAM ( LUNDRV, 'PTH', FAIL, ERRMSG )
        ELSE
          CALL SKPSEC ( LUNDRV, WRNMSG, FAIL, ERRMSG )
        END IF  
      ELSE IF ( KEYCHK .EQ. '*rad' ) THEN
        IF ( GOTRAD ) RETURN 
        GOTRAD = .TRUE.
        IF ( RADFLG ) THEN
          CALL INPNAM ( LUNDRV, 'RAD', FAIL, ERRMSG )
        ELSE
          CALL SKPSEC ( LUNDRV, WRNMSG, FAIL, ERRMSG )
        END IF  
      ELSE IF ( KEYCHK .EQ. '*rej' ) THEN
        IF ( GOTREJ ) RETURN 
        GOTREJ = .TRUE.
        IF ( REJFLG ) THEN
          CALL INPREJ ( LUNDRV, FAIL, ERRMSG )
        ELSE
          CALL SKPSEC ( LUNDRV, WRNMSG, FAIL, ERRMSG )
        END IF  
      ELSE IF ( KEYCHK .EQ. '*rjt' ) THEN
        IF ( GOTRJT ) RETURN 
        GOTRJT = .TRUE.
        IF ( RJTFLG ) THEN
          CALL INPNAM ( LUNDRV, 'RJT', FAIL, ERRMSG )
        ELSE
          CALL SKPSEC ( LUNDRV, WRNMSG, FAIL, ERRMSG )
        END IF  
      ELSE IF ( KEYCHK .EQ. '*sfc' ) THEN
        IF ( GOTSFC ) RETURN 
        GOTSFC = .TRUE.
        IF ( SFCFLG ) THEN
          CALL INPSFC ( LUNDRV, LUNNXT, FAIL, ERRMSG )
        ELSE
          CALL SKPSEC ( LUNDRV, WRNMSG, FAIL, ERRMSG )
        END IF  
      ELSE IF ( KEYCHK .EQ. '*shp' ) THEN
        IF ( GOTSHP ) RETURN 
        GOTSHP = .TRUE.
        IF ( SHPFLG ) THEN
          CALL INPSHP ( LUNDRV, FAIL, ERRMSG )
        ELSE
          CALL SKPSEC ( LUNDRV, WRNMSG, FAIL, ERRMSG )
        END IF  
      ELSE IF ( KEYCHK .EQ. '*tab' ) THEN
        IF ( GOTTAB ) RETURN 
        GOTTAB = .TRUE.
        IF ( TABFLG ) THEN
          CALL INPNAM ( LUNDRV, 'TAB', FAIL, ERRMSG )
        ELSE
          CALL SKPSEC ( LUNDRV, WRNMSG, FAIL, ERRMSG )
        END IF  
      ELSE IF ( KEYCHK .EQ. '*tps' ) THEN
        IF ( GOTTPS ) RETURN 
        GOTTPS = .TRUE.
        IF ( TPSFLG ) THEN
          CALL INPTPS ( LUNDRV, LUNNXT, FAIL, ERRMSG )
        ELSE
          CALL SKPSEC ( LUNDRV, WRNMSG, FAIL, ERRMSG )
        END IF  
      ELSE IF ( KEYCHK .EQ. '*tra' ) THEN
        IF ( GOTTRA ) RETURN 
        GOTTRA = .TRUE.
        IF ( TRAFLG ) THEN
          CALL INPNAM ( LUNDRV, 'TRA', FAIL, ERRMSG )
        ELSE
          CALL SKPSEC ( LUNDRV, WRNMSG, FAIL, ERRMSG )
        END IF  
      ELSE IF ( KEYCHK .EQ. '*wid' ) THEN
        IF ( GOTWID ) RETURN 
        GOTWID = .TRUE.
        IF ( WIDFLG ) THEN
          CALL INPNAM ( LUNDRV, 'WID', FAIL, ERRMSG )
        ELSE
          CALL SKPSEC ( LUNDRV, WRNMSG, FAIL, ERRMSG )
        END IF  
      ELSE IF ( KEYCHK .EQ. '*xsc' ) THEN
        IF ( GOTXSC ) RETURN
        GOTXSC = .TRUE.
        IF ( GETXSC ) THEN
          CALL INPXSC ( LUNDRV, LUNNXT, FAIL, ERRMSG )
        ELSE
          CALL SKPSEC ( LUNDRV, WRNMSG, FAIL, ERRMSG )
        END IF
      ELSE IF ( KEYCHK .EQ. '*end' ) THEN
C
C End of Driver File reached. Check all required optional files read in.
C 
        CLOSE ( LUNDRV, IOSTAT=IOS, ERR=900 )
        FAIL = .TRUE.
        IF ( CRVFLG .AND. .NOT. GOTCRV ) THEN
          ERRMSG = 'F-RFMINP: No *CRV data supplied'
        ELSE IF ( FINFLG .AND. .NOT. GOTFIN ) THEN
          ERRMSG = 'F-RFMINP: No *FIN data supplied'
        ELSE IF ( FOVFLG .AND. .NOT. GOTFOV ) THEN
          ERRMSG = 'F-RFMINP: No *FOV data supplied'
        ELSE IF ( GRDFLG .AND. .NOT. GOTGRD ) THEN
          ERRMSG = 'F-RFMINP: No *GRD data supplied'
        ELSE IF ( GETHIT .AND. .NOT. GOTHIT .AND. .NOT. LUTFLG ) THEN
          ERRMSG = 'F-RFMINP: No *HIT data supplied'
        ELSE IF ( ILSFLG .AND. .NOT. GOTILS ) THEN
          ERRMSG = 'F-RFMINP: No *ILS data supplied'
        ELSE IF ( JACFLG .AND. .NOT. GOTJAC ) THEN
          ERRMSG = 'F-RFMINP: No *JAC data supplied'
        ELSE IF ( LEVFLG .AND. .NOT. GOTLEV ) THEN
          ERRMSG = 'F-RFMINP: No *LEV data supplied'
        ELSE IF ( LUTFLG .AND. .NOT. GOTLUT ) THEN
          ERRMSG = 'F-RFMINP: No *LUT data supplied'
        ELSE IF ( GETNTE .AND. .NOT. GOTNTE ) THEN
          ERRMSG = 'F-RFMINP: No *NTE data supplied'
        ELSE IF ( OBSFLG .AND. .NOT. GOTOBS ) THEN
          ERRMSG = 'F-RFMINP: No *OBS data supplied'
        ELSE IF ( REJFLG .AND. .NOT. GOTREJ ) THEN
          ERRMSG = 'F-RFMINP: No *REJ data supplied'
        ELSE IF ( SFCFLG .AND. .NOT. GOTSFC ) THEN
          ERRMSG = 'F-RFMINP: No *SFC data supplied'
        ELSE IF ( SHPFLG .AND. .NOT. GOTSHP ) THEN
          ERRMSG = 'F-RFMINP: No *SHP data supplied'
        ELSE IF ( TPSFLG .AND. .NOT. GOTTPS ) THEN
          ERRMSG = 'F-RFMINP: No *TPS data supplied'
        ELSE IF ( GETXSC .AND. .NOT. GOTXSC .AND. .NOT. LUTFLG ) THEN
          ERRMSG = 'F-RFMINP: No *XSC data supplied'
        ELSE
          CALL INPCHK ( FAIL, ERRMSG )    ! Sets FAIL=FALSE for exit if all OK
        END IF
        RETURN
      ELSE
        FAIL = .TRUE.
        IF ( KEYCHK .EQ. '*hdr' .OR. KEYCHK .EQ. '*flg' .OR.
     &       KEYCHK .EQ. '*spc' .OR. KEYCHK .EQ. '*gas' .OR.
     &       KEYCHK .EQ. '*atm' .OR. KEYCHK .EQ. '*tan' .OR.
     &       KEYCHK .EQ. '*dim' .OR. KEYCHK .EQ. '*ele' .OR.
     &       KEYCHK .EQ. '*geo' .OR. KEYCHK .EQ. '*len' .OR.
     &       KEYCHK .EQ. '*sec'                              ) THEN 
          ERRMSG = 'F-RFMINP: Duplicate Primary Section:'//KEYCHK
        ELSE 
          ERRMSG = 'F-RFMINP: Unrecognised Key:'//RECORD(1:20)
        END IF
        RETURN
      END IF
      IF ( FAIL ) RETURN
      GOTO 100                      ! Read next record
C
C Exit with Fatal I/O Error 
C 
  900 WRITE ( ERRMSG, '(A,I11)' ) 
     &  'F-RFMINP: I/O Failure on Driver File. IOSTAT=', IOS
      FAIL = .TRUE.
      END
