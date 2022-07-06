      SUBROUTINE RFMXSC ( FAIL, ERRMSG )
C
C VERSION
C     08-AUG-13  AD  New molecule indices - check for old values still in use
C                    No temperature extrapolation unless V42 flag enabled.
C     01-DEC-10  AD  Keep separate MIN/MAX WNO for each table.
C     19-AUG-02  AD  SAVE everything
C     17-DEC-01  AD  Initialise & Save NXSP=0 to avoid warning messges
C     17-APR-01  AD  Move LXSC from local to to xsccom.inc 
C     29-DEC-01  AD  Adapted to cope with HITRAN 2000 format
C     07-SEP-99  AD  Allow for p,T tabulation instead of just T
C     09-AUG-99  AD  Change DELXSW from Real to D.P.
C     10-JUN-99  AD  Check that XSC data not superseded by LUT.   
C     03-MAR-97  AD  Version 3.
C     28-DEC-96  AD  Use QTINT rather than extrapolation beyond tabul.Temps.
C     27-DEC-96  AD  Correction: Set IXST when interpolating Temp X/S data.
C     01-OCT-96  AD  Version 2.
C     28-SEP-96  AD  Add Aerosol - avoid calling QTINT
C     01-SEP-96  AD  Version 1.
C 
C DESCRIPTION    
C     Calculate the absorption due to cross-section species 
C     Called by RFM for each wide mesh interval.
C     First read in external cross-section data for the heavy molecules. 
C     The cross-section data is interpolated onto the fine wavenumber grid and 
C     the absorption calculated for the relevant paths by triangulation on the
C     (irregular) tabulation p,T grid.
C     The cross-section data files are assumed to be set out: 
C                 
C                 *** MOL NXST              
C                 (A3,I7,I10)
C                 CHARID WNOMIN WNOMAX NXSW TEMXST PREXST HEADER
C                 (A10,F10.4,F10.4,I10,F10.1,F10.4,A) 
C                               or  ...F10.4,F10.1...
C                 XSCTAB(I),I=1,NXSW (plus offset for each pT table)
C                 (1P10E10.3)    
C                 
C                 ***    EACH SECTION STARTS WITH THIS FLAG
C                 MOL    IS THE MOLECULAR ID (>32). FOR THE SAME 
C                        MOLECULE THERE MAY BE MORE THAN ONE XSEC SET 
C                        FOR EACH DIFFERENT BAND. EACH BAND HAS A 
C                        SEPARATE HEADER.
C                 NXST   NO OF PRESSURE,TEMPERATURE TABLES
C                 CHARID MOLECULE CHARACTER IDENTIFIER (CHEMICAL NAME)
C                 WNOMIN WAVENUMBER OF FIRST CROSS-SECTION DATA POINT
C                 WNOMAX WAVENUMBER OF LAST CROSS-SECTION DATA POINT
C                 NXSW   NUMBER OF CROSS-SECTION DATA POINTS
C                 TEMXST TEMPERATURE AT WHICH DATA WERE TAKEN (K)
C                 PREXST PRESSURE AT WHICH DATA WERE TAKEN (TORR)
C                 HEADER WHICH INCLUDES INFORMATION ON THE
C                        LABORATORY AT WHICH DATA WAS TAKEN AND
C                        MOLECULE COMMON NAME
C                 XSCTAB  CROSS-SECTION DATA IN UNITS (CM2/MOLEC)
C
C                  
      IMPLICIT NONE
      SAVE
C
C ARGUMENTS      
      LOGICAL      FAIL    !  O  T=a fatal error has been detected
      CHARACTER*80 ERRMSG  !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  IDGNEW ! Convert old RFM index for .xsc data to new value
     &, QTINT  ! Calculate total internal partition sums for tabulated data.   
     &, TRIANG ! Interpolation of irregular 2D grid by triangulation.
      INTEGER IDGNEW
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
      INCLUDE 'shpcon.inc' ! Line-shape codes
C
C COMMON VARIABLES
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'idxcon.inc' ! HITRAN molecular indices
      INCLUDE 'pthcom.inc' ! RFM path data
      INCLUDE 'xflcom.inc' ! X-section data files
      INCLUDE 'xsccom.inc' ! X-section data 
C
C LOCAL CONSTANTS
      REAL TEMTOL          ! Temp. tolerance [K] before applying QTINT factor
        PARAMETER ( TEMTOL = 0.1 ) 
C
C LOCAL VARIABLES
      INTEGER  IDGOLD    ! Old RFM index for some gases (for V42FLG)
      INTEGER  IDXTRI(3,MAXCLC) ! Indices for triangulation points 
      INTEGER  IFIN      ! Fine mesh grid point counter (0:NFIN)
      INTEGER  IGAS      ! Gas counter
      INTEGER  IOFXST(MAXXST) ! Offset index for each (p,T) dataset
      INTEGER  IOS       ! Saved value of IOSTAT for error messages
      INTEGER  IPTH      ! Path counter
      INTEGER  IREC      ! Record counter
      INTEGER  IXSC      ! Counter for X/S data sets
      INTEGER  IXSP,JXSP ! Index for X/S data points
      INTEGER  IXST      ! Counter for X/S data temperature tabulation points
      INTEGER  IXFL      ! Index of X/S file
      INTEGER  LUN       ! Logical Unit Number for X/S file
      INTEGER  MOL       ! X/S Molecule ID code (51:63), 50=aerosol
      INTEGER  NXSP      ! No. of X/S data points
      INTEGER  NXST      ! No. of p,T tabulations
      INTEGER  NXSW      ! No. of wavenumber tabulations
      REAL     AMT       ! Gas amount in current path [ moles ? ]
      REAL     DXSP      ! Fracs. of Waveno.tabul.intvl above FM Grid Pt.
      REAL     PRELIM    ! Interpolated pressre [Torr] from table
      REAL     PREXST(MAXXST)        ! List of X/S data Pres.tab. pts [Torr]
      REAL     PTORR     ! Pressure [Torr] (760 Torr = 1 atm)
      REAL     SQ        ! Ratio of Tot.Partit.Sum at 296K to sum at path temp.
      REAL     SQLIM     ! Ratio of TIPS at 296K to TIPS at tabulation limit
      REAL     TEM       ! Temperature of current path
      REAL     TEMFAC(MAXCLC) ! Temperature extrapolation factor
      REAL     TEMLIM    ! Interpolated temperature from table
      REAL     TEMXST(MAXXST)        ! List of X/S data temp.tab. pts [K]
      REAL     WGTTRI(3,MAXCLC)      ! Triangulation weights for each path
      REAL     XSC       ! Value of interpolated X/S data to path,FM Grid pt.
      REAL     XSCFIN(MAXXST) ! X/S data tabulation interp.to FM grid.
      REAL     XSCTAB(MAXXSP) ! X/S data tabulation
      REAL     XXSP      ! Position of FM Grid pt in wavenumber tabulation axis
      LOGICAL  H2KFMT    ! T=HITRAN 2000 format, F=Earlier format
      LOGICAL  INIT      ! T=initialise triangulation routine
      CHARACTER*1      CDUMMY ! Dummy read character
      CHARACTER*10     FMT    ! File Format identifier
      CHARACTER*80     RECORD ! Text record read from X/S file
      DOUBLE PRECISION DELXSW(MAXXST) ! Intvl [/cm] in X/S waveno. tabulation
      DOUBLE PRECISION WMNXST(MAXXST) ! Min.Wno for each X/S tabulation
      DOUBLE PRECISION WMXXST(MAXXST) ! Max.Wno for each X/S tabulation
      DOUBLE PRECISION WNO    ! Waveno of current Fine Mesh grid point
      DOUBLE PRECISION WNOMIN,WNOMAX ! Min,Max Wavenumbers for X/S data
C
C DATA STATEMENTS
      DATA NXSP / 0 / 
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Loop over cross-section data file
C
      DO IXSC = 1, NXSC
        IGAS = IGSXSC(IXSC)
        IF ( SHPGAS(IGAS) .EQ. SHPXSC .AND.          ! Not replaced by LUT
     &       WN2XSC(IXSC) .GE. WN1FIN .AND.
     &       WN1XSC(IXSC) .LE. WN2FIN ) THEN
          IF ( IXSC .NE. LXSC ) THEN                 ! Load data from file
            IXFL = IFLXSC(IXSC)
            LUN = LUNXFL(IXFL)            
            REWIND ( LUN, IOSTAT=IOS, ERR=900 )
            DO IREC = 1, IRCXSC(IXSC)-1
              READ ( LUN, '(A1)', IOSTAT=IOS, ERR=900 ) CDUMMY
            END DO
            READ ( LUN, '(3X,I7,I10,A10)', IOSTAT=IOS, ERR=900 ) 
     &        MOL, NXST, FMT
            H2KFMT = FMT .NE. ' '
C Convert old xsc molecular index to new, if required
            MOL = IDGNEW ( MOL ) 
C For V42FLG also require old molecule ID for QTINT for 5 specific molecules
            IF ( V42FLG ) THEN
              IDGOLD = MOL
              IF ( MOL .EQ. 111 ) IDGOLD = 51    ! f11
              IF ( MOL .EQ. 122 ) IDGOLD = 56    ! f22
              IF ( MOL .EQ. 104 ) IDGOLD = 60    ! ccl4
              IF ( MOL .EQ. 101 ) IDGOLD = 61    ! clono2
              IF ( MOL .EQ. 111 ) IDGOLD = 62    ! n2o5
            END IF
            IF (MOL .NE. IDXGAS(IGAS)) STOP 'F-RFMXSC: Logical Error#1'
            NXSP = 0
            DO IXST = 1, NXST
C Due to the different formats for expressing tabulated temperature and 
C pressure, and the file record beginning with a un-quoted character string, 
C it is necessary to read the data into a character string, and then perform a
C free-format read on the part of the string known to contain the appropriate 
C information.
              READ ( LUN, '(A)', IOSTAT=IOS, ERR=900 ) RECORD
              IF ( H2KFMT ) THEN    ! Ignore characters 1-20
                READ ( RECORD(21:60), *, IOSTAT=IOS, ERR=900 ) 
     &            WNOMIN, WNOMAX, NXSW, TEMXST(IXST), PREXST(IXST)
              ELSE                  ! Ignore characters 1-10
                READ ( RECORD(11:60), *, IOSTAT=IOS, ERR=900 ) 
     &            WNOMIN, WNOMAX, NXSW, TEMXST(IXST), PREXST(IXST)
              END IF
              DELXSW(IXST) = ( WNOMAX - WNOMIN ) / ( NXSW - 1 ) 
              IOFXST(IXST) = NXSP 
              WMNXST(IXST) = WNOMIN
              WMXXST(IXST) = WNOMAX
              READ ( LUN, *, ERR=900, IOSTAT=IOS )
     &          ( XSCTAB(IXSP), IXSP = NXSP + 1, NXSP + NXSW )
              NXSP = NXSP + NXSW
            END DO
            INIT = .TRUE.                             ! Set to FALSE by TRIANG
            DO IPTH = 1, NCLC  
              IF ( IGSPTH(IPTH) .EQ. IGAS ) THEN
                TEM = TEMPTH(IPTH)
C NB: Triangulation assumes temperature and pressure axes are appropriately
C scaled, so keep P in Torr so variation in each axis is a few hundred units.
                PTORR = PREPTH(IPTH) * 760.0          ! Convert atm to Torr
                CALL TRIANG ( TEM, PTORR, NXST, TEMXST, PREXST, INIT,
     &                        IDXTRI(1,IPTH), WGTTRI(1,IPTH),  
     &                        TEMLIM, PRELIM, FAIL, ERRMSG  )
                IF ( FAIL ) RETURN ! Shouldn't happen since tested earlier
C For old (v4.2) only ..
C If actual temperature (TEM) is beyond tabulation scale (TEMLIM), 
C extrapolate by applying partition function ratio SQ appropriate for each 
C molecule (not for aerosol, though)
                IF ( V42FLG .AND.             
     &               MOL .NE. IDXAER .AND.       
     &               ABS ( TEM - TEMLIM ) .GT. TEMTOL ) THEN
                  CALL QTINT ( IDGOLD, TEMLIM, SQLIM )
                  CALL QTINT ( IDGOLD, TEM,    SQ    )
                  TEMFAC(IPTH) = SQ / SQLIM
                ELSE                     ! no T extrapolation in new RFM
                  TEMFAC(IPTH) = 1.0
                END IF
              END IF
            END DO
            LXSC = IXSC
          END IF
C
C Interpolate the cross-section data onto the fine mesh
C
          IXSFIN(IGAS) = 0         ! info for WID diagnostics
          DO IFIN = 1, NFIN
            WNO = WNOFIN(IFIN) 
            DO IXST = 1, NXST
              IF ( WNO .LT. WMNXST(IXST) .OR. 
     &             WNO .GE. WMXXST(IXST)      ) THEN
                XSCFIN(IXST) = 0.0
              ELSE
                IXSFIN(IGAS) = 1
                XXSP = 1.0 + SNGL ( ( WNO - WMNXST(IXST) ) 
     &                                      / DELXSW(IXST) )
                IXSP = INT ( XXSP ) + IOFXST(IXST)
                JXSP = MIN ( 1 + IXSP, NXSP )  ! (if WNO=WNOMAX) limit to NXST
                DXSP = MOD ( XXSP, 1.0 )
                XSCFIN(IXST) = XSCTAB(IXSP) * ( 1.0 - DXSP ) + 
     &                         XSCTAB(JXSP) * DXSP
              END IF
            END DO
C 
            DO IPTH = 1, NCLC                                 ! Loop over paths
              IF ( IGAS .EQ. IGSPTH(IPTH) ) THEN
                AMT = AVOG * AMTPTH(IPTH)
C
C Interpolate between temperature data sets, calculate absorption
C
                XSC = WGTTRI(1,IPTH) * XSCFIN(IDXTRI(1,IPTH)) +
     &                WGTTRI(2,IPTH) * XSCFIN(IDXTRI(2,IPTH)) +
     &                WGTTRI(3,IPTH) * XSCFIN(IDXTRI(3,IPTH))
                ABSFIN(IFIN,IPTH) = ABSFIN(IFIN,IPTH) + 
     &                              XSC * AMT * TEMFAC(IPTH)
                CNTFIN(IFIN,IPTH) = ABSFIN(IFIN,IPTH)
              END IF
            END DO
          END DO
        END IF
      END DO
C 
      FAIL = .FALSE.
      RETURN
C
  900 FAIL = .TRUE. 
      WRITE ( ERRMSG, '(A,I3,A,I11)' ) 
     &  'F-RFMXSC: I/O failure on X/S file#', IXFL, ' IOSTAT=', IOS
C
      END
