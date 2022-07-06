      SUBROUTINE TPSFIL ( LUNTPS, NAMTPS, FAIL, ERRMSG )
C
C VERSION
C     17OCT13 AD Change to DO WHILE ... structure
C     13FEB03 AD Original.
C
C DESCRIPTION
C     Open TIPS data file, check and load contents.
C     Called by INPTPS for each field in *TPS section of driver table.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER       LUNTPS  !  I  Temporary LUN for TPS File
      CHARACTER*(*) NAMTPS  !  I  Name of TPS file
      LOGICAL       FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  NXTREC ! Load next record from RFM input file.
     &, OPNFIL ! Open ASCII input file, and skip/Log any initial comments.
     &, RFMLOG ! Write text message to RFM log file.
C
C GLOBAL CONSTANTS
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES 
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'tpscom.inc' ! Tabulated TIPS data
C
C LOCAL VARIABLES
      INTEGER      IOS      ! Saved value of IOSTAT for error messages
      INTEGER      IDGAS    ! HITRAN Index of molecule read from file
      INTEGER      IGAS     ! Index of molecule in *GAS variables
      INTEGER      ISO      ! Isotope number read from file (HITRAN numbering)
      INTEGER      IPT      ! Counter for tabulated data points
      INTEGER      ITPS     ! Index of molec/isotope in *TPS variables
      INTEGER      JTPS     ! Counter for molec/isotope in *TPS variables
      INTEGER      NPT      ! No. temperature tabulation points
      REAL         DUMMY    ! Dummy variable for skipping ignored data
      REAL         PT1      ! First temperature tabulation point [K]
      REAL         PT2      ! Last temperature tabulation point [K]
      REAL         PTD      ! Temperature tabulation increment [K]
      REAL         TPSREF   ! TIPS function evaluated at Ref.Temp.
      REAL         TPSTEM(MAXTPT) ! Values read directly from file
      REAL         XPT      ! Position of Ref.Temp in temperature axis intvl.
      REAL         XREF     ! Position of Ref.Temp in temperature axis
      LOGICAL      ENDSEC   ! Set TRUE when end of file reached
      CHARACTER*80 RECORD   ! Text record read from driver table
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      CALL OPNFIL ( LUNTPS, NAMTPS, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
C Loop here for start of each molecule/isotope entry
      CALL NXTREC ( LUNTPS, RECORD, ENDSEC, FAIL, ERRMSG ) 
      IF ( FAIL ) RETURN
      DO WHILE ( .NOT. ENDSEC ) 
        READ ( RECORD, *, IOSTAT=IOS, ERR=900 ) 
     &    IDGAS, ISO, NPT, PT1, PTD
        IF ( IDGAS .LT. 1 .OR. IDGAS .GT. MAXHLN ) THEN
          FAIL = .TRUE.
          WRITE ( ERRMSG, '(A,I11)' )
     &      'F-TPSFIL: Illegal value of IDGAS=', IDGAS
          RETURN
        END IF
        IF ( IGSMOL(IDGAS) .GT. 0 ) THEN       ! This molecule is used, so load
          IGAS = IGSMOL(IDGAS)
          IF ( ISO .LT. 1 .OR. ISO .GT. NISGAS(IGAS) ) THEN
            FAIL = .TRUE.
            WRITE ( ERRMSG, '(A,I3,A,I11)' )
     &        'F-TPSFIL: Gas ID=', IDGAS, ' has illegal Isotope#=', ISO
            RETURN
          END IF
          IF ( NPT .GT. MAXTPT ) THEN
            FAIL = .TRUE.
            WRITE ( ERRMSG, '(A,A,A,I3,A,I2,A,I11,A,I9)' )
     &        'F-TPSFIL: Gas ', CODGAS(IGAS), ' ID#', IDGAS, ' Iso#', 
     &        ISO, ' tabulated for ', NPT, ' Tems, >MAXTPT=', MAXTPT
            RETURN
          END IF   
C Check if data already loaded
          DO JTPS = 1, NTPS
            IF ( IDGAS .EQ. IDGTPS(JTPS) .AND. 
     &             ISO .EQ. ISOTPS(JTPS)       ) THEN
              WRITE ( RECORD, '(A,A,A,I3,A,I2,A)' )
     &          'W-TPSFIL: TIPS data for ', CODGAS(IGAS), ' Gas ID#', 
     &          IDGAS, ' Iso#', ISO, ' superseded'
              CALL RFMLOG ( RECORD, FAIL, ERRMSG )
              IF ( FAIL ) RETURN
              ITPS = JTPS
              GOTO 100
            END IF
          END DO
C If not, check space for additional data
          IF ( NTPS .EQ. MAXTPS ) THEN
            FAIL = .TRUE.
            WRITE ( ERRMSG, '(A,I11)' ) 
     &        'F-TPSFIL: No. TIPS tabulations > MAXTPS =', MAXTPS
            RETURN
          END IF
          WRITE ( RECORD, '(A,A,A,I3,A,I2)' )
     &      'I-TPSFIL: Loading TIPS data for ', CODGAS(IGAS), 
     &      ' Gas ID#', IDGAS, ' Iso#', ISO
          CALL RFMLOG ( RECORD, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
          NTPS = NTPS + 1
          ITPS = NTPS
  100     CONTINUE
C Load/Reload TIPS data
          READ ( LUNTPS, *, IOSTAT=IOS, ERR=900 ) 
     &      ( TPSTEM(IPT), IPT = 1, NPT )
C Find value corresponding to Reference Temperature (296K)
          PT2 = PT1 + ( NPT - 1 ) * PTD
          IF ( ( TEMREF .LE. PT1 .AND. TEMREF .LE. PT2 ) .OR.
     &         ( TEMREF .GE. PT1 .AND. TEMREF .GE. PT2 )      ) THEN
            FAIL = .TRUE.
            WRITE ( ERRMSG, '(A,F6.1)' )
     &        'F-TPSFIL: Tabulation does not span T.Ref=', TEMREF
            RETURN
          END IF
          XREF = 1.0 + ( TEMREF - PT1 ) / PTD
          IF ( XREF .LT. 1.0 .OR. XREF .GT. FLOAT ( NPT ) ) 
     &      STOP 'F-TPSFIL: Logical Error'
          IPT = INT ( XREF )
          XPT = XREF - FLOAT ( IPT ) 
          TPSREF = TPSTEM(IPT) * ( 1.0 - XPT ) + TPSTEM(IPT+1) * XPT
C Insert Q(296)/Q(T) into QFCTPS array
          DO IPT = 1, NPT
            QFCTPS(IPT,ITPS) = TPSREF / TPSTEM(IPT)
          END DO
          IDGTPS(ITPS) = IDGAS
          ISOTPS(ITPS) = ISO
          NPTTPS(ITPS) = NPT
          PT1TPS(ITPS) = PT1
          PTDTPS(ITPS) = PTD 
          IDXTPS(ISO,IDGAS) = ITPS
        ELSE
          READ ( LUNTPS, *, IOSTAT=IOS, ERR=900 ) (DUMMY, IPT = 1, NPT)
        END IF
C
C Get next molecule/isotope data from file
        CALL NXTREC ( LUNTPS, RECORD, ENDSEC, FAIL, ERRMSG ) 
        IF ( FAIL ) RETURN
      END DO
C
C Normal exit
      CLOSE ( LUNTPS, ERR=900, IOSTAT=IOS )
      FAIL = .FALSE.
      RETURN
C
C Exit with fatal I/O errors
 900  CONTINUE
      WRITE ( ERRMSG, '(A,I11)' ) 
     &  'F-TPSFIL: I/O error reading TIPS file, IOSTAT=', IOS
      FAIL = .TRUE.
C
      END
