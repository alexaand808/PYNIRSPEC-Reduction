      SUBROUTINE FILPRF ( LUNATM, DIRECT, NLEV, IPRF, FAIL, ERRMSG )
C
C VERSION
C     08OCT13 AD Bug#82: Add CHKPRF to check input profiles
C     23JUL09 AD Correction: Add VIBATM as local array
C     08JAN08 AD Allow for Vib.Tem profiles
C     11NOV03 AD Local flag to determine aerosol interpolation
C     20JAN03 AD Set CHGATM = TRUE on reading first profile
C     12DEC01 AD Save HGTLEV
C     06JUN01 AD Check heights increase monotonically.
C                    Remove call to PRFCHK
C     04NOV99 AD Set RETURN if FAIL=TRUE before PRFCHK
C     21JUN99 AD Add PRFCHK subroutine
C     24JAN99 AD Replace OFMFLG by LINFLG
C     25JUL97 AD Increase MAXLEV from 121 to 201
C     03MAR97 AD Version 3.
C     13JAN97 AD Correction to error message generated if NATM>MAXATM
C                    Increase MAXLEV from 100 to 121
C     01OCT96 AD Version 2.
C     01SEP96 AD Version 1.
C
C DESCRIPTION
C     Read atmospheric profile from file into /ATMCOM/.
C     Called by ATMFIL for each requested profile in file.
C     If IPRF .ge. 100000 then assume this is a Vib.Tem-KT profile.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER       LUNATM  !  I  LUN for Atmospheric profiles 
      LOGICAL       DIRECT  !  I  T=Load directly, F=interpolate to alts.
      INTEGER       NLEV    !  I  No. of profile levels to use
      INTEGER       IPRF    !  I  Code identifying type of profile
      LOGICAL       FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  ATMNTE ! Load vibrational temperature profile from .atm file
     &, CHKPRF ! Check atmospheric profile on input.
     &, INTERP ! General interpolation of arrays
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C LOCAL CONSTANTS
      INTEGER       MAXLEV         ! Max.no levels in supplied profiles
        PARAMETER ( MAXLEV = 201 ) ! Should be comparable or .GT. MAXATM
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data 
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
C
C LOCAL VARIABLES
      INTEGER  IATM            ! Surface counter for used profile (1:NATM)
      INTEGER  ILEV            ! Surface counter for supplied profile (1:NLEV)
      INTEGER  IOS             ! Saved value of IOSTAT for error message
      LOGICAL  AERLIN          ! T=linearly interpolate aerosol,F=use LINFLG
      REAL     HGTLEV(MAXLEV)  ! Altitudes of supplied profiles
      REAL     PRFLEV(MAXLEV)  ! Data values of supplied profiles
      REAL     VIBATM(MAXATM)  ! Vib.Tem profile interpolated to .atm levels
      SAVE HGTLEV
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Set AERLIN = true if aerosol always to be interpolated linearly with alt
C Set AERLIN = false if handled the same way as VMR (ie log unless LINFLG set)
      AERLIN = .TRUE.
C
C If this is the first file read (DIRECT=TRUE) load directly into /ATMCOM/
C
      IF ( DIRECT ) THEN        
        IF ( IPRF .EQ. -2 ) THEN                          ! Altitude [km]
          NATM = NLEV
          FAIL = .TRUE.
          IF ( NATM .GT. MAXATM ) THEN
            WRITE ( ERRMSG, '(A,I11,A,I6)' )
     &        'F-FILPRF: No.profile levs=', NATM,
     &        ' > RFMSIZ.INC dimension MAXATM =', MAXATM
            RETURN
          ELSE IF ( NATM .EQ. 1 .AND. .NOT. HOMFLG ) THEN
            ERRMSG = 'F-FILPRF: Single atmos.level specified but '//
     &               'RFM not flagged for Homog. path calc.'
            RETURN
          ELSE IF ( HOMFLG .AND. NATM .NE. 1 ) THEN
            ERRMSG = 'F-FILPRF: RFM set for Homog. path calc. but'//
     &               ' multi-layer atmos.profile supplied'
            RETURN
          ELSE
            READ ( LUNATM, *, IOSTAT=IOS, ERR=900 ) 
     &          ( HGTATM(IATM), IATM = 1, NATM )
            CHGATM(1) = .TRUE.
            CALL CHKPRF ( 'HGT', NATM, HGTATM, FAIL, ERRMSG )
            IF ( FAIL ) RETURN
            DO IATM = 2, NATM
              CHGATM(IATM) = .TRUE.
            END DO
          END IF
        ELSE IF ( IPRF .EQ. -1 ) THEN                     ! Pressure [mb]
          READ ( LUNATM, *, IOSTAT=IOS, ERR=900 ) 
     &        ( PREATM(IATM), IATM = 1, NATM )
          CALL CHKPRF ( 'PRE', NATM, PREATM, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        ELSE IF ( IPRF .EQ. 0 ) THEN                      ! Temperature [K]
          READ ( LUNATM, *, IOSTAT=IOS, ERR=900 ) 
     &        ( TEMATM(IATM), IATM = 1, NATM )
          CALL CHKPRF ( 'TEM', NATM, TEMATM, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        ELSE IF ( IPRF .GE. 1000000 ) THEN                ! VT-KT [K]
          READ ( LUNATM, *, IOSTAT=IOS, ERR=900 ) 
     &        ( VIBATM(IATM), IATM = 1, NATM )            ! +ve or -ve, no check
          CALL ATMNTE ( IPRF, VIBATM, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        ELSE                                               ! VMR [ppmv]
          READ ( LUNATM, *, IOSTAT=IOS, ERR=900 )
     &        ( VMRATM(IATM,IPRF), IATM = 1, NATM )       
          CALL CHKPRF ( CODGAS(IPRF), NATM, VMRATM(1,IPRF), FAIL,ERRMSG)
          IF ( FAIL ) RETURN
          DO IATM = 1, NATM
            VMRATM(IATM,IPRF) = VMRATM(IATM,IPRF)*1.0E-6  ! convert ppmv to ppv
          END DO
        END IF
C
C If not first file, then interpolate from altitudes HGTLEV to alts HGTATM
C
      ELSE
        IF ( IPRF .EQ. -2 ) THEN                 ! Set up altitudes for file
          IF ( NLEV .GT. MAXLEV ) THEN
            FAIL = .TRUE.
            WRITE ( ERRMSG, '(A,I11,A,I6)' )
     &        'F-FILPRF: No.profile levs=', NLEV,
     &        ' > Local dimension MAXLEV=', MAXLEV
            RETURN
          END IF
          READ ( LUNATM, *, IOSTAT=IOS, ERR=900 ) 
     &        ( HGTLEV(ILEV), ILEV = 1, NLEV )
          CALL CHKPRF ( 'HGT', NLEV, HGTLEV, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        ELSE IF ( IPRF .EQ. -1 ) THEN            
          READ ( LUNATM, *, IOSTAT=IOS, ERR=900 )          ! Read pressure
     &      ( PRFLEV(ILEV), ILEV = 1, NLEV )
          CALL CHKPRF ( 'PRE', NLEV, PRFLEV, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
          CALL INTERP ( NLEV, NATM, HGTLEV, HGTATM,        ! Log interpolation
     &                  PRFLEV, .TRUE., .TRUE., PREATM )   ! plus extrapolation
        ELSE IF ( IPRF .EQ. 0 ) THEN              
          READ ( LUNATM, *, IOSTAT=IOS, ERR=900 )          ! Read temperature 
     &      ( PRFLEV(ILEV), ILEV = 1, NLEV )
          CALL CHKPRF ( 'TEM', NLEV, PRFLEV, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
          CALL INTERP ( NLEV, NATM, HGTLEV, HGTATM,        ! Linear interp
     &                  PRFLEV, .FALSE., .FALSE., TEMATM ) ! no extrapolation
        ELSE IF ( IPRF .GE. 1000000 ) THEN                             
          READ ( LUNATM, *, IOSTAT=IOS, ERR=900 )          ! Read Vib. Temp.
     &      ( PRFLEV(ILEV), ILEV = 1, NLEV )
          CALL INTERP ( NLEV, NATM, HGTLEV, HGTATM,        ! Linear interp
     &                  PRFLEV, .FALSE., .FALSE., VIBATM ) ! no extrapolation
          CALL ATMNTE ( IPRF, VIBATM, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        ELSE 
          READ ( LUNATM, *, IOSTAT=IOS, ERR=900 )          ! Read vmr
     &      ( PRFLEV(ILEV), ILEV = 1, NLEV )
          CALL CHKPRF ( CODGAS(IPRF), NLEV, PRFLEV, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
          IF ( AERLIN .AND. IPRF .EQ. IAXGAS ) THEN        ! Lin.Interp aerosol
            CALL INTERP ( NLEV, NATM, HGTLEV, HGTATM, PRFLEV, 
     &                    .FALSE., .FALSE., VMRATM(1,IPRF) ) ! No Extrap.
          ELSE
            CALL INTERP ( NLEV, NATM, HGTLEV, HGTATM, PRFLEV, ! Log/Lin Interp
     &                  (.NOT. LINFLG), .FALSE., VMRATM(1,IPRF) ) ! No Extrap.
          END IF
          DO IATM = 1, NATM
            VMRATM(IATM,IPRF) = VMRATM(IATM,IPRF)*1.0E-6  ! convert ppmv to ppv
          END DO
        END IF               
      END IF               
C
      FAIL = .FALSE.
      RETURN
C
  900 WRITE ( ERRMSG, '(A,I11)' )
     &  'F-FILPRF: I/O failure reading Atm.Profile file. IOSTAT=', IOS
      FAIL = .TRUE.
C
      END
