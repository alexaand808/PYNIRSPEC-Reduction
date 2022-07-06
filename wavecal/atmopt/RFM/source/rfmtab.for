      SUBROUTINE RFMTAB ( FAIL, ERRMSG )
C
C VERSION
C     07-JUN-00  AD  Allow for irregular grids
C     09-AUG-99  AD  Make XV = SNGL( ) 
C     04-JUN-99  AD  Allow for Binary files as well
C     04-MAY-99  AD  Allow for NP or NT = 1 in look-up tables
C     28-APR-99  AD  Allow for -lnp,Tem axis increment scaling factors
C                    Also allow for temperature offset profiles
C     19-JUN-98  AD  Correction: ensure buffer entries separated correctly 
C                    when two or more TAB.LUT files are used.
C     07-JAN-98  AD  Use ARGMIN from PHYCON.INC instead of 1.0E-38
C     01-NOV-97  AD  Original.
C 
C DESCRIPTION    
C     Calculate the absorption using external Abs.Coeff. Tables.
C     Called by RFM for each wide mesh interval.
C     Note that different tables are stored in the same two buffers using
C     different IX values, so that only part of these buffers are written with
C     new data from each record.
C                  
      IMPLICIT NONE
C
C ARGUMENTS
      LOGICAL       FAIL      !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG    !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  TABAXS ! Resample (p,T) axis from TAB input file.
     &, TEMOFF ! Temperature offset for pretabulated (p,T) data.
      REAL TEMOFF
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'shpcon.inc' ! Line-shape codes
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
C
C COMMON VARIABLES
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'pthcom.inc' ! RFM path data
      INCLUDE 'lutcom.inc' ! Look-Up Table data
C
C LOCAL VARIABLES
      INTEGER  IBUF      ! Buffer pointer for upper wavenumber
      INTEGER  IFIN      ! Fine mesh grid point counter (0:NFIN)
      INTEGER  IG        ! Wavenumber record index for current TAB file
      INTEGER  IGAS      ! Gas counter
      INTEGER  IGL       ! Wavenumber record index in IGLUT array
      INTEGER  II,IJ,JI,JJ ! Indices of 2-d interpolation points in K matrix
      INTEGER  IOS       ! Saved value of IOSTAT for error message
      INTEGER  IP        ! Index of lower interp. point in -lnp tabulation axis
      INTEGER  IPTH      ! Path counter
      INTEGER  ITI, ITJ  ! Index of lower interp. point in Tem tabulation axis
      INTEGER  ILUT      ! Counter for LUT data sets
      INTEGER  IX        ! Counter for p,T points
      INTEGER  IX0,IX1,IXN ! Counter for 0th, 1st and Nth elements in NX axis
      INTEGER  LBUF      ! Buffer pointer for lower wavenumber
      INTEGER  LUN       ! LUN of Abs.Coeff file
      INTEGER  NP        ! Local value of number of points in -lnp axis
      INTEGER  NT        ! Local value of number of points in Temp axis
      REAL     ABS       ! Single gas optical depth for path segment 
      REAL     DP        ! Local value of -lnp axis increment [p in mb]
      REAL     DT        ! Local value of Temp axis increment [K]
      REAL     EP        ! Fraction of -lnp Tab.interval above path value
      REAL     ETI,ETJ   ! Fraction of Temp Tab.interval above path value
      REAL     EV        ! Fraction of Wno  Tab.interval above Fine Mesh value
      REAL  KVHI, KVLO   ! Interp.Abs.Coeffs at high, low wavenumbers
      REAL     LNP       ! -lnp [p in mb] of path segment
      REAL     P1        ! Local value of 1st point in -lnp axis [p in mb]
      REAL  TABBUF(MAXLUX,2) ! Storage for 2 consec.wno records of TAB file
      REAL     T1        ! Local value of 1st point in Temp axis [K]
      REAL     TEM       ! Temperature of current path
      REAL  WII,WIJ,WJI,WJJ ! Weights of 2-d interpolation points in K matrix
      REAL     XP        ! Position of path -lnp in tabulation axis
      REAL     XT        ! Position of path temp in tabulation axis
      DOUBLE PRECISION WNO ! Waveno of current Fine Mesh grid point
      DOUBLE PRECISION WNOLOW ! Waveno of lower buffer contents
      DOUBLE PRECISION WNOUPP ! Waveno of upper buffer contents
      LOGICAL  BINFIL    ! T=binary file, F=ASCII file
      LOGICAL  USEAXS    ! T=Use resampled (p,T) axes
      LOGICAL  USELIN    ! T=Tabulated function is Abs.Coeff. k
      LOGICAL  USE4RT    ! T=Tabulated function is SQRT(SQRT(k))
C
C EXECUTABLE CODE --------------------------------------------------------------
C
C Loop over all LUT for this spectral range
      DO ILUT = 1, NLUT
        IF ( LUNLUT(ILUT) .NE. 0 ) THEN      ! Flag for Abs.Coeff. table
          LUN = LUNLUT(ILUT)
          BINFIL = BINLUT(ILUT)
          IX0 = IXOFF(ILUT)
          IX1 = IX0 + 1
          IXN = IX0 + NPLUT(ILUT) * NTLUT(ILUT)
          NP = NPLUT(ILUT)
          P1 = P1LUT(ILUT)
          DP = DPLUT(ILUT)
          NT = NTLUT(ILUT)
          T1 = T1LUT(ILUT)
          DT = DTLUT(ILUT)
          USELIN = TABLUT(ILUT) .EQ. 'lin'
          USE4RT = TABLUT(ILUT) .EQ. '4rt'
          USEAXS = NDPLUT(ILUT) .NE. 1 .OR. NDTLUT(ILUT) .NE. 1
C
C Set up wavenumbers for contents of upper and lower buffers for this LUT
C For new spectral range, IRCLUT=0 
          IG = IRCLUT(ILUT)
          IGL = MAX ( IG, 2 ) + IGOFF(ILUT) ! MAX ensures IGL valid if IG=0
          IBUF = 1 + MOD ( IG, 2 ) 
          LBUF = 3 - IBUF
          WNOLOW = V1LUT(ILUT) + ( IGLUT(IGL-1) - 1 ) * DVLUT(ILUT)
          WNOUPP = V1LUT(ILUT) + ( IGLUT(IGL) - 1 ) * DVLUT(ILUT)
C
C Loop over fine mesh grid points
          DO IFIN = 1, NFIN
            WNO = WNOFIN(IFIN)
C Advance buffers until WNOUPP larger than fine mesh point or two recs.loaded
            DO WHILE ( ( WNO .GT. WNOUPP .AND. IG .LT. NGLUT(ILUT) ) 
     &                 .OR. IG .LT. 2 )
              IG = IG + 1
              IRCLUT(ILUT) = IG
              LBUF = IBUF
              IBUF = 3 - IBUF
              IGL = IG + IGOFF(ILUT)
              WNOLOW = WNOUPP
              WNOUPP = V1LUT(ILUT) + ( IGLUT(IGL) - 1 ) * DVLUT(ILUT)
              IF ( BINFIL ) THEN
                READ ( LUN, IOSTAT=IOS, ERR=900 ) 
     &            ( TABBUF(IX,IBUF), IX = IX1, IXN )
              ELSE 
                READ ( LUN, *, IOSTAT=IOS, ERR=900 ) 
     &            ( TABBUF(IX,IBUF), IX = IX1, IXN )
              END IF
C Convert contents of TABBUF to ln(k)
              IF ( USE4RT ) THEN
                DO IX = IX1, IXN
                  TABBUF(IX,IBUF) = 
     &              4.0 * LOG ( MAX ( ARGMIN, TABBUF(IX,IBUF) ) )
                END DO
              ELSE IF ( USELIN ) THEN
                DO IX = IX1, IXN
                  TABBUF(IX,IBUF) = 
     &              LOG ( MAX ( ARGMIN, TABBUF(IX,IBUF) ) )
                END DO
              END IF
C Undersample p,T axes if required
              IF ( USEAXS ) CALL TABAXS ( ILUT, 
     &          NP, P1, DP, NT, T1, DT, TABBUF(IX1,IBUF) )
            END DO
C Now have two wavenumbers loaded into buffers.
C If these bracket the fine grid wavenumber, then interpolate 
            EV = SNGL ( ( WNO - WNOLOW ) / ( WNOUPP - WNOLOW ) )
            IF ( EV .GE. 0.0 .AND. EV .LE. 1.0 ) THEN
              DO IPTH = 1, NCLC
                IGAS = IGSPTH(IPTH) 
                IF ( ILUT .EQ. SHPGAS(IGAS) - SHPLUT ) THEN
                  LNP = - LOG ( ATMB * PREPTH(IPTH) )
                  TEM = TEMPTH(IPTH)
C Interpolation points and weights in -lnp axis
                  XP = ( LNP - P1 ) / DP + 1.0
                  XP = MIN ( MAX ( 1.0, XP ), FLOAT(NP) )        !limit 1:NP
                  IP = MIN ( INT ( XP ), NP-1 )                  !limit 1:NP-1
                  EP = XP - FLOAT ( IP )
C Interpolation points and weights in Temperature axis
                  IF ( OFFLUT(ILUT) ) THEN
                    LNP = P1 + (IP-1) * DP                       ! -lnp at IP
                    TEM = TEMPTH(IPTH) - TEMOFF(LNP)
                  END IF
                  XT = ( TEM - T1 ) / DT + 1.0
                  XT = MIN ( MAX ( 1.0, XT ), FLOAT(NT) )        !limit 1:NT
                  ITI = MIN ( INT ( XT ), NT-1 )                 !limit 1:NT-1
                  ETI = XT - FLOAT ( ITI )
                  IF ( OFFLUT(ILUT) ) THEN
                    LNP = LNP + DP                               !-lnp at IP+1
                    TEM = TEMPTH(IPTH) - TEMOFF(LNP)
                    XT = ( TEM - T1 ) / DT + 1.0
                    XT = MIN ( MAX ( 1.0, XT ), FLOAT(NT) )      !limit 1:NT
                    ITJ = MIN ( INT ( XT ), NT-1 )               !limit 1:NT-1
                    ETJ = XT - FLOAT ( ITJ )
                  ELSE
                    ITJ = ITI
                    ETJ = ETI
                  END IF
C Indices and weights in 1:NX dimension for (p,T) interpolation
                  II = IX0 + IP + NP * ( ITI- 1 )        ! low -lnp, low  tem
                  JI = II + 1 + ( ITJ - ITI ) * NP       ! high -lnp, low  tem
                  IJ = II + NP                           ! low  -lnp, high tem
                  JJ = IJ + 1                            ! high -lnp, high tem
C Might end up with II, IJ or JI (not JJ) .LE.0 if NP or NT=1, in which case
C the associated weight W** should also be zero, but reset index=1 to ensure
C no problems accessing array elements
                  II = MAX ( 1, II )
                  IJ = MAX ( 1, IJ )
                  JI = MAX ( 1, JI )
                  WII = ( 1.0 - EP ) * ( 1.0 - ETI )
                  WJI = EP * ( 1.0 - ETJ )
                  WIJ = ( 1.0 - EP ) * ETI
                  WJJ = EP * ETJ
C
                  KVLO = WII*TABBUF(II,LBUF) + WIJ*TABBUF(IJ,LBUF) 
     &                 + WJI*TABBUF(JI,LBUF) + WJJ*TABBUF(JJ,LBUF)
                  KVHI = WII*TABBUF(II,IBUF) + WIJ*TABBUF(IJ,IBUF) 
     &                 + WJI*TABBUF(JI,IBUF) + WJJ*TABBUF(JJ,IBUF)
C
C Interpolate absorption coeff array to fine mesh and scale by absorber amounts
C for single-gas optical depths. Note that ABSFIN may contain continuum data
C already, so add on this term.
C Factor 1.0e7 to convert AMT [kmoles/cm^2]) to [moles/m^2] for compatibility
C with pretabulated K coefficients.
C
                  ABS = EXP ( (1.0 - EV) * KVLO + EV * KVHI ) * 
     &                  AMTPTH(IPTH) * 1.0E7
                  ABSFIN(IFIN,IPTH) = ABSFIN(IFIN,IPTH) + ABS  
                  CNTFIN(IFIN,IPTH) = ABSFIN(IFIN,IPTH)
                END IF
              END DO
            END IF
          END DO
        END IF
      END DO
      FAIL = .FALSE.
      RETURN
C
  900 WRITE ( ERRMSG, '(A,I11)')
     &  'F-RFMTAB: Input failure on TAB file. IOSTAT=', IOS
      FAIL = .TRUE.
C
      END
