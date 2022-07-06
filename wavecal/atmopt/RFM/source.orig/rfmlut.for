      SUBROUTINE RFMLUT 
C
C VERSION
C     08-JUN-00  AD  Allow for LUTs on irregular grids
C     09-AUG-99  AD  Use explicit XV = SNGL( ) 
C     04-MAY-99  AD  Allow for NP or NT = 1 in look-up tables
C     28-APR-99  AD  Allow use of offset temperature profile
C     07-JAN-98  AD  Use ARGMIN from PHYCON.INC instead of KMIN.
C     07-NOV-97  AD  Check for -ve k before taking ln(k) for interpolation.
C     01-NOV-97  AD  Interpolate for (p,T) in Log(K) whatever tabulation used
C                    Also check for use of external Abs.Coeff tables
C     18-JUL-97  AD  Test for tabulation function TABLUT
C     08-JUL-97  AD  correction to table-edge handling: XP, XT
C     02-JUL-97  AD  Original.
C 
C DESCRIPTION    
C     Calculate the absorption using Look-Up Tables.
C     Called by RFM for each wide mesh interval.
C                  
      IMPLICIT NONE
C
      EXTERNAL 
     &  TEMOFF ! Temperature offset for pretabulated (p,T) data.
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
      INTEGER  IFIN      ! Fine mesh grid point counter (0:NFIN)
      INTEGER  IGAS      ! Gas counter
      INTEGER  II,IJ,JI,JJ ! Indices of 2-d interpolation points in K matrix
      INTEGER  IL        ! Counter for basis vectors
      INTEGER  IP        ! Index of lower interp. point in -lnp tabulation axis
      INTEGER  IPTH      ! Path counter
      INTEGER  ITI,ITJ   ! Indices of lower interp. point in Tem tabulation axis
      INTEGER  ILUT      ! Counter for LUT data sets
      INTEGER  IV        ! Index for Tab.grid for particular LUT
      INTEGER  IVL       ! Absolute index for IV in IVLUT array 
      INTEGER  IVLMAX    ! Index for IVLUT array element corresp. to WNOMAX
      INTEGER  IVLMIN    ! Index for IVLUT array element corresp. to WNOMIN
      INTEGER  IVMAX(MAXLUT) ! Tab.grid pointers correspoding to WNOMAX
      INTEGER  IVMIN(MAXLUT) ! Tab.grid pointers correspoding to WNOMIN
      INTEGER  IVUMTX    ! Index for U-matrix wavenumber axis
      REAL     ABS       ! Single gas optical depth for path segment 
      REAL     DP        ! Fraction of -lnp Tab.interval above path value
      REAL     DTI,DTJ   ! Fraction of Temp Tab.interval above path value
      REAL     DV        ! Fraction of Wno  Tab.interval above Fine Mesh value
      REAL     KABS(MAXLUV) ! Reconstructed LUT Absorption Coefficient.
      REAL  KII,KIJ,KJI,KJJ ! Reconstructed LUT tabulated function
      REAL     LNP       ! -lnp [p in mb] of path segment
      REAL     TEM       ! Temperature of current path
      REAL  WII,WIJ,WJI,WJJ ! Weights of 2-d interpolation points in K matrix
      REAL     XP        ! Position of path -lnp in tabulation axis
      REAL     XT        ! Position of path temp in tabulation axis
      DOUBLE PRECISION WNO ! Waveno of current Fine Mesh grid point
      DOUBLE PRECISION WNO2ND ! Next tabulated wno above WNOMIN
      DOUBLE PRECISION WNOLOW ! Tabulated wno immediately below fine grid point
      DOUBLE PRECISION WNOMAX ! Upper tabulated wno spanning WM interval
      DOUBLE PRECISION WNOMIN ! Lower tabulated wno spanning WM interval
      DOUBLE PRECISION WNOUPP ! Tabulated wno immediately above fine grid point
      LOGICAL  USELOG    ! T=Tabulated function is ln(k)
      LOGICAL  USE4RT    ! T=Tabulated function is SQRT(SQRT(k))
C
C DATA STATEMENTS
      DATA IVMIN / MAXLUT * 1 /
      DATA IVMAX / MAXLUT * 1 /
      SAVE IVMIN, IVMAX
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      DO IPTH = 1, NCLC
        IGAS = IGSPTH(IPTH) 
        IF ( SHPGAS(IGAS) .GT. SHPLUT ) THEN        ! Marked for use of LUT
          ILUT = SHPGAS(IGAS) - SHPLUT
          IF ( LUNLUT(ILUT) .EQ. 0 ) THEN           ! Use SVD-compressed table
            USELOG = TABLUT(ILUT) .EQ. 'log'
            USE4RT = TABLUT(ILUT) .EQ. '4rt'
            LNP = - LOG ( ATMB * PREPTH(IPTH) )
            TEM = TEMPTH(IPTH) 
C Interpolation points and weights in -lnp axis
            XP = ( LNP - P1LUT(ILUT) ) / DPLUT(ILUT) + 1.0
            XP = MIN ( MAX ( 1.0, XP ), FLOAT(NPLUT(ILUT)) )  ! limit to 1:NP
            IP = MIN ( INT ( XP ), NPLUT(ILUT)-1 )            ! limit to 1:NP-1
            DP = XP - FLOAT ( IP )
C Interpolation points and weights in Temperature axis
            IF ( OFFLUT(ILUT) ) THEN
              LNP = P1LUT(ILUT) + (IP-1) * DPLUT(ILUT)        ! -lnp at IP
              TEM = TEMPTH(IPTH) - TEMOFF(LNP)
            END IF
            XT = ( TEM - T1LUT(ILUT) ) / DTLUT(ILUT) + 1.0
            XT = MIN ( MAX ( 1.0, XT ), FLOAT(NTLUT(ILUT)) )  ! limit to 1:NT
            ITI = MIN ( INT ( XT ), NTLUT(ILUT)-1 )           ! limit to 1:NT-1
            DTI = XT - FLOAT ( ITI )
            IF ( OFFLUT(ILUT) ) THEN
              LNP = LNP + DPLUT(ILUT)
              TEM = TEMPTH(IPTH) - TEMOFF(LNP)
              XT = ( TEM - T1LUT(ILUT) ) / DTLUT(ILUT) + 1.0
              XT = MIN ( MAX ( 1.0, XT ), FLOAT(NTLUT(ILUT)) ) ! limit to 1:NT
              ITJ = MIN ( INT ( XT ), NTLUT(ILUT)-1 )         ! limit to 1:NT-1
              DTJ = XT - FLOAT ( ITJ ) 
            ELSE
              ITJ = ITI
              DTJ = DTI
            END IF
C Indices and weights in 1:NX dimension for (p,T) interpolation
            II = IXOFF(ILUT) + IP + NPLUT(ILUT) * ( ITI - 1 ) ! lo -lnp, lo tem
            JI = II + 1 + ( ITJ - ITI ) * NPLUT(ILUT)         ! hi -lnp, lo tem
            IJ = II + NPLUT(ILUT)                             ! lo -lnp, hi tem
            JJ = IJ + 1 + ( ITJ - ITI ) * NPLUT(ILUT)         ! hi -lnp, hi tem
C Might end up with II, IJ or JI (not JJ) .LE.0 if NP or NT=1, in which case
C the associated weight W** should also be zero, but reset index=1 to ensure
C no problems accessing array elements
            II = MAX ( 1, II )
            IJ = MAX ( 1, IJ )
            JI = MAX ( 1, JI )
            WII = ( 1.0 - DP ) * ( 1.0 - DTI )
            WJI = DP * ( 1.0 - DTJ )
            WIJ = ( 1.0 - DP ) * DTI
            WJJ = DP * DTJ
C
C Multiply out basis vectors to construct Absorption Coefficient array
C
C Find index range of LUT wavenumber axes that span this interval
C In the following, the tabulated points are WNOMIN,WNO2ND,...,WNOMAX
C and the aim is to set the tabulation pointers so that the end points of the 
C current WM interval (WN1FIN:WN2FIN) are each bracketed by the end two tabulated
C points: WNOMIN < WN1FIN < WNO2ND and (penultimate pt) < WN2FIN < WNOMAX 
C
C Use previously saved value for WNOMIN for this LUT
            IVLMIN = IVMIN(ILUT) + IVOFF(ILUT)
            WNOMIN = V1LUT(ILUT) + ( IVLUT(IVLMIN) - 1 ) * DVLUT(ILUT)
C
C If WM intvl starts below current lower tabulated point, this could just be 
C because the spectral range is still below the lowest tabulated point, or it
C could be the start of a new spectral range. Either way: reset pointer to
C lowest tabulated points.
            IF ( WNOMIN .GT. WN1FIN ) THEN  
              IVMIN(ILUT) = 1
              IVLMIN = 1 + IVOFF(ILUT)
              WNOMIN = V1LUT(ILUT) 
              IVMAX(ILUT) = 1
              IVLMAX = 1 + IVOFF(ILUT)
              WNOMAX = V1LUT(ILUT) 
            END IF
C
C Advance WNOMIN,WNO2ND pointers until WN1FIN < WNO2ND
            WNO2ND = V1LUT(ILUT) + ( IVLUT(IVLMIN+1)-1 ) * DVLUT(ILUT)
            DO WHILE ( WNO2ND .LT. WN1FIN .AND. 
     &                 IVMIN(ILUT) .LT. NVLUT(ILUT)-1 )
              IVMIN(ILUT) = IVMIN(ILUT) + 1
              IVLMIN = IVLMIN + 1
              WNOMIN = WNO2ND
              WNO2ND = V1LUT(ILUT) + ( IVLUT(IVLMIN+1)-1 ) * DVLUT(ILUT)
            END DO
C
C Adance WNOMAX pointer until WN2FIN < WNOMAX
            IVLMAX = IVMAX(ILUT) + IVOFF(ILUT)
            WNOMAX = V1LUT(ILUT) + ( IVLUT(IVLMAX) - 1 ) * DVLUT(ILUT)
            DO WHILE ( WNOMAX .LT. WN2FIN .AND. 
     &                 IVMAX(ILUT) .LT. NVLUT(ILUT) )
              IVMAX(ILUT) = IVMAX(ILUT) + 1
              IVLMAX = IVLMAX + 1
              WNOMAX= V1LUT(ILUT) + ( IVLUT(IVLMAX) - 1 ) * DVLUT(ILUT)
            END DO
C
C Expand basis vectors to fill KABS array for required wavenumber range to span
C this widemesh interval
            DO IV = IVMIN(ILUT), IVMAX(ILUT)
              IVUMTX = IV + IVOFF(ILUT)
              KABS(IV) = 0.0                           ! Initialise to zero 
              KII = 0.0
              KIJ = 0.0
              KJI = 0.0
              KJJ = 0.0
              DO IL = 1, NLLUT(ILUT)                ! loop over all basis vector
                KII = KII + U(IVUMTX,IL) * K(IL,II) ! multiply U * K(p,T) at the
                KIJ = KIJ + U(IVUMTX,IL) * K(IL,IJ) ! four points required for
                KJI = KJI + U(IVUMTX,IL) * K(IL,JI) ! interpolation
                KJJ = KJJ + U(IVUMTX,IL) * K(IL,JJ) 
              END DO
              IF ( USELOG ) THEN
                KABS(IV) = WII*KII + WIJ*KIJ + WJI*KJI + WJJ*KJJ 
              ELSE
                KABS(IV) = WII * LOG ( MAX ( KII, ARGMIN ) ) + 
     &                     WIJ * LOG ( MAX ( KIJ, ARGMIN ) ) + 
     &                     WJI * LOG ( MAX ( KJI, ARGMIN ) ) +
     &                     WJJ * LOG ( MAX ( KJJ, ARGMIN ) ) 
                IF ( USE4RT ) KABS(IV) = 4.0*KABS(IV)
              END IF
            END DO
C
C Interpolate absorption coeff array to fine mesh and scale by absorber amounts
C for single-gas optical depths. Note that ABSFIN may contain continuum data
C already, so add on this term.
C Factor 1.0e7 to convert AMT [kmoles/cm^2]) to [moles/m^2] for compatibility
C with pretabulated K coefficients.
C
            IVL   = IVLMIN
            IV    = IVLMIN - IVOFF(ILUT)
            WNOLOW = WNOMIN
            WNOUPP = WNO2ND
            DO IFIN = 1, NFIN
              WNO = WNOFIN(IFIN)
              DO WHILE ( WNO .GT. WNOUPP .AND. IV .LT. NVLUT(ILUT)-1 ) 
                IVL   = IVL + 1
                IV    = IV + 1
                WNOLOW = WNOUPP
                WNOUPP = V1LUT(ILUT) + ( IVLUT(IVL+1)-1 ) * DVLUT(ILUT)
              END DO
              DV = SNGL ( ( WNO - WNOLOW ) / ( WNOUPP - WNOLOW ) )
              IF ( DV .GE. 0.0 .AND. DV .LE. 1.0 ) THEN
                ABS = EXP ( ( 1.0 - DV ) * KABS(IV) + DV * KABS(IV+1) )  
     &                * AMTPTH(IPTH) * 1.0E7
                ABSFIN(IFIN,IPTH) = ABSFIN(IFIN,IPTH) + ABS  
                CNTFIN(IFIN,IPTH) = ABSFIN(IFIN,IPTH)
              END IF
            END DO
          END IF
        END IF
      END DO
C
      END
