      SUBROUTINE RFMNAD 
C
C VERSION
C     20-JAN-11  AD  Bug#77 Test DLPMIN to avoid instabilities for thin layers
C     15-JAN-04  AD  Use IDIR=0 in argument to IDXPTH
C     27-AUG-03  AD  Add comments on how to remove altitude scale dependence
C                    and correct definitions of units of USUM
C     16-AUG-01  AD  Rewrite, allowing for LIN flag
C     28-APR-00  AD  Set PSIPTH for ITAN=2,3 etc
C     18-APR-00  AD  Allow for FLX. Correct PSIPTH.
C     04-JAN-00  AD  Set PSIPTH.
C                    Add extra argument to IDXPTH
C     11-AUG-99  AD  Allow for special values of x.
C     17-DEC-98  AD  Bug fix: allow for VMR=0 when taking Log(VMR)
C     10-DEC-98  AD  Bug fix: calculate SECANG=1/SIN instead of 1/COS
C     29-JUL-98  AD  Get angle from USRTAN instead of HGTTAN
C     14-JUL-98  AD  Comment change only
C     03-MAR-97  AD  Version 3.
C     15-JAN-97  AD  Add zenith as well as nadir viewing
C                    Interpret 'HGTTAN' as airmass factor.
C     25-OCT-96  AD  Revised to work in pressure coordinates (not height)
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Calculate nadir/zenith path parameters through atmosphere.
C     Called by RFM if NAD or ZEN flags enabled.
C     Calculate CG quantities for each path = absorber & layer
C
C     Use notation Dx = x_l - x_u where x_l,x_u are x evaluated at the lower
C     and upper boundaries respectively.
C     v = vmr, p = pressure, T = temperature, e=partial pressure, u=abs.amount
C     Assumes T varies linearly with altitude in all cases.
C
C     LIN flag disabled: 
C       Assumes v varies within a layer as v = v_l(p/p_l)^b
C       where b=ln(v_u/v_l)/ln(p_u/p_l)
C          u =  1/(gM(1+b)) Dpv         
C         pu =  1/(gM(2+b)) Dp^2v      
C         eu = 1/(2gM(1+b)) Dp^2v^2    
C         Tu = 1/(gM(1+b)) ( DpTv - DTDpv/((1+b)Dlnp) ) 
C       special cases 
C         b=-1,   u = (1/gM)(v_l p_l Dlnp)     
C                eu = (1/gM)(v_l^2 p_l^2 Dlnp )
C                Tu = (1/gM)(v_l p_l (T_u+T_l)/2 Dlnp )
C         b=-2,  pu = (1/gM)(v_l p_l^2 Dlnp )
C     LIN flag enabled
C       Assumes v varies within a layer as v = v_l + a * (z - z_l)
C          u = (1/gM)(Dpv - DpDv/Dlnp)
C         pu = (1/gM)(Dp^2v - Dp^2Dv/(2Dlnp))
C         eu = (1/gM)(Dp^2v^2 - Dp^2vDv/Dlnp + Dp^2(Dv)^2/(2Dlnp)^2) )
C         Tu = (1/gM)(DpTv - (DvDpT+DTDpv)/(Dlnp) + 2DpDTDv/(Dlnp)^2 )
C
C     To avoid making assumptions about g or M, for hydrostatic equilbm assume
C                     M.g = R <T> ln(P2/P1) / ( Z2-Z1 )
C     However, to make the calculation independent of altitude scale the value
C     of M.g can be explicitly inserted (see comments inserted in code below
C     for how and where to do this)
C
C     It is assumed that the CG pressure, part.press and temperature do not
C     vary with airmass (ie zenith angle) so calculated for first angle listed
C     in *TAN section only.
C
      IMPLICIT NONE
C
      EXTERNAL
     &  IDXPTH ! Indexing function for Pth# from Tan.pt#, Atm.Layer# and Gas#
      INTEGER IDXPTH
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
C
C LOCAL CONSTANTS
      REAL     DLPMIN      ! Minimum value for dlnp for explicit differencing
        PARAMETER ( DLPMIN = 1.0E-3 ) ! 0.1%, about 7m min layer thickness
C Note: thin layers can cause problems in evaluating DPV etc, so avoid
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data  
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'tancom.inc' ! Tangent heights
      INCLUDE 'pthcom.inc' ! RFM path data
C
C LOCAL VARIABLES
      INTEGER  IATM       ! Atmos. layer# for path
      INTEGER  IDIR       ! Direction (away from obs): +1=zen, -1=nad
      INTEGER  IGAS       ! Gas# for path
      INTEGER  IPTH       ! Path counter
      INTEGER  ITAN       ! Tangent path# for path
      INTEGER  JATM       ! Index of profile point at top of segment
      INTEGER  JPTH       ! Secondary path counter
      REAL     B          ! VMR profile parameter: VMR = VMR1*(P/P1)^X
      REAL     DLNP       ! ln(p_l/p_u)
      REAL     DP         ! p_l - p_u
      REAL     DP2        ! p_l^2 - p_u^2
      REAL     DP2V       ! p_l^2 v_l - p_u^2 v_u
      REAL     DP2V2      ! p_l^2 v_l^2 - p_u^2 v_u^2
      REAL     DPT        ! p_l T_l - p_u T_u 
      REAL     DPTV       ! p_l T_l v_l - p_u T_u v_u
      REAL     DPV        ! p_l v_l - p_u v_l
      REAL     DT         ! T_l - T_u
      REAL     DV         ! v_l - v_u
      REAL     ECG        ! Curtis-Godson partial pressure [mb]
      REAL     HGT1, HGT2 ! Hgts [km] at lower,upper segment bounds
      REAL     PCG        ! Curtis-Godson pressure [mb]
      REAL     PRE1, PRE2 ! Pressure [mb] at upp,low segment bounds
      REAL     SECANG     ! Sec(zen.ang), = airmass factor
      REAL     TAVG       ! Mid-point temperature [K] of segment 
      REAL     TCG        ! Curtis-Godson temperature [K]
      REAL     TEM1, TEM2 ! Temperatures [K] at upp,low seg.bounds
      REAL     USUM       ! Absorber Amount in layer [10^5*kmole/m^2]
      REAL     VMR1, VMR2 ! Vmr [ppv] at upp,low seg.bounds
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      IF ( ZENFLG ) THEN
        IDIR = 1
      ELSE
        IDIR = -1
      END IF
C
C Loop over all gases
      DO IGAS = 1, NGAS
        HGT2 = HGTATM(1)
        VMR2 = VMRATM(1,IGAS)
        PRE2 = PREATM(1)
        TEM2 = TEMATM(1)
        DO IATM = 1, NATM-1
          JATM = IATM + 1
          IPTH = IDXPTH ( 1, IATM, IGAS, 0 )        ! 1 = ITAN, 0 = IDIR
C Copy previous segment upper value to current segment lower values
          VMR1 = VMR2
          HGT1 = HGT2
          PRE1 = PRE2
          TEM1 = TEM2
C Set current segment upper values from profile layer values
          VMR2 = VMRATM(JATM,IGAS) 
          HGT2 = HGTATM(JATM)
          PRE2 = PREATM(JATM)
          TEM2 = TEMATM(JATM) 
C Calculate summation terms for CG integrals (allow for VMR1=0 and/or VMR2=0)
          DLNP  = LOG(PRE1/PRE2)
          DP2V  = PRE1*PRE1*VMR1 - PRE2*PRE2*VMR2
          DP2V2 = PRE1*PRE1*VMR1*VMR1 - PRE2*PRE2*VMR2*VMR2
          DPTV  = PRE1*TEM1*VMR1 - PRE2*TEM2*VMR2
          DPV   = PRE1*VMR1 - PRE2*VMR2
          DT    = TEM1 - TEM2
          TAVG = 0.5 * ( TEM1 + TEM2 )
          IF ( LINFLG .OR. VMR1 .LT. ARGMIN 
     &                .OR. VMR2 .LT. ARGMIN ) THEN
            DP    = PRE1 - PRE2
            DP2   = PRE1*PRE1 - PRE2*PRE2
            DPT   = PRE1*TEM1 - PRE2*TEM2
            DV    = VMR1 - VMR2
            USUM  = DPV - DP * DV / DLNP
            IF ( USUM .NE. 0.0 .AND. DLNP .GE. DLPMIN ) THEN
              PCG = 0.5 * ( DP2V - 0.5*DP2*DV/DLNP ) / USUM
              ECG = 0.5 * ( DP2V2 - 
     &              ( DP2V*DV - 0.5*DP2*DV*DV/DLNP ) / DLNP ) / USUM
              TCG = ( DPTV - 
     &        ( DV*DPT + DT*DPV - 2.0*DP*DT*DV/DLNP ) / DLNP ) / USUM
            ELSE 
              PCG = SQRT ( PRE1 * PRE2 )
              TCG = TAVG
              ECG = 0.0
            END IF
          ELSE
            B = LOG(VMR1/VMR2) / DLNP
            IF ( ABS ( 1.0 + B ) .LE. 1.0E-3 
     &                .OR. DLNP .LT. DLPMIN ) THEN
              USUM = PRE1*VMR1*DLNP
              PCG  = DP2V / ( 2.0 + B ) / USUM
              ECG  = PRE1*VMR1
              TCG  = TAVG 
            ELSE 
              USUM = DPV / ( 1.0 + B )
              IF ( ABS ( 2.0 + B ) .LE. 1.0E-3 ) THEN
                PCG = PRE1*PRE1*VMR1*DLNP
              ELSE
                PCG = DP2V / ( 2.0 + B ) / USUM
              END IF
              ECG = 0.5 * DP2V2 / DPV
              TCG =  DPTV / DPV - DT / ( 1.0 + B ) / DLNP 
            END IF
          END IF
C So far USUM is (absorber amount)*(gM) and is in [mb]. 
C Required units are kmole.cm^-2. 
C So multiply by 100 (convert mb to Pa) for SI units, then divide by 10000
C to convert m^-2 to cm^-2, net factor 0.01
C
C Remove gM factor from USUM - see notes at top of subroutine.
C Factor 10 represents net factor 0.01 * 1000 m/km (convert HGT to SI)
          USUM = USUM * ABS ( HGT2 - HGT1 ) / RGAS / TAVG / DLNP * 10.0
C NB to make nadir path calculations independent of altitude scale replace
C above scaling of USUM with following line
c          USUM = USUM / ( 9.80 * 28.964 ) * 0.01
C
C Calculate required CG integrals for layer, convert p from mb to atm.
          PREPTH(IPTH) = PCG / ATMB
          PPAPTH(IPTH) = ECG / ATMB
          TEMPTH(IPTH) = TCG
C Set abs. amounts for different airmass (=sec(theta)) values stored as USRTAN
          IF ( FLXFLG ) THEN        ! FLX flag uses only single vertical path
            SECANG = 1.0
          ELSE IF ( USRELE ) THEN 
            SECANG = 1.0 / SIN ( ABS(USRTAN(1)) * DTORAD )
          ELSE
            SECANG = USRTAN(1)
          END IF
          AMTPTH(IPTH) = USUM * SECANG
          RAYPTH(IPTH) = ( HGT2 - HGT1 ) * SECANG
C For non-GRA cases, PSIPTH contain zenith angle of ray
          PSIPTH(IPTH) = ACOS ( 1.0/SECANG ) / DTORAD
          IF ( ZENFLG ) PSIPTH(IPTH) = 180.0 - PSIPTH(IPTH)
          IF ( .NOT. FLXFLG ) THEN       ! FLX only uses ITAN=1
            DO ITAN = 2, NTAN
              IF ( USRELE ) THEN 
                SECANG = 1.0 / SIN ( ABS(USRTAN(ITAN)) * DTORAD )
              ELSE
                SECANG = USRTAN(ITAN)
              END IF
              JPTH = IDXPTH ( ITAN, IATM, IGAS, 0 )          
              AMTPTH(JPTH) = USUM * SECANG
              RAYPTH(JPTH) = ( HGT2 - HGT1 ) * SECANG
              TEMPTH(JPTH) = TEMPTH(IPTH)
              PREPTH(JPTH) = PREPTH(IPTH)
              PPAPTH(JPTH) = PPAPTH(IPTH)
              PSIPTH(JPTH) = ACOS ( 1.0/SECANG ) / DTORAD
              IF ( ZENFLG ) PSIPTH(JPTH) = 180.0 - PSIPTH(JPTH)
            END DO
          END IF
        END DO
C
C For MTX+TRA flag, create extra paths for calculation of absorption coeff.
C at selected output levels (marked by IDIR=2).
        IF ( MTXFLG .AND. TRAFLG ) THEN
          DO ITAN = 1, NTAN 
            IATM = IATTAN(ITAN)
            IPTH = IDXPTH ( 1, IATM, IGAS, 2 )     ! 1=ITAN, 2=IDIR
            IF ( IPTH .EQ. 0 ) STOP 'F-RFMNAD: Logical Error'
            AMTPTH(IPTH) = VMRATM(IATM,IGAS) * DNSATM(IATM) ! set amt=densy/cm3
            PREPTH(IPTH) = PREATM(IATM) / ATMB
            PPAPTH(IPTH) = VMRATM(IATM,IGAS) * PREPTH(IPTH)  
            TEMPTH(IPTH) = TEMATM(IATM)
            PSIPTH(IPTH) = 0.0
            RAYPTH(IPTH) = 1.0E-5                 ! set = 1cm
          END DO
        END IF
      END DO
C
      END
