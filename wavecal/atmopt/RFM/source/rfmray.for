      SUBROUTINE RFMRAY
C
C VERSION
C     24-OCT-13  AD  Bg#85 fix: Check for DELS12,23=0.0 for tan.ht = top of atm.
C     08-MAR-05  AD  Bug fix: Rename DYS to DY1 and initialise as 1-sin not 0
C     25-OCT-03  AD  Remove CHGATM test if OBS flag is true
C     20-JAN-03  AD  Add CHGATM test
C     28-NOV-01  AD  Check for DZ3.GE.ZTOP-ZMIN instead of .GE.ZTOP
C                    in order to avoid too thin a layer
C                    Rename DZMAX to DSMAX
C     23-APR-01  AD  Add check for -ve VMRs when calculating CG values
C     13-APR-01  AD  Check VMR value before interpolating log(VMR)
C     22-MAR-01  AD  Test CHGGAS to see if paths have changed
C     27-APR-00  AD  Correct calculation of DSUMPP
C     04-JAN-00  AD  Use more accurate weights for integration
C                    Convert RADCRV from SP to DP
C                    Add extra argument to IDXPTH
C     04-NOV-99  AD  Add MAX(*,0.0) before taking DCOS1=SQRT(*) for new layer
C     09-AUG-99  AD  Add explicit SNGL ( ) and DBLE ( )
C                    Convert DSH and DRFRn from S.P. to D.P.
C     11-MAY-99  AD  Limit DZMAX to at least 1 metre.
C     24-JAN-99  AD  Replace OFMFLG by LINFLG. Remove GL2FLG
C     29-JUL-98  AD  Initialise DSIN1 with SZNTAN for each ray path
C     14-JUL-98  AD  Comment change only.
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     17-SEP-96  AD  Use RADCRV rather than REARTH for loc.radius of curvature
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Construct ray paths through atmosphere.
C     Called once by RFM.
C     Calculate CG quantities for each path, defined as a combination of a 
C     single absorber and single layer.
C     Profile layers are sub-divided into segments (according to NP) purely 
C     for the ray-tracing within this module.
C
C     Each segment has a lower, mid.pt (by altitude) and upper boundary denoted
C     by 1,2,3 respectively. Further quarter-points, denoted by M, may be 
C     required near the tangent point for calculating COS term.
C
C     The algorithm is derived from GENLN2, which is itself taken from FASCODE.
C     The CG integration is new, and derived from integrating a quadratic
C     fitted to three points (y1,y2,y3) at (f,0,g).
C     The integral is:
C        I = (w1*y1 + w2*y2 + w3*y3)*h    (h is the total width, h=g-f).
C     where  w1 = (g+2*f)/(6*f), 
C            w3 = (f+2*g)/(6*g)
C            w2 = -h^2/(6*f*g)
C     NB: if f=-g, ie interval h is split equally, this reduces to 
C     (w1,w2,w3) = (1,4,1)/6, which is Simpson's rule. However, it is the 
C     altitude which is split in half, not the path length, so integrating  
C     along the path the split is not at the mid-point.
C
C REFERENCE
C     GALLERY, W.O, F.X. KNEIZYS and S.A. CLOUGH
C     Air Mass Computer Program for Atmospheric Transmittance/Radiance 
C     Calculation: FSCATM.
C     Report: AFGL-TR-83-0065, Hanscom AFB, Mass. 01731.
C
      IMPLICIT NONE
      SAVE            ! seems to be required to work under linux g77
C
      EXTERNAL
     &  IDXPTH ! Indexing function for Pth# from Tan.pt#, Atm.Layer# and Gas#
      INTEGER IDXPTH
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'atmcom.inc' ! Atmospheric profile data  
      INCLUDE 'crvcom.inc' ! Local radius of curvature
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'tancom.inc' ! Tangent heights
      INCLUDE 'pthcom.inc' ! RFM path data
C
C LOCAL CONSTANTS
      INTEGER       NP          ! Min. no.of path segments for integration [NP]
        PARAMETER ( NP = 10 )
C
C LOCAL VARIABLES
      INTEGER  IATM                  ! Atmos. layer# for path
      INTEGER  IGAS                  ! Gas# for path
      INTEGER  I,J                   ! Lower,Upper indices in interpolation
      INTEGER  IPTH                  ! Path counter
      INTEGER  ITAN                  ! Tangent path# for path
      REAL     DNS1,DNS2,DNS3        ! Air molecular density [/cm^3]
      REAL     FG1,FG2,FG3           ! Gas molecular density [/cm^3]
      REAL     FI,FJ                 ! Factors for I,J interp. in *ATM arrays
      REAL     FP1,FP2,FP3           ! Gas molec.density * Pres [mb/cm^3]
      REAL     FPP1,FPP2,FPP3        ! Gas molec.density * Part.Pres [mb/cm^3]
      REAL     FT1,FT2,FT3           ! Gas molec.density * Temp [K/cm^3]
      REAL     PRE1,PRE2,PRE3        ! Pressure [mb]
      REAL     TEM1,TEM2,TEM3        ! Temperatures [K]
      REAL     VMR1,VMR2,VMR3        ! Absorber vmrs [ppv]
      REAL     WGT1, WGT2, WGT3      ! Weights for CG integration
      REAL     ZBOT,ZTOP             ! Low,High alts.[km] of path in layer
      REAL     Z1, Z2, Z3, ZM        ! Altitudes [km] segment bounds & midpts
  
      DOUBLE PRECISION DCC                     ! n.r.sin(theta)
      DOUBLE PRECISION DCOS1,DCOS2,DCOS3       ! Cos of zenith angle 
      DOUBLE PRECISION DELS12, DELS23, DELS13  ! Path lengths 1-2,2-3,1-3 [km]
      DOUBLE PRECISION DELZ                    ! Alt.thickness of segment [km]
      DOUBLE PRECISION DRAD1,DRAD2,DRAD3,DRADM ! Rad.of DZn from centre earth
      DOUBLE PRECISION DRAT1,DRAT2,DRAT3,DRATM ! R/(n/n') (R=DRAD)
      DOUBLE PRECISION DRFR1,DRFR2,DRFR3,DRFRM ! Refractivity
      DOUBLE PRECISION DSDX1,DSDX2,DSDX3       ! dS/dX, where S=path length
      DOUBLE PRECISION DSH           ! Scale height for density or refrac. [km]
      DOUBLE PRECISION DSIN1,DSIN2,DSIN3,DSINM ! Sine of zenith angle 
      DOUBLE PRECISION DSMAX                   ! Max length of path segment[km]
      DOUBLE PRECISION DSRAY                   ! Length of ray
      DOUBLE PRECISION DSUMU,DSUMTU,DSUMPU,DSUMPP ! Sums for C-G integrals
      DOUBLE PRECISION DELX12,DELX23           ! DX increment from 1-2, and 2-3
      DOUBLE PRECISION DZ1,DZ2,DZ3,DZM         ! Altitudes [km] 
      DOUBLE PRECISION DZMIN                   ! Min segment thickness [km]
      DOUBLE PRECISION DX1,DX2,DX3             ! Intermediate variables
      DOUBLE PRECISION DY1,DY2,DY3             ! Integrals for eval.cos near tp
C
      LOGICAL FINAL     ! TRUE = last segment in current path 
      LOGICAL FIRST     ! TRUE = first segment in current path
C
C These DATA statements are purely to satisfy FLINT - actual values never used
        DATA FG3 / 0.0 /
        DATA FP3 / 0.0 /
        DATA FPP3 / 0.0 /
        DATA FT3 / 0.0 /
        DATA Z3  / 0.0 /
        DATA DCOS3 / 0.0 /
        DATA DRAD3 / 0.0 /
        DATA DRAT3 / 0.0 / 
        DATA DSDX3 / 0.0 /
        DATA DSIN3 / 0.0 /
        DATA DZ2 / 0.0 / 
        DATA DZ3 / 0.0 /
        DATA DX3 / 0.0 /
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Loop over all tangent paths
C
      DO ITAN = 1, MTAN
        DO IGAS = 1, NGAS
          IF ( .NOT. CHGGAS(IGAS) ) GOTO 200   ! Skip if paths unchanged
C
C Establish lower and upper altitudes of current path from atmosp.layer
C
          ZTOP = HGTTAN(ITAN)
          DSIN1 = SZNTAN(ITAN)       ! Initial zenith angle=1 at tangent point
          DO IATM = IATTAN(ITAN), NATM - 1
C If atmospheric layer unchanged only need to initialise DSIN1 and ZTOP to
C the values that would have been calculated at upper boundary
C NB cannot use this formula for observer within the atmosphere since GEOTAN
C is no longer the tangent point projected from conditions of n=1
            IF ( .NOT. CHGATM(IATM) .AND. .NOT. OBSFLG ) THEN
              ZTOP = HGTATM(IATM+1)
              DSIN1 = ( DBLE ( GEOTAN(ITAN) ) + RADCRV ) /
     &                ( DBLE ( RFRATM(IATM+1) ) + 1.0D0 ) /
     &                ( DBLE ( ZTOP ) + RADCRV )
              GOTO 100
            END IF
            IPTH = IDXPTH ( ITAN, IATM, IGAS, 0 )
            I = IATM
            J = IATM + 1 
            ZBOT = ZTOP                 ! Lower layer hgt or Tan.Pt.
            Z1 = ZBOT                   ! Initial lower boundary altitude
            DZ1 = DBLE ( Z1 )           ! DP version
            ZTOP = HGTATM(IATM+1)
            DSH = DBLE ( DSHATM(IATM) ) ! Scale Ht for refractivity [km]
C
            FINAL = .FALSE.
            FIRST = .TRUE.
C
C Distance from centre of earth
            DRAD1 = DZ1 + RADCRV
C
C Nominal path segment length based on layer thickness/NP 
C (Limit to .GE. 1 metre in case t.p. ZBOT close to ZTOP )
            DSMAX = DBLE ( MAX ( ( ZTOP - ZBOT ) / NP, 0.001 ) ) 
C
C Set a minimum vertical thickness of a segment based on a horizontal
C distance DSMAX through a layer (from ds^2 = 2.R.dz approximately)
            DZMIN = 0.5D0 * DSMAX**2 / DRAD1
C
C Interpolate required quantities at Z1
C
            FJ = ( Z1 - HGTATM(I) ) / ( HGTATM(J) - HGTATM(I) )
            FI = 1.0 - FJ
            DRFR1 = DBLE ( EXP ( FI * LNRATM(I) + FJ * LNRATM(J) ) )
C
C Calculate snells law constant R = - r/(n/n'), =0 if refractivity=0
C
            IF ( GEOFLG ) THEN
              DRAT1 = 0.0D0
              DCC = DSIN1 * DRAD1
            ELSE
              DRAT1 = DRAD1 * DRFR1 / ( DSH * (1.D0 + DRFR1) ) 
              DCC = DSIN1 * DRAD1 * ( DRFR1 + 1.D0 )
            END IF
C
            DNS1 = EXP ( FI * LNDATM(I) + FJ * LNDATM(J) )
            IF ( LINFLG .OR. VMRATM(I,IGAS) .LT. ARGMIN 
     &                  .OR. VMRATM(J,IGAS) .LT. ARGMIN ) THEN
              VMR1 = FI * VMRATM(I,IGAS) + FJ * VMRATM(J,IGAS) 
            ELSE
              VMR1 = EXP ( FI * LNVATM(I,IGAS) + FJ * LNVATM(J,IGAS) )
            END IF
            PRE1 = EXP ( FI * LNPATM(I) + FJ * LNPATM(J) )
            TEM1 = FI * TEMATM(I) + FJ * TEMATM(J) 
C
C Initialise layer lower boundary weighted quantities for integration
C
            DSUMU = 0.D0                        ! Path gas amount
            DSUMTU = 0.D0                       ! Path average temperature
            DSUMPU = 0.D0                       ! Path average pressure
            DSUMPP = 0.D0                       ! Path average partial pressure
            FG1  = DNS1 * VMR1                     
            FT1  = FG1 * TEM1
            FP1  = FG1 * PRE1
            FPP1 = FP1 * VMR1 
C
C Initialize lower boundary values (NB DSIN1 may not have been accurately
C calculated if previous layer always used the near t.p. approximation so set
C MAX(*,0.0) to be safe.
C Cos of initial zenith angle
            DCOS1 = DSQRT ( MAX ( 1.D0 - DSIN1**2, 0.D0 ) ) 
            DY1 = MAX ( 1.0D0 - DSIN1, 0.D0 )        ! Matches DCOS1 definition
            DX1 = DRAD1 * DCOS1   ! Change to intermediate variable X
            DSDX1 = 1.D0/(1.D0 - DRAT1 * DSIN1**2 ) ! ds/dx
            DSRAY = 0.D0          ! Ray length is zero at this stage
C
C Loop over sub-layers
C Calculate extent of this sub-layer 
C
C Allow for cos-->0 at the tangent height
C
            DO WHILE ( .NOT. FINAL )
C
C Update quantites at new layer level 1 from old layer level 3
C
              IF ( .NOT. FIRST ) THEN
                DZ1 = DZ3
                Z1 = Z3
                DRAD1 = DRAD3
                DSIN1 = DSIN3
                DCOS1 = DCOS3
                DX1 = DX3
                DRAT1 = DRAT3
                DSDX1 = DSDX3
                FG1 = FG3
                FT1 = FT3
                FP1 = FP3
                FPP1 = FPP3
              END IF
              FIRST = .FALSE.
C
C Calculate altitude increment for next segment
C
              DELZ = MAX ( DSMAX * DCOS1, DZMIN )
              DZ3 = DZ1 + DELZ
              IF ( DZ3 .GE. DBLE(ZTOP) - DZMIN ) THEN ! top of layer reached
                DZ3 = DBLE ( ZTOP )
                DELZ = DZ3 - DZ1
                FINAL = .TRUE.
              END IF
              DZ2 = 0.5D0 * ( DZ1 + DZ3 )
              Z2 = SNGL ( DZ2 )
              Z3 = SNGL ( DZ3 )
C
C Calculate refracted path parameters at Z2
C
              DRAD2 = RADCRV + DZ2
              FJ = ( Z2 - HGTATM(I) ) / ( HGTATM(J) - HGTATM(I) )
              FI = 1.0 - FJ
              DRFR2 = DBLE ( EXP ( FI * LNRATM(I) + FJ * LNRATM(J) ) )
              IF ( GEOFLG ) THEN
                DRAT2 = 0.0D0
                DSIN2 = DCC / DRAD2
              ELSE
                DRAT2 = DRAD2 * DRFR2 / ( DSH * (1.D0 + DRFR2) ) 
                DSIN2 = DCC / ( (DRFR2 + 1.D0) * DRAD2 )
              END IF
C
C Interpolate required quantities at Z2 and evaluate required functions
C
              DNS2 = EXP ( FI * LNDATM(I) + FJ * LNDATM(J) )
              IF ( LINFLG .OR. VMRATM(I,IGAS) .LT. ARGMIN 
     &                    .OR. VMRATM(J,IGAS) .LT. ARGMIN ) THEN 
                VMR2 = FI * VMRATM(I,IGAS) + FJ * VMRATM(J,IGAS) 
              ELSE
                VMR2 = EXP ( FI * LNVATM(I,IGAS) + FJ * LNVATM(J,IGAS) )
              END IF
              PRE2 = EXP ( FI * LNPATM(I) + FJ * LNPATM(J) )
              TEM2 = FI * TEMATM(I) + FJ * TEMATM(J) 
              FG2 = DNS2 * VMR2 
              FT2 = FG2 * TEM2
              FP2 = FG2 * PRE2
              FPP2 = FP2 * VMR2
C
C Calculate refracted path parameters at Z3
C
              DRAD3 = RADCRV + DZ3
              FJ = ( Z3 - HGTATM(I) ) / ( HGTATM(J) - HGTATM(I) )
              FI = 1.0 - FJ
              DRFR3 = DBLE ( EXP ( FI * LNRATM(I) + FJ * LNRATM(J) ) )
              IF ( GEOFLG ) THEN                 ! No refraction
                DRAT3 = 0.0D0
                DSIN3 = DCC / DRAD3
              ELSE
                DRAT3 = DRAD3 * DRFR3 / ( DSH * (1.D0 + DRFR3) ) 
                DSIN3 = DCC / ( (DRFR3 + 1.D0) * DRAD3 )
              END IF
C
C Interpolate required quantities at Z3 and evaluate required functions
C
              DNS3 = EXP ( FI * LNDATM(I) + FJ * LNDATM(J) )
              IF ( LINFLG .OR. VMRATM(I,IGAS) .LT. ARGMIN
     &                    .OR. VMRATM(J,IGAS) .LT. ARGMIN ) THEN
                VMR3 = FI * VMRATM(I,IGAS) + FJ * VMRATM(J,IGAS) 
              ELSE
                VMR3 = EXP ( FI * LNVATM(I,IGAS) + FJ * LNVATM(J,IGAS) )
              END IF
              PRE3 = EXP ( FI * LNPATM(I) + FJ * LNPATM(J) )
              TEM3 = FI * TEMATM(I) + FJ * TEMATM(J) 
              FG3 = DNS3 * VMR3 
              FT3 = FG3 * TEM3 
              FP3 = FG3 * PRE3 
              FPP3 = FP3 * VMR3
C
C Need the cosine of the angle, leads to precision problems. near the tangent 
C height
C
              IF ( (1.D0 - DSIN2) .GT. 1.0D-5 ) THEN    ! Not too close to t.p.
                DCOS2 = DSQRT ( 1.D0 - DSIN2**2 )
                DCOS3 = DSQRT ( 1.D0 - DSIN3**2 )
              ELSE                                      ! Close to t.p.
                DZM = 0.5D0 * ( DZ1 + DZ2 )
                ZM = SNGL ( DZM )
                DRADM = DZM + RADCRV
                IF ( GEOFLG ) THEN
                  DRATM = 0.D0
                  DSINM = DCC / DRADM
                ELSE
                  FJ = ( ZM - HGTATM(I) ) / ( HGTATM(J) - HGTATM(I) )
                  FI = 1.0 - FJ
                  DRFRM = DBLE ( EXP ( FI * LNRATM(I) + FJ * LNRATM(J)))
                  DRATM = DRADM * DRFRM / ( DSH * (1.D0 + DRFRM) ) 
                  DSINM = DCC / ( (DRFRM + 1.D0) * DRADM )
                END IF
                DY2 = DY1 + ( DZM - DZ1 ) / 3.D0 * ( 
     &                  ( 1.D0 - DRAT1) * DSIN1 / DRAD1 +
     &                  4.D0 * (1.D0 - DRATM) * DSINM / DRADM + 
     &                  (1.D0 - DRAT2) * DSIN2 / DRAD2          )
                DCOS2 = DSQRT ( 2.D0*DY2 - DY2**2 )
                DZM = 0.5D0 * ( DZ2 + DZ3 )
                ZM = SNGL ( DZM ) 
                DRADM = DZM + RADCRV
                IF ( GEOFLG ) THEN
                  DRFRM = 0.0D0
                  DRATM = 0.0D0
                ELSE
                  FJ = ( ZM - HGTATM(I) ) / ( HGTATM(J) - HGTATM(I) )
                  FI = 1.0 - FJ
                  DRFRM = DBLE ( EXP ( FI * LNRATM(I) + FJ * LNRATM(J)))
                  DRATM = DRADM * DRFRM / ( DSH * (1.D0 + DRFRM) ) 
                ENDIF
                DSINM = DCC / ( (DRFRM + 1.D0) * DRADM )
                DY3 = DY2 + ( DZM - DZ2 ) / 3.D0 * ( 
     &                  ( 1.D0 - DRAT2 ) * DSIN2 / DRAD2 +
     &                  4.D0 * (1.D0 - DRATM) * DSINM / DRADM + 
     &                  ( 1.D0 - DRAT3 ) * DSIN3 / DRAD3        )
                DCOS3 = DSQRT ( 2.D0 * DY3 - DY3**2 )
                DY1 = DY3
              ENDIF            ! End special case for calculating Cos
C
              DX2 = DRAD2 * DCOS2
              DX3 = DRAD3 * DCOS3
              DELX12 = DX2 - DX1
              DELX23 = DX3 - DX2
C
C Ray path length
C
              DSDX2 = 1.D0 / ( 1.D0 - DRAT2 * DSIN2**2 )
              DSDX3 = 1.D0 / ( 1.D0 - DRAT3 * DSIN3**2 )
              DELS12 = 0.5D0 * DELX12 * ( DSDX1 + DSDX2 )
              DELS23 = 0.5D0 * DELX23 * ( DSDX2 + DSDX3 )
              DELS13 = DELS12 + DELS23
              DSRAY = DSRAY + DELS13
C
C Path gas amount, Av.Temp, Av. Pres, and Av.P.Pres 
C The following integral is for a quadratic fitted to three points
C (f=-DELS12, g=DELS23 h=DELS13 from notes at top).
C
C If tangent height is top of atmosphere, DELS12,23=0.0 so avoid division by 0
              IF ( DELS12 .EQ. 0.0 .OR. DELS23 .EQ. 0.0 ) THEN
                WGT1 = 0.0
                WGT2 = 0.0
                WGT3 = 0.0
              ELSE
                WGT1 = SNGL(DELS13*(2.0D0*DELS12-DELS23)/(6.0D0*DELS12))
                WGT3 = SNGL(DELS13*(2.0D0*DELS23-DELS12)/(6.0D0*DELS23))
                WGT2 = SNGL(DELS13**3/(6.0D0*DELS12*DELS23))
              END IF
C
              DSUMU = DSUMU  + DBLE ( WGT1*FG1 + WGT2*FG2 + WGT3*FG3 ) 
              DSUMTU= DSUMTU + DBLE ( WGT1*FT1 + WGT2*FT2 + WGT3*FT3 )
              DSUMPU= DSUMPU + DBLE ( WGT1*FP1 + WGT2*FP2 + WGT3*FP3 )
              DSUMPP= DSUMPP + DBLE ( WGT1*FPP1+ WGT2*FPP2+ WGT3*FPP3)
            END DO
C
C Form final integrated quantites
C
            IF ( DSUMU .NE. 0.0D0 .AND. VMRATM(I,IGAS) .GE. 0.0
     &           .AND. VMRATM(J,IGAS) .GE. 0.0 ) THEN
              TEMPTH(IPTH) = SNGL ( DSUMTU / DSUMU )
              PREPTH(IPTH) = SNGL ( DSUMPU / ( DSUMU * DBLE(ATMB) ) )
              PPAPTH(IPTH) = SNGL ( DSUMPP / ( DSUMU * DBLE(ATMB) ) )
C
C If gas amount is zero return averages for pressure and temperature
C Factor 1e5 in Amt calculation = cm/km
            ELSE
              TEMPTH(IPTH) = 0.5 * ( TEMATM(I) + TEMATM(J) )
              PREPTH(IPTH) = SQRT ( PREATM(I) * PREATM(J) ) / ATMB
              PPAPTH(IPTH) = 0.0
            ENDIF
            AMTPTH(IPTH) = 1.0E5 * SNGL(DSUMU) / AVOG 
            RAYPTH(IPTH) = SNGL ( DSRAY )
C
C Save path parameter DSIN1 required for next layer 
C ZTOP also required, everything else is recalculated using the FIRST flag.
            DSIN1 = DSIN3
            PSIPTH(IPTH) = SNGL ( ASIN ( DSIN1 ) ) / DTORAD
C
  100       CONTINUE
          END DO
  200     CONTINUE
        END DO
      END DO
C
      END
