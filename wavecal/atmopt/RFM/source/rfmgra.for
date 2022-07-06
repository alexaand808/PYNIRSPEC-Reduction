      SUBROUTINE RFMGRA
C
C VERSION
C     23-APR-01  AD  Avoid problems with -ve VMR when calc. CG P, T
C     27-APR-00  AD  Use DSHALF as argument to RAYGRA instead of DS/2.0
C     03-JAN-00  AD  Original. Based on RFMRAY
C
C DESCRIPTION
C     Construct ray paths through 2D atmosphere.
C     Called once by RFM if GRAFLG enabled.
C     Calculate CG quantities for each path, defined as a combination of a 
C     single absorber and single layer.
C     Profile layers are sub-divided into segments (according to NP) purely for 
C     the ray-tracing within this module.
C
      IMPLICIT NONE
C
      EXTERNAL
     &  IDXPTH ! Indexing function for Pth# from Tan.pt#, Atm.Layer# and Gas#
     &, RAYGRA ! Ray-tracing in a 2-D atmosphere (z,psi).
     &, VALGRA ! Function interpolate value from 2D atmospheric field.
      INTEGER IDXPTH
      REAL    VALGRA
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data  
      INCLUDE 'crvcom.inc' ! Local radius of curvature
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'tancom.inc' ! Tangent heights
      INCLUDE 'pthcom.inc' ! RFM path data
C
C LOCAL CONSTANTS
      INTEGER       NP       ! Max. no.of path segments/layer for integration 
        PARAMETER ( NP = 10 )
C
C LOCAL VARIABLES
      INTEGER  IATM         ! Atmos. layer# for path
      INTEGER  IDIR         ! Direction counter
      INTEGER  IGAS         ! Gas# for path
      INTEGER  IPTH         ! Path counter
      INTEGER  ITAN         ! Tangent path# for path
      LOGICAL  FINISH       ! TRUE = last segment in current layer
      REAL     AMT1,AMT2,AMTM ! Absorber amounts at segment bounds & mid-pt
      REAL     DNS1,DNS2,DNSM ! Air molecular density [/cm^3]
      REAL     DS           ! Path length [km] of segment
      REAL     DSHALF       ! Half of DS
      REAL     DZMAX        ! Max segment thickness [km]
      REAL     P1, P2, PM   ! Horiz.angles [deg] of segment bounds & mid-pt
      REAL     PRE1,PRE2,PREM ! Pressure [mb]
      REAL     PNEAR        ! Horiz.angle [deg] at near boundary of layer
      REAL     T1, T2, TM   ! Zenith angles [deg] of segment bounds & mid-pt
      REAL     TEM1,TEM2,TEMM ! Temperatures [K]
      REAL     VMR1(MAXGAS) ! VMR at lower end of segment
      REAL     VMR2(MAXGAS) ! VMR at upper end of segment
      REAL     VMRM(MAXGAS) ! VMR at mid-pt of segment
      REAL     ZBOT,ZTOP    ! Low,High alts.[km] of layer
      REAL     Z1, Z2, ZM   ! Altitudes [km] of segment bounds & mid-pt
      REAL     ZTAN         ! Tangent height [km] or Geom.Tan height
C
      DOUBLE PRECISION DSUMDS         ! Cumulative path length
      DOUBLE PRECISION DSUMU(MAXGAS)  ! Sum of absorber amount U in layer
      DOUBLE PRECISION DSUMTU(MAXGAS) ! Sum of Temp*U in layer
      DOUBLE PRECISION DSUMPU(MAXGAS) ! Sum of Pres*U in layer
      DOUBLE PRECISION DSUMPP(MAXGAS) ! Sum of P.Pres*U in layer
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      PNEAR = 0.0                  ! Avoid warning about uninitialised variable
C
      DO ITAN = 1, MTAN                           ! Loop over all tangent paths
C ZTAN used for path length estimate DS. For surface intersecting rays, this
C needs to be based on geometric tangent height since HGTTAN=0.
        ZTAN = MIN ( HGTTAN(ITAN), GEOTAN(ITAN) )
        DO IDIR = -1, 1, 2                        ! Loop over each half
          Z2 = HGTTAN(ITAN)                       ! Initialise at tangent pt.
          P2 = PSITAN(ITAN)
C If IDIR=-1, ray from t.p. towards observer, so T is 0:90 (clockwise from zen)
C If IDIR=+1, ray from t.p. away from observer, so T is 0:-90 (a/c from zen)
          T2 = - SIGN ( SNGL(ASIN(SZNTAN(ITAN)))/DTORAD, FLOAT(IDIR))
          DNS2 = VALGRA ( Z2, P2, 'DNS', 0 )
          PRE2 = VALGRA ( Z2, P2, 'PRE', 0 ) 
          TEM2 = VALGRA ( Z2, P2, 'TEM', 0 ) 
          DO IGAS = 1, NGAS
            VMR2(IGAS) = VALGRA ( Z2, P2, 'VMR', IGAS )
          END DO
C
          DO IATM = IATTAN(ITAN), NATM - 1        ! Loop over atmos.layers
C
            IF ( IDIR .EQ. 1 ) PNEAR = P2  ! Near=lower boundary for up path
C Nominal path length DS based on approx NP intervals through each layer
C (in the case of tangent layer: NP intervals if t.p. were at base of layer)
            ZTOP  = HGTATM(IATM+1)
            ZBOT  = HGTATM(IATM)     
            DZMAX = ( ZTOP - ZBOT ) / NP
            IF ( ZTAN .GT. ZBOT ) THEN                 ! actual tangent layer
              DS = SQRT ( 2.0 * SNGL(RADCRV) * DZMAX ) 
            ELSE
              DS = SQRT ( 2.0 * SNGL(RADCRV) ) * 
     &              ( SQRT (ZBOT+DZMAX-ZTAN) - SQRT (ZBOT-ZTAN) )
            END IF
C
C Initialise CG integrations for this layer
            DSUMDS = 0.0D0
            DO IGAS = 1, NGAS
              DSUMU(IGAS)  = 0.0D0               ! Path gas amount
              DSUMTU(IGAS) = 0.0D0               ! Path average temperature
              DSUMPU(IGAS) = 0.0D0               ! Path average pressure
              DSUMPP(IGAS) = 0.0D0               ! Path average partial pressure
            END DO
C
            FINISH = .FALSE.
C
            DO WHILE ( .NOT. FINISH )        ! Repeat for each path segment
              Z1 = Z2
              P1 = P2
              T1 = T2
              DNS1 = DNS2
              PRE1 = PRE2
              TEM1 = TEM2
              DO IGAS = 1, NGAS
                VMR1(IGAS) = VMR2(IGAS)
              END DO
              CALL RAYGRA ( Z1, P1, T1, 0, DS, Z2, P2, T2 )
              IF ( Z2 .GT. ZTOP ) THEN      ! Stop at top of layer
                Z2 = ZTOP
                CALL RAYGRA ( Z1, P1, T1, 1, DS, Z2, P2, T2 ) ! DS changed
                FINISH = .TRUE.
              END IF
              DNS2 = VALGRA ( Z2, P2, 'DNS', 0 )
              PRE2 = VALGRA ( Z2, P2, 'PRE', 0 ) 
              TEM2 = VALGRA ( Z2, P2, 'TEM', 0 ) 
              DO IGAS = 1, NGAS
                VMR2(IGAS) = VALGRA ( Z2, P2, 'VMR', IGAS )
              END DO
              DSHALF = DS/2.0
              CALL RAYGRA ( Z1, P1, T1, 0, DSHALF, ZM, PM, TM )
              DNSM = VALGRA ( ZM, PM, 'DNS', 0 )
              PREM = VALGRA ( ZM, PM, 'PRE', 0 ) 
              TEMM = VALGRA ( ZM, PM, 'TEM', 0 ) 
              DO IGAS = 1, NGAS
                VMRM(IGAS) = VALGRA ( ZM, PM, 'VMR', IGAS )
              END DO
C
C Add to CG integrals using Simpson's rule (1,4,1)/6
              DSUMDS = DSUMDS + DBLE ( DS ) 
              DO IGAS = 1, NGAS
                AMT1 = DS * DNS1 * VMR1(IGAS) / 6.0
                AMT2 = DS * DNS2 * VMR2(IGAS) / 6.0
                AMTM = 4.0 * DS * DNSM * VMRM(IGAS) / 6.0
C For trapezoidal rule, uncomment the next 3 lines
C	amt1 = 0.5 * ds * dns1 * vmr1(igas)                
C	amt2 = 0.5 * ds * dns2 * vmr2(igas)                 
C	amtm = 0.0 
                DSUMU(IGAS)  = DSUMU(IGAS)  + DBLE(AMT1+AMTM+AMT2)
                DSUMTU(IGAS) = DSUMTU(IGAS) + 
     &            DBLE ( AMT1*TEM1 + AMTM*TEMM + AMT2*TEM2 )
                DSUMPU(IGAS) = DSUMPU(IGAS) + 
     &            DBLE ( AMT1*PRE1 + AMTM*PREM + AMT2*PRE2 )
                DSUMPP(IGAS) = DSUMPP(IGAS) + 
     &            DBLE ( AMT1*PRE1*VMR1(IGAS) + AMTM*PREM*VMRM(IGAS) 
     &                                        + AMT2*PRE2*VMR2(IGAS) )
              END DO
            END DO
C
C Form final integrated quantites
            DO IGAS = 1, NGAS
              IPTH = IDXPTH ( ITAN, IATM, IGAS, IDIR )
              IF ( DSUMU(IGAS) .NE. 0.0D0 .AND. 
     &             VMRATM(IATM,IGAS) .GE. 0.0 .AND.
     &             VMRATM(IATM+1,IGAS) .GE. 0.0 ) THEN
                TEMPTH(IPTH) = SNGL ( DSUMTU(IGAS) / DSUMU(IGAS) )
                PREPTH(IPTH) = SNGL ( DSUMPU(IGAS) / DSUMU(IGAS) )/ATMB
                PPAPTH(IPTH) = SNGL ( DSUMPP(IGAS) / DSUMU(IGAS) )/ATMB
C If gas amount is zero return averages for pressure and temperature
C Factor 1e5 in Amt calculation = cm/km
              ELSE
                TEMPTH(IPTH) = 0.5 * ( TEMATM(IATM) + TEMATM(IATM+1) )
                PREPTH(IPTH) = SQRT(PREATM(IATM)*PREATM(IATM+1))/ ATMB
                PPAPTH(IPTH) = 0.0
              ENDIF
              AMTPTH(IPTH) = 1.0E5 * SNGL ( DSUMU(IGAS) ) / AVOG 
              RAYPTH(IPTH) = SNGL ( DSUMDS )
              IF ( IDIR .EQ. -1 ) PNEAR = P2 ! Near=upper boundary for down path
              PSIPTH(IPTH) = PNEAR
            END DO
          END DO
        END DO
      END DO
C
      END
