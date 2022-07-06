      SUBROUTINE RFMFLX 
C
C VERSION
C     23-APR-12  AD  Allow for spectrally varying surface emissivity
C     26-SEP-08  AD  Bug#72: SAVE ITNATM,WGTATM between calls
C     03-DEC-07  AD  Bug#67: SAVE NQAD,XQAD,WQAD between calls
C     03-MAY-05  AD  Add jdxcon.inc
C     15-JAN-03  AD  Add FLXEFN
C     10-MAR-03  AD  Set PIFAC correctly
C     17-DEC-01  AD  Initialise TEMPTB and PIFAC
C     04-NOV-01  AD  Limit calculations to IATTAN(1) & upwards for ZEN flag.
C                    Allow OPTFLG to store optical depths
C     31-OCT-01  AD  Add VRTFLG
C                    Correction: avoid accumulating RADTAN if downward path
C                    calculation above a reflecting surface
C     28-MAY-00  AD  Convert RADTAN and ABSTAN from S.P. to D.P.
C                    And use these directly for fluxes and cooling rates
C                    Also convert PLANCK to return D.P. BBFSFC
C     27-APR-00  AD  Original.
C
C DESCRIPTION
C     Perform RFM flux calculation
C     Called by RFM for each widemesh interval if FLX flag is selected.
C
C     Cooling rates are calculated from dT/dt = dF/dz*(1/rho.cp)
C     where rho = atmospheric density
C     dF/dz is calculated from a quadratic fit to F(z) at 3 successive levels
C     in the atmospheric profile.
C     F is in nW/cm2..., z in km, rho in molec/cm3, cp in J/K/kmol=W.s/K/kmol, 
C     so dT/dt*(1e-9W/nW)*(1e-5km/cm)*(Na molec/kmol)*(864000s/day)= K/day
C
      IMPLICIT NONE
C
      EXTERNAL
     &  FLXEFN ! Absorption Coeff. calc. for Escape Fn normalisation
     &, FLXTRA ! Radiative transfer calculation for flux calculations.
     &, GAUQAD ! Values and Weights for Gaussian First Moment Quadrature
     &, PLANCK ! Calculate Planck emission spectrum at specified temperature
     &, PTBATM ! Perturb/unperturb atmospheric profiles for Jacobian calc.
C
C GLOBAL CONSTANTS
      INCLUDE 'jdxcon.inc' ! Indices for non-VMR Jacobians
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
      INCLUDE 'rfmcon.inc' ! Constants for RFM calculations
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data 
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'jaccom.inc' ! Jacobian data
      INCLUDE 'sfccom.inc' ! Surface parameters
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL CONSTANTS
      INTEGER       MAXWGT       ! Max.no o/p level cooling rates for each IATM
        PARAMETER ( MAXWGT = 5 ) ! Quadratic fit implies a maximum=5
      INTEGER       MAXQAD       ! No.points required for quadrature
        PARAMETER ( MAXQAD = 4 ) ! Must be .LE. number listed in GAUQAD.FOR
C
C LOCAL VARIABLES
      INTEGER  IATM       ! Atmospheric level/layer counter
      INTEGER  IATM1      ! Lowest layer# required
      INTEGER  IFIN       ! Fine mesh grid point counter
      INTEGER  IJAC       ! Jacobian counter
      INTEGER  IQAD       ! Counter for quadrature points
      INTEGER  ITAN       ! Output level index (1:NTAN)
      INTEGER  IOPATM(MAXATM) ! Index of o/p level (if any) for each atm level
      INTEGER  IPTB       ! Perturbed layer/secondary level for MTX calc.
      INTEGER  ITNATM(MAXWGT,MAXATM) ! o/p level for each weight
      INTEGER  IWGT       ! Counter over weights for each ATM level
      INTEGER  JATM       ! Atmospheric level/layer counter
      INTEGER  JTAN       ! Index of unptb/jac ray in calc.paths (1:NJAC+1)
      INTEGER  KATM       ! Atmospheric level/layer counter
      INTEGER  NQAD       ! No. quadrature points used (either=MAXQAD or 1)
      INTEGER  NWTATM(MAXATM) ! No.weights for each ATM level (1:MAXWGT)
      LOGICAL  FIRST      ! T=first call, F=subsequent call
      LOGICAL  TEMPTB     ! T=current temperature profile is perturbed.
      LOGICAL  SETOPT     ! T=calculate optical depths
      LOGICAL  SETTRA     ! T=calculate flux transmittances
      REAL     COOFAC     ! Factor for cooling rates
      REAL     SECFAC     ! Scale factor for absorber amount
      REAL     SRCFAC     ! Perturbation applied to source functions
      REAL     WGTATM(MAXWGT,MAXATM) ! Weights for cooling rates
      REAL     XIJ,XJK,XKI ! Altitude separation between atm.levels
      DOUBLE PRECISION ABSOFF(MAXFIN) ! Offset optical depth  MTX+TRA/ABS calcs
      DOUBLE PRECISION ABSQAD(MAXFIN) ! Integrated optical depth along path
      DOUBLE PRECISION BBFSFC(MAXFIN) ! Planck function for surface temperature
      DOUBLE PRECISION EMSFIN(MAXFIN) ! Surface emissivity
      DOUBLE PRECISION PIFAC          ! Either pi or 1 (VRTFLG)
      DOUBLE PRECISION RADQAD(MAXFIN) ! Integrated radiance along path
      DOUBLE PRECISION RADSFC(MAXFIN) ! Radiance from surface (incl.reflected)
      DOUBLE PRECISION XQAD(MAXQAD)   ! Quadrature points = cos(zenith angle)
      DOUBLE PRECISION WQAD(MAXQAD)   ! Weights for quadrature points
C
C DATA STATEMENTS
      DATA FIRST  / .TRUE. /
      DATA IOPATM / MAXATM * 0 /
      DATA NWTATM / MAXATM * 0 /
      DATA PIFAC / 0.0D0 /  ! Avoids warning message over unitialized variable
      DATA XQAD / MAXQAD * 0.0D0 /
      DATA WQAD / MAXQAD * 0.0D0 /
      SAVE FIRST, IOPATM, NQAD, NWTATM, PIFAC, XQAD, WQAD
      SAVE ITNATM, WGTATM
C       
C EXECUTABLE CODE -------------------------------------------------------------
C
C On first call of this subroutine...
      IF ( FIRST ) THEN
        IF ( VRTFLG ) THEN
          NQAD = 1
          XQAD(1) = 1.0D0
          PIFAC   = 1.0D0
          WQAD(1) = 1.0D0
        ELSE
          NQAD = MAXQAD
          PIFAC = DBLE ( PI )
          CALL GAUQAD ( NQAD, XQAD, WQAD )    ! Set quadrature points
          DO IQAD = 1, NQAD
            WQAD(IQAD) = WQAD(IQAD) * 2.0D0 * PIFAC 
          END DO
        END IF
        DO ITAN = 1,  NTAN                   ! Find atm layers for each o/p lev
          JATM = IATTAN(ITAN)
          IOPATM(JATM) = ITAN                  !   atm layer at output level
C For cooling rates, fit quadratic to atm.levels i,j,k and calculate gradient
C at level j. Coefficients do not depend on sequence i,j,k so can also be used
C if j is one of the end points (if JATM=1 or NATM).
          IF ( COOFLG ) THEN
            IATM = JATM-1                        ! sequence i,j,k
            KATM = JATM+1
            IF ( IATM .LT. 1 ) IATM = 3          ! sequence j,k,i 
            IF ( KATM .GT. NATM ) KATM = NATM-2  ! sequence k,i,j
            NWTATM(IATM) = NWTATM(IATM) + 1
            NWTATM(JATM) = NWTATM(JATM) + 1
            NWTATM(KATM) = NWTATM(KATM) + 1
            ITNATM(NWTATM(IATM),IATM) = ITAN
            ITNATM(NWTATM(JATM),JATM) = ITAN
            ITNATM(NWTATM(KATM),KATM) = ITAN
            XIJ = HGTATM(IATM) - HGTATM(JATM)
            XJK = HGTATM(JATM) - HGTATM(KATM)
            XKI = HGTATM(KATM) - HGTATM(IATM)
            COOFAC = (AVOG * 86400.0E-14) / DNSATM(JATM) / CPMOL
            WGTATM(NWTATM(IATM),IATM) = - COOFAC * XJK / XIJ / XKI
            WGTATM(NWTATM(JATM),JATM) = COOFAC * ( 1.0/XJK - 1.0/XIJ ) 
            WGTATM(NWTATM(KATM),KATM) = COOFAC * XIJ / XJK / XKI
          END IF
        END DO  
        FIRST = .FALSE.
      END IF
C
      TEMPTB = .FALSE.
      IPTB = 0
      IJAC = 0
      CALL SFCEMS ( NFIN, WNOFIN, EMSFIN ) ! Calculate sfc.emiss. on fine grid
 100  CONTINUE    ! Repeat from here with different JTAN,IJAC if JACobians reqd
C
C Initialise output arrays to zero (contributions from each layer summed here)
      DO ITAN = 1, NTAN              ! 1:NTAN first time, NTAN+1:2*NTAN IJAC=1
        JTAN = ITAN                                   ! unperturbed paths
        IF ( IJAC .GT. 0 ) JTAN = ITNJAC(ITAN,IJAC)   ! perturbed paths
        DO IFIN = 1, NFIN
          RADTAN(IFIN,JTAN) = 0.0D0 
          ABSTAN(IFIN,JTAN) = 0.0D0  ! Used to store cooling rates or opt.depth
          TRATAN(IFIN,JTAN) = 0.0D0  ! NB initialise to 0 for summation, not 1
        END DO
      END DO     
      ITAN = IJAC + 1  ! Index for path calculations (not output levels)
C
C If ZEN option, don't need to integrate down to the surface
      IF ( ZENFLG ) THEN             ! Matches logic in FLXPTH
        IATM1 = IATTAN(1) 
        IF ( COOFLG .AND. IATM1 .GT. 1 ) IATM1 = IATM1 - 1
      ELSE
C Set up surface emission contribution to upward flux
C (surface parameters are undefined if ZENFLG is TRUE)
        IATM1 = 1
        CALL PLANCK ( TEMSFC, NFIN, WNOFIN, BBFSFC ) 
        DO IFIN = 1, NFIN           ! Reduce by emissivity factor
          RADSFC(IFIN) = BBFSFC(IFIN) * EMSFIN(IFIN)
        END DO
      END IF
      SETTRA = ( TRAFLG .OR. ABSFLG )
      SETOPT = OPTFLG
C If there is no downward contribution, skip the first part
      IF ( NADFLG .AND. .NOT. RFLSFC ) GOTO 200  
C
C Downward radiance flux calculation
      DO IQAD = 1, NQAD       ! Loop over slant path for each quadrature point 
C
        SECFAC = 1.0 / SNGL ( XQAD(IQAD) )      ! Sec(theta) scale factor
        DO IFIN = 1, NFIN                ! Initialise calc. ray path
          RADQAD(IFIN) = 0.0D0           ! Stores integrated radiance
          ABSQAD(IFIN) = 0.0D0           ! Stores integrated optical depth
          ABSOFF(IFIN) = 0.0D0            ! Initialise for IPTB=NATM
        END DO
C
C MTX: only calculate transmittances from level IPTB downwards
        IF ( MTXFLG ) SETTRA = .FALSE.
        IF ( MTXFLG ) SETOPT = .FALSE.
        DO IATM = NATM, IATM1, -1
          IF ( IATM .EQ. IPTB ) THEN
            SRCFAC = ( 1.0 + PTBSRC )
          ELSE
            SRCFAC = 1.0
          END IF
C If constructing flux transmittance matrix, require ABSOFF to be set for  
C secondary level before calculating transmittances to lower levels.
          IF ( IATM .LT. NATM ) 
     &       CALL FLXTRA ( RADQAD, ABSQAD, ITAN, IATM, -1, 
     &                                      SECFAC, SRCFAC ) 
          IF ( IATM .EQ. IPTB ) THEN   ! Only if MTXFLG since IPTB=0 otherwise
            DO IFIN = 1, NFIN
              ABSOFF(IFIN) = ABSQAD(IFIN)
            END DO
            SETTRA = ( TRAFLG .OR. ABSFLG ) 
            SETOPT = OPTFLG
          END IF
          IF ( IOPATM(IATM) .NE. 0 .AND. .NOT. NADFLG ) THEN
            JTAN = IOPATM(IATM)          ! ATM level coincides with output ITAN
            IF ( IJAC .GT. 0 ) JTAN = ITNJAC(JTAN,IJAC)
            DO IFIN = 1, NFIN 
              RADTAN(IFIN,JTAN) = RADTAN(IFIN,JTAN) 
     &                          - RADQAD(IFIN) * WQAD(IQAD) 
            END DO
            IF ( SETTRA ) THEN
              DO IFIN = 1, NFIN
                TRATAN(IFIN,JTAN) = TRATAN(IFIN,JTAN) + 
     &            EXP ( - ( ABSQAD(IFIN) - ABSOFF(IFIN) ) ) 
     &            * WQAD(IQAD) / PIFAC
              END DO
            END IF          
            IF ( SETOPT ) THEN
              DO IFIN = 1, NFIN
                ABSTAN(IFIN,JTAN) = ABSTAN(IFIN,JTAN) + 
     &            ( ABSQAD(IFIN) - ABSOFF(IFIN) ) 
     &            * WQAD(IQAD) / PIFAC
              END DO
            END IF          
          END IF          
          IF ( COOFLG ) THEN
            DO IWGT = 1, NWTATM(IATM)
              JTAN = ITNATM(IWGT,IATM) 
              IF ( IJAC .GT. 0 ) JTAN = ITNJAC(JTAN,IJAC)
              DO IFIN = 1, NFIN 
                ABSTAN(IFIN,JTAN) = ABSTAN(IFIN,JTAN) 
     &             - RADQAD(IFIN) * WQAD(IQAD) * DBLE(WGTATM(IWGT,IATM))
              END DO
            END DO
          END IF          
        END DO
C
C Contribution from Surface term. 
C Reflected *radiance* is downwelling *irradiance* divided by pi
        IF ( RFLSFC ) THEN            ! Allow for reflected component
          DO IFIN = 1, NFIN
            RADSFC(IFIN) = RADSFC(IFIN) + 
     &        WQAD(IQAD) * RADQAD(IFIN) * (1.0D0-EMSFIN(IFIN)) / PIFAC
          END DO
        END IF
      END DO
C
      IF ( ZENFLG ) GOTO 300                
 200  CONTINUE
C
C Upward radiance flux calculation
      DO IQAD = 1, NQAD
C
        SECFAC = 1.0 / SNGL ( XQAD(IQAD) )    ! Sec(theta) scale factor
        DO IFIN = 1, NFIN                ! Initialise calc. ray path
          RADQAD(IFIN) = RADSFC(IFIN)    ! Stores integrated radiance
          ABSQAD(IFIN) = 0.0D0           ! Stores integrated optical depth
          ABSOFF(IFIN) = 0.0D0           ! Initialise for IPTB=1
        END DO
C
C NB: IATM layer indices differ from level indices by one on the upward 
C path, so call FLXTRA at end of IATM loop rather than beginning.
C
        IF ( MTXFLG ) SETTRA = .FALSE.
        IF ( MTXFLG ) SETOPT = .FALSE.
        DO IATM = 1, NATM
          IF ( IATM .EQ. IPTB ) THEN
            SRCFAC = ( 1.0 + PTBSRC )
            DO IFIN = 1, NFIN
              ABSOFF(IFIN) = ABSQAD(IFIN)
            END DO
          ELSE
            SRCFAC = 1.0
          END IF
          IF ( IOPATM(IATM) .NE. 0 ) THEN
            JTAN = IOPATM(IATM)          ! ATM level coincides with output ITAN
            IF ( IJAC .GT. 0 ) JTAN = ITNJAC(JTAN,IJAC)
            DO IFIN = 1, NFIN 
              RADTAN(IFIN,JTAN) = RADTAN(IFIN,JTAN) 
     &                          + RADQAD(IFIN) * WQAD(IQAD) 
            END DO
            IF ( SETTRA ) THEN 
              DO IFIN = 1, NFIN
                TRATAN(IFIN,JTAN) = TRATAN(IFIN,JTAN) +
     &            EXP ( - ( ABSQAD(IFIN) - ABSOFF(IFIN) ) ) 
     &            * WQAD(IQAD) / PIFAC
              END DO
            END IF
            IF ( SETOPT ) THEN 
              DO IFIN = 1, NFIN
                ABSTAN(IFIN,JTAN) = ABSTAN(IFIN,JTAN) +
     &            ( ABSQAD(IFIN) - ABSOFF(IFIN) ) 
     &            * WQAD(IQAD) / PIFAC
              END DO
            END IF
          END IF          
          IF ( COOFLG ) THEN
            DO IWGT = 1, NWTATM(IATM)
              JTAN = ITNATM(IWGT,IATM) 
              IF ( IJAC .GT. 0 ) JTAN = ITNJAC(JTAN,IJAC)
              DO IFIN = 1, NFIN 
                ABSTAN(IFIN,JTAN) = ABSTAN(IFIN,JTAN) 
     &             + RADQAD(IFIN) * WQAD(IQAD) * DBLE(WGTATM(IWGT,IATM))
              END DO
            END DO
          END IF          
          IF ( IATM .LT. NATM ) 
     &      CALL FLXTRA ( RADQAD, ABSQAD, ITAN, IATM, 1, SECFAC, SRCFAC)
          IF ( IATM .EQ. IPTB ) SETTRA = ( TRAFLG .OR. ABSFLG )
          IF ( IATM .EQ. IPTB ) SETOPT = OPTFLG
        END DO
C
      END DO
C
 300  CONTINUE
C
C Last part sets up perturbed profiles for repeated calcs for Jacobians
C (each call undoes the perturbations set up by the previous call)
      IF ( JACFLG ) THEN
        IF ( IJAC .LT. NJAC ) THEN
          IJAC = IJAC + 1
          IF ( IGSJAC(IJAC) .EQ. MAXGAS+JDXTEM .AND. BFXFLG ) THEN  ! =Tem.Ptb.
            CALL PTBATM ( IJAC )
            TEMPTB = .TRUE.
          ELSE IF ( TEMPTB ) THEN
            CALL PTBATM ( NJAC+1 )       ! Undo last temperature perturbation
          END IF
          GOTO 100                       ! Repeat calculation for IJAC
C If last ptb was a temperature perturbation, this needs to be undone since
C actual temperatures required for Planck functions in radiative transfer
C with BFX flag = TRUE. Other (VMR, PRE) perturbations can be left.
        END IF
        IF ( TEMPTB ) CALL PTBATM ( NJAC+1 ) ! Undo temperature ptb.
      ELSE IF ( MTXFLG ) THEN
        IF ( IJAC .LT. NJAC ) THEN
          IJAC = IJAC + 1
          IPTB = IATJAC(IJAC)
C Write absorption coefficient to TRATAN(*,IATM) for each output level IATM
          IF ( TRAFLG ) CALL FLXEFN ( TRATAN(1,IJAC), IPTB )
          GOTO 100
        END IF
      END IF
C
      END
