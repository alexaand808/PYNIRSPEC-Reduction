      SUBROUTINE RADTRA ( ITAN, ISTA, IEND, IDIR, INIT )
C
C VERSION
C     03-MAY-05  AD  Add jdxcon.inc
C     13-JUN-04  AD  Save DSTTAN
C     15-JAN-04  AD  Add IDXDIR instead of IDIR as argument to IDXPTH,JDXPTH
C     01-JAN-04  AD  Add LEV flag.
C     17-DEC-01  AD  Initialise TEMPTB, IPTH, TANHGT
C     05-NOV-01  AD  Add JDXPTH to speed up calculations
C     11-APR-01  AD  Calculate FACTOR allowing for -ve absorber amounts
C     26-MAY-00  AD  Convert to D.P.
C     17-APR-00  AD  Comment change only.
C     02-JAN-00  AD  Allow for GRA option.
C                    Make RADCRV explicitly SNGL. 
C                    Add extra argument to IDXPTH
C     11-AUG-99  AD  Remove redundant SAVE
C     09-AUG-99  AD  Use explicit DBLE ( ) in RADTAN calculations
C     11-APR-99  AD  Use CG weak-limit for BFX option
C     01-FEB-99  AD  Change to use normal path indexing for Jacobian calcs.
C     02-JAN-99  AD  Original. Extracted from RFMRAD.
C
C DESCRIPTION
C     Radiative transfer calculation for single ray
C     Called by RFMRAD for up and down parts for each ray path.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER ITAN  !  I  Index of ray path
      INTEGER ISTA  !  I  Index of starting layer in atmosphere
      INTEGER IEND  !  I  Index of finishing layer in atmosphere
      INTEGER IDIR  !  I  Direction of calculation: -1=down, +1=upwards
      LOGICAL INIT  ! I/O T on input=initialise, set F on exit
C
      EXTERNAL
     &  BBFOPT ! Calculate interp. factor for Planck Function across layer
     &, IDXPTH ! Indexing function for Pth# from Tan.pt#, Atm.Layer# and Gas#
     &, JDXPTH ! Indexing function for Pth# from Tan.pt#, Atm.Layer# and Gas#
     &, LEVUPD ! Set up calculations for intermediate output levels
     &, PLANCK ! Calculate Planck emission spectrum at specified temperature 
     &, PTBATM ! Perturb/unperturb atmospheric profiles for Jacobian calc.
     &, VALGRA ! Function interpolate value from 2D atmospheric field.
      INTEGER IDXPTH, JDXPTH
      REAL    VALGRA
C
C GLOBAL CONSTANTS
      INCLUDE 'jdxcon.inc' ! Indices for non-VMR Jacobians
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data 
      INCLUDE 'crvcom.inc' ! Local radius of curvature
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'jaccom.inc' ! Jacobian data
      INCLUDE 'pthcom.inc' ! RFM path data
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL CONSTANTS
      REAL          SMALL              ! Small no. to avoid divide-by-zero
        PARAMETER ( SMALL = 1.0E-35 )  ! Same as used in GENLN2 EMISSN.FOR
      DOUBLE PRECISION  DSMALL             ! D.P. version of SMALL
        PARAMETER ( DSMALL = 1.0D-35 ) 
C
C LOCAL VARIABLES
      INTEGER   IATM     !  Index of current atmospheric layer on downward path
      INTEGER   IDXDIR   !  Direction index if relevant to path calculation
      INTEGER   IFIN     !  Fine mesh grid point counter (1:NFIN)
      INTEGER   IGAS     !  Absorber counter
      INTEGER   IJAC     !  Counter for Jacobian elements
      INTEGER   INXT     !  Offset Index of 'next' layer for BFX calculations
      INTEGER   IPRV     !  Offset Index of 'prev' layer for BFX calculations
      INTEGER   IPTH     !  Path counter
      INTEGER   JPTH     !  Index of path to be used for scaling absorption
      INTEGER   JTAN     !  Index of original plus Jacobian tan.paths
      INTEGER   KTAN     !  Index of tan.paths used for intermed.level o/p
      INTEGER   KTAN1, KTAN2 ! Range of KTAN values used for current ITAN
      REAL      DSTPRV   !  Saved value of previous DSTTAN
      REAL      DSTTAN   !  Distance from layer boundary to tangent point
      REAL      FACTOR   !  Path scaling factor for Absorption
      REAL      TANHGT   !  Projected tangent altitude [km]
      REAL      TEMLEV   ! Temperature at "near" boundary of layer
      REAL      TEMTAN   ! Tangent point temperature [K]
      LOGICAL   LEVINI   ! T=initialise intermediate output level calculations
      LOGICAL   TEMPTB   ! T=current temperature profile is perturbed.
      DOUBLE PRECISION ABSLAY(MAXFIN) ! Layer optical depth
      DOUBLE PRECISION BBFFIN(MAXFIN) ! Planck emission
      DOUBLE PRECISION BBFLEV(MAXFIN) ! Saved value of previous BBFFIN
      DOUBLE PRECISION BBFWGT(MAXFIN) ! Planck Wgt factors for BFX calculations
      DOUBLE PRECISION SRCLAY(MAXFIN) ! Absorption-wgtd source func. for layer
C
C DATA STATEMENTS
      DATA DSTTAN / 0.0 /     ! Value shouldn't actually be used
      SAVE DSTTAN
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Unnecessary, but avoids warnings over use of uninitialised variables
      IPTH = 0        
      TANHGT = 0.0  

C Note that ISTA=IEND can either mean a homogeneous path calculation or 
C part of either an upward or downward calculation so use IDIR 
      IF ( IDIR .EQ. 1 ) THEN        ! upward calculation (or homog.)
        INXT = 1
        IPRV = 0
      ELSE                           ! downward calculation 
        INXT = 0
        IPRV = 1
      END IF
      IF ( GRAFLG ) THEN
        IDXDIR = IDIR
      ELSE
        IDXDIR = 0
      END IF
C
      TEMPTB = .FALSE.
      IJAC = 0
      JTAN = ITAN
C
C Repeat from here with different JTAN values if performing Jacobian calcs.
  100 CONTINUE
      IF ( INIT ) THEN
        DO IFIN = 1, NFIN
          ABSTAN(IFIN,JTAN) = 0.0D0
          RADTAN(IFIN,JTAN) = 0.0D0
          TRATAN(IFIN,JTAN) = 1.0D0
        END DO
        IF ( LEVFLG ) THEN
          DO IJAC = 1, NJAC
            KTAN = ITNJAC(JTAN,IJAC) 
            IF ( KTAN .NE. 0 ) THEN
              DO IFIN = 1, NFIN
                ABSTAN(IFIN,KTAN) = 0.0D0
                RADTAN(IFIN,KTAN) = 0.0D0
                TRATAN(IFIN,KTAN) = 1.0D0
              END DO
            END IF
          END DO
          LEVINI = .TRUE.        ! Reset by LEVUPD
        END IF
C
C If modelling BBF v. Optical depth, initialise (NB BFXFLG excludes HOMFLG)
        IF ( BFXFLG ) THEN
          IF ( NADFLG .OR. ZENFLG ) THEN
            TANHGT = 0.0
            DSTTAN = 0.0
          ELSE 
            TANHGT = MIN ( HGTTAN(ITAN), GEOTAN(ITAN) )
            DSTTAN = SQRT ( 2.0 * ( SNGL(RADCRV) + TANHGT ) 
     &                 * MAX ( ( HGTATM(ISTA+IPRV) - TANHGT ), 0.0 ) ) 
          END IF 
        END IF
      END IF
C
      DO IATM = ISTA, IEND, IDIR
        DO IFIN = 1, NFIN
          SRCLAY(IFIN) = 0.0D0
          ABSLAY(IFIN) = 0.0D0
        END DO
        DO IGAS = 1, NGAS
          IPTH = JDXPTH ( JTAN, IATM, IGAS, IDXDIR )
C If there isn't a Jacobian-perturbed path use original unperturbed path
          IF ( IPTH .EQ. 0 ) IPTH = IDXPTH ( ITAN, IATM, IGAS, IDXDIR )
          IF ( CLCPTH(IPTH) ) THEN
            JPTH = IPTH
            FACTOR = 1.0
          ELSE
            JPTH = ISCPTH(IPTH)
            FACTOR = AMTPTH(IPTH) / 
     &        SIGN ( MAX ( ABS(AMTPTH(JPTH)), SMALL ), AMTPTH(JPTH) )
          END IF
          IF ( NTEFLG ) THEN
            CALL PLANCK ( TEMPTH(IPTH), NFIN, WNOFIN, BBFFIN )
            DO IFIN = 1, NFIN
              SRCLAY(IFIN) = SRCLAY(IFIN) + BBFFIN(IFIN) * 
     &                         DBLE ( CNTFIN(IFIN,JPTH) * FACTOR )
              ABSLAY(IFIN) = ABSLAY(IFIN) + 
     &                         DBLE ( ABSFIN(IFIN,JPTH) * FACTOR )
            END DO
          ELSE
            CALL PLANCK ( TEMPTH(IPTH), NFIN, WNOFIN, BBFFIN )
            DO IFIN = 1, NFIN  
              SRCLAY(IFIN) = SRCLAY(IFIN) + BBFFIN(IFIN) * 
     &                         DBLE ( ABSFIN(IFIN,JPTH) * FACTOR )
              ABSLAY(IFIN) = ABSLAY(IFIN) + 
     &                         DBLE ( ABSFIN(IFIN,JPTH) * FACTOR )
            END DO
          END IF
        END DO
C
        IF ( BFXFLG ) THEN
C In tangent layer, need to interpolate for tangent point temperature
C IATM=IATTAN only happens for limb-viewing since IATTAN=0 for ZEN/NAD
C For surface-intersecting rays, IATTAN=1 and HGTTAN = HGTATM(1) so FACTOR=0
C and TEMATM(1) gets used anyway
          IF ( IATM .EQ. IATTAN(ITAN) .AND. IDIR .EQ. +1 ) THEN
            IF ( GRAFLG ) THEN
              TEMTAN = VALGRA ( HGTTAN(ITAN), PSITAN(ITAN), 'TEM', 0 )
            ELSE
              FACTOR = ( HGTTAN(ITAN) - HGTATM(IATM) ) / 
     &                 ( HGTATM(IATM+1)- HGTATM(IATM) )
              TEMTAN = TEMATM(IATM) + FACTOR * 
     &                 ( TEMATM(IATM+1) - TEMATM(IATM) ) 
            END IF
            CALL PLANCK ( TEMTAN, NFIN, WNOFIN, BBFLEV )
          ELSE 
            IF ( GRAFLG ) THEN          ! PSIPTH contains "near" value of psi
              TEMLEV = 
     &          VALGRA ( HGTATM(IATM+IPRV), PSIPTH(IPTH), 'TEM', 0 )
            ELSE
              TEMLEV = TEMATM(IATM+IPRV)
            END IF
            CALL PLANCK ( TEMLEV, NFIN, WNOFIN, BBFLEV )
          END IF
          DSTPRV = DSTTAN
          IF ( NADFLG .OR. ZENFLG ) THEN
            DSTTAN = 0.0
          ELSE
            DSTTAN = SQRT ( 2.0 * ( SNGL(RADCRV) + TANHGT ) 
     &               * MAX ( ( HGTATM(IATM+INXT) - TANHGT ), 0.0 ) ) 
          END IF
          CALL BBFOPT ( DSTPRV, DSTTAN, NFIN, ABSLAY, BBFWGT )
          DO IFIN = 1, NFIN
            SRCLAY(IFIN) = SRCLAY(IFIN) /    
     &        SIGN ( MAX ( ABS(ABSLAY(IFIN)), DSMALL ), ABSLAY(IFIN)) 
            RADTAN(IFIN,JTAN) = RADTAN(IFIN,JTAN) + 
     &        TRATAN(IFIN,JTAN) 
     &        *  ( BBFLEV(IFIN) + 
     &             BBFWGT(IFIN) * ( SRCLAY(IFIN) - BBFLEV(IFIN) ) ) 
     &        * ( 1.0D0 - EXP ( -ABSLAY(IFIN) ) )
            ABSTAN(IFIN,JTAN) = ABSTAN(IFIN,JTAN) + ABSLAY(IFIN)
            TRATAN(IFIN,JTAN) = EXP ( -ABSTAN(IFIN,JTAN) )
          END DO
        ELSE
          DO IFIN = 1, NFIN
            SRCLAY(IFIN) = SRCLAY(IFIN) /    
     &        SIGN ( MAX ( ABS(ABSLAY(IFIN)), DSMALL ), ABSLAY(IFIN)) 
            RADTAN(IFIN,JTAN) = RADTAN(IFIN,JTAN) +
     &        TRATAN(IFIN,JTAN) * SRCLAY(IFIN)
     &        * ( 1.0D0 - EXP ( - ABSLAY(IFIN) ) ) 
            ABSTAN(IFIN,JTAN) = ABSTAN(IFIN,JTAN) + ABSLAY(IFIN)
            TRATAN(IFIN,JTAN) = EXP ( - ABSTAN(IFIN,JTAN) )
          END DO
        END IF
C Update summations for intermediate output levels starting at various points
C along tangent path
        IF ( LEVFLG ) THEN
          CALL LEVUPD ( LEVINI, ITAN, IATM, IDIR, KTAN1, KTAN2 )
          DO KTAN = KTAN1, KTAN2
            DO IFIN = 1, NFIN
              IF ( BFXFLG ) THEN
                RADTAN(IFIN,KTAN) = RADTAN(IFIN,KTAN) + 
     &          TRATAN(IFIN,KTAN) 
     &          *  ( BBFLEV(IFIN) + 
     &                BBFWGT(IFIN) * ( SRCLAY(IFIN) - BBFLEV(IFIN) ) ) 
     &          * ( 1.0D0 - EXP ( -ABSLAY(IFIN) ) )
              ELSE
                RADTAN(IFIN,KTAN) = RADTAN(IFIN,KTAN) +
     &            TRATAN(IFIN,KTAN) * SRCLAY(IFIN) 
     &            * ( 1.0D0 - EXP ( - ABSLAY(IFIN) ) ) 
              END IF
              ABSTAN(IFIN,KTAN) = ABSTAN(IFIN,KTAN) + ABSLAY(IFIN)
              TRATAN(IFIN,KTAN) = EXP ( -ABSTAN(IFIN,KTAN) )
            END DO
          END DO
        END IF
C
      END DO
C
C Last part sets up perturbed profiles for repeated calcs for Jacobians
C (each call undoes the perturbations set up by the previous call)
C If JTAN=0 this means zero Jacobian for ITAN,IJAC combination, so skip.
      IF ( JACFLG ) THEN 
        DO WHILE ( IJAC .LT. NJAC )
          IJAC = IJAC + 1
          JTAN = ITNJAC(ITAN,IJAC)
          IF ( JTAN .NE. 0 ) THEN
            IF ( IGSJAC(IJAC) .EQ. MAXGAS+JDXTEM .AND. BFXFLG ) THEN
              CALL PTBATM ( IJAC )
              TEMPTB = .TRUE.
            ELSE IF ( TEMPTB ) THEN
              CALL PTBATM ( NJAC+1 )       ! Undo last temperature perturbation
            END IF
            GOTO 100                       ! Repeat calculation for ITAN=JTAN
          END IF
        END DO
C If last ptb was a temperature perturbation, this needs to be undone since
C actual temperatures required for Planck functions in radiative transfer.
C Other (VMR, PRE) perturbations can be left until next call of PTBATM.
        IF ( TEMPTB ) CALL PTBATM ( NJAC+1 ) ! Undo temperature ptb.
      END IF
C
      INIT = .FALSE.
C
      END 
