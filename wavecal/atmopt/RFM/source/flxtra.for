      SUBROUTINE FLXTRA ( RADFLX, ABSFLX, ITAN, IATM, IDIR, 
     &                    SECFAC, SRCFAC )
C
C VERSION
C     17-JAN-04  AD  Use 0 as direction argument to IDXPTH
C     26-MAY-00  AD  Convert from S.P. to D.P.
C     23-APR-00  AD  Original.
C
C DESCRIPTION
C     Radiative transfer calculations for flux mode.
C     Called by RFMFLX for each layer in atmosphere.
C     Assumes that vertical path information is loaded into pthcom.inc.
C
      IMPLICIT NONE
C
C ARGUMENTS
      DOUBLE PRECISION RADFLX(*) ! I/O Radiance
      DOUBLE PRECISION ABSFLX(*) ! I/O Integrated optical thickness
      INTEGER IATM      !  I  Index of atmospheric layer 
      INTEGER ITAN      !  I  Index of unperturbed/perturbed ray path
      INTEGER IDIR      !  I  Direction of calculation: -1=down, +1=upwards
      REAL    SECFAC    !  I  Sec(theta) scale factor for absorber amount
      REAL    SRCFAC    !  I  Scale factor for source function
C
      EXTERNAL
     &  BBFOPT ! Calculate interp. factor for Planck Function across layer
     &, IDXPTH ! Indexing function for Pth# from Tan.pt#, Atm.Layer# and Gas#
     &, PLANCK ! Calculate Planck emission spectrum at specified temperature 
      INTEGER IDXPTH
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data 
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'pthcom.inc' ! RFM path data
C
C LOCAL CONSTANTS
      REAL          SMALL             ! Small no. to avoid divide-by-zero
        PARAMETER ( SMALL = 1.0E-35 ) ! Same as used in GENLN2 EMISSN.FOR
      DOUBLE PRECISION DSMALL         ! D.P. version of SMALL
        PARAMETER ( DSMALL = 1.0D-35 ) 
C
C LOCAL VARIABLES
      INTEGER   IFIN     !  Fine mesh grid point counter (1:NFIN)
      INTEGER   IGAS     !  Absorber counter
      INTEGER   INEAR    !  level index offset for "near" boundary (BFX only)
      INTEGER   IPTH     !  Path counter
      INTEGER   JPTH     !  Secondary Path counter
      REAL      AMTFAC   !  Net scaling factor for absorber amount
      DOUBLE PRECISION BBFLEV(MAXFIN) ! Saved value of previous BBFFIN
      DOUBLE PRECISION BBFWGT(MAXFIN) ! Planck Wgt factors for BFX calculations
      DOUBLE PRECISION ABSLAY(MAXFIN) ! Layer optical depth
      DOUBLE PRECISION BBFFIN(MAXFIN) ! Planck emission
      DOUBLE PRECISION SRCLAY(MAXFIN) ! Absorption-weighted source function for layer
      REAL      TEMLEV   ! Temperature at "near" boundary of layer
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      IF ( IDIR .EQ. 1 ) THEN        ! upward flux calculation 
        INEAR = 1
      ELSE                           ! downward flux calculation 
        INEAR = 0
      END IF
C
C Initialise layer optical thickness and source functions
      DO IFIN = 1, NFIN
        SRCLAY(IFIN) = 0.0D0
        ABSLAY(IFIN) = 0.0D0
      END DO
C
      DO IGAS = 1, NGAS
        IPTH = IDXPTH ( ITAN, IATM, IGAS, 0 )
C If there isn't a Jacobian-perturbed path use original unperturbed path
        IF ( IPTH .EQ. 0 ) IPTH = IDXPTH ( 1, IATM, IGAS, 0 ) 
C All unperturbed paths should be explicitly calculated, but some Jacobian
C paths may be scaled by the absorber amount+perturbation
        IF ( CLCPTH(IPTH) ) THEN
          AMTFAC = SECFAC
        ELSE
          JPTH = IPTH
          IPTH = ISCPTH(JPTH)
          AMTFAC = SECFAC * AMTPTH(JPTH)/(MAX(AMTPTH(IPTH),SMALL))
        END IF
C
        CALL PLANCK ( TEMPTH(IPTH), NFIN, WNOFIN, BBFFIN )
        IF ( NTEFLG ) THEN
          DO IFIN = 1, NFIN
            SRCLAY(IFIN) = SRCLAY(IFIN) + BBFFIN(IFIN) *
     &                       DBLE ( CNTFIN(IFIN,IPTH) * AMTFAC )
            ABSLAY(IFIN) = ABSLAY(IFIN) + 
     &                       DBLE ( ABSFIN(IFIN,IPTH) * AMTFAC )
          END DO
        ELSE
          DO IFIN = 1, NFIN  
            SRCLAY(IFIN) = SRCLAY(IFIN) + BBFFIN(IFIN) *
     &                       DBLE ( ABSFIN(IFIN,IPTH) * AMTFAC )
            ABSLAY(IFIN) = ABSLAY(IFIN) + 
     &                       DBLE ( ABSFIN(IFIN,IPTH) * AMTFAC )
          END DO
        END IF
      END DO
      IF ( BFXFLG ) THEN
        TEMLEV = TEMATM(IATM+INEAR)
        CALL PLANCK ( TEMLEV, NFIN, WNOFIN, BBFLEV )
        CALL BBFOPT ( 0.0, 0.0, NFIN, ABSLAY, BBFWGT )
        DO IFIN = 1, NFIN
          SRCLAY(IFIN) = SRCLAY(IFIN) * DBLE ( SRCFAC ) /    
     &      SIGN ( MAX ( ABS(ABSLAY(IFIN)), DSMALL ), ABSLAY(IFIN)) 
          RADFLX(IFIN) = RADFLX(IFIN) * EXP ( - ABSLAY(IFIN) ) + 
     &      ( BBFLEV(IFIN) + 
     &          BBFWGT(IFIN) * ( SRCLAY(IFIN) - BBFLEV(IFIN) ) )
     &      * ( 1.0D0 - EXP ( - ABSLAY(IFIN) ) )
          ABSFLX(IFIN) = ABSFLX(IFIN) + ABSLAY(IFIN) 
        END DO
      ELSE
        DO IFIN = 1, NFIN
          SRCLAY(IFIN) = SRCLAY(IFIN) * DBLE ( SRCFAC ) /    
     &      SIGN ( MAX ( ABS(ABSLAY(IFIN)), DSMALL ), ABSLAY(IFIN)) 
          RADFLX(IFIN) = RADFLX(IFIN) * EXP ( - ABSLAY(IFIN) ) +
     &      SRCLAY(IFIN) * ( 1.0D0 - EXP ( - ABSLAY(IFIN) ) )
          ABSFLX(IFIN) = ABSFLX(IFIN) + ABSLAY(IFIN)
        END DO
      END IF
C
      END 
