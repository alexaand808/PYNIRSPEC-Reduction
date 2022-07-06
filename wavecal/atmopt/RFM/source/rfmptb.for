      SUBROUTINE RFMPTB ( FAIL, ERRMSG )
C
C VERSION
C     03-MAY-05  AD  Allow for vibrational temperature Jacobians.
C                    Add jdxcon.inc
C     03-SEP-03  AD  Allow for surface parameter Jacobians
C     19-FEB-03  AD  Bug#42 Use KDIR
C     04-JAN-02  AD  Bug#34 Allow for JPTH itself being scaled
C     08-JUN-01  AD  Check LTAN v MAXTAN allowing for LOS flag
C     12-APR-01  AD  Check for ABS(AMTPTH) when comparing changes
C     22-MAR-01  AD  Add PPASAV, RAYSAV
C     21-APR-00  AD  Allow for FLX flag. 
C                    Use scaled paths where possible for extra Jac.paths
C     30-DEC-99  AD  Add IDRPTH and remove unused IDIPTH
C     31-JUL-99  AD  Remove unused variable SIGDIF
C     04-FEB-99  AD  Store perturbations in *PTH rather than *PTB.
C     03-JAN-99  AD  Original.
C
C DESCRIPTION
C     Set up list of pertubations required for RFM Jacobian calc.
C     Called once by RFM if JAC option enabled.
C
      IMPLICIT NONE
C
C ARGUMENTS
      LOGICAL      FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG !  O  Error message written if FAIL is TRUE 
C
      EXTERNAL
     &  HOMPTH ! Set up homogeneous paths for RFM calculations.
     &, IDXPTH ! Indexing function for Pth# from Tan.pt#, Atm.Layer# and Gas#
     &, NTEPTH ! Determine paths affected by Non-LTE Jacobian.
     &, PTBATM ! Perturb/unperturb atmospheric profiles for Jacobian calc.
     &, PTHSRT ! Sort Calculated paths to start of /PTHCOM/ list.
     &, RFMLOG ! Write text message to RFM log file.
     &, RFMNAD ! Calculate nadir/zenith path parameters through atmosphere
     &, RFMRAY ! Construct ray paths through atmosphere
      INTEGER IDXPTH
C
C GLOBAL CONSTANTS
      INCLUDE 'jdxcon.inc' ! Indices for non-VMR Jacobians
      INCLUDE 'rfmcon.inc' ! Constants for RFM calculations
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'jaccom.inc' ! Jacobian data
      INCLUDE 'pthcom.inc' ! RFM path data
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL CONSTANTS
      REAL          SIGCLC ! Fraction of p,T ptb requiring new abs.coeff.calc
        PARAMETER ( SIGCLC = 0.1 )
      REAL          SIGPTH ! Fraction of p,T,v ptb requiring new path segment.
        PARAMETER ( SIGPTH = 0.01 )
C
C LOCAL VARIABLES
      INTEGER IJAC   ! Jacobian element counter
      INTEGER IPTH, JPTH, KPTH   ! Path indices
      INTEGER ITAN, JTAN, KTAN   ! Tangent ray indices
      INTEGER KDIR   ! Direction of scaled path
      INTEGER LCLC   ! Saved value of NCLC from entry
      INTEGER LPTH   ! Saved value of NPTH from entry
      LOGICAL NTEJAC ! T=Vib.Temp Jacobian
      LOGICAL SFCJAC ! T=Surface parameter jacobian, F=atmospheric parameter
      REAL    AMTSAV(MAXPTH) ! Saved values of absorber amounts in paths
      REAL    PPASAV(MAXPTH) ! Saved values of pressures in paths
      REAL    PRESAV(MAXPTH) ! Saved values of pressures in paths
      REAL    RAYSAV(MAXPTH) ! Saved values of pressures in paths
      REAL    TEMSAV(MAXPTH) ! Saved values of CG path temperatures
      CHARACTER*80 MESSGE    ! Text message for log file
C
C EXECUTABLE CODE ------------------------------------------------------------
C
      LCLC = NCLC
      LPTH = NPTH
      DO IPTH = 1, NPTH
        AMTSAV(IPTH) = AMTPTH(IPTH)
        PPASAV(IPTH) = PPAPTH(IPTH)
        PRESAV(IPTH) = PREPTH(IPTH)
        RAYSAV(IPTH) = RAYPTH(IPTH)
        TEMSAV(IPTH) = TEMPTH(IPTH)
      END DO
C
C At some stage, either here (No FOV), or after RFMFOV, Jacobian paths will
C be assigned sequentially to tan.paths in /TANCOM/ - check enough space.
C Allow for any additional paths required for LOS flag (only assigned by
C RFMLOS immediately prior to output, so not included in LTAN)
      LTAN = NTAN * ( 1 + NJAC )
      IF ( LTAN + LOSTAN .GT. MAXTAN ) THEN
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,I10,A,I8)' ) 
     &    'F-RFMPTB: Jac.Calcs require', LTAN + LOSTAN, 
     &    ' tan.pths > MAXTAN (RFMSIZ.INC)=', MAXTAN
        RETURN
      END IF
C 
C For each Jacobian element, perturb atmosphere and recalculate paths
C For limbview with FOV, tangent paths are assigned later if any significant
C perturbation in atmospheric path is found from Jacobians. However, surface
C parameter perturbations won't show up so assign spare paths here.
      JTAN = MTAN
      DO IJAC = 1, NJAC
        IF ( FOVFLG ) THEN     ! No spare tan.paths assigned for IJAC calc.yet
          SFCJAC = IGSJAC(IJAC) .EQ. MAXGAS + JDXSFT .OR. ! Surface temperature
     &             IGSJAC(IJAC) .EQ. MAXGAS + JDXSFE      ! Surface emissivity
          DO ITAN = 1, MTAN       
            IF ( SFCJAC .AND. SFCTAN(ITAN) ) THEN ! ray intersects surface
              IF ( JTAN .EQ. MAXTAN ) THEN
                FAIL = .TRUE.
                WRITE ( ERRMSG, '(A,I8)' )
     &            'F-RFMPTB: Jac.Calcs require > MAXTAN '//
     &            '(RFMSIZ.INC) Tan.paths, MAXTAN=', MAXTAN
                RETURN
              END IF
              JTAN = JTAN + 1
              ITNJAC(ITAN,IJAC) = JTAN
              ITNJAC(JTAN,IJAC) = ITAN
            ELSE 
              ITNJAC(ITAN,IJAC) = 0      ! assign later if required
            END IF
          END DO
        ELSE
          DO ITAN = 1, MTAN    ! Assign tan.paths for each ITAN,IJAC 
            JTAN = JTAN + 1    ! Should end up with JTAN=LTAN
            ITNJAC(ITAN,IJAC) = JTAN
            ITNJAC(JTAN,IJAC) = ITAN
          END DO
        END IF
C
        NTEJAC = IGSJAC(IJAC) .GT. MAXGAS + JDXNTE
        IF ( NTEJAC ) CALL NTEPTH ( IJAC, LPTH )
C
        CALL PTBATM ( IJAC )      
        IF ( HOMFLG ) THEN
          CALL HOMPTH ( FAIL, ERRMSG )   ! Previously tested, so shouldn't fail
          IF ( FAIL ) STOP 'F-RFMPTB: Logical Error#1'
        ELSE IF ( NADFLG .OR. ZENFLG .OR. FLXFLG ) THEN
          CALL RFMNAD
        ELSE IF ( GRAFLG ) THEN
          CALL RFMGRA
        ELSE 
          CALL RFMRAY
        END IF
C
C Check over all paths to see which ones have changed significantly
C
        DO IPTH = 1, LPTH
          IF ( ( ABS ( PREPTH(IPTH) - PRESAV(IPTH) )            ! New PTH reqd.
     &           .GT. SIGPTH * PTBPRE * PRESAV(IPTH) ) .OR. 
     &         ( ABS ( TEMPTH(IPTH) - TEMSAV(IPTH) ) 
     &           .GT. SIGPTH * PTBTEM ) .OR.
     &         ( ABS ( AMTPTH(IPTH) - AMTSAV(IPTH) ) 
     &           .GT. SIGPTH * PTBVMR * ABS(AMTSAV(IPTH)) ) .OR.
     &         ( NTEJAC .AND. IVJPTH(IPTH) .GT. 0 ) ) THEN 
            IF ( NPTH .EQ. MAXPTH ) THEN
              FAIL = .TRUE.
              WRITE ( ERRMSG, '(A,I11)' )
     &          'F-RFMPTB: No. Paths for Jacobians > '//
     &          'MAXPTH in RFMSIZ.INC, =', MAXPTH
              RETURN
            END IF
            NPTH = NPTH + 1
            IATPTH(NPTH) = IATPTH(IPTH)
            IDRPTH(NPTH) = IDRPTH(IPTH)
            IGSPTH(NPTH) = IGSPTH(IPTH)
            IVJPTH(NPTH) = IVJPTH(IPTH)
            NLNPTH(NPTH) = NLNPTH(IPTH)
            AMTPTH(NPTH) = AMTPTH(IPTH)
            PPAPTH(NPTH) = PPAPTH(IPTH)
            PREPTH(NPTH) = PREPTH(IPTH)
            PSIPTH(NPTH) = PSIPTH(IPTH)
            RAYPTH(NPTH) = RAYPTH(IPTH)
            TEMPTH(NPTH) = TEMPTH(IPTH)
C If this is a new path, check if tangent path exists for ITAN,IJAC radiative 
C transfer calculations. If not, establish one as JTAN.
            ITAN = ITNPTH(IPTH)       ! Should only return 1 with FLX flag
            IF ( ITNJAC(ITAN,IJAC) .EQ. 0 ) THEN   ! should only be if FOV flag
              IF ( .NOT. FOVFLG ) STOP 'F-RFMPTB: Logical Error#2'
              IF ( JTAN .EQ. MAXTAN ) THEN
                FAIL = .TRUE.
                WRITE ( ERRMSG, '(A,I8)' )
     &            'F-RFMPTB: Jac.Calcs require > MAXTAN '//
     &            '(RFMSIZ.INC) Tan.paths, MAXTAN=', MAXTAN
                RETURN
              END IF
              JTAN = JTAN + 1
              ITNJAC(ITAN,IJAC) = JTAN
              ITNJAC(JTAN,IJAC) = ITAN
            END IF
C Flux calculations use separate ITAN assignments for output paths and 
C calculated/Jacobian paths
            IF ( FLXFLG ) THEN
              ITNPTH(NPTH) = 1 + IJAC
            ELSE
              ITNPTH(NPTH) = ITNJAC(ITAN,IJAC)
            END IF
C Check if this can be scaled to unperturbed path (ie only absorber amount
C changes significantly) - doesn't apply for VT Jacobians since pT unchanged
            IF ( ABS ( PREPTH(IPTH) - PRESAV(IPTH) ) .LT.
     &                SIGCLC * PTBPRE * PRESAV(IPTH) .AND.
     &           ABS ( TEMPTH(IPTH) - TEMSAV(IPTH) ) .LT.
     &                             SIGCLC * PTBTEM   .AND.
     &           .NOT. NTEJAC ) THEN 
              CLCPTH(NPTH) = .FALSE.
              ISCPTH(NPTH) = ISCPTH(IPTH)
            ELSE IF ( .NOT. CLCPTH(IPTH) ) THEN
C If not, then if unperturbed path was scaled, scale to equivalent ptb path.
C Also applies to VT Jacobians
              KPTH = ISCPTH(IPTH)       ! Index of unperturbed calc. path
              KTAN = ITNPTH(KPTH)       ! Index of KPTH tangent ray
              KDIR = IDRPTH(KPTH)
              JPTH = IDXPTH ( ITNJAC(KTAN,IJAC), IATPTH(IPTH), 
     &                        IGSPTH(IPTH), KDIR )  ! ptb calc.path
              IF ( JPTH .EQ. 0 ) STOP 'F-RFMPTB: Logical Error#3'
              ISCPTH(NPTH) = ISCPTH(JPTH)  ! Allow for JPTH itself being scaled
              CLCPTH(NPTH) = .FALSE.
C If not, then set up a new calculated path
            ELSE
              IF ( NCLC .EQ. MAXCLC ) THEN
                FAIL = .TRUE.
                WRITE ( ERRMSG, '(A,I11)' )
     &            'F-RFMPTB: No.Calc.Paths for Jacobians > '//
     &            'MAXCLC in RFMSIZ.INC, =', MAXCLC
                RETURN
              END IF
              NCLC = NCLC + 1
              ISCPTH(NPTH) = NPTH
              CLCPTH(NPTH) = .TRUE.
            END IF
          END IF
          AMTPTH(IPTH) = AMTSAV(IPTH)
          PPAPTH(IPTH) = PPASAV(IPTH)
          PREPTH(IPTH) = PRESAV(IPTH)
          RAYPTH(IPTH) = RAYSAV(IPTH)
          TEMPTH(IPTH) = TEMSAV(IPTH)
          IVJPTH(IPTH) = 0
        END DO
      END DO
C
C Construct info message detailing no.tangent paths required
      WRITE ( MESSGE, '(A,I7,A,I7)' ) 
     &  'I-RFMPTB: ', MAX ( LTAN, JTAN ) - MTAN, 
     &  ' extra tan.paths reqd for Jac.Calcs. New total=',MAX(LTAN,JTAN)
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      WRITE ( MESSGE, '(A,I7,A,I7)' ) 'I-RFMPTB: ', NCLC - LCLC,
     &  ' extra Clc.paths reqd for Jac.Calcs. New total=', NCLC
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      WRITE ( MESSGE, '(A,I7,A,I7)' ) 'I-RFMPTB: ', NPTH - LPTH,
     &  ' extra C-G.paths reqd for Jac.Calcs. New total=', NPTH
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
C Reset atmosphere to unperturbed values
      CALL PTBATM ( NJAC+1 )
C
C Resort paths to set 1:NCLC as calculated paths
      CALL PTHSRT
C
      FAIL = .FALSE.
C
      END
