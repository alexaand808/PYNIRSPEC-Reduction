      SUBROUTINE WRTOUT ( TYPE, FAIL, ERRMSG )
C
C VERSION
C     14-SEP-11  AD  Add RJT option.
C     12-OCT-05  AD  Bug fix: correct BBT Jacobians
C     09-AUG-03  AD  Add BBT option
C     20-JUL-02  AD  Revise handling of LUNFLG option
C     17-DEC-01  AD  Initialise LUN
C     08-JUN-01  AD  Output 1:LTAN+LOSTAN instead of 1:LTAN 
C     16-JUN-00  AD  Add factor 1e-5 to convert flux rads from nW/cm2 to W/m2
C     09-JUN-00  AD  Rewritten.
C                    Allow for DBL flag
C                    Change DP TRA output from (X,3F20.16) to *
C     22-APR-00  AD  Allow for COO flag
C     11-AUG-99  AD  Remove redundant fincom.inc
C     21-JUN-99  AD  Bug in use of LUNFLG - writing to wrong files
C     02-FEB-99  AD  Add LUNOUT, FILOUT
C     02-JAN-99  AD  Change 1:NTAN to 1:LTAN to allow for any Jacobian spectra.
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Write RFM output spectra for current widemesh interval.
C     Called by RFM for each widemesh interval.
C
C     NB: When used with the LUN Flag, this routine reopens the existing files
C     to APPEND the new data, which is not standard F77 
C     the form of the open-statement may vary from one system to another.
C
      IMPLICIT NONE
C
C ARGUMENTS
      CHARACTER*3  TYPE    !  I  Type of spectrum 'ABS','COO',etc
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error occurs
      CHARACTER*80 ERRMSG  !  O  Error message written if FAIL is TRUE
C
C GLOBAL CONSTANTS
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'outcom.inc' ! RFM output file data
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL VARIABLES
      INTEGER      IFIL    ! Output file counter
      INTEGER      IOUT    ! Output point counter
      INTEGER      IOS     ! Value of IOSTAT saved for error messages.
      INTEGER      ITAN    ! Tangent height counter
      INTEGER      LUN     ! LUN for transmission output files
      INTEGER      NOMTAN  ! Nominal tangent path (of Jacobian tan.paths)
      INTEGER      NOUT    ! No.pts written for each WM interval
      CHARACTER*11 FMTTYP  ! Format type, 'FORMATTED' or 'UNFORMATTED'
      DOUBLE PRECISION BBTNOM  ! BBT of nominal (unperturbed) radiance
      DOUBLE PRECISION BBTOUT(MAXFIN) ! Brightness temperature spectrum 
      DOUBLE PRECISION BBTPTB  ! BBT of perturbed radiance RADPTB
      DOUBLE PRECISION RADFAC  ! Factor to convert nW/cm2 to W/m2 for Flux rad
      DOUBLE PRECISION RADPTB  ! Total radiance including Jacobian perturbation
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      FAIL = .FALSE.
      NOUT = IOUT2 - IOUT1 + 1
      IF ( BINFLG ) THEN
        FMTTYP = 'UNFORMATTED'
      ELSE
        FMTTYP = 'FORMATTED'
      END IF
C
      IF ( FLXFLG .AND. .NOT. VRTFLG ) THEN
        RADFAC = 1.0D-5
      ELSE
        RADFAC = 1.0D0
      END IF
C
      LUN = 0
      IF ( TYPE .EQ. 'ABS' ) LUN = LUNABS
      IF ( TYPE .EQ. 'BBT' ) LUN = LUNBBT
      IF ( TYPE .EQ. 'COO' ) LUN = LUNCOO
      IF ( TYPE .EQ. 'OPT' ) LUN = LUNOPT
      IF ( TYPE .EQ. 'RAD' ) LUN = LUNRAD
      IF ( TYPE .EQ. 'RJT' ) LUN = LUNRJT
      IF ( TYPE .EQ. 'TRA' ) LUN = LUNTRA
      IF ( LUN .EQ. 0 ) STOP 'F-WRTOUT: Logical Error#1'
      IF ( LUNFLG ) THEN
        IF ( TYPE .EQ. 'ABS' ) IFIL = IFLABS
        IF ( TYPE .EQ. 'BBT' ) IFIL = IFLBBT
        IF ( TYPE .EQ. 'COO' ) IFIL = IFLCOO
        IF ( TYPE .EQ. 'OPT' ) IFIL = IFLOPT
        IF ( TYPE .EQ. 'RAD' ) IFIL = IFLRAD
        IF ( TYPE .EQ. 'RJT' ) IFIL = IFLRJT
        IF ( TYPE .EQ. 'TRA' ) IFIL = IFLTRA
      END IF
C
      DO ITAN = 1, LTAN + LOSTAN
        IF ( LUNFLG ) THEN
          IFIL = IFIL + 1
          OPEN ( UNIT=LUN, FILE=FILOUT(IFIL), STATUS='OLD',
     &               ACCESS='APPEND',                ! non-standard F77
C    &               POSITION='APPEND',              ! Standard F90
     &               FORM=FMTTYP, IOSTAT=IOS, ERR=900 )
        END IF
C
        IF ( TYPE .EQ. 'ABS' ) THEN
          IF ( BINFLG ) THEN
            IF ( DBLFLG ) THEN
              WRITE ( LUN, IOSTAT=IOS, ERR=900 ) NOUT, 
     &          ( 1.D0 - TRATAN(IOUT,ITAN), IOUT = IOUT1, IOUT2 )
            ELSE
              WRITE ( LUN, IOSTAT=IOS, ERR=900 ) NOUT, 
     &          ( SNGL(1.D0-TRATAN(IOUT,ITAN)), IOUT = IOUT1, IOUT2 )
            END IF
          ELSE   
            IF ( DBLFLG ) THEN
              WRITE ( LUN, *, IOSTAT=IOS, ERR=900 )
     &          ( 1.D0 - TRATAN(IOUT,ITAN), IOUT = IOUT1, IOUT2 )
            ELSE
              WRITE ( LUN, *, IOSTAT=IOS, ERR=900 )
     &          ( SNGL(1.D0-TRATAN(IOUT,ITAN)), IOUT = IOUT1, IOUT2 )
            END IF
          END IF
        ELSE IF ( TYPE .EQ. 'BBT' ) THEN ! Convert radiance to Bright.Temp.
C If calculating Jacobians (ITAN>NTAN) then need to convert radiance 
C differences back to absolute radiance before calculating BBT
          IF ( ITAN .LE. NTAN ) THEN          ! Nominal spectra
            DO IOUT = IOUT1, IOUT2
              IF ( RADTAN(IOUT,ITAN) .GT. ARGMIN ) THEN
                BBTOUT(IOUT) = C2 * WNOFIN(IOUT) /  ! C1,C2 in phycon.inc
     &            LOG ( 1.0D0 + C1*WNOFIN(IOUT)**3 / RADTAN(IOUT,ITAN) ) 
              ELSE
                BBTOUT(IOUT) = 0.0D0
              END IF
            END DO
          ELSE                                ! Jacobian spectra
            NOMTAN = MOD ( ITAN, NTAN )             
            IF ( NOMTAN .EQ. 0 ) NOMTAN = NTAN
            DO IOUT = IOUT1, IOUT2
              IF ( RADTAN(IOUT,NOMTAN) .GT. ARGMIN ) THEN
                RADPTB = RADTAN(IOUT,ITAN) + RADTAN(IOUT,NOMTAN)
                BBTNOM = C2 * WNOFIN(IOUT) / 
     &            LOG ( 1.0D0 + C1*WNOFIN(IOUT)**3/RADTAN(IOUT,NOMTAN)) 
                BBTPTB = C2 * WNOFIN(IOUT) / 
     &            LOG ( 1.0D0 + C1*WNOFIN(IOUT)**3 / RADPTB ) 
                BBTOUT(IOUT) = BBTPTB - BBTNOM
              ELSE
                BBTOUT(IOUT) = 0.0D0
              END IF
            END DO
          END IF
          IF ( BINFLG ) THEN
            IF ( DBLFLG ) THEN
              WRITE ( LUN, IOSTAT=IOS, ERR=900 ) NOUT, 
     &          ( BBTOUT(IOUT), IOUT = IOUT1, IOUT2 )
            ELSE
              WRITE ( LUN, IOSTAT=IOS, ERR=900 ) NOUT, 
     &          ( SNGL(BBTOUT(IOUT)), IOUT = IOUT1, IOUT2 )
            END IF
          ELSE   
            IF ( DBLFLG ) THEN
              WRITE ( LUN, *, IOSTAT=IOS, ERR=900 )
     &          ( BBTOUT(IOUT), IOUT = IOUT1, IOUT2 )
            ELSE
              WRITE ( LUN, *, IOSTAT=IOS, ERR=900 )
     &          ( SNGL(BBTOUT(IOUT)), IOUT = IOUT1, IOUT2 )
            END IF
          END IF
        ELSE IF ( TYPE .EQ. 'COO' .OR. TYPE .EQ. 'OPT' ) THEN    
C ! NB: Cooling rates stored in ABSTAN instead of optical depth
          IF ( BINFLG ) THEN
            IF ( DBLFLG ) THEN
              WRITE ( LUN, IOSTAT=IOS, ERR=900 ) NOUT, 
     &          ( ABSTAN(IOUT,ITAN), IOUT = IOUT1, IOUT2 )
            ELSE
              WRITE ( LUN, IOSTAT=IOS, ERR=900 ) NOUT, 
     &          ( SNGL(ABSTAN(IOUT,ITAN)), IOUT = IOUT1, IOUT2 )
            END IF
          ELSE   
            IF ( DBLFLG ) THEN
              WRITE ( LUN, *, IOSTAT=IOS, ERR=900 )
     &        ( ABSTAN(IOUT,ITAN), IOUT = IOUT1, IOUT2 )
            ELSE
              WRITE ( LUN, *, IOSTAT=IOS, ERR=900 )
     &          ( SNGL(ABSTAN(IOUT,ITAN)), IOUT = IOUT1, IOUT2 )
            END IF
          END IF
        ELSE IF ( TYPE .EQ. 'RAD' ) THEN
          IF ( BINFLG ) THEN
            IF ( DBLFLG ) THEN
              WRITE ( LUN, IOSTAT=IOS, ERR=900 ) NOUT, 
     &          ( RADTAN(IOUT,ITAN)*RADFAC, IOUT = IOUT1, IOUT2 )
            ELSE
              WRITE ( LUN, IOSTAT=IOS, ERR=900 ) NOUT, 
     &          ( SNGL(RADTAN(IOUT,ITAN)*RADFAC), IOUT = IOUT1, IOUT2 )
            END IF
          ELSE   
            IF ( DBLFLG ) THEN
              WRITE ( LUN, *, IOSTAT=IOS, ERR=900 )
     &          ( RADTAN(IOUT,ITAN)*RADFAC, IOUT = IOUT1, IOUT2 )
            ELSE
              WRITE ( LUN, *, IOSTAT=IOS, ERR=900 )
     &          ( SNGL(RADTAN(IOUT,ITAN)*RADFAC), IOUT = IOUT1, IOUT2 )
            END IF
          END IF
        ELSE IF ( TYPE .EQ. 'RJT' ) THEN ! Convert radiance to Bright.Temp.
          IF ( BINFLG ) THEN
            IF ( DBLFLG ) THEN
              WRITE ( LUN, IOSTAT=IOS, ERR=900 ) NOUT, 
     &          ( C2 * RADTAN(IOUT,ITAN) / ( C1 * WNOFIN(IOUT)**2 ), 
     &            IOUT = IOUT1, IOUT2 )
            ELSE
              WRITE ( LUN, IOSTAT=IOS, ERR=900 ) NOUT, 
     &          ( SNGL ( C2 * RADTAN(IOUT,ITAN) / 
     &            ( C1 * WNOFIN(IOUT)**2) ), IOUT = IOUT1, IOUT2 )
            END IF
          ELSE   
            IF ( DBLFLG ) THEN
              WRITE ( LUN, *, IOSTAT=IOS, ERR=900 )
     &          C2 * RADTAN(IOUT,ITAN) / ( C1 * WNOFIN(IOUT)**2 )
            ELSE
              WRITE ( LUN, *, IOSTAT=IOS, ERR=900 )
     &          ( SNGL ( C2 * RADTAN(IOUT,ITAN) / 
     &                 ( C1 * WNOFIN(IOUT)**2 ) ), IOUT = IOUT1, IOUT2 )
            END IF
          END IF
        ELSE IF ( TYPE .EQ. 'TRA' ) THEN
          IF ( BINFLG ) THEN
            IF ( DBLFLG ) THEN
              WRITE ( LUN, IOSTAT=IOS, ERR=900 ) NOUT, 
     &           ( TRATAN(IOUT,ITAN), IOUT = IOUT1, IOUT2 )
            ELSE
              WRITE ( LUN, IOSTAT=IOS, ERR=900 ) NOUT, 
     &          ( SNGL(TRATAN(IOUT,ITAN)), IOUT = IOUT1, IOUT2 )
            END IF
          ELSE   
            IF ( DBLFLG ) THEN 
              WRITE ( LUN, *, IOSTAT=IOS, ERR=900 )    ! was: '(X,3F20.16)'
     &          ( TRATAN(IOUT,ITAN), IOUT = IOUT1, IOUT2 )
            ELSE
              WRITE ( LUN, *, IOSTAT=IOS, ERR=900 )
     &          ( SNGL(TRATAN(IOUT,ITAN)), IOUT = IOUT1, IOUT2 )
            END IF
          END IF
        ELSE
          STOP 'F-WRTOUT: Logical Error#2'
        END IF
C
        IF ( LUNFLG ) THEN
          CLOSE ( LUN, IOSTAT=IOS, ERR=900 )
        ELSE
          LUN = LUN + 1
        END IF
      END DO
C
      FAIL = .FALSE.
      RETURN
C
  900 FAIL = .TRUE.
      WRITE ( ERRMSG, '(A,I11)' ) 
     &  'F-WRTOUT: I/O failure on '//TYPE//' output file. IOSTAT=', IOS
      END
