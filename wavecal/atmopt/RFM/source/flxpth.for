      SUBROUTINE FLXPTH ( FAIL, ERRMSG )
C
C VERSION
C     07-AUG-13  AD  Remove redundant local variable MESSGE
C     16-JAN-04  AD  Add paths for storing absorption coeff if MTX+TRA
C     30-DEC-03  AD  Replace OBSATM with ATMLEV
C     04-NOV-01  AD  Limit paths to lowest output level if ZEN flag selected
C     27-APR-00  AD  Original.
C
C DESCRIPTION
C     Set up paths for flux calculations.
C     Called once by INPCHK if FLX flag is enabled.
C     Each level specified in the *LEV (=*TAN) section of the driver table
C     is associated with 1:NTAN. 
C     Rays at various angles for calculation are set up as NTAN+1:MTAN
C
      IMPLICIT NONE 
C
C ARGUMENTS
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  ATMLEV ! Find/insert atmospheric level for given altitude.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmcon.inc' ! Constants for RFM calculations
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES      
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'jaccom.inc' ! Jacobian data
      INCLUDE 'pthcom.inc' ! RFM path data
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL VARIABLES
      INTEGER      IATM    ! Counter for atmospheric profile levels
      INTEGER      IATM1   ! Index of layer required
      INTEGER      IGAS    ! Counter for absorbers
      INTEGER      IJAC    ! Counter for Jacobian elements (MTX option)
      INTEGER      IPTH    ! Counter for path segments
      INTEGER      ITAN    ! Counter for output levels
      INTEGER      JTAN    ! Secondary counter for output levels
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Assign atmospheric index to IATTAN for each output level 1:NTAN, inserting 
C extra profile levels if required. Since USRTAN is arranged in increasing 
C altitude, this should not affect the IATTAN index of previous levels.
      DO ITAN = 1, NTAN
        CALL ATMLEV ( USRTAN(ITAN), IATM, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        IATTAN(ITAN) = IATM
        CLCTAN(ITAN) = .FALSE.
      END DO
C
C If ZENith viewing, only require layers from top of the atmosphere down to
C the lowest output level, but include an extra level if cooling rates are
C required since this uses a quadratic fit over profile levels
      IF ( ZENFLG ) THEN                   ! matches logic in FLXPTH
        IATM1 = IATTAN(1)
        IF ( COOFLG .AND. IATM1 .GT. 1 ) IATM1 = IATM1 - 1
      ELSE
        IATM1 = 1
      END IF
C
C Determine number of calculated path segments required
      NCLC = NGAS * ( NATM - IATM1 )
C Extra paths for storing abs.coeff at output levels
      IF ( MTXFLG .AND. TRAFLG ) NCLC = NCLC + ( NTAN * NGAS ) 
      IF ( NCLC .GT. MAXCLC ) THEN
        WRITE ( ERRMSG, '(A,I11,A,I7)' )
     &    'F-FLXPTH: No.paths for line-by-line calc=', NCLC,
     &    ' >MAXCLC in RFMSIZ=', MAXCLC
        FAIL = .TRUE.
        RETURN
      END IF
C
C Determine total number path segments required
      NPTH = NCLC 
      IF ( NPTH .GT. MAXPTH ) THEN
        WRITE ( ERRMSG, '(A,I11,A,I7)' )
     &    'F-FLXPTH: Total no.paths required calc=', NPTH,
     &    ' >MAXPTH in RFMSIZ=', MAXPTH
        FAIL = .TRUE.
        RETURN
      END IF
C
C Assign path segments nominally to ITAN=1 although TAN arrays not actually
C used for calculations
      IPTH = 0
      DO IGAS = 1, NGAS
        DO IATM = IATM1, NATM-1
          IPTH = IPTH + 1
          IGSPTH(IPTH) = IGAS
          ITNPTH(IPTH) = 1
          IATPTH(IPTH) = IATM
          CLCPTH(IPTH) = .TRUE.
          ISCPTH(IPTH) = IPTH
          IDRPTH(IPTH) = 0
          NLNPTH(IPTH) = 0
        END DO
      END DO
C
C For combination of MTX+TRA flags set up extra paths for calculating C
C absorption coefficients at each output level
      IF ( MTXFLG .AND. TRAFLG ) THEN
        DO IGAS = 1, NGAS
          DO ITAN = 1, NTAN
            IPTH = IPTH + 1
            IGSPTH(IPTH) = IGAS
            ITNPTH(IPTH) = 1
            IATPTH(IPTH) = IATTAN(ITAN)
            CLCPTH(IPTH) = .TRUE.
            ISCPTH(IPTH) = IPTH
            IDRPTH(IPTH) = 2       ! Special flag to indicate atmos.level calc.
            NLNPTH(IPTH) = 0
          END DO
        END DO
      END IF
C
C If MTX option assign extra TAN arrays to hold outputs for each combination 
C of levels. This uses the Jacobian arrays since JAC and MTX incompatible.
      IF ( MTXFLG ) THEN
        IF ( NTAN .EQ. 1 ) THEN
          ERRMSG = 
     &      'F-FLXPTH: MTX option requires at least 2 output levels'
          FAIL = .TRUE.
          RETURN
        END IF
        NJAC = NTAN
        IF ( NJAC .GT. MAXJAC ) THEN
          WRITE ( ERRMSG, '(A,I10,A,I10)' )
     &      'F-FLXPTH: MTX Clc requires', NJAC,
     &      ' Jac.elements>MAXJAC (RFMSIZ.INC)=', MAXJAC
          FAIL = .TRUE.
          RETURN
        END IF
        LTAN = NTAN * ( 1 + NTAN ) 
        IF ( LTAN .GT. MAXTAN ) THEN
          WRITE ( ERRMSG, '(A,I10,A,I8)' )
     &      'F-FLXPTH: MTX Calcs require', LTAN,
     &      ' tan.pths > MAXTAN (RFMSIZ.INC)=', MAXTAN
          FAIL = .TRUE.
          RETURN
        END IF
C
        IJAC = 0
        LTAN = NTAN                    ! should end up as LTAN=NTAN+NJAC
        DO ITAN = 1, NTAN              ! output levels
          DO JTAN = 1, NTAN            ! secondary levels/ptb layers
            IJAC = JTAN
            LTAN = LTAN + 1
            IATJAC(IJAC) = IATTAN(JTAN)
            ITNJAC(ITAN,IJAC) = LTAN
            ITNJAC(LTAN,IJAC) = ITAN
            IGSJAC(IJAC) = 0           ! Value 0 used by JACNAM
          END DO
        END DO
C
      END IF
C
      FAIL = .FALSE.
C
      END
