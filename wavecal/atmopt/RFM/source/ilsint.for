      SUBROUTINE ILSINT ( ISPC, ILIM, FNC, ILOW, IUPP, FAIL, ERRMSG )
C
C VERSION
C     24-APR-12  AD  Change arguments to DLOOKP
C     04-APR-06  AD  Change DFLOAT to DBLE
C     27-APR-00  AD  Change ILS Fn to DP
C     04-JAN-99  AD  Add AVG flag option 
C     05-AUG-99  AD  Use ILOW, IUPP =  NINT( ) and WNO = SNGL( ) 
C                    Allow for PTD being positive or negative
C     23-APR-98  AD  Documentation change: not called by ILSCHK.
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Interpolate ILS function to fine-grid spacing.
C     Called by RFMILS for each spectral range.
C     This requires NFIN, WNRFIN for current spectral range to be in /FINCOM/
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      ISPC            !  I  Index of Spectral range for calc.
      INTEGER      ILIM            !  I  Limits on FNC array
      DOUBLE PRECISION FNC(-ILIM:ILIM) !  O  Array for interpolated ILS
      INTEGER      ILOW            !  O  Index of lowest non-zero pt in FNC
      INTEGER      IUPP            !  O  Index of highest non-zero pt in FNC
      LOGICAL      FAIL            !  O  Set TRUE if a fatal error occurs
      CHARACTER*80 ERRMSG          !  O  Error message written if FAIL=TRUE
C
      EXTERNAL
     &  DLOOKP  !  General purpose interpolation routine
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'ilscom.inc' ! Instrument Lineshape functions
      INCLUDE 'outcom.inc' ! RFM output file data
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
C 
C LOCAL VARIABLES
      INTEGER  I    ! Counter for points in interpolated ILS
      INTEGER  IPT  ! Counter for tabulated ILS points
      INTEGER  ILS  ! Index of Tabulated Instrument Line Shape
      INTEGER  J    ! Guessed value for interpolation
      INTEGER  NPT  ! No. of points in Tabulated ILS
      DOUBLE PRECISION PT1  ! Starting wavenumber of Tabulated ILS [/cm]
      DOUBLE PRECISION PTD  ! Wavenumber increment of Tabulated ILS
      DOUBLE PRECISION SUM  ! For normalising interpolated ILS
      DOUBLE PRECISION WNO  ! Wavenumber of each interpolated point
      DOUBLE PRECISION WNOILS(MAXILP) ! Wavenumber array for tabulated ILS 
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      IF ( AVGFLG ) THEN
        IUPP = NINT ( WNROUT / WNRFIN )
        ILOW = - IUPP
        IF ( IUPP .GT. ILIM ) THEN
          FAIL = .TRUE.
          WRITE ( ERRMSG, '(A,I11)' ) 
     &      'F-ILSINT: MAXILF in rfmsiz.inc too small. '//
     &      'Required value .GE.', IUPP
          RETURN
        END IF
c        FNC(0) = 1.0D0 / DFLOAT(IUPP) 
        FNC(0) = 1.0D0 / DBLE(IUPP) 
        DO I = 1, IUPP
c          FNC(I) = ( 1.0D0 - DFLOAT(I)/DFLOAT(IUPP) ) / DFLOAT(IUPP)
          FNC(I) = ( 1.0D0 - DBLE(I)/DBLE(IUPP) ) / DBLE(IUPP)
          FNC(-I) = FNC(I)
        END DO
        RETURN
      END IF
C
      ILS = ILSSPC(ISPC)
      PT1 = PT1ILS(ILS)
      PTD = PTDILS(ILS)
      NPT = NPTILS(ILS)
C
      ILOW = NINT ( PT1 / WNRFIN )          ! PT1 is usually negative
      IUPP = NINT ( ( PT1 + ( NPT - 1 ) * PTD ) / WNRFIN )
C
C Shouldn't fail here since MAXILF limit already checked by ILSCHK 
      IF ( MIN ( ILOW, IUPP ) .LT. -ILIM .OR. 
     &     MAX ( ILOW, IUPP ) .GT. +ILIM      ) THEN
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,2I11)' ) 
     &    'F-ILSINT: interpolated ILS too large. Reqd MIN,MAX=', 
     &    ILOW, IUPP
        RETURN
      END IF
C
C Construct wavenumbers array for tabulated ILS
C
      DO IPT = 1, NPT 
        WNOILS(IPT) = PT1 + ( IPT - 1 ) * PTD 
      END DO
C
      J = 1
      SUM = 0.0D0
      DO I = ILOW, IUPP
        WNO = WNRFIN * I  
        CALL DLOOKP ( WNO, FNC(I), WNOILS, FNCILS(1,ILS), NPT, J )
        SUM = SUM + FNC(I) 
      END DO
C
      DO I = ILOW, IUPP
        FNC(I) = FNC(I) / SUM
      END DO
C
      FAIL = .FALSE.
      END
