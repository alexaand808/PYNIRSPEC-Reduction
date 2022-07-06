      SUBROUTINE QTFCT ( IDGAS, ISO, TEM, SQ )
C
C VERSION
C     01OCT13 AD Add new molecules and change warnings using QTWARN
C     08AUG13 AD New isotopes and molecules. Assume QTEM values now 
C                already scaled by Q(296K).
C                Limit the number of warning messages issues
C                If V42 flag, call old version QTFCT_V42 instead
C     03APR13 AD Extra isotopes: print warning and use main isotope TIPS
C     23JUL09 AD Correction to calculation of max ISPE 
C     14APR08 AD Rewritten using 4pt Lagrangian interpolation
C                Add include file tpsdat.inc.
C                Externally the same as previous version
C     23FEB03 AD Add BrO as ID=40, move GeH4 (was ID=40) to new ID=46
C     13FEB03 AD Updated to use TIPS97.FOR coefficients plus 
C                extra values supplied by Bianca Dinelli (BMD).
C                Correction to handling of minor isotopes.
C     23JAN01 AD Local version, new HNO3 coefficients
C     16OCT00 AD Check for IDGAS >38 as well
C     13SEP99 AD Add extra isotope for C2H4 (ID=38)
C     09AUG99 AD Remove surplus brackets around QCOEFH(ISPC,4)
C     08JUL98 AD Add new HNO3 coefficients (commented out at the moment)
C     03NOV97 AD Remove duplicate NO high-temp data (QCOEFH)
C     25AUG97 AD Add extra data QCOEFH to allow for T=500-1500K 
C     03MAR97 AD Version 3.
C     12DEC96 AD Also check for C2H6 (no TIPS coefficients).
C     23NOV96 AD Revised for HITRAN'96 (except for HNO3 - keep HITRAN'92).
C     01OCT96 AD Version 2.
C     01SEP96 AD Version 1.
C
C DESCRIPTION
C     Calculate total internal partition sums.
C     Called by ADJUST, NTECLC and MX2CLC.
C     The ratio of total internal partition sums, Q(296K)/Q(Path Temp) is 
C     calculated for a given molecule and isotopic species. The line strength 
C     can then be adjusted by this factor.
C
C     Based on program TIPS_2003.FOR by R.R.GAMACHE
C     The original version used 4th order polynomials to span entire 
C     temperature ranges, however the 2003 version uses local 4th order 
C     Lagrangian interpolation which should be more accurate and spans a 
C     greater temperature range overall
C     
C     In general, Lagrangian interpolation of y(x) is given by
C         y(x) = sum_ijkl  yi(x-xj)(x-xk)(x-xl)/[(xi-xj)(xi-xk)(xi-xl)]
C     where ijkl are the four different points used for interpolation.
C
C     In this code (as in the Gamache version) these points are chosen such 
C     that x lies between j and k, except at the extreme ends of the tabulated
C     temperature where the points used are the edge points of the tabulation.
C     In this code it it assumed that the xi are equally spaced (25K intervals)
C     (although this is not checked) which allows for significant simplification
C     compared to the more general Gamache subroutine "AtoB"
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER IDGAS !  I  HITRAN gas ID
      INTEGER ISO   !  I  HITRAN isotope ID 
      REAL    TEM   !  I  Path temperature [K]
      REAL    SQ    !  O  Ratio of tot. part. sum at 296K to tps at pth temp
C
      EXTERNAL
     &  QTFCT_V42 ! Old version of QTFCT used in RFM v4.2
     &, QTWARN    ! Handle warning messages from subroutine QTFCT
     &, TPSLKP    ! Look up tabulated TIPS factor
      REAL TPSLKP
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM array sizes
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
C
C LOCAL CONSTANTS
      INTEGER MAXSPE              ! Total number of isotope species (all gases)
        PARAMETER ( MAXSPE = 133 )
      INTEGER MAXIDG              ! Highest HITRAN index for TIPS data
        PARAMETER ( MAXIDG = 48 ) ! = 1 + No of defined molecules in tpsdat.inc
      INTEGER NT                  ! No. of temperatures for TIPS tabulations
        PARAMETER ( NT = 119 )    
C The following parameters are actually obtainable from the TEMVAL data but
C set as parameters to speed up code
      REAL    DTEM                ! Temperature increment in tabulated data
        PARAMETER ( DTEM = 25.0 )
      REAL    TEM1                ! Lowest temperature in tabulated data
        PARAMETER ( TEM1 = 60.0 )
      REAL    TEM2                ! Highest temperature in tabulated data
        PARAMETER ( TEM2 = 3010.0 )   ! NT = (TEM2-TEM1)/DTEM + 1 
      REAL    DTI3                ! 1/(DTEM)**3
        PARAMETER ( DTI3 = 6.4E-5 )  ! 6.4E-5 = (1/25)^3
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'tpscom.inc' ! Tabulated TIPS data
C
C LOCAL VARIABLES
      INTEGER I               ! Counter for DATA statements
      INTEGER IT,JT,KT,LT     ! Indices of 4 consecutive tabulation points
      INTEGER ISPE            ! Counter for isotope species (all gases)
      INTEGER ISPOFF(MAXIDG)  ! Offset for each gas in the QTEM arrays
      REAL    DTI,DTJ,DTK,DTL ! Position of TEM relative to tabulation points
      REAL    QTEM(NT,MAXSPE) ! Tabulated TIPS values for each ISPE, T
      REAL    TEMVAL(NT)      ! Tabulation temperatures [K]
      REAL    QT              ! Interpolated TIPS value at T
C
C DATA STATEMENTS
      include 'tpsdat.inc'    ! Data for TIPS Lagrangian interpolation
      SAVE ISPOFF, TEMVAL, QTEM
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      IF ( V42FLG .AND. IDGAS .LE. 38 ) THEN ! Old QTFCT has same IDs 1:38 only
        CALL QTFCT_V42 ( IDGAS, ISO, TEM, SQ )
        RETURN
      END IF
C
      IF ( IDGAS .GE. MAXIDG .OR. ! Assume extra molecule
     &     IDGAS .EQ. 39 .OR.     ! no data yet for CH3OH
     &     IDGAS .EQ. 43 ) THEN   ! no data yet for SO3
        SQ = ( TEMREF / TEM )**1.5
        CALL QTWARN ( IDGAS, -1, TEM )   ! -1 flags this as no TIPS for molec
        RETURN
      ELSE IF ( TPSFLG .AND. IDXTPS(ISO,IDGAS) .GT. 0 ) THEN
        SQ = TPSLKP ( IDGAS, ISO, TEM )
      ELSE
        ISPE = ISPOFF(IDGAS) + ISO     ! Position in array of molecular isotope
C Check for isotopes not included in TIPS data
C Note that ISPOFF array extends to max IDGAS plus 1 to allow for this check
        IF ( ISPE .GT. ISPOFF(IDGAS+1) ) THEN
          CALL QTWARN ( IDGAS, ISO, TEM ) 
C Just use main isotope value instead 
          ISPE = ISPOFF(IDGAS) + 1
        END IF 
C Find indices of 4 consecutive points for interpolation, IT,JT,KT,LT, aiming
C to bracket TEM between JT LE TEM LT KT, but ensuring that the 4 chosen points
C remain within the tabulation if TEM is near the edges
        JT = INT ( ( TEM - TEM1 ) / DTEM ) + 1
        IF ( JT .LT. 2 ) THEN 
          IF ( TEM .LT. TEM1 ) CALL QTWARN ( IDGAS, ISO, TEM )
          JT = 2
        ELSE IF ( JT .GT. NT-2 ) THEN
          IF ( TEM .GT. TEM2 ) CALL QTWARN ( IDGAS, ISO, TEM )
          JT = NT - 2
        END IF
        IT = JT - 1
        KT = JT + 1
        LT = JT + 2
C
        DTI = TEM - TEMVAL(IT)
        DTJ = TEM - TEMVAL(JT)
        DTK = TEM - TEMVAL(KT)
        DTL = TEM - TEMVAL(LT)
C
        QT = DTI3 * ( -1.0/6.0*DTJ*DTK*DTL * QTEM(IT,ISPE)
     &                +1.0/2.0*DTI*DTK*DTL * QTEM(JT,ISPE)
     &                -1.0/2.0*DTI*DTJ*DTL * QTEM(KT,ISPE)
     &                +1.0/6.0*DTI*DTJ*DTK * QTEM(LT,ISPE) )
C
        SQ = 1.0D0 / QT
      END IF
C
      END
