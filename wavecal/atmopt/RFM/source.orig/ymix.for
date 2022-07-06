      REAL FUNCTION YMIX ( TEM, PRE, PPA )
C
C VERSION
C     11-JAN-05  AD  Only check part of BSLQ string since some versions of
C                    HITRAN have additional information
C     11-AUG-99  AD  Add SNGL ( ) in calc. of YMIX. Remove surplus brackets
C     17-AUG-98  AD  Remove general SAVE statement
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-OCT-96  AD  Change to use GENLN2 v4 coefficients.
C     19-SEP-96  AD  Comment changes only
C     01-SEP-96  AD  Version 1.
C 
C DESCRIPTION
C     Calculate line-mixing y-coefficient
C     Called by ADJUST if MIXFLG enabled.
C                      CO2 Q-branch mixing coefficient data for:
C                  1   618 cm-1  3  2, 626 (Q50 - Q 2,  615 - 619 cm-1)
C                  2   648 cm-1  2  1, 636 (Q 2 - Q50,  648 - 652 cm-1)
C (currently dummy)3   662 cm-1  2  1, 628 (Q 1 - Q50,  662 - 665 cm-1)
C                  4   667 cm-1  2  1, 626 (Q 2 - Q50,  667 - 671 cm-1)
C                  5   668 cm-1  4  2, 626 (Q 2 - Q81,  667 - 675 cm-1)
C                  6   720 cm-1  5  2, 626 (Q50 - Q 2,  718 - 721 cm-1)
C                  7   741 cm-1  8  4, 626 (Q81 - Q 2,  733 - 742 cm-1)
C                  8   791 cm-1  8  3, 626 (Q 2 - Q50,  791 - 794 cm-1)
C                  9  1932 cm-1  6  1, 626 (Q 2 - Q50, 1932 -1937 cm-1)
C                 10  2076 cm-1  8  1, 626 (Q 2 - Q50, 2076 -2080 cm-1)
C                 11  2093 cm-1 14  2, 626 (Q 2 - Q70, 2093 -2095 cm-1)
C                 12  2129 cm-1 15  2, 626 (Q50 - Q 2, 2128 -2130 cm-1)
C     Refer to Tobin & Strow for details of quantum numbers.
C     To be used with CO2MIX lines with LABEL '101293'
C
C     1st order line mixing = Y(T)*(200/T)^0.75, with
C     Y(T) = COMIX(1) + COMIX(2)*(T - 200) +
C     COMIX(3)*(T - 200)**2 + COMIX(4)*(T - 200)**3
C
C     Data relating to the ranges of sets is incorporated in SPCWID.FOR
C
C REFERENCE
C     TOBIN, D. and L.L. STROW
C     A Compilation of First-order Line-mixing Coefficients for CO2 Q-branches
C     J.Q.S.R.T., 52, 281 (1994).
C
C
      IMPLICIT NONE
C
C ARGUMENTS
      REAL  TEM   !  I  CG temperature [K]
      REAL  PRE   !  I  CG Pressure [atm]
      REAL  PPA   !  I  CG Partial Pressure [atm]
C
C LOCAL CONSTANTS
      INTEGER MAXSET                 ! Max.no of sets 
        PARAMETER ( MAXSET = 12 )
      INTEGER MAXMLN                 ! Max.no.lines in any set
        PARAMETER ( MAXMLN = 80 )
      INTEGER MAXMCO                 ! No. temperature polynomial coefficients
        PARAMETER ( MAXMCO = 4 ) 
      REAL    TEMREF                 ! Reference temperature for calculations
        PARAMETER ( TEMREF = 200.0 )
C
C COMMON VARIABLES
      INCLUDE 'hitcom.inc' ! HITRAN line database variables
C  
C LOCAL VARIABLES
      INTEGER  IMCO    ! Counter for temperature polynomial coefficients [IC]
      INTEGER  IMLN    ! Counter for lines in each set [IQ]
      INTEGER  ISET    ! Counter for number of sets [IS]
      INTEGER  IMMOL(MAXSET)  ! GAS IDs of each set       
      INTEGER  IMISO(MAXSET)  ! Isotope IDs of each set       
      INTEGER  MGQL(MAXSET)   ! AFGL global quanta IDs of lower states
      INTEGER  MGQU(MAXSET)   ! AFGL global quanta IDs of upper states
      INTEGER  NMLN(MAXSET)   ! No lines in each set [NSNUM]
      CHARACTER*9 MLQL(MAXMLN,MAXSET) ! AFGL local quanta ID of lower states
      DOUBLE PRECISION TEMDIF  ! Difference from calc.ref. temperature
      DOUBLE PRECISION YMIXS  ! Self-broadened mixing
      DOUBLE PRECISION YMIXF  ! Foreign-broadened mixing
      DOUBLE PRECISION COMIXS(MAXMCO,MAXMLN,MAXSET) ! Self-broad. Mix. Coeffs.
      DOUBLE PRECISION COMIXF(MAXMCO,MAXMLN,MAXSET) ! For.-broad. Mix.Coeffs
         DATA IMMOL / MAXSET * 2 /        ! 2 is CO2
         DATA IMISO / 1, 2, 3, 9 * 1 /    ! 1 is 626, 2 is 636, 3 is 628
         DATA MGQL  / 2, 1, 1, 1, 2, 2, 4, 3, 1, 1, 2, 2 /
         DATA MGQU  / 3, 2, 2, 2, 4, 5, 8, 8, 6, 8, 14, 15 /
         DATA NMLN  / 25, 25, 50, 25, 80, 25, 80, 25, 25, 25, 69, 25 /
         SAVE IMMOL, IMISO, MGQL, MGQU, NMLN
      INCLUDE 'mixdat.inc'  ! Data for MLQL,COMIXS,COMIXF
         SAVE MLQL, COMIXS, COMIXF
C
C EXECUTABLE CODE -------------------------------------------------------------
C      
      DO ISET = 1, MAXSET
        IF ( IMMOL(ISET) .EQ. IDGAS .AND. IMISO(ISET) .EQ. ISO .AND.
     &       IUSGQ .EQ. MGQU(ISET) .AND. ILSGQ .EQ. MGQL(ISET) ) THEN
          DO IMLN = 1, NMLN(ISET)
            IF ( BSLQ(5:9) .EQ. MLQL(IMLN,ISET)(5:9) ) THEN
              TEMDIF = DBLE ( TEM - TEMREF )
              YMIXS =   COMIXS(1,IMLN,ISET) + TEMDIF *
     &                ( COMIXS(2,IMLN,ISET) + TEMDIF * 
     &                ( COMIXS(3,IMLN,ISET) + TEMDIF * 
     &                  COMIXS(4,IMLN,ISET)            ))
              YMIXF =   COMIXF(1,IMLN,ISET) + TEMDIF *
     &                ( COMIXF(2,IMLN,ISET) + TEMDIF *
     &                ( COMIXF(3,IMLN,ISET) + TEMDIF * 
     &                  COMIXF(4,IMLN,ISET)            ))
              YMIX = ( TEMREF / TEM ) ** 0.75 *  
     &          SNGL ( DBLE (PPA) * YMIXS + DBLE (PRE-PPA) * YMIXF )
              RETURN
            END IF
          END DO
        END IF 
      END DO
      YMIX = 0.0
      END
