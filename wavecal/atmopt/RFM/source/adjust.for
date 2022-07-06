      SUBROUTINE ADJUST ( IPTH )
C
C VERSION
C     20-JUL-05  AD  Set YMXADJ=0 to be safe.
C     05-MAY-05  AD  Add PTBVIB: Allow for vibrational temperature Jacobians
C     10-AUG-03  AD  Change VLIGHT from SNGL to DBLE 
C     02-AUG-99  AD  Add some DBLE ( ) and SNGL ( ) to make explicit
C     21-OCT-97  AD  Change ANLTE,CNLTE to DOUBLE PRECISION
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     28-SEP-96  AD  Add pressure shift to WNOADJ
C                    Remove references to QTINT (only called for line data)
C     20-SEP-96  AD  Change name of IUSVIB,ILSVIB to IUSNTE,ILSNTE
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Adjust line parameters for path conditions
C     Called by RFMFIN and RFMWID.
C     The line strength and width at 296K (TEMREF) and 1 atm (PREF) read from 
C     the HITRAN line data file are adjusted for atmospheric path conditions. 
C     The doppler halfwidth is also calculated. These parameters are then used 
C     in the lineshape formulation. 
C
C     Adjusts data in /HITCOM/ and loads into /ADJCOM/
C   
C     Line mixing data and temperature dependence interpolation due to 
C     L.L.Strow (private communication with Dave Edwards)
C
C     Line strengths for different isotopes of a gas on the HITRAN data base 
C     weighted by atmospheric abundance of the particular isotope. Absolute 
C     strengths may be obtained by dividing by this abundance i.e. 
C     STREN/ABUN(IDGAS,ISO) where IDGAS is the gas ID and ISO the isotope ID. 
C     This will be important when performing calculations for other planetary 
C     atmospheres.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER IPTH !  I  Index of current path
C
      EXTERNAL
     &  NTECLC ! Calculate various non-LTE parameters for line 
     &, PTBVIB ! Perturb/unperturb VT profiles for Jacobian calc.
     &, QTFCT  ! Calculate total internal partition sums.
     &, YMIX   ! Calculate line-mixing y-coefficient
C
      REAL  YMIX
C
C GLOBAL CONSTANTS
      INCLUDE 'idxcon.inc' ! HITRAN molecular indices
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
C
C COMMON VARIABLES
      INCLUDE 'adjcom.inc' ! Path-Adjusted line data 
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'hitcom.inc' ! HITRAN line database variables
      INCLUDE 'pthcom.inc' ! RFM path data
C
C LOCAL VARIABLES
      INTEGER          IGAS      ! Array index for gas
      REAL             BR        ! Broadening factor at Standard Temp.
      REAL             PRE       ! CG Pressure [mb]
      REAL             PPA       ! CG Partial pressure [mb]
      REAL             SQ        ! Ratio of tps@296K/tps@path_temp
      REAL             TEM       ! CG Temperature [K]
      DOUBLE PRECISION ANLTE     ! Non-LTE Correction factor for k absorption
      DOUBLE PRECISION CNLTE     ! Non-LTE Correction factor for c absorption
      DOUBLE PRECISION SB        ! exp( -hcE_l/kT_path ) / exp( -hcE_l/kT_ref )
      DOUBLE PRECISION SE        ! Ration of stimulated emission @path/@ref
      DOUBLE PRECISION GAMMA     ! exp ( -hcv/kT )
      DOUBLE PRECISION GAMREF    ! exp ( -hcv/kT_ref )
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      PRE = PREPTH(IPTH)
      PPA = PPAPTH(IPTH)
      TEM = TEMPTH(IPTH)
      IGAS = IGSPTH(IPTH)
      WNOADJ = WNUM + DBLE ( PREPTH(IPTH) * TSP ) ! Pressure shift (normally 0)
C
C Convert for line width in cm-1 at 296K and 1atm.
C  
      IF ( SBROAD .EQ. 0.0 ) THEN                 ! If self-broad.HW is zero...
        IF ( IDGAS .EQ. IDXH2O ) THEN                   ! ID=1 = Water vapour 
          SBROAD = 5.0 * ABROAD
          BR = ( ABROAD * ( PRE - PPA ) + SBROAD * PPA ) / PREREF
        ELSE                     ! SBROAD = ABROAD for want of something better
          BR = ABROAD * PRE / PREREF  
        ENDIF
      ELSE
        BR = ( ABROAD * ( PRE - PPA ) + SBROAD * PPA ) / PREREF
      ENDIF
      WIDADJ = BR * ( TEMREF / TEM )**ABCOEF
C
C Calculate Doppler half-width at half-max HWHM in cm-1. 
C
      DOPADJ = SNGL ( WNUM / VLIGHT ) * 
     &                    SQRT ( R2 * TEM / WGTGAS(ISO,IGAS) )
C
C Calculate the line mixing y coefficient (only CO2 lines at present)
C
      IF ( MIXFLG .AND. IDGAS .EQ. IDXCO2 ) THEN          ! ID= 2 = CO2
        YMXADJ = YMIX ( TEM, PRE, PPA )
      ELSE
        YMXADJ = 0.0
      ENDIF
C
C Convert for line strength in cm-1.(mol.cm-2)-1 at 296K.
C
C Boltzman factor for lower state energy
C
      SB = DEXP ( DBLE(ELS) * C2 * DBLE(TEM-TEMREF)/DBLE(TEM*TEMREF) )
C
C Stimulated emission 
C
      GAMMA = DEXP ( -C2 * WNUM / DBLE ( TEM ) )
      GAMREF = DEXP ( -C2 * WNUM / DBLE ( TEMREF ) )
      SE = ( 1.D0 - GAMMA ) / ( 1.D0 - GAMREF )
C
C Nonlte calculation of absorption coefficient modifiers
C
      SQ = 1.0                                       ! just to keep FLINT happy
      IF ( IUSNTE .NE. 0 .OR. ILSNTE .NE. 0 ) THEN
        IF ( IVJPTH(IPTH) .GT. 0 ) CALL PTBVIB ( IVJPTH(IPTH) ) 
        CALL NTECLC ( PRE, TEM, GAMMA, ANLTE, CNLTE, SQ )
        IF ( IVJPTH(IPTH) .GT. 0 ) CALL PTBVIB ( 0 ) 
      ELSE
        ANLTE = 1.0D0
        CNLTE = 1.0D0       
        CALL QTFCT ( IDGAS, ISO, TEM, SQ )
      ENDIF
C
      STRADJ = AMTPTH(IPTH) * STREN * SNGL ( SB * SE ) * SQ
      ANTADJ = ANLTE
      CNTADJ = CNLTE
C
      END            
