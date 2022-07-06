C $Id: rfrgl2.for,v 1.9 1997/03/03 17:17:53 adarby Exp $
      REAL FUNCTION REFRAC ( IATM )
C
C VERSION
C     23-MAR-00  AD  Revised to use Kaye & Laby form
C     09-AUG-99  AD  Use local variable HGT in LOOKUP 
C     24-JAN-99  AD  Original. Based on RFRGL2.
C
C DESCRIPTION
C     Calculate refractivity.
C     Called for each profile level by ATMAUX if GL2 option selected.
C     Basically the GENLN2 algorithm but with corrected value of reference
C     temperature (GENLN2 uses TEMREF=296.0, correct value TREF=288.16)
C     Uses any loaded water vapour profile, otherwise US Std Atmos. profile.
C     Ref: Tables of Physical and Chemical Constants (16th Ed), 
C          G. W. C. Kaye and T. H. Laby
C          Published by Longman, 1995.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER  IATM      !  I  Index of atmospheric level
C
      EXTERNAL
     &  LOOKUP ! General purpose interpolation routine
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
      INCLUDE 'idxcon.inc' ! HITRAN molecular indices
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
C
C LOCAL CONSTANTS
      INTEGER       MAXUSA        ! No. elements in US Std Atm. profiles
        PARAMETER ( MAXUSA = 50 )
      REAL          TREF          ! Ref.Temperature for Edlen refrac. calc.
        PARAMETER ( TREF = 288.16 )
C
C LOCAL VARIABLES
      INTEGER  IGSH2O ! Index of H2O in *GAS if stored as a profile, else=0
      INTEGER  JGUESS ! Guess value for interpolation point
      LOGICAL  FIRST  ! T=first time this function gets called
      REAL     HGT    ! Height required for Look-up
      REAL     WPP    ! Water Vapour Partial Pressure [mb]
      REAL     HGTUSA(MAXUSA) ! US Std Atm. altitudes [km] 
      REAL     H2OUSA(MAXUSA) ! US Std Atm. Water Vapour [ppmv]
      REAL     SIGMA2 ! (Wavenumber[/micron])^2
        DATA IGSH2O / 0 /        
        DATA JGUESS / 1 / 
        DATA FIRST / .TRUE. /
        SAVE IGSH2O, JGUESS, FIRST
        DATA HGTUSA / 
     &    0.0,  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0, 
     &   10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0,
     &   20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 27.5, 30.0, 32.5, 35.0,
     &   37.5, 40.0, 42.5, 45.0, 47.5, 50.0, 55.0, 60.0, 65.0, 70.0,
     &   75.0, 80.0, 85.0, 90.0, 95.0,100.0,105.0,110.0,115.0,120.0 /
        SAVE HGTUSA
        DATA H2OUSA / 
     &   7.745E+03, 6.071E+03, 4.631E+03, 3.182E+03, 2.158E+03,
     &   1.397E+03, 9.254E+02, 5.720E+02, 3.667E+02, 1.583E+02,
     &   6.996E+01, 3.613E+01, 1.906E+01, 1.085E+01, 5.927E+00,
     &   5.000E+00, 3.950E+00, 3.850E+00, 3.825E+00, 3.850E+00,
     &   3.900E+00, 3.975E+00, 4.065E+00, 4.200E+00, 4.300E+00,
     &   4.425E+00, 4.575E+00, 4.725E+00, 4.825E+00, 4.900E+00,
     &   4.950E+00, 5.025E+00, 5.150E+00, 5.225E+00, 5.250E+00,
     &   5.225E+00, 5.100E+00, 4.750E+00, 4.200E+00, 3.500E+00,
     &   2.825E+00, 2.050E+00, 1.330E+00, 8.500E-01, 5.400E-01,
     &   4.000E-01, 3.400E-01, 2.800E-01, 2.400E-01, 2.000E-01 /
        SAVE H2OUSA
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      IF ( FIRST ) THEN 
        IGSH2O = IGSMOL(IDXH2O)
        FIRST = .FALSE.
      END IF
C
      IF ( IGSH2O .EQ. 0 ) THEN
        HGT = HGTATM(IATM)
        CALL LOOKUP ( HGT, WPP, HGTUSA, H2OUSA, MAXUSA,
     &                JGUESS, 1 )
        WPP = WPP * PREATM(IATM) * 1.0E-6
      ELSE
        WPP = VMRATM(IATM,IGSH2O) * PREATM(IATM)
      END IF
C
C K&L formula expressed in sigma^2, where sigma=1/wavelength[microns] 
      SIGMA2 = (SNGL(WMDSPC)*0.0001)**2
      REFRAC = 1.0E-8 * PREATM(IATM) / ATMB * TREF / TEMATM(IATM) *
     &         (    8342.54 +
     &           2406147.0 / ( 130.0 - SIGMA2 ) +
     &             15998.0 / (  38.9 - SIGMA2 )    ) 
c K&L Water Vapour correction only appropriate for 405-644nm, so comment out
c     &       - 1.0E-8 * WPP * ( 3.7345 - 0.0401 * SIGMA2 ) 
C
      END
C-------------------------------------------------------------------------------
C                                NOTICE
C
C     This software module is part of the MIPAS Reference Forward Model supplied
C     to the European Space Agency under ESTEC contract no. 11886/96/NL/GS.
C        
C     All rights, title and interest in and to the copyright
C     and any other proprietary right, is held by
C     The University Corporation for Atmospheric Research, having a
C     correspondence address of P.O. Box 3000, Boulder, Colorado 80307, USA.
C
C     However, note that all inquiries concerning the MIPAS
C     Reference Forward Model should be submitted to the The Manager (Projects),
C     AOPP,  Clarendon Laboratory, Parks Road, Oxford, OX1 3PU, UK.
C     (Tel: +44-8165-272900,    Fax: +44-1865-272924).  
C
C-------------------------------------------------------------------------------
