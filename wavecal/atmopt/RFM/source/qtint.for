C $Id: qtint.for,v 1.8 1997/03/03 17:17:39 adarby Exp $
      SUBROUTINE QTINT ( IDGAS, TEM, SQ )
C
C VERSION
C      17-AUG-98  AD  Remove general SAVE statement.
C      03-MAR-97  AD  Version 3.
C      01-OCT-96  AD  Version 2.
C      01-SEP-96  AD  Version 1.
C 
C DESCRIPTION 
C     Calculate total internal partition sums for tabulated data. 
C     Called by ADJUST, NTECLC and RFMXSC.
C     Primarily for use with the heavy molecule cross-section data when there 
C     are no temperature dependent measurements. The total internal partition 
C     sums for the first 32 HITRAN molecules are calculated by routine QTFCT.
C     Data from Aaron Goldman at NCAR taken from GENLN2 BLOCK DATA QSHEV
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER IDGAS !  I  HITRAN gas ID
      REAL    TEM   !  I  Path temperature [K]
      REAL    SQ    !  O  Ratio of tot. partition sum at 296K to tps at pth temp
C
      EXTERNAL
     &  LOOKUP ! General purpose interpolation routine
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
C
C LOCAL CONSTANTS 
      INTEGER       MAXITP         ! No of temp.pts for qvib interpolation
        PARAMETER ( MAXITP = 18 )
C
C LOCAL VARIABLES
      INTEGER J    ! Index for accessing lookup interpolation
      REAL    SPR
      REAL    SPT
      REAL    SPV
      REAL    TEMITP(MAXITP)  ! Vibrational partition function data
      REAL    VP51(MAXITP)    ! Vib.Part. Fn for CFCl3 (F11) 
      REAL    VP56(MAXITP)    ! ID=56: CHClF2 (F22)
      REAL    VP60(MAXITP)    ! ID=60: CCl4
      REAL    VP61(MAXITP)    ! ID=61: ClONO2
      REAL    VP62(MAXITP)    ! ID=62: N2O5
C
        DATA  J / 1 /
        SAVE  J
C
        DATA  TEMITP /
     & 170.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0,
     & 260.0, 270.0, 280.0, 290.0, 296.0, 300.0, 310.0, 320.0, 330.0 /
       SAVE TEMITP
C
        DATA VP51 /
     & 1.51031,  1.61072,  1.72276,  1.84728,  1.98524,  2.13770,
     & 2.30584,  2.49093,  2.69438,  2.91772,  3.16262,  3.43089,
     & 3.72449,  3.91368,  4.04554,  4.39633,  4.77933,  5.19722/
        DATA VP56 /
     & 1.08491,  1.10448,  1.12629,  1.15033,  1.17665,  1.20528,
     & 1.23625,  1.26961,  1.30541,  1.34372,  1.38461,  1.42816,
     & 1.47446,  1.50359,  1.52359,  1.57567,  1.63079,  1.68909/
        DATA VP60 / 
     & 1.82149,  1.98259,  2.16373,  2.36682,  2.59396,  2.84752,
     & 3.13009,  3.44455,  3.79404,  4.18202,  4.61229,  5.08897,
     & 5.61659,  5.95964,  6.20006,  6.84474,  7.55644,  8.34148/
        DATA VP61 / 
     & 1.80702,  1.91774,  2.03739,  2.16660,  2.30604,  2.45642,
     & 2.61853,  2.79320,  2.98129,  3.18376,  3.40161,  3.63590,
     & 3.88777,  4.04783,  4.15843,  4.44916,  4.76131,  5.09634/
        DATA VP62 /
     & 1.13198,  1.16234,  1.19659,  1.23494,  1.27764,  1.32494,
     & 1.37716,  1.43460,  1.49765,  1.56668,  1.64214,  1.72450,
     & 1.81429,  1.87196,  1.91207,  2.01846,  2.13414,  2.25985/
       SAVE VP51, VP56, VP60, VP61, VP62
C
C EXECUTABLE CODE --------------------------------------------------------------
C
C Rotational partition function - assume non-linear molecules
C
      SPR = ( TEMREF / TEM )**1.5
C
C Vibrational partition function by linear interpolation/extrapolation
C
      IF ( IDGAS .EQ. 51 ) THEN 
        CALL LOOKUP ( TEM, SPT, TEMITP, VP51, MAXITP, J, 1 )
        SPV = VP51(14) / SPT
      ELSE IF ( IDGAS .EQ. 56 ) THEN
        CALL LOOKUP ( TEM, SPT, TEMITP, VP56, MAXITP, J, 1 )
        SPV = VP56(14) / SPT
      ELSE IF ( IDGAS .EQ. 60 ) THEN
        CALL LOOKUP ( TEM, SPT, TEMITP, VP60, MAXITP, J, 1 )
        SPV = VP60(14) / SPT
      ELSE IF ( IDGAS .EQ. 61 ) THEN
        CALL LOOKUP ( TEM, SPT, TEMITP, VP61, MAXITP, J, 1 )
        SPV = VP61(14) / SPT
      ELSE IF ( IDGAS .EQ. 62 ) THEN
        CALL LOOKUP ( TEM, SPT, TEMITP, VP62, MAXITP, J, 1 )
        SPV = VP62(14) / SPT
      ELSE 
        SPV = 1.0
      END IF
C
      SQ = SPR*SPV
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
