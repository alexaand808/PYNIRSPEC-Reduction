C $Id: rfmfov.for,v 1.8 1997/03/03 17:17:43 adarby Exp $
      SUBROUTINE RFMFOV 
C
C VERSION
C     26-MAY-00  AD  Convert from S.P. to D.P.
C     09-AUG-99  AD  Use explicit DBLE( ) in transmittance calc.
C     03-JAN-99  AD  Allow for Jacobians
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Apply field-of-view convolution.
C     Called by RFM for each widemesh interval.
C     With Jacobian calculations the FOV convolved Jacobians end up in
C     NTAN+1:2*NTAN (IJAC=1);  2*NTAN+1:3*NTAN (IJAC=2) ... etc.
C     Routine RFMPTB checks that there will be enough space for this.
C
      IMPLICIT NONE
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'fovcom.inc' ! FOV Data
      INCLUDE 'jaccom.inc' ! Jacobian data
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL VARIABLES
      INTEGER          IFIN           ! Fine mesh counter
      INTEGER          IFOV           ! Field of view point counter
      INTEGER          IJAC           ! Jacobian element counter
      INTEGER          ITAN           ! Tangent height counter
      INTEGER          JTAN           ! Index of Tan.Hgt contributing to FOV
      INTEGER          NXTTAN         ! Index of next o/p tan.path 
      DOUBLE PRECISION RADTMP(MAXTAN) ! Temporary Radiance array
      DOUBLE PRECISION TRATMP(MAXTAN) ! Temporary transmission array
C
C EXECUTABLE CODE ------------------------------------------------------------
C
      DO IFIN = 1, NFIN
        DO ITAN = 1, NTAN
          RADTMP(ITAN) = 0.0D0
          TRATMP(ITAN) = 0.0D0
          DO IFOV = 1, NFOV
            JTAN = IFVTAN(IFOV,ITAN)
            RADTMP(ITAN) = RADTMP(ITAN) + FNCFOV(IFOV)*RADTAN(IFIN,JTAN)
            TRATMP(ITAN) = TRATMP(ITAN) + FNCFOV(IFOV)*TRATAN(IFIN,JTAN)
          END DO
        END DO
        NXTTAN = NTAN
        IF ( JACFLG ) THEN
          DO IJAC = 1, NJAC
            DO ITAN = 1, NTAN
              NXTTAN = NXTTAN + 1
              RADTMP(NXTTAN) = 0.0D0
              TRATMP(NXTTAN) = 0.0D0
              DO IFOV = 1, NFOV
                JTAN = IFVTAN(IFOV,ITAN)
                JTAN = ITNJAC(JTAN,IJAC)
                IF ( JTAN .NE. 0 ) THEN     ! 0=tan.path insensitive to ptb
                  RADTMP(NXTTAN) = RADTMP(NXTTAN) + 
     &                             FNCFOV(IFOV) * RADTAN(IFIN,JTAN)
                  TRATMP(NXTTAN) = TRATMP(NXTTAN) + 
     &                             FNCFOV(IFOV) * TRATAN(IFIN,JTAN)
                END IF
              END DO
            END DO
          END DO
        END IF
        DO ITAN = 1, LTAN
          RADTAN(IFIN,ITAN) = RADTMP(ITAN)
          TRATAN(IFIN,ITAN) = TRATMP(ITAN)
        END DO
      END DO
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
