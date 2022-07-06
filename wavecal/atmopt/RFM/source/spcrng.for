C $Id: spcrng.for,v 1.8 1997/03/03 17:17:58 adarby Exp $
      SUBROUTINE SPCRNG ( ISPC, RECORD, FAIL, ERRMSG )
C
C VERSION
C      23-APR-99  AD  Adapt for C*8 label instead of C*6
C      03-MAR-97  AD  Version 3.
C      01-OCT-96  AD  Version 2.
C      01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Check Spectral Range & Resln parameters and update table.
C
      IMPLICIT NONE 
C
C ARGUMENTS 
      INTEGER       ISPC    !  I  Index in table (range 1:NSPC)
      CHARACTER*(*) RECORD  !  I  Text record containing range & resln params
      LOGICAL       FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG  !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  RFMLOG ! Write text message to RFM log file.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmcon.inc' ! Constants for RFM calculations
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
C 
C LOCAL VARIABLES
      INTEGER          IOS    ! Saved value of IOSTAT for error message
      DOUBLE PRECISION WNOLOW ! Lower limit of Spc.range read from file [cm]
      DOUBLE PRECISION WNOUPP ! Upper limit of Spc.range read from file [cm]
      DOUBLE PRECISION WNORES ! Resolution of Spc.range read from file [cm]
      CHARACTER*8      LABEL  ! Spc.range label read from file
      CHARACTER*80     MESSGE ! Text for messages printed to Log file
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      LABEL = LABSPC(ISPC)
C
C Check that no calculation yet defined for this spectral range (if a non-zero
C resolution is specified then the range is already marked for use in this run)
C
      IF ( WNRSPC(ISPC) .NE. 0.0 ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-SPCRNG: Calculation already defined for '//
     &           'spectral range, Label: '//LABEL
        RETURN
      END IF
C
C Check that 3 numbers can be read successfully from record (actual values are
C checked later by SPCCHK)
C       
      READ ( RECORD, *, IOSTAT=IOS ) WNOLOW, WNOUPP, WNORES
      IF ( IOS .NE. 0 ) THEN
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,I11)' ) 
     &   'F-SPCRNG: Error reading Rng/Res params for Label='//
     &    LABEL//'. IOSTAT=', IOS
        RETURN
      END IF
C
C If range already defined - check if now redefining range and issue warning
C NB: OK to redefine ranges for which no calculation yet specified (WNRSPC=0.0)
C
      IF ( WNLSPC(ISPC) .NE. 0.0 .AND. WNUSPC(ISPC) .NE. 0.0 ) THEN  ! defined
        IF ( WNOLOW .NE. WNLSPC(ISPC) .OR. 
     &       WNOUPP .NE. WNUSPC(ISPC)       ) THEN
          MESSGE = 'W-SPCRNG: Redefining Spectral Range for Label: '//
     &      LABEL
          CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        END IF 
      END IF 
C
      WNLSPC(ISPC) = WNOLOW
      WNUSPC(ISPC) = WNOUPP
      WNRSPC(ISPC) = WNORES
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
