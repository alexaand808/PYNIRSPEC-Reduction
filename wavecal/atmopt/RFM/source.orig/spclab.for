C $Id: spclab.for,v 1.8 1997/03/03 17:17:56 adarby Exp $
      SUBROUTINE SPCLAB ( LABEL, IDXSPC, NEW, FAIL, ERRMSG )
C
C VERSION
C      22-JAN-98  AD  Split concatenation so FORTRAN to C converter can cope
C      03-MAR-97  AD  Version 3.
C      01-OCT-96  AD  Version 2.
C      01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Check Spc.range label find index in table.
C
      IMPLICIT NONE 
C
C ARGUMENTS 
      CHARACTER*(*) LABEL   !  I  Spectral range label.
      INTEGER       IDXSPC  !  O  Table index for label
      LOGICAL       NEW     !  O  T=this is a new label added to table
      LOGICAL       FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG  !  O  Error message written if FAIL is TRUE
C
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmcon.inc' ! Constants for RFM calculations
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
C 
C LOCAL VARIABLES
      INTEGER    ILIM      ! Pointer to last part of field written out
      INTEGER    ISPC      ! Counter for elements in Spc.range table
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      FAIL = .FALSE.
C
C Check Label within permitted length (LENLAB)
C
      IF ( LEN ( LABEL ) .GT. LENLAB ) THEN
        FAIL = .TRUE.
        ILIM = MIN ( LEN ( LABEL ), 24 )      ! 24 is max accomodated in ERRMSG
        WRITE ( ERRMSG, '(A,A,A,I2)' )
     &    'F-SPCLAB: Spc.range Label: ', LABEL(1:ILIM),
     &    ' exceeds max.length LENLAB=', LENLAB
          RETURN
      END IF
C  
C Find if label already listed - and either overwrite or add to list
C
      DO ISPC = 1, NSPC
        IF ( LABEL .EQ. LABSPC(ISPC) ) THEN 
          IDXSPC = ISPC
          NEW = .FALSE.
          RETURN
        END IF
      END DO
C
      NEW = .TRUE.
      IF ( NSPC .EQ. MAXSPC ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-SPCLAB: Array size MAXSPC in RFMSIZ.INC too small '
     &           //'to add new label: '//LABEL
      ELSE 
        NSPC = NSPC + 1
        IDXSPC = NSPC
        LABSPC(NSPC) = LABEL
        WNLSPC(NSPC) = 0.D0 
        WNUSPC(NSPC) = 0.D0 
        WNRSPC(NSPC) = 0.D0 
      END IF
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
