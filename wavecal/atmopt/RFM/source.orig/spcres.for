C $Id: spcres.for,v 1.7 1997/03/03 17:17:57 adarby Exp $
      SUBROUTINE SPCRES ( RECORD, FAIL, ERRMSG )
C
C VERSION
C      23-APR-99  AD  Adapt for C*8 LABEL instead of C*6
C      03-MAR-97  AD  Version 3.
C      01-OCT-96  AD  Version 2.
C      01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Set calculation resolution for labelled spectral range.
C     Called by INPSPC for each record in *SPC section with 2 fields.
C     
      IMPLICIT NONE
C
C ARGUMENTS
      CHARACTER*80 RECORD  !  I  Record read from Driver Table
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  SPCLAB ! Check Spc.range label find index in table.
     &, TXTFLD ! Identify start and end points of text field in record
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
C
C LOCAL VARIABLES
      INTEGER          IDXSPC ! Table index of required spectral range 
      INTEGER          IEND   ! Pointer to end of field in text record
      INTEGER          IOS    ! Saved value of IOSTAT   
      INTEGER          ISTA   ! Pointer to start of field in text record
      LOGICAL          NEW    ! Set TRUE if label not previously listed
      DOUBLE PRECISION WNORES ! Spectral resolution for calculations
      CHARACTER*8      LABEL  ! Checked Label for required spectral range 
C
C EXECUTABLE CODE --------------------------------------------------------------
C
C Locate first field in record - should be Label - and check that a range has
C already been tabulated for it (ie not a NEW label)
C 
      CALL TXTFLD ( 80, RECORD, 1, ISTA, IEND )
      CALL SPCLAB ( RECORD(ISTA:IEND), IDXSPC, NEW, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      LABEL = LABSPC(IDXSPC)
      IF ( NEW ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-SPCRES: No pretabulated Spectral Range for Label: '
     &           //LABEL
        RETURN
      END IF
C
C Check that no calculation yet defined for this spectral range (if a non-zero
C resolution is specified then the range is already marked for use in this run)
C
      IF ( WNRSPC(IDXSPC) .NE. 0.0 ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-SPCRES: Resolution already defined for '//
     &           'Spectral Range, Label: '//LABEL
        RETURN
      END IF
C
C Read resolution from remainder of RECORD - value will be checked later
C
      READ ( RECORD(IEND+1:), *, IOSTAT=IOS ) WNORES
      IF ( IOS .NE. 0 ) THEN
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,I11)' ) 
     &    'F-SPCRES: Failed to read Resln for Spc.range '//
     &    'label='//LABEL//'. IOSTAT=', IOS
        RETURN
      END IF
C
      WNRSPC(IDXSPC) = WNORES
      FAIL = .FALSE.
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
