C $Id: spcfil.for,v 1.9 1997/03/03 17:17:55 adarby Exp $
      SUBROUTINE SPCFIL ( LUNSPC, NAMSPC, FAIL, ERRMSG )
C
C VERSION
C      22-JAN-98  AD  Remove redundant EXTERNAL declaration: RFMLOG
C      03-MAR-97  AD  Version 3.
C      01-OCT-96  AD  Version 2.
C      01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Open file and load pretabulated spectral range data.
C     Called by INPSPC for each record in *SPC section with 1 field.
C
      IMPLICIT NONE 
C
C ARGUMENTS 
      INTEGER       LUNSPC  !  I  Logical unit number for opening file
      CHARACTER*(*) NAMSPC  !  I  Name of spectral range data file
      LOGICAL       FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG  !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  NXTREC ! Load next record from RFM input file
     &, OPNFIL ! Open ASCII input file, and skip/Log any initial comments.
     &, SPCLAB ! Check Spc.range label find index in table.
     &, SPCRNG ! Check Spectral Range & Resln parameters and update table.
     &, TXTFLD ! Identify start and end points of text field in record
C
C LOCAL VARIABLES
      INTEGER      IDXSPC    ! Index of spectral range in tabulated data
      INTEGER      ISTA,IEND ! Pointers to start,end of field in RECORD
      LOGICAL      ENDSEC    !  TRUE = end of section or file
      LOGICAL      LDUMMY    ! TRUE = new spectral range added (ignored here)
      CHARACTER*80 RECORD    ! Record read from Spc.range file
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      CALL OPNFIL ( LUNSPC, NAMSPC, FAIL, ERRMSG )
      IF ( FAIL ) RETURN                    
C
C Load each record from Spc.range file
C
  100 CALL NXTREC ( LUNSPC, RECORD, ENDSEC, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      IF ( .NOT. ENDSEC ) THEN            
C
C Extract Label and index in internal table of spectral ranges 
C
        CALL TXTFLD ( 80, RECORD, 1, ISTA, IEND )
        CALL SPCLAB ( RECORD(ISTA:IEND), IDXSPC, LDUMMY, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
C
C Pass remainder of record for extraction of range and resolution parameters
C
        CALL SPCRNG ( IDXSPC, RECORD(IEND+1:), FAIL, ERRMSG )
        IF ( FAIL ) RETURN
C
C Get next record from Spectral Range file
C
        GOTO 100 
      END IF
C
      CLOSE ( LUNSPC )
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
