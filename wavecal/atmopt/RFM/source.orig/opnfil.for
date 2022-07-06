C $Id: opnfil.for,v 1.13 1997/03/03 17:17:35 adarby Exp $
      SUBROUTINE OPNFIL ( LUNFIL, NAMFIL, FAIL, ERRMSG )
C
C VERSION
C     05-DEC-00  AD  Comment out READONLY
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     18-SEP-96  AD  Don't convert filename to lower case before opening.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Open ASCII input file, and skip/Log any initial comments.
C     General purpose routine.
C     First record written to Log file if it begins with '!'
C     All subsequent records beginning with '!' are skipped
C     Sets file pointer ready to read 1st record not beginning with '!'.
C
      IMPLICIT NONE 
C
C ARGUMENTS
      INTEGER       LUNFIL !  I  LUN for opening/reading file
      CHARACTER*(*) NAMFIL ! I/O File name (converted to lower case internally)
      LOGICAL       FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  NXTREC ! Load next record from RFM input file.
     &, RFMLOG ! Write text message to RFM log file.
C
C LOCAL VARIABLES
      INTEGER       IOS    ! Saved value of IOSTAT for error messages
      LOGICAL       LDUMMY ! Set TRUE if end-of-section/file reached (ignored)
      CHARACTER*101 MESSGE ! Text message sent to LOG file
      CHARACTER*80  RECORD ! Text record read from file
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      MESSGE = 'I-OPNFIL: Opening file: '//NAMFIL   ! Max 80 chars for NAMFIL
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
      OPEN ( UNIT=LUNFIL, FILE=NAMFIL, STATUS='OLD', 
c            readonly,                                         ! VMS
     &       IOSTAT=IOS, ERR=900 )
C
C Print first record of file to Log file (should be a comment anyway)
C
      READ ( LUNFIL, '(A80)', IOSTAT=IOS, ERR=900 ) RECORD
      CALL RFMLOG ( '  '//RECORD, FAIL, ERRMSG )
      BACKSPACE ( LUNFIL, IOSTAT=IOS, ERR=900 ) 
C
C Set pointer to first non-comment, non-blank record
C
      CALL NXTREC ( LUNFIL, RECORD, LDUMMY, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
      BACKSPACE ( LUNFIL, IOSTAT=IOS, ERR=900 )
C
      RETURN
C
  900 FAIL = .TRUE.
      WRITE ( ERRMSG, '(A,I11)' )       
     &  'F-OPNFIL: I/O failure on file. IOSTAT=', IOS 
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
