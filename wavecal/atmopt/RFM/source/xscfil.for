C $Id: xscfil.for,v 1.9 1997/03/03 17:18:06 adarby Exp $
      SUBROUTINE XSCFIL ( LUNXSC, NAMXSC, FAIL, ERRMSG )
C
C VERSION
C     05-DEC-00  AD  Comment out READONLY
C     03-MAR-97  AD  Version 3.
C     24-DEC-96  AD  Add readonly to OPEN statement for VAX use
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Open X-Section data file read from RFM driver table.
C     Called by INPXSC for each file specified in *XSC section of driver table.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER       LUNXSC  ! I/O LUN for XSC File/next free LUN
      CHARACTER*(*) NAMXSC  !  I  Name of XSC file
      LOGICAL       FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  FILXFL ! Check X/S data input file and store in /XFLCOM/
     &, RFMLOG ! Write text message to RFM log file.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C LOCAL VARIABLES
      INTEGER       IOS         ! Saved value of IOSTAT for error messages
      LOGICAL       REJECT      ! T=ignore file - no useful data
      CHARACTER*132 MESSGE      ! Text message sent to LOG file
C
C EXECUTABLE CODE --------------------------------------------------------------
C
C Send message to LOG file saying which file is about to be opened
C
      MESSGE = 'I-XSCFIL: Opening X-Secn File: '//NAMXSC
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
C Open File 
C
      OPEN ( UNIT=LUNXSC, FILE=NAMXSC, STATUS='OLD', 
c             readonly,                                           ! VMS
     &        IOSTAT=IOS, ERR=900 )
C
C Check file is useful, and store files in *XFL 
C
      CALL FILXFL ( LUNXSC, REJECT, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
      IF ( REJECT ) THEN
        CLOSE ( LUNXSC )
      ELSE
        LUNXSC = LUNXSC + 1
      END IF
C
      FAIL = .FALSE.
      RETURN
C
  900 WRITE ( ERRMSG, '(A,I11)') 
     &  'F-XSCFIL: I/O failure on X/S file. IOSTAT=', IOS
      FAIL = .TRUE.
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
