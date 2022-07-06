C $Id: gasfil.for,v 1.7 1997/03/03 17:17:10 adarby Exp $
      SUBROUTINE GASFIL ( LUNGAS, NAMGAS, FAIL, ERRMSG )
C
C VERSION
C      03-MAR-97  AD  Version 3.
C      01-OCT-96  AD  Version 2.
C      01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Open file and load pretabulated gases.
C     Called by INPGAS if field in *GAS section not recognised as a gas.
C     I.e the field is presumed to be a filename.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER       LUNGAS  !  I  LUN for Gas File
      CHARACTER*(*) NAMGAS  !  I  Name of Gas file
      LOGICAL       FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  GASCHK ! Check Gas name and set indices.
     &, NXTREC ! Load next record from RFM input file.
     &, OPNFIL ! Open ASCII input file, and skip/Log any initial comments.
     &, TXTFLD ! Identify start and end points of text field in record
C
C LOCAL VARIABLES
      INTEGER      IFLD           ! Counter for fields in record
      INTEGER      IEND           ! End location for text identifying gas
      INTEGER      ISTA           ! Starting location for text identifying gas
      LOGICAL      ENDSEC         ! Set TRUE when end of section or file reached
      LOGICAL      NOTGAS         ! Set TRUE if file gas not recognised
      CHARACTER*80 RECORD         ! Text record read from driver table
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      CALL OPNFIL ( LUNGAS, NAMGAS, FAIL, ERRMSG )
      IF ( FAIL ) THEN
        IEND = MIN ( LEN(NAMGAS), 29 )  ! Prevents error message overflow
        ERRMSG = 'F-GASFIL: Unidentified gas/file:'//NAMGAS(1:IEND)
     &          //ERRMSG(62:80)
        RETURN
      END IF
C     
  100 CALL NXTREC ( LUNGAS, RECORD, ENDSEC, FAIL, ERRMSG ) 
      IF ( FAIL ) RETURN
      IF ( .NOT. ENDSEC ) THEN 
C
C Loop here for each absorbing species in record
C
        IFLD = 0
  200   IFLD = IFLD + 1
        CALL TXTFLD ( 80, RECORD, IFLD, ISTA, IEND )    ! Isolate gas string
        IF ( ISTA .NE. 0 ) THEN                         ! Found next gas
          CALL GASCHK ( RECORD(ISTA:IEND), NOTGAS, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
          IF ( NOTGAS ) THEN
            FAIL = .TRUE.
            IEND = MIN ( IEND, ISTA+37 )
            ERRMSG = 'F-GASFIL: File contained unrecognised gas:'//
     &                       RECORD(ISTA:IEND)
            RETURN
          END IF 
          GOTO 200
        ELSE
          GOTO 100
        END IF
      END IF
C
      CLOSE ( LUNGAS )
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
