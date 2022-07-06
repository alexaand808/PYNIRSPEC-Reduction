C $Id: fovfil.for,v 1.8 1997/03/03 17:17:09 adarby Exp $
      SUBROUTINE FOVFIL ( LUNFOV, NAMFOV, FAIL, ERRMSG )
C
C VERSION
C     26-MAY-00  AD  Convert FNCFOV from SP to DP
C     18-JUL-98  AD  Interpret -ve NFOV as angular FOV data
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Open, read FOV data file and close.
C     Called by INPFOV if filename listed in *FOV section of Driver Table.
C
      IMPLICIT NONE 
C
C ARGUMENTS
      INTEGER       LUNFOV  !  I  LUN for FOV data file
      CHARACTER*(*) NAMFOV  !  I  Name of FOV data file
      LOGICAL       FAIL    !  O  T=A fatal error was detected
      CHARACTER*80  ERRMSG  !  O  Error message written if FAIL is TRUE
C      
      EXTERNAL
     &  OPNFIL ! Open ASCII input file, and skip/Log any initial comments.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'fovcom.inc' ! FOV Data
C
C LOCAL VARIABLES      
      INTEGER IFOV           ! Counter for FOV tabulation points
      INTEGER IOS            ! Saved value of IOSTAT for error message
      REAL    TABFOV(MAXFOV) ! Tabulated FOV response function
      DOUBLE PRECISION SUM   ! Summation for normalising FOV function
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      CALL OPNFIL ( LUNFOV, NAMFOV, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C      
      READ ( LUNFOV, *, IOSTAT=IOS, ERR=900 ) NFOV
      ELEFOV = NFOV .LT. 0
      NFOV = ABS ( NFOV )
C
      FAIL = .TRUE.
      IF ( NFOV .LE. 2 ) THEN
        WRITE ( ERRMSG, '(A,I2,A)' )
     &    'F-FOVFIL: No.tabulated FOV pts=', NFOV,
     &    ' - should be at least 3'
      ELSE IF ( NFOV .GT. MAXFOV ) THEN
        WRITE ( ERRMSG, '(A,I11,A,I11,A)' )
     &    'F-FOVFIL: No.tabulated FOV pts=', NFOV,
     &    ' > MAXFOV=', MAXFOV, ' in RFMSIZ.INC'
      ELSE IF ( ELEFOV .AND. .NOT. OBSFLG ) THEN
        ERRMSG = 'F-FOVFIL: Require Observation altitude (OBS Flag)'//
     &           ' to use angular FOV tabulation'
      ELSE IF ( .NOT. ELEFOV .AND. OBSFLG ) THEN
        ERRMSG ='F-FOVFIL: Require angular FOV tabulation with OBS flag'
      ELSE
        FAIL = .FALSE.
      END IF
      IF ( FAIL ) RETURN
C
      READ ( LUNFOV, *, IOSTAT=IOS, ERR=900 ) 
     &  ( ALTFOV(IFOV), IFOV = 1, NFOV )
      READ ( LUNFOV, *, IOSTAT=IOS, ERR=900 ) 
     &  ( TABFOV(IFOV), IFOV = 1, NFOV )
C
      CLOSE ( LUNFOV, IOSTAT=IOS, ERR=900 )
C
      IF ( TABFOV(1) .NE. 0.0 .OR. TABFOV(NFOV) .NE. 0.0 ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-FOVFIL: FOV Function non-zero at one or both ends'
        RETURN
      END IF
C
      FNCFOV(1) = DBLE ( ( ALTFOV(2) - ALTFOV(1) ) * TABFOV(2) )
      FNCFOV(NFOV) = DBLE ( ( ALTFOV(NFOV) - ALTFOV(NFOV-1) ) 
     &                      * TABFOV(NFOV-1) )
      SUM = FNCFOV(1) + FNCFOV(NFOV)
      DO IFOV = 2, NFOV - 1
        FNCFOV(IFOV) = 
     &    DBLE (   ( ALTFOV(IFOV) - ALTFOV(IFOV-1) ) 
     &           * ( TABFOV(IFOV-1) + 2 * TABFOV(IFOV) ) 
     &           + ( ALTFOV(IFOV+1) - ALTFOV(IFOV) ) 
     &           * ( TABFOV(IFOV+1) + 2 * TABFOV(IFOV) ) )
        SUM = SUM + FNCFOV(IFOV)
      END DO
C
      DO IFOV = 1, NFOV
        FNCFOV(IFOV) = FNCFOV(IFOV) / SUM
      END DO
C
      FAIL = .FALSE.
      RETURN
C
  900 FAIL = .TRUE.
      WRITE ( ERRMSG, '(A,I11)' ) 
     &  'F-FOVFIL: I/O failure on FOV file. IOSTAT=', IOS
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
