      SUBROUTINE WRTTAB ( FAIL, ERRMSG )
C
C VERSION
C      31-OCT-97  AD  Original.
C
C DESCRIPTION
C     Write RFM tabulated absorption coefficients.
C     Called by RFM for each widemesh interval if TAB flag selected.
C
      IMPLICIT NONE
C
C ARGUMENTS
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error occurs
      CHARACTER*80 ERRMSG  !  O  Error message written if FAIL is TRUE
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'outcom.inc' ! RFM output file data
      INCLUDE 'tabcom.inc' ! Table axes for tabulated absorb.coeffs.
C
C LOCAL VARIABLES
      INTEGER       ICLC   ! Calculated path counter
      INTEGER       IFIN   ! Fine mesh grid-point counter
      INTEGER       IGAS   ! Counter for absorbers
      INTEGER       IOFF   ! Path offset index for different gases
      INTEGER       IOS    ! Value of IOSTAT saved for error messages.
      INTEGER       LUN    ! Local value of Logical Unit Number.
      INTEGER       NTAB   ! No.tabulated (p,T) points for each wavenumber
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      FAIL = .FALSE.
      NTAB = NLPTAB * NTMTAB
C
      DO IGAS = 1, NGAS
        LUN = LUNTAB + IGAS - 1 
        IOFF = ( IGAS - 1 ) * NTAB
        DO IFIN = IOUT1, IOUT2
          IF ( BINFLG ) THEN
            WRITE ( LUN, IOSTAT=IOS, ERR=900 ) 
     &        ( ABSFIN(IFIN,ICLC), ICLC = IOFF+1, IOFF+NTAB )
          ELSE
            WRITE ( LUN, *, IOSTAT=IOS, ERR=900 )
     &        ( ABSFIN(IFIN,ICLC), ICLC = IOFF+1, IOFF+NTAB )
          END IF
        END DO
      END DO
C
      FAIL = .FALSE.
      RETURN
C
  900 FAIL = .TRUE.
      WRITE ( ERRMSG, '(A,I11)' ) 
     &  'F-WRTTAB: I/O failure on output file. IOSTAT=', IOS
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
