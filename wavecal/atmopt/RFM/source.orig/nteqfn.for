C $Id: nteqfn.for,v 1.5 1997/03/03 17:17:33 adarby Exp $
      SUBROUTINE NTEQFN ( IDMOL, IDISO, IATMLO, 
     &                    NLEV, LNPLEV, QFNLEV, FAIL, ERRMSG )
C
C VERSION
C     20-FEB-01  AD  Duplicte VTs to top of atmospheric profile if necessary
C     01-JUL-98  AD  Remove IATMHI argument. Assume NTE profile always extends
C                    to top of ATM profile
C     04-APR-97  AD  Remove FIRST argument - assume NQFN initialised elsewhere.
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     21-SEP-96  AD  Original.
C
C DESCRIPTION
C     Store Vib.Partition Function data read from NTE file. 
C     Called by NTEFIL for each partitition function in each file.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER       IDMOL      !  I  HITRAN ID for molecule 
      INTEGER       IDISO      !  I  Isotope#
      INTEGER       IATMLO     !  I  Lower index of ATM profile for interpolatn.
      INTEGER       NLEV       !  I  No.of profile levels
      REAL          LNPLEV(*)  !  I  ln(p/mb) profile
      REAL          QFNLEV(*)  !  I  Partition Function profile
      LOGICAL       FAIL       !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG     !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  LOOKUP !  General purpose interpolation routine
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'qfncom.inc' ! non-LTE Partition Functions 
C
C LOCAL VARIABLES
      INTEGER      IATM    ! Counter for atmospheric profile layers
      INTEGER      IGAS    ! Index of GAS in RFM list /GASCOM/
      INTEGER      IQFN    ! Counter for Stored Partition Function profiles
      INTEGER      J       ! Saved index for interpolation LOOKUP    
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      FAIL   = .FALSE.
C
C Check if molecule actually required
C
      IF ( IGSMOL(IDMOL) .EQ. 0 ) RETURN
C
C Check partition function for molecule,isotope not already loaded
C 
      DO IQFN = 1, NQFN
        IF ( IDGQFN(IQFN) .EQ. IDMOL .AND.
     &       ISOQFN(IQFN) .EQ. IDISO       ) THEN
          FAIL = .TRUE.
          IGAS = IGSMOL(IDMOL)
          WRITE ( ERRMSG, '(A,I11)' )
     &      'F-NTEQFN: Repeated NTE Vib.Partition Fn for Gas='//
     &      CODGAS(IGAS)//' Iso#=', IDISO
          RETURN
        END IF
      END DO
C
C Check sufficient space to add to list
C
      IF ( NQFN .EQ. MAXQFN ) THEN
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,I11)' ) 
     &    'F-NTEQFN: No.Different Vib.Part.Fns > '//
     &    'MAXQFN in RFMSIZ.INC,=', MAXQFN
        RETURN
      END IF
C
C Increment NQFN and add molecule,isotope to list 
C
      NQFN = NQFN + 1
      IDGQFN(NQFN) = IDMOL
      ISOQFN(NQFN) = IDISO
C
      J = 1
      DO IATM = 1, IATMLO-1
        PRFQFN(IATM,NQFN) = 1.0
      END DO
      DO IATM = IATMLO, NATM
        IF ( LNPATM(IATM) .LT. LNPLEV(NLEV) ) THEN
          PRFQFN(IATM,NQFN) = QFNLEV(NLEV)
        ELSE
          CALL LOOKUP ( LNPATM(IATM), PRFQFN(IATM,NQFN), LNPLEV,
     &                        QFNLEV, NLEV, J, 1 )
        END IF
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
