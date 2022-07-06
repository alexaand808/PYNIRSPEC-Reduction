      SUBROUTINE TABPTH ( FAIL, ERRMSG )
C
C VERSION
C      29-APR-99  AD  Allow for temperature offset
C      14-JUL-98  AD  Comment change only.
C      22-OCT-97  AD  Original. 
C
C DESCRIPTION
C     Set up table-entry paths for RFM calculations.
C     Called by INPTAN if TAB flag enabled.
C     A path is defined for each absorber, pressure and temperature.
C
      IMPLICIT NONE
C
C ARGUMENTS
      LOGICAL      FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG !  O  Contains error message if FAIL is TRUE.
C
      EXTERNAL
     &  LOOKUP ! General purpose interpolation routine
     &, TEMOFF ! Temperature offset for pretabulated (p,T) data.
      REAL TEMOFF
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data 
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'pthcom.inc' ! RFM path data
      INCLUDE 'tabcom.inc' ! Table axes for tabulated absorb.coeffs.
C
C LOCAL VARIABLES
      INTEGER  IGAS   !  Gas counter
      INTEGER  IPRE   !  Pressure tabulation domain counter
      INTEGER  IPTH   !  Path segment counter
      INTEGER  ITEM   !  Temperature tabulation domain counter
      INTEGER  J      !  Starting guess for LOOKUP 
      REAL     LNP    !  Ln(p/mb)
      REAL     VMR    !  Volume mixing ratio [ppv]
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      NPTH = NGAS * NTMTAB * NLPTAB
      IF ( NPTH .GT. MAXPTH ) THEN
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,I11,A,I7)' )
     &    'F-TABPTH: No.paths reqd=', NPTH,
     &    ' >MAXPTH in RFMSIZ=', MAXPTH
        RETURN
      END IF
C
      NCLC = NPTH
      IF ( NCLC .GT. MAXCLC ) THEN
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,I11,A,I7)' )
     &    'F-TABPTH: No.paths for line-by-line calc=', NCLC,
     &    ' >MAXCLC in RFMSIZ=', MAXCLC
        RETURN
      END IF
C
      IPTH = 0
      J = 1
      DO IGAS = 1, NGAS
        DO ITEM = 1, NTMTAB
          DO IPRE = 1, NLPTAB
            IPTH = IPTH + 1
            IGSPTH(IPTH) = IGAS
            PREPTH(IPTH) = EXP ( -LP1TAB -(IPRE-1) * LPDTAB ) / ATMB
            TEMPTH(IPTH) = TM1TAB + ( ITEM - 1 ) * TMDTAB
            IF ( OFFTAB ) THEN
              LNP = - LP1TAB - (IPRE-1) * LPDTAB 
              TEMPTH(IPTH) = TEMPTH(IPTH) + TEMOFF ( LNP )
            END IF
            IF ( HOMFLG ) THEN
              PPAPTH(IPTH) = VMRATM(1,IGAS) * PREPTH(IPTH)
            ELSE
              LNP = LOG ( ATMB * PREPTH(IPTH) )
              CALL LOOKUP ( LNP, VMR, LNPATM, VMRATM(1,IGAS), NATM, J,1)
              PPAPTH(IPTH) = VMR * PREPTH(IPTH)
            END IF
            AMTPTH(IPTH) = 1.0E-7 ! ! 10^4 for cm^2 to m^2, 10^3 kmoles to moles
            CLCPTH(IPTH) = .TRUE.
          END DO
        END DO
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
