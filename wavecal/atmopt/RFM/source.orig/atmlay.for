C $Id: atmlay.for,v 1.8 1997/03/03 17:17:00 adarby Exp $
      SUBROUTINE ATMLAY ( FAIL, ERRMSG )
C
C VERSION
C     31-DEC-99  AD  Insert extra layers into horiz.gradient profiles as well
C     24-JAN-99  AD  Replace OFMFLG by LINFLG.
C                    Remove GL2FLG and module LAYGL2
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Subdivide atmospheric profile layers.
C     Called by INPATM once if LAYFLG option selected.
C     Three layering algorithms are available according to whether the GL2 or
C     OFM or neither option is chosen.
C
      IMPLICIT NONE
C
C ARGUMENTS
      LOGICAL      FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  BRAKET ! Return lower index I of the two points in array ARRAY(N)
     &, LAYDEF ! Sub-layer profile using default RFM algorithm
     &, MOVGRA ! Exchange profiles between /ATMCOM/ and /GRACOM/ arrays.
     &, RFMLOG ! Write text message to RFM log file.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes      
      INCLUDE 'phycon.inc' ! Math and Physical constants
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'gracom.inc' ! Atmospheric 2-D field data
C
C LOCAL VARIABLES
      INTEGER IGAS           ! Absorber counter
      INTEGER IGRA           ! 2-D profile location counter
      INTEGER INEW           ! New profile layer counter
      INTEGER IOLD           ! Original profile layer counter
      INTEGER JOLD           ! Index of Layer immediately above IOLD
      INTEGER NNEW           ! No.layers in new profiles
      INTEGER NOLD           ! No.layers in original profiles
      REAL    FI, FJ         ! Interpolation factors for values at IOLD, JOLD
      REAL    HGTNEW(MAXATM) ! New profile altitudes [km]
      REAL    HGTOLD(MAXATM) ! Original profile altitudes [km]
      REAL    PREOLD(MAXATM) ! Original profile pressures [mb]
      REAL    TEMOLD(MAXATM) ! Original profile temperatures [K]
      REAL    VMROLD(MAXATM,MAXGAS) ! Original profile vmrs [ppv]
      CHARACTER*80 MESSGE    ! Text message sent to Log file
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Save original profiles of height, pressure, temperature and vmr.
C
      NOLD = NATM
      IGRA = 0
C Do Height separately since don't want HGTOLD modified when routine repeats
C from 100 for other profiles in horizontal gradient option.
      DO IOLD = 1, NOLD
        HGTOLD(IOLD) = HGTATM(IOLD)
      END DO

 100  CONTINUE
      DO IOLD = 1, NOLD 
        PREOLD(IOLD) = PREATM(IOLD)
        TEMOLD(IOLD) = TEMATM(IOLD)
        DO IGAS = 1, NGAS
          VMROLD(IOLD,IGAS) = VMRATM(IOLD,IGAS)
        END DO
      END DO
C
C Perform sub-layering  
C
      IF ( IGRA .EQ. 0 ) THEN
        CALL LAYDEF ( NNEW, HGTNEW, FAIL, ERRMSG )
        IF ( FAIL ) RETURN        ! only fail if reqd NNEW > MAXATM
      END IF
C
      IF ( NNEW .EQ. NATM ) THEN
        MESSGE = 'I-ATMLAY: Atmospheric profile layers unchanged'
      ELSE 
        WRITE ( MESSGE, '(A,I11,A,I11)' )
     &    'I-ATMLAY: Number of Atmospheric layers changed from ', NATM,
     &    ' to ', NNEW
      END IF
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
      IF ( FAIL ) RETURN        
      IF ( NNEW .EQ. NATM ) RETURN
C
C Interpolate old profiles on to new altitude levels
C
      DO INEW = 2, NNEW         ! NB all layering leaves element 1 unchanged
        HGTATM(INEW) = HGTNEW(INEW)
        CHGATM(INEW) = .TRUE.
        CALL BRAKET ( HGTOLD, NOLD, HGTNEW(INEW), IOLD )
        IOLD = MIN ( NOLD-1, IOLD )
        JOLD = IOLD + 1
        FJ = ( HGTNEW(INEW) - HGTOLD(IOLD) ) / 
     &       ( HGTOLD(JOLD) - HGTOLD(IOLD) )
        FI = 1.0 - FJ
        TEMATM(INEW) = FI * TEMOLD(IOLD) + FJ * TEMOLD(JOLD)
        PREATM(INEW) = EXP ( FI * LOG ( PREOLD(IOLD) ) +
     &                       FJ * LOG ( PREOLD(JOLD) )   )
        DO IGAS = 1, NGAS
          IF ( LINFLG .OR. VMROLD(IOLD,IGAS) .LT. ARGMIN 
     &                .OR. VMROLD(JOLD,IGAS) .LT. ARGMIN ) THEN ! Interp. VMR
            VMRATM(INEW,IGAS) = FI * VMROLD(IOLD,IGAS) +
     &                          FJ * VMROLD(JOLD,IGAS)
          ELSE                      ! Linear interpolation in Log(vmr)
            VMRATM(INEW,IGAS) = EXP ( FI * LOG ( VMROLD(IOLD,IGAS) ) +
     &                                FJ * LOG ( VMROLD(JOLD,IGAS) )   )
          END IF
        END DO
      END DO
C
C If working with horizontal gradients relayer each profile in turn
      IF ( GRAFLG ) THEN
C First copy current profile in atmcom back to gracom
        NATM = NNEW                       ! move NNEW levels
        IF ( IGRA .EQ. 0 ) THEN           ! represents ref (psi=0) profile
          CALL MOVGRA ( 0, NOMGRA )
        ELSE                              ! other profiles
          CALL MOVGRA ( 0, IGRA )
        END IF
        IGRA = IGRA + 1                   ! Get next profile
        IF ( IGRA .EQ. NOMGRA ) IGRA = IGRA + 1  ! Already done ref. profile
C Then load next profile from gracom into atmcom.
        IF ( IGRA .LE. NGRA ) THEN        
          NATM = NOLD                     ! move NOLD levels
          CALL MOVGRA ( IGRA, 0 )
          GOTO 100
        ELSE
          NATM = NNEW
          CALL MOVGRA ( NOMGRA, 0 )       ! Reload ref. profile into /ATMCOM/
        END IF
      END IF
C
      NATM = NNEW                         ! Set NATM if GRAFLG not enabled
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
