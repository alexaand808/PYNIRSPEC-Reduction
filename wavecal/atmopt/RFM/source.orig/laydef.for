C $Id: laydef.for,v 1.8 1997/03/03 17:17:28 adarby Exp $
      SUBROUTINE LAYDEF ( NNEW, HGTNEW, FAIL, ERRMSG )
C
C VERSION
C      03-MAR-97  AD  Version 3.
C      01-OCT-96  AD  Version 2.
C      01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Subdivide atmospheric profile layers.
C     Called by ATMLAY if neither OFM nor GL2 flags selected.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      NNEW         !  O  Number of layers in new profile
      REAL         HGTNEW(*)    !  O  Base Altitudes [km] of new profile layers
      LOGICAL      FAIL         !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG       !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  WIDGL2 ! Calculate typical average half-width of line 
      REAL WIDGL2
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes      
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
C
C LOCAL CONSTANTS
      REAL          TDIFB         ! Max temp variation across base layer [K]
        PARAMETER ( TDIFB = 5.0 )
      REAL          TDIFT         ! Max temp variation across base layer [K]
        PARAMETER ( TDIFT = 20.0 )
      REAL          WLIM          ! Max Voigt hw variation across layer [/cm]
        PARAMETER ( WLIM  = 1.5 )
C
C LOCAL VARIABLES
      INTEGER INEW           ! New profile layer counter
      INTEGER IATM           ! Atmospheric layer counter
      INTEGER JATM           ! Index of layer above IATM
      INTEGER JNEW           ! Secondary counter for shuffling new layers
      REAL    FI, FJ         ! Interpolation factors from layers IATM, JATM
      REAL    HB,HT          ! Alt [km] of bottom,top boundaries of new layer 
      REAL    HGTBOT         ! Altitude at bottom of old (&new) profiles [km]
      REAL    HGTTOP         ! Altitude at top of old (&new) profiles [km]
      REAL    HGTSPN         ! Altitude range spanned by old (&new) profiles[km]
      REAL    PT             ! Press [mb] of top boundary of new layer 
      REAL    TB,TT          ! Temp [mb] of bottom,top boundaries of new layer 
      REAL    TLIM           ! Maximum allowed temperature variation
      REAL    WB,WT          ! Voigt hwdth [cm] of bot,top bounds of new layer 
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      HGTBOT = HGTATM(1)
      HGTTOP = HGTATM(NATM)
      HGTSPN = HGTTOP - HGTBOT
C
      DO IATM = 1, NATM
        HGTNEW(IATM) = HGTATM(IATM)
      END DO
      NNEW = NATM
C
      HT = HGTATM(1) 
      PT = PREATM(1)
      TT = TEMATM(1)
      WT = WIDGL2 ( PT, TT ) 
C
      INEW = 1
      DO IATM = 1, NATM - 1
        JATM = IATM+1
        DO WHILE ( HGTNEW(INEW) .LT. HGTATM(JATM) )
          INEW = INEW + 1
          HB = HT
          TB = TT
          WB = WT
          TLIM = TDIFB * ( TDIFT / TDIFB )**( ( HB - HGTBOT ) / HGTSPN )
C
  200     CONTINUE
          HT = HGTNEW(INEW)
          FJ = ( HT - HGTATM(IATM) ) / ( HGTATM(JATM) - HGTATM(IATM) )
          FI = 1.0 - FJ
          TT = TEMATM(IATM) * FI + TEMATM(JATM) * FJ
          PT = EXP ( FI * LOG ( PREATM(IATM) ) + 
     &               FJ * LOG ( PREATM(JATM) )   )
          WT = WIDGL2 ( PT, TT ) 
C
          IF ( ABS ( TT - TB ) .GT. TLIM .OR.
     &         WT/WB .GT. WLIM .OR. WB/WT .GT. WLIM ) THEN
            IF ( NNEW .GE. MAXATM ) THEN
              FAIL = .TRUE.
              WRITE ( ERRMSG, '(A,I11)' )
     &          'F-LAYDEF: No profile levels reqd for sublayering '//
     &          '> MAXATM in RFMSIZ,=', MAXATM
              RETURN
            END IF
C 
            DO JNEW = NNEW, INEW, -1
              HGTNEW(JNEW+1) = HGTNEW(JNEW)
            END DO
            HGTNEW(INEW) = 0.5 * ( HT + HB )
            NNEW = NNEW + 1
            GOTO 200
          END IF 
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
