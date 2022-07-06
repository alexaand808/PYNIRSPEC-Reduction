      SUBROUTINE CHKRES ( WNRMAX, FAIL, ERRMSG )
C
C VERSION
C     07-AUG-13  AD  Change FORMAT descriptors to avoid ifort warnings
C     18-MAR-91  AD  Original. Replaces CHKWID
C
C DESCRIPTION
C     Check tangent pt line-widths compatible with spectral resln.
C     Called by FINCHK for highest tangent height or hom.path
C     NB should only be called if IATTAN(NTAN) is defined and non-zero.
C
      IMPLICIT NONE
C
C ARGUMENTS
      REAL         WNRMAX !  I  Maximum spacing of fine mesh resolution
      LOGICAL      FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  RFMLOG ! Write text message to RFM log file.
     &, WIDGL2 ! Calculate typical average half-width of line
      REAL  WIDGL2
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL CONSTANTS
      REAL          WIDLIM ! Min No. calc.pts/half width for adequate resoln
        PARAMETER ( WIDLIM = 1.5 )
C
C LOCAL VARIABLES
      INTEGER      IATM   ! Index of Atmospheric layer with highest tan.pt.
      REAL         PRE    ! Pressure [mb] corresponding to bottom of layer
      REAL         TEM    ! Temperature [K] corresponding to bottom of layer
      REAL         WID    ! Typical line half-width [cm-1]
      CHARACTER*80 MESSGE ! Warning message sent to runlog file
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      FAIL = .FALSE.
C
      IF ( NTAN .EQ. 0 ) STOP 'F-CHKRES: Logical Error#1'
      IF ( IATTAN(NTAN) .EQ. 0 ) STOP 'F-CHKRES: Logical Error#2'
C
      IATM = IATTAN(NTAN)
      TEM = TEMATM(IATM)
      PRE = PREATM(IATM)
      WID = WIDGL2 ( PRE, TEM )
      IF ( WID * WIDLIM .LT. WNRMAX ) THEN
        IF ( NATM .GT. 1 ) THEN
          WRITE ( MESSGE, '(A,F7.3,A,F7.4,A)' ) 
     &      'W-CHKRES: linewidth at Tan.Ht=', HGTATM(IATM),
     &      ' is ~', WID, ' cm-1 - will not be resolved'
        ELSE 
          WRITE ( MESSGE, '(A,F7.4,A)' ) 
     &      'W-CHKRES: linewidth is ~', WID, 
     &      ' cm-1 - will not be resolved'
        END IF
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG ) 
      END IF    
C
      END
