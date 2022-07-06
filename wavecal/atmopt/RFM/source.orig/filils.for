      SUBROUTINE FILILS ( LUNILS, WNOLOW, WNOUPP, FAIL, ERRMSG )
C
C VERSION
C     09-AUG-03  AD  Allow for GHz rather than Wno 
C     15-MAY-02  AD  Remove redundant rfmcon.inc
C     27-APR-00  AD  Convert PT2ILS to DP
C     07-OCT-99  AD  Remove checks on ILS width v. WIDILS (now in ILSCHK)
C                    Keep default ILS as last in array
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     30-SEP-96  AD  Correction to format of error messages
C     29-SEP-96  AD  Check width of user-spec ILS cf. WIDILS
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Check ILS Function and load /ILSCOM/
C     Called by ILSFIL.
C
      IMPLICIT NONE 
C
C ARGUMENTS 
      INTEGER          LUNILS !  I  Lun for ILS file
      DOUBLE PRECISION WNOLOW !  I  Lower wavenumber for ILS application
      DOUBLE PRECISION WNOUPP !  I  Upper wavenumber for ILS application
      LOGICAL          FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80     ERRMSG !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  RFMLOG ! Write text message to RFM log file
C
C GLOBAL CONSTANTS
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C LOCAL CONSTANTS
      DOUBLE PRECISION WNOFAC ! GHz to Wno conversion factor
        PARAMETER ( WNOFAC = 1.0D7 / VLIGHT )   ! VLIGHT in phycon.inc
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'ilscom.inc' ! Instrument Lineshape functions.
C 
C LOCAL VARIABLES
      INTEGER      IILS    ! Index for ILS arrays
      INTEGER      IOS     ! IOSTAT saved for error message
      INTEGER      IPT     ! Counter for reading Tabulation points
      LOGICAL      ANYDEF  ! Default ILS function loaded (can only be one)
      CHARACTER*80 WRNMSG  ! Warning message sent to log file
C
C DATA STATEMENTS
      DATA ANYDEF / .FALSE. /
      SAVE ANYDEF
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      FAIL = .FALSE.
C
C Check enough space for extra table
C 
      IF ( NILS .EQ. MAXILS ) THEN
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,I11)' ) 
     &    'F-FILILS: No space for extra ILS Function. '//
     &    'MAXILS in RFMSIZ.INC=', MAXILS
        RETURN
      END IF
      NILS = NILS + 1
C
C If default loaded, then keep this as last in array
      IF ( ANYDEF ) THEN
        NPTILS(NILS) = NPTILS(NILS-1)
        PT1ILS(NILS) = PT1ILS(NILS-1)
        PTDILS(NILS) = PTDILS(NILS-1)
        WNUILS(NILS) = WNUILS(NILS-1)       ! Should be 0.0
        WNLILS(NILS) = WNLILS(NILS-1)       ! Should be 0.0
        DO IPT = 1, NPTILS(NILS)
          FNCILS(IPT,NILS) = FNCILS(IPT,NILS-1)
        END DO
        IILS = NILS-1      ! Load next data into penultimate element
      ELSE
        IILS = NILS        ! Load next data into last element
      END IF
C
C Check default wavenumber range
C 
      IF ( WNOLOW .EQ. 0.0D0 .AND. WNOUPP .EQ. 0.0D0 ) THEN
        IF ( ANYDEF ) THEN
          FAIL = .TRUE. 
          ERRMSG = 'F-FILILS: No range specified for ILS Fn.'
     &      //' but default ILS Fn already loaded.'
          RETURN
        ELSE
          ANYDEF = .TRUE.
        END IF
      END IF
      WNLILS(IILS) = WNOLOW
      WNUILS(IILS) = WNOUPP
C
C Read No.pts, initial point and wavenumber spacing for ILS function
C
      READ ( LUNILS, *, IOSTAT=IOS, ERR=900 ) 
     &  NPTILS(IILS), PT1ILS(IILS), PTDILS(IILS)
C
C If NPts is negative, convert from GHz to Wno
      IF ( NPTILS(IILS) .LT. 0 ) THEN
        NPTILS(IILS) = -NPTILS(IILS)
        WNLILS(IILS) = WNLILS(IILS) * WNOFAC
        WNUILS(IILS) = WNUILS(IILS) * WNOFAC
        PT1ILS(IILS) = PT1ILS(IILS) * WNOFAC
        PTDILS(IILS) = PTDILS(IILS) * WNOFAC
        IF ( .NOT. GHZFLG ) THEN 
          WRNMSG = 'W-FILILS: ILS specified in GHz - converting to Wno'
          CALL RFMLOG ( WRNMSG, FAIL, ERRMSG )
        END IF
      ELSE      ! If NPt positive, ILS in Wno, so warn if RFM calc in GHz
        IF ( GHZFLG ) THEN 
          WRNMSG = 'W-FILILS: ILS specified in Wno - converting to GHz'
          CALL RFMLOG ( WRNMSG, FAIL, ERRMSG )
        END IF
      END IF
      IF ( FAIL ) RETURN
C
      PT2ILS(IILS) = PT1ILS(IILS) + ( NPTILS(IILS) - 1 ) * PTDILS(IILS)
      IF ( NPTILS(IILS) .GT. MAXILP ) THEN
        FAIL = .TRUE.
        WRITE ( ERRMSG, '(A,I11,A,I11)' ) 
     &    'F-FILILS: ILS Fn. has ', NPTILS(IILS),
     &    ' Pts, Max in RFMSIZ.INC is MAXILP=', MAXILP
        RETURN
      END IF
C
      READ ( LUNILS, *, IOSTAT=IOS, ERR=900 ) 
     &  ( FNCILS(IPT,IILS), IPT = 1, NPTILS(IILS) )
      FAIL = .FALSE.
      RETURN
C
  900 FAIL = .TRUE.
      WRITE ( ERRMSG, '(A,I11)' )
     &  'F-FILILS: I/O error on ILS data. IOSTAT=', IOS
      END
