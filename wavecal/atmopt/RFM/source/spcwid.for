      SUBROUTINE SPCWID ( ISPC )
C
C VERSION
C     28-NOV-01  AD  Correction to setting line-mixing limits.
C     16-MAR-01  AD  Add 0.50001 instead of 0.5 in calc of NWID to ensure 
C                    extra wm interval used in borderline cases
C     05-JAN-99  AD  Allow for AVG option.
C     07-OCT-99  AD  Set DW1ILS,DW2ILS=DELWID rather than WIDILS
C     11-AUG-99  AD  Make D.P. operations explicit
C     17-DEC-97  AD  Use actual ILS width to set widemesh limits
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-OCT-96  AD  Update with new ranges for CO2 Line mixing bands
C     19-SEP-96  AD  Add checks to ensure entire line-mixing sets are included
C                    Only add ILS margin if ILS option selected
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Set Wide Mesh grid parameters from spectral range/resln.
C     Called by SPCCHK and RFMSPC for each required spectral range.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      ISPC    !  I  Index of spectral range 
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'rfmcon.inc' ! Constants for RFM calculations
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'ilscom.inc' ! Instrument Lineshape functions.
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
      INCLUDE 'widcom.inc' ! Wide mesh data
C
C LOCAL CONSTANTS
      INTEGER       MAXSET         ! No. of sets of line-mixing data
        PARAMETER ( MAXSET = 12 )   ! Should match value in YMIX.FOR
C
C LOCAL VARIABLES
      INTEGER          IILS    ! Index for ILS function
      INTEGER          ISET    ! Mixing set counter
      DOUBLE PRECISION DW1ILS  ! Lower Wno adjustment [cm-1] for ILS convol.
      DOUBLE PRECISION DW2ILS  ! Upper Wno adjustment [cm-1] for ILS convol.
      REAL             WNLSET(MAXSET)   ! Lower Wno boundaries of each set
      REAL             WNUSET(MAXSET)   ! Upper Wno boundaries of each set
C
      DATA WNLSET / 615.9, 648.4,  662.3,  667.3,  667.7,  718.2,
     &              733.9, 791.4, 1932.4, 2076.8, 2093.3, 2128.3 /
      DATA WNUSET / 618.1, 651.1,  664.9,  670.1,  674.5,  720.8,
     &              741.8, 793.9, 1936.1, 2079.8, 2094.6, 2129.8 /
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Calculation of DELWID should be the same as in ILSCHK.
      IF ( WNRSPC(ISPC) .GT. 1.0D0 ) THEN        ! Treat as N.Pts/Wno
        DELWID = 1.D0
      ELSE
        DELWID = WNRSPC(ISPC) * NINT ( 1.0D0 / WNRSPC(ISPC) ) 
      END IF
      IF ( ILSFLG ) THEN
        IILS = ILSSPC(ISPC)
        IF ( IILS .EQ. 0 ) THEN   ! data not yet loaded, use default widths
          DW1ILS = -DELWID
          DW2ILS =  DELWID
        ELSE                      ! ILS data now loaded, so use actual widths
          DW1ILS = PT1ILS(IILS) 
          DW2ILS = PT2ILS(IILS)
        END IF
      ELSE IF ( AVGFLG ) THEN    ! Increase by one output grid interval
        IF ( WNRSPC(ISPC) .GT. 1.0D0 ) THEN
          DW1ILS = -1.0D0/WNRSPC(ISPC)
        ELSE
          DW1ILS = -WNRSPC(ISPC)
        END IF
        DW2ILS = -DW1ILS
      ELSE
        DW1ILS = 0.D0
        DW2ILS = 0.D0
      END IF
C
C Determine no. wide mesh grid points required (NB: IWID = 0,1 .... NWID so
C NWID+1 grid points, NWID intervals).
C
      WN1WID = INT ( WNLSPC(ISPC) + DW1ILS )
      NWID   = NINT ( ( WNUSPC(ISPC) + DW2ILS  
     &                  - WN1WID ) / DELWID + 0.50001D0 ) 
      WN2WID = WN1WID + NWID * DELWID
C
C Set limits for including lines in widemesh calculations
C
      WNLWID = WN1WID - FWIND
      WNUWID = WN2WID + FWIND
C
C If line-mixing used, it is important that all lines within a set are 
C included, so check that lower,upper limits don't split one of the sets - 
C if so, increase lower,upper margin by 1 cm-1 and retest.
      IF ( MIXFLG ) THEN
 100    CONTINUE        ! Repeat from here if line-mixing margins are increased
        DO ISET = 1, MAXSET
          IF ( WNLWID .GE. WNLSET(ISET) .AND. 
     &         WNLWID .LE. WNUSET(ISET) ) THEN
            WNLWID = WNLWID - 1.0D0
            GOTO 100
          ELSE IF ( WNUWID .GE. WNLSET(ISET) .AND. 
     &              WNUWID .LE. WNUSET(ISET) ) THEN
            WNUWID = WNUWID + 1.0D0
            GOTO 100
          END IF
        END DO
      END IF
C
      END
