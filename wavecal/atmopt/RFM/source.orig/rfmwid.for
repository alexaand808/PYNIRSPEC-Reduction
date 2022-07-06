      SUBROUTINE RFMWID ( FAIL, ERRMSG )
C
C VERSION
C     20-JUL-05  AD  Only call MIXSHP for CO2 lines.
C     13-MAR-04  AD  Set SUBWNG according to default isotope
C     25-FEB-04  AD  MORSE: Add check on CLCPTH
C     11-AUG-03  AD  Add VVWSHP, VVWCOR
C     12-JUN-02  AD  Simplify: assume only one HITRAN file
C     16-FEB-02  AD  Allow for different isotopes
C     22-AUG-01  AD  Change error message for unrecognised lineshape value
C     13-APR-00  AD  Add counts for NILWID and NOLWID
C     04-NOV-99  AD  Check for CTMGAS rather than CTMFLG before SUBWNG
C     09-AUG-99  AD  Use explicit DBLE, SNGL in CNTADJ,ANTADJ summations.
C     24-JAN-99  AD  Remove GL2FLG
C     03-MAR-97  AD  Version 3.
C     18-DEC-96  AD  Subtract offsets if CKD H2O Continuum applied
C     18-OCT-96  AD  Set absorption=0.0 beyond 25cm-1 from line centre 
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION  
C     Perform wide mesh absorption calculations.
C     Called once by RFM.
C     Make a wide pass over the spectrum frequency grid. Each line is read in 
C     turn and the appropriate paths identified. The absorption in freqency 
C     interval some distance off due to the line wings is calculated at three
C     points within the interval for future quadratic interpolation. 
C     Whether or not the line wing absorption in a given frequency interval is
C     treated in this way is determined by the frequency windows in operation. 
C
      IMPLICIT NONE
C
C ARGUMENTS      
      LOGICAL      FAIL      !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG    !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  ADJUST ! Adjust line data parameters for path conditions
     &, AWIDTH ! Calculate average half-width of line
     &, CHISHP ! Calculate Voigt Line shape allowing for chi-factor.
     &, DOPSHP ! Calculate Doppler Line shape.
     &, INIHFL ! Initialise the HITRAN line data file
     &, LORSHP ! Calculate Lorentz line shape
     &, MIXSHP ! Calculate Voigt Line shape allowing for line-mixing
     &, REAHIT ! Read HITRAN line data file for RFM 
     &, SUB25W ! Subtract value at 25cm-1 from abs.coeff.
     &, VOISHP ! Calculate Voigt line shape
     &, VVWCOR ! Apply Van Vleck-Weisskopf correction to line shape
     &, VVWSHP ! Calculate Van Vleck-Weisskopf Line shape.
      REAL    AWIDTH
C
C GLOBAL CONSTANTS
      INCLUDE 'idxcon.inc' ! HITRAN molecular indices
      INCLUDE 'rfmsiz.inc' ! MIPAS RFM Array sizes
      INCLUDE 'rfmcon.inc' ! Constants for RFM calculations
      INCLUDE 'shpcon.inc' ! Line Shape codes
C
C COMMON VARIABLES
      INCLUDE 'adjcom.inc' ! Path-Adjusted line data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'hitcom.inc' ! HITRAN line database variables
      INCLUDE 'outcom.inc' ! RFM output file data
      INCLUDE 'pthcom.inc' ! RFM path data
      INCLUDE 'rejcom.inc' ! Minimum Line strength limits
      INCLUDE 'widcom.inc' ! Widemesh pass data
C
C LOCAL VARIABLES
      INTEGER     IDGPTH      ! HITRAN Gas ID for path
      INTEGER     IGAS        ! Absorber counter
      INTEGER     IPTH        ! Path counter
      INTEGER     IQAD        ! Parabolic point counter for each WM interval
      INTEGER     ISHP        ! Lineshape index for gas in path
      INTEGER     IWD2        ! Half-Mesh grid point index
      INTEGER     IWD2L,IWD2U ! Low/Upp Half-Mesh points incl.line contrib.
      INTEGER     JEXCL,JEXCU ! Closest Low/Upp Wide Mesh Intvl including line 
      INTEGER     JWID        ! Wide mesh interval counter
      INTEGER     JWIDL,JWIDU ! Low/Upp Wide Mesh Intvls incl.line contrib.
      INTEGER     NWD2        ! No. Half-Mesh points spanned by line
      LOGICAL     EOF         ! T = end-of-file found
      LOGICAL     COUNTD(MAXWID) ! Set TRUE if line already counted for WMItvl.
      LOGICAL     COUNTL      ! Set TRUE if line already counted for total.
      LOGICAL     SUBWNG      ! T = Subtract abs.coeff at 25cm-1
      REAL        AVWDTH      ! Average line width
      DOUBLE PRECISION WNOLOW ! Lower wavenumber of region of line influence 
      DOUBLE PRECISION WNOUPP ! Upper wavenumber of region of line influence 
      REAL        ABSORP(0:MAXWD2) ! Absorption
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Initiate reading of line data
C 
      CALL INIHFL ( WNLWID, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
C
  100 CALL REAHIT ( EOF, FAIL, ERRMSG )       ! Get next line from file
      IF ( FAIL ) RETURN
      IF ( WNUM .LE. WNUWID .AND. .NOT. EOF ) THEN
        JWIDL = MAX (    1, INT( 1.0D0+(WNUM-FWIND-WN1WID)/DELWID )  )
        JWIDU = MIN ( NWID, INT( 1.0D0+(WNUM+FWIND-WN1WID)/DELWID )  )
        JEXCL = MIN ( NWID, INT( 1.0D0+(WNUM-FEXC-WN1WID)/DELWID ) -1)
        JEXCU = MAX (    1, INT( 1.0D0+(WNUM+FEXC-WN1WID)/DELWID ) +1)
        IWD2L = 2 * ( JWIDL - 1 )
        IWD2U = 2 * JWIDU 
        NWD2  = IWD2U - IWD2L + 1
        COUNTL = .FALSE.
        DO JWID = 1, NWID
          COUNTD(JWID) = .FALSE.
        END DO
        DO IPTH = 1, NCLC                         ! Loop over paths 
C
C Does this path use this line
          IGAS = IGSPTH(IPTH)
          IDGPTH = IDXGAS(IGAS)
          IF ( CLCPTH(IPTH) .AND. IDGPTH .EQ. IDGAS .AND. 
     &         ISOGAS(ISO,IGAS) .EQ. IGAS ) THEN
C
C Calculate path adjusted linshape parameters 
            SUBWNG = CTMGAS(ISOGAS(0,IGAS)) .AND. ( IDGPTH .EQ. IDXH2O ) 
            ISHP = SHPGAS(IGAS)
            CALL ADJUST ( IPTH )
            IF ( STRADJ .GE. WIDREJ ) THEN        ! Line not too weak
C
C Loop over coarse freq. intervals
              IF ( WNGFLG ) THEN
                AVWDTH = AWIDTH ( IPTH, WNUM )
                WNOLOW = WNUM - MIN ( FWIND, DBLE(NWDCUT*AVWDTH) )
                WNOUPP = WNUM + MIN ( FWIND, DBLE(NWDCUT*AVWDTH) )
                JWIDL = MAX ( 1,    1+INT( (WNOLOW-WN1WID) / DELWID ))
                JWIDU = MIN ( NWID, 1+INT( (WNOUPP-WN1WID) / DELWID ))
                IWD2L = 2 * ( JWIDL - 1 )
                IWD2U = 2 * JWIDU 
                NWD2  = IWD2U - IWD2L + 1
                IF ( NWD2 .LT. 3 ) GOTO 100
              END IF
              IF ( ISHP .EQ. SHPVOI ) THEN
                IF ( MIXFLG .AND. IDGPTH .EQ. IDXCO2 ) THEN
                  CALL MIXSHP ( NWD2, WNOWD2(IWD2L), ABSORP(IWD2L) )
                ELSE
                  CALL VOISHP ( NWD2, WNOWD2(IWD2L), ABSORP(IWD2L) )
                END IF
              ELSE IF ( ISHP .EQ. SHPLOR ) THEN
                CALL LORSHP ( NWD2, WNOWD2(IWD2L), ABSORP(IWD2L) )
              ELSE IF ( ISHP .EQ. SHPDOP ) THEN
                CALL DOPSHP ( NWD2, WNOWD2(IWD2L), ABSORP(IWD2L) )
              ELSE IF ( ISHP .EQ. SHPCHI ) THEN
                CALL CHISHP ( IPTH, NWD2,WNOWD2(IWD2L),ABSORP(IWD2L) )
              ELSE IF ( ISHP .EQ. SHPVVW ) THEN
                CALL VVWSHP ( NWD2, WNOWD2(IWD2L), ABSORP(IWD2L) )
              ELSE
                STOP 'F-RFMWID: Unrecognised lineshape value'
              ENDIF
              IF ( VVWFLG .AND. ISHP .NE. SHPVVW ) 
     &          CALL VVWCOR ( NWD2, WNOWD2(IWD2L), ABSORP(IWD2L) )
C If applying CKD H2O continuum need to subtract abs.coeff at 25cm
              IF ( SUBWNG ) 
     &          CALL SUB25W ( NWD2, ABSORP(IWD2L) )
C
C Set absorption=0.0 beyond 25cm-1 from line centre 
              IF ( WNOWD2(IWD2L) .LT. WNUM-FWIND ) 
     &          ABSORP(IWD2L) = 0.0
              IF ( WNOWD2(IWD2L+1) .LT. WNUM-FWIND ) 
     &          ABSORP(IWD2L+1) = 0.0
              IF ( WNOWD2(IWD2U-1) .GT. WNUM+FWIND ) 
     &          ABSORP(IWD2U-1) = 0.0
              IF ( WNOWD2(IWD2U) .GT. WNUM+FWIND ) 
     &          ABSORP(IWD2U) = 0.0
C
              DO JWID = JWIDL, JEXCL 
                DO IQAD = 1, 3
                  IWD2 = 2 * JWID + IQAD - 3
                  ABSWID(IQAD,JWID,IPTH) = ABSWID(IQAD,JWID,IPTH) + 
     &                        SNGL ( ANTADJ * DBLE ( ABSORP(IWD2) ) )
                  CNTWID(IQAD,JWID,IPTH) = CNTWID(IQAD,JWID,IPTH) + 
     &                        SNGL ( CNTADJ * DBLE ( ABSORP(IWD2) ) )
                END DO
                IF ( .NOT. COUNTD(JWID) ) THEN
                  NLNWID(JWID,IGAS) = NLNWID(JWID,IGAS) + 1
                  COUNTD(JWID) = .TRUE.
                END IF
              END DO
              DO JWID = JEXCU, JWIDU
                DO IQAD = 1, 3
                  IWD2 = 2 * JWID + IQAD - 3
                  ABSWID(IQAD,JWID,IPTH) = ABSWID(IQAD,JWID,IPTH) + 
     &                        SNGL ( ANTADJ * DBLE ( ABSORP(IWD2) ) )
                  CNTWID(IQAD,JWID,IPTH) = CNTWID(IQAD,JWID,IPTH) + 
     &                        SNGL ( CNTADJ * DBLE ( ABSORP(IWD2) ) )
                END DO
                IF ( .NOT. COUNTD(JWID) ) THEN
                  NLNWID(JWID,IGAS) = NLNWID(JWID,IGAS) + 1
                  COUNTD(JWID) = .TRUE.
                END IF
              END DO
              IF ( .NOT. COUNTL ) THEN ! Update line ct stats for spec.rng.
                IF ( WNUM .LT. WNLOUT .OR. WNUM .GT. WNUOUT ) THEN ! outside
                  NOLWID(IGAS) = NOLWID(IGAS) + 1
                ELSE                                               ! inside
                   NILWID(IGAS) = NILWID(IGAS) + 1                    
                END IF
                COUNTL = .TRUE.
              END IF
            END IF      ! End case of line strong enough for inclusion
          END IF        ! End case of a correct line for path
        END DO          ! End loop over paths
        GOTO 100        ! End loop over lines
      END IF            ! Carry on from here if EOF or last WNUM read > reqd.
C
      END       
