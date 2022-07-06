      SUBROUTINE RFMSPC ( ISPC, FAIL, ERRMSG )
C
C VERSION
C     17-AUG-06  AD  Calculate WNDFIN
C     24-JAN-03  AD  Calculate IOFOUT
C     02-JAN-03  AD  Move SPCGRD,SPCLUT to end of subroutine
C     13-APR-00  AD  Initialise line count stats NLNWID, NILWID, NOLWID
C     09-AUG-99  AD  Change WNOWID from Real to D.P.
C     11-DEC-97  AD  Add SPCGRD
C     03-JUL-97  AD  Add SPCLUT, FLGCOM.INC, arguments FAIL, ERRMSG 
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Initialise RFM for each new spectral range.
C     Called by RFM once for each Spectral range.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      ISPC    !  I  Current spectral range
      LOGICAL      FAIL    !  O  T=a fatal error has been detected
      CHARACTER*80 ERRMSG  !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  SPCFIN ! Set Fine Mesh grid parameters from spectral range/resln.
     &, SPCGRD ! Initialise GRD data for each new spectral range.
     &, SPCLUT ! Initialise LUT data for each new spectral range.
     &, SPCWID ! Set Wide Mesh grid parameters from spectral range/resln.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'rfmcon.inc' ! Constants for RFM calculations
C
C COMMON VARIABLES
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'outcom.inc' ! RFM output file data
      INCLUDE 'pthcom.inc' ! RFM path data
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
      INCLUDE 'widcom.inc' ! Wide mesh data
C
C LOCAL VARIABLES
      INTEGER IGAS         ! Absorber counter
      INTEGER IPTH         ! Path counter (1:CLC) (Calc.paths only)
      INTEGER IQAD         ! Parabolic point counter (1:3)
      INTEGER IWID         ! Widemesh boundary counter (0:NWID)
      INTEGER JWID         ! Widemesh interval counter (1:NWID)
      INTEGER MFIN         ! No.fine grid points per output grid point
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      CALL SPCWID ( ISPC )  ! Sets NWID, DELWID, WNLWID, WNUWID, WN1WID, WN2WID
      CALL SPCFIN ( ISPC )  ! Sets NFIN, WNRFIN
C
C Set Wavenumbers at widemesh grid points and half-grid points
C
      DO IWID = 0, NWID
        WNOWID(IWID) = WN1WID + IWID * DELWID
        WNOWD2(2*IWID) = WNOWID(IWID)             ! WNOWD2(0:2*MAXWID+1)
        WNOWD2(2*IWID+1) = WNOWID(IWID) + 0.5D0 * DELWID 
      END DO
C
C Initialise parabolic absorption values
C
      DO IPTH = 1, NCLC 
        DO JWID = 1, NWID
          DO IQAD = 1, 3
            ABSWID(IQAD,JWID,IPTH) = 0.0
            CNTWID(IQAD,JWID,IPTH) = 0.0
          END DO
        END DO
      END DO
C
C Initialise line count stats
      DO IGAS = 1, NGAS
        DO JWID = 1, NWID
          NLNWID(JWID,IGAS) = 0
        END DO
        NILWID(IGAS) = 0
        NOLWID(IGAS) = 0
      END DO
C
C Save upper and lower wavenumbers required for output spectra
C
      WNLOUT = WNLSPC(ISPC)
      WNUOUT = WNUSPC(ISPC)
      IF ( WNRSPC(ISPC) .GT. 1.0D0 ) THEN
        WNROUT = 1.0D0 / WNRSPC(ISPC)
      ELSE
        WNROUT = WNRSPC(ISPC)
      END IF
      NPWOUT = NINT ( DELWID / WNROUT )
      MFIN = NINT ( WNROUT / WNRFIN ) 
C Shift fine grid by fraction of WNRFIN if required so that it matches o/p grid
      WNDFIN = MOD ( WNLOUT - WNOWID(0), WNRFIN )
      IOFOUT = NINT ( ( WNLOUT - WNOWID(0) - WNDFIN ) / WNRFIN )
      IOFOUT = MOD ( IOFOUT, MFIN )
      LABOUT = LABSPC(ISPC)
C
C Calculate total number of points required for output spectra
C
      NPTOUT = NINT ( ( WNUOUT - WNLOUT ) / WNROUT ) + 1
C
      IF ( GRDFLG ) THEN
        CALL SPCGRD ( ISPC, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
C
      IF ( LUTFLG ) THEN
        CALL SPCLUT ( ISPC, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
C
      END
