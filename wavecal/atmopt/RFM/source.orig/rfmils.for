      SUBROUTINE RFMILS ( IWID, ISPC, FAIL, ERRMSG )
C
C VERSION
C     17-AUG-06  AD  Correction: set IOPOFF allowing for offset grids
C     01-JUN-04  AD  Correction: Check IG1BUF LE NGRD (for IWID=NWID+1)
C     09-AUG-03  AD  Overwrite WNOFIN with output grid
C     24-JAN-03  AD  Use IOFOUT instead of calculating IDXOFF
C     13-JAN-03  AD  Allow for OUTGRD = TRUE
C     02-JAN-03  AD  Change logic to try and speed up
C     04-DEC-01  AD  Correction: Save IDXBUF as well
C     26-MAY-00  AD  Convert RADTAN, RADBUF from S.P. to D.P.
C     27-APR-00  AD  Convert ILS to DP
C     07-OCT-99  AD  Remove MAXILF as a local constant, now in RFMSIZ.INC
C     09-AUG-99  AD  Set TRATAN = 0.0D0 rather than = 0.0
C     02-JAN-99  AD  Use LTAN rather than NTAN to allow for Jacobians.
C     09-NOV-98  AD  Bug fix: add IDXOFF to allow for offset o/p grids.
C     17-AUG-98  AD  Remove general SAVE statement
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Convolve RFM spectra with ILS.
C     Called by RFM for each spectral range if ILS flag enabled.
C
      IMPLICIT NONE
      SAVE
C
C ARGUMENTS
      INTEGER      ISPC    !  I  Index of current Spectral range 
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error occurs
      CHARACTER*80 ERRMSG  !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  ILSINT ! Interpolate ILS function to fine-grid spacing.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'grdcom.inc' ! Irregular grid
      INCLUDE 'outcom.inc' ! RFM output file data
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL CONSTANTS
      INTEGER MAXIL2               ! Used to initialise ILS DATA statement
        PARAMETER ( MAXIL2 = 2 * MAXILF + 1 )
      INTEGER NBUF                 ! No wide-mesh buffers for storing FM data
        PARAMETER ( NBUF = 3 )     ! 3=central mesh plus one either side
C
C LOCAL VARIABLES
      INTEGER I                 ! Index of point in ILS function
      INTEGER IGRD              ! Irregular grid point counter
      INTEGER I1,I2             ! Indices marking useful range of F(I1:I2)
      INTEGER IMIN,IMAX         ! range of I1:I2 within each buffer
      INTEGER IBUF              ! Buffer index 
      INTEGER IDXBUF(NBUF)      ! Fine Grid Index of first point in each buffer
      INTEGER IDXFIN            ! Fine Grid Index of each fine mesh point
      INTEGER IDXMAX            ! Max fine grid index spanned by ILS
      INTEGER IDXMIN            ! Min fine grid index spanned by ILS
      INTEGER IDXOUT            ! Fine Grid Index of Output Point 
      INTEGER IFIN              ! Fine grid point counter
      INTEGER IG1BUF(NBUF)      ! Index of first irreg.grid point in buffer
      INTEGER IG2BUF(NBUF)      ! Index of last  irreg.grid point in buffer
      INTEGER IOPG              ! Counter for o/p pts affected by i.g.point
      INTEGER IOPOFF            ! Output grid point offset for current intvl.
      INTEGER IOUT              ! Output grid point counter
      INTEGER ITAN              ! Tangent height counter
      INTEGER IWID              ! Widemesh interval counter
      INTEGER JBUF              ! Buffer counter
      INTEGER LBUF              ! Index of previous widemesh buffer
      INTEGER MFIN              ! No.fine grid points per output grid point
      DOUBLE PRECISION ILS(-MAXILF:MAXILF) ! ILS interp.to Fine grid spacing
      DOUBLE PRECISION RADBUF(MAXFIN,MAXTAN,NBUF) ! Radiance on fine grid
      DOUBLE PRECISION RADSUM   ! Temporary store for radiance summation
      DOUBLE PRECISION TRABUF(MAXFIN,MAXTAN,NBUF) ! Transmission on fine grid
      DOUBLE PRECISION TRASUM   ! Temporary store for transmittance summation
      DOUBLE PRECISION WNOOFF   ! Wno of first point on o/p grid 
C
C DATA STATEMENTS
        DATA IBUF / 0 /             ! To satisfy FLINT
        DATA MFIN / 0 /             ! To satisfy FLINT
        DATA ILS / MAXIL2 * 0.0D0 / ! To satisfy FLINT
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C For first wide-mesh interval in each spectral range, just initialise and
C copy into buffer - no ILS convolution yet.
C
C If GRDFLG is TRUE, NFIN is the number of irregular grid points,
C If GRDFLG is FALSE, NFIN is the number of fine grid points 
      IF ( IWID .EQ. 1 ) THEN
        IBUF = 1
        IDXBUF(IBUF) = 1
        IG1BUF(IBUF) = 1
        IG2BUF(IBUF) = NFIN
        IOPOFF = - ( ( WNLOUT - WN1FIN + 0.1*WNRFIN ) / WNROUT ) ! integer div.
        WNOOFF = WN1FIN + IOFOUT * WNRFIN 
        DO JBUF = 2, NBUF
          IDXBUF(JBUF) = 0
          IG1BUF(JBUF) = 0     ! not necessary to initialise IG2BUF
        END DO
        IF ( .NOT. ( GRDFLG .AND. OUTGRD ) ) THEN  ! using full fine grid
          MFIN = NFIN / NPWOUT                     ! Should be exact multiple
          CALL ILSINT ( ISPC, MAXILF, ILS, I1, I2, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
          IF ( I1 .LT. -NFIN .OR. I2 .GT. NFIN ) THEN  
            FAIL = .TRUE.
            ERRMSG = 'F-RFMILS: ILS width > assumed max of 1 widemesh'
            RETURN
          END IF
        END IF
        DO ITAN = 1, LTAN
          DO IFIN = 1, NFIN
            RADBUF(IFIN,ITAN,IBUF) = RADTAN(IFIN,ITAN)
            TRABUF(IFIN,ITAN,IBUF) = TRATAN(IFIN,ITAN)
          END DO
        END DO
        RETURN
      END IF
C
C For subsequent widemesh intervals, copy current calculations into buffer
C
      LBUF = IBUF
      IBUF = IBUF + 1
      IF ( IBUF .GT. NBUF ) IBUF = 1
      IDXBUF(IBUF) = IDXBUF(LBUF) + NFIN 
      DO ITAN = 1, LTAN
        DO IFIN = 1, NFIN
          RADBUF(IFIN,ITAN,IBUF) = RADTAN(IFIN,ITAN)
          TRABUF(IFIN,ITAN,IBUF) = TRATAN(IFIN,ITAN)
        END DO
      END DO
C
C Overwrite WNOFIN with wavenumber grid for output interval
C WN1FIN is current fine mesh, output interval lags by 1,
C IOFOUT is offset of output grid relative to WN1FIN in finemesh intervals
      DO IOUT = 1, NPWOUT
        WNOFIN(IOUT) = WNOOFF + ( IOUT - 1 ) * WNROUT
      END DO
      WNOOFF = WN1FIN + IOFOUT * WNRFIN    ! Set for next output interval
C
C Perform direct interpolation of fine grid to output grid
      IF ( GRDFLG .AND. OUTGRD ) THEN
C Set pointers for current buffer IBUF 
        IG1BUF(IBUF) = IG2BUF(LBUF) + 1
        IG2BUF(IBUF) = IG2BUF(LBUF) + NFIN
C Initialise output grid (centred on previous buffer)
        DO IOUT = 1, NPWOUT     
          DO ITAN = 1, LTAN
            RADTAN(IOUT,ITAN) = 0.0D0
            TRATAN(IOUT,ITAN) = 0.0D0
          END DO
        END DO
C Loop over buffers storing rad and tra calc. at irreg grid pts
C (note that for IWID=NWID+1 NFIN will be still be at previous value)
        DO JBUF = 1, NBUF            ! check if buffer loaded 
          IF ( IG1BUF(JBUF) .GT. 0 .AND. IG1BUF(JBUF) .LE. NGRD ) THEN  
            IFIN = 0                 
            DO IGRD = IG1BUF(JBUF), IG2BUF(JBUF) ! Loop over i.g. within buffer
              IFIN = IFIN + 1        ! Counter over i.g. points within buffer
              DO IOPG = 1, NOPGRD(IGRD)     ! Affected o/p points for each ig
                IOUT = IOPGRD(IOPG,IGRD) - IOPOFF  ! O/p posn in widemesh
                IF ( IOUT .GE. 1 .AND. IOUT .LE. NPWOUT ) THEN ! Within intvl
                  DO ITAN = 1, LTAN
                    RADTAN(IOUT,ITAN) = RADTAN(IOUT,ITAN) +
     &                RADBUF(IFIN,ITAN,JBUF) * WGTGRD(IOPG,IGRD)
                    TRATAN(IOUT,ITAN) = TRATAN(IOUT,ITAN) +
     &                TRABUF(IFIN,ITAN,JBUF) * WGTGRD(IOPG,IGRD)
                  END DO
                END IF
              END DO
            END DO
          END IF
        END DO
        IOPOFF = IOPOFF + NPWOUT
        FAIL = .FALSE.
        RETURN
      END IF
C
C Perform ILS convolution for previous widemesh interval starting at IDXOUT
      IDXOUT = 1 + ( IWID - 2 ) * NFIN + IOFOUT
      DO IOUT = 1, NPWOUT     ! Output grid = MFIN fine mesh intvls.
        DO ITAN = 1, LTAN
          RADTAN(IOUT,ITAN) = 0.0D0
          TRATAN(IOUT,ITAN) = 0.0D0
        END DO
        IDXMIN = IDXOUT + I1  ! Establish fine grid index covered by ILS
        IDXMAX = IDXOUT + I2
        DO JBUF = 1, NBUF
          IF ( IDXBUF(JBUF) .NE. 0 .AND.           ! If buffer written
     &         IDXBUF(JBUF) .LE. IDXMAX .AND.      ! and overlaps ILS
     &         IDXBUF(JBUF) + NFIN .GT. IDXMIN ) THEN  ! NB buffer is 1:NFIN-1
C Restrict I1,I2 to buffer range IDXBUF(1:NFIN-1)
            IMIN = I1 + MAX ( 0, IDXBUF(JBUF) - IDXMIN)
            IMAX = I2 + MIN ( 0, IDXBUF(JBUF)+NFIN-1 - IDXMAX )
            IDXFIN = IDXOUT + IMIN
            DO ITAN = 1, LTAN
              RADSUM = 0.0D0
              TRASUM = 0.0D0
              IFIN = IDXFIN - IDXBUF(JBUF) 
              DO I = IMIN, IMAX
                IFIN = IFIN + 1 
                RADSUM = RADSUM + RADBUF(IFIN,ITAN,JBUF) * ILS(I)
                TRASUM = TRASUM + TRABUF(IFIN,ITAN,JBUF) * ILS(I)
              END DO
              RADTAN(IOUT,ITAN) = RADTAN(IOUT,ITAN) + RADSUM
              TRATAN(IOUT,ITAN) = TRATAN(IOUT,ITAN) + TRASUM
            END DO
          END IF
        END DO
        IDXOUT = IDXOUT + MFIN
      END DO
C
      FAIL = .FALSE.
      END
