      SUBROUTINE GRDILS ( ISPC, FAIL, ERRMSG )
C
C VERSION
C     07-AUG-13  AD  Remove redundant local variable WNOOFF
C     17-AUG-06  AD  Add WNDFIN adjustment for additional grid points
C     17-JAN-03  AD  Original.
C
C DESCRIPTION
C     Calc. Irreg. grid direct contributions to convolved output grid.
C     Called by SPCGRD for each spec.range if GRD and (ILS or AVG) flags set.
C     If there is insufficient storage space (MAXOPG) OUTGRD will be set
C     false and RFM then interpolates irreg.grid to full grid followed by ILS
C     convolution (and a warning message sent to the runlog file)
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      ISPC    !  I  Index of current Spectral range
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error occurs
      CHARACTER*80 ERRMSG  !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  DBRAKT ! Return lower index I of the two points in array ARRAY(N)
     &, ILSINT ! Interpolate ILS function to fine-grid spacing
     &, RFMLOG ! Write text message to RFM log file.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C 
C COMMON VARIABLES
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'grdcom.inc' ! Irregular grid 
      INCLUDE 'outcom.inc' ! RFM output file data
      INCLUDE 'widcom.inc' ! Wide mesh data
C
C LOCAL CONSTANTS
      INTEGER MAXIL2       ! Used to initialise ILS DATA statement
        PARAMETER ( MAXIL2 = 2 * MAXILF + 1 )
      DOUBLE PRECISION WNOTOL     ! Tolerance as fraction of fine grid spacing
        PARAMETER ( WNOTOL = 0.01D0 )
C
C LOCAL VARIABLES
      INTEGER I            ! Counter for ILS fine grid intervals
      INTEGER I1, I2       ! Indices marking useful range of ILS(I1:I2)
      INTEGER I1OFF,I2OFF  ! Offsets to confine I1,I2 to within irreg.grid
      INTEGER IEND         ! 1=lower,2=upper boundary of widemesh interval
      INTEGER IGRD         ! Counter for irregular grid points
      INTEGER IOUT         ! Counter for output grid points
      INTEGER IWID         ! Widemesh interval counter
      INTEGER LGRD         ! Store original value of NGRD
      INTEGER MGRD         ! Index for inserting new irreg.grid point
      INTEGER NOP1         ! No.o/p points for current lower irr.grd.pt
      INTEGER NOP2         ! No.o/p points for current upper irr.grd.pt
      CHARACTER*80 MESSGE  ! Text message for runlog
      DOUBLE PRECISION ILS(-MAXILF:MAXILF) ! ILS interp.to Fine grid spacing
      DOUBLE PRECISION WEIGHT ! Weight for current lower irr.grd.pt
      DOUBLE PRECISION WN1GRD ! Lower boundary of current irreg.grid interval
      DOUBLE PRECISION WN2GRD ! Upper boundary of current irreg.grid interval
      DOUBLE PRECISION WNODEL ! Tolerance in wavenumbers
      DOUBLE PRECISION WNOEND ! Wavenumber at either boundary of widemesh intvl
      DOUBLE PRECISION WNOILS ! Wavenumber for ILS convolution
      DOUBLE PRECISION WNOOUT ! Wavenumber of each output grid point
C
C DATA STATEMENTS
      DATA ILS / MAXIL2 * 0.0D0 /     ! Not actually necessary
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Currently only set up for linearly interpolated irregular grids
      IF ( FNCGRD .NE. 'lin' ) THEN
        OUTGRD = .FALSE.
        MESSGE = 'W-GRDILS: direct irr.grd to o/p '//
     &           'only available for linear interp.'
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        RETURN
      END IF
C
      WNODEL = WNOTOL * WNRFIN
C
C To simulate exactly the results from interpolation to fine grid followed by
C convolution is necessary to insert additional grid points at the widemesh
C boundaries (normally done in RFMGRD but bypassed if OUTGRD is set TRUE).
C This isn't necessary with direct interpolation, but just added to provide
C a direct comparison
C Comment this next section out if boundary grid points are NOT required
c      goto 100
      LGRD = NGRD
      DO IWID = 1, NWID
        DO IEND = 1, 2
          IF ( IEND .EQ. 1 ) THEN            ! lower wno boundary
            WNOEND = WNOWID(IWID-1) + WNDFIN
          ELSE
            WNOEND = WNOWID(IWID) - WNRFIN + WNDFIN  ! upper wno boundary
          END IF
          CALL DBRAKT ( WNOGRD, NGRD, WNOEND+WNODEL, MGRD ) 
          IF ( MGRD .EQ. 0 .OR. 
     &         ABS ( WNOGRD ( MAX ( MGRD,1 ) ) - WNOEND ) 
     &         .GT. WNODEL ) THEN ! add point
            IF ( NGRD .EQ. MAXGRD ) THEN
              FAIL = .TRUE.
              ERRMSG = 
     &          'F-GRDILS: MAXGRD too small to allow extra points'
              RETURN
            END IF
            NGRD = NGRD + 1
            DO IGRD = NGRD, MGRD+2, -1
              WNOGRD(IGRD) = WNOGRD(IGRD-1)
            END DO
            WNOGRD(MGRD+1) = WNOEND
          END IF
        END DO
      END DO
      WRITE ( MESSGE, '(A,I5,A)' ) 'I-GRDILS: ', NGRD-LGRD, 
     &  ' extra irreg grid points added at widemesh boundaries'
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C End of section adding boundary grid points
C
 100  continue
      DO IGRD = 1, NGRD
        NOPGRD(IGRD) = 0
        DO IOUT = 1, MAXOPG
          WGTGRD(IOUT,IGRD) = 0.0D0
        END DO
      END DO
C
      CALL ILSINT ( ISPC, MAXILF, ILS, I1, I2, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
c     do i = i1, i2
c         ils(i) = 0.0d0
c     end do
c      ils(0) = 1.0d0
C
C
      WNOOUT = WNLOUT                            ! Lower WNO for output
      DO IOUT = 1, NPTOUT                      ! Loop over all output points
C
C Check highest wavenumber required for ILS convolution overlaps irreg.grid
        WNOILS = WNOOUT + I2 * WNRFIN            ! Higher WNO for ILS
        I2OFF = 0
        DO WHILE ( WNOILS .GT. WNOGRD(NGRD) + WNODEL )       
          WNOILS = WNOILS - WNRFIN
          I2OFF = I2OFF + 1
        END DO
C
C Check lowest wavenumber required for ILS convolution overlaps irreg.grid
        WNOILS = WNOOUT + I1 * WNRFIN            ! Lower WNO for ILS
        I1OFF = 0
        DO WHILE ( WNOILS .LT. WNOGRD(1) - WNODEL )       
          WNOILS = WNOILS + WNRFIN
          I1OFF = I1OFF + 1
        END DO
C
C Find irregular grid points IGRD:IGRD+1 which span initial low end WNOILS
C Small rounding errors may cause values IGRD=NGRD
        CALL DBRAKT ( WNOGRD, NGRD, WNOILS, IGRD ) 
C Rounding errors may also cause IGRD value one below desired value
        IF ( ABS (  WNOGRD(IGRD+1) - WNOILS ) .LE. WNODEL ) 
     &    IGRD = IGRD + 1
        IF ( IGRD .EQ. 0    ) STOP 'F-GRDILS: Logical Error#1'
        IF ( IGRD .EQ. NGRD ) IGRD = NGRD - 1
C
C Enter/increment output point in list of outputs to which IGRD contributes
        WN2GRD = WNOGRD(IGRD)        ! Copied to WN1GRD on first iter of I.
        NOP2 = NOPGRD(IGRD)
        IF ( NOP2 .GT. 0 ) THEN
          IF ( IOPGRD(NOP2,IGRD) .NE. IOUT ) THEN
            IF ( NOP2 .EQ. MAXOPG ) GOTO 900
            NOP2 = NOP2 + 1
            IOPGRD(NOP2,IGRD) = IOUT
            NOPGRD(IGRD) = NOP2
          END IF
        ELSE
          IOPGRD(1,IGRD) = IOUT
          NOP2 = 1
          NOPGRD(IGRD) = NOP2
        END IF
C
C Loop over ILS fine grid within irregular grid
        IGRD = IGRD - 1       ! Force upper grid point to be loaded on Iter#1
        DO I = I1 + I1OFF, I2 - I2OFF
C
C Move IGRD:IGRD+1 to next irr.grid interval if WNOILS is above IGRD+1
C (ILS and irregular grid both on the same fine grid, so only check once)
C True on Iter#1
          IF ( WNOILS .GT. WNOGRD(IGRD+1) - WNODEL 
     &                        .AND. IGRD+1 .LT. NGRD ) THEN   
            WN1GRD = WN2GRD
            NOP1 = NOP2
            IGRD = IGRD + 1
            IF ( IGRD .GE. NGRD ) STOP 'F-GRDILS: Logical Error#3'
            WN2GRD = WNOGRD(IGRD+1)
            NOP2 = NOPGRD(IGRD+1)
            IF ( NOP2 .GT. 0 ) THEN
              IF ( IOPGRD(NOP2,IGRD+1) .NE. IOUT ) THEN
                IF ( NOP2 .EQ. MAXOPG ) GOTO 900
                NOP2 = NOP2 + 1
                IOPGRD(NOP2,IGRD+1) = IOUT
                NOPGRD(IGRD+1) = NOP2
              END IF
            ELSE
              IOPGRD(1,IGRD+1) = IOUT
              NOP2 = 1
              NOPGRD(IGRD+1) = NOP2
            END IF
          END IF
C
          WEIGHT = ( WN2GRD - WNOILS ) / ( WN2GRD - WN1GRD )
          WGTGRD(NOP1,IGRD) = WGTGRD(NOP1,IGRD) + WEIGHT * ILS(I) 
          WGTGRD(NOP2,IGRD+1) = WGTGRD(NOP2,IGRD+1) + 
     &      ( 1.0D0 - WEIGHT ) * ILS(I)
          WNOILS = WNOILS + WNRFIN
        END DO
        WNOOUT = WNOOUT + WNROUT
      END DO
C
C Normal exit
      OUTGRD = .TRUE.
      RETURN
C 
C Exit without direct mapping - go via interpolation & convolution instead.
 900  CONTINUE
      OUTGRD = .FALSE.
      MESSGE = 'W-GRDILS: MAXOPG too small for direct irr.grd to o/p'
      CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
C
      END
