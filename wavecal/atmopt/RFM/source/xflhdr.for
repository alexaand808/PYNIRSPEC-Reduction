      SUBROUTINE XFLHDR ( LUNXFL, IDGAS, WN1SPC, WN2SPC, 
     &                    NXST, NXSP, NREC, HEADER, FAIL, ERRMSG )
C
C VERSION
C     29-AUG-13  AD  Original. Extracted from parts of FILXFL.
C
C DESCRIPTION
C     Extract header info for next spectral range from .xsc files.
C
      IMPLICIT NONE 
C
C ARGUMENTS
      INTEGER          LUNXFL !  I  LUN for current X/S file
      INTEGER          IDGAS  !  O  ID of gas
      DOUBLE PRECISION WN1SPC !  O  Min. wavenumber [cm-1] for spectral range 
      DOUBLE PRECISION WN2SPC !  O  Max. wavenumber [cm-1] for spectral range 
      INTEGER          NXST   !  O  No. of (p,T) tables X/S data set 
      INTEGER          NXSP   !  O  Total No. data points in spectral range
      INTEGER          NREC   !  O  No. records read from file
      CHARACTER*(*)    HEADER !  O  RFM header for next spectral range
      LOGICAL          FAIL   !  O  T=A fatal error was detected
      CHARACTER*80     ERRMSG !  O  Error message written if FAIL is TRUE
C      
C LOCAL VARIABLES      
      INTEGER          IOS      ! Saved value of IOSTAT for error message
      INTEGER          IREC     ! Record counter
      INTEGER          IXST     ! Counter for (p,T) tables
      INTEGER          IXSW     ! Counter for data points
      INTEGER          NDREC    ! No. of data records in (p,T) table
      INTEGER          NXSW     ! No. of data points in each (p,T) table
      LOGICAL          GOTRNG   ! T=found spectral range in RFM header
      LOGICAL          H2KFMT   ! T=HITRAN 2000 format, F=Earlier format
      REAL             RDUMMY   ! Dummy variable for reading data values
      DOUBLE PRECISION WN1XST   ! Min. wavenumber [cm-1] for (pT) table
      DOUBLE PRECISION WN2XST   ! Max. wavenumber [cm-1] for (pT) table
      CHARACTER*10     FMT      ! File Format identifier
      CHARACTER*100    RECORD   ! Record read from X/S data file
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      NREC = 0            ! value returned if End-of-file is reached
C
      READ ( LUNXFL, '(A)', IOSTAT=IOS, ERR=900, END=100 ) HEADER
      NREC = 1
      IF ( HEADER(1:3) .NE. '***' ) THEN
        FAIL = .TRUE.
        ERRMSG = 'F-XFLHDR: Unexpected record in file, not RFM header'
        RETURN
      END IF
C
C First see if wavenumber range is included in header (true for more recent
C files created by hitxsc but not older data )
      WN1SPC = 0.0
      WN2SPC = 0.0 
      READ ( HEADER, '(3X,I7,I10,A10,2F10.4)', IOSTAT=IOS ) 
     &  IDGAS, NXST, FMT, WN1SPC, WN2SPC
      GOTRNG = IOS .EQ. 0 .AND. ( WN2SPC .NE. 0.0 )  ! got wno range from header
      IF ( .NOT. GOTRNG ) THEN
        READ ( HEADER, '(3X,I7,I10,A10)', IOSTAT=IOS, ERR=900 ) 
     &    IDGAS, NXST, FMT
        WN1SPC = 1.0D6
        WN2SPC = 0.0
      END IF
C
      H2KFMT = FMT .NE. ' '
C
      NXSP = 0
      DO IXST = 1, NXST
        READ ( LUNXFL, '(A)', IOSTAT=IOS, ERR=900) RECORD
        NREC = NREC + 1
        IF ( H2KFMT ) THEN
          READ ( RECORD(21:60), *, IOSTAT=IOS, ERR=900 ) 
     &      WN1XST, WN2XST, NXSW
        ELSE
          READ ( RECORD(11:60), *, IOSTAT=IOS, ERR=900 ) 
     &      WN1XST, WN2XST, NXSW
        END IF
        NXSP = NXSP + NXSW
C If spectral range not in header, construct from individual pT table ranges
        IF ( .NOT. GOTRNG ) THEN
          WN1SPC = MIN ( WN1SPC, WN1XST )
          WN2SPC = MAX ( WN2SPC, WN2XST )
        END IF
C
C No. of data records of assuming 10 data values per record
        NDREC = ( NXSW - 1 ) / 10 + 1    ! integer division
        DO IREC = 1, NDREC-1
          READ ( LUNXFL, '(10F10.3)' ) ( RDUMMY, IXSW = 1, 10 )
        END DO
        READ ( LUNXFL, '(10F10.3)' ) 
     &    ( RDUMMY, IXSW = ( NDREC - 1 ) * 10 + 1, NXSW ) 
        NREC = NREC + NDREC

      END DO
 100  CONTINUE         ! jump here if end-of-file is reached
      FAIL = .FALSE.
      RETURN
C
  900 FAIL = .TRUE.
      WRITE ( ERRMSG, '(A,I11)' ) 
     &  'F-XFLHDR: I/O failure on X/S file. IOSTAT=', IOS
C
      END
