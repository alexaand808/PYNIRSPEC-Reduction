      SUBROUTINE LUTGRD ( LUNLUT, BINFIL, MAXGRD, 
     &                    NGRD, IDXGRD, FAIL, ERRMSG )
C
C VERSION
C     14-DEC-01  AD  Ensure ICHAR and ITWO initialised.
C     08-JUN-00  AD  Original.
C
C DESCRIPTION
C     Read coded irregular grid data and turn into integer indices.
C     Called by LUTFIL for each LUT/TAB file
C     NB: if NGRD>MAXGRD, only indices 1:MAXGRD will be written to IDXGRD. 
C     Also permissible to have NFIN and NGRD as same external variable.
C     Leaves file pointer at start of next record after coded grid.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      LUNLUT         !  I  LUN for GRD file
      LOGICAL      BINFIL         !  I  T=Binary file, F=ASCII
      INTEGER      MAXGRD         !  I  External array size of IDXGRD 
      INTEGER      NGRD           ! I/O I=No.Full grid, O=No.Irreg grid points
      INTEGER      IDXGRD(MAXGRD) !  O  Array of irreg.grid indices
      LOGICAL      FAIL           !  O  T=Fatal error, F=OK
      CHARACTER*80 ERRMSG         !  O  Error message written if FAIL is TRUE
C
C LOCAL VARIABLES
      INTEGER      ICHAR !  Character# in coded irreg.grid data record
      INTEGER      IFIN  !  Regular fine grid index
      INTEGER      IOS   !  Status of I/O Error 
      INTEGER      ITWO  !  Integer corresponding to digit within character.
      INTEGER      IVAL  !  Integer value of hexadecimal character
      INTEGER      NFIN  !  Local version of NFUL
      CHARACTER*50 REC50 !  Record containing coded irregular grid data
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      NFIN = NGRD
      NGRD = 0
      ICHAR = 0          ! Not necessary, but avoids "unitialized" warnings
      ITWO = 0           ! Not necessary, but avoids "unitialized" warnings
C
      DO IFIN = 1, NFIN
        IF ( MOD ( IFIN, 200 ) .EQ. 1 ) THEN
          IF ( BINFIL ) THEN
            READ ( LUNLUT, ERR=100, IOSTAT=IOS ) REC50
          ELSE
            READ ( LUNLUT, '(A)', ERR=100, IOSTAT=IOS ) REC50
          END IF
          ICHAR = 0
        END IF
        IF ( MOD ( IFIN, 4 ) .EQ. 1 ) THEN
          ICHAR = ICHAR + 1
          READ ( REC50(ICHAR:ICHAR), '(Z1)', ERR=100, IOSTAT=IOS ) IVAL
          ITWO = 16
        END IF
        ITWO = ITWO/2
        IF ( IVAL .GE. ITWO ) THEN
          NGRD = NGRD + 1
          IF ( NGRD .LE. MAXGRD ) IDXGRD(NGRD) = IFIN
          IVAL = IVAL - ITWO
        END IF
      END DO
C
      FAIL = .FALSE.
      RETURN
C
 100  CONTINUE
      FAIL = .TRUE.
      WRITE ( ERRMSG, '(A,I11)' ) 
     &  'F-LUTGRD: I/O failure reading coded grid data. IOSTAT=', IOS
C
      END
