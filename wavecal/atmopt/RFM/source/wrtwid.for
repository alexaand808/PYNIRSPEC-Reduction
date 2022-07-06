      SUBROUTINE WRTWID ( IWID, FAIL, ERRMSG )
C
C VERSION
C     06-APR-06  AD  Remove leading 'X' in output format
C     27-APR-00  AD  Initialise NTOEXT, NTOINT
C     13-APR-00  AD  Add total internal/external line counts   
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Write RFM widemesh line count diagnostics.
C     Called by RFM for each widemesh interval if WID flag selected.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      IWID    !  I  Index of current wide-mesh interval
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error occurs
      CHARACTER*80 ERRMSG  !  O  Error message written if FAIL is TRUE
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
      INCLUDE 'rfmcon.inc' ! Constants for RFM calculations
C
C COMMON VARIABLES
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'outcom.inc' ! RFM output file data
      INCLUDE 'cyccom.inc' ! Cyclic line data buffers
      INCLUDE 'fincom.inc' ! Fine mesh data
      INCLUDE 'widcom.inc' ! Wide mesh data
C
C LOCAL VARIABLES
      INTEGER       ICYC   ! Pointer for cyclic buffer
      INTEGER       IGAS   ! Counter for absorbers
      INTEGER       ILIN   ! Line counter for cyclic buffer
      INTEGER       IOS    ! Value of IOSTAT saved for error messages.
      INTEGER       LUN    ! Local value of Logical Unit Number.
      INTEGER       NLNFIN(MAXGAS) ! No.lines used for fine mesh calculations
      INTEGER       NLNITV(MAXGAS) ! No.lines within wide mesh interval
      INTEGER       NTOCTM ! Total number of continuum absorbers in WM interval
      INTEGER       NTOEXT ! Total number of external lines
      INTEGER       NTOFIN ! Total number of lines used for fine mesh calc.
      INTEGER       NTOINT ! Total number of internal lines
      INTEGER       NTOITV ! Total number of lines within wide mesh interval
      INTEGER       NTOWID ! Total number of lines used for wide mesh calc.
      INTEGER       NTOXSC ! Total number of X/Sect absorbers in WM interval
      DOUBLE PRECISION WNO ! Wavenumber [/cm]
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      DO IGAS = 1, NGAS
        NLNFIN(IGAS) = 0
        NLNITV(IGAS) = 0
      END DO
C
      DO ILIN = 1, NLIN
        ICYC = MOD ( ILIN + ICYC1 - 2, MAXCYC ) + 1
        WNO = WNOCYC(ICYC)
        IF ( WNO .GE. WNLFIN .AND. WNO .LE. WNUFIN ) THEN
          IGAS = IGSMOL(IDGCYC(ICYC))
          NLNFIN(IGAS) = NLNFIN(IGAS) + 1
          IF ( WNO .GE. WN1FIN .AND. WNO .LT. WN2FIN ) THEN
            NLNITV(IGAS) = NLNITV(IGAS) + 1
          END IF
        END IF
      END DO
      NTOWID = 0
      NTOCTM = 0
      NTOXSC = 0
      NTOITV = 0
      NTOFIN = 0
      NTOEXT = 0
      NTOINT = 0
      IF ( NGAS .EQ. 1 ) THEN  ! If one gas, there is no total lines file
        LUN = LUNWID - 1
      ELSE                     ! If >1 gas, LUNWID is total lines file
        LUN = LUNWID
      END IF
      DO IGAS = 1, NGAS
        LUN = LUN + 1
        WRITE ( LUN, '(I5,F12.5,6I10)', IOSTAT=IOS, ERR=900 ) IWID,
     &    WNOWID(IWID-1), NLNWID(IWID,IGAS)+NLNFIN(IGAS), 
     &    NLNWID(IWID,IGAS), NLNFIN(IGAS),  NLNITV(IGAS),  IXSFIN(IGAS),
     &    ICTWID(IWID,IGAS)
        NTOWID = NTOWID + NLNWID(IWID,IGAS)
        NTOCTM = NTOCTM + ICTWID(IWID,IGAS)
        NTOXSC = NTOXSC + IXSFIN(IGAS)
        NTOFIN = NTOFIN + NLNFIN(IGAS)
        NTOITV = NTOITV + NLNITV(IGAS)
        IF ( IWID .EQ. NWID ) THEN 
          WRITE ( LUN, '(2I11,A)' ) NILWID(IGAS), NOLWID(IGAS),
     &      ' = Internal, external lines used for whole range'
          NTOINT = NTOINT + NILWID(IGAS)
          NTOEXT = NTOEXT + NOLWID(IGAS)
        END IF          
      END DO
C
      IF ( NGAS .GT. 1 ) THEN
        LUN = LUNWID 
        WRITE ( LUN, '(I5,F12.5,6I10)', IOSTAT=IOS, ERR=900 ) IWID,
     &    WNOWID(IWID-1), NTOWID+NTOFIN, NTOWID, NTOFIN, NTOITV,  
     &    NTOXSC, NTOCTM
        IF ( IWID .EQ. NWID ) THEN 
          WRITE ( LUN, '(2I11,A)' ) NTOINT, NTOEXT, 
     &      ' = Internal, external lines used for whole range'
        END IF          
      END IF 
C
      FAIL = .FALSE.
      RETURN
C
  900 FAIL = .TRUE.
      WRITE ( ERRMSG, '(A,I11)' ) 
     &  'F-WRTWID: I/O failure on output file. IOSTAT=', IOS
      END
