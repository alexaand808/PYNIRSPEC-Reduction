      SUBROUTINE RFMPRF ( LUNPRF, FAIL, ERRMSG )
C
C VERSION
C     29-AUG-13  AD  Remove redundant local variable IDXPRF
C                    Add MOVGRA to EXTERNAL list
C     26-APR-11  AD  Bug#79: Correct handling of isotopic profiles
C     24-MAR-11  AD  Bug#78: Multiply aerosol by 1E6 compared to previous
C     20-JUN-07  AD  Bug#63: Check ISOMOL
C     06-JAN-04  AD  Original.
C
C DESCRIPTION
C     Open PRF files and write out RFM internal profile.
C     Called once by RFM if PRF flag selected.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      LUNPRF  !  I  Next available LUN
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error occurs
      CHARACTER*80 ERRMSG  !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  MAKNAM ! Construct filename for RFM output files.
     &, MOVGRA ! Exchange profiles between /ATMCOM/ and /GRACOM/ arrays.
     &, RFMLOG ! Write text message to RFM log file.
C
C GLOBAL CONSTANTS
      INCLUDE 'idxcon.inc' ! HITRAN molecular indices
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'gracom.inc' ! Atmospheric 2-D field data
      INCLUDE 'outcom.inc' ! RFM output file data
C
C LOCAL VARIABLES
      INTEGER       IATM   ! Atmospheric layer counter
      INTEGER       IEND   ! Index of last non-blank character in CODGAS
      INTEGER       IGAS   ! Absorber counter
      INTEGER       IGRA   ! Horizontal angle counter
      INTEGER       IMOL   ! HITRAN index of species
      INTEGER       IOS    ! Value of IOSTAT saved for error messages.
      INTEGER       ISO    ! Isotope counter
      INTEGER       MGRA   ! Either 1 (1D atm) or NGRA (2D atm)
      CHARACTER*132 MESSGE ! Text message sent to Log file.
      CHARACTER*80  NAMFIL ! Name of file actually opened (incl. RUNID)
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      IF ( GRAFLG ) THEN
        MGRA = NGRA
      ELSE
        MGRA = 1
      END IF
C
      DO IGRA = 1, MGRA
C
C Construct filename and open file
        IF ( GRAFLG ) THEN       ! Insert angle into output filename
          CALL MAKNAM ( 0, 0, 0, IGRA, NAMPRF, NAMFIL )
        ELSE                     ! No angle required
          CALL MAKNAM ( 0, 0, 0, 0, NAMPRF, NAMFIL )
        END IF
        MESSGE = 'I-RFMPRF: Opening output file: '//NAMFIL
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
C
        IF ( NEWFLG ) THEN
          OPEN ( UNIT=LUNPRF, FILE=NAMFIL, STATUS='NEW', 
     &           IOSTAT=IOS, ERR=900 )
        ELSE
          OPEN ( UNIT=LUNPRF, FILE=NAMFIL, STATUS='UNKNOWN', 
     &           IOSTAT=IOS, ERR=900 )
        END IF
C
C Write File Header 
        WRITE ( LUNPRF, '(A,F8.3,A)', IOSTAT=IOS, ERR=900 ) 
     &      '! RFM internal profile written by RFM v.'//VERSN
        WRITE ( LUNPRF, '(A,A)', IOSTAT=IOS, ERR=900 ) '!',TXTHDR
C        
        WRITE ( LUNPRF, '(I6,A)', IOSTAT=IOS, ERR=900 ) 
     &    NATM, ' ! No. profile levels'
        WRITE ( LUNPRF, '(A)', IOSTAT=IOS, ERR=900 ) '*HGT [km]'
        WRITE ( LUNPRF, *, IOSTAT=IOS, ERR=900 ) 
     &    ( HGTATM(IATM), IATM = 1, NATM )
C
        IF ( GRAFLG ) CALL MOVGRA ( IGRA, 0 )
C
        WRITE ( LUNPRF, '(A)', IOSTAT=IOS, ERR=900 ) '*PRE [mb]'
        WRITE ( LUNPRF, *, IOSTAT=IOS, ERR=900 ) 
     &    ( PREATM(IATM), IATM = 1, NATM )
C
        WRITE ( LUNPRF, '(A)', IOSTAT=IOS, ERR=900 ) '*TEM [K]'
        WRITE ( LUNPRF, *, IOSTAT=IOS, ERR=900 ) 
     &    ( TEMATM(IATM), IATM = 1, NATM )
C
        DO IGAS = 1, NGAS
C          write (*,*) igas, codgas(igas), nisgas(igas)
          IMOL = IDXGAS(IGAS)
          IF ( IMOL .EQ. IDXAER ) THEN
            WRITE ( LUNPRF, '(A)', IOSTAT=IOS, ERR=900 ) 
     &        '*' // CODGAS(IGAS) // ' [km-1]'
            WRITE ( LUNPRF, *, IOSTAT=IOS, ERR=900 ) 
     &        ( VMRATM(IATM,IGAS) * DNSATM(IATM)*1.0E5, IATM = 1, NATM )
          ELSE
            IEND = INDEX ( CODGAS(IGAS), ' ' ) - 1
            IF ( IEND .EQ. -1 ) IEND = LEN ( CODGAS(IGAS) )
            IF ( .NOT. ISOMOL(IMOL) ) THEN   ! No isotope information
              WRITE ( LUNPRF, '(3A)', IOSTAT=IOS, ERR=900 ) 
     &          '*', CODGAS(IGAS)(1:IEND), ' [ppmv]'
            ELSE                                   ! Particular isotope
c              write (*,*) 'igas=', igas
c              do iso = 0, nisgas(igas)
c                write (*,*) 'iso=', iso, isogas(iso,igas)
c              enddo
              IF ( ISOGAS(0,IGAS) .EQ. IGAS ) THEN
                WRITE ( LUNPRF, '(3A,I1,A)', IOSTAT=IOS, ERR=900 ) 
     &            '*', CODGAS(IGAS)(1:IEND), ' [ppmv]'
              ELSE
                DO ISO = 1, NISGAS(IGAS)
                  IF ( ISOGAS(ISO,IGAS) .EQ. IGAS ) THEN
                    WRITE ( LUNPRF, '(3A,I1,A)', IOSTAT=IOS, ERR=900 ) 
     &                '*', CODGAS(IGAS)(1:IEND), '(', ISO, ') [ppmv]'
                    GOTO 100
                  END IF
                END DO
                STOP 'F-RFMPRF: Logical error'
              ENDIF
 100          CONTINUE
            END IF
            WRITE ( LUNPRF, *, IOSTAT=IOS, ERR=900 ) 
     &        ( VMRATM(IATM,IGAS) * 1.0E6, IATM = 1, NATM )
          END IF
        END DO
C
        WRITE ( LUNPRF, '(A)', IOSTAT=IOS, ERR=900 ) '*END'
C
      END DO
C
      CLOSE ( LUNPRF, IOSTAT=IOS, ERR=900 )
      IF ( GRAFLG ) CALL MOVGRA ( NOMGRA, 0 )          
      FAIL = .FALSE.
      RETURN
C
  900 FAIL = .TRUE.
      WRITE ( ERRMSG, '(A,I11)' ) 
     &  'F-RFMPRF: I/O failure on output file. IOSTAT=', IOS
      END
