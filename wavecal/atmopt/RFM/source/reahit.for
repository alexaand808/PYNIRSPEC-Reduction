      SUBROUTINE REAHIT ( EOF, FAIL, ERRMSG ) 
C
C VERSION
C     07AUG13 AD Combine AI and REF as local variable SPARE9 
C     09JAN08 AD Remove special case for HDO
C     29APR04 AD Add check on SHPGAS of isotopes
C     13JUN02 AD Allow for case no HITRAN file being used
C     12JUN02 AD Simplify: assume single HITRAN file. Remove IHFL argument
C     30APR02 AD Allow for isotopes
C     16FEB02 AD Allow for HDO
C     01DEC97 AD Check for SHPGAS = SHPCTM before setting record pointer
C     02OCT97 AD Change handling of Qualifier data (...QAL variables)
C     25JUL97 AD Add check for isotope# = 0
C     22JUL97 AD Add checks for line data qualifiers
C     03JUL97 AD Add test for SHPLUT
C     02APR97 AD Correction to allow for LSTAT values .LE. -10 
C                    (matters when HITRAN data contains superseded lines).
C     03MAR97 AD Version 3.
C     01OCT96 AD Version 2.
C     20SEP96 AD Rename VIBFLG to NTEFLG, SETVIB to SETNTE
C     01SEP96 AD Version 1.
C
C DESCRIPTION
C     Read HITRAN line data file for RFM.
C     Called by REACYC, RFMWID.
C
      IMPLICIT NONE
C
C ARGUMENTS      
      LOGICAL      EOF    !  O  Set TRUE if end-of-file reached
      LOGICAL      FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  RFMLOG ! Write text message to RFM log file.
     &, SETNTE ! Set IUSNTE and ILSNTE in /HITCOM/ if line affected by non-LTE.
C
C GLOBAL CONSTANTS
      INCLUDE 'idxcon.inc' ! HITRAN molecular indices
      INCLUDE 'rfmsiz.inc'  !  RFM array sizes
C
C LOCAL CONSTANTS
      INTEGER MAXIGI       ! Max no different warnings for unknown isotopes
        PARAMETER ( MAXIGI = 10 )
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'hflcom.inc' ! HITRAN line data file
      INCLUDE 'hitcom.inc' ! HITRAN line database variables
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'qalcom.inc' ! List of band/isotope qualifiers for line data.
      INCLUDE 'shpcon.inc' ! Line-shape codes
C
C LOCAL VARIABLES
      INTEGER      IBND    ! Selected band counter
      INTEGER      IGAS    ! Counter for gases to be used in RFM
      INTEGER      IGI     ! Combination of IGAS, ISO for warning messages
      INTEGER      IOS     ! Saved value of IOSTAT for error messages
      INTEGER      IREC    ! Record# 
      INTEGER      IIGI    ! Counter for warning messages
      INTEGER      NIGI    ! Number of warning messages so far
      INTEGER      IGILST(MAXIGI) ! List of unrecognised gas,isotope combs.
      LOGICAL      USELIN  ! T=Use this line 
      CHARACTER*80 MESSGE  ! Text of message sent to log file
      CHARACTER*9  SPARE9  ! Spare characters in binary record
C
C DATA STATEMENTS
      DATA NIGI / 0 /
      DATA IGILST / MAXIGI * 0 /
      SAVE NIGI, IGILST      
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      FAIL = .FALSE.
      EOF = LUNHFL .EQ. 0           ! Flag for no HITRAN data file being used
      IF ( EOF ) RETURN 
C
C Find next usefule record in file. Consider only gases for which SHPGAS is
C one of the line shapes, and not SHPCTM, SHPXSC or SHPLUT.
      DO WHILE ( .NOT. EOF ) 
        IREC = IR2HFL
        DO IGAS = 1, NGAS
C Assume SHPLUT > SHPXSC > SHPCTM and IFPHFL=0 for isotopes
          IF ( SHPGAS(IGAS) .LT. SHPCTM .AND. IFPHFL(IGAS) .GT. 0 ) 
     &      IREC = MIN ( IREC, IFPHFL(IGAS) )
        END DO
        EOF = IREC .EQ. IR2HFL
C
        IF ( .NOT. EOF ) THEN
          READ ( LUNHFL, REC=IREC, IOSTAT=IOS, ERR=900 ) LSTAT
C
          IF ( ABS ( LSTAT ) .GE. 10 ) THEN               ! Found a line record
            READ ( LUNHFL, REC=IREC, IOSTAT=IOS, ERR=900 )      ! Load /HITCOM/
     &        LSTAT, IDGAS, ISO, WNUM, STREN, TPROB, ABROAD, SBROAD,
     &        ELS, ABCOEF, TSP, IUSGQ, ILSGQ, USLQ, BSLQ, SPARE9, IFWDPT
C
C See if this line belongs to a required gas (should do) 
            IF ( IGSMOL(IDGAS) .GT. 0 ) THEN
C
C Set/Update forward pointer for the next call 
              IGAS = IGSMOL(IDGAS) 
              IFPHFL(IGAS) = IFWDPT + IREC
C
C Check if this line is from latest version of data in HITRAN file
              USELIN = ( LSTAT .GE. 10 )
C
C Check that isotope# is recognised
              IF ( ISO .LT. 1 .OR. ISO .GT. NISGAS(IGAS) ) THEN
                USELIN = .FALSE.
                IF ( NIGI .LT. MAXIGI ) THEN
                  IGI = IGAS*100 + ISO
                  DO IIGI = 1, NIGI
                    IF ( IGI .EQ. IGILST(IIGI) ) GOTO 100
                  END DO
                  NIGI = NIGI + 1
                  IGILST(NIGI) = IGI
                  WRITE ( MESSGE, '(A,F10.4,A,I6)' )
     &              'W-REAHIT: '//CODGAS(IGAS)//' Line at', WNUM,
     &              ' ignored - unrecognised Isotope#=', ISO
                  CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
                  IF ( FAIL ) RETURN
                END IF
 100            CONTINUE
              END IF
C
C Check isotope is also marked for line-by-line calculation
              IF ( USELIN ) 
     &          USELIN = SHPGAS(ISOGAS(ISO,IGAS)) .LT. SHPCTM
C
C If any qualifiers apply to this molecule, check it meets selection criteria
              IF ( USELIN .AND. USEQAL(IGAS) ) THEN
                DO IBND = 1, NBND(ISO,IGAS)
                  IF ( IUSQAL(IBND,ISO,IGAS) .EQ. 0 .AND. 
     &                 ILSQAL(IBND,ISO,IGAS) .EQ. 0 ) GOTO 200  ! selected 
                  IF ( IUSGQ .EQ. IUSQAL(IBND,ISO,IGAS) .AND.
     &                 ILSGQ .EQ. ILSQAL(IBND,ISO,IGAS) ) GOTO 200 
                  IF ( ILSQAL(IBND,ISO,IGAS) .EQ. 0 .AND.
     &                 IUSGQ .EQ. IUSQAL(IBND,ISO,IGAS) ) GOTO 200
                  IF ( IUSQAL(IBND,ISO,IGAS) .EQ. 0 .AND.
     &                 ILSGQ .EQ. ILSQAL(IBND,ISO,IGAS) ) GOTO 200
                END DO
                USELIN = .FALSE.                           ! Band not selected
  200           CONTINUE
              END IF                                  
C
              IF ( USELIN ) THEN
                IF ( ISO .GE. 1 .AND. ISO .LE. NISGAS(IGAS) ) THEN
                  IRCHFL = IREC
                  IF ( NTEFLG ) CALL SETNTE
                  RETURN             ! Exit with line loaded in /HITCOM/ buffer
                ELSE
                  WRITE ( MESSGE, '(A,F12.6,A,I11)' )
     &              'W-REAHIT: '//CODGAS(IGAS)//' Line at', WNUM,
     &              ' ignored - unlisted Isotope#=', ISO
                  CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
                  IF ( FAIL ) RETURN
                END IF
              END IF
            END IF
          ELSE IF ( LSTAT .EQ. -2 ) THEN
            EOF = .TRUE.
          END IF
        END IF
      END DO
C
C Return with EOF = TRUE and no usable data in /HITCOM/
C
      RETURN
C
  900 FAIL = .TRUE.
      WRITE ( ERRMSG, '(A,I11,A,I11)' )
     &  'F-REAHIT: Failed to read Rec#', IREC, 
     &  ' in HITRAN File. IOSTAT=', IOS
C
      END
