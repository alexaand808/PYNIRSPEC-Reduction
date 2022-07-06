      SUBROUTINE GASQAL ( QALSTR, FAIL, ERRMSG )
C
C VERSION
C     31-JAN-04  AD  Correct format statement for error message
C     22-JAN-98  AD  Split concatenation so FORTRAN to C converter can cope
C     02-OCT-97  AD  Allow band-isotope combinations.
C     22-JUL-97  AD  Original
C
C DESCRIPTION
C     Check & load isotope/band qualifiers for each molecule 
C     Called by GASCHK for each molecule with appended qualifier string.
C     This module assumes that qualifiers refer to last gas loaded into
C     /GASCOM/ with index NGAS.
C     The layout of qualifiers is, eg for CO2,
C       CO2(3:4)(1:2)(2)(1)(2:3)(3)(1:3)(2:4)(4)(5)
C     which will (a) select bands (3:4), (1:2) for all isotopes
C                (b) select band (2:3) for isotopes 1,2
C                (c) select bands (1:3), (2:3) for isotope 3
C                (d) select all bands for isotopes 4,5
C     The structure of this subroutine is:
C       Initialise ISOLST with ALL isotopes for gas
C       For each Qualifier ...
C         If Qualifier is a Band ...
C           Set band for all isotopes currently stored in ISOLST
C         Else Qualifier is an Isotope
C           If previous Qualifier was a band then erase current ISOLST
C           Add Isotope to ISOLST
C         End if
C       End loop over Qualifiers
C       If last Qualifier was an Isotope then set all bands for current ISOLST
C
      IMPLICIT NONE
C
C ARGUMENTS
      CHARACTER*(*) QALSTR  !  I  String containing list of qualifiers
      LOGICAL       FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     & RFMLOG ! Write text message to RFM log file.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'qalcom.inc' ! List of band/isotope qualifiers for line data.
C
C LOCAL VARIABLES
      INTEGER   IBND ! Selected band counter
      INTEGER   IISO ! Isotope counter
      INTEGER   ILS  ! Selected lower state Global Quantum number
      INTEGER   IOS  ! Value of IOSTAT saved for error messages
      INTEGER   IPT  ! Pointer to '(' element of current qualifier 
      INTEGER   IQAL ! Qualifier counter
      INTEGER   ISO  ! Isotope counter/index
      INTEGER   ISOLST(MAXISO) ! List of isotopes using current band
      INTEGER   IUS  ! Selected upper state Global Quantum number
      INTEGER   JPT  ! Pointer to ':' element (if any) of current qualifier 
      INTEGER   KPT  ! Pointer to ')' element of current qualifier 
      INTEGER   LPT  ! Pointer to element of QALSTR
      INTEGER   MPT  ! Pointer to limit length of string written to ERRMSG
      INTEGER   NISO ! No.of isotopes in ISOLST
      LOGICAL   LASBND    ! T=last qualifier read was a band (F=an isotope)
      CHARACTER*30 ERRSTR ! Substring for error messages
      CHARACTER*80 LOGMSG ! Message sent to RFM log file
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      LPT = LEN ( QALSTR )
C
C Initialise for each new gas
      IQAL = 0
      IPT = 1
      KPT = 0
      LASBND = .TRUE. 
      NISO = NISGAS(NGAS)
      DO ISO = 1, NISO
        ISOLST(ISO) = ISO
        NBND(ISO,NGAS) = 0
      END DO
      WRITE ( LOGMSG, '(A)' )
     &  'I-GASQAL: Gas='//CODGAS(NGAS)//
     &      ' List of selected isotopes and/or bands ... '
      CALL RFMLOG ( LOGMSG, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
C Repeat until all qualifiers decoded from QALSTR
      DO WHILE ( KPT .LT. LPT )
        IQAL = IQAL + 1
C
C Set up substring for any error messages applicable to this gas/qualifier
        WRITE ( ERRSTR, '(A,I2,A)' )
     &    'F-GASQAL: Gas='//CODGAS(NGAS)//' Qual#', IQAL, ':'
C
C Determine limits IPT:KPT of substring containing qualifier
        IF ( QALSTR(IPT:IPT) .NE. '(' ) THEN
          FAIL = .TRUE.
          WRITE ( ERRMSG, '(A,A,A)' ) ERRSTR//
     &      'read error. Expected ''('',  found ''',
     &      QALSTR(IPT:IPT), ''''
          RETURN
        END IF
        KPT = INDEX ( QALSTR(IPT+1:LPT), ')' )
        IF ( KPT .EQ. 0 ) THEN
          MPT = MIN ( LPT, IPT+9 )
          FAIL = .TRUE.
          WRITE ( ERRMSG, '(A,A,A)' ) ERRSTR//
     &      'no '')'' found after ''', QALSTR(IPT:MPT), '...'''
          RETURN
        END IF
        KPT = KPT + IPT
C
C Read qualifier contents - start by looking for ':' implying band selection
        JPT = INDEX ( QALSTR(IPT+1:KPT-1), ':' ) 
        IF ( JPT .NE. 0 ) THEN                               
          LASBND = .TRUE.
          JPT = JPT + IPT
          IF ( QALSTR(IPT+1:JPT-1) .EQ. '*' ) THEN
            ILS = 0
          ELSE
            READ ( QALSTR(IPT+1:JPT-1), *, ERR=900, IOSTAT=IOS ) ILS
            IF ( ILS .LE. 0 .OR. ILS .GT. 999 ) THEN
              FAIL = .TRUE.
              WRITE ( ERRMSG, '(A,I11)' ) ERRSTR//
     &        'ILSGQ out of range 1:999, value=', ILS
              RETURN
            END IF
          END IF
          IF ( QALSTR(JPT+1:KPT-1) .EQ. '*' ) THEN
            IUS = 0
          ELSE
            READ ( QALSTR(JPT+1:KPT-1), *, ERR=900, IOSTAT=IOS ) IUS
            IF ( IUS .LE. 0 .OR. IUS .GT. 999 ) THEN
              FAIL = .TRUE.
              WRITE ( ERRMSG, '(A,I11)' ) ERRSTR//
     &        'IUSGQ out of range 1:999, value=', IUS
              RETURN
            END IF
          END IF
          DO IISO = 1, NISO
            ISO = ISOLST(IISO)
C
            DO IBND = 1, NBND(ISO,NGAS)
C Check this isn't a repeated band,isotope combination
              IF ( ILSQAL(IBND,ISO,NGAS) .EQ. ILS .AND.
     &             IUSQAL(IBND,ISO,NGAS) .EQ. IUS ) THEN
                FAIL = .TRUE.
                WRITE ( ERRMSG, '(A,A,A,I2)' ) ERRSTR//
     &          'repeated band,isotope ', QALSTR(IPT:KPT), ',', ISO
                RETURN
              END IF
              IF ( IUS .EQ. IUSQAL(IBND,ISO,NGAS) .AND. 
     &          ( ILS .EQ. 0 .OR. ILSQAL(IBND,ISO,NGAS) .EQ. 0 ) ) THEN
                FAIL = .TRUE.
                WRITE ( ERRMSG, '(A,I2)' ) ERRSTR//
     &           'Inconsistent use of (*: ), Isotope#', ISO 
                RETURN
              END IF
              IF ( ILS .EQ. ILSQAL(IBND,ISO,NGAS) .AND. 
     &          ( IUS .EQ. 0 .OR. IUSQAL(IBND,ISO,NGAS) .EQ. 0 ) ) THEN
                FAIL = .TRUE.
                WRITE ( ERRMSG, '(A,I2)' ) ERRSTR//
     &           'Inconsistent use of ( :*), Isotope#', ISO 
                RETURN
              END IF
            END DO
C
C Checks complete - load selected band and print info message to log file
            IF ( NBND(ISO,NGAS) .EQ. MAXBND ) THEN
              FAIL = .TRUE.
              WRITE ( ERRMSG, '(A,I11)' ) ERRSTR//
     &        'No.selected bands > MAXBND in RFMSIZ, =', MAXBND
              RETURN
            ELSE
              NBND(ISO,NGAS) = NBND(ISO,NGAS) + 1
              ILSQAL(IBND,ISO,NGAS) = ILS
              IUSQAL(IBND,ISO,NGAS) = IUS
            END IF
          END DO
          WRITE ( LOGMSG, '(A,A)' )
     &      'I-GASQAL:     Select lines from vib.band low:high=',
     &      QALSTR(IPT:KPT)
          CALL RFMLOG ( LOGMSG, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
C
C If no `:' found in qualifier, assume it contains an isotope number
        ELSE                                                      
          READ ( QALSTR(IPT+1:KPT-1), *, ERR=901, IOSTAT=IOS ) ISO
C
C Check selected isotope within valid range for molecule
          IF ( ISO .LE. 0 .OR. ISO .GT. NISGAS(NGAS) ) THEN
            FAIL = .TRUE.
            WRITE ( ERRMSG, '(A,I11,A,I4,A)' ) ERRSTR//
     &        'bad isotope#=', ISO, '  (range=1:', NISGAS(NGAS), ')'
            RETURN
C
C Selected isotope OK, so print info message to RFM log file
          ELSE
            WRITE ( LOGMSG, '(A,I4)' ) 
     &        'I-GASQAL:   Select isotope#', ISO
            CALL RFMLOG ( LOGMSG, FAIL, ERRMSG )
            IF ( FAIL ) RETURN
            IF ( LASBND ) THEN
              NISO = 0
              LASBND = .FALSE.
            END IF
            DO IISO = 1, NISO
              IF ( ISO .EQ. ISOLST(IISO) ) THEN
                FAIL = .TRUE.
                WRITE ( ERRMSG, '(A,I2)' ) ERRSTR//
     &          'repeated isotope#', ISO
                RETURN
              END IF
            END DO            
            NISO = NISO + 1
            ISOLST(NISO) = ISO
          END IF
        END IF
C
C Assume next character after ')' is start of next qualifier: '(' 
        IPT = KPT + 1
      END DO                           
C
C Successfully decoded entire QALSTR 
C
C Check if finished with an isotope list, in which case set ILS,IUS wildcards
C
      IF ( .NOT. LASBND ) THEN
        DO IISO = 1, NISO
          ISO = ISOLST(IISO)
          IF ( NBND(ISO,NGAS) .NE. 0 ) THEN
            FAIL = .TRUE.
            WRITE ( ERRMSG, '(A,I2)' ) ERRSTR//
     &      'bands already selected for isotope#', ISOLST(ISO)
            RETURN
          END IF
          NBND(ISO,NGAS) = 1
          ILSQAL(1,ISO,NGAS) = 0
          IUSQAL(1,ISO,NGAS) = 0
        END DO
        IF ( NISO .EQ. NISGAS(NGAS) ) THEN
          WRITE ( LOGMSG, '(A)' ) 
     &      'W-GASQAL: all isotopes selected'
          CALL RFMLOG ( LOGMSG, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        ELSE
          WRITE ( LOGMSG, '(A,I4)' ) 
     &      'I-GASQAL:     Select all bands'
          CALL RFMLOG ( LOGMSG, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        END IF
C
      END IF
C
      FAIL = .FALSE.
      RETURN
C
  900 WRITE ( ERRMSG, '(A,I11)' ) ERRSTR//
     &  'failed to read quantum numbers. IOSTAT=', IOS
      FAIL = .TRUE.
      RETURN
C
  901 WRITE ( ERRMSG, '(A,I11)' ) ERRSTR//
     &  'failed to read isotope#, IOSTAT=', IOS
      FAIL = .TRUE.
C
      END 
