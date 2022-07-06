      SUBROUTINE GASALL ( ALLSTR, FINAL, FAIL, ERRMSG )
C
C VERSION
C     15OCT13 AD SAVE IMIN
C     14AUG13 AD Add IMIN argument to LSTABS
C                    Use MOLIDX rather than local array of molecule names
C                    Remove checks for clono2/q sf6/q etc
C                    Write CODMOL without any leading blanks
C     23JUL09 AD Additional molecules (as in chkgas.for)
C     06MAY05 AD Save IMIN
C     06MAR02 AD Add SF6 as molecule#64, and HDO as molecule#39
C     14DEC01 AD Ensure JABS initialised
C     07FEB01 AD Add extra CFCs/HFCs ID=70-81, plus molecules 40-45
C     21DEC98 AD Original.
C
C DESCRIPTION
C     Add wildcard gases to list of absorbers
C     Called twice by INPGAS if a wildcard symbol '*' is present in *GAS 
C     section.
C
      IMPLICIT NONE
C
C ARGUMENTS 
      CHARACTER*(*) ALLSTR !  I  Text string from *GAS section
      LOGICAL       FINAL  !  I  True = Final Call to this module
      LOGICAL       FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG !  O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  GASCHK ! Check Gas name and set indices.
     &, LSTABS ! Construct ordered list of absorbers for given spec.range
     &, MOLIDX !  Give molecule name for HITRAN/RFM index, or vice-versa
     &, RFMLOG ! Write text message to RFM log file.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
C
C LOCAL CONSTANTS
      INTEGER       MAXLST ! Max.no. molecules returned by LSTABS
        PARAMETER ( MAXLST = 45 )  ! warning issued if this is limiting
      INTEGER       MAXSTR ! Max.length of string listing all abs.species
        PARAMETER ( MAXSTR = 400 ) ! Guess OK 
C
C LOCAL VARIABLES
      INTEGER       IABS   ! Absorber counter
      INTEGER       IOS    ! Saved value of IOSTAT for error messages
      INTEGER       IPT    ! Pointer to location within text string
      INTEGER       IMIN   ! Minimum strength to use
      INTEGER       IMOL   ! HITRAN index of each molecule from LSTABS
      INTEGER       ISPC   ! Spectral range counter
      INTEGER       JPT    ! Pointer to location within text string
      INTEGER       LPT    ! No.of characters in ALLSTR
      INTEGER       MOLLST(MAXLST) ! List of molecule IDs
      INTEGER       NABS   ! No. of molecules required/found by LSTABS
      REAL          OPTLST(MAXLST) ! List of "optical strengths"
      REAL          WNOH   ! Upper wavenumber of each spectral range 
      REAL          WNOL   ! Lower wavenumber of each spectral range 
      CHARACTER*7   CODMOL ! Molecules name or index
      CHARACTER*80  LOGMSG ! Information message for log file
      CHARACTER*(MAXSTR) ABSSTR ! Text string listing all absorbers
      LOGICAL       FIRST  ! TRUE=first call of this routine
      LOGICAL       NOTGAS ! TRUE=string supplied to GASCHK not recognised gas
C
C DATA STATEMENTS
      DATA FIRST / .TRUE. /
      DATA IMIN / 0 /  ! default minimum if not specified via quantifier

      SAVE FIRST, IMIN
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C The first part (FINAL=FALSE) is executed for any occurrence of a text field
C starting with '*' in the *GAS section of the driver table (should only be 
C one at most).
      IF ( .NOT. FINAL ) THEN
        IF ( .NOT. FIRST ) THEN
          FAIL = .TRUE.
          ERRMSG = 'F-GASALL: Repeated wildcard ''*'' in *GAS section'
          RETURN
        END IF
        FIRST = .FALSE.
        LPT = LEN ( ALLSTR )
        IF ( LPT .GT. 5 ) THEN
          FAIL = .TRUE.
          LPT = MIN ( LPT, 10 )
          ERRMSG = 'F-GASALL: Invalid wildcard in *GAS section '//
     &             '(string too long): '//ALLSTR(1:LPT)
        END IF
        LOGMSG = 'I-GASALL: Include all absorbers'
        IF ( LPT .GT. 1 ) THEN
          IPT = INDEX ( ALLSTR, '(' )
          JPT = INDEX ( ALLSTR, ')' ) 
          IF ( IPT .NE. 2 .OR. JPT .LT. 4 .OR. JPT .GT. 5 ) THEN
            FAIL = .TRUE.
            ERRMSG = 'F-GASALL: Invalid wildcard in *GAS section '//
     &               '(Brackets not as expected): '//ALLSTR(1:LPT)
            RETURN
          END IF
          IF ( JPT .EQ. 4 ) THEN
            READ ( ALLSTR(3:4), '(I1)', ERR=900, IOSTAT=IOS ) IMIN
          ELSE 
            READ ( ALLSTR(3:5), '(I2)', ERR=900, IOSTAT=IOS ) IMIN
          END IF
          IF ( IMIN .LT. 0 ) THEN
            FAIL = .TRUE.
            ERRMSG = 'F-GASALL: Invalid wildcard in *GAS section '//
     &             '(-ve qualifier not allowed): '//ALLSTR(1:LPT)
            RETURN
          END IF
          WRITE ( LOGMSG, '(A,I2)' ) 'I-GASALL: Include all '//
     &      'absorbers with ''optical strength'' .GE. ', IMIN  
        END IF
        CALL RFMLOG ( LOGMSG, FAIL, ERRMSG )
        RETURN
      END IF
C
C Following is to be executed (FINAL=TRUE) after all explicit gases have been
C read from *GAS section
      IF ( FIRST ) THEN          ! No extra "wildcard" gases to be added
        FAIL = .FALSE.
        RETURN                   
      END IF
C
      DO ISPC = 1, NSPC
        WNOL = SNGL ( WNLSPC(ISPC) )
        WNOH = SNGL ( WNUSPC(ISPC) )
        NABS = MAXLST
        CALL LSTABS ( WNOL, WNOH, IMIN, NABS, MOLLST, OPTLST )
C
C Check if no.absorbers returned by LSTABS limited by MAXLST array size
        IF ( NABS .EQ. MAXLST ) THEN
          WRITE ( LOGMSG, '(A,I3,A)' ) 
     &      'W-GASALL: No.absorbers found limited to ',  MAXLST,
     &      '.Increase local parameter MAXLST?'
          CALL RFMLOG ( LOGMSG, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        END IF
C
C Establish all required absorbers for this spectral range and write to log file
        WRITE ( LOGMSG, '(A,I2,A,A)' ) 
     &    'I-GASALL: Absorbers(Strength .GE.', IMIN, 
     &    ') for spectral range: ', LABSPC(ISPC)
        CALL RFMLOG ( LOGMSG, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        ABSSTR = ' '
        IPT = 5
        DO IABS = 1, NABS
          IMOL = MOLLST(IABS)
          CALL MOLIDX ( IMOL, CODMOL, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
          JPT = IPT + INDEX ( CODMOL, ' ' ) - 2
          IF ( JPT .LT. IPT ) JPT = IPT + 6
          IF ( JPT .GT. MAXSTR-4 ) THEN ! can't fit anymore text into ABSSTR
            IPT = MAXSTR-3              ! should ensure no more text attempted
            JPT = MAXSTR
            ABSSTR(IPT:JPT) = 'etc.'
            WRITE ( *, * ) 'w-gasall: MAXSTR too small' ! increase MAXSTR
          ELSE
            ABSSTR(IPT:JPT) = CODMOL
            IF ( OPTLST(IABS) .GE. 10 ) THEN
              WRITE ( ABSSTR(JPT+1:JPT+4), '(A,I2,A)' ) 
     &          '(', INT(OPTLST(IABS)), ')' 
              IPT = JPT + 6
            ELSE
              WRITE ( ABSSTR(JPT+1:JPT+3), '(A,I1,A)' )
     &          '(', INT(OPTLST(IABS)), ')' 
              IPT = JPT + 5
            END IF
          END IF
        END DO
        CALL RFMLOG ( ABSSTR(1:IPT), FAIL, ERRMSG )
        IF ( FAIL ) RETURN
C
C For each absorber not already included, add to list
        DO IABS = 1, NABS
          IMOL = MOLLST(IABS)
          IF ( IGSMOL(IMOL) .EQ. 0 ) THEN    ! Not already specified by user
            IF ( IMOL .LE. 9 ) THEN  ! Add to list (without any leading blanks)
              WRITE ( CODMOL, '(I1)' ) IMOL       
            ELSE IF ( IMOL .LE. 99 ) THEN
              WRITE ( CODMOL, '(I2)' ) IMOL
            ELSE 
              WRITE ( CODMOL, '(I3)' ) IMOL
            END IF
            CALL GASCHK ( CODMOL, NOTGAS, FAIL, ERRMSG )
C Since the only molecule indices sent from this routine to GASCHK are from the
C LSTABS module, they should all be recognised, so NOTGAS=T shouldn't happen.
            IF ( NOTGAS ) STOP 'F-GASALL: Logical Error'
            IF ( FAIL ) RETURN
         END IF
        END DO           ! End loop over potential new absorbers
      END DO             ! End loop over spectral ranges
C
      RETURN
  900 CONTINUE
      FAIL = .TRUE.
      WRITE ( ERRMSG, '(A,I11)' )
     &  'F-GASALL: Error reading * qualifier in *GAS section,'//
     &  'IOSTAT=', IOS
      END
