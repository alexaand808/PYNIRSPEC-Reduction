      SUBROUTINE ATMCHK ( GOTTEM, GOTPRE, GOTGAS, FAIL, ERRMSG )
C
C VERSION
C     03-APR-13  AD  Change sum VMR > 1 from fatal to warning
C     17-JAN-13  AD  Bug#81: write warnings for -ve aerosol extinction
C     20-DEC-01  AD  Replace "PRFCHK" with "ATMCHK" in error messages
C     06-JUN-01  AD  Original. Replaces old subroutine PRFCHK.
C
C DESCRIPTION
C     Check atmospheric profiles are reasonable.
C     Called by ATMFIL for every set of profiles loaded from a file.
C
      IMPLICIT NONE
C
C ARGUMENTS
      LOGICAL       GOTTEM    !  I  Set TRUE if Temperature profile read in
      LOGICAL       GOTPRE    !  I  Set TRUE if Pressure profile read in
      LOGICAL       GOTGAS(*) !  I  Set TRUE if vmr profile read in
      LOGICAL       FAIL      !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG    !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  RFMLOG ! Write text message to RFM log file.
C
C GLOBAL CONSTANTS
      INCLUDE 'idxcon.inc' ! HITRAN molecular indices
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'gascom.inc' ! Molecular and isotope data
C
C LOCAL VARIABLES
      INTEGER IATM         ! Atmospheric level index
      INTEGER IGAS         ! Absorber index
      REAL    SUM          ! Summation for composition at each level
      CHARACTER*80 WRNMSG  ! Warning message
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Check that pressures are +ve and decrease monotonically
      IF ( GOTPRE ) THEN       ! PRE profile loaded
        DO IATM = 1, NATM
          IF ( PREATM(IATM) .LE. 0.0 ) THEN
            WRITE ( ERRMSG, '(A,1PG9.2,A,G9.2,A)' )
     &        'F-ATMCHK: Illegal PRE profile value =', PREATM(IATM),
     &        ' at HGT ', HGTATM(IATM), ' km'
            FAIL = .TRUE.
            RETURN
          END IF
        END DO
        DO IATM = 2, NATM
          IF ( PREATM(IATM) .GE. PREATM(IATM-1) ) THEN
            ERRMSG = 
     &       'F-ATMCHK: PRE  profile doesn''t decrease monotonically'
            FAIL = .TRUE.
            RETURN
          END IF
        END DO
      END IF
C
C Check that temperatures are +ve
      IF ( GOTTEM ) THEN      ! TEM profile
        DO IATM = 1, NATM            
          IF ( TEMATM(IATM) .LE. 0.0 ) THEN
            WRITE ( ERRMSG, '(A,1PG9.2,A,G9.2,A)' )
     &        'F-ATMCHK: Illegal TEM profile value =', TEMATM(IATM),
     &        ' at HGT ', HGTATM(IATM), ' km'
            FAIL = .TRUE.
            RETURN
          END IF
        END DO
      END IF

C Check that total VMR doesn't exceed 1ppv.
      DO IATM = 1, NATM
        SUM = 0.0
        DO IGAS = 1, NGAS 
          IF ( GOTGAS(IGAS) ) THEN
            IF ( VMRATM(IATM,IGAS) .LT. 0.0 ) THEN
C Allow negative aerosol values (with warning) but not VMR
              IF ( IDXGAS(IGAS) .EQ. IDXAER ) THEN          
                WRITE ( WRNMSG, '(A,1PG9.2,A,G9.2,A)' )
     &            'W-ATMCHK: Negative '//CODGAS(IGAS)//
     &            ' profile value =', VMRATM(IATM,IGAS),
     &            ' at HGT ', HGTATM(IATM), ' km'
                CALL RFMLOG ( WRNMSG, FAIL, ERRMSG )
                IF ( FAIL ) RETURN
              ELSE
                WRITE ( ERRMSG, '(A,1PG9.2,A,G9.2,A)' )
     &          'F-ATMCHK: Negative '//CODGAS(IGAS)//
     &          ' profile value =', VMRATM(IATM,IGAS),
     &          ' at HGT ', HGTATM(IATM), ' km'
                FAIL = .TRUE.
                RETURN
              END IF
            END IF
            IF ( IDXGAS(IGAS) .NE. IDXAER ) THEN          
              SUM = SUM + VMRATM(IATM,IGAS) 
              IF ( SUM .GT. 1.0 ) THEN
                WRITE ( WRNMSG, '(A,A7,A,G9.2,A)' )
     &            'W-ATMCHK: Total VMR > 1ppv adding GAS=', 
     &            CODGAS(IGAS), ' at HGT=', HGTATM(IATM), ' km'
                CALL RFMLOG ( WRNMSG, FAIL, ERRMSG )  
                IF ( FAIL ) RETURN
                RETURN
              END IF
            END IF
          END IF
        END DO
      END DO
C
      FAIL = .FALSE.
      END
