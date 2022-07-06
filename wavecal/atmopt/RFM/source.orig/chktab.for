      SUBROUTINE CHKTAB ( VALTAB, FAIL, ERRMSG )
C
C VERSION
C     19-FEB-03  AD  Comment change only
C     11-AUG-99  AD  Use explicit SAVE of local variables rather than general
C     28-APR-99  AD  Modified to include temperature offset profile
C     23-JUN-98  AD  Correction to Bug Fix: set LP1TAB correctly
C     24-APR-98  AD  Bug Fix: ensure this works for NPRE or NTEM = 1
C     22-OCT-97  AD  Original.
C
C DESCRIPTION
C     Check Tabulation parameters.
C     Called by TANCHK for each of 6 fields in *TAN / *DIM section.
C
      IMPLICIT NONE
C
      EXTERNAL
     &  RFMLOG ! Write text message to RFM log file
C
C ARGUMENTS
      REAL          VALTAB  !  I  Tabulation parameter to be tested 
      LOGICAL       FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG  !  O  Error message returned if FAIL is TRUE
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'tancom.inc' ! Tangent heights
      INCLUDE 'tabcom.inc' ! Table axes for tabulated absorb.coeffs.
C
C LOCAL VARIABLES
      REAL    PRE1, PRE2  ! Pressure limits [mb] for tabulation
      REAL    TEM2        ! 2nd temperature limit [K] for tabulation
      CHARACTER*80 LOGMSG ! Text string written to Log file.
      SAVE    PRE1, PRE2, TEM2
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C NTAN is just used to keep track of the 6 required tabulation parameters
      NTAN = NTAN + 1
C
C Check if more than 6 tabulation parameters have been supplied
C
      IF ( NTAN .GT. 6 ) THEN
        WRITE ( ERRMSG, '(A,G12.3)' ) 'F-CHKTAB: '//
     &  'Unexpected extra number in *TAN/*DIM section, value=', VALTAB
        FAIL = .TRUE.
        RETURN
      END IF
C
C All tabulation parameters must be +ve to be valid, except for temperature
C limits: -ve, assume the offset profile is to be applied
C 
      IF ( VALTAB .LE. 0.0 .AND. ( NTAN .NE. 4 .AND. NTAN .NE. 5 )) THEN
        WRITE ( ERRMSG, '(A,I1,A,G12.3)' ) 
     &    'F-CHKTAB: Tabulation parameter#', NTAN, 
     &    ' in *TAN/*DIM section is not +ve, =', VALTAB
        FAIL = .TRUE.
        RETURN
      END IF
C
C First and second parameters are high,low pressure limits (or low,high)
C
      IF ( NTAN .EQ. 1 .OR. NTAN .EQ. 2 ) THEN         
        IF ( VALTAB .GT. 1014.0 .OR. VALTAB .LT. 0.0001 ) THEN
          WRITE ( LOGMSG, '(A,G12.3)' ) 'W-CHKTAB: '//
     &    'Dubious value for pressure [mb] limit: ', VALTAB
          CALL RFMLOG ( LOGMSG, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        END IF
        IF ( .NOT. HOMFLG ) THEN
          IF ( VALTAB .GT. PREATM(1) .OR. 
     &         VALTAB .LT. PREATM(NATM)   ) THEN
            WRITE ( ERRMSG, '(A,G12.3)' ) 'F-CHKTAB: Press.Tabulation'
     &        //' limit outside atmospheric profile, value=', VALTAB
            FAIL = .TRUE.
            RETURN
          END IF
        END IF
        IF ( NTAN .EQ. 1 ) THEN 
          PRE1 = VALTAB 
          LP1TAB = - LOG ( PRE1 )
        ELSE
          PRE2 = VALTAB 
        END IF
C
C Third parameter is number of points in pressure domain
C
      ELSE IF ( NTAN .EQ. 3 ) THEN
        IF ( MOD ( VALTAB, 1.0 ) .NE. 0.0 ) THEN
          WRITE ( ERRMSG, '(A,G12.3)' ) 'F-CHKTAB: '//
     &    'No.Press.pts is not an integer, value=', VALTAB
          FAIL = .TRUE.
          RETURN 
        END IF
        NLPTAB = NINT ( VALTAB )
        IF ( NLPTAB .EQ. 1 ) THEN                ! Single pressure point
          IF ( PRE1 .NE. PRE2 ) THEN
            ERRMSG = 'F-CHKTAB: Single pressure point '//
     &               'but Lower & Upper pressures are not identical'
            FAIL = .TRUE.
            RETURN
          END IF
          WRITE ( LOGMSG, '(A,G12.3,A)' ) 
     &      'I-CHKTAB: Pressure tabulation: ', PRE1, 
     &      ' [mb], 1 point only'
          LPDTAB = 0.0
        ELSE                                     ! 2 or more pressure points
          WRITE ( LOGMSG, '(A,G12.3,A,G12.3,A,I6,A)' ) 
     &      'I-CHKTAB: Pressure tabulation: ', PRE1, ' to ',
     &      PRE2, ' [mb], ', NLPTAB, ' points'
          LPDTAB = ( - LOG ( PRE2 ) - LP1TAB ) / ( NLPTAB - 1 )
        END IF
        CALL RFMLOG ( LOGMSG, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
C
C Fourth and Fifth parameters are high,low temperature limits (or low,high)
C
      ELSE IF ( NTAN .EQ. 4 .OR. NTAN .EQ. 5 ) THEN    
        IF ( NTAN .EQ. 4 .AND. VALTAB .LT. 0 ) THEN
          OFFTAB = .TRUE.
        ELSE IF ( NTAN .EQ. 5 .AND. VALTAB .LT. 0.0 ) THEN
          IF ( OFFTAB ) THEN
            ERRMSG = 'F-CHKTAB: both temperature limits in '//
     &               '*TAN/*DIM section are -ve'
            FAIL = .TRUE.
            RETURN
          ELSE
            OFFTAB = .TRUE.
          END IF
        END IF
C
        IF ( ( OFFTAB .AND. ABS(VALTAB) .GT. 100.0 ) .OR.
     &       ( .NOT. OFFTAB .AND. VALTAB .GT. 350.0 ) .OR. 
     &       ( .NOT. OFFTAB .AND. VALTAB .LT. 120.0 )      ) THEN 
          WRITE ( LOGMSG, '(A,G12.3)' ) 'W-CHKTAB: '//
     &    'Dubious value for temperature [K] limit: ', VALTAB
          CALL RFMLOG ( LOGMSG, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        END IF
        IF ( NTAN .EQ. 4 ) THEN 
          TM1TAB = VALTAB 
        ELSE
          TEM2 = VALTAB 
        END IF
C
C Sixth parameter is number of points in temperature domain
C
      ELSE IF ( NTAN .EQ. 6 ) THEN
        IF ( MOD ( VALTAB, 1.0 ) .NE. 0.0 ) THEN
          WRITE ( ERRMSG, '(A,G12.3)' ) 'F-CHKTAB: '//
     &    'No.Temp.pts is not an integer, value=', VALTAB
          FAIL = .TRUE.
          RETURN 
        END IF
        NTMTAB = NINT ( VALTAB )
        IF ( NTMTAB .EQ. 1 ) THEN                ! Single temperature point
          IF ( TM1TAB .NE. TEM2 ) THEN
            ERRMSG = 'F-CHKTAB: Single temperature point '//
     &               'but Lower & Upper temps. are not identical'
            FAIL = .TRUE.
            RETURN
          END IF
          WRITE ( LOGMSG, '(A,F10.3,A)' ) 
     &      'I-CHKTAB: Temperature tabulation: ', TM1TAB, 
     &      ' [K], 1 point only'
          TMDTAB = 0.0
        ELSE                                      ! 2 or more temperature points
          WRITE ( LOGMSG, '(A,F10.3,A,F10.3,A,I6,A)' ) 
     &      'I-CHKTAB: Temperature tabulation: ', TM1TAB, ' to ',
     &      TEM2, ' [K], ', NTMTAB, ' points'
          TMDTAB = ( TEM2 - TM1TAB ) / ( NTMTAB - 1 )
        END IF
        CALL RFMLOG ( LOGMSG, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
        IF ( OFFTAB ) THEN
          LOGMSG = 'I-CHKTAB: Assuming temperature offset profile'
          CALL RFMLOG ( LOGMSG, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
        END IF
      END IF
      FAIL = .FALSE.
C
      END
C-------------------------------------------------------------------------------
C                                NOTICE
C
C     This software module is part of the MIPAS Reference Forward Model supplied
C     to the European Space Agency under ESTEC contract no. 11886/96/NL/GS.
C        
C     All rights, title and interest in and to the copyright
C     and any other proprietary right, is held by
C     The University Corporation for Atmospheric Research, having a
C     correspondence address of P.O. Box 3000, Boulder, Colorado 80307, USA.
C
C     However, note that all inquiries concerning the MIPAS
C     Reference Forward Model should be submitted to the The Manager (Projects),
C     AOPP,  Clarendon Laboratory, Parks Road, Oxford, OX1 3PU, UK.
C     (Tel: +44-8165-272900,    Fax: +44-1865-272924).  
C
C-------------------------------------------------------------------------------
