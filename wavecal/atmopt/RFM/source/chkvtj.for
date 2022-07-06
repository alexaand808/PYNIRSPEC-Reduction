      SUBROUTINE CHKVTJ ( FAIL, ERRMSG )
C
C VERSION
C     03-MAY-05  AD  Original.
C
C DESCRIPTION
C     Check vibrational temperature Jacobians
C     Called once by INPCHK if both JAC and NTE flags enabled.
C     Checks that vibrational temperature profiles have been loaded for all
C     the VT Jacobians specified in the *JAC section and resets IGSJAC to
C     point to the appropriate profile.
C
      IMPLICIT NONE
C
C ARGUMENTS
      LOGICAL       FAIL       !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG     !  O  Error message returned if FAIL is TRUE
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'jdxcon.inc' ! Indices for non-VMR Jacobians
C
C COMMON VARIABLES
      INCLUDE 'jaccom.inc' ! Jacobian data
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'ntecom.inc' ! Non-LTE data 
C
C LOCAL VARIABLES
      INTEGER      IJAC    ! Counter for Jacobians
      INTEGER      IDG     ! HITRAN species index for VT Jacobian 
      INTEGER      ISO     ! HITRAN isotope index for VT Jacobian 
      INTEGER      IGQ     ! HITRAN Global Quantum index for VT Jacobian
      INTEGER      INTE    ! Counter for NTE Vib.Temp profiles
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      DO IJAC = 1, NJAC
        IF ( IGSJAC(IJAC) .GT. 1000000 ) THEN    ! Decode into gasid,iso,gq
C Decoding needs to match coding in subroutine JACISO
          IDG = IGSJAC(IJAC) / 1000000
          ISO = MOD ( IGSJAC(IJAC), 1000000 ) / 1000
          IGQ = MOD ( IGSJAC(IJAC), 1000 )
          DO INTE = 1, NNTE
            IF ( IDG .EQ. IDGNTE(INTE) .AND. ISO .EQ. ISONTE(INTE) 
     &           .AND. IGQ .EQ. IGQNTE(INTE) ) THEN
C Recoding needs to match decoding in subroutine PTBNTE
              IGSJAC(IJAC) = MAXGAS + JDXNTE + INTE
              GOTO 100
            END IF
          END DO
          FAIL = .TRUE.
          WRITE ( ERRMSG, '(A,A,A,I3,A,I3)' )
     &      'F-CHKVTJ: No VT profile loaded for gas=', 
     &      CODGAS(IGSMOL(IDG)), ' isotope=', ISO, ' GQIndex=', IGQ
          RETURN
        END IF
 100    CONTINUE
      END DO
C
      FAIL   = .FALSE.
      END      
