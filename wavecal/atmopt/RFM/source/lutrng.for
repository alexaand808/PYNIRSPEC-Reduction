      SUBROUTINE LUTRNG ( ILUT, IGAS, FAIL, ERRMSG )
C
C VERSION
C     09-SEP-99  AD  Suppress warnings if NP or NT=1, also >max -lnp edge
C     01-AUG-97  AD  Allow for DPLUT, DTLUT being negative as well.
C     03-JUL-97  AD  Original.
C
C DESCRIPTION
C     Check p,T range of Look-Up Table data.
C     Called by SPCLUT for LUT file used.
C     Checks that p,T range of LUT encompasses all path p,T conditions required
C     for particular gas. If not, a warning is issued but calculation proceeds
C     (RFMLUT sets outside values to value at nearest edge of table)
C     Note that local calculation is in terms of Max/Min Pressure instead of
C     Min/Max -lnp so Max/Min appear the wrong way around.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER       ILUT    !  I  Index of Look-Up Table to be tested
      INTEGER       IGAS    !  I  Index of Gas for LUT file
      LOGICAL       FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  RFMLOG ! Write text message to RFM log file.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
C
C COMMON VARIABLES
      INCLUDE 'lutcom.inc' ! Look-Up Table data
      INCLUDE 'pthcom.inc' ! Look-Up Table file data
C
C LOCAL VARIABLES
      INTEGER   IPTH   ! Path counter
      INTEGER   NPMAX  ! No.times Max -lnp (Min Press) is exceeded 
      INTEGER   NPMIN  ! No.times Min -lnp (Max Press) is exceeded 
      INTEGER   NTMAX  ! No.times Max Temp is exceeded 
      INTEGER   NTMIN  ! No.times Min Temp is exceeded 
      INTEGER   NTOTAL ! Total no.calculated paths using this gas
      REAL      LNPMAX ! Max -lnp [p in mb] in LUT tabulation
      REAL      LNPMIN ! Min -lnp [p in mb] in LUT tabulation
      REAL      PREMAX ! Max Pressure [atm] in LUT tabulation
      REAL      PREMIN ! Min Pressure [atm] in LUT tabulation
      REAL      TEMMAX ! Max Temperature [K] in LUT tabulation
      REAL      TEMMIN ! Min Temperature [K] in LUT tabulation
      CHARACTER*80 WRNMSG ! Warning message sent to Log file
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      NPMIN = 0
      NPMAX = 0
      NTMIN = 0
      NTMAX = 0
      NTOTAL = 0
C
      IF ( DPLUT(ILUT) .GT. 0 ) THEN
        LNPMIN = P1LUT(ILUT)
        LNPMAX = P1LUT(ILUT) + ( NPLUT(ILUT) - 1 ) * DPLUT(ILUT)
      ELSE 
        LNPMAX = P1LUT(ILUT)
        LNPMIN = P1LUT(ILUT) + ( NPLUT(ILUT) - 1 ) * DPLUT(ILUT)
      END IF
      PREMAX = EXP ( -LNPMIN ) / ATMB
      PREMIN = EXP ( -LNPMAX ) / ATMB
      IF ( DTLUT(ILUT) .GT. 0 ) THEN
        TEMMIN = T1LUT(ILUT)
        TEMMAX = T1LUT(ILUT) + ( NTLUT(ILUT) - 1 ) * DTLUT(ILUT)
      ELSE
        TEMMAX = T1LUT(ILUT)
        TEMMIN = T1LUT(ILUT) + ( NTLUT(ILUT) - 1 ) * DTLUT(ILUT)
      END IF
C
      DO IPTH = 1, NCLC
        IF ( IGAS .EQ. IGSPTH(IPTH) ) THEN
          NTOTAL = NTOTAL + 1
          IF ( PREPTH(IPTH) .GT. PREMAX ) NPMIN = NPMIN + 1     
          IF ( PREPTH(IPTH) .LT. PREMIN ) NPMAX = NPMAX + 1     
          IF ( TEMPTH(IPTH) .GT. TEMMAX ) NTMAX = NTMAX + 1
          IF ( TEMPTH(IPTH) .LT. TEMMIN ) NTMIN = NTMIN + 1
        END IF
      END DO 
C suppress warnings for exceeding low-pressure (=high alt) edge of table
c      IF ( NPMAX .GT. 0 .AND. NPLUT(ILUT) .GT. 1 ) THEN
c        WRITE ( WRNMSG, '(A,I6,A1,I7,A,F5.1,A2)' )
c     &    'W-LUTRNG: Max -lnp edge of LUT exceeded for', NPMAX,
c     &    '/', NTOTAL, ' calc paths (', 
c     &    100.0*FLOAT(NPMAX)/FLOAT(NTOTAL), '%)'
c        CALL RFMLOG ( WRNMSG, FAIL, ERRMSG )
c        IF ( FAIL ) RETURN
c      END IF
C
      IF ( NPMIN .GT. 0 .AND. NPLUT(ILUT) .GT. 1 ) THEN
        WRITE ( WRNMSG, '(A,I6,A1,I7,A,F5.1,A2)' )
     &    'W-LUTRNG: Min -lnp edge of LUT exceeded for', NPMIN,
     &    '/', NTOTAL, ' calc paths (', 
     &    100.0*FLOAT(NPMIN)/FLOAT(NTOTAL), '%)'
        CALL RFMLOG ( WRNMSG, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
C
      IF ( NTMAX .GT. 0 .AND. NTLUT(ILUT) .GT. 1 ) THEN
        WRITE ( WRNMSG, '(A,I6,A1,I7,A,F5.1,A2)' )
     &    'W-LUTRNG: Max Temp edge of LUT exceeded for', NTMAX,
     &    '/', NTOTAL, ' calc paths (', 
     &    100.0*FLOAT(NTMAX)/FLOAT(NTOTAL), '%)'
        CALL RFMLOG ( WRNMSG, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
C
      IF ( NTMIN .GT. 0 .AND. NTLUT(ILUT) .GT. 1 ) THEN
        WRITE ( WRNMSG, '(A,I6,A1,I7,A,F5.1,A2)' )
     &    'W-LUTRNG: Min Temp edge of LUT exceeded for', NTMIN,
     &    '/', NTOTAL, ' calc paths (', 
     &    100.0*FLOAT(NTMIN)/FLOAT(NTOTAL), '%)'
        CALL RFMLOG ( WRNMSG, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
C
      END
