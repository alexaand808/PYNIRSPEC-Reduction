C $Id: braket.for,v 1.8 1997/03/07 10:01:58 adarby Exp $
      SUBROUTINE BRAKET ( ARRAY, N, X, I )
C
C VERSION
C      03-MAR-97  AD  Version 3.
C      01-OCT-96  AD  Version 2.
C      01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Return lower index I of the two points in array ARRAY(N)
C     whose values bracket X.  
C     i.e. L:  ARRAY(I)>X>ARRAY(I+1)  or  ARRAY(I)<X<ARRAY(I+1)
C     Special cases:
C         If X<ARRAY(1)<ARRAY(N)  or  X>ARRAY(1)>ARRAY(N), then I = 0
C         If X<ARRAY(N)<ARRAY(1)  or  X>ARRAY(N)>ARRAY(1), then I = N
C    Routine does not check parameters or if the array is monotonic
C    Note: there is also a DP version of this subroutine: DBRAKT.FOR
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER N          !  I  Array length
      REAL    ARRAY(N)   !  I  Data Array
      REAL    X          !  I  Value to be bracketted
      INTEGER I          !  O  Index of array element below X
C
C LOCAL VARIABLES
      INTEGER J, K, M
      LOGICAL A                                    ! True if an increasing array
C
CEXECUTABLE CODE ---------------------------------------------------------------
C
      K = 1
      A = ARRAY(1) .LT. ARRAY(N)
      IF ( I.LE.0 .OR. I.GT.N ) I = N/2
C
      K = 1
      IF ( X.GE.ARRAY(I) .EQV. A ) THEN
    1   J = I + K
        IF ( J .GT. N ) THEN
          J = N + 1
        ELSE IF ( X.GE.ARRAY(J) .EQV. A ) THEN
          I = J
          K = K + K
          GO TO 1
        ENDIF
      ELSE
        J = I
    2   I = J - K
        IF ( I .LT. 1 ) THEN
          I = 0
        ELSE IF ( X.LT.ARRAY(I) .EQV. A ) THEN
          J = I
          K = K + K
          GO TO 2
        ENDIF
      ENDIF
C
    3 CONTINUE
      IF ( J-I .NE. 1 ) THEN
        M = (J+I)/2
        IF ( X.GT.ARRAY(M) .EQV. A ) THEN
          I = M
        ELSE
          J = M
        ENDIF
        GO TO 3
      ENDIF
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
