C $Id: lookup.for,v 1.9 1997/03/03 17:17:29 adarby Exp $
      SUBROUTINE LOOKUP ( X, Y, XTAB, YTAB, N, JGUESS, MORD )
C
C VERSION
C      03-MAR-97  AD  Version 3.
C      01-OCT-96  AD  Version 2.
C      01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     General purpose interpolation routine.
C     General purpose
C     Given arrays of N corresponding values in arrays XTAB and YTAB the routine
C     will return the interpolated value (Y) at the point X using a polynomial
C     of order MORD.  XTAB must be monotonic.
C     If the routine is called with JGUESS such that X is approximately 
C     X(JGUESS) the routine will work faster but JGUESS can take any value on
C     input.  On output JGUESS is set to the highest of the MORD element indices
C     of the input table that was used for the interpolation; this is a
C     value appropriate for making another entry in the same part of the table.
C     Note: also a DP version of this subroutine: DLOOKP.FOR
C
      IMPLICIT NONE
C
C ARGUMENTS
      REAL X          ! I  Point at which interpolation is required
      REAL Y          ! O  Answer
      REAL XTAB(*)    ! I  Tabulation of N x coordinates
      REAL YTAB(*)    ! I  Tabulation of N y coordinates
      INTEGER N       ! I  Number of points in table
      INTEGER JGUESS  !I/O Entry point to table
      INTEGER MORD    ! I  Order of interpolating polynomial
C
      EXTERNAL 
     &  BRAKET ! Return lower index I of the two points in array ARRAY(N)
C
C LOCAL VARIABLES
      INTEGER I, J, NN 
      REAL T
C
C- EXECUTABLE CODE -------------------------------------------------------------
C
      J = JGUESS
      CALL BRAKET ( XTAB, N, X, J ) 
      NN     = MAX ( 1, J-(MORD-1)/2 )
      JGUESS = MIN ( N, NN+MORD )
      NN     = JGUESS - MORD
C
      Y = 0.0
      DO J = NN, JGUESS
        T = YTAB(J)
        DO I = NN, JGUESS
         IF ( I .NE. J ) T = T * ( X - XTAB(I) ) / ( XTAB(J) - XTAB(I) )
        ENDDO
        Y = Y + T
      ENDDO
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
