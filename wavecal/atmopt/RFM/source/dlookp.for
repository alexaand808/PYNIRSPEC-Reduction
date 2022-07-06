      SUBROUTINE DLOOKP ( X, Y, XTAB, YTAB, N, IGUESS )
C
C VERSION
C     24-APR-12  AD  Remove MORD argument, and simplify internal logic.
C     27-APR-00  AD  Original. DP version of LOOKUP.
C
C DESCRIPTION
C     General purpose DP interpolation routine.
C     Last modification limits Y to end-values of YTAB.
C     Given arrays of N corresponding values in arrays XTAB and YTAB the routine
C     will return the linearly interpolated value (Y) at the point X.
C     XTAB must be monotonic.
C     If the routine is called with JGUESS such that X is approximately 
C     X(JGUESS) the routine will work faster but JGUESS can take any value on
C     input.  On output JGUESS is set to the lower of the two indices
C     of the input table that was used for the interpolation; this is a
C     value appropriate for making another entry in the same part of the table.
C
      IMPLICIT NONE
C
C ARGUMENTS
      DOUBLE PRECISION X       !  I  Point at which interpolation is required
      DOUBLE PRECISION Y       !  O  Answer
      DOUBLE PRECISION XTAB(*) !  I  Tabulation of N x coordinates
      DOUBLE PRECISION YTAB(*) !  I  Tabulation of N y coordinates
      INTEGER          N       !  I  Number of points in table
      INTEGER          IGUESS  ! I/O Entry point to table
C
      EXTERNAL 
     &  DBRAKT ! Return lower index I of the two points in array ARRAY(N)
C
C LOCAL VARIABLES
      INTEGER I
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      I = IGUESS
      CALL DBRAKT ( XTAB, N, X, I ) 
      IF ( I .EQ. 0 ) THEN
        Y = YTAB(1)
        IGUESS = 1
      ELSE IF ( I .EQ. N ) THEN
        Y = YTAB(N)
        IGUESS = N
      ELSE 
        Y = YTAB(I) + ( X - XTAB(I) ) / ( XTAB(I+1) - XTAB(I) ) *
     &      ( YTAB(I+1) - YTAB(I) )
        IGUESS = I
      END IF
C
      END
