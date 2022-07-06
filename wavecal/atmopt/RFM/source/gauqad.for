      SUBROUTINE GAUQAD ( N, X, W )
C
C VERSION
C     16-APR-00  AD  Original.
C
C DESCRIPTION
C     Values and Weights for Gaussian First Moment Quadrature.
C     General purpose module.
C     Suitable for integrating x.f(x).dx in interval x=0:1
C     Numbers taken from Abramowitz and Stegun, p.921.
C     Values NQAD=1:4 are identical to those listed in Clough et al 1992.
C     Digits here have been checked!
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER          N    !  I  Degree of quadrature required (1:5)
      DOUBLE PRECISION X(*) !  O  List of N x-values (abscissa) (0:1)
      DOUBLE PRECISION W(*) !  O  List of N weights (0:1)
C
C NB: arrays X,W are dimensioned as (*) rather than (N) to avoid crashing the
C     code if N=0 or some negative value.
C
C REFERENCES
C     Handbook of Mathematical Functions 
C     M.Abramowitz and I.A.Stegun (Eds)
C     9th Dover printing, New York, 1972.
C
C     Clough et al, J.Geophys.Res. 97, 157610-15785 (1992).
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      IF ( N .EQ. 1 ) THEN
        X(1) = 0.6666666667D0
        W(1) = 0.5D0
      ELSE IF ( N .EQ. 2 ) THEN
        X(1) = 0.3550510257D0
        X(2) = 0.8449489743D0
        W(1) = 0.1819586183D0
        W(2) = 0.3180413817D0
      ELSE IF ( N .EQ. 3 ) THEN
        X(1) = 0.2123405382D0
        X(2) = 0.5905331356D0
        X(3) = 0.9114120405D0
        W(1) = 0.0698269799D0
        W(2) = 0.2292411064D0
        W(3) = 0.2009319137D0
      ELSE IF ( N .EQ. 4 ) THEN
        X(1) = 0.1397598643D0
        X(2) = 0.4164095676D0
        X(3) = 0.7231569864D0
        X(4) = 0.9428958039D0
        W(1) = 0.0311809710D0
        W(2) = 0.1298475476D0
        W(3) = 0.2034645680D0
        W(4) = 0.1355069134D0
      ELSE IF ( N .EQ. 5 ) THEN
        X(1) = 0.0985350858D0
        X(2) = 0.3045357266D0
        X(3) = 0.5620251898D0
        X(4) = 0.8019865821D0
        X(5) = 0.9601901429D0
        W(1) = 0.0157479145D0
        W(2) = 0.0739088701D0
        W(3) = 0.1463869871D0
        W(4) = 0.1671746381D0
        W(5) = 0.0967815902D0
      ELSE
        STOP 'F-GAUQAD: Argument N out of range'
      END IF
C
      END
