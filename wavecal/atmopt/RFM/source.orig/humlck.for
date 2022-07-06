C $Id: humlck.for,v 1.7 1997/03/03 17:17:13 adarby Exp $
      SUBROUTINE HUMLCK ( N, X, Y, V )
C
C VERSION
C      03-MAR-97  AD  Version 3.
C      01-OCT-96  AD  Version 2.
C      01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Calculate Humlicek complex prob.function for Voigt Line shape.
C
C     The Voigt lineshape formulation:
C
C                 g(X,Y) = S * g0 * K(X,Y)
C                 g0 = 1/Ad * SQRT(ln2/pi)
C                 X = (nu - nu0)/Ad *SQRT(ln2)
C                 Y = Al/Ad *SQRT(ln2)
C                 K(X,Y) = Y/pi * 
C                 INT^(+infty)_(-infty){exp(-t**2)/[Y**2 + (X-t)**2]}dt
C
C     This routine calculates the complex probability function using a 
C     vectorized version of the Humlicek JQSRT V27 437 1982 paper. 
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER          N             !  I  No. wavenumber points to be evaluted
      REAL             X(N)          !  I  X array 
      REAL             Y             !  I  Y value 
      COMPLEX          V(N)          !  O  Absorption 
C
C LOCAL VARIABLES
      INTEGER  I           
      REAL     S1V, S2V
      COMPLEX  U,T
C
C EXECUTABLE CODE --------------------------------------------------------------
C
C Sort the (x,y) pairs into the 4 regions of the Humlicek expressions of the 
C Voigt line profile.
C
      DO I = 1, N 
        S1V = ABS(X(I)) + Y
        S2V = ( 0.195 * ABS(X(I)) ) - 0.176
        T = CMPLX ( Y, -X(I) )
C
C Region 1 of Humlicek
C
        IF ( S1V .GE. 15.0 ) THEN
          V(I) = T * 0.5641896 / ( 0.5 + (T*T) )
C
C Region 2 OF Humlicek
C
        ELSEIF ( S1V .GE. 5.5 ) THEN
          U = T * T
          V(I) = T * ( 1.410474 + U * 0.5641896 ) / 
     &               ( 0.75 + U * ( 3.0 + U ) )
C
C Region 3 of Humlicek
C
        ELSE IF ( Y .GE. S2V ) THEN
          V(I) = ( 16.4955 + T * ( 20.20933 + T * ( 11.96482 + 
     &                T * ( 3.778987 + T * 0.5642236 ) ) ) ) /
     &              ( 16.4955 + T * ( 38.82363 + T * ( 39.27121 +
     &                T * ( 21.69274 + T * ( 6.699398 + T ) ) ) ) )
C
C Region 4 of Humlicek
C
        ELSE
          U = T * T
          V(I) = CEXP ( U ) - T * ( 36183.31 - U * ( 3321.9905 -
     &            U * ( 1540.787 - U * ( 219.0313 - U * ( 35.76683 - 
     &            U * ( 1.320522 - U * 0.56419 ) ) ) ) ) ) /
     &            ( 32066.6 - U * ( 24322.84 - U * ( 9022.228 - 
     &            U * ( 2186.181 - U * ( 364.2191 - U * ( 61.57037 - 
     &            U * ( 1.841439-U ) ) ) ) ) ) )
        END IF
      END DO
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
