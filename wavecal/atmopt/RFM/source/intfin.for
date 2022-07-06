      SUBROUTINE INTFIN ( IFIN, II, IJ, IK, IL )
C
C VERSION
C     04-MAY-01  AD  Correction: loop to LTAN, not NTAN
C     26-MAY-00  AD  Convert RADTAN and ABSTAN from S.P. to D.P.
C     11-AUG-99  AD  Remove redundant fincom.inc
C     05-AUG-99  AD  Use DBLE ( ) in calculations of TRATAN 
C     19-DEC-97  AD  Original.
C
C DESCRIPTION    
C     Apply irregular grid interpolation function.
C     Called by RFMINT for each finemesh interval to be interpolated.
C     This module applies interpolation over all tangent heights for the 
C     arrays RADTAN (radiance), TRATAN (transmission) and ABSTAN (opt.depth)
C     containing the path-integrated values at the calculated grid points.
C     Four indices are input (II,IJ,IK,IL), where IJ,IK are the indices of the
C     calculated values at wavenumbers below and above the required 
C     interpolation point (IFIN), and (II,IL) are the next points below and
C     above - although these are not used for the two-point interpolations.
C     This module is not called for IFIN equal to any of the calc.points.
C                  
      IMPLICIT NONE
C
C ARGUMENTS      
      INTEGER IFIN !  I  Index of required interpolation point
      INTEGER II   !  I  Index of calc.value two below required point
      INTEGER IJ   !  I  Index of calc.value below required point
      INTEGER IK   !  I  Index of calc.value above required point
      INTEGER IL   !  I  Index of calc.value two above required point
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'grdcom.inc' ! Irregular grid
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL CONSTANTS
      DOUBLE PRECISION SMALL  ! Min.val of rad,abs,tra for inverse or log fns.
        PARAMETER ( SMALL = 1.0D-32 )
C
C LOCAL VARIABLES
      INTEGER I,J,K,L ! Local values of II,IJ,IK,IL
      INTEGER M    ! Either II or IL for 3-point interpolations
      INTEGER ITAN ! Tangent path counter
      INTEGER N    ! Local value of IFIN
      DOUBLE PRECISION AI  ! Weight for point two below required point
      DOUBLE PRECISION AJ  ! Weight for point below required point
      DOUBLE PRECISION AK  ! Weight for point above required point
      DOUBLE PRECISION AL  ! Weight for point two above required point
      DOUBLE PRECISION AM  ! Weight for point two above or below reqd point
      DOUBLE PRECISION Y   ! Intermediate value of required function
C
C EXECUTABLE CODE -------------------------------------------------------------
C 
      I = II
      J = IJ
      K = IK
      L = IL
      N = IFIN
      IF ( TWOGRD ) THEN
        AJ = DBLE ( N - K ) / DBLE ( J - K )            
        AK = DBLE ( N - J ) / DBLE ( K - J )
        IF ( FNCGRD .EQ. 'lin' ) THEN                      ! ax + b
          DO ITAN = 1, LTAN
            RADTAN(N,ITAN) = AJ*RADTAN(J,ITAN) + AK*RADTAN(K,ITAN)
            TRATAN(N,ITAN) = AJ*TRATAN(J,ITAN) + AK*TRATAN(K,ITAN)
            ABSTAN(N,ITAN) = AJ*ABSTAN(J,ITAN) + AK*ABSTAN(K,ITAN)
          END DO
        ELSE IF ( FNCGRD .EQ. '1li' ) THEN                 ! 1/(ax+b)
          DO ITAN = 1, LTAN
            Y = AJ / MAX ( SMALL, RADTAN(J,ITAN) ) + 
     &          AK / MAX ( SMALL, RADTAN(K,ITAN) )
            RADTAN(N,ITAN) = 1.0D0 / Y
            Y = AJ / MAX ( SMALL, TRATAN(J,ITAN) ) + 
     &          AK / MAX ( SMALL, TRATAN(K,ITAN) )
            TRATAN(N,ITAN) = 1.0D0 / Y
            Y = AJ / MAX ( SMALL, ABSTAN(J,ITAN) ) + 
     &          AK / MAX ( SMALL, ABSTAN(K,ITAN) )
            ABSTAN(N,ITAN) = 1.0D0 / Y
          END DO
        ELSE IF ( FNCGRD .EQ. '1sq' ) THEN                 ! 1/(ax+b)^2
          DO ITAN = 1, LTAN
            Y = AJ / SQRT ( MAX ( SMALL, RADTAN(J,ITAN) ) ) + 
     &          AK / SQRT ( MAX ( SMALL, RADTAN(K,ITAN) ) )
            RADTAN(N,ITAN) = 1.0D0 / MAX ( SMALL, Y**2 )
            Y = AJ / SQRT ( MAX ( SMALL, TRATAN(J,ITAN) ) ) + 
     &          AK / SQRT ( MAX ( SMALL, TRATAN(K,ITAN) ) )
            TRATAN(N,ITAN) = 1.0D0 / MAX ( SMALL, Y**2 )
            Y = AJ / SQRT ( MAX ( SMALL, ABSTAN(J,ITAN) ) ) + 
     &          AK / SQRT ( MAX ( SMALL, ABSTAN(K,ITAN) ) )
            ABSTAN(N,ITAN) = 1.0D0 / MAX ( SMALL, Y**2 )
          END DO
        ELSE IF ( FNCGRD .EQ. 'lor' ) THEN                 ! a/(x^2+b)
          AJ = DBLE ( (K-N)*(K+N) ) / DBLE ( (K-J)*(K+J) )
          AK = DBLE ( (J-N)*(J+N) ) / DBLE ( (J-K)*(J+K) )          
          DO ITAN = 1, LTAN
            Y = AJ / MAX ( SMALL, RADTAN(J,ITAN) ) + 
     &          AK / MAX ( SMALL, RADTAN(K,ITAN) )
            RADTAN(N,ITAN) = 1.0D0 / MAX ( SMALL, Y )
            Y = AJ / MAX ( SMALL, TRATAN(J,ITAN) ) + 
     &          AK / MAX ( SMALL, TRATAN(K,ITAN) )
            TRATAN(N,ITAN) = 1.0D0 / MAX ( SMALL, Y )
            Y = AJ / MAX ( SMALL, ABSTAN(J,ITAN) ) + 
     &          AK / MAX ( SMALL, ABSTAN(K,ITAN) )
            ABSTAN(N,ITAN) = 1.0D0 / MAX ( SMALL, Y )
          END DO
        ELSE IF ( FNCGRD .EQ. 'lnl' ) THEN                 ! exp(ax+b)
          DO ITAN = 1, LTAN
            Y = AJ * LOG ( MAX ( SMALL, RADTAN(J,ITAN) ) ) +
     &          AK * LOG ( MAX ( SMALL, RADTAN(K,ITAN) ) ) 
            RADTAN(N,ITAN) = EXP ( Y )
            Y = AJ * LOG ( MAX ( SMALL, TRATAN(J,ITAN) ) ) +
     &          AK * LOG ( MAX ( SMALL, TRATAN(K,ITAN) ) ) 
            TRATAN(N,ITAN) = EXP ( Y )
            Y = AJ * LOG ( MAX ( SMALL, ABSTAN(J,ITAN) ) ) +
     &          AK * LOG ( MAX ( SMALL, ABSTAN(K,ITAN) ) ) 
            ABSTAN(N,ITAN) = EXP ( Y )
          END DO
        ELSE
          STOP 'F-INTFIN: Unrecognised function (TWOGRD=T)'
        END IF
      ELSE IF ( FNCGRD .EQ. 'qad' .OR. FNCGRD .EQ. '1qa' ) THEN  
        IF ( N - I .LT. L - N ) THEN
          M = I
        ELSE
          M = L
        END IF
        AM = DBLE (           (N-J)*(N-K) ) /
     &       DBLE (             (M-J)*  (M-K) )
        AJ = DBLE ( (N-M)          *(N-K) ) /
     &       DBLE (   (J-M)          *  (J-K) )
        AK = DBLE ( (N-M)*(N-J)           ) /
     &       DBLE (   (K-M)*  (K-J)           )
        IF ( FNCGRD .EQ. 'qad' ) THEN                ! ax^2 + bx + c
          DO ITAN = 1, LTAN
            RADTAN(N,ITAN) =  AJ*RADTAN(J,ITAN) + AK*RADTAN(K,ITAN)
     &                         + AM*RADTAN(M,ITAN)          
            TRATAN(N,ITAN) =  AJ*TRATAN(J,ITAN) + AK*TRATAN(K,ITAN)
     &                         + AM*TRATAN(M,ITAN)          
            ABSTAN(N,ITAN) =  AJ*ABSTAN(J,ITAN) + AK*ABSTAN(K,ITAN)
     &                         + AM*ABSTAN(M,ITAN)          
          END DO
        ELSE                                  ! '1qa'  1/(ax^2 + bx + c)
          DO ITAN = 1, LTAN
            Y = AJ / MAX ( SMALL, RADTAN(J,ITAN) ) + 
     &          AK / MAX ( SMALL, RADTAN(K,ITAN) ) +
     &          AM / MAX ( SMALL, RADTAN(M,ITAN) ) 
            RADTAN(N,ITAN) = 1.0D0 / MAX ( Y, SMALL )
            Y = AJ / MAX ( SMALL, TRATAN(J,ITAN) ) + 
     &          AK / MAX ( SMALL, TRATAN(K,ITAN) ) +
     &          AM / MAX ( SMALL, TRATAN(M,ITAN) ) 
            TRATAN(N,ITAN) = 1.0D0 / MAX ( Y, SMALL )
            Y = AJ / MAX ( SMALL, ABSTAN(J,ITAN) ) + 
     &          AK / MAX ( SMALL, ABSTAN(K,ITAN) ) +
     &          AM / MAX ( SMALL, ABSTAN(M,ITAN) ) 
            ABSTAN(N,ITAN) = 1.0D0 / MAX ( Y, SMALL )
          END DO
        END IF
      ELSE
        AI = DBLE (           (N-J)*(N-K)*(N-L) ) /
     &       DBLE (             (I-J)*  (I-K)*  (I-L) )
        AJ = DBLE ( (N-I)          *(N-K)*(N-L) ) /
     &       DBLE (   (J-I)          *  (J-K)*  (J-L) )
        AK = DBLE ( (N-I)*(N-J)          *(N-L) ) /
     &       DBLE (   (K-I)*  (K-J)          *  (K-L) )
        AL = DBLE ( (N-I)*(N-J)*(N-K) ) /
     &       DBLE (   (L-I)*  (L-J)*  (L-K) )
        IF ( FNCGRD .EQ. 'cub' ) THEN
          DO ITAN = 1, LTAN
            RADTAN(N,ITAN) = AI*RADTAN(I,ITAN) + AJ*RADTAN(J,ITAN) 
     &                        + AK*RADTAN(K,ITAN) + AL*RADTAN(L,ITAN)
            TRATAN(N,ITAN) = AI*TRATAN(I,ITAN) + AJ*TRATAN(J,ITAN) 
     &                        + AK*TRATAN(K,ITAN) + AL*TRATAN(L,ITAN)
            ABSTAN(N,ITAN) = AI*ABSTAN(I,ITAN) + AJ*ABSTAN(J,ITAN) 
     &                        + AK*ABSTAN(K,ITAN) + AL*ABSTAN(L,ITAN)
          END DO
        ELSE IF ( FNCGRD .EQ. '1cu' ) THEN
          DO ITAN = 1, LTAN
            Y = AI / MAX ( SMALL, RADTAN(I,ITAN) ) +
     &          AJ / MAX ( SMALL, RADTAN(J,ITAN) ) +
     &          AK / MAX ( SMALL, RADTAN(K,ITAN) ) +
     &          AL / MAX ( SMALL, RADTAN(L,ITAN) ) 
            RADTAN(N,ITAN) = 1.0D0 / Y
            Y = AI / MAX ( SMALL, TRATAN(I,ITAN) ) +
     &          AJ / MAX ( SMALL, TRATAN(J,ITAN) ) +
     &          AK / MAX ( SMALL, TRATAN(K,ITAN) ) +
     &          AL / MAX ( SMALL, TRATAN(L,ITAN) ) 
            TRATAN(N,ITAN) = 1.0D0 / Y
            Y = AI / MAX ( SMALL, ABSTAN(I,ITAN) ) +
     &          AJ / MAX ( SMALL, ABSTAN(J,ITAN) ) +
     &          AK / MAX ( SMALL, ABSTAN(K,ITAN) ) +
     &          AL / MAX ( SMALL, ABSTAN(L,ITAN) ) 
            ABSTAN(N,ITAN) = 1.0D0 / Y
          END DO
        ELSE IF ( FNCGRD .EQ. 'lnc' ) THEN
          DO ITAN = 1, LTAN
            Y = AI * LOG ( MAX ( SMALL, RADTAN(I,ITAN) ) ) +
     &          AJ * LOG ( MAX ( SMALL, RADTAN(J,ITAN) ) ) +
     &          AK * LOG ( MAX ( SMALL, RADTAN(K,ITAN) ) ) +
     &          AL * LOG ( MAX ( SMALL, RADTAN(L,ITAN) ) ) 
            RADTAN(N,ITAN) = EXP ( Y ) 
            Y = AI * LOG ( MAX ( SMALL, TRATAN(I,ITAN) ) ) +
     &          AJ * LOG ( MAX ( SMALL, TRATAN(J,ITAN) ) ) +
     &          AK * LOG ( MAX ( SMALL, TRATAN(K,ITAN) ) ) +
     &          AL * LOG ( MAX ( SMALL, TRATAN(L,ITAN) ) ) 
            TRATAN(N,ITAN) = EXP ( Y ) 
            Y = AI * LOG ( MAX ( SMALL, ABSTAN(I,ITAN) ) ) +
     &          AJ * LOG ( MAX ( SMALL, ABSTAN(J,ITAN) ) ) +
     &          AK * LOG ( MAX ( SMALL, ABSTAN(K,ITAN) ) ) +
     &          AL * LOG ( MAX ( SMALL, ABSTAN(L,ITAN) ) ) 
            ABSTAN(N,ITAN) = EXP ( Y ) 
          END DO
        ELSE
          STOP 'F-INTFIN: Unrecognised function (TWOGRD=F)'
        END IF
      END IF
C
      END 
