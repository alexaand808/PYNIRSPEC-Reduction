      SUBROUTINE TABAXS ( ILUT, NP, P1, DP, NT, T1, DT, X )
C
C VERSION
C      04-MAY-99  AD  Original.
C 
C DESCRIPTION    
C     Resample (p,T) axis from TAB input file.
C     Called by RFMTAB for each record read from TAB file if qualifiers used
C     to resample (p,T) axes.
C     Notes: interpretation of resampling increments NDP, NDT
C        +ve multiply DP, DT by positive integer factor
C        -ve reverse axes (changing P1,T1) before multiplying DP,DT by 
C            negative integer factor
C        0   Take single point nearest middle of axes (changes P1,T1).
C                  
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER ILUT   !  I  Index of LUT in LUTCOM.INC
      INTEGER NP     !  O  Resampled number of points on -lnp axis
      REAL    P1     !  O  Resampled value of first -lnp point
      REAL    DP     !  O  Resampled -lnp axis increment
      INTEGER NT     !  O  Resampled number of points on -lnp axis
      REAL    T1     !  O  Resampled value of first Temp point
      REAL    DT     !  O  Resampled Temp axis increment
      REAL    X(*)   ! I/O Original/resampled combined (p,T) axis (DIM=NP*NT)
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'lutcom.inc' ! Look-Up Table data
C
C LOCAL VARIABLES
      INTEGER IP        ! DO loop counter for original -lnp axis
      INTEGER IP1, IP2  ! DO loop limits for original -lnp axis
      INTEGER IPSTEP    ! Step size for loop over original -lnp axis      
      INTEGER IT        ! DO loop counter for original Temp axis
      INTEGER IT1, IT2  ! DO loop limits for original Temp axis
      INTEGER ITSTEP    ! Step size for loop over original Temp axis
      INTEGER IX        ! Counter for original axis
      INTEGER IX0       ! Offset for IX indexing
      INTEGER IY        ! Counter for new axis
      REAL    Y(MAXLUX) ! Work space
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      IY = 0
      IF ( NDPLUT(ILUT) .GT. 0 ) THEN       
        IP1 = 1
        IP2 = NPLUT(ILUT)
        IPSTEP = NDPLUT(ILUT)
        P1  = P1LUT(ILUT)
        DP  = DPLUT(ILUT) * NDPLUT(ILUT)
        NP = 1 + ( NPLUT(ILUT) - 1 ) / ABS( NDPLUT(ILUT) )
      ELSE IF ( NDPLUT(ILUT) .LT. 0 ) THEN
        IP1 = NPLUT(ILUT)
        IP2 = 1
        IPSTEP = NDPLUT(ILUT)
        P1  = P1LUT(ILUT) + ( NPLUT(ILUT) - 1 ) * DPLUT(ILUT)
        DP  = DPLUT(ILUT) * NDPLUT(ILUT)
        NP = 1 + ( NPLUT(ILUT) - 1 ) / ABS( NDPLUT(ILUT) )
      ELSE
        IP1 = 1 + INT ( ( NPLUT(ILUT) - 1 ) / 2 ) ! NP=1,2,3,4,5 :IP1=1,1,2,2,3
        IP2 = IP1
        IPSTEP = 1
        P1  = P1LUT(ILUT) + ( IP1 - 1 ) * DPLUT(ILUT)
        DP  = DPLUT(ILUT)
        NP  = 1
      END IF
      IF ( NDTLUT(ILUT) .GT. 0 ) THEN
        IT1 = 1
        IT2 = NTLUT(ILUT)
        ITSTEP = NDTLUT(ILUT)
        T1  = T1LUT(ILUT)
        DT  = DTLUT(ILUT) * NDTLUT(ILUT)
        NT = 1 + ( NTLUT(ILUT) - 1 ) / ABS( NDTLUT(ILUT) )
      ELSE IF ( NDTLUT(ILUT) .LT. 0 ) THEN
        IT1 = NTLUT(ILUT)
        IT2 = 1
        ITSTEP = NDTLUT(ILUT)
        T1  = T1LUT(ILUT) + ( NTLUT(ILUT) - 1 ) * DTLUT(ILUT)
        DT  = DTLUT(ILUT) * NDTLUT(ILUT)
        NT = 1 + ( NTLUT(ILUT) - 1 ) / ABS( NDTLUT(ILUT) )
      ELSE
        IT1 = 1 + INT ( ( NTLUT(ILUT) - 1 ) / 2 )  
        IT2 = IT1
        ITSTEP = 1
        T1  = T1LUT(ILUT) + ( IT1 - 1 ) * DTLUT(ILUT)
        DT  = DTLUT(ILUT) 
        NT  = 1
      END IF
C
      IY = 0
      DO IT = IT1, IT2, ITSTEP
        IX0 = NPLUT(ILUT) * ( IT - 1 ) 
        DO IP = IP1, IP2, IPSTEP
          IX = IX0 + IP
          IY = IY + 1
          Y(IY) = X(IX)
        END DO
      END DO
C
      DO IY = 1, NP * NT
        X(IY) = Y(IY)
      END DO
C
      END
