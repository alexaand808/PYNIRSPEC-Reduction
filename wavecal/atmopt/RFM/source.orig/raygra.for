      SUBROUTINE RAYGRA ( HGTINI, PSIINI, ZENINI, IFIN, S, 
     &                    HGTFIN, PSIFIN, ZENFIN )
C
C VERSION
C     17-DEC-01  AD  Initialise DFIN to avoid warning messages
C     23-MAR-01  AD  Add more intelligent adjustment of DS for convergence.
C     27-APR-00  AD  Remove redundant variables REFRAC and RFRIDX.
C     02-JAN-00  AD  Original.
C
C DESCRIPTION
C     Ray-tracing in a 2-D atmosphere (z,psi).
C
C     If IFIN=0, terminates after one step DS
C
      IMPLICIT NONE
C
C ARGUMENTS
      REAL     HGTINI  !  I  Initial altitude [km]
      REAL     PSIINI  !  I  Initial LOS angle [deg]
      REAL     ZENINI  !  I  Initial Zenith Angle [deg]
      INTEGER  IFIN    !  I  Termination parameter (1=HGT,2=PSI,3=ZEN)
      REAL     S       ! I/O Init Path increment (I) or Cumulative path (O)[km]
      REAL     HGTFIN  ! I/O Final height [km]
      REAL     PSIFIN  ! I/O Final LOS angle [deg]
      REAL     ZENFIN  ! I/O Final Zenith Angle [deg]
C
      EXTERNAL
     &  GRADVS ! Derivatives for 2-D ray-tracing.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'phycon.inc' ! Physical and Mathematical constants
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'crvcom.inc' ! Local radius of curvature
C
C LOCAL CONSTANTS
      DOUBLE PRECISION TTOL
        PARAMETER ( TTOL = 1.0D-6 )
      DOUBLE PRECISION PTOL
        PARAMETER ( PTOL = 1.0D-6 )
      DOUBLE PRECISION RTOL
        PARAMETER ( RTOL = 1.0D-3 )
C
C LOCAL VARIABLES
      LOGICAL        SURFCE ! T=stopped at surface, F=stopped at reqd.limit
      DOUBLE PRECISION DFIN ! Difference from target value of P,R,T
      DOUBLE PRECISION DPDS, DPDS1, DPDS2, DPDS3, DPDS4 ! d(P)/d(S)
      DOUBLE PRECISION DRDS, DRDS1, DRDS2, DRDS3, DRDS4 ! d(R)/d(S)
      DOUBLE PRECISION DTDS, DTDS1, DTDS2, DTDS3, DTDS4 ! d(T)/d(S)
      DOUBLE PRECISION DS   ! Path length increment [km]
      DOUBLE PRECISION DS2  ! DS/2
      DOUBLE PRECISION DSFC ! Distance above surface [km]
      DOUBLE PRECISION DSSUM ! Cumulatived path [km]
      DOUBLE PRECISION P    ! D.P. horizontal angle [rad]
      DOUBLE PRECISION PFIN ! Termination value of P (if IFIN=2)
      DOUBLE PRECISION PINT   ! P evaluated at intermediate points
      DOUBLE PRECISION R    ! Radius coordinate of path [km]
      DOUBLE PRECISION RFIN ! Termination value of R (if IFIN=1)
      DOUBLE PRECISION RINT   ! R evaluated at intermediate points
      DOUBLE PRECISION RSFC ! Radius of curvature at earth's surface
      DOUBLE PRECISION T    ! D.P. zenith angle [rad]
      DOUBLE PRECISION TFIN ! Termination value of T (if IFIN=3)
      DOUBLE PRECISION TINT   ! T evaluated at intermediate points
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      SURFCE = .FALSE.
      RSFC = DBLE ( HGTATM(1) ) + RADCRV 
      DS = DBLE ( S ) 
      DSSUM = 0.0D0
      DFIN =0.0D0            ! Not necessary but avoids warnings
C
C Look up atmospheric parameters for (z,psi)_INI
C
      R = RADCRV + DBLE ( HGTINI )
      P = DBLE ( PSIINI * DTORAD )
      T = DBLE ( ZENINI * DTORAD )
C
C Repeat from here for each path increment
  100 CONTINUE
      DSFC = R - RSFC 
      IF ( IFIN .NE. 0 .AND. ABS ( DSFC ) .LE. RTOL ) THEN
        SURFCE = .TRUE.
        GOTO 300
      ELSE IF ( IFIN .EQ. 1 ) THEN
        RFIN = RADCRV + DBLE ( HGTFIN )
        DFIN = R - RFIN 
        IF ( ABS ( DFIN ) .LE. RTOL ) GOTO 300
      ELSE IF ( IFIN .EQ. 2 ) THEN
        PFIN = DBLE ( PSIFIN * DTORAD )
        DFIN = P - PFIN 
        IF ( ABS ( DFIN ) .LE. PTOL ) GOTO 300
      ELSE IF ( IFIN .EQ. 3 ) THEN
        TFIN = DBLE ( ZENFIN * DTORAD )
        DFIN = T - TFIN 
        IF ( ABS ( DFIN ) .LE. TTOL ) GOTO 300
      END IF
C 
C Repeat from here if modifying step size
  200 CONTINUE
      DS2 = DS * 0.5D0
C Runge-Kutta 4th order numerical integration increments y according to 
C  y_(i+1) = y_i + ( dy/dx_1 + 2*dy/dx_2 + 2*dy/dx_3 + dy/dx_4 )*dx/6
C where derivatives dy/dx_i are evaluated at: 
C  dy/dx_1: y_i 
C  dy/dx_2: y_i + (dy/dx_1)*dx/2
C  dy/dx_3: y_i + (dy/dx_2)*dx/2 
C  dy/dx_4: y_i + (dy/dx_3)*dx
      CALL GRADVS ( P, R, T, DPDS1, DRDS1, DTDS1 )      
      PINT = P + DPDS1 * DS2
      RINT = R + DRDS1 * DS2
      TINT = T + DTDS1 * DS2
      CALL GRADVS ( PINT, RINT, TINT, DPDS2, DRDS2, DTDS2 )   
      PINT = P + DPDS2 * DS2              
      RINT = R + DRDS2 * DS2 
      TINT = T + DTDS2 * DS2 
      CALL GRADVS ( PINT, RINT, TINT, DPDS3, DRDS3, DTDS3 )   
      PINT = P + DPDS3 * DS
      RINT = R + DRDS3 * DS
      TINT = T + DTDS3 * DS
      CALL GRADVS ( PINT, RINT, TINT, DPDS4, DRDS4, DTDS4 )      
      DPDS = ( DPDS1 + DPDS4 + 2.0D0 * ( DPDS2 + DPDS3 ) ) / 6.0D0
      DRDS = ( DRDS1 + DRDS4 + 2.0D0 * ( DRDS2 + DRDS3 ) ) / 6.0D0
      DTDS = ( DTDS1 + DTDS4 + 2.0D0 * ( DTDS2 + DTDS3 ) ) / 6.0D0
C      
C Putting these lines in gives 1st order (Euler) integration
C	dpds = dpds1
C	drds = drds1
C	dtds = dtds1
C
C If moving further under surface or further away from termination criterion, 
C set smaller stepsize DS
      IF ( IFIN .NE. 0 .AND. DSFC + DRDS * DS .LT. - RTOL ) THEN
        DS = - DSFC / DRDS 
        GOTO 200
      ELSE IF ( IFIN .EQ. 1 .AND. 
     &  ABS ( DFIN + DRDS * DS ) .GT. ABS ( DFIN ) ) THEN
        DS = -DS * DFIN / ( DFIN + DRDS * DS )
        GOTO 200
      ELSE IF ( IFIN .EQ. 2 .AND. 
     &  ABS ( DFIN + DPDS * DS ) .GT. ABS ( DFIN ) ) THEN
        DS = -DS * DFIN / ( DFIN + DPDS * DS )
        GOTO 200
      ELSE IF ( IFIN .EQ. 3 .AND. 
     &  ABS ( DFIN + DTDS * DS ) .GT. ABS ( DFIN ) ) THEN
        DS = -DS * DFIN / ( DFIN + DTDS * DS )
        GOTO 200
      ELSE                               ! Add increment with same step size
        P = P + DPDS * DS
        R = R + DRDS * DS
        T = T + DTDS * DS
        DSSUM = DSSUM + DS
        IF ( IFIN .NE. 0 ) GOTO 100      ! IFIN=0 means one iteration only
      END IF
C
C Set final values to last iteration unless original term.criterion reached.
  300 CONTINUE

      IF ( IFIN .EQ. 0 ) THEN
        HGTFIN = SNGL ( R - RADCRV )
        PSIFIN = SNGL ( P ) / DTORAD 
        ZENFIN = SNGL ( T ) / DTORAD
      ELSE IF ( SURFCE ) THEN
        HGTFIN = HGTATM(1)
        PSIFIN = SNGL ( P ) / DTORAD 
        ZENFIN = SNGL ( T ) / DTORAD
      ELSE 
        IF ( IFIN .NE. 1 ) HGTFIN = SNGL ( R - RADCRV )
        IF ( IFIN .NE. 2 ) PSIFIN = SNGL ( P ) / DTORAD 
        IF ( IFIN .NE. 3 ) ZENFIN = SNGL ( T ) / DTORAD
      END IF
      S = SNGL ( DSSUM )
C
      END
