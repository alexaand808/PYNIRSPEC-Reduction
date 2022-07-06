      SUBROUTINE TRIANG ( XO, YO, N, X, Y, INIT, 
     &                    IFIT, WFIT, XFIT, YFIT, FAIL, ERRMSG ) 
C
C VERSION
C     03-MAY-05  AD  Correction: Save NTRI
C     17-DEC-01  AD  Save FINITE
C     03-SEP-99  AD  Original.
C
C DESCRIPTION
C     2-D interpolation of irregular grid using triangulation.
C
C     A plane is fitted to the three triangulation points and the interpolated
C     value on this plane is expressed as a set of weights applied to the 
C     value of the function at the triangulation points.
C     If the interpolation point lies outside the triangulated field, a value
C     is returned from the closest point along the boundary. In this case 
C     the arguments XFIT, YFIT will be different from the input coordinates 
C     XO,YO.
C
C     NB: the output indices of the interpolated points may be repeated 
C     (eg if only one point is chosen, the index will be repeated three times
C     and the corresponding weights will be 1.0, 0.0 and 0.0).
C
      IMPLICIT NONE
C
C ARGUMENTS
      REAL    XO          !  I  x-coordinate of interpolated point
      REAL    YO          !  I  y-coordinate of interpolated point
      INTEGER N           !  I  No. of given points in field (must be .GE. 1) 
      REAL    X(*)        !  I  x-coordinates of points in field, dim (N)
      REAL    Y(*)        !  I  y-coordinates of points in field, dim (N)
      LOGICAL INIT        ! I/O T=Initialise new grid (set F on exit)
      INTEGER IFIT(3)     !  O  Indices of points used for interpolation
      REAL    WFIT(3)     !  O  Weights of points used for interpolation 
      REAL    XFIT        !  O  x-coordinate of point actually fitted 
      REAL    YFIT        !  O  y-coordinate of point actually fitted
      LOGICAL FAIL        !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG !  O  Error message written if FAIL is TRUE 
C
C LOCAL CONSTANTS
      INTEGER       MAXN   !  Maximum expected size of N (for local arrays)
        PARAMETER ( MAXN = 1000 )
      REAL          SMALL  !  Minimum significant difference in cosine. 
        PARAMETER ( SMALL = 1.0E-5 ) 
C
C LOCAL VARIABLES
      INTEGER I,J,K,L      ! Indices of original points X,Y
      INTEGER II,JJ,KK,LL  ! Saved values of indices
      INTEGER IDXPER(MAXN) ! Index of perimeter points
      INTEGER IDXTRI(MAXN,3) ! Indices of triangle vertices
      INTEGER ISID         ! Side counter (of triangle)
      INTEGER ITRI,JTRI,KTRI,LTRI ! Triangle counters
      INTEGER IXMAX        ! (Original) Index of maximum value of x-coordinate
      INTEGER NPER         ! Number of perimeter points
      INTEGER NSAVE        ! Saved value of N from previous call 
      INTEGER NTRI         ! Number of triangles constructed
      REAL    CIJK,CIJL,CJKL,CKIL ! Cross products of vectors
      REAL    CIJO,CJKO,CKIO ! Cross products of vectors
      REAL    CIJTOL, CJKTOL, CKITOL ! Cross-product tolerances
      REAL    CKLI,CLIJ    ! Cross products of vectors 
      REAL    COSIJK, COSKLI ! Cosine of angles IJK and KLI
      REAL    COSINE       ! Cosine of angle between two vectors
      REAL    COSMAX       ! Maximum value of cosine
      REAL    DIJK, DKLI   ! Dot products of vectors
      REAL    DIJO, DJIO   ! Dot products of vectors
      REAL    DLIJ, DLJK, DLKI  ! Dot products of vectors
      REAL    FACTX,FACTY,FACTO ! Factors in computation of triangulation wgts
      REAL    FLIP         ! Benefit (+ve reduced distance) from diagonal flip
      REAL    FLIPMX       ! Maximum value of FLIP found so far
      REAL    R            ! Distance between points
      REAL    RMIN         ! Minimum distance from interp.point 
      REAL    RTOL         ! Tolerance in comparing R
      REAL    S            ! Length of vector
      REAL    SMIN         ! Minimum length of vector
      REAL    SIJ, SJK, SKL, SLI ! Distances between points
      REAL    SJP2         ! Square of length of vector JP 
      REAL    SINE         ! Sine of angle
      REAL    UIJ,UIL,UKI  ! x-coordinate of unit vectors along perimeter
      REAL    VIJ,VIL,VKI  ! y-coordinate of unit vectors along perimeter
      REAL    XMAX         ! Maximum value of x-coordinate
      LOGICAL FINITE       ! T=space has finite area, F=straight line
      LOGICAL INSIDE       ! T=inside triangle, F=outside
      LOGICAL PERIM(MAXN)  ! T=perimeter point, F=internal point
C
C DATA STATEMENTS
      DATA NSAVE  / 0 / 
      SAVE NSAVE, IDXTRI, FINITE, NTRI
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C First check that argument N is within valid range  (1:MAXN)
      FAIL = .TRUE.
      IF ( N .LE. 0 ) THEN
        WRITE ( ERRMSG, '(A,I11)' ) 
     &    'F-TRIANG: invalid argument N (ie not >0), value =', N
        RETURN
      ELSE IF ( N .GT. MAXN ) THEN
        WRITE ( ERRMSG, '(A,I11,A,I11)' ) 
     &    'F-TRIANG: No.points (N)=', N, 
     &    ' > Local array dimension MAXN=', MAXN
        RETURN
      ELSE IF ( .NOT. INIT .AND. N .NE. NSAVE ) THEN
        ERRMSG ='F-TRIANG: Not initialised for new X,Y array'
        RETURN
      END IF
      FAIL = .FALSE.
      NSAVE = N
C
C Deal with special cases of N=1 and N=2
      IF ( N .EQ. 1 ) THEN  ! Simply take value of single grid point
        IFIT(1) = 1
        IFIT(2) = 1
        IFIT(3) = 1
        WFIT(1) = 1.0
        WFIT(2) = 0.0
        WFIT(3) = 0.0
        XFIT    = X(1)
        YFIT    = Y(1)
        INIT    = .FALSE.
        RETURN
      ELSE IF ( N .EQ. 2 ) THEN ! Interpolate normally from line IJ
        I = 1
        J = 2
        IFIT(3) = 1
        WFIT(3) = 0.0
        DIJO = (X(I)-X(J)) * (XO-X(J)) + (Y(I)-Y(J)) * (YO-Y(J))
        DJIO = (X(J)-X(I)) * (XO-X(I)) + (Y(J)-Y(I)) * (YO-Y(I))
        IF ( DIJO .GT. 0.0 .AND. DJIO .GT. 0.0 ) THEN  ! Interpolate
          SJP2 = DIJO**2 / ( (X(J)-X(I))**2 + (Y(J)-Y(I))**2 )
          IFIT(1) = I
          IFIT(2) = J
          WFIT(1) = SJP2 / DIJO
          WFIT(2) = 1.0 - WFIT(1)
          XFIT = X(I)*WFIT(1) + X(J)*WFIT(2)
          YFIT = Y(I)*WFIT(1) + Y(J)*WFIT(2)
        ELSE                                 ! Duplicate closest point
          IFIT(2) = 1
          WFIT(2) = 0.0
          WFIT(1) = 1.0
          IF ( (X(I)-XO)**2 + (Y(I)-YO)**2 .LT.
     &         (X(J)-XO)**2 + (Y(J)-YO)**2      ) THEN  ! I is closest point
            IFIT(1) = I
            XFIT    = X(I)
            YFIT    = Y(I)
          ELSE                                          ! J is closest point
            IFIT(1) = J
            XFIT    = X(J)
            YFIT    = Y(J)
          END IF
        END IF
        INIT = .FALSE.
        RETURN
      END IF            
C
C Continue from here for N > 2, ie if triangulation possible
      IF ( .NOT. INIT ) GOTO 300
      INIT = .FALSE.
C
C Initialise for new array X,Y
C Step#1: establish point with the maximum x-coordinate 
      XMAX = X(1)
      IXMAX = 1
      PERIM(1) = .FALSE.
      DO I = 2, N
        PERIM(I) = .FALSE.
        IF ( X(I) .GT. XMAX ) THEN
          IXMAX = I
          XMAX = X(I)
        END IF
      END DO
C
C Step#2: Determine the array of points that form the perimeter
C Given a known point along the perimeter I ...
C for each unit vector KI along the perimeter (leading to point (I)), we want 
C to find the next point (J) whose unit vector IJ lies at the smallest angle 
C relative to KI (because we are proceeding anti-c/w around the perimeter, all
C angles are between 0 and 180 degrees by definition). 
C Cos(angle) between vectors is a suitable parameter to maximise since this
C decreases monotonically from 1 (0 degrees) to -1 (180 degrees).
C
      NPER = 0
      PERIM(IXMAX) = .TRUE.      ! Prevent IXMAX selected first time
      I = IXMAX       ! Start at a known perimeter point, (X(IXMAX),Y(IXMAX))
      UKI = 0.0       ! Start with unit vector KI in +y direction
      VKI = 1.0  
      J = 0
      DO WHILE ( J .NE. IXMAX )        ! Stop when returns to point IXMAX
        RMIN = 0.0
        COSMAX = -2.0                  ! Initialise search for next perim.point
        DO L = 1, N                    ! Try all points in field (L)
          IF ( .NOT. PERIM(L) ) THEN    ! avoiding previous perimeter points
            R  = SQRT ( (X(L)-X(I))**2 + (Y(L)-Y(I))**2 )
            UIL = ( X(L) - X(I) ) / R       ! Unit vector in dir. (I) to (J)
            VIL = ( Y(L) - Y(I) ) / R
            COSINE = UKI * UIL + VKI * VIL  ! Dot-product of KI.IL = cos(angle)
            IF ( COSINE .GT. COSMAX+SMALL ) THEN  ! New direction
              RMIN   = R
              J      = L
              UIJ    = UIL
              VIJ    = VIL
              COSMAX = COSINE
            ELSE IF ( COSINE .GT. COSMAX-SMALL .AND.        ! Same direction
     &                   R .LT. RMIN               ) THEN ! But closer
              J = L
              RMIN = R
            END IF
          END IF
        END DO   
        NPER = NPER + 1                ! Add point I to perimeter list
        IDXPER(NPER) = J
        PERIM(J) = .TRUE.
        I   = J                        ! Set current perimeter point I
        UKI = UIJ                      ! Unit vector along last perim.section
        VKI = VIJ
        PERIM(IXMAX) = .FALSE.         ! Allow point IXMAX to be selected 
      END DO
      PERIM(IXMAX) = .TRUE.
C
C Step#3: Form triangles by drawing lines from IXMAX to every other perimeter
C point. Since there are NPER perimeter points, there must be NPER-2 triangles.
C Also check if perimeter section subtends a finite angle at IXMAX: if it 
C does then the grid has a finite area (otherwise it is a straight line)
      NTRI = NPER - 2
      FINITE = .FALSE.
      DO ITRI = 1, NTRI
        I = IXMAX
        J = IDXPER(ITRI)
        K = IDXPER(ITRI+1)
        IDXTRI(ITRI,1) = I
        IDXTRI(ITRI,2) = J
        IDXTRI(ITRI,3) = K
        IF ( .NOT. FINITE ) THEN
          COSINE = ((X(J)-X(I))*(X(K)-X(I)) + (Y(J)-Y(I))*(Y(K)-Y(I)))
     &             / SQRT ( (X(J)-X(I))**2 + (Y(J)-Y(I))**2 ) 
     &             / SQRT ( (X(K)-X(I))**2 + (Y(K)-Y(I))**2 )
          FINITE = ( 1.0 - ABS ( COSINE ) ) .GT. SMALL
        END IF
      END DO
C
C Step#4: For each non-perimeter point L, find which existing triangle IJK it 
C lies inside, and sub-divide that triangle into three smaller triangles by 
C drawing lines from point to each vertex: IJL, JKL, KIL
      DO L = 1, N
        IF ( .NOT. PERIM(L) ) THEN
          DO ITRI = 1, NTRI
            I = IDXTRI(ITRI,1)
            J = IDXTRI(ITRI,2)
            K = IDXTRI(ITRI,3)
C If the point L is enclosed, then the cross-products of the vectors formed
C by each of the sides in turn (IJ, JK, KI) and the lines from each vertex to L
C (JL, KL, IL) must all be of the same sign. 
            CIJL = (X(J)-X(I))*(Y(L)-Y(J)) - (Y(J)-Y(I))*(X(L)-X(J))
            CJKL = (X(K)-X(J))*(Y(L)-Y(K)) - (Y(K)-Y(J))*(X(L)-X(K))
            CKIL = (X(I)-X(K))*(Y(L)-Y(I)) - (Y(I)-Y(K))*(X(L)-X(I))
            INSIDE = (CIJL.GE.0.0 .AND. CJKL.GE.0.0 .AND. CKIL.GE.0.0)
     &          .OR. (CIJL.LE.0.0 .AND. CJKL.LE.0.0 .AND. CKIL.LE.0.0) 
C However, if all three cross-products are zero, point lies on same straight
C line as the other three points, so need to determine if interpolation 
C ("inside") or extrapolation ("outside"). Do this by looking for at least one 
C negative dot product between LI.LJ, LJ.LK and LK.LI
            IF ( CIJL.EQ.0.0 .AND. CJKL.EQ.0.0 .AND. CKIL.EQ.0.0 ) THEN
              DLIJ = (X(I)-X(L))*(X(J)-X(L)) + (Y(I)-Y(L))*(Y(J)-Y(L))
              DLJK = (X(J)-X(L))*(X(K)-X(L)) + (Y(J)-Y(L))*(Y(K)-Y(L))
              DLKI = (X(K)-X(L))*(X(I)-X(L)) + (Y(K)-Y(L))*(Y(I)-Y(L))
              INSIDE = DLIJ.LT.0.0 .OR. DLJK.LT.0.0 .OR. DLKI.LT.0.0
            END IF
            IF ( INSIDE ) THEN
              IDXTRI(ITRI,3) = L        ! ITRI becomes IJL
              IDXTRI(NTRI+1,1) = I      ! NTRI+1 becomes IKL
              IDXTRI(NTRI+1,2) = K
              IDXTRI(NTRI+1,3) = L
              IDXTRI(NTRI+2,1) = J      ! NTRI+1 becomes JKL
              IDXTRI(NTRI+2,2) = K
              IDXTRI(NTRI+2,3) = L
              NTRI = NTRI + 2           ! 2 more triangles added (net effect)
              GOTO 100
            END IF
          END DO
          STOP 'F-TRIANG: logical error#1'  ! All points must lie in a triangle
 100      CONTINUE
        END IF
      END DO
C
C At this point the entire 2D field has been triangulated 
C
C Step#5: Optimising the triangulation (try to make triangles more equilateral)
C Look for pairs of triangles IJK, KLI which share a side IK
 200  CONTINUE
      FLIPMX = 0.0
      DO JTRI = 1, NTRI-1                       ! Loop over all triangles
        DO ISID = 1, 3                          ! Loop over each side of JTRI
          I = IDXTRI(JTRI,ISID)                 ! I = 1, 2, 3
          J = IDXTRI(JTRI,1+MOD(ISID,3))        ! J = 2, 3, 1
          K = IDXTRI(JTRI,1+MOD(ISID+1,3))      ! K = 3, 1, 2
          DO LTRI = JTRI+1, NTRI                ! Compare with other triangles
            IF ( ( I .EQ. IDXTRI(LTRI,1) .OR.
     &             I .EQ. IDXTRI(LTRI,2) .OR.
     &             I .EQ. IDXTRI(LTRI,3)      ) .AND.
     &           ( K .EQ. IDXTRI(LTRI,1) .OR.
     &             K .EQ. IDXTRI(LTRI,2) .OR.
     &             K .EQ. IDXTRI(LTRI,3)      ) ) THEN  ! Found shared boundary
              L = IDXTRI(LTRI,1)                        ! Establish 4th point, L
              IF ( L .EQ. I .OR. L .EQ. K ) THEN
                L = IDXTRI(LTRI,2)
                IF ( L .EQ. I .OR. L .EQ. K ) THEN
                  L = IDXTRI(LTRI,3)
                END IF
              END IF
C We now have a quadrilateral IJKL. To allow the diagonal IK to be flipped to 
C JL it is necessary that each corner is convex (inside angle up to 180deg)
C For this, the sine (ie cross-product) between vectors representing 
C successive sides IJ, JK, JL, LI must all be the same sign.
              CIJK = (X(J)-X(I))*(Y(K)-Y(J))-(X(K)-X(J))*(Y(J)-Y(I))
              CJKL = (X(K)-X(J))*(Y(L)-Y(K))-(X(L)-X(K))*(Y(K)-Y(J))
              CKLI = (X(L)-X(K))*(Y(I)-Y(L))-(X(I)-X(L))*(Y(L)-Y(K))
              CLIJ = (X(I)-X(L))*(Y(J)-Y(I))-(X(J)-X(I))*(Y(I)-Y(L))
              IF ( ( CIJK .GE. 0.0 .AND. CJKL .GE. 0.0 .AND. 
     &               CKLI .GE. 0.0 .AND. CLIJ .GE. 0.0       ) .OR.
     &             ( CIJK .LE. 0.0 .AND. CJKL .LE. 0.0 .AND. 
     &               CKLI .LE. 0.0 .AND. CLIJ .LE. 0.0       ) ) THEN
C
C Flip if the circle through one triangle includes the fourth point.
C This is true if the sum of angles IJK+KLI > 180. If this is true, then it
C will not be true for the other pair of corners (since sum all angles=360),
C so "benefit" is increase in sin(sum opposite corners).
C
                DIJK = (X(I)-X(J))*(X(K)-X(J))+(Y(I)-Y(J))*(Y(K)-Y(J))
                DKLI = (X(K)-X(L))*(X(I)-X(L))+(Y(K)-Y(L))*(Y(I)-Y(L))
                SIJ = SQRT ( (X(J)-X(I))**2 + (Y(J)-Y(I))**2 )
                SJK = SQRT ( (X(K)-X(J))**2 + (Y(K)-Y(J))**2 )
                SKL = SQRT ( (X(L)-X(K))**2 + (Y(L)-Y(K))**2 )             
                SLI = SQRT ( (X(I)-X(L))**2 + (Y(I)-Y(L))**2 )
                COSIJK = DIJK / ( SIJ * SJK )
                COSKLI = DKLI / ( SKL * SLI )
                SINE = COSIJK * SQRT ( MAX ( 0.0, 1.0 - COSKLI**2 ) ) +
     &                 COSKLI * SQRT ( MAX ( 0.0, 1.0 - COSIJK**2 ) )
                FLIP = - SINE
C
C If two "flat" triangles (so cos=-1), swap diagonal to shorten boundary
                IF (COSIJK+1.0 .LE. SMALL .AND. COSKLI+1.0 .LE. SMALL) 
     &            FLIP = 2.0
C
                IF ( FLIP .GT. FLIPMX ) THEN      ! Found best so far
                  FLIPMX = FLIP
                  II = I
                  JJ = J
                  KK = K
                  LL = L
                  KTRI = JTRI
                  ITRI = LTRI
                END IF     ! end test for best flip so far
              END IF       ! end test for flippable diagonals
            END IF         ! end test for shared boundary
          END DO           ! end loop over LTRI
        END DO             ! end loop over sides of JTRI
      END DO               ! end loop over JTRI
C
C Flip diagonal for pair of triangles which gives most improvement
C So triangles IJK, KLI become LIJ, JKL
      IF ( FLIPMX .GT. SMALL ) THEN
        IDXTRI(ITRI,1) = LL
        IDXTRI(ITRI,2) = II
        IDXTRI(ITRI,3) = JJ
        IDXTRI(KTRI,1) = JJ
        IDXTRI(KTRI,2) = KK
        IDXTRI(KTRI,3) = LL
        GOTO 200               ! Repeat until no more improvements found
      END IF
C
C Jump here if triangulation already completed for this grid
 300  CONTINUE
C 
C First, find closest grid point to interpolation point and set to return value
C at this point if no other location is found.
      I = 1
      RMIN = ( XO - X(1) )**2 + ( YO - Y(1) )**2
      DO J = 2, N
        R = ( XO - X(J) )**2 + ( YO - Y(J) )**2
        IF ( R .LT. RMIN ) THEN
          I = J
          RMIN = R
        END IF
      END DO
      IFIT(1) = I
      IFIT(2) = 1
      IFIT(3) = 1
      WFIT(1) = 1.0
      WFIT(2) = 0.0
      WFIT(3) = 0.0
      XFIT = X(I)
      YFIT = Y(I)
      IF ( RMIN .EQ. 0.0 ) RETURN
      SMIN = 0.0
C
C Look for existing triangle enclosing interpolated point.
C Enclosing defined by cross-products, as before, but normalise by square of
C length of side to avoid choosing "flat" triangles
      DO ITRI = 1, NTRI
        IF ( FINITE ) THEN
          I = IDXTRI(ITRI,1)
          J = IDXTRI(ITRI,2)
          K = IDXTRI(ITRI,3)
          CIJO = ( X(J)-X(I))*(YO-Y(J)) - (Y(J)-Y(I))*(XO-X(J) ) 
          CIJTOL = SMALL * SQRT (( (X(J)-X(I))**2 + (Y(J)-Y(I))**2 )* 
     &                           ( ( XO -X(J))**2 + ( YO -Y(J))**2 ))
          CJKO = ( X(K)-X(J))*(YO-Y(K)) - (Y(K)-Y(J))*(XO-X(K) ) 
          CJKTOL = SMALL * SQRT (( (X(K)-X(J))**2 + (Y(K)-Y(J))**2 )*
     &                           ( ( XO -X(K))**2 + ( YO -Y(K))**2 ))
          CKIO = ( X(I)-X(K))*(YO-Y(I)) - (Y(I)-Y(K))*(XO-X(I) ) 
          CKITOL = SMALL * SQRT (( (X(I)-X(K))**2 + (Y(I)-Y(K))**2 )*
     &                           ( ( XO -X(I))**2 + ( YO -Y(I))**2 ))
          IF ( ( CIJO .GT. CIJTOL .AND. 
     &           CJKO .GT. CJKTOL .AND. CKIO .GT. CKITOL ) .OR. 
     &         ( CIJO .LT. -CIJTOL .AND. 
     &           CJKO .LT. -CJKTOL .AND. CKIO .LT. -CKITOL )    ) THEN
C
C Fit plane to points I,J,K and interpolate plane to value at (XO,YO).
C This next bit of code effectively does both by inverting the 3x3 matrix to
C establish the weights (WFIT) for I,J,K in the final interpolation directly
            IFIT(1) = I
            IFIT(2) = J
            IFIT(3) = K
            FACTX = 1.0/( ( X(I)-X(J) ) * ( Y(J)-Y(K) ) - 
     &                    ( X(J)-X(K) ) * ( Y(I)-Y(J) )   )
            FACTY = 1.0/( ( Y(I)-Y(J) ) * ( X(J)-X(K) ) - 
     &                    ( Y(J)-Y(K) ) * ( X(I)-X(J) )   )
            FACTO = 1.0/( X(K)*Y(J) - X(J)*Y(K) + 
     &                    X(I)*Y(K) - X(K)*Y(I) + 
     &                    X(J)*Y(I) - X(I)*Y(J)   )
            WFIT(1) = XO * ( Y(J)-Y(K) ) * FACTX + 
     &                YO * ( X(J)-X(K) ) * FACTY + 
     &                ( Y(J)*X(K) - X(J)*Y(K) ) * FACTO
            WFIT(2) = XO * ( Y(K)-Y(I) ) * FACTX +
     &                YO * ( X(K)-X(I) ) * FACTY +
     &                ( Y(K)*X(I) - X(K)*Y(I) ) * FACTO
            WFIT(3) = XO * ( Y(I)-Y(J) ) * FACTX +
     &                YO * ( X(I)-X(J) ) * FACTY +
     &                ( Y(I)*X(J) - X(I)*Y(J) ) * FACTO
            XFIT = XO
            YFIT = YO
            RETURN             ! Exit with triangulation completed
          END IF
        END IF
C For each side of triangle in turn, find if the interpolation point lies 
C perpendicular ie if the triangle OIJ has acute angles at corners I and J
C so that both dot products are > 0. Note that the previous point-by-point
C test will have already found the case of dot product=0.
        DO ISID = 1, 3
          I = IDXTRI(ITRI,ISID)                ! I = 1, 2, 3
          J = IDXTRI(ITRI,1+MOD(ISID,3))       ! J = 2, 3, 1
          DIJO = (X(I)-X(J)) * (XO-X(J)) + (Y(I)-Y(J)) * (YO-Y(J))
          DJIO = (X(J)-X(I)) * (XO-X(I)) + (Y(J)-Y(I)) * (YO-Y(I))
          IF (DIJO .GT. 0.0 .AND. DJIO .GT. 0.0 )  THEN
C Project interpolation point O to point P along line IJ. 
C IF OP < RMIN, save this point (P) as potentially closest point.
C If OP = RMIN, only save this point if length of IJ (=S) is shorter than
C previous (representing a linear interpolation between two closer points).
C The distance JP = JO.cos(IJO) = DIJO/IJ since DIJO = IJ.JO.cos(IJO)
C So from right-angle triang JPO, distance OP is given by OP^2 = JO^2 - JP^2
C But DIJO = IJ.JO.cos(IJO) so cos(IJO) = DIJO/(IJ.JO)
C Also: DIJO + DJIO = IJ.(JO.cos(IJO)+IO.cos(JIO)) = IJ^2. 
            SJP2 = DIJO**2 / ( (X(J)-X(I))**2 + (Y(J)-Y(I))**2 )  
            R = (X(J)-XO)**2 + (Y(J)-YO)**2 - SJP2
            S = DIJO + DJIO
            RTOL = SJP2*SMALL + S*SMALL
            IF ( R .LT. RMIN-RTOL .OR. 
     &           ( R .LT. RMIN+RTOL .AND. S.LT.SMIN) ) THEN
C
C Weight of point I is then JP/IJ. But JP = DIJO/IJ, so Wgt=JP^2/DIJO
C But DIJO = IJ.JO.cos(IJO), so PJ/IJ=DIJO/(IJ)**2 
C NB: DIJO=0 shouldn't arise since then R=RMIN from previous loop looking for
C closest point (J) and SMIN was set to 0.0
              IFIT(1) = I
              IFIT(2) = J
              IFIT(3) = 1
              WFIT(1) = SJP2 / DIJO 
              WFIT(2) = 1.0 - WFIT(1)
              WFIT(3) = 0.0
              XFIT = X(I)*WFIT(1) + X(J)*WFIT(2)
              YFIT = Y(I)*WFIT(1) + Y(J)*WFIT(2)
              RMIN = R
              SMIN = DIJO + DIJO
            END IF
          END IF
        END DO
      END DO
C
      END
