      SUBROUTINE RFMSCA ( FAIL, ERRMSG )
C
C VERSION
C     15-OCT-04  AD  Initialise FAIL=FALSE
C     16-JAN-04  AD  Set IDIR=0 in IDXPTH argument unless GRAFLG is TRUE
C     26-JAN-01  AD  Check NPTH > 0 before writing diagnostic message
C     30-DEC-99  AD  Add IDIR argument to IDXPTH
C                    Add check for temperature compatibility
C     02-FEB-99  AD  Set ISCPTH to path index if not scaled (was 0)
C     03-MAR-97  AD  Version 3.
C     20-DEC-96  AD  Correction: Error when for NCLC .GT. MAXCLC (was .GE.)
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Determine which path absorption calculations can be scaled.
C     Called by RFM in limb-viewing mode unless CLC option selected.
C
      IMPLICIT NONE
C
C ARGUMENTS
      LOGICAL      FAIL   ! O  T=A fatal error has occurred
      CHARACTER*80 ERRMSG ! O  Error message written if FAIL is TRUE
C
      EXTERNAL
     &  IDXPTH ! Indexing function for Pth# from Tan.pt#, Atm.Layer# and Gas#
     &, PTHSRT ! Sort Calculated paths to start of /PTHCOM/ list
     &, RFMLOG ! Write text message to RFM log file.
C
      INTEGER IDXPTH
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C LOCAL CONSTANTS
      INTEGER NEXTRA               ! No.layers above tangent layer for explicit
        PARAMETER ( NEXTRA = 1 )   ! line-by-line calculation
      REAL PCGMAX                  ! Allowed difference in ln(p) for scaling
        PARAMETER ( PCGMAX = 0.01 ) ! 0.01 = 1% 
      REAL TCGMAX                  ! Allowed difference in T for scaling
        PARAMETER ( TCGMAX = 1.0 ) ! 1.0 = 1K
C
C COMMON VARIABLES
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'pthcom.inc' ! RFM path data
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL VARIABLES
      INTEGER      IATM    ! Atmospheric profile layer
      INTEGER      IDIR    ! Direction pointer
      INTEGER      IGAS    ! Gas counter
      INTEGER      IPTH    ! Path counter
      INTEGER      ITAN    ! Tangent path counter
      INTEGER      JPTH    ! Path counter
      INTEGER      NDIR    ! 0=direction irrelevant, 1=both up and down reqd.
      REAL         PCGDIF  ! Difference in Ln(C-G Press) between two paths
      REAL         TCGDIF  ! Difference in (C-G Temp) between two paths
      CHARACTER*80 MESSGE  ! Info message sent to Log file
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      FAIL = .FALSE.
C Initially mark all paths for no line-by-line calculation
C
      DO IPTH = 1, NPTH
        CLCPTH(IPTH) = .FALSE.
      END DO
      NCLC = 0
C
      IF ( GRAFLG ) THEN 
        NDIR = 1
      ELSE
        NDIR = -1
      END IF
C
C Set all paths in tangent layer + NEXTRA layers for line-by-line calculation
C
      DO ITAN = 1, MTAN
        IF ( CLCTAN(ITAN) ) THEN                 ! Tangent path calc. required
          DO IATM = IATTAN(ITAN), MIN ( NATM-1, IATTAN(ITAN)+NEXTRA)
            DO IGAS = 1, NGAS
              DO IDIR = -1, NDIR, 2
                IF ( GRAFLG ) THEN
                  IPTH = IDXPTH ( ITAN, IATM, IGAS, IDIR )
                ELSE
                  IPTH = IDXPTH ( ITAN, IATM, IGAS, 0 )
                END IF
                ISCPTH(IPTH) = IPTH
                CLCPTH(IPTH) = .TRUE.
                NCLC = NCLC + 1
              END DO
            END DO
          END DO
        END IF
      END DO
C 
C For all higher paths, see if there is a similar path marked for calculation
C which can be scaled, otherwise flag for calculation (needs to start from
C tangent level to catch any paths not marked for calculation by CLCTAN=T.)
C
      DO ITAN = 1, MTAN
        DO IATM = IATTAN(ITAN), NATM-1
          DO IGAS = 1, NGAS
            DO IDIR = -1, NDIR, 2
              IF ( GRAFLG ) THEN
                IPTH = IDXPTH ( ITAN, IATM, IGAS, IDIR )
              ELSE
                IPTH = IDXPTH ( ITAN, IATM, IGAS, 0 )
              END IF
              IF ( .NOT. CLCPTH(IPTH) ) THEN
                DO JPTH = 1, NPTH
                  IF ( CLCPTH(JPTH) .AND. IATM .EQ. IATPTH(JPTH) .AND.
     &                 IGAS .EQ. IGSPTH(JPTH) ) THEN
                    PCGDIF = ABS ( LOG ( PREPTH(IPTH)/PREPTH(JPTH) ) )
                    TCGDIF = ABS ( TEMPTH(IPTH) - TEMPTH(JPTH) )
                    IF ( PCGDIF .LT. PCGMAX .AND. 
     &                   TCGDIF .LT. TCGMAX       ) THEN
                      ISCPTH(IPTH) = JPTH
                      GOTO 100
                    END IF
                  END IF
                END DO
C
C If the end of the loop is reached no similar path has been found 
C so flag for calculation
C
                CLCPTH(IPTH) = .TRUE.
                ISCPTH(IPTH) = IPTH
                NCLC = NCLC + 1
  100           CONTINUE
              END IF
            END DO
          END DO
        END DO
      END DO
C
      IF ( NPTH .GT. 0 ) THEN
        WRITE ( MESSGE, '(A,I8,A,I8,A,I3,A)' )
     &    'I-RFMSCA: LBL Absorption will be calc. for', NCLC, 
     &    ' paths out of', NPTH, '(=', NINT(100.0*NCLC/FLOAT(NPTH)),'%)'
        CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
        IF ( FAIL ) RETURN
      END IF
C
      IF ( NCLC .GT. MAXCLC ) THEN
        WRITE ( ERRMSG, '(A,I11,A,I7)' )
     &    'F-RFMSCA: No.paths for line-by-line calc=', NCLC,
     &    ' >MAXCLC in RFMSIZ=', MAXCLC
        FAIL = .TRUE.
      ELSE
        CALL PTHSRT
      END IF
C
      END
