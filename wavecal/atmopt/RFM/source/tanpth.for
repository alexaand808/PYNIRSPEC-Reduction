      SUBROUTINE TANPTH ( FAIL, ERRMSG )
C
C VERSION
C     16-APR-10  AD  Increase length of character field in MAXPTH error message
C     17-FEB-04  AD  Initialise NLNPTH
C     16-JAN-04  AD  Set IDRPTH=0 unless GRAFLG = T
C     19-AUG-02  AD  Initialise FAIL to FALSE
C     30-DEC-99  AD  Set IDIR to distinguish up and down paths (GRAFLG)
C     02-FEB-99  AD  Set ISCPTH to current path index (was 0)
C                    Use CLCFLG instead of .NOT. SCAFLG
C     30-DEC-98  AD  Comment change only.
C     14-JUL-98  AD  Comment change only.
C     03-MAR-97  AD  Version 3.
C     20-DEC-96  AD  Change I4 to I5 in error message to avoid overflows.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Construct list of RFM paths from tang.hts and atmos.profiles.
C     Called once by INPCHK.
C     A path is defined for each combination of absorber and atmos layer.
C
      IMPLICIT NONE
C
C ARGUMENTS
      LOGICAL      FAIL   !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG !  O  Contains error message if FAIL is TRUE.
C
      EXTERNAL
     &  PTHSRT ! Sort Calculated paths to start of /PTHCOM/ list.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'atmcom.inc' ! Atmospheric profile data  
      INCLUDE 'tancom.inc' ! Tangent heights
      INCLUDE 'pthcom.inc' ! RFM path data
C
C LOCAL VARIABLES
      INTEGER  IGAS   !  Gas counter
      INTEGER  ITAN   !  Tangent path counter
      INTEGER  IATM   !  Atmospheric layer counter
      INTEGER  IDIR   !  Direction pointer for path 
      INTEGER  JATM   !  Index of layer containing tangent point
      INTEGER  NDIR   !  -1=no horizontal gradient, +1=use gradient
      INTEGER  NEXTRA !  No.of extra PTHs required for this tangent height
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      FAIL = .FALSE.
      NPTH = 0
      NCLC = 0
      IF ( GRAFLG ) THEN
        NDIR = 1
      ELSE
        NDIR = -1
      END IF
      DO ITAN = 1, MTAN           
        JATM = NATM - 1
        DO WHILE ( HGTTAN(ITAN) .LT. HGTATM(JATM) .AND. JATM .GT. 1 )
          JATM = JATM - 1
        END DO
        NEXTRA = NGAS * ( NATM - JATM ) 
        IF ( GRAFLG ) NEXTRA = NEXTRA * 2
        IF ( NPTH + NEXTRA .GT. MAXPTH ) THEN
          WRITE ( ERRMSG, '(A,F8.3,A,I6,A,I6)' )
     &      'F-TANPTH: Up to Tan.Ht=', HGTTAN(ITAN),
     &      'km No.Paths reqd NPTH=', NPTH + NEXTRA,
     &      ' > MAXPTH=', MAXPTH
          FAIL = .TRUE.
          RETURN
        END IF
C
        DO IGAS = 1, NGAS
          DO IATM = JATM, NATM-1
            DO IDIR = -1, NDIR, 2
              NPTH = NPTH + 1
              IGSPTH(NPTH) = IGAS
              ITNPTH(NPTH) = ITAN
              IATPTH(NPTH) = IATM
              IF ( GRAFLG ) THEN
                IDRPTH(NPTH) = IDIR
              ELSE
                IDRPTH(NPTH) = 0
              END IF
              CLCPTH(NPTH) = CLCTAN(ITAN)
              ISCPTH(NPTH) = NPTH
              NLNPTH(NPTH) = 0
              IF ( CLCPTH(NPTH) ) NCLC = NCLC + 1
            END DO
          END DO
        END DO
        IATTAN(ITAN) = JATM
      END DO
C
C If not using scaled paths, then current value of NCLC is the number of 
C calculated paths that will be required. If using scaled paths, this will be
C checked later in RFMSCA.
C
      IF ( CLCFLG ) THEN
        IF ( NCLC .GT. MAXCLC ) THEN
          WRITE ( ERRMSG, '(A,I11,A,I7)' )
     &      'F-TANPTH: No.paths for line-by-line calc=', NCLC,
     &      ' >MAXCLC in RFMSIZ=', MAXCLC
          FAIL = .TRUE.
        ELSE
          CALL PTHSRT
        END IF
      END IF
C
      END
