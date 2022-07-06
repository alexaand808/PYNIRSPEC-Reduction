      SUBROUTINE LOSNAM ( ISPC, ITAN, INPNAM, OUTNAM, HEADRJ )
C
C VERSION
C     07-AUG-13  AD  Remove redundant local varbls IGAS,IATM,IHGT,GASSTR,HGTSTR
C     08-JUN-01  AD  Original. Based on JACNAM.
C
C DESCRIPTION
C     Construct filename for RFM LOS-Jacobian files.
C     Called by OPNOUT for each tangent height (if LOS flag enabled).
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER       ISPC    !  I  Spectral range number
      INTEGER       ITAN    !  I  Tangent height number
      CHARACTER*(*) INPNAM  !  I  Filename template
      CHARACTER*(*) OUTNAM  !  O  Jacobian filename
      CHARACTER*(*) HEADRJ  !  O  Header record giving Jacobian info
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES 
      INCLUDE 'outcom.inc' ! RFM output file data
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL CONSTANTS
      CHARACTER*1   CMINUS ! Character used for negative altitudes
        PARAMETER ( CMINUS = '-' )
C
C LOCAL VARIABLES
      INTEGER      IEND   ! Length of filename so far
      INTEGER      IWLD   ! Pointer to wildcard in input filename
      INTEGER      L      ! Total length of INSERT
      INTEGER      LENGTH ! Length of substrings
      CHARACTER*80 INSERT ! Inserted component of filename
      CHARACTER*6  TANSTR ! Tangent height text
C
C EXECUTABLE CODE ------------------------------------------------------------
C
C Start by assuming output name = input name, removing all trailing spaces
C
      IEND = INDEX ( INPNAM, ' ' ) - 1            ! point to last character
      IF ( IEND .EQ. -1 ) IEND = LEN ( INPNAM )   ! no trailing spaces
C
C Remove wildcard from the input name (should always exist for Jacobian files)
C
      IWLD = INDEX ( INPNAM, '*' )
      IF ( IWLD .EQ. 0 ) STOP 'F-LOSNAM: Logical Error#1'
C
C Construct string to be inserted at wildcard position
      INSERT = ' '
      L = 0
C
C Insert any spectral range label
C
      IF ( ISPC .NE. 0 ) THEN
        LENGTH = INDEX ( LABSPC(ISPC)//' ', ' ' ) - 1
        IF ( LENGTH .GT. 0 )
     &      INSERT(L+1:L+LENGTH) = LABSPC(ISPC)(1:LENGTH)
          L = L + LENGTH
      END IF
C
C Insert tangent height string
      IF ( ITAN .EQ. 0 ) STOP 'F-LOSNAM: Logical Error#2'
      IF ( STRTAN(ITAN)(1:1) .EQ. '-' ) THEN
        TANSTR = CMINUS//STRTAN(ITAN)(2:6)
        LENGTH = 6
      ELSE
        TANSTR = STRTAN(ITAN)(1:5)
        LENGTH = 5
      END IF
      INSERT(L+1:L+LENGTH) = TANSTR(1:LENGTH)
      L = L + LENGTH
C
C Insert '_' character
      INSERT(L+1:L+4) = '_los'
      L = L + 4
C
C Construct Jacobian filename
      IF ( IWLD .EQ. 1 ) THEN
        OUTNAM = INSERT(1:L)//INPNAM(2:IEND)
      ELSE IF ( IWLD .EQ. IEND ) THEN
        OUTNAM = INPNAM(1:IEND-1)//INSERT(1:L)
      ELSE
        OUTNAM = INPNAM(1:IWLD-1)//INSERT(1:L)//INPNAM(IWLD+1:IEND)
      END IF
      IEND = IEND + L - 1   ! Substract 1 since also removing '*' char.
C
      IF ( USRELE ) THEN
        HEADRJ =
     &  '! Pointing Jacobian spectrum for 0.001 deg perturbation'
      ELSE 
        HEADRJ =
     &  '! Pointing Jacobian spectrum for 0.001 km perturbation'
      END IF
C
C If a Run ID is specified, append that to filename.
C
      IF ( RUNID .NE. ' ' ) THEN
        LENGTH = INDEX ( RUNID//' ', ' ' ) - 1  
        OUTNAM = OUTNAM(1:IEND)//RUNID(1:LENGTH)
        IEND = IEND + LENGTH
      END IF
C
      END
