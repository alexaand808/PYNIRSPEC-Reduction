      SUBROUTINE MAKNAM ( ISPC, IGAS, ITAN, IGRA, INPNAM, OUTNAM )
C
C VERSION
C     09-MAY-04  AD  Allow for isotopic species
C     06-JAN-04  AD  Add IGRA argument
C     05-DEC-00  AD  Avoid concatenation of variable character strings
C     08-DEC-99  AD  Rewritten to avoid non-standard FORTRAN usage
C     18-JUL-98  AD  Use STRTAN for tan.hgt identifier. Add CMINUS
C     31-OCT-97  AD  Change sequence to write GAS after SPC label.      
C     22-OCT-97  AD  Only write gas name if wildcard exists.
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     18-SEP-96  AD  Remove redundant LOCASE
C                     Always insert tangent height if there is a wildcard
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Construct filename for RFM output files.
C     General purpose.
C     CMINUS can be changed to some other character if operating system does 
C     not allow '-' sign in filenames.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER       ISPC   !  I  Spectral range number
      INTEGER       ITAN   !  I  Tangent height number
      INTEGER       IGAS   !  I  Absorber number
      INTEGER       IGRA   !  I  Horizontal (psi) angle number
      CHARACTER*(*) INPNAM !  I  Input (basic) filename
      CHARACTER*(*) OUTNAM !  O  Output filename
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES 
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'gracom.inc' ! Atmospheric 2-D field data
      INCLUDE 'outcom.inc' ! RFM output file data
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL CONSTANTS
      CHARACTER*1  CMINUS  ! Character used for negative tan.hgts, psi angles
        PARAMETER ( CMINUS = '-' )
      CHARACTER*1   CISO   ! Character used before Isotope# in filename
        PARAMETER ( CISO = 'i' )
C
C LOCAL VARIABLES
      INTEGER      IEND    ! Length of filename so far
      INTEGER      IWLD    ! Index of wildcard in input filename
      INTEGER      L       ! Total length of INSERT
      INTEGER      LENGTH  ! Length of sub-string to be inserted
      CHARACTER*80 INSERT  ! Inserted component of filename
      CHARACTER*2  ISOSTR  ! Isotope component of filename
      CHARACTER*6  TANSTR  ! Tangent height text 
      CHARACTER*6  GRASTR  ! Horizontal angle text 
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Start by assuming output name = input name, removing all trailing spaces
      IEND = INDEX ( INPNAM, ' ' ) - 1          ! point to last character
      IF ( IEND .EQ. -1 ) IEND = LEN ( INPNAM ) ! No trailing spaces in INPNAM 
C
C Find if there is a wildcard in the input name. If there is, note its position
C and move it to the front of the filename temporarily (or leave it there if
C that's where it is already).
C
      IWLD = INDEX ( INPNAM, '*' )
C
      IF ( IWLD .EQ. 0 ) THEN      ! No wildcard, so keep same filename
        OUTNAM = INPNAM(1:IEND)
      ELSE                         ! Wildcard, so insert strings (if any)
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
C Insert any absorber name 
C
        IF ( IGAS .NE. 0 ) THEN
          LENGTH = INDEX ( CODGAS(IGAS)//' ', ' ' ) - 1
          IF ( LENGTH .GT. 0 )
     &      INSERT(L+1:L+LENGTH) = CODGAS(IGAS)(1:LENGTH)
          L = L + LENGTH
C
C Insert isotopic information (if relevant)
          IF ( ISOMOL(IDXGAS(IGAS)) ) THEN               ! Split by isotope
            WRITE ( ISOSTR, '(A,I1)' ) CISO, IDIGAS(IGAS)
            LENGTH = LEN ( ISOSTR )
            INSERT(L+1:L+LENGTH) = ISOSTR
            L = L + LENGTH
          END IF
        END IF
C
C Insert any tangent height
C
        IF ( ITAN .NE. 0 ) THEN
          IF ( STRTAN(ITAN)(1:1) .EQ. '-' ) THEN
            TANSTR = CMINUS//STRTAN(ITAN)(2:6)
            LENGTH = 6
          ELSE
            TANSTR = STRTAN(ITAN)(1:5)
            LENGTH = 5
          END IF
          INSERT(L+1:L+LENGTH) = TANSTR(1:LENGTH)
          L = L + LENGTH
        END IF
C
C Insert any horizontal (psi) angle 
C
        IF ( IGRA .NE. 0 ) THEN
          IF ( PSIGRA(IGRA) .GE. 0.0 ) THEN     ! Psi angles should be <99deg
            WRITE ( GRASTR, '(I5.5)' ) NINT ( PSIGRA(IGRA) * 1000.0 )
            LENGTH = 5
          ELSE
            WRITE ( GRASTR, '(A,I5.5)' ) 
     &        CMINUS, -NINT ( PSIGRA(IGRA) * 1000.0 )
            LENGTH = 6
          END IF
          INSERT(L+1:L+LENGTH) = GRASTR(1:LENGTH)
          L = L + LENGTH
        END IF
C
C Construct new filename name including any insertion
C
        IF ( L .GT. 0 ) THEN
          IF ( IWLD .EQ. 1 ) THEN
            OUTNAM = INSERT(1:L)//INPNAM(2:IEND)
          ELSE IF ( IWLD .EQ. IEND ) THEN
            OUTNAM = INPNAM(1:IEND-1)//INSERT(1:L)
          ELSE
            OUTNAM = INPNAM(1:IWLD-1)//INSERT(1:L)//INPNAM(IWLD+1:IEND)
          END IF
        ELSE
          IF ( IWLD .EQ. 1 ) THEN
            OUTNAM = INPNAM(2:IEND)
          ELSE IF ( IWLD .EQ. IEND ) THEN
            OUTNAM = INPNAM(1:IEND-1)
          ELSE
            OUTNAM = INPNAM(1:IWLD-1)//INPNAM(IWLD+1:IEND)
          END IF
        END IF
C
C Update IEND to point to last character in OUTNAM 
C
        IEND = IEND + L - 1  ! NB subtract 1 since '*' character is removed
C
      END IF
C
C If a Run ID is specified, append that to filename.
C
      IF ( RUNID .NE. ' ' ) THEN
        LENGTH = INDEX ( RUNID//' ', ' ' ) - 1  
        IF ( LENGTH .GT. 0 ) THEN
          OUTNAM(IEND+1:IEND+LENGTH) = RUNID(1:LENGTH)
          IEND = IEND + LENGTH
        END IF
      END IF
C
      END
