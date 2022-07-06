      SUBROUTINE JACNAM ( ISPC, ITAN, IJAC, INPNAM, OUTNAM, HEADRJ )
C
C VERSION
C     29-AUG-13  AD  Remove redundant local variable IISO
C     17-MAR-09  AD  Bug#73: Avoid accessing HGTATM(IATM) if IATM=0
C     03-MAY-05  AD  Add jdxcon.inc. Allow for VT Jacobians
C     18-MAR-04  AD  Change to use 'iN' rather than '(N)' for Isotope#N
C     01-JAN-04  AD  Define 'down' an 'up' for use with LEV flag
C                    Allow for isotopic Jacobians
C     03-SEP-03  AD  Only write ptb altitude if IATJAC > 0
C     05-DEC-00  AD  Avoid concatenating variable length strings
C     23-APR-00  AD  Allow for MTXFLG
C     08-DEC-99  AD  Rewritten to avoid non-standard FORTRAN usage
C     31-JUL-99  AD  Remove unused variables ISTA, LENEND, STREND
C     06-APR-99  AD  Correction to calculation of LENGAS.
C     02-JAN-99  AD  Original. Based on MAKNAM
C
C DESCRIPTION
C     Construct filename for RFM Jacobian files.
C     Called by OPNOUT for each Jacobian file.
C     Also called if MTX option selected.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER       ISPC    !  I  Spectral range number
      INTEGER       ITAN    !  I  Tangent height number
      INTEGER       IJAC    !  I  Index of Jacobian element
      CHARACTER*(*) INPNAM  !  I  Filename template
      CHARACTER*(*) OUTNAM  !  O  Jacobian filename
      CHARACTER*(*) HEADRJ  !  O  Header record giving Jacobian info
C
C GLOBAL CONSTANTS
      INCLUDE 'jdxcon.inc' ! Indices for non-VMR Jacobians
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
C
C COMMON VARIABLES 
      INCLUDE 'atmcom.inc' ! Atmospheric profile data
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'jaccom.inc' ! Jacobian data
      INCLUDE 'ntecom.inc' ! Non-LTE data
      INCLUDE 'outcom.inc' ! RFM output file data
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
      INCLUDE 'tancom.inc' ! Tangent heights
C
C LOCAL CONSTANTS
      CHARACTER*1   CMINUS ! Character used for negative altitudes
        PARAMETER ( CMINUS = '-' )
      CHARACTER*1   CISO   ! Character used before Isotope# in filename
        PARAMETER ( CISO = 'i' )
      CHARACTER*1   CVIB   ! Character used before Vib.Glob.Quant# in filename
        PARAMETER ( CVIB = 'v' )
C
C LOCAL VARIABLES
      INTEGER      IATM   ! Atmos.profile level of Jacobian element
      INTEGER      IDG    ! HITRAN index of gas
      INTEGER      IEND   ! Length of filename so far
      INTEGER      IGAS   ! Index of Jacobian species
      INTEGER      IHGT   ! Height [m or km] of Jacobian element
      INTEGER      INTE   ! Index of vibrational level for VT Jacobians
      INTEGER      IWLD   ! Pointer to wildcard in input filename
      INTEGER      L      ! Total length of INSERT
      INTEGER      LENGTH ! Length of substrings
      CHARACTER*7  GASSTR ! Code for Jacobian species, eg 'o3', 'tem'
      CHARACTER*6  HGTSTR ! Perturbation height text 
      CHARACTER*80 INSERT ! Inserted component of filename
      CHARACTER*2  ISOSTR ! Isotope component of filename
      CHARACTER*3  VIBSTR ! Global Quantum# component of filename
      CHARACTER*6  TANSTR ! Tangent height text
C
C EXECUTABLE CODE ------------------------------------------------------------
C
C Start by assuming output name = input name, removing all trailing spaces
C
c      IEND = INDEX ( INPNAM//' ', ' ' ) - 1
      IEND = INDEX ( INPNAM, ' ' ) - 1            ! point to last character
      IF ( IEND .EQ. -1 ) IEND = LEN ( INPNAM )   ! no trailing spaces
C
C Remove wildcard from the input name (should always exist for Jacobian files)
C
      IWLD = INDEX ( INPNAM, '*' )
      IF ( IWLD .EQ. 0 ) STOP 'F-JACNAM: Logical Error#1'
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
      IF ( ITAN .EQ. 0 ) STOP 'F-JACNAM: Logical Error#2'
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
      INSERT(L+1:L+1) = '_'
      L = L + 1
C
C Insert perturbed parameter string (not with MTX flag)
      IGAS = IGSJAC(IJAC)
      IF ( IGAS .GT. 0 ) THEN       ! JAC,LEV flags define IGSJAC, MTX leaves=0
        IF ( IGAS .EQ. MAXGAS+JDXTEM ) THEN
          GASSTR = 'tem'
        ELSE IF ( IGAS .EQ. MAXGAS+JDXPRE ) THEN
          GASSTR = 'pre'
        ELSE IF ( IGAS .EQ. MAXGAS+JDXSFT ) THEN
          GASSTR = 'sfctem'
        ELSE IF ( IGAS .EQ. MAXGAS+JDXSFE ) THEN
          GASSTR = 'sfcems'
C LEV flag: Values for 'up' and 'down' are assigned in subroutine LEVTAN
        ELSE IF ( IGAS .EQ. MAXGAS+JDXLVU ) THEN   ! Used with LEV flag
          GASSTR = 'up'
        ELSE IF ( IGAS .EQ. MAXGAS+JDXLVD ) THEN   ! Used with LEV flag
          GASSTR = 'down'
C VT Jacobians require gas, isotope and global quantum number
        ELSE IF ( IGAS .GT. MAXGAS+JDXNTE ) THEN   ! Just gas code for now
          INTE = IGAS - ( MAXGAS + JDXNTE )
          IDG  = IDGNTE(INTE)
          GASSTR = CODGAS(IGSMOL(IDG))
        ELSE
          GASSTR = CODGAS(IGAS)
        END IF
        LENGTH = INDEX ( GASSTR//' ', ' ' ) - 1
        INSERT(L+1:L+LENGTH) = GASSTR(1:LENGTH)
        L = L + LENGTH
      END IF
C
C Insert isotopic information (if relevant)
      IF ( IGAS .GT. 0 .AND. IGAS .LE. MAXGAS ) THEN  ! Only gas jacobians
        IF ( IDIGAS(IGAS) .NE. 0 ) THEN                         ! Isotope
          WRITE ( ISOSTR, '(A,I1)' ) CISO, IDIGAS(IGAS)
          LENGTH = LEN ( ISOSTR )
          INSERT(L+1:L+LENGTH) = ISOSTR
          L = L + LENGTH
        END IF
      END IF          
C
C Insert isotopic and global quantum number for VT Jacobians
      IF ( IGAS .GT. MAXGAS+JDXNTE ) THEN
        WRITE ( ISOSTR, '(A,I1)' ) CISO, ISONTE(INTE)
        LENGTH = LEN ( ISOSTR )
        INSERT(L+1:L+LENGTH) = ISOSTR
        L = L + LENGTH
        WRITE ( VIBSTR, '(A,I2.2)' ) CVIB, IGQNTE(INTE)
        LENGTH = LEN ( VIBSTR )
        INSERT(L+1:L+LENGTH) = VIBSTR
        L = L + LENGTH
      END IF        
C
C Insert perturbation altitude string
      IATM = IATJAC(IJAC) 
      IF ( IATM .GT. 0 ) THEN
        IHGT = NINT ( HGTATM(IATM)*1000.0 )
        IF ( ABS ( IHGT ) .GT. 99999 ) IHGT = NINT ( HGTATM(IATM) )
        IF ( IHGT .GE. 0 ) THEN
          WRITE ( HGTSTR, '(I5.5)' ) IHGT
          LENGTH = 5
        ELSE
          WRITE ( HGTSTR, '(A,I5.5)' ) CMINUS, ABS(IHGT)
          LENGTH = 6
        END IF
        INSERT(L+1:L+LENGTH) = HGTSTR(1:LENGTH)
        L = L + LENGTH
      ENDIF
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
      IF ( JACFLG ) THEN
        IF ( IATM .GT. 0 ) THEN
          WRITE ( HEADRJ, '(A,A7,A,F10.3,A)' ) 
     &      '! Jacobian spectrum for ', GASSTR, 
     &      ' perturbed at altitude = ', HGTATM(IATM), ' [km]'
        ELSE                   ! Surface or column perturbation
          WRITE ( HEADRJ, '(A,A7)' ) 
     &      '! Jacobian spectrum for ', GASSTR 
        END IF
      ELSE IF ( MTXFLG ) THEN
        WRITE ( HEADRJ, '(A,F10.3,A,F10.3,A)' )
     &    '! Flux matrix spectrum between levels ', HGTTAN(ITAN), 
     &    ' - ', HGTATM(IATM), ' km'
      ELSE                                           ! LEV flag
        WRITE ( HEADRJ, '(A,F10.3,A,A,A)' )
     &    '! Intermediate spectrum at level ', HGTATM(IATM), 
     &    ' km in ', GASSTR, ' path'
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
