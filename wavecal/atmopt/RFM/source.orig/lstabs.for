      SUBROUTINE LSTABS ( WNOL, WNOH, IMIN, NREQ, MOLLST, OPTLST )
C
C VERSION
C     07-AUG-13  AD  Add IMIN argument
C                    Updated for new HITRAN/RFM indices, 
C                    Increase MAXABS to 153, MAXOPT to 7858
C     06-MAR-02  AD  Increase MAXABS from 63 to 64, MAXOPT 6097 to 6106 (sf6)
C     16-JUN-01  AD  Change tabulation for revised optdat.inc
C     20-DEC-98  AD  Original.
C
C DESCRIPTION
C     Construct ordered list of absorbers for given spec.range.
C     On entry, NGAS specifies the number of absorbers required (value .GE. 45
C     will be guaranteed to return all absorbing species) and the arrays 
C     MOLLST and OPTLST should be externally dimensioned at least this size.
C     On exit, NGAS contains the number of absorbing species found, MOLLST
C     contains their HITRAN molecular ID's and OPTLST an absorption parameter
C     (largest value of OPTLST is first element, corresponding to maximum
C     absorption).      
C
      IMPLICIT NONE
C
C ARGUMENTS
      REAL     WNOL         !  I  Lower Wavenumber [cm-1] of required range
      REAL     WNOH         !  I  Upper Wavenumber [cm-1] of required range
      INTEGER  IMIN         !  I  Threshold for min. optical depth parameter
      INTEGER  NREQ         ! I/O I=Max No. of molecules, O=number found
      INTEGER  MOLLST(NREQ) !  O  Ordered list of HITRAN ID's for molecules
      REAL     OPTLST(NREQ) !  O  Value of optical depth parameter
C
C LOCAL CONSTANTS
      INTEGER       MAXOPT         ! No.points in OPTDAT.INC arrays
        PARAMETER ( MAXOPT = 7858 )
      INTEGER       MAXABS         ! Highest index of HITRAN-defined molecules.
        PARAMETER ( MAXABS = 153 )   
C
C LOCAL VARIABLES
      INTEGER   I      ! Counter for DATA statements in include file
      INTEGER   IDX    ! Offset index for molecule in data arrays OLOG10,WAVENO
      INTEGER   IGAS   ! Counter for absorbing species within range
      INTEGER   ILOG   ! Current value of OLOG10 array
      INTEGER   IMAX   ! Maximum value of ILOG for current molecule
      INTEGER   IMOL   ! Counter for all HITRAN species (1:MAXABS)
      INTEGER   ISUM   ! Summation of opt.dep.values for each molecule
      INTEGER   IWNO   ! Lower wavenumber of each molecules opt.dep. data range
      INTEGER   IWNOH  ! Integer wavenumber corresponding to WNOH
      INTEGER   IWNOL  ! Integer wavenumber corresponding to WNOL 
      INTEGER   JGAS   ! Position to insert current molecule into list
      INTEGER   JWNO   ! Upper wavenumber of each molecules opt.dep data range
      INTEGER   NGAS   ! No.absorbers found (output value of NREQ)
      REAL      FACTOR ! Factor for 'averaging' optical depth
      REAL      OPT    ! Combined maximum & average optical.depths
      INTEGER   OLOG10(MAXOPT) ! Optical depth spectra
      INTEGER   WAVENO(MAXOPT) ! Wavenumber indices for OLOG10
      INTEGER   IDXMOL(MAXABS) ! Offset indices for each molecule
C
C DATA STATEMENTS
      INCLUDE 'optdat.inc'   ! Spectral Optical depth data OLOG10,WAVENO,IDXMOL
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      NGAS  = 0                          ! Initialise no.of absorbers found
      IWNOL = INT ( WNOL )
      IWNOH = INT ( WNOH )
      FACTOR = 0.01 / FLOAT ( IWNOH - IWNOL + 1 ) 
      DO IMOL = 1, MAXABS
        IF ( IDXMOL(IMOL) .NE. 0 ) THEN  ! Exclude undefined molecular indices
          IDX = IDXMOL(IMOL)
          JWNO = WAVENO(IDX)             ! Upper limit of first range
          DO WHILE ( JWNO .LT. IWNOL )   ! Step over data ranges until overlap
            IDX = IDX + 1
            JWNO = JWNO + WAVENO(IDX)    ! Upper limit of next data range
          END DO    
          IWNO = JWNO - WAVENO(IDX) + 1  ! Lower wno of current data range   
          ILOG = OLOG10(IDX)             ! Value of Opt.Dep at WNOL
          IMAX = ILOG                    ! Initialise IMAX for this molecule
          ISUM = ( MIN(JWNO,IWNOH) - IWNOL + 1 ) * ILOG 
          DO WHILE ( JWNO .LT. IWNOH )   ! Include any further mol.data ranges
            IDX = IDX + 1
            ILOG = OLOG10(IDX)
            IMAX = MAX ( IMAX, ILOG )           
            IWNO = JWNO + 1              ! Update lower wno of new data range
            JWNO = JWNO + WAVENO(IDX)    ! & Upper wno
            ISUM = ISUM + ( MIN(JWNO,IWNOH) - IWNO + 1 ) * ILOG  ! & Summation
          END DO
          OPT = FLOAT(IMAX) + FACTOR * FLOAT(ISUM) ! Combine Max & Avg opt.dep
          IF ( IMAX .GT. 0 .AND.                   ! Finite absorption
     &         OPT .GT. FLOAT(IMIN) ) THEN         ! Above reqd threshold
            JGAS = NGAS + 1                        ! Find position in list
            DO IGAS = NGAS, 1, -1                
              IF ( OPTLST(IGAS) .LT. OPT ) THEN    ! Found smaller absorber
                IF ( IGAS .LT. NREQ ) THEN         ! Unless at end of list ...
                  OPTLST(IGAS+1) = OPTLST(IGAS)    ! Shuffle down in list
                  MOLLST(IGAS+1) = MOLLST(IGAS)
                END IF
                JGAS = IGAS                        ! Save insertion location 
              END IF
            END DO
            IF ( JGAS .LE. NREQ ) THEN             ! If not beyond end of list
              OPTLST(JGAS) = OPT                   ! Insert current into list
              MOLLST(JGAS) = IMOL
            END IF
            IF ( NGAS .LT. NREQ ) NGAS = NGAS + 1  ! Increase number stored
          END IF                                   ! End case of sig.absorber
        END IF                                     ! End case of defined mol.
      END DO                                       ! End loop over molecules
C
      NREQ = NGAS
C
      END
