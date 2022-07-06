      SUBROUTINE INPSPC ( LUNDRV, LUNSPC, FAIL, ERRMSG )
C
C VERSION
C     01-APR-08  AD  Bug#70: use limit ISTA:IEND in argument to SPCFIL 
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Read RFM wavenumber range and resolution from driver table 
C     Called by RFMINP once following *SPC marker
C     Requested parameters checked against values specified in RFM requirements
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      LUNDRV  !  I  LUN for Driver File
      INTEGER      LUNSPC  !  I  LUN for spectral range tables
      LOGICAL      FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80 ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  NXTREC ! Load next record from RFM input file.
     &, SPCFIL ! Open file and load pretabulated spectral range data.
     &, SPCLAB ! Check Spc.range label find index in table.
     &, SPCRES ! Set calculation resolution for labelled spectral range.
     &, SPCRNG ! Check Spectral Range & Resln parameters and update table.
     &, TXTFLD ! Identify start and end points of text field in record
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'rfmfil.inc' ! Standard filenames of RFM I/O files
C
C COMMON VARIABLES
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
C
C LOCAL VARIABLES
      INTEGER      IDXSPC  ! Index of spectral range in table
      INTEGER      IEND    ! Pointer to end of field in text record
      INTEGER      ISTA    ! Pointer to start of field in text record
      LOGICAL      ENDSEC  ! Set TRUE when last record in section is found
      LOGICAL      LDUMMY  ! Set TRUE if new label found (ignored)
      CHARACTER*80 RECORD  ! Text record containing spec.range data
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C Initialise tabulated spectral ranges
C
      NSPC = 0
C
C Load next driver table record containing data
C
  100 CALL NXTREC ( LUNDRV, RECORD, ENDSEC, FAIL, ERRMSG )
      IF ( FAIL ) RETURN
C
C Find how many fields in driver table record and call appropriate routine
C 
      IF ( .NOT. ENDSEC ) THEN
        CALL TXTFLD ( 80, RECORD, 2, ISTA, IEND )
        IF ( ISTA .EQ. 0 ) THEN       ! 1 field only: assume Spc.Range filename
          CALL TXTFLD ( 80, RECORD, 1, ISTA, IEND ) ! Extract .spc filename
          CALL SPCFIL ( LUNSPC, RECORD(ISTA:IEND), FAIL, ERRMSG )  
          IF ( FAIL ) RETURN
        ELSE
          CALL TXTFLD ( 80, RECORD, 3, ISTA, IEND )              
          IF ( ISTA .EQ. 0 ) THEN     ! 2 fields - assume Label + Resln
            CALL SPCRES ( RECORD, FAIL, ERRMSG )                 ! 2 Fields
            IF ( FAIL ) RETURN
          ELSE
            CALL TXTFLD ( 80, RECORD, 4, ISTA, IEND )
            IF ( ISTA .EQ. 0 ) THEN   ! 3 fields - assume unlabelled Spc.range
              CALL SPCLAB ( ' ', IDXSPC, LDUMMY, FAIL, ERRMSG )
              IF ( FAIL ) RETURN
              CALL SPCRNG ( IDXSPC, RECORD, FAIL, ERRMSG ) 
              IF ( FAIL ) RETURN
            ELSE                      ! 4 fields - assume labelled Spc.range
              CALL TXTFLD ( 80, RECORD, 1, ISTA, IEND )
              CALL SPCLAB ( RECORD(ISTA:IEND), IDXSPC, LDUMMY, 
     &                      FAIL,ERRMSG)
              IF ( FAIL ) RETURN
              CALL SPCRNG ( IDXSPC, RECORD(IEND+1:), FAIL, ERRMSG )
              IF ( FAIL ) RETURN
            END IF
          END IF
        END IF          
        GOTO 100
      END IF
C
C Check range and resolution parameters are reasonable
C
      CALL SPCCHK ( FAIL, ERRMSG )
C
      END
