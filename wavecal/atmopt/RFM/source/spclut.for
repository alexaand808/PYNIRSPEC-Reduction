      SUBROUTINE SPCLUT ( ISPC, FAIL, ERRMSG )
C
C VERSION
C     06-APR-06  AD  Change 'X' to '1X' in format
C     09-JUN-05  AD  Remove test for old format LUT files.
C     09-MAY-04  AD  Allow for isotopic LUTs
C     02-JAN-03  AD  Add outcom.inc since *OUT variables now set beforehand
C     08-JUN-00  AD  Allow for Irregular Grid LUT files
C     27-APR-00  AD  Convert ILS to DP
C     05-JAN-00  AD  Allow for AVG option.
C     11-AUG-99  AD  Change variable names and make D.P. explicit
C                    Initialise and SAVE SHPSAV
C     04-JUN-99  AD  Allow for binary LUT files as well. 
C     24-MAY-99  AD  Set DV,DP,DT=1 if NP,NV,NT=1
C     28-APR-99  AD  Allow for NDPLFL, NDTLFL: -lnp,Tem axis increment factors
C                    Also set OFFLUT 
C     23-APR-98  AD  Adapt error messages for C*8 LABSPC and MWCODE
C                    Also check for C*6 MWCODE on input file
C     19-JUN-98  AD  Separate TAB.LUT and Cmp.LUT storage.
C     20-JAN-98  AD  Read number of singular vectors determined by NSVLFL
C     17-DEC-97  AD  Correction to calculation of WNOMIN with ILS function
C     01-NOV-97  AD  Set NLUT, LUNLUT
C     03-JUL-97  AD  Original.
C
C DESCRIPTION
C     Initialise LUT data for each new spectral range.
C     Called RFMSPC once for each spectral range.
C     Note that all array dimensions have been checked previously by LUTCHK.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER       ISPC    !  I  Current spectral range
      LOGICAL       FAIL    !  O  Set TRUE if a fatal error is detected
      CHARACTER*80  ERRMSG  !  O  Error message returned if FAIL is TRUE
C
      EXTERNAL
     &  LUTGRD ! Read coded irregular grid data and turn into integer
     &, LUTRNG ! Check p,T range of Look-Up Table data.
     &, RFMLOG ! Write text message to RFM log file.
C
C GLOBAL CONSTANTS
      INCLUDE 'rfmsiz.inc' ! RFM Array sizes
      INCLUDE 'shpcon.inc' ! Line-shape codes
C
C COMMON VARIABLES
      INCLUDE 'flgcom.inc' ! RFM option flags
      INCLUDE 'gascom.inc' ! Molecular and isotope data
      INCLUDE 'ilscom.inc' ! Instrument Lineshape functions.
      INCLUDE 'lflcom.inc' ! Look-Up Table file data
      INCLUDE 'lutcom.inc' ! Look-Up Table data
      INCLUDE 'outcom.inc' ! RFM output file data
      INCLUDE 'spccom.inc' ! Spectral Range & Resolution data
C
C LOCAL VARIABLES
      INTEGER IDG            ! HITRAN ID of gas, read from LUT file 
      INTEGER IDI            ! Isotopic ID of gas, read from LUT file
      INTEGER IG             ! Index for Wno dimension of TAB files
      INTEGER IG0            ! Pointer for Wno dimension offset
      INTEGER IGAS           ! Index of gas in LUT file
      INTEGER IHEAD          ! Counter for header records to be skipped
      INTEGER IILS           ! Index of ILS data
      INTEGER IL             ! Basis vector counter
      INTEGER ILFL           ! LUT file counter
      INTEGER IMOL           ! HITRAN ID of molecule
      INTEGER IOS            ! Saved value of IOSTAT for error messages
      INTEGER IV             ! Relative pointer for Wno dimension 
      INTEGER IVMAX,IVMIN    ! Max/Min values of IV index req to span Range+ILS
      INTEGER IV0            ! Pointer for Wno dimension offset
      INTEGER IX             ! Counter for p,T dimension
      INTEGER IX0            ! Pointer for Cmp. p,T dimension offset
      INTEGER IY0            ! Pointer for Tab.file p,T dimension offset
      INTEGER L              ! Length of molecule name
      INTEGER LUN            ! LUN for LUT data
      INTEGER ILUT           ! Index of LUTs within Spc.range 
      INTEGER NHEAD          ! No.of header records to be skipped
      INTEGER NSPACE         ! No.of free array elements remaining
      INTEGER NSV            ! Number of singular vectors required from file
      INTEGER NX             ! Number of p,T points in LUT
      INTEGER SHPSAV(MAXGAS) ! Saved values of SHPGAS
      LOGICAL BINFIL         ! T=Binary file, F=ASCII file
      LOGICAL IRRFIL         ! T=Irreg.Grid file, F=Full grid file
      REAL    RDUMMY         ! Dummy variable for skipping excluded SV data
      CHARACTER*8  MWCODE    ! MW code read from LUT file
      CHARACTER*80 RECORD    ! Record read from LUT file
      CHARACTER*3  CDUMMY    ! Tabulated function
      CHARACTER*80 MESSGE    ! Message for RFM log file
      DOUBLE PRECISION WNOMAX ! Highest Wno [/cm] required from LUT 
      DOUBLE PRECISION WNOMIN !  Lowest Wno [/cm] required from LUT 
C
      DATA SHPSAV / MAXGAS * 0 /       ! Purely to satisfy ftnchek
      SAVE SHPSAV
C
C EXECUTABLE CODE -------------------------------------------------------------
C
C For each spectral range, any gas using LUT data has SHPGAS overwritten with 
C info identifying LUT, so this needs to be replaced from previous Spc.range
C
      DO ILFL = 1, NLFL
        IF ( SPCLFL(ILFL) .EQ. ISPC-1 ) THEN
          IGAS = GASLFL(ILFL)
          SHPGAS(IGAS) = SHPSAV(IGAS)
        END IF
      END DO
C
      IG0 = 0      
      IV0 = 0      
      IX0 = 0
      IY0 = 0
      ILUT = 0
C
      DO ILFL = 1, NLFL
        IF ( SPCLFL(ILFL) .EQ. ISPC ) THEN
          IGAS = GASLFL(ILFL)
          IMOL = IDXGAS(IGAS)
          ILUT = ILUT + 1
          SHPSAV(IGAS) = SHPGAS(IGAS) 
          SHPGAS(IGAS) = SHPLUT + ILUT
          LUN = LUNLFL(ILFL)
          TABLUT(ILUT) = TABLFL(ILFL)
          NSV = NSVLFL(ILFL)
          NHEAD = NHDLFL(ILFL)
          BINFIL = BINLFL(ILFL)
          BINLUT(ILUT) = BINFIL
          IRRFIL = IRRLFL(ILFL)
          IRRLUT(ILUT) = IRRFIL
          REWIND ( LUN, IOSTAT=IOS, ERR=900 )

          IF ( BINFIL ) THEN
            DO IHEAD = 1, NHEAD
              READ ( LUN, ERR=900, IOSTAT=IOS ) RECORD(1:1)
            END DO
            READ ( LUN, ERR=900, IOSTAT=IOS ) RECORD(1:15)
          ELSE
            DO IHEAD = 1, NHEAD
              READ ( LUN, '(A)', ERR=900, IOSTAT=IOS ) RECORD(1:1)
            END DO
            READ ( LUN, '(A)', ERR=900, IOSTAT=IOS ) RECORD
          END IF
          IF ( ISOMOL(IMOL) ) THEN
            READ (RECORD, '(A8,1X,I2,1X,I1,1X,A3)', IOSTAT=IOS, ERR=900) 
     &        MWCODE, IDG, IDI, CDUMMY
            L = INDEX ( CODGAS(IGAS)//' ', ' ' ) -1
            WRITE ( MESSGE, '(A,A,A,A,A,A,A,I1)' )
     &        'I-SPCLUT: Spc.Rng=', LABSPC(ISPC),
     &        ' Loading LUT data, MWCODE=', MWCODE,
     &        ', Gas=', CODGAS(IGAS)(1:L), ', Iso=', IDIGAS(IGAS)
            IF ( IDG .NE. IMOL ) STOP 'F-SPCLUT: Logical Error#1'
            IF ( IDI .NE. IDIGAS(IGAS) ) 
     &        STOP 'F-SPCLUT: Logical Error#2'
          ELSE
            READ ( RECORD, '(A8,1X,I2,1X,A3)', IOSTAT=IOS, ERR=900 ) 
     &        MWCODE, IDG, CDUMMY
            WRITE ( MESSGE, '(A,A,A,A,A,A)' )
     &        'I-SPCLUT: Spc.Rng=', LABSPC(ISPC),
     &        ' Loading LUT data, MWCODE=', MWCODE,
     &        ', Gas=', CODGAS(IGAS)
            IF ( IDG .NE. IMOL ) STOP 'F-SPCLUT: Logical Error#3'
          END IF
          CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
C
          IF ( BINFIL ) THEN
            READ ( LUN, IOSTAT=IOS, ERR=900 ) 
     &        NLLUT(ILUT), NVLUT(ILUT), V1LUT(ILUT), DVLUT(ILUT),
     &        NPLUT(ILUT), P1LUT(ILUT), DPLUT(ILUT),
     &        NTLUT(ILUT), T1LUT(ILUT), DTLUT(ILUT)
          ELSE
            READ ( LUN, *, IOSTAT=IOS, ERR=900 ) 
     &        NLLUT(ILUT), NVLUT(ILUT), V1LUT(ILUT), DVLUT(ILUT),
     &        NPLUT(ILUT), P1LUT(ILUT), DPLUT(ILUT),
     &        NTLUT(ILUT), T1LUT(ILUT), DTLUT(ILUT)
          END IF
          NVLUT(ILUT) = ABS ( NVLUT(ILUT) )   ! -ve value = Irreg.grid (IRRFIL)
          IF ( NVLUT(ILUT) .EQ. 1. ) DVLUT(ILUT) = 1.0D0
          IF ( NPLUT(ILUT) .EQ. 1. ) DPLUT(ILUT) = 1.0
          IF ( NTLUT(ILUT) .EQ. 1. ) DTLUT(ILUT) = 1.0
C
          CALL LUTRNG ( ILUT, IGAS, FAIL, ERRMSG )
          IF ( FAIL ) RETURN
C
          WNOMAX = WNUSPC(ISPC) 
          WNOMIN = WNLSPC(ISPC) 
          IF ( ILSFLG ) THEN
            IILS = ILSSPC(ISPC)
            WNOMIN = WNOMIN + PT1ILS(IILS) 
            WNOMAX = WNOMAX + PT2ILS(IILS)
          ELSE IF ( AVGFLG ) THEN   ! Increment by WNROUT 
            WNOMIN = WNOMIN - WNROUT
            WNOMAX = WNOMAX + WNROUT
          END IF
C
          IVMIN = 1 + NINT ( ( WNOMIN - V1LUT(ILUT) ) / DVLUT(ILUT) ) 
          IVMAX = 1 + NINT ( ( WNOMAX - V1LUT(ILUT) ) / DVLUT(ILUT) ) 
          IF ( IVMIN .LT. 1 .OR. IVMAX .GT. NVLUT(ILUT) ) THEN
            WRITE ( MESSGE, '(A,A6,A)' )
     &        'W-SPCLUT: ILS wider than LUT, Spc.Rng=', 
     &        LABSPC(ISPC), '. Assume zero absorp.outside LUT'
            CALL RFMLOG ( MESSGE, FAIL, ERRMSG )
            IF ( FAIL ) RETURN
          END IF
C
C Apply temperature offset profile if either temperature limit is negative
          OFFLUT(ILUT) = ( T1LUT(ILUT) .LT. 0.0 ) .OR.
     &      ( T1LUT(ILUT) + (NTLUT(ILUT)-1) * DTLUT(ILUT) .LT. 0.0 )
C
          IF ( NLLUT(ILUT) .GT. 0 ) THEN        ! Using LUT file (not TAB file)
            LUNLUT(ILUT) = 0
            IVOFF(ILUT) = IV0
C
C If irregular grid, read from file and load indices into IVLUT using same 
C wavenumber indexing as U matrix
            IF ( IRRFIL ) THEN
              NSPACE = MAXLUV - IV0
              CALL LUTGRD ( LUN, BINFIL, NSPACE, NVLUT(ILUT),
     &                      IVLUT(IV0+1), FAIL, ERRMSG ) 
              IF ( FAIL ) RETURN
            ELSE                                ! Full grid
              DO IV = 1, NVLUT(ILUT)
                IVLUT(IV0+IV) = IV
              END DO
            END IF

            IF ( BINFIL ) THEN
              DO IV = IV0 + 1, IV0 + NVLUT(ILUT)
                READ ( LUN, IOSTAT=IOS, ERR=900 ) 
     &            ( U(IV,IL), IL = 1, NSV ), 
     &            ( RDUMMY, IL = NSV+1, NLLUT(ILUT) )
              END DO
            ELSE
              DO IV = IV0 + 1, IV0 + NVLUT(ILUT)
                READ ( LUN, *, IOSTAT=IOS, ERR=900 ) 
     &            ( U(IV,IL), IL = 1, NSV ), 
     &            ( RDUMMY, IL = NSV+1, NLLUT(ILUT) )
              END DO
            END IF
            IV0 = IV0 + NVLUT(ILUT)
            IXOFF(ILUT) = IX0
            NX = NPLUT(ILUT) * NTLUT(ILUT)
            IF ( BINFIL ) THEN
              DO IX = IX0 + 1, IX0 + NX
                READ ( LUN, IOSTAT=IOS, ERR=900 ) 
     &            ( K(IL,IX), IL = 1, NSV ),   
     &            ( RDUMMY, IL = NSV+1, NLLUT(ILUT) )  
              END DO
            ELSE
              DO IX = IX0 + 1, IX0 + NX
                READ ( LUN, *, IOSTAT=IOS, ERR=900 ) 
     &            ( K(IL,IX), IL = 1, NSV ),   
     &            ( RDUMMY, IL = NSV+1, NLLUT(ILUT) )  
              END DO
            END IF
            IX0 = IX0 + NX
            NLLUT(ILUT) = NSV
          ELSE                                  ! Using TAB file (not LUT file)
            NGLUT(ILUT) = NVLUT(ILUT)
            IF ( IRRFIL ) THEN
              NSPACE = MAXLUG - IG0
              CALL LUTGRD ( LUN, BINFIL, NSPACE, NGLUT(ILUT),
     &                      IGLUT(IG0+1), FAIL, ERRMSG ) 
              IF ( FAIL ) RETURN
              IGOFF(ILUT) = IG0
              IG0 = IG0 + NGLUT(ILUT)
            ELSE
              DO IG = 1, NGLUT(ILUT)
                IGLUT(IG0+IG) = IG
              END DO
            END IF
            IXOFF(ILUT)  = IY0
            IY0 = IY0 + NPLUT(ILUT) * NTLUT(ILUT)
            LUNLUT(ILUT) = LUN           
            NDPLUT(ILUT) = NDPLFL(ILFL)
            NDTLUT(ILUT) = NDTLFL(ILFL)
            IRCLUT(ILUT) = 0
          END IF
C
        END IF
      END DO
C
      NLUT = ILUT
      FAIL = .FALSE.
      RETURN
C
  900 WRITE ( ERRMSG, '(A,I11)')
     &  'F-SPCLUT: I/O failure on LUT file. IOSTAT=', IOS
      FAIL = .TRUE.
C
      END
