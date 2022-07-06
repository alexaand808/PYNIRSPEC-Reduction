C GRACOM.INC -----------------------------------------------------------------
C
      COMMON / GRACOM / NGRA, NOMGRA, FLPGRA, FLTGRA, FLVGRA, 
     &                  DNSGRA, PREGRA, PSIGRA, RFRGRA, TEMGRA, VMRGRA
C
C VERSION
C     30-DEC-99  AD  Original.
C
C DESCRIPTION
C     Atmospheric 2-D Field data.
C     Mostly set by PRFGRA, missing values filled by ATMGRA.
C     DNSGRA and RFRGRA set by ATMAUX.
C
      INTEGER NGRA                          ! No. of profile locations 
      INTEGER NOMGRA                        ! Index of ref (psi=0)  profile
      LOGICAL FLPGRA(MAXGRA)                ! Flag T if Press.profile set
      LOGICAL FLTGRA(MAXGRA)                ! Flag T if Temp.profile set
      LOGICAL FLVGRA(MAXGAS,MAXGRA)         ! Flag T if VMR profile set
      REAL    DNSGRA(MAXATM,MAXGRA)         ! 2-D Number Density field [/cm^3]
      REAL    PREGRA(MAXATM,MAXGRA)         ! 2-D Pressure field
      REAL    PSIGRA(MAXGRA)                ! Psi values of each location
      REAL    RFRGRA(MAXATM,MAXGRA)         ! 2-D Refractivity field
      REAL    TEMGRA(MAXATM,MAXGRA)         ! 2-D Temperature field
      REAL    VMRGRA(MAXATM,MAXGAS,MAXGRA)  ! 2-D VMR field for each species
C
C -----------------------------------------------------------------------------
