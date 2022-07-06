C RFMSIZ.INC ------------------------------------------------------------------
C 
C VERSION
C     10-JAN-13  AD  Add GEISA molecules & Bro: change MAXHLN from 47 to 56
C                    Remove redundant MAXHFS
C     01-AUG-13  AD  Changes for HITRAN 2012 compatibility:
C                    Increase MAXMOL from 85 to 160
C                    Increase MAXHLN from 46 to 47
C                    Increase MAXISO from 9 to 11
C     23-APR-12  AD  Add MAXSFC
C     14-APR-08  AD  Increase MAXMOL from 81 to 85
C     18-NOV-04  AD  Remove MAXPNT
C     23-FEB-03  AD  Increase MAXHLN from 45 to 46
C     13-FEB-03  AD  Add MAXTPS, MAXTPT
C     03-JAN-03  AD  Add MAXOPG
C     11-JUN-02  AD  Remove MAXHFL
C     12-DEC-01  AD  Increase MAXPTR from 42 to 56
C     07-FEB-01  AD  Increase MAXMOL from 69 to 81
C     06-JUN-00  AD  Add MAXLUG
C     30-DEC-99  AD  Add MAXGRA
C     07-OCT-99  AD  Add MAXILF
C     07-AUG-99  AD  Redefine MAXXSW as MAXXSP     
C     28-JUL-99  RJW Replace MAX function in PARAMETER statement for MAXPNT
C     28-APR-99  AD  Remove MAXLUY (redundant)
C     02-FEB-99  AD  Add MAXOUT. Remove MAXPTB
C     03-JAN-99  AD  Add MAXJAC, MAXPTB
C     21-DEC-98  AD  Remove MAXFIL (redundant), 
C                    Set all MAXLU* arrays to 1 as default
C                    Increase MAXMOL from 63 to 69
C     15-OCT-98  AD  Comment regarding MAXLFL relative to MAXLUT
C     18-JUN-98  AD  Add MAXLUY
C     10-DEC-97  AD  Add MAXGFL, initially set to 10
C                    Add MAXGRD, initially set to 1000                      
C     21-JUL-97  AD  Add MAXBND, initially set to 5
C     02-JUL-97  AD  Add MAXLFL, MAXLUL, MAXLUT, MAXLUX, MAXLUV
C     03-MAR-97  AD  Version 3.
C     01-OCT-96  AD  Version 2.
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     Parameters defining RFM array sizes
C     Apart from the parameters which depend on other parameters (also indicated
C     by the comments "don't change this!") all of these values can be safely
C     altered to any positive integer (1,2,...) - the RFM code will check that
C     enough array space is allocated before using the arrays.
C     MAXCLC, MAXWID and MAXFIN are usually the critical parameters in defining 
C     total RFM memory requirements.
C
      INTEGER MAXATM               ! Max no. of profile data points 
        PARAMETER ( MAXATM = 521 ) ! 
      INTEGER MAXBND               ! Max no. of selected vibration bands 
        PARAMETER ( MAXBND = 5 )   ! 
C
C Note: MAXCLC is a major factor in determining total memory (MAXCLC*MAXWID and
C MAXCLC*MAXFIN). Its size determines the total number of path segments that 
C can be calculated within a single run of the RFM (product of the number of
C tangent heights, incl.FOV convolution, gases, and atmospheric layers). You 
C probably want to set this as large as your machine can handle.
      INTEGER MAXCLC               ! Max no. of paths with line-by-line calc.
        PARAMETER ( MAXCLC = 10000 ) ! 1260 500 
C
      INTEGER MAXCYC               ! Max no. of lines stored in cyclic buffer
        PARAMETER ( MAXCYC = 10000 )! 10000 adequate for standard HITRAN 96

C Note: MAXFIN is a major factor in determining total memory (MAXFIN*MAXCLC)
C Its size determines the maximum resolution (pts/wavenumber) that can be 
C calculated with the RFM. Probably best to change MAXCLC or MAXWID first.
      INTEGER MAXFIN               ! Max no. of fine mesh grid pts per widemesh
        PARAMETER ( MAXFIN = 2016 ) ! 2000 required to resolve Doppler widths 
C                                   ! in the middle atmosphere
      INTEGER MAXFOV               ! Max no. of FOV tabulation points
        PARAMETER ( MAXFOV = 21 )  ! 
      INTEGER MAXGAS               ! Max no. gases read from line/xsec dataset
        PARAMETER ( MAXGAS = 25 )
      INTEGER MAXGFL               ! Max no. of GRD input files
        PARAMETER ( MAXGFL = 2 )  ! 
      INTEGER MAXGRA               ! Max no. of horizontal gradient profiles
        PARAMETER ( MAXGRA = 3 )   
      INTEGER MAXGRD               ! Max no. of Irreg. GRD points (.GE. 2)
        PARAMETER ( MAXGRD = 5000 )   ! 2 is the minimum value allowed
C note: if MAXHLN > 56 no of fwd ptr blocks in .hit file has to be increased
      INTEGER MAXHLN               ! No. species with line data
        PARAMETER ( MAXHLN = 56 )  ! 56 is highest defined species

      INTEGER MAXILF               ! Max.no fine-grid pts in either wing of ILS
        PARAMETER ( MAXILF = 2501 )! 2501 allows +/- 1.0cm-1 @ 0.0004cm-1 spacing
      INTEGER MAXILS               ! Max.no ILS functions that can be stored
        PARAMETER ( MAXILS = 6 )   ! 6 = 1 for each passband, plus 1 general
      INTEGER MAXILP               ! Max.no tabulation points for ILS functions
        PARAMETER ( MAXILP = 5001 )! 5001 allows +/- 1.0cm-1 @ 0.0004cm-1 spacing
      INTEGER MAXISO               ! Highest molecular isotope id in use [9]
        PARAMETER ( MAXISO = 11 )   ! Dimension of DATA statement in GASISO.
C
      INTEGER MAXJAC               ! Max. no. Jacobian retrieval elements
        PARAMETER ( MAXJAC = 300 )  ! 37 allows for 2 species 5,8...56 km +1
C
C Note: the MAXL** dimensions refer to arrays which are only used when running
C the RFM with Look-Up Tables (LUT flag) and may be reduced to 1 otherwise
      INTEGER MAXLFL               ! Max no. LUT files used in single RFM run
        PARAMETER ( MAXLFL = 10 )   ! Separate LUN reserved for each.
      INTEGER MAXLUG               ! Max no. Irreg.Grid points in TAB input files
        PARAMETER ( MAXLUG = 100)   !
      INTEGER MAXLUL               ! Max no. basis vectors in any LUT
        PARAMETER ( MAXLUL = 30 )   !
      INTEGER MAXLUT               ! Max no. LUTs used in single spectral range
        PARAMETER ( MAXLUT = 10 )   ! Should be .LE. MAXLFL
      INTEGER MAXLUV               ! Max no. LUT Wno pts in single spec.range
        PARAMETER ( MAXLUV = 1 )   !
      INTEGER MAXLUX               ! Max no. Cmp.LUT p,T pts in single spc.range
        PARAMETER ( MAXLUX = 150 )   !
C
      INTEGER MAXMOL               ! Highest molecular id in use [81]
        PARAMETER ( MAXMOL = 160 )  ! NB: Used as dimension for DATA statements
      INTEGER MAXNTE               ! Max no. of non-LTE vibrational Tem.profiles
        PARAMETER ( MAXNTE = 500 )
C
      INTEGER MAXOPG               ! Max no.output points for each irr.grd pt
        PARAMETER ( MAXOPG = 19 )  ! Set to 1 unless using GRD + ILS/AVG options
      INTEGER MAXOUT
        PARAMETER ( MAXOUT = 3000 ) ! (approx) Max no.output files written
C
C Note: MAXPTH refers to total paths as opposed to calculated paths (MAXCLC)
C Using the SCA flag typically means only 10% of the total are calculated 
      INTEGER MAXPTH               ! Max no. of paths to be used
        PARAMETER ( MAXPTH = 120000 )  
C
      INTEGER MAXPTR               ! Max no.Forward Pointers in HITLIN file
        PARAMETER ( MAXPTR = 56 )  ! 56 allows for up to 4 records of 14
      INTEGER MAXQFN               ! Max no. of non-LTE Partition profs.stored
        PARAMETER ( MAXQFN = 10 )  ! ie maximum no.usable from any one NTE file.
      INTEGER MAXSFC               ! Max no. tab.points for surface emissivity
        PARAMETER ( MAXSFC = 2 ) 
      INTEGER MAXSPC               ! Max no. of spectral ranges stored 
        PARAMETER ( MAXSPC = 30 )  
      INTEGER MAXTAN               ! Max no. of tangent points to be calculated
        PARAMETER ( MAXTAN = 4000 )  
C
C The following can be set to 1 unless the TPS flag is used
      INTEGER MAXTPS               ! Max no. of external TIPS tabulations
        PARAMETER ( MAXTPS = 10 ) 
      INTEGER MAXTPT               ! Max no. temperatures in TIPS tabulations
        PARAMETER ( MAXTPT = 2931 )
C
C Note: MAXWID is a major factor in determining total memory (MAXWID*MAXCLC)
C Its size determines the widest spectrum (in wavenumbers) that can be 
C calculated with a single run of the RFM. 
      INTEGER MAXWID               ! Max. no. of widemesh intervals
        PARAMETER ( MAXWID = 1900 ) ! 1900 = 600-2500 in steps of 1.0
      INTEGER MAXWD2               ! 2 * MAXWID + 1 (+1 allows quick loading) 
        PARAMETER ( MAXWD2 = 2 * MAXWID + 1 ) ! Don't change this!
C
C Array sizes associated with .xsc files
      INTEGER MAXXFL               ! Max no. of X/S line data input files
        PARAMETER ( MAXXFL = 10 )   ! 
      INTEGER MAXXSC               ! Max.no. X/S data sets to be used
        PARAMETER ( MAXXSC = 20 )
      INTEGER MAXXST               ! Max no. of X/S p,T tables per molecule
        PARAMETER ( MAXXST = 62 )  ! 55 is max.no sets in HITRAN2K X/S data 
      INTEGER MAXXSP               ! Max no. of X/S waveno.pts for any molecule
        PARAMETER ( MAXXSP = 3500000 ) ! 400k req for F11 in HITRAN 96 X/S
C
C -----------------------------------------------------------------------------
