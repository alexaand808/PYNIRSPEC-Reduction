C FLGCOM.INC ------------------------------------------------------------------
      COMMON / FLGCOM / ABSFLG, AVGFLG, BBTFLG, BFXFLG, BINFLG, CLCFLG, 
     &                  COOFLG, CRVFLG, CTMFLG, DBLFLG, FINFLG, FLXFLG, 
     &                  FOVFLG, FVZFLG, GEOFLG, GHZFLG, GRAFLG, GRDFLG, 
     &                  HOMFLG, ILSFLG, JACFLG, LAYFLG, LEVFLG, LINFLG, 
     &                  LOSFLG, LUNFLG, LUTFLG, MIXFLG, MTXFLG, NADFLG, 
     &                  NEWFLG, NTEFLG, OBSFLG, OPTFLG, PRFFLG, PTHFLG, 
     &                  QADFLG, RADFLG, REJFLG, RJTFLG, SFCFLG, SHPFLG, 
     &                  TABFLG, TPSFLG, TRAFLG, V42FLG, VRTFLG, VVWFLG, 
     &                  WIDFLG, WNGFLG, ZENFLG
C
C VERSION
C     08-AUG-13  AD  Add V42FLG
C     13-SEP-11  AD  Add RJTFLG
C     06-JAN-04  AD  Add PRFFLG
C     29-DEC-03  AD  Add LEVFLG
C     11-AUG-03  AD  Add VVWFLG
C     08-AUG-03  AD  Add BBTFLG, GHZFLG
C     12-FEB-03  AD  Add TPSFLG
C     31-OCT-01  AD  Add VRTFLG
C     07-JUN-01  AD  Add LOSFLG
C     26-MAY-00  AD  Add DBLFLG
C     08-APR-00  AD  Add COOFLG, FLXFLG, MTXFLG
C     04-JAN-00  AD  Add GRAFLG, AVGFLG
C     02-FEB-99  AD  Add CLCFLG, LUNFLG. Remove SCAFLG.
C     24-JAN-99  AD  Remove GL2FLG, OFMFLG, IQDFLG. 
C                    Add LINFLG, FVZFLG, QADFLG.
C     30-DEC-98  AD  Add JACFLG
C     29-DEC-98  AD  Add BFXFLG
C     13-JUL-98  AD  Add SFCFLG
C     09-DEC-97  AD  Add GRDFLG
C     21-OCT-97  AD  Add TABFLG
C     03-JUL-97  AD  Add LUTFLG
C     03-MAR-97  AD  Version 3.
C     16-JAN-97  AD  Add OBSFLG and ZENFLG
C     01-OCT-96  AD  Version 2.
C     28-SEP-96  AD  Add IQDFLG
C     20-SEP-96  AD  Add CRV and FIN flags
C                    Rename VIBFLG to NTEFLG
C     01-SEP-96  AD  Version 1.
C
C DESCRIPTION
C     RFM option flags, loaded by INPFLG.
C
      LOGICAL ABSFLG ! T = output absorption spectra
      LOGICAL AVGFLG ! T = spectrally average output
      LOGICAL BBTFLG ! T = output brightness temperatures
      LOGICAL BFXFLG ! T = use Planck Fn variation with optical depth in layer
      LOGICAL BINFLG ! T = write output spectra to binary files
      LOGICAL CLCFLG ! T = Calculate absorption coefficients for each path
      LOGICAL COOFLG ! T = Calculate cooling rates
      LOGICAL CRVFLG ! T = user-defined local radius of curvature
      LOGICAL CTMFLG ! T = use continuum data
      LOGICAL DBLFLG ! T = write output spectra with Double Precision
      LOGICAL FINFLG ! T = user-defined fine mesh prior to ILS convolution
      LOGICAL FLXFLG ! T = apply solid angle integration for flux calculation
      LOGICAL FOVFLG ! T = apply FOV convolution 
      LOGICAL FVZFLG ! T = FOV function defined as refracted tangent hts.
      LOGICAL GEOFLG ! T = use geometric (ie non-refractive) ray paths
      LOGICAL GHZFLG ! T = use GHz rather than wavenumber for spectral axis
      LOGICAL GRAFLG ! T = use horizontal gradients
      LOGICAL GRDFLG ! T = use user-supplied fine grid 
      LOGICAL HOMFLG ! T = use homogeneous path 
      LOGICAL ILSFLG ! T = use instrument lineshape convolution
      LOGICAL JACFLG ! T = calculate Jacobian spectra
      LOGICAL LAYFLG ! T = use sub-layering 
      LOGICAL LEVFLG ! T = intermediate outputs at specified altitude levels.
      LOGICAL LINFLG ! T = Assume VMR varies linearly with altitude
      LOGICAL LOSFLG ! T = calculate elevation pointing Jacobian spectra
      LOGICAL LUNFLG ! T = Use same logical unit number for all output spectra
      LOGICAL LUTFLG ! T = use Look-Up Tables for absorption coefficients
      LOGICAL MIXFLG ! T = use line-mixing
      LOGICAL MTXFLG ! T = produce inter-level matrix of flux calculations
      LOGICAL NADFLG ! T = nadir-viewing
      LOGICAL NEWFLG ! T = open new files STATUS='NEW' rather than 'UNKNOWN'
      LOGICAL NTEFLG ! T = use non-LTE vibrational temperatures 
      LOGICAL OBSFLG ! T = observer located within atmosphere
      LOGICAL OPTFLG ! T = output optical depth spectra
      LOGICAL PRFFLG ! T = output internal atmospheric profile
      LOGICAL PTHFLG ! T = output Path diagnostics
      LOGICAL QADFLG ! T = use simple quadratic fit to line wings
      LOGICAL RADFLG ! T = output Radiance spectra
      LOGICAL REJFLG ! T = implement line rejection
      LOGICAL RJTFLG ! T = output Rayleigh-Jeans brightness temperatures
      LOGICAL SFCFLG ! T = Allow for opaque surface under atmosphere
      LOGICAL SHPFLG ! T = use non-Voigt line shapes
      LOGICAL TABFLG ! T = output tabulated absorption coefficients
      LOGICAL TPSFLG ! T = use tabulated Total Internal Partition Sum data
      LOGICAL TRAFLG ! T = output Transmission spectra
      LOGICAL V42FLG ! T = reproduce RFM v4.28 behaviour
      LOGICAL VRTFLG ! T = flux-type calculations with single vertical path
      LOGICAL VVWFLG ! T = apply Van Vleck-Weisskopf correction to lineshapes
      LOGICAL WIDFLG ! T = output wide-mesh calc diagnostics
      LOGICAL WNGFLG ! T = ignore lines beyond NWDCUT Voigt half-widths
      LOGICAL ZENFLG ! T = zenith-viewing
C
C -----------------------------------------------------------------------------
