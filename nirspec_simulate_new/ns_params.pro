;**********************************************************
; NAME:                      NS_PARAMS
;
; PURPOSE:
; NS_PARAMS
;     This file contains the NIRSPEC parameters that are needed
;     by both the EFS and the DRP.  It is designed to be included
;     directly into an IDL *.pro file.  If parameters like grating
;     groove densities change, then this is the file that should
;     be updated.
;
; MODIFICATION HISTORY:
;
;     Started: 29 OCT 1998  - By James E. Larkin
;     modified 09 DEC 2000  - JEL modified integration times
;                             for flat field lamps
;     modified 26 Dec 2001  - GMH adjusted flat lamp times
;     modified 24 Sep 2009  - jlyke AO pupil stop replaces
;                             BrGamma filter in FW2
;     modified 16 Jul 2012  - gdoppmann,  OPEN replaces 0.576x24 slit
;
;     modified Nov 2018     - GDoppmann,
;                             PINHOLE replaces 0.720x24 slit
;                             pix, npix set to H2RG params
;                             xarc, yarc made smaller by factor of 0.666
;**********************************************************


;--------------------------------
; Set up echelle parameters
;--------------------------------
echsigma = 1./(23.21*1.00396)   ; echelle groove density
echgamma = (5.0/180.0)*!PI      ; the QLM angle
echalpha = 63.0*!dtor           ; nominal echelle angle
;--------------------------------
; Some cross disperser parameters
;--------------------------------
grasigma = 1./(75.465*1.00396)  ; cross disperser groove density
grahaang = 25.00*!pi/180.       ; the half opening angle
grablz = 10.0*!dtor
;--------------------------------
; The TMA 
;--------------------------------
few = 465.                       ; effective focal length in slit width
feh = 406.                       ; effective focal length in slit height
dist_few = fltarr(4)
dist_feh = fltarr(4)
dist_few = [486.682, 444.050, 451.472, 494.105] ; effective fl in slit wid. for corners
dist_feh = [421.442, 404.167, 405.943, 423.218] ; effective fl in slit wid. for corners
;dist_few = [488, 445, 447, 488] ; effective fl in slit wid. for corners
;dist_feh = [423, 401, 403, 424] ; effective fl in slit wid. for corners

;--------------------------------
; The Detector
;--------------------------------
; MODIFIED - GWD  Nov 2018
;            
;            new YARC and XARC values have a factor of 0.666 applied for
;            the finer plate scale

;pix = 0.027			; pixel size in mm
;npix=1024			; size of array
;yarc = 0.2			; arcseconds per y pixel in high or lowres
;xarc = 0.16			; arcseconds per x pixel in high res
pix = 0.018			; pixel size in mm
npix=2048			; size of array
yarc = 0.133			; arcseconds per y pixel in high or lowres
xarc = 0.107   		        ; arcseconds per x pixel in high res
;chipang = 4.657*!dtor           ; Rotation angle of chip.
;chipang = 4.657                 ; Rotation angle of chip.
chipang = 5.079                 ; Rotation angle of chip.
dx0 = 0                         ; detector offset in x 
dy0 = 0                         ; detector offset in y
satval = 32000                  ; saturation 
;---------------------
; slit parameters
srot0 = 0.                      ; slit rotation angle
fcol = 1200                     ; focal length mm
slitscale = 0.4854              ; mm/arcsec for slit
;----------------------------
; The filter wheels
numfil = 19 
A = {filters, name: 'J', l1:1.7, l2:2.0, graorder:3, aspect: 0.8, $
     yaspect: 0.5, wheel: 1, position: 1, low_flat: 1.0, high_flat: 1.0 }
filterarr = REPLICATE({filters},numfil)
filterarr(0) = {filters, 'NIRSPEC-1', 0.000947, 0.001121, 4, 0.9, 0.9, 2, 0, $
                29,59}
filterarr(1) = {filters, 'NIRSPEC-2', 0.001089, 0.001293, 4, 1.1, 1.1, 2, 1, $
               29,59}
filterarr(2) = {filters, 'NIRSPEC-3(J)', 0.001143, 0.001375, 3, 1.1, 1.0,2,2, $
               1.9, 39}
filterarr(3) = {filters, 'NIRSPEC-4', 0.001241, 0.001593, 3, 1.1, 0.85, 2, 3, $
               1.3, 17}
filterarr(4) = {filters, 'NIRSPEC-5(H)', 0.001431, 0.001808, 3, 1.1,0.8,2, $
               4, 1.5, 17}
;filterarr(5) = {filters, 'NIRSPEC-5(2nd)', 0.001431, 0.001808, 2,1.1, 0.8, $
;               2, 4, 2, 20}
filterarr(5) = {filters, 'NIRSPEC-6', 0.001558, 0.002315, 2, 1.0, 0.6, 2, 5, $
               0.25, 1.5}
filterarr(6) = {filters, 'NIRSPEC-7', 0.001839, 0.002630, 2, 1.3, 0.6, 2, 6, $
               0.25, 3.5}
;filterarr(7) = {filters, 'BR-GAMMA(2.16)', 0.002155, 0.002175,$ 
;                2,35.0,22.0,2,7,0.25,2}
filterarr(7) = {filters, 'CO(2.3)', 0.0022810, 0.002305, 2, 37.0, 20, 2,8, $
               0.25, 2}
filterarr(8) = {filters, 'K', 0.001996, 0.002382, 2, 2.5, 1.2, 2, 9, $
               0.25, 3.5}
filterarr(9) = {filters, 'K-PRIME', 0.00195, 0.002295, 2, 2.5, 1.2, 2, 10, $
               0.25, 3.5}
filterarr(10) = {filters, 'L-PRIME', 0.00342, 0.00412, 1, 4.0, 1.3, 1, 3, $
               0.25, 1}
filterarr(11) = {filters, 'M-PRIME', 0.00457, 0.00481, 1, 14, 4, 1, 4, $
                0.25, 1}
filterarr(12) = {filters, 'KL', 0.002134, 0.004228, 1, 0.8, 0.4, 1, 5, $
                0.25, 1}
filterarr(13) = {filters, 'HeI(1.08)', 0.0010776, 0.0010884, 4, 20.0,20.0,1,6,$
                10, 120}
filterarr(14) = {filters, 'PA-BETA(1.28)', 0.0012757, 0.0012888,4,18.0,18,1,7,$
                10, 120}
filterarr(15) = {filters, 'FeII(1.64)', 0.001639, 0.001654, 3, 18., 18., 1, 8,$
                1, 12}
filterarr(16) = {filters, 'H2(2.12)', 0.002110, 0.002129, 2, 35., 22., 1, 9,$
                0.25, 2}
filterarr(17) = {filters, 'M-WIDE', 0.004420, 0.005530, 1, 3.0, 0.75, 1, 10,$
                0.25, 1}
filterarr(18) = {filters, 'BLANK', 1, 2, 0, 1, 1, 1, 11, 1, 1}

;----------------------------
; Slit wheel
numslit = 12
A = {slits, name: '0.48x12', length: 12.0, width: 0.48, clock:0.1945}
slitlist = REPLICATE({slits},numslit)
slitlist(0) = {slits, '0.144x12', 12.0, 0.144, 0.24435}
slitlist(1) = {slits, '0.288x12', 12.0, 0.288, 0.23562}
slitlist(2) = {slits, '0.432x12', 12.0, 0.432, 0.25307}
slitlist(3) = {slits, '0.576x12', 12.0, 0.576, 0.24435}
slitlist(4) = {slits, '0.720x12', 12.0, 0.720, 0.23562}
slitlist(5) = {slits, '0.288x24', 24.0, 0.288, 0.17977}
slitlist(6) = {slits, '0.432x24', 24.0, 0.432, 0.23911}
;slitlist(7) = {slits, '0.576x24', 24.0, 0.576, 0.23562}
slitlist(7) = {slits, 'OPEN', 24.0, 0.576, 0.23562}
;slitlist(8) = {slits, '0.720x24', 24.0, 0.720, 0.24435}
slitlist(8) = {slits, 'PINHOLE', 24.0, 0.720, 0.24435}
slitlist(9) = {slits, '42x0.380', 0.38, 42.0, -0.04712}
slitlist(10) = {slits, '42x0.570', 0.57, 42.0, -0.04887}
slitlist(11) = {slits, '42x0.760', 0.76, 42.0, -0.04887}

;----------------------------
; Global magnification
xmag = 1.0355
ymag = 1.0086

