PRO chg_atm,pwv_new=pwv_new,pre_new=pre_new

atm = '../ATM/ngt.atm'
reaprf, 'HGT', atm, hgt
reaprf, 'PRE', atm, pre
reaprf, 'TEM', atm, tem
;reaprf, 'N2', atm, N2
reaprf, 'O2', atm, O2
reaprf, 'CO2', atm, CO2
reaprf, 'O3', atm, O3
reaprf, 'H2O', atm, H2O
reaprf, 'CH4', atm, CH4
 reaprf, 'N2O', atm, N2O
;reaprf, 'HNO3', atm, HNO3
reaprf, 'CO', atm, CO
;reaprf, 'NO2', atm, NO2
;reaprf, 'NO', atm, NO
;reaprf, 'NH3', atm, NH3
;reaprf, 'OCS', atm, OCS
;reaprf, 'SO2', atm, SO2

IF KEYWORD_SET(pre_new) THEN BEGIN
   pre_scl = pre_new/pre[0]
   pre = pre*pre_scl
ENDIF

low = 2.6
up  = 120.
nint = 1000.
hgt_int = (up-low)*findgen(nint)/(nint-1)+low
pre_int = INTERPOL(pre,hgt,hgt_int)
h2o_int = INTERPOL(h2o,hgt,hgt_int)
tem_int = INTERPOL(tem,hgt,hgt_int)

RR = 8.314472d7
h2o_amu = 18.01528
n_h2o = pre_int*1d3*h2o_int*1d-6/(RR*tem_int)*h2o_amu ;in g/cm^3
pwv_g_cm2 = int_tabulated(hgt_int*1d5,n_h2o,/sort, /double)
pwv_mm  = pwv_g_cm2*10.
PRINT, 'PWV [mm]: ',pwv_mm, '  from an altitude of: ',low, ' km'
PRINT, 'Pressure at ', low, ' km: ', interpol(pre,hgt,low), ' mb'

IF KEYWORD_SET(pwv_new) THEN BEGIN
   h2o_scl = pwv_new/pwv_mm
   h2o = h2o*h2o_scl
ENDIF

atm_new = '../ATM/tmp.atm'
wriprf, 'HGT', atm_new, hgt
wriprf, 'PRE', atm_new, pre
wriprf, 'TEM', atm_new, tem
;wriprf, 'N2', atm_new, N2
wriprf, 'O2', atm_new, O2
wriprf, 'CO2', atm_new, CO2
wriprf, 'O3', atm_new, O3
wriprf, 'H2O', atm_new, H2O
wriprf, 'CH4', atm_new, CH4
wriprf, 'N2O', atm_new, N2O
;wriprf, 'HNO3', atm_new, HNO3
wriprf, 'CO', atm_new, CO
;wriprf, 'NO2', atm_new, NO2
;wriprf, 'NO', atm_new, NO
;wriprf, 'NH3', atm_new, NH3
;wriprf, 'OCS', atm_new, OCS
;wriprf, 'SO2', atm_new, SO2

stop


END
