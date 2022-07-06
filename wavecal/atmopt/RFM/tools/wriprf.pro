PRO wriprf, gas, file, prof

sgas = '*' + STRUPCASE(gas) + ' '
l    = STRLEN(sgas)
nlev = N_ELEMENTS(prof)
IF gas EQ 'HGT' THEN BEGIN
   OPENW, lun, file,/get_lun
   PRINTF, lun, nlev
   PRINTF, lun, sgas
ENDIF ELSE BEGIN
   OPENW, lun, file, /get_lun, /APPEND
   PRINTF, lun, sgas
ENDELSE


PRINTF, lun, prof
CLOSE, lun
FREE_LUN, lun

END
