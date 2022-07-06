PRO reaprf, gas, file, prof
;
; 28-NOV-00  AD
; IDL procedure for reading in profile from RFM .atm file
;
;   gas  (in)  string    name of species (eg 'CH4', 'HGT')
;   file (in)  string    name of .atm file (eg 'std.atm')
;   prof (out) real arr  profile
;
sgas = '*' + STRUPCASE(gas) + ' '
l = STRLEN(sgas)
OPENR, 1, file
header = '!'
nlev = 0
WHILE STRMID(header,0,1) EQ '!' DO READF, 1, header
READS, header, nlev
prof = FLTARR(nlev)
WHILE STRMID(header,0,l) NE sgas DO BEGIN
  READF, 1, header
  IF STRLEN(header) LT l THEN header = header + ' ' 
  header = STRUPCASE(header)
ENDWHILE
READF, 1, prof
CLOSE, 1
END

	
