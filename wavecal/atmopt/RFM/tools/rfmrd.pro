PRO rfmrd, filename, wnoarr, datarr
;
; 25-AUG-03  AD  Change NWIDEMESH from INTARR to LONARR
; 20-DEC-00  VJ  Add capability to read binary files
; 28-NOV-00  AD  Original
;
; IDL procedure for reading asc or bin RFM spectra
; Binary files are identified by presence of '.bin' string in filename
;
;   filename (in)  string   name of RFM  spec. file (eg 'rad_01000.asc')
;   wnoarr   (out) dbl arr  wavenumber array
;   datarr   (out) dbl arr  data array
;

wno1  = 0.D0
wnod  = 0.D0   
nwno  = long(0)
commnt = ''

; ASCII or binary?

IF STRPOS(filename,'.bin') GT 0 THEN BEGIN    ; binary file

nwidemesh = LONARR (60000L)   ; max. no. of nwidemesh points

OPENR, 1, filename, /F77_UNFORMATTED
READU, 1, commnt
READU, 1, commnt
READU, 1, commnt
READU, 1, nwno, wno1, wnod
datarr   = DBLARR (nwno)
wnoarr = wno1 + FINDGEN(nwno) * wnod

; Read number of elements per record

i 	  = 0L
n	  = LONG(0)
WHILE (NOT EOF (1)) DO BEGIN
  READU, 1, n
;   print ,'n',n
  nwidemesh (i) = n
  i = i + 1
ENDWHILE
nrecord = i 
CLOSE, 1


; Read datarr
x = ''
i       = 0L
j       = 0L

OPENR, 1, filename, /F77_UNFORMATTED
FOR m = 0, 3 DO READU, 1, x		; skip first 4 records
FOR m = 0L, long(nrecord - 1) DO BEGIN
  recarray = FLTARR (nwidemesh (m))
  IF (NOT EOF (1)) THEN READU, 1, n, recarray
  FOR i = j, j + nwidemesh (m) - 1 DO BEGIN
    datarr (i) = recarray (i - j)
  ENDFOR
  j = j + nwidemesh (m)
ENDFOR
CLOSE, 1

ENDIF ELSE BEGIN    ; ASCII file

OPENR, 1, filename
READF, 1, commnt     ; skip header records (assume always 3)
READF, 1, commnt
READF, 1, commnt

READF, 1, nwno, wno1, wnod
WNOARR = wno1 + FINDGEN(nwno) * wnod
DATARR = DBLARR(nwno)
READF, 1, datarr
CLOSE, 1

ENDELSE
END
