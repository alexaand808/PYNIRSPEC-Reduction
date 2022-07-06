pro read_data, the_file, x, y, num
; Generic IDL routine to read in two column data.  It will ignore
; lines that begin with a # sign in the first column.  Other character
; data will cause an error.  You specify the filename as the first argument
; and you give two arrays x and y which have enough elements to hold the
; entire data set.  The routine returns num = to the number of elements
; successfully read.

openr, 1, the_file

line = ''
i = 0
while not EOF(1) do begin
  readf, 1, line
  a = strmid(line, 0, 1)
  if ( a ne "#") then begin
    line = strcompress(line)
    parts = str_sep(line, ' ')
    x(i) = float(parts(0))
    y(i) = float(parts(1))
    i = i + 1
  endif
endwhile
num = i
print, 'Read in ',num,' values.'
close, 1

end ;program
