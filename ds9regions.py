# ---------------------------------------------------------------------------------------
# Check initialization file yranges 
# July 5, 2022
# ---------------------------------------------------------------------------------------
import astropy.io.fits as fits
import numpy as np
import sys
sys.path.insert(0, '/home/cross/nirspec/pipeline/reduction/')

# Input initialization file and the locations of A and B nod files
# Copy output into terminal to draw regions on fits file in ds9

inifile = "inifiles/nirspec_17jan2022_COp.ini"

### if reduction has already run
#fits_file = "2021_11_23/FZTau_COice/PLOTS/A-B-after-bpRemoval.fits"

### if reduction has NOT yet run, uncomment this section
### ----------------------------------------
data_path = '/export/nobackup1/nirspec/data/2022_01_17/spec/'
abeam = "nspec220117_0215.fits"
bbeam = "nspec220117_0216.fits"

def diff(abeam,bbeam):
    a = fits.open(data_path+abeam)[0].data
    b = fits.open(data_path+bbeam)[0].data
    difference = np.rot90((a - b),k=3)
    newfile = abeam[:-9] + 'diff.fits'
    newfile = 'nspec_diff.fits'
    fits.writeto(newfile, difference, overwrite=True)
    #Image.fromarray(difference).save("your_file.png")
    print('saved: ' + newfile)
    return newfile
fits_file = diff(abeam,bbeam)
# ### ----------------------------------------

settingsfile = open(inifile,'r')
Lines = settingsfile.readlines()
yranges = {}
for line in Lines:
    if "yrange" in line and not "#" in line:
        #print(line)
        line=line.replace(' ','')
        line=line.split('=')
        name=line[0]
        line=line[1].replace('[','')
        line=line.replace(']','')
        ilist=line.split(',')
        ilist=map(int,ilist)
        yranges[name] = list(ilist)

print("\ncopy into command line:")
print("ds9",fits_file,"-zscale",end=" ")
# box x y width height
for onum in range(len(yranges)):
    coords = yranges["yrange"+str(onum+1)]
    y1 = coords[0]
    y2 = coords[1]
    ymid = (y1+y2)/2
    height = y2-y1
    print("-regions command 'box 1024",str(ymid),"2048",str(height),"'",end=" ")
    print("-regions command 'line 0",str(ymid),"2048",str(ymid),"'",end=" ")
print()