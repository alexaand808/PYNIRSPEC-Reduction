#-------------------------------------------------------------------------------------------------------
#	
#	Script to gather the header information from a directory of fits files
#   Maria Camarca
#   Oct 2020
# 
#-------------------------------------------------------------------------------------------------------
import os
import astropy.io.fits as fits
import glob
from pathlib import Path

def logfile(dir, filename = None, print_terminal = True):
    '''
    purpose: list the contents of a spectral directory in an output file
    params: 
    dir: the directory of the fits files to be read.
    filename: the name of the output file. default is "[date]-logfile.txt", where
    date is read from the second to last part of the pathname (based on current structure for Geoff Blake's
    nirspec data on angela)
    '''
    if filename is None:
        p = Path(dir)
        #this lil chunk is supposed to pull out the folder with the date and do a completely arbitrary replacement
        # of underscores with dashes (maria personal preference). just delete the .replace("_","-") if unwanted
        filename = p.parts[-2].replace("_","-") + '-logfile.txt'
    
    file_list = sorted(glob.glob(dir + '/*.fits'))
    print(file_list)
    f      = open(filename, "w")
    hd     = ['File','UTC','Target Name','Object Name','Airmass','Time','Read','Slit','Filter','Echelle / CDG']
    string = '{:<7}{:<14}{:<16}{:<20}{:<9}{:<.3f}/{:<5}{:<7}{:11}{:<15}{} / {}  {}'
    hd_str = '{:<7}{:<14}{:<16}{:<20}{:<9}{:<11}{:<7}{:<11}{:<15}{}'
    space  = '-' *110
    
    f.write('Directory: '+ dir + "\n")
    f.write(space + "\n" + hd_str.format(*hd) + "\n" + space + "\n")
    
    if print_terminal: 
        print(space + "\n" + hd_str.format(*hd) + "\n" + space)

    
    for i in file_list:
        entry = fits.getheader(i, 0)
        #flats check: mirror is in and the halogen lamp is on
        if ''.join(entry['Halogen'].split()).lower() == 'on' and ''.join(entry['CALMPOS'].split()).lower() == 'in':
            flat = 'Flat = ON  Lamp = ' + str(entry['HALOLEV'])
        elif 'flat' in entry['OBJECT'].lower(): 
            flat = entry['OBJECT']
        else: flat = ''
       
        var = [entry['FRAMENUM'], entry['UTC'], entry['TARGNAME'], entry['OBJECT'],entry['AIRMASS'], entry['TRUITIME'], entry['COADDS'],
               entry['SAMPMODE'], entry['SLITNAME'], entry['FILTER'], entry['ECHLPOS'], entry['DISPPOS'], flat]
        
        if print_terminal: 
            print(string.format(*var))
        f.write(string.format(*var) + "\n")
    
    f.close()
    print('\n' + '-----> Saving '+ filename)

filepath = input('Enter directory: ') 
logfile(filepath)