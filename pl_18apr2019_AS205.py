#-------------------------------------------------------------------------------------------------------
#	PYNIRSPEC REDUCTION PIPELINE -- Initialization File
#	July 5, 2022
#-------------------------------------------------------------------------------------------------------

from asyncio import base_tasks
import matplotlib.pylab as plt
import sys
sys.path.insert(0, '/home/cross/nirspec/pipeline/reduction/')
import pynirspec.pynirspec_python3 as pn

#-------------------------------------------------------------------------------------------------------
#	Enter input/output locations:
#-------------------------------------------------------------------------------------------------------
path = '/export/nobackup1/nirspec/data/2019_04_18/spec/'
base 		= '190418'
ut_date 	= '2019_04_18'
output_path = '/home/cross/nirspec/pipeline/reduction/'


#-------------------------------------------------------------------------------------------------------
#	Enter data ranges for flats and darks:
#-------------------------------------------------------------------------------------------------------
Flat_M_COice	= [(158,168)] 		# CO_ice 61.16/36.94
Flat_M_COp		= [(147,157)]		# CO_p 63.35/37.19
FlatDark    	= [(182,186)]  
ObsDark     	= [(200,204)]
#ObsDark = FlatDark 


#-------------------------------------------------------------------------------------------------------
#	AS205 -- CO ICE
#-------------------------------------------------------------------------------------------------------
# Enter the file numbers for the disk (science) and star (standard) observations.
# Choose whether to reduce individual nod pairs or all the frames coadded together.
# If doing the star, change the sci ranges to the star's, and make the std None.
Flat_MWide 		= Flat_M_COice
SettingsFile	= 'inifiles/nirspec_18apr2019_as205.ini'
sci_tname		= 'AS205_COice'
std_tname		= 'HR5984'

SciRanges 		= [(39,46)] 					
#SciRanges   = [(k, k+1) for k in range(39,46,2)]		
StdRanges 		= [(96,99)]


NRanges = len(SciRanges)
for i in range(0,NRanges):
	Red = pn.Reduction(flat_range=Flat_MWide,flat_dark_range=FlatDark, dark_range=ObsDark,
		sci_range=SciRanges[i], std_range=StdRanges, path=path, output_path=output_path, base=base,
		shift=0.0, dtau=-0.0, level1=True, level2=False, SettingsFile=SettingsFile, ut_date=ut_date,
		sci_tname=sci_tname, std_tname=std_tname, nirspec=2, hold_plots=False)
assert(1==0)




#-------------------------------------------------------------------------------------------------------
#	HR5984 - standard -- CO ICE
#-------------------------------------------------------------------------------------------------------

SciRanges = [(96,99)]

## use 'nirspec_18apr2019_hr5984_am2.ini' for second set
NRanges = len(SciRanges)
for i in range(0,NRanges):
	Red = pn.Reduction(flat_range=Flat_MWide, flat_dark_range=FlatDark, dark_range=ObsDark,
			   sci_range=SciRanges[i], std_range=None, path=path, output_path=output_path, base=base,
			   shift=0.0, dtau=-0.0,level1=True,level2=False, SettingsFile=SettingsFile,ut_date=ut_date,
			   sci_tname='HR5984', std_tname='', nirspec=2, hold_plots=False)
assert(1==0)







#-------------------------------------------------------------------------------------------------------
#	AS205 -- CO P
#-------------------------------------------------------------------------------------------------------
Flat_MWide 		= Flat_M_COp
SettingsFile	= 'inifiles/nirspec_18apr2019_as205_B.ini'
sci_tname		= 'AS205_COp'
std_tname		= 'HR5984'

# Disk and star ranges
SciRanges = [(23,38)]
#SciRanges   = [(k, k+1) for k in range(23,38,2)]	
StdRanges = [(69,72)]

NRanges = len(SciRanges)
for i in range(0,NRanges):
	Red = pn.Reduction(flat_range=Flat_MWide, flat_dark_range=FlatDark, dark_range=ObsDark,
			   sci_range=SciRanges[i], std_range=StdRanges, path=path, output_path=output_path, base=base,
			   shift=0.0, dtau=-0.0,level1=True,level2=False, SettingsFile=SettingsFile,ut_date=ut_date,
			   sci_tname=sci_tname, std_tname=std_tname, nirspec=2, hold_plots=False)

assert(1==0)


#-------------------------------------------------------------------------------------------------------
#	HR5984 - standard -- CO P
#-------------------------------------------------------------------------------------------------------

#SciRanges = [(69, 70), (71, 72)]
SciRanges = [(69,72)]

NRanges = len(SciRanges)
for i in range(0,NRanges):
	Red = pn.Reduction(flat_range=Flat_MWide, flat_dark_range=FlatDark, dark_range=ObsDark,
			   sci_range=SciRanges[i], std_range=None, path=path, output_path=output_path, base=base,
			   shift=0.0, dtau=-0.0,level1=True,level2=False, SettingsFile=SettingsFile,ut_date=ut_date,
			   sci_tname='HR5984', std_tname='', nirspec=2, hold_plots=False)
assert(1==0)






