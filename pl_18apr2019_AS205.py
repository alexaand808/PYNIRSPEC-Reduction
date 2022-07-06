#-------------------------------------------------------------------------------------------------------
#	PYNIRSPEC REDUCTION PIPELINE -- Initialization File
#	July 5, 2022
#-------------------------------------------------------------------------------------------------------

import matplotlib.pylab as plt
import sys
#sys.path.insert(0, '/home/cross/nirspec/pipeline/reduction/')
import pynirspec.pynirspec_python3 as pn

#-------------------------------------------------------------------------------------------------------
#	Enter input/output locations:
#-------------------------------------------------------------------------------------------------------
path = '/export/nobackup1/nirspec/data/2019_04_18/spec/'
ut_date = '2019_04_18'
output_path = '/home/cross/nirspec/pipeline/reduction/'


#-------------------------------------------------------------------------------------------------------
#	Enter data ranges for flats and darks:
#-------------------------------------------------------------------------------------------------------
Flat_M_COIce	= [(158,168)] 		# CO_ice 61.16/36.94
Flat_M_COP		= [(147,157)]		# CO_p 63.35/37.19
FlatDark    	= [(182,186)]  
ObsDark     	= [(200,204)]
#ObsDark = FlatDark 


#-------------------------------------------------------------------------------------------------------
#	AS205 -- CO ICE
#-------------------------------------------------------------------------------------------------------

# Enter the file numbers for the disk (science) and star (standard) observations.
# Choose whether to reduce individual nod pairs or all the frames coadded together.
# If doing the star, change the sci ranges to the star's, and make the std None.
SciRanges 		= [(39,46)] 					
#SciRanges   = [(k, k+1) for k in range(39,46,2)]		
StdRanges 		= [(96,99)]


NRanges = len(SciRanges)
for i in range(0,NRanges):
	Red = pn.Reduction(flat_range=Flat_M_COIce,flat_dark_range=FlatDark, dark_range=ObsDark,
		sci_range=SciRanges[i], std_range=StdRanges, path=path, output_path=output_path, base='190418',
		shift=0.0, dtau=-0.0, level1=True, level2=False, SettingsFile='nirspec_18apr2019_as205.ini', ut_date=ut_date,
		sci_tname='AS205', std_tname='HR5984', nirspec=2, hold_plots=False)
assert(1==0)




#-------------------------------------------------------------------------------------------------------
#	HR5984 -- CO ICE
#-------------------------------------------------------------------------------------------------------

SciRanges = [(96,99)]

## use 'nirspec_18apr2019_hr5984_am2.ini' for second set
NRanges = len(SciRanges)
for i in range(0,NRanges):
	Red = pn.Reduction(flat_range=Flat_M_COIce, flat_dark_range=FlatDark, dark_range=ObsDark,
			   sci_range=SciRanges[i], std_range=None, path=path, output_path=output_path, base='190418',
			   shift=0.0, dtau=-0.0,level1=True,level2=False, SettingsFile='nirspec_18apr2019_as205.ini',ut_date=ut_date,
			   sci_tname='HR5984', std_tname='', nirspec=2, hold_plots=False)
assert(1==0)







#-------------------------------------------------------------------------------------------------------
#	AS205 -- CO P
#-------------------------------------------------------------------------------------------------------

# Disk and star ranges
# ==================================
#SciRanges   = [(k, k+1) for k in range(23,38,2)]	
SciRanges = [(23,38)]
StdRanges = [(69,72)]
#settings_file = 'AS205-M-co-p.ini'

NRanges = len(SciRanges)
for i in range(0,NRanges):
	Red = pn.Reduction(flat_range=Flat_M_COP, flat_dark_range=FlatDark, dark_range=ObsDark,
			   sci_range=SciRanges[i], std_range=StdRanges, path=path, output_path=output_path, base='190418',
			   shift=0.0, dtau=-0.0,level1=False,level2=True, SettingsFile='nirspec_18apr2019_as205_B.ini',ut_date=ut_date,
			   sci_tname='AS205', std_tname='HR5984', nirspec=2, hold_plots=False)

assert(1==0)


#-------------------------------------------------------------------------------------------------------
#	HR5984 - standard -- CO P
#-------------------------------------------------------------------------------------------------------

#SciRanges = [(69, 70), (71, 72)]
SciRanges = [(69,72)]

NRanges = len(SciRanges)
for i in range(0,NRanges):
	Red = pn.Reduction(flat_range=Flat_M_COP, flat_dark_range=FlatDark, dark_range=ObsDark,
			   sci_range=SciRanges[i], std_range=None, path=path, output_path=output_path, base='190418',
			   shift=0.0, dtau=-0.0,level1=True,level2=False, SettingsFile='nirspec_18apr2019_hr5984_2_am2.ini',ut_date=ut_date,
			   sci_tname='HR5984', std_tname='', nirspec=2, hold_plots=False)
assert(1==0)






