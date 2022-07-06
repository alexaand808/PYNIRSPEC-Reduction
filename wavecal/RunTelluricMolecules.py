import WaveCal_orig_Lband_NIRSPEC2 as tell
import matplotlib.pyplot as plt
import numpy as np


res = 35000
minimum = 2.82366 #3.62446
maximum = 2.87938	#3.69897

mollist = ['H2O', 'CO2']#, 'CO2']#, 'N2O', 'OCS', 'O3']
scale = [0.2,1.2]#,1]#,0.2,0.2,0.2]
cal_full = tell.WaveCal_molecule(wave_min=minimum, wave_max=maximum,resolution=res,scale=scale,mollist=mollist)
tellwave, tellflux = cal_full.WaveMod, cal_full.SkyRadMod 
	
for num, mol in enumerate(mollist):	
	plt.plot(tellwave, tellflux[num], label=mol)


#plt.plot(cal_full.FullWaveMod, cal_full.FullSkyRadMod,color='k',label='Full')
plt.legend()	
plt.show()	

