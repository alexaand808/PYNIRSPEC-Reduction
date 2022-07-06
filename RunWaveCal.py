#-------------------------------------------------------------------------------------------------------
#	Main run file for wavelength calibration
#-------------------------------------------------------------------------------------------------------

import matplotlib.pylab as plt
import tkinter
import numpy as np
import importlib	
#import Fringe_class as fringe ## 



#-------------------------------------------------------------------------------------------------------
### CHOOSE NIRSPEC ###
#-------------------------------------------------------------------------------------------------------
#nirspec = 1
nirspec = 2


#-------------------------------------------------------------------------------------------------------
### CHOOSE BAND ###
#-------------------------------------------------------------------------------------------------------
#band = 'L'
#band = 'K'
band = 'M_COice'
#band = 'M_COp'

if nirspec == 1:
	import wavecal.WaveCal_orig_Lband as WaveCal

	n_order_K = 6
	n_order_L = 5
	n_order_M = 2

	wave_min_K = np.array([2.34474, 2.27431, 2.20846, 2.14638, 2.08678, 2.03186])
	wave_max_K = np.array([2.38127, 2.30955, 2.24218, 2.17881, 2.11799, 2.06194])

	### for K-long
	#wave_min_K = np.array([2.38157,2.31,2.24245,2.17894,2.11878,2.0617])
	#wave_max_K = np.array([2.41566,2.34284,2.27485,2.20861,2.14639,2.08703])

	wave_min_L = np.array([3.40170, 3.25510, 3.11947, 2.99529, 2.88066])
	wave_max_L = np.array([3.45520, 3.30559, 3.16903, 3.04284, 2.92635])

	wave_min_M = np.array([4.95145, 4.64382])
	wave_max_M = np.array([5.03151, 4.71965])

if nirspec == 2:
	import wavecal.WaveCal_orig_Lband_NIRSPEC2 as WaveCal

	n_order_L = 7
	n_order_M = 4
	
	wave_min_L = np.array([3.62763, 3.46095, 3.3101, 3.17376, 3.04799, 2.93166, 2.82366])
	wave_max_L = np.array([3.69645, 3.53032, 3.37729, 3.23805, 3.10872, 2.98951, 2.87938])
	
	wave_min_Mice = np.array([5.3095, 4.94263, 4.64594, 4.36894])
	wave_max_Mice = np.array([5.42629, 5.06949, 4.74741, 4.46483])
	
	wave_min_Mp = np.array([5.40199, 5.04451, 4.73953, 4.45493])
	wave_max_Mp = np.array([5.50795, 5.15746, 4.83772, 4.54388])
	
	wave_min_K = []
	wave_max_K = []
			
importlib.reload(WaveCal)		## python 3.4 and up


if band == 'K' : 
	n_order = n_order_K
	wave_min = wave_min_K
	wave_max = wave_max_K
	#modelfile = '/home/cbuzard/Code/PCA/StellarModels/stellarmodel_HD187123_K_ordersonly.dat'
if band == 'L' : 
	n_order = n_order_L
	wave_min = wave_min_L
	wave_max = wave_max_L
if band == 'M_COice' : 
	n_order = n_order_M
	wave_min = wave_min_Mice
	wave_max = wave_max_Mice
if band == 'M_COp' : 
	n_order = n_order_M
	wave_min = wave_min_Mp
	wave_max = wave_max_Mp


	
#-------------------------------------------------------------------------------------------------------
##### FIT DISPERSION SOLUTIONS #####
# input path to WAVE directory and change # orders
#-------------------------------------------------------------------------------------------------------
ObsFilename=['']*n_order
#wave_path = '2019_04_18/AS205/WAVE/AS205_20190418_115311_wave'
wave_path = '2019_04_18/AS205/WAVE/HR5984_20190418_141954_wave'

## ORDER 1
ObsFilename[0] = wave_path+'1.fits'		# M_ice

## ORDER 2
ObsFilename[1] = wave_path+'2.fits'		# M_ice

## ORDER 3
ObsFilename[2] = wave_path+'3.fits'	        # M_ice

## ORDER 4
ObsFilename[3] = wave_path+'4.fits'	        # M_ice

### ORDER 5
#ObsFilename[4] = 'reducs/2019_04_08/HD187123/WAVE/HD187123_20190408_1429_wave5.fits' 

### ORDER 6
#ObsFilename[5] = 'reducs/2019_04_08/HD187123/WAVE/HD187123_20190408_1429_wave6.fits'

### ORDER 7
#ObsFilename[6] = 'reducs/2019_04_08/HD187123/WAVE/HD187123_20190408_1429_wave7.fits'


if band=='L':
	WaterAbundance = 0.2 #0.5
	CarbonMonoxideAbundance = 1 
	MethaneAbundance = 1
	N2OAbund = 0
	O3Abund = 0
	CO2Abund = 1.2
	
if band=='K':
	# need vbary and vrad for K band wavelength solutions
	# from PyAstronomy import pyasl
	# RA = 360*13/24.+47/60.+16.04/3600 # tau boo
	# Dec = 17.+27/60.+24.39/3600		# tau boo
	# heli, bary = pyasl.baryvel(JD, deq=2000.0)
	# vh, vb = pyasl.baryCorr(JD, RA, Dec, deq=2000.0)

	vrad = -16.965
	#vbary = 15.8769330581	 # 4/6/2015 13:22
	vbary = 17.0942194754	# 5/8/2015 10:30

	#2015apr06
	WaterAbundance = 0.48455072339
	MethaneAbundance = 1.07139188017
	CarbonDioxideAbundance = 1.42491674
	CarbonMonoxideAbundance = 0.835038289258 


if band == 'M_COp' or band == 'M_COice':
	#apr08 abundances # tboo  2015
	WaterAbundance = 0.1
	CarbonDioxideAbundance = 1.2
	CarbonMonoxideAbundance = 1.2 
	#N2OAbundance = 0.2
	#OCSAbundance = 	0.2
	#O3Abundance = 0.2

#-------------------------------------------------------------------------------------------------------
### FitMethod ###
#-------------------------------------------------------------------------------------------------------
#FittingMethod = 'Polynomial'		### initial guess set by first 2 coefs
FittingMethod = 'SetPoints'			## initial guess set by poly fit to (up to 5) set points

#-------------------------------------------------------------------------------------------------------
### SetPointsName ###
#-------------------------------------------------------------------------------------------------------
date 	= '20190418'
target 	= 'AS205'

SetPointsName = []
for order,filename in enumerate(ObsFilename):
        filename = date+'_'+str(order+1)+'_'+target+'.dat' 
        filepath =  'WaveCalSetPoints/'+filename
        with open(filepath) as f: # create SetPointsFile
                if not 'Pos 0 Neg 0' in f.read():
                        f.write("Pos 0 Neg 0")
        SetPointsName.append(filepath)
print(SetPointsName)

# SetPointsName = ['WaveCalSetPoints/20190418_1_AS205.dat',
# 					'WaveCalSetPoints/20190418_2_AS205.dat',
# 					'WaveCalSetPoints/20190418_3_AS205.dat',
# 					'WaveCalSetPoints/20190418_4_AS205.dat']

#-------------------------------------------------------------------------------------------------------
### Wavelength Calibration ###
# Choose order(s) to analyze
#-------------------------------------------------------------------------------------------------------

for order,filename in enumerate(ObsFilename) :
        #if order == 1 : continue		## skip this
	#if order == 0 : continue
#	if order != 0: continue 		## do only this
	print ('-----------------------------------------------------')
	print ('-----------------------------------------------------')
	print ('FILE', ObsFilename[order] )
	print ('ORDER', order+1, 'ORDER', order+1, 'ORDER', order+1, 'ORDER', order+1, 'ORDER', order+1, 'ORDER', order+1 )
	print ('-----------------------------------------------------')
	
	######## default is H2O=0.5, CO=1, have to define all others
	
	if band == 'K' : 
		cal = WaveCal.WaveCal_st(ObsFilename=filename, ModelFilename=modelfile,wave_min=wave_min[order], wave_max=wave_max[order], vrad=vrad, vbary=vbary,
							 H2O_scale=WaterAbundance, CH4_scale=MethaneAbundance, CO2_scale=CarbonDioxideAbundance, 
							 CO_scale = CarbonMonoxideAbundance)
	if band == 'L' : 
		cal = WaveCal.WaveCal(ObsFilename=filename, wave_min=wave_min[order], wave_max=wave_max[order],FittingMethod=FittingMethod,
									SetPointsFile=SetPointsName[order],H2O_scale=WaterAbundance, CH4_scale=MethaneAbundance, CO2_scale=CO2Abund,
									CO_scale = CarbonMonoxideAbundance, N2O_scale=N2OAbund,O3_scale=O3Abund)

	if band == 'M_COice' or band == 'M_COp': 
		cal = WaveCal.WaveCal(ObsFilename=filename, wave_min=wave_min[order], wave_max=wave_max[order],FittingMethod=FittingMethod,SetPointsFile=SetPointsName[order],
							H2O_scale=WaterAbundance, CO2_scale=CarbonDioxideAbundance,CO_scale = CarbonMonoxideAbundance) 
							#N2O_scale=N2OAbundance,OCS_scale=OCSAbundance,O3_scale=O3Abundance)	 

print (cal)

##### QUICK REDUCTION TO GET 1D WAVE SPECTRA FOR FITTING #####

#path = '/home/gablakers/NirspecData/2015nov21/spec/'
#Flat_KL = [(239,366)] 
#FlatDark_KL = [(129,133)]	
#ObsDark = [(230,232)] #Darks to not match observations
#
### 51 Peg
##SciRanges = [(354,357)] #397)]
##StdRanges = [(398,401)]
#
### HD 88133
#SciRanges = [(142,145)]
#StdRanges = [(138,141)]
#
#qr = WaveCal.QuickRed(flat_range=Flat_KL, flat_dark_range=FlatDark_KL, dark_range=ObsDark,
#					  sci_range=[SciRanges[0]], std_range=[StdRanges[0]], path=path, base='nov21',
#					  shift=0.0, dtau=-0.0,level1=True,level2=True, SettingsFile='nirspec_21nov2015.ini',
#					  sci_tname='HD88133', std_tname='HR4024', hold_plots=True)
#
