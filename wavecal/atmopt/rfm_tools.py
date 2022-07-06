import re
import json
import warnings
import os

import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt
import subprocess as sp
import scipy.integrate as integrate

#import pyfits as pf
import astropy.io.fits as pf
from lmfit import minimize, Parameters, Parameter

import tempfile
import os

class Atmosphere():
	'''
	Module for reading, writing and manipulating model atmospheres in the RFM format.

	Arguments:
		filename = string: name of an RFM-compliant atmosphere file. Currently does not accept 
		comma-separated data fields, though. 
	'''
	def __init__(self, filename):

		
		self.profiles = {}		  
		self.read(filename)

		#initial scale factors are 1
		self.scale_facs = {}
		for key in self.profiles.keys():
			self.scale_facs[key] = 1.

		self.profiles_scl = self.profiles.copy()



	def read(self, filename):
		hcond = re.compile('\*')
		ccond = re.compile('\!')

		file  = open(filename, 'r')
		lines = file.readlines()

		for line in lines:
			#Is it a comment line?
			if not ccond.search(line):
				#Is it a header line?
				if hcond.match(line):
					current_header = line.split()[0][1:]
					self.profiles[current_header] = np.array([])
					
				else: 
					self.profiles[current_header] = \
						np.append(self.profiles[current_header],
								  [float(s) for s in line.split()] )

		file.close()
		
		return self.profiles

	def writeAll(self, filename):
		
		file = open(filename,'w')
		file.write('! This atmosphere file is created by the rfm tools module.\n')

		length = self.profiles['HGT'].shape
		file.write('		'+str(length[0]))
		file.write(' ! Profile Levels\n')

		self.writeProf(file,'HGT')
		self.writeProf(file,'PRE')
		self.writeProf(file,'TEM')
		
		for type in self.profiles.keys():
			if type!='HGT' and type!='PRE' and type!='TEM' and type!='END':
				self.writeProf(file,type)

		self.writeProf(file,'END')
		file.close()

	def getPWV(self,alt=0.):
		RR		= 8.314472e7
		h2o_amu = 18.01528
		
		hgt = self.profiles_scl['HGT']
		pre = self.profiles_scl['PRE']
		tem = self.profiles_scl['TEM']
		h2o = self.profiles_scl['H2O']

		hgt_int = np.arange(alt,hgt.max(),0.1)
		pre_int = np.interp(hgt_int,hgt,pre)
		tem_int = np.interp(hgt_int,hgt,tem)
		h2o_int = np.interp(hgt_int,hgt,h2o)

		n_h2o	  = pre_int*1e3*h2o_int*1e-6/(RR*tem_int)*h2o_amu #in g/cm^3
		pwv_g_cm2 = integrate.simps(n_h2o,hgt_int*1e5)
		pwv_mm	  = pwv_g_cm2*10.

		return pwv_mm

	def writeProf(self, file, type):
			if type=='HGT':
				unit = '[km]'
			elif type=='PRE':
				unit = '[mb]'
			elif type=='TEM':
				unit = '[K]'
			else:
				unit = '[ppmv]'
			file.write('*'+type+' '+unit+'\n')

			self.profiles_scl[type].tofile(file, sep="	 ")

			file.write('\n')
		
	def scaleProf(self,type,factor):
		self.scale_facs[type]	= factor
		self.profiles_scl[type] = self.profiles[type]*self.scale_facs[type]

	def subGDAS(self,GDAS_file):
		'''
		Substitute the lower 26 km with a GDAS sounding.
		'''
		profiles = self.profiles
		
		gdas = self._readGDAS(GDAS_file)
		gdas_maxhgt = gdas['HGT'].max()
		gsubs = np.where(profiles['HGT']<gdas_maxhgt)
		ghgts = profiles['HGT'][gsubs]
		keys = ['PRE','TEM', 'H2O']

		for key in keys:
			profiles[key][gsubs] = np.interp(ghgts,gdas['HGT'],gdas[key])
		

	def _readGDAS(self,GDAS_file):
		file = open(GDAS_file, 'r')
		lines = file.readlines()
		data_lines = lines[39:58]
		pre = np.array([])
		hgt = np.array([]) 
		tem = np.array([]) 
		dew = np.array([])
		hum = np.array([])

		for line in data_lines:
			data = [float(s) for s in line.split()]
			pre = np.append(pre,data[0])
			hgt = np.append(hgt,data[1]/1e3)
			tem = np.append(tem,data[2]+273.15)
			dew = np.append(dew,data[3])
			
		data_lines = lines[94:113]
		for line in data_lines:
			data = [float(s) for s in line.split()]
			hum = np.append(hum,data[1])

		pmmv = [self.Rel2AbsHum(hum[i],tem[i],pre[i]) for i in np.arange(hum.size)]
		
		return {'HGT':hgt,'PRE':pre,'TEM':tem,'DEW':dew,'H2O':pmmv,'HUM':hum}
		
	def Rel2AbsHum(self,RH,T,P):
		'''
		arguments:
		RH: Relative humidity (unitless)
		T: Temperature in Kelvin
		P: Pressure in millibar == hPa
		'''
				
		#http://www.vaisala.com/Vaisala%20Documents/Application%20notes/Humidity_Conversion_Formulas_B210973EN-D.pdf
		if T>=273.15:
			Tc = 647.096 #K
			Pc = 220640. #hPa
			C1 = -7.85951783
			C2 = 1.84408259
			C3 = -11.7866497
			C4 = 22.6807411
			C5 = -15.9618719
			C6 = 1.80122502

			V  = 1-T/Tc
			Pw = Pc*np.exp(Tc/T*(C1*V+C2*V**1.5+C3*V**3+C4*V**3.5+C5*V**4+C6*V**7.5))
		elif T<273.15:
			Pn = 6.11657 #hPa
			Tn = 273.16	 #K
			a0 = -13.928169
			a1 = 34.707823
			
			S  = T/Tn
			Pw = np.exp(a0*(1.-S**-1.5)+a1*(1-S**-1.25))

		pmmv = Pw/(P-Pw) * 1e6

		return pmmv

class Model():
	def __init__(self,**kwargs):
		self.Envi = Environment()
		self.rfm_exe = self.Envi.exe_file
		self.outname = 'rfm.tra'

	def RFM(self,**kwargs):
		self.makeDriverFile('rfm.drv',outname=self.outname,**kwargs)
		output = sp.check_output([self.rfm_exe])
		self.spectrum = self.readFromFile(self.outname)
		self.radiance = self.readFromFile('rfm.rad')
		return 1
		
	def makeDriverFile(self,filename,pwv=2.,am=1.05,wmin=2200.,wmax=2210.,
					   sampling=0.01,alt=4.15,
					   mols=['H2O','CO2','O3','CO','CH4'],
					   outname='rfm.tra',atm_files=None):

		if atm_files is None:
			atm_files = [self.Envi.atm_file]
		hitran_path = self.Envi.hit_file
		
		file = open(filename, 'w')

		file.write('*HDR\n		 ')
		file.write('RFM Driver file generated by RFM tools.\n')
		file.write('*FLG\n		 ')
		file.write('OBS TRA RAD ZEN DBL		 !\n')
		file.write('*SPC\n		 ')
		file.write(str(wmin)+' '+str(wmax)+' '+str(sampling)+'\n')
		file.write('*GAS\n		 ')
		for mol in mols:
			file.write(mol+' ')
		file.write('\n')
		file.write('*ATM\n		')
		for atm_file in atm_files:
			file.write(atm_file+'\n')
		file.write('*TAN\n		')
		file.write(str(am)+'\n')
		file.write('*HIT\n		')
		file.write(hitran_path+'\n')
		file.write('*OBS\n		')
		file.write(str(alt)+'\n')
		file.write('*TRA\n		')
		file.write(self.outname+'\n')
		file.write('*RAD\n		')
		file.write('rfm.rad\n')
		file.write('*END\n')

		file.close()

	def readFromFile(self,filename):
		
		ccond = re.compile('\!')
		
		file = open(filename,'r')
		lines = file.readlines()

		first_line_read = False

		data = np.array([])
		for line in lines:
			#Is it a comment line?
			if not ccond.search(line):
				#Is it the first line?
				if not first_line_read:
					first_line = line.split()
					npoints	 = int(first_line[0])
					wmin	 = float(first_line[1])
					sampling = float(first_line[2])
					wmax	 = float(first_line[3])
					label	 = first_line[4]
					label_len = len(label)
					#Remove inverted commas
					label = label[1:label_len-1]

					first_line_read = True
				else:
					data = np.append(data,[float(s) for s in line.split()])
 
		file.close()

		wave = np.linspace(wmin,wmax,npoints)

		#Everything in order of increasing wavelength
		return {'wavelength':1e4/wave[::-1],'wavenumber':wave[::-1],label:data[::-1]}
		
	def blur(self,rpower=1e5):
		wave = self.spectrum['wavelength']
		tran = self.spectrum['Transmission']
		radi = self.radiance['Radiance']
		
		wmin = wave.min()
		wmax = wave.max()
		
		nx = wave.size
		x  = np.arange(nx)
		
		A = wmin
		B = np.log(wmax/wmin)/nx
		wave_constfwhm = A*np.exp(B*x)

		tran_constfwhm = np.interp(wave_constfwhm, wave, tran)		  
		radi_constfwhm = np.interp(wave_constfwhm, wave, radi)	
		
		#Sanity check that the FWHM is indeed constant
		dwdx_constfwhm = np.diff(wave_constfwhm)
		fwhm_pix = wave_constfwhm[1:]/rpower/dwdx_constfwhm
		if (fwhm_pix.min()-fwhm_pix.max())/fwhm_pix.max() > 1e-5: 
			#This should not happen
			warnings.warn('The FWHM is not constant in units of pixels.')
		
		fwhm_pix  = fwhm_pix[0]
		sigma_pix = fwhm_pix/2.3548
		kx = np.arange(nx)-(nx-1)/2.
		kernel = 1./(sigma_pix*np.sqrt(2.*np.pi))*np.exp(-kx**2/(2.*sigma_pix**2)) 
		
		tran_conv = fft.ifft(fft.fft(tran_constfwhm)*np.conj(fft.fft(kernel)))
		tran_conv = fft.fftshift(tran_conv).real
		tran_oldsampling = np.interp(wave,wave_constfwhm,tran_conv)

		radi_conv = fft.ifft(fft.fft(radi_constfwhm)*np.conj(fft.fft(kernel)))
		radi_conv = fft.fftshift(radi_conv).real
		radi_oldsampling = np.interp(wave,wave_constfwhm,radi_conv)

		self.obs_spectrum = {'wavelength':wave,'wavenumber':1e4/wave, 
							 'Transmission':tran_oldsampling,'Radiance':radi_oldsampling}

		return self.obs_spectrum

	def plot(self,unit='cm-1'):
		if unit=='micron':
			plt.plot(self.obs_spectrum['wavelength'],self.obs_spectrum['Transmission'])
		elif unit=='cm-1':
			plt.plot(self.obs_spectrum['wavenumber'],self.obs_spectrum['Transmission'])
		else:
			raise Exception("Unknown unit in model plot")

class ObsSpectrum():
	def __init__(self, filename,input_unit='nm',**kwargs):
		try:
			self._readCRIRES(filename,**kwargs)
		except:
			self._readNIRSPEC(filename,**kwargs)

	def _readNIRSPEC(self,filename,input_unit='micron',beam='pos',**kwargs):
		data = pf.getdata(filename)
		self.header = pf.getheader(filename)
		self.beam = beam
		if beam=='pos':
			self.wave = data['wave_pos']
			self.flux = data['flux_pos']
			self.error = data['uflux_pos']
			self.sky = data['sky_pos']
			self.usky = data['usky_pos']
		if beam=='neg':
			self.wave = data['wave_neg']
			self.flux = data['flux_neg']
			self.error = data['uflux_neg']
			self.sky = data['sky_neg']
			self.usky = data['usky_neg']
			
		self.wave_unit=input_unit
	
	def _readCRIRES(self,filename,input_unit='nm',ext=1,**kwargs):
		data = pf.getdata(filename,ext)
		self.wave  = np.float64(data['wavelength'].squeeze())
		self.flux  = np.float64(data['flux_opt'].squeeze())
		self.error = np.float64(data['uncertainty_opt'].squeeze())
		
		if input_unit=='nm':
			self.wave = self.wave/1e3
			self.wave_unit = 'micron'

	def normalize(self):
		self.flux = self.flux/self.flux.max()

	def getMinMaxWave(self):
		return (self.wave.min(),self.wave.max())

	def plot(self,scale=1,range=None):
		plt.plot(self.wave,self.flux*scale)

class Optimize():
	def __init__(self,spec_file, H2O_scale=0.5, CO_scale=1.0, am=1.,rpower=1e5,cull=50,**kwargs):

		self.cull = cull
		
		self.am = am
		self.Envi = Environment()

		Obs = ObsSpectrum(spec_file,**kwargs)
		Obs.normalize()
		self.wave_min_max = Obs.getMinMaxWave()

		A = Atmosphere(filename=self.Envi.atm_file)
		if 'gdas' in kwargs:
			A.subGDAS(kwargs['gdas'])
		A.scaleProf('H2O',H2O_scale)
		A.scaleProf('CO',CO_scale)
		tf=tempfile.NamedTemporaryFile()
		
		A.writeAll(tf.name+'.atm')
		
		M = Model()
		M.RFM(atm_files=[tf.name+'.atm'],wmin=1e4/self.wave_min_max[1],wmax=1e4/self.wave_min_max[0],
			  am=am)
		os.remove(tf.name+'.atm')
		M.blur()

		self.A	 = A
		self.Obs = Obs
		self.M	 = M

		self.bestfit = self.wavecal(rpower=rpower)

		self.Obs.wave = self.bestfit['wave']
		
	def wavecal(self,rpower=1e5):

		size = self.Obs.wave.size
		try:
			Aini = self.Obs.header['WAVA_'+self.Obs.beam]
			Bini = self.Obs.header['WAVB_'+self.Obs.beam]
			Cini = self.Obs.header['WAVC_'+self.Obs.beam]
			rpower = self.Obs.header['RPOW_'+self.Obs.beam]
		except:
			Aini = self.Obs.wave[0]
			Bini = np.abs(self.Obs.wave[size-1]-self.Obs.wave[0])/size
			Cini = 0.

		params = Parameters()
		params.add('contA', value=1.0)
		params.add('contB', value=0.,vary=False)
		params.add('R',value=rpower, min=rpower*0.8,max=rpower*1.2,vary=True)
		params.add('A',value=Aini,vary=True) 
		params.add('B',value=Bini,vary=True) 
		params.add('C',value=Cini*1e6,vary=True)
		params.add('D',value=1e-6,vary=True)
		params.add('E',value=1e-6,vary=True) 

		## NRC EDITED
		#params = Parameters()
		#params.add('contA', value=1.0)
		#params.add('contB', value=0., vary=False)
		#params.add('R',value=rpower, vary=True, min=rpower*0.8, max=rpower*1.2)
		#params.add('A',value=Aini,vary=True, min=0.9*Aini, max=1.1*Aini) #NRC (05 SEPT 2014)
		#params.add('B',value=Bini,vary=True, min=0.9*Bini, max=1.1*Bini) #NRC (05 SEPT 2015)
		#params.add('C',value=Cini*1e6,vary=True, min=0.9*Cini*1e6, max=1.1*Cini*1e6) #NRC (05 SEPT 2015)
		#params.add('D',value=1.0e-6,vary=True)
		#params.add('E',value=1.0e-6,vary=True) 
		#########
		#self.Obs.sky = self.Obs.sky + 20000.0

		
		leastsq_kws = {'factor':1e-6}
		result = minimize(self._wavecalResidual, params, args=(self.M,self.Obs),method='powell')#,nan_policy='omit')#,**leastsq_kws)
		params = result.params		  
		residual = self._wavecalResidual #result.residual

		#NRC (25 July 2014) - CODE NOW ALSO PRINTs A,B, AND C
		print ('')
		print (self.Obs.beam)
		print ('-----')
		print ('A = ', params['A'])
		print ('B = ', params['B']) 
		print ('C = ', params['C'])
		print ('D = ', params['D'])
		print ('E = ', params['E'])
		print ('R = ', params['R'])
		#print 'RedChi / Success = ', result.redchi, ' / ', result.success
		print ('')
		print ('=====')

		index	 = np.arange(size)
		wave_fit = params['A'].value+params['B'].value*index+params['C'].value*(index/1000.)**2.+ \
			params['D'].value*(index/1000.)**3.+params['E'].value*(index/1000.)**4.
		transmission = np.interp(wave_fit,self.M.obs_spectrum['wavelength'],self.M.obs_spectrum['Transmission'])
		radiance = np.interp(wave_fit,self.M.obs_spectrum['wavelength'],self.M.obs_spectrum['Radiance'])*params['contA'].value
		
		return {'wave':wave_fit,'Transmission':transmission,
				'Radiance':radiance,
				'flux':self.Obs.flux,'sky':self.Obs.sky,
				'residual':residual, 'R_fit':params['R'], 'A_fit':params['A'],'B_fit':params['B'],
				'D_fit':params['D'],'E_fit':params['E'],'C_fit':params['C'],'contA':params['contA'],'contB':params['contB']}

	def _wavecalResidual(self, params, Model, Obs):
		contA  = params['contA'].value
		contB  = params['contB'].value
		Aval   = params['A'].value
		Bval   = params['B'].value
		Cval   = params['C'].value
		Dval   = params['D'].value
		Eval   = params['E'].value
		Rval   = params['R'].value

		index  = np.arange(Obs.wave.size)
		wave   = Aval+Bval*index+Cval*(index/1000.)**2.+Dval*(index/1000.)**3.+Eval*(index/1000.)**4.
		
		cont   = contA#+contB*index

		Model.blur(rpower=Rval)

		model_int = np.interp(wave,Model.obs_spectrum['wavelength'],Model.obs_spectrum['Radiance'])
		model_scl = model_int*cont
#		 import pdb;pdb.set_trace()
		gsubs = np.where(Obs.sky/Obs.usky>20.)
		residual = (Obs.sky-model_scl)[gsubs]#[self.cull:Obs.sky.size-self.cull]
		#/Obs.error[self.cull:Obs.flux.size-self.cull]
				
		return 	residual

	def scaleSpecies(self, Rval,H2O_scale):
		params = Parameters()
		params.add('H2O_scale', value=H2O_scale,min=0.01,max=3.,vary=True)
		params.add('CO2_scale', value=1.,vary=False)
		params.add('CO_scale', value=1.,vary=False)
		params.add('O3_scale', value=1.,vary=False)
		params.add('CH4_scale', value=1.,vary=False)
		params.add('AM', value = self.am, vary=False)
		params.add('scale', value=1., vary=True)
		params.add('R', value=1e5, vary=True)
		
		result = minimize(self._scaleSpeciesResidual, params, 
						  args=(self.M, self.Obs, self.A,self.wave_min_max))
		params = result.params
	   
		return {'wave':self.Obs.wave,'Transmission':self.Obs.flux-result.residual}

	def _scaleSpeciesResidual(self, params, Model, Obs, Atm, wave_min_max):
		
		H2O_scale = params['H2O_scale'].value
		CO2_scale = params['CO2_scale'].value
		CO_scale  = params['CO_scale'].value
		O3_scale  = params['O3_scale'].value
		CH4_scale = params['CH4_scale'].value
		Rval	  = params['R'].value
		AMval	  = params['AM'].value
		scale = params['scale'].value

		Atm.scaleProf('H2O',H2O_scale)
		Atm.scaleProf('CO2',CO2_scale)
		Atm.scaleProf('CO',CO_scale)
		Atm.scaleProf('O3',O3_scale)
		Atm.scaleProf('CH4',CH4_scale)
		Atm.writeAll('tmp.atm')

		Model.RFM(atm_files=['tmp.atm'],wmin=1e4/wave_min_max[1],wmax=1e4/wave_min_max[0],
				  am=AMval)
		Model.blur(rpower=Rval)

		model_int = np.interp(Obs.wave,Model.obs_spectrum['wavelength'],Model.obs_spectrum['Transmission'])
		
		gsubs = np.where(Obs.flux>0.5)
		return (model_int*scale-Obs.flux)

	def getWave(self):
		return self.Obs.wave

	def getModel(self):
		return {'wave':self.bestfit['wave'],'Transmission':self.bestfit['Transmission'],
				'flux':self.bestfit['flux'],'Radiance':self.bestfit['Radiance'],'sky':self.bestfit['sky']}
	
	def plotFit(self):
		plt.plot(self.Obs.wave,self.Obs.flux,drawstyle='steps-mid')
		plt.plot(self.M.spectrum['wavelength'],self.M.obs_spectrum['Transmission'])
		plt.show()
	   
class Environment():
	def __init__(self,exe_file='RFM/bin/rfm', atm_file='RFM/ATM/equ.atm', hit_file='RFM/HIT/hitran08.hit'):
		sys_dir = self.getSysPath()
		
		self.exe_file = os.path.join(sys_dir,exe_file)
		self.atm_file = os.path.join(sys_dir,atm_file)
		self.hit_file = os.path.join(sys_dir,hit_file)

	def getSysPath(self):
		sys_dir, this_filename = os.path.split(__file__)
		return sys_dir

		
