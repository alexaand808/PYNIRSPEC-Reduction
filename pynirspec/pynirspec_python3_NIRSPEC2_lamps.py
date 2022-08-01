import warnings
import json
import os
#import ConfigParser as cp
import configparser as cp
import numpy as np
import numpy.ma as ma
import astropy.io.fits as pf
import scipy.fftpack as fp
from scipy.stats import tmean, tvar, sigmaclip
from scipy.ndimage.filters import median_filter, maximum_filter1d, maximum_filter, gaussian_filter1d
from scipy.optimize import curve_fit
from scipy import constants
import scipy.signal
from scipy.signal import medfilt
import matplotlib.pylab as plt
from collections import Counter
import sys
import inpaint as inpaint
sys.path.insert(0, '../')
import wavecal.atmopt.rfm_tools as rfm
from PyAstronomy import pyasl
from lmfit import minimize, Parameters, Parameter
from itertools import groupby, count
from operator import itemgetter
from scipy.interpolate import UnivariateSpline


class Environment():
	'''
	Class to encapsulate global environment parameters.

	Parameters
	----------
	config_file: string
		ConfigParser compliant file containing global parameters.

	Methods
	-------
	
	Attributes
	----------
	pars: SafeConfigParser object
		Contains all global parameters

	'''
	def __init__(self,settings_file=None,detpars_file=None,nirspec=1):

		sys_dir = self._getSysPath()
		#print ("sys_dir: " + sys_dir)

		#self.settings = cp.SafeConfigParser()
		self.settings = cp.ConfigParser()
		if settings_file is None:
			#use the system settings
			self.settings.read(sys_dir+'/'+'nirspec.ini')
		else:
			self.settings.read(settings_file)
			
		#self.detpars = cp.SafeConfigParser()
		self.detpars = cp.ConfigParser()
		if nirspec == 1:
			if detpars_file is None:
				self.detpars.read(sys_dir+'/'+'detector_nirspec1.ini')
			else:
				self.detpars.read(detpars_file)
		if nirspec == 2:
			if detpars_file is None:
				self.detpars.read(sys_dir+'/'+'detector_nirspec2.ini')
			else:
				self.detpars.read(detpars_file)



	def _getSysPath(self):
		sys_dir, this_filename = os.path.split(__file__)
		return sys_dir

	def getItems(self,option):
		return [self.settings.get(section,option) for section in self.settings.sections()]

	def getSections(self):
		return self.settings.sections()

	def getWaveRange(self,setting,onum):
		range_str = self.settings.get(setting,'wrange'+str(int(onum)))
		range = json.loads(range_str)
		return range
		
	def getYRange(self,setting,onum):
		range_str = self.settings.get(setting,'yrange'+str(int(onum)))
		range = json.loads(range_str)
		return range

	def getDispersion(self,setting,onum):
		A_str = self.settings.get(setting,'A'+str(int(onum)))
		B_str = self.settings.get(setting,'B'+str(int(onum)))
		C_str = self.settings.get(setting,'C'+str(int(onum)))
		R_str = self.settings.get(setting,'R'+str(int(onum)))
		As = json.loads(A_str)
		Bs = json.loads(B_str)
		Cs = json.loads(C_str)
		Rs = json.loads(R_str)
		return {'A':As,'B':Bs,'C':Cs,'R':Rs}

	# Returns a dict containing detector gain, read noise, dark current
	def getDetPars(self):
		gain = self.detpars.getfloat('Detector','gain')
		rn	 = self.detpars.getfloat('Detector','rn')
		dc	 = self.detpars.getfloat('Detector','dc')
		return {'gain':gain,'rn':rn,'dc':dc}

	# How many orders of spectra
	def getNOrders(self,setting):
		return self.settings.getint(setting,'norders')

class Observation():
	'''
	Private object containing a NIRSPEC observation - that is, all exposures related to a single type of activity.

	Any specific activity (Darks, Flats, Science, etc.) are modeled as classes derived off the Observation class.

	Parameters
	----------
	filelist: List of strings
		List of data (.fits) files associated with the observation.
	type: string
		type of observation (e.g., Dark). 

	Attributes
	----------
	type	# darks/flats/sci etc.
	Envi	# instrument settings and observation details (as contained in the .ini files)
	flist	# file list
	planes	# list of opened fits files: e.g. planes[0] is the first file, and planes[0][0] is the PrimaryHDU of the fits file
	header	# header of a specific fits file

	Methods
	-------
	getSetting
	getNOrders
	getTargetName
	subtractFromStack
	divideInStack
	writeImage
   
	'''
	def __init__(self,filelist,type='image', SettingsFile=None,
				 tname=None,out_dir=None,nirspec=2):
		self.type = type
		self.Envi = Environment(settings_file=SettingsFile,nirspec=nirspec)
		self.flist = filelist
		self.nirspec = nirspec
		self._openList(tname)
		self._makeHeader(tname)
		self.sci_tname = tname
		self.output_path = out_dir


	# Open the files and get exp and det parameters from header and .ini file respectively
	def _openList(self, tname):
		warnings.resetwarnings()
		warnings.filterwarnings('ignore', category=UserWarning, append=True)

		# Construct a list of opened fits files (HDU lists)
		self.planes = []
		for file in self.flist:
			plane = pf.open(file,ignore_missing_end=True)
			self.planes.append(plane)
		self._makeHeader(tname)

		# Get exposure parameters from header
		self.exp_pars = self.getExpPars()
		# Get detector parameters from .ini files
		self.det_pars = self.Envi.getDetPars()

	# Extract header info from the first file
	def _makeHeader(self, tname=None):
		self.header = self.planes[0][0].header 
		
		#The NIRSPEC headers have a couple of illegally formatted keywords
		try:
			del self.header['GAIN.SPE']
		except:
			print ('GAIN.SPE already deleted from hdr.')
		try:
			del self.header['FREQ.SPE']
		except:
			print ('FREQ.SPE already deleted from hdr.')

		# Replace object name in header with custom supplied name
		if (tname is not None):
			self.header['OBJECT'] = tname.replace(' ','')

	# Returns a dict of exposure info
	def getExpPars(self):
		coadds	 = self.getKeyword('COADDS')[0]
		sampmode = self.getKeyword('SAMPMODE')[0]
		if self.nirspec == 1:
			nreads	  = self.getKeyword('MULTISPE')[0]
			itime	 = self.getKeyword('ITIME')[0]
		if self.nirspec == 2:	
			nreads	 = self.getKeyword('NUMREADS')[0]
			itime	 = self.getKeyword('TRUITIME')[0]
		# Number of images
		nexp	 = len(self.planes)
		return {'coadds':coadds,'sampmode':sampmode,'nreads':nreads,'itime':itime,'nexp':nexp}

	# Return correct setting name, echelle pos, and crossdisp pos
	def getSetting(self):
		echelle	  = self.getKeyword('ECHLPOS')
		crossdisp = self.getKeyword('DISPPOS')

		# Check that all exposures are taken with the same echelle position and cross disperser position
		assert len([e for e in echelle if e==echelle[0]])==len(echelle), \
			'All exposures must be taken with the same setting!'
		assert len([c for c in crossdisp if c==crossdisp[0]])==len(crossdisp), \
			'All exposures must be taken with the same setting!'

		echelle	  = echelle[0]
		crossdisp = crossdisp[0]

		echelles   = self.Envi.getItems('echelle')
		crossdisps = self.Envi.getItems('crossdisp')

		# Note: the following steps are unnecessary when user supplies a custom .ini file
		# Get indices of echelle position in .ini file that match that in the header
		setsub1 = [i for i,v in enumerate(echelles) if float(v)==echelle]

		# Get indices of crossdisp position in .ini file that match that in the header
		setsub2 = [i for i,v in enumerate(crossdisps) if float(v)==crossdisp]

		# Determine the correct section (configuration) as
		# the one containing both the correct echelle pos and crossdisp pos
		sub = [i for i in setsub1 if i in setsub2]
		return self.Envi.getSections()[sub[0]],echelle,crossdisp

	# Return a list of airmasses, converting faulty values to the mean
	def getAirmass(self):
		airmasses = self.getKeyword('AIRMASS') 
		airmasses = np.array([airmass if type(airmass) is not str else -1 for airmass in airmasses])
		mean = np.mean(airmasses[np.where(airmasses!=-1)])
		airmasses[np.where(airmasses==-1)] = mean
		return airmasses

	# Get number of orders
	def getNOrders(self):
		setting,echelle,crossdisp = self.getSetting()
		return self.Envi.getNOrders(setting)

	def getTargetName(self):
		target_name = self.header['OBJECT']
		return target_name

	# Prepare data stack in units of e- counts, and an error stack
	def _getStack(self):
		
		nexp = self.exp_pars['nexp']
		nx = self.planes[0][0].header['NAXIS1']
		ny = self.planes[0][0].header['NAXIS2']
		
		stack = np.zeros((nx,ny,nexp))
		ustack = np.zeros((nx,ny,nexp))
		
		for i,plane in enumerate(self.planes):
			# data stack (counts*gain = e- count)
			stacki = plane[0].data*self.det_pars['gain'] #convert everything to e-
			#if self.nirspec == 2:
			#	stacki = np.transpose(stacki)	## change to orientation of NIRSPEC1
			#	stacki = np.flip(stacki,axis=1)
			stack[:,:,i]  = stacki
			# error stack
			ustacki = self._error(plane[0].data)
			#if self.nirspec == 2:
			#	ustacki = np.transpose(ustacki)
			#	ustacki = np.flip(ustacki,axis=1)
			ustack[:,:,i] = ustacki

		return stack,ustack

	# Compute uncertainty as sqrt of variance
	def _error(self,data):
		var_data = np.abs(data+self.exp_pars['itime']*self.det_pars['dc']+
						  self.det_pars['rn']**2/self.exp_pars['nreads'])
		return np.sqrt(var_data)

	# Method to subtract darks from flats
	def subtractFromStack(self,Obs):
		warnings.resetwarnings()
		warnings.filterwarnings('ignore', category=RuntimeWarning, append=True)
		try:
			self.uimage = np.sqrt(self.uimage**2+Obs.uimage**2)
			self.image -= Obs.image
		except:
			print ('Subtraction failed - no image calculated')

	# Flat-field correction
	def divideInStack(self,Obs,lamp=None):
		warnings.resetwarnings()
		warnings.filterwarnings('ignore', category=RuntimeWarning, append=True)
		
		### identify where bad pixels are in flats
		#img_medfilt = medfilt(Obs.image, 3)
		#sigma = np.median(Obs.image)
		#diff = np.abs((Obs.image - img_medfilt))
		#sigma = np.std(diff)
		#plt.figure()
		#plt.imshow(np.abs((Obs.image - img_medfilt)) / sigma, cmap='gray')
		#plt.colorbar()
		#plt.figure()
		#plt.imshow(np.abs((Obs.image - img_medfilt)), cmap='gray')
		#plt.colorbar()
		#badimg = np.abs((Obs.image - img_medfilt)) / sigma > 3.0
		#print(np.mean(diff),np.median(diff),np.std(diff))
		#plt.show()
		#assert(1==0)
		#badimg = np.abs((Obs.image - img_medfilt)) / sigma == 1.0
		#ones = np.ones(np.shape(Obs.image))
		#ones[badimg] = 0
		#plt.imshow(ones, cmap='gray')
		#plt.colorbar
		#plt.show()
		#assert(1==0)
		#Obs.image[badimg] = img_medfilt[badimg]
		
		if not lamp:
			stacktouse, ustacktouse = self.stack, self.ustack
		else:	
			stacktouse, ustacktouse = lamp.image, lamp.uimage
		
		cube = np.zeros((stacktouse.shape[2], stacktouse.shape[0], stacktouse.shape[1]))
		cube2 = np.zeros((stacktouse.shape[2], stacktouse.shape[0], stacktouse.shape[1]))
		for i in np.arange(stacktouse.shape[2]):
			plane = stacktouse[:,:,i]
			### the orders shifted when we switched settings (on 4/18/19) 
			## so check where flats are relative to sci images
			cube[i,:,:] = plane
			
			uplane = ustacktouse[:,:,i]
			uplane = np.sqrt((uplane/plane)**2+(Obs.uimage/Obs.image)**2)
			# Divide data by normalized flats
			plane /= Obs.image
			## replace bad pixels identified in flats
			#img_maxfilt = maximum_filter1d(np.abs(plane), 3, axis=0)
			#plane[badimg] = img_maxfilt[badimg]
			
			uplane *= np.abs(plane)
			
			### the orders shifted when we switched settings (on 4/18/19) 
			## so check where flats are relative to sci images
			cube2[i,:,:] = Obs.image

			stacktouse[:,:,i]  = plane
			ustacktouse[:,:,i] = uplane
		
		#pf.writeto(self.plots_dir + '/sciimagethatwillbedivided.fits', cube, overwrite=True)
		#pf.writeto(self.plots_dir + '/flattodivide.fits', cube2, overwrite=True)
		#assert(1==0)	
		
		if not lamp:
			self.stack, self.ustack = stacktouse, ustacktouse
		else:	
			lamp.image, lamp.uimage = stacktouse, ustacktouse
		
		
	# Collapse darks and flats sequences into single images, with a weighted average (assuming Gaussian error)
	def _collapseStack(self,stack=None,ustack=None,method='SigClip',sig=50.):
		'''
		If called without the stack keyword set, this will collapse the entire stack.
		However, the internal stack is overridden if a different stack is passed.
		For instance, this could be a stack of nod pairs.
		'''
		if stack is None:
			stack,ustack = self.stack,self.ustack

		# Mask invalid points so that errors do not occur in math operations
		masked_stack = ma.masked_invalid(stack)
		masked_ustack = ma.masked_invalid(ustack)

		# Take weighted average of a sequence of images
		image = ma.average(masked_stack,2,weights=1./masked_ustack**2)
		uimage = np.sqrt(ma.mean(masked_ustack**2,2)/ma.count(masked_ustack,2))

		return image, uimage

	# Save image (data + uncertainty)
	def writeImage(self,filename=None):
		if filename is None:
			filename = self.type+'.fits'
			
		hdu = pf.PrimaryHDU(self.image.data)
		uhdu = pf.ImageHDU(self.uimage.data)
		hdulist = pf.HDUList([hdu,uhdu])
		hdulist.writeto(filename,overwrite=True)

	# Returns list of values for a specific header key
	def getKeyword(self,keyword):
		try:
			klist = [plane[0].header[keyword] for plane in self.planes]
			return klist
		except ValueError:
			print ("Invalid header keyword")


## Flat class, extends the class Observation
class Lamp(Observation):
	def __init__(self,filelist,dark=None,norm_thres=5000.,save=False,sci_tname=None,out_dir=None,nirspec=2,**kwargs):

		# Inherits __init__ function from Observation class
		Observation.__init__(self,filelist,tname=sci_tname,out_dir=out_dir,nirspec=nirspec,**kwargs)
		self.type = 'lamp'
		# Convert images to units of e- count
		self.stack,self.ustack = self._getStack()
		# Number of flats
		self.nplanes = self.stack.shape[2]
		# Frame size of flats
		self.height = self.stack[:,:,0].shape[0]
		# Get setting parameters
		self.setting,self.echelle,self.crossdisp = self.getSetting()
		## Combine flats through a weighted average
		#self.image,self.uimage = self._collapseStack()
		self.image,self.uimage = self.stack,self.ustack
		self.out_dir = out_dir
		self.nirspec = nirspec

		# Whether to subtract flat darks from flats
		if dark:
			self.subtractFromStack(dark)
		# Normalize flats by median
		self._normalize(norm_thres)

		# Where the flat field is faulty, it's set to 1 to avoid divide by zeros.
		#self.image[np.where(self.image<0.1)] = 1
		#if save:
		#	self.writeImage(filename=self.out_dir + self.type + '.fits')
		
		## make the right shape (1,2048,2048)
		#self.image = np.array([self.image])
		#self.uimage = np.array([self.uimage])
		
	# Normalize flats with median value (larger than a given threshold)
	def _normalize(self,norm_thres):
		flux = np.median(self.image[np.where(self.image>norm_thres)])
		self.image = self.image/flux
		self.uimage = self.uimage/flux
	# Not used
	def makeMask(self):
		return np.where(self.image<0.1,0,1)


## Flat class, extends the class Observation
class Flat(Observation):
	def __init__(self,filelist,dark=None,norm_thres=5000.,save=False,sci_tname=None,out_dir=None,nirspec=2,**kwargs):

		# Inherits __init__ function from Observation class
		Observation.__init__(self,filelist,tname=sci_tname,out_dir=out_dir,nirspec=nirspec,**kwargs)
		self.type = 'flat'
		# Convert images to units of e- count
		self.stack,self.ustack = self._getStack()
		# Number of flats
		self.nplanes = self.stack.shape[2]
		# Frame size of flats
		self.height = self.stack[:,:,0].shape[0]
		# Get setting parameters
		self.setting,self.echelle,self.crossdisp = self.getSetting()
		# Combine flats through a weighted average
		self.image,self.uimage = self._collapseStack()
		self.out_dir = out_dir
		self.nirspec = nirspec

		# Whether to subtract flat darks from flats
		if dark:
			self.subtractFromStack(dark)
		# Normalize flats by median
		self._normalize(norm_thres)

		# Where the flat field is faulty, it's set to 1 to avoid divide by zeros.
		self.image[np.where(self.image<0.1)] = 1
		if save:
			self.writeImage(filename=self.out_dir + self.type + '.fits')

	# Normalize flats with median value (larger than a given threshold)
	def _normalize(self,norm_thres):
		flux = np.median(self.image[np.where(self.image>norm_thres)])
		self.image = self.image/flux
		self.uimage = self.uimage/flux
	# Not used
	def makeMask(self):
		return np.where(self.image<0.1,0,1)

## Dark class, extends the class Observation
class Dark(Observation):
	def __init__(self,filelist,save=False,sci_tname=None,out_dir=None,nirspec=2,**kwargs):
		# Inherits __init__ function from Observation class
		Observation.__init__(self,filelist,tname=sci_tname,out_dir=out_dir,nirspec=nirspec,**kwargs)
		self.type = 'dark'
		self.out_dir = out_dir
		# Convert images to units of e- count
		self.stack,self.ustack = self._getStack()
		# Number of darks
		self.nplanes = self.stack.shape[2]
		# Frame size of darks
		self.height = self.stack[:,:,0].shape[0]
		# Combine darks using a weighted average
		self.image,self.uimage = self._collapseStack()
		# Make bad pixel map (see function below)
		self._badPixMap(filename=self.out_dir + '/badpix.dmp')
		# Whether to save combined dark frame
		if save:
			self.writeImage(filename=self.out_dir + 'flat-' + self.type + '.fits')
		self.nirspec = nirspec

	# Create bad pixel map
	def _badPixMap(self,clip=30,filename='badpix.dmp'):
		median = np.median(self.image)
		# Variance is computed for values between -100 and 100 - why a range???
		var	 = tvar(self.image,(-100,100))
		# Bad pixels are defined as those where the count-median is greater than 30*Sqrt(variance)
		self.badpix = ma.masked_greater(self.image-median,clip*np.sqrt(var))

		# Save the bad pixel map as a numpy masked array
		if filename is not None:
			self.badpix.dump(filename)

## Class for science images (both A and B nods). Prepares A-B images, as well as a sky image
class Nod(Observation):
	def __init__(self,filelist,dark=None,flat=None,lamp=None,badpix='badpix.dmp',plots_dir=None,tname=None,nirspec=2,**kwargs):

		# Inherits Observation
		Observation.__init__(self,filelist,nirspec=nirspec,**kwargs)
		self.type = 'nod'				 
		self.setting,self.echelle,self.crossdisp = self.getSetting()
		# Get list of airmasses through the science sequence
		self.airmasses = self.getAirmass()
		# Mean airmass
		self.airmass = np.mean(self.airmasses)
		self.plots_dir = plots_dir
		self.tname = tname
		self.nirspec = nirspec
		
		# List of RA, DEC, and file numbers for science images
		RAs	 = self.getKeyword('RA')
		## convert RAs from h:m:s to deg
		RAs_deg = []
		for ra in RAs:
			h,m,s = ra.split(':')
			ra_deg = pyasl.hmsToDeg(float(h),float(m),float(s))
			RAs_deg.append(ra_deg)
		RAs = RAs_deg

		DECs = self.getKeyword('DEC')
		## convert Decs from deg:m:s to deg
		DECs_deg = []
		for dec in DECs:
			d,m,s = dec.split(':')
			dec_deg = pyasl.dmsToDeg(float(d),float(m),float(s))
			DECs_deg.append(dec_deg)
		DECs = DECs_deg	  
		
		if self.nirspec == 1:
			FileNums = self.getKeyword('FILENUM')
		if self.nirspec == 2:
			FileNums = self.getKeyword('FRAMENUM')		  

		#  Get indices of each AB nod pair, e.g. [(0, 1), (3, 2), (4, 5), (7, 6), (8, 9), (11, 10)]
		pairs = self._getPairs(RAs,DECs,FileNums)
		# Get sequence of A-B images from nod pairs
		self.stack,self.ustack = self._makePairStack(pairs)

		# Frame size
		self.height = self.stack[:,:,0].shape[0]

		# Save each A-B pair as fits cube
		cube = np.zeros((self.stack.shape[2], self.stack.shape[0], self.stack.shape[1]))
		for slice in range(self.stack.shape[2]):
			im = self.stack[:, :, slice]
			cube[slice, :, :] = im
		print(np.shape(cube))
		pf.writeto(self.plots_dir+'/A-B-before-flattening.fits', cube, overwrite=True)

		## stop here to see where shift happens
		#assert(1==0)

		# Divide by normalized flats
		if flat:
			self.divideInStack(flat)
			## divide lamps by flat
			self.divideInStack(flat,lamp=lamp)		
		
		#print(lamp_stack)
		#print(np.shape(lamp_stack))

		cube = np.zeros((self.stack.shape[2], self.stack.shape[0], self.stack.shape[1]))
		for slice in range(self.stack.shape[2]):
			im = self.stack[:, :, slice]
			cube[slice, :, :] = im
		pf.writeto(self.plots_dir + '/A-B-before-bpRemovel-after-flattening.fits', cube, overwrite=True)
		
		#cube = np.zeros((self.ustack.shape[2], self.ustack.shape[0], self.ustack.shape[1]))
		#for slice in range(self.ustack.shape[2]):
		#	im = self.ustack[:, :, slice]
		#	cube[slice, :, :] = im
		#pf.writeto(self.plots_dir + '/A-B-before-bpRemovel-after-flattening-unc.fits', cube, overwrite=True)

		cube = np.zeros((lamp.image.shape[2], lamp.image.shape[0], lamp.image.shape[1]))
		for slice in range(lamp.image.shape[2]):
			im = lamp.image[:, :, slice]
			cube[slice, :, :] = im
		pf.writeto(self.plots_dir + '/lamps-after-flattening.fits', cube, overwrite=True)


		# Correct bad pixels through a weighted local mean
		if badpix:
			if dark:
				badmask = np.load(badpix,allow_pickle=True)
				self._correctBadPix(badmask)
				self._correctBadPix(badmask,lamp=lamp)

		cube = np.zeros((self.stack.shape[2], self.stack.shape[0], self.stack.shape[1]))
		for slice in range(self.stack.shape[2]):
			im = self.stack[:, :, slice]
			cube[slice, :, :] = im
		pf.writeto(self.plots_dir + '/A-B-after-bpRemoval.fits', cube, overwrite=True)
		
		cube = np.zeros((lamp.image.shape[2], lamp.image.shape[0], lamp.image.shape[1]))
		for slice in range(lamp.image.shape[2]):
			im = lamp.image[:, :, slice]
			cube[slice, :, :] = im
		pf.writeto(self.plots_dir + '/lamps-after-bpRemoval.fits', cube, overwrite=True)
		
		
		# Reassign variable names for A-B images
		self.TargetStack, self.UTargetStack = self.stack, self.ustack
		self.LampStack, self.ULampStack = lamp.image, lamp.uimage

		# An averaged sky frame is constructed using the off-beam pixels
		# Separated stacks for A images and B images
		stackA,ustackA,stackB,ustackB = self._makeSingleBeamStacks(pairs)
		beamStacks = [(stackA,ustackA),(stackB,ustackB)] 

		# Creates preprocessed images for A image and B image separately)
		beam_sky_stacks, beam_usky_stacks = [], []
		count = 0
		for beamStack in beamStacks:
			self.stack = beamStack[0]
			self.ustack = beamStack[1]

			if count ==0:
				fname1 = self.plots_dir + '/Anod-after-darksub.fits'
				fname2 = self.plots_dir + '/Anod-after-flattening.fits'
				fname3 = self.plots_dir + '/Anod-after-bpRemoval.fits'
			elif count ==1:
				fname1 = self.plots_dir + '/Bnod-after-darksub.fits'
				fname2 = self.plots_dir + '/Bnod-after-flattening.fits'
				fname3 = self.plots_dir + '/Bnod-after-bpRemoval.fits'

			cube = np.zeros((self.stack.shape[2], self.stack.shape[0], self.stack.shape[1]))
			for slice in range(self.stack.shape[2]):
				im = self.stack[:, :, slice]
				uim = self.ustack[:, :, slice]
				self.image = im
				self.uimage = uim
				### Subtract dark to test destriping
				if dark:
					self.subtractFromStack(dark)
				# Update the stack after dark subtracting for one image
				self.stack[:,:,slice] = self.image
				self.ustack[:,:,slice] = self.uimage
				cube[slice, :, :] = self.image
			print(np.shape(cube))
			pf.writeto(fname1, cube, overwrite=True)
			#assert(1==0)

			if flat:
				self.divideInStack(flat)
			cube = np.zeros((self.stack.shape[2], self.stack.shape[0], self.stack.shape[1]))
			for slice in range(self.stack.shape[2]):
				im = self.stack[:, :, slice]
				cube[slice, :, :] = im
			pf.writeto(fname2, cube, overwrite=True)

			if badpix:
				if dark:
					badmask = np.load(badpix,allow_pickle=True)
					self._correctBadPix(badmask)

			cube = np.zeros((self.stack.shape[2], self.stack.shape[0], self.stack.shape[1]))
			for slice in range(self.stack.shape[2]):
				im = self.stack[:, :, slice]
				cube[slice, :, :] = im
			pf.writeto(fname3, cube, overwrite=True)

			beam_sky_stacks.append(self.stack)
			beam_usky_stacks.append(self.ustack)
			count =+ 1

		self.beamSkyStacks	= beam_sky_stacks
		self.beamUSkyStacks = beam_usky_stacks
		#assert(1==0)
		
	# Correct bad pixels using an iterative local mean method from inpaint.py
	def _correctBadPix(self,badmask,lamp=None):
		
		if not lamp:
			stacktouse, ustacktouse = self.stack, self.ustack
		else:
			stacktouse, ustacktouse = lamp.image, lamp.uimage
	
		for i in np.arange(stacktouse.shape[2]):
			plane = stacktouse[:,:,i]			 
			maskedImage = np.ma.array(plane, mask=badmask.mask)
			NANMask = maskedImage.filled(np.NaN)
			stacktouse[:,:,i] = inpaint.replace_nans(NANMask, method='localmean')		   
			
			plane = ustacktouse[:,:,i]			  
			maskedImage = np.ma.array(plane, mask=badmask.mask)
			NANMask = maskedImage.filled(np.NaN)			
			ustacktouse[:,:,i] = inpaint.replace_nans(NANMask,method='localmean')		   
		
		if not lamp:
			self.stack, self.ustack = stacktouse, ustacktouse
		else:
			lamp.image, lamp.uimage = stacktouse, ustacktouse 
			
	# Compute A-B image for each AB nod pair
	def _makePairStack(self,pairs):
		# Number of AB pairs
		npairs = len(pairs)
		# x and y frame sizes
		nx = self.planes[0][0].header['NAXIS1']
		ny = self.planes[0][0].header['NAXIS2']

		# Convert sci images to units of e- counts
		stack,ustack = self._getStack()
		pair_stack	= np.zeros((nx,ny,npairs))
		pair_ustack = np.zeros((nx,ny,npairs))

		# Iterate over all pairs to compute A-B for each pair
		for i,pair in enumerate(pairs):
			pair_stack[:,:,i] = stack[:,:,pair[0]] - stack[:,:,pair[1]]
			pair_ustack[:,:,i] = np.sqrt(ustack[:,:,pair[0]]**2 + ustack[:,:,pair[1]]**2)

		
		return pair_stack,pair_ustack

	# Returns stacks for A nods and B nods
	def _makeSingleBeamStacks(self,pairs):
		npairs = len(pairs)
		nx = self.planes[0][0].header['NAXIS1']
		ny = self.planes[0][0].header['NAXIS2']

		stack,ustack = self._getStack()
		stackA	= np.zeros((nx,ny,npairs))
		ustackA = np.zeros((nx,ny,npairs))
		stackB	= np.zeros((nx,ny,npairs))
		ustackB = np.zeros((nx,ny,npairs))
		
		for i,pair in enumerate(pairs):
			stackA[:,:,i] = stack[:,:,pair[0]]
			ustackA[:,:,i] = ustack[:,:,pair[0]]
			stackB[:,:,i] = stack[:,:,pair[1]]
			ustackB[:,:,i] = ustack[:,:,pair[1]]
		return stackA,ustackA,stackB,ustackB

	# Determine which images are A nods, and which are B nods
	# Returns a list of arrays, where each array is a AB pair
	# e.g. [(0, 1), (3, 2), (4, 5), (7, 6), (8, 9), (11, 10)] for 12 images in a standard ABBA pattern
	def _getPairs(self,RAs,DECs,FileNums):
		nexp = len(RAs)				   
		assert nexp % 2 ==0, "There must be an even number of exposures"

		AorB = []
		RA_A  = RAs[0]
		DEC_A = DECs[0]

		RA_B  = RAs[1]
		DEC_B = DECs[1]
		
		dist_a = np.sqrt(((np.array(RAs)-RA_A)*3600)**2+((np.array(DECs)-DEC_A)*3600)**2)
		dist_b = np.sqrt(((np.array(RAs)-RA_B)*3600)**2+((np.array(DECs)-DEC_B)*3600)**2)
		
		for i in range(0,nexp):
			if (dist_a[i] < dist_b[i]):
				AorB.append('A')
			else:
				AorB.append('B')

		print('')
		print('%7s %12s %12s %5s' % ('Exp.', 'RA', 'DEC', 'AorB'))
		print('%7s %12s %12s %5s' % (str('-'*7), str('-'*12), str('-'*12), str('-'*5)))
		for i in range(0,len(RAs)):
			print('%7i %12f %12f %5s' % (FileNums[i], RAs[i], DECs[i], AorB[i]))
		print('\n')

		#LOOK FOR PROBLEMS IN RA AND DEC FROM HEADERS
		problem = False
		if (0.0 in RAs or 0.0 in DECs):
			problem = True

		TotAorB = Counter((np.array(AorB))) 
		TotA, TotB = TotAorB['A'], TotAorB['B']
		if (TotA != TotB):
			problem = True

		if (problem):
			print ('Warning!! - There may be a problem with the pointing in these exposures.')

			while (True):
				yn_correct = input('Is the AB pattern shown above correct? [yes/no]: ').lower()
				if (yn_correct == 'no'):
					do_fix = True
					break
				elif (yn_correct == 'yes'):
					do_fix = False
					break
				else:
					print("Please type 'yes' or 'no'.")

			while (do_fix):
				yn_abba = input('Assume an ABBA pattern? [yes/no]: ').lower()
				if (yn_abba == 'yes'):
					if (nexp % 4 == 0):
						AorB = 'ABBA'*int(nexp/4)
						do_fix = False
					else:
						print('The number of exposures must be a multiple of 4 to assume an ABBA pattern.')
				elif (yn_abba == 'no'):
					do_setnod = True
					while (do_setnod):
						yn_pattern = input('Set nod pattern? [yes/no]: ').lower()
						if (yn_pattern == 'yes'):
							nod_pattern = input('Import pattern, e.g. ABBA BA ABBA \n (case or spaces do not matter): ').upper().replace(" ", "")
							AorB = list(nod_pattern)
							do_setnod, do_fix = False, False
						elif (yn_pattern == 'no'):
							print('Please fix image headers and run the pipeline again.')
							sys.exit()
						else:
							print("Please type 'yes' or 'no'.")
				else:
					print("Please type 'yes' or 'no'.")

			print('')
			print ('## Revised AB Pairs ##')
			print ('======================')
			print ('%7s %12s %12s %5s' % ('Exp.', 'RA', 'DEC', 'AorB'))
			print ('%7s %12s %12s %5s' % (str('-'*7), str('-'*12), str('-'*12), str('-'*5)))	
			for i in range(0,len(RAs)):
				print  ('%7i %12f %12f %5s' % (FileNums[i], RAs[i], DECs[i], AorB[i]))
			print ('\n')
		
		ii = 0
		pairs = []
		while ii<nexp:
			if AorB[ii]!=AorB[ii+1]:
				if AorB[ii]=='A':
					pair = (ii,ii+1)
				else:
					pair = (ii+1,ii)
				pairs.append(pair)
				ii+=2
			else:
				ii+=1
		#print (pairs)
		return pairs

## Class to prepare a single order, resulting in a 2D spectral order, saved in SPEC2D/
## 1) create A-B images for a specific order, 2) shift images to match pos of the first image,
## 3) combine images, and 4) rectify the combined image by fitting a polynomial
class Order():
	def __init__(self,Nod,onum=1,trace=None,write_path=None, targetType = 'bright'):
		self.type = 'order'
		self.header	 = Nod.header
		self.setting = Nod.setting
		self.echelle = Nod.echelle
		self.crossdisp = Nod.crossdisp
		self.airmass = Nod.airmass
		self.Envi	 = Nod.Envi
		self.tname = Nod.tname
		self.targetType = targetType
		# Order being processed
		self.onum	 = onum
		# Spatial range of order in units of pixels
		self.yrange = self.Envi.getYRange(self.setting,onum)
		
		# Crop A-B images in the spatial direction to extract a single order
		self.stack	= Nod.TargetStack[self.yrange[0]:self.yrange[1],:,:]
		self.ustack = Nod.UTargetStack[self.yrange[0]:self.yrange[1],:,:]
		nexp = len(self.stack[0,0,:])
		self.write_path = write_path
		
		self.lampstack	= Nod.LampStack[self.yrange[0]:self.yrange[1],:,:]
		self.ulampstack = Nod.ULampStack[self.yrange[0]:self.yrange[1],:,:]
		
		# CB ADD - from MLB below
		## are there nans in original ustack? yes
		print(np.isnan(self.stack).any())
		print(np.isnan(self.ustack).any())
		if np.isnan(self.ustack).any():
			## make median filter to replace nans in ustack
			img_medfilt = medfilt(self.ustack, 3)
			badimg = np.where(np.isnan(self.ustack))
			img_rectc = self.ustack
			# Replace bad pixels with median filtered results
			img_rectc[badimg] = img_medfilt[badimg]
			self.ustack = img_rectc
			
			#sigma = np.median(self.uimage_rect)
			#badimg = np.abs((self.image_rect - img_medfilt)) / sigma > 7.0
		
		
		# Save intermediate products to plots folder
		plots_path = self.write_path.replace('SPEC2D', 'PLOTS/') + 'A-B-order' + str(onum) + '.fits'
		cube = np.zeros((self.stack.shape[2], self.stack.shape[0], self.stack.shape[1]))
		for slice in range(self.stack.shape[2]):
			im = self.stack[:, :, slice]
			cube[slice, :, :] = im
		pf.writeto(plots_path, cube, overwrite=True)

		# Save intermediate products to plots folder
		plots_path = self.write_path.replace('SPEC2D', 'PLOTS/') + 'lamps-order' + str(onum) + '.fits'
		cube = np.zeros((self.lampstack.shape[2], self.lampstack.shape[0], self.lampstack.shape[1]))
		for slice in range(self.lampstack.shape[2]):
			im = self.lampstack[:, :, slice]
			cube[slice, :, :] = im
		pf.writeto(plots_path, cube, overwrite=True)
		
		
		# Get y offsets between the same order from each image
		offsets1, offsets2 = self._findYOffsets()
		# Shift based on offsets, so that when combining the images, the order from different images fit each other
		self.stack1	 = self._yShift(offsets1,self.stack)
		self.ustack1 = self._yShift(offsets1,self.ustack)
		self.stack2	 = self._yShift(offsets2,self.stack)
		self.ustack2 = self._yShift(offsets2,self.ustack)

		# Save intermediate products to plots folder
		plots_path = self.write_path.replace('SPEC2D', 'PLOTS/') + 'yShifted-posTrace-order'+str(onum)+'.fits'
		cube = np.zeros((self.stack1.shape[2], self.stack1.shape[0], self.stack1.shape[1]))
		for slice in range(self.stack1.shape[2]):
			im = self.stack1[:, :, slice]
			cube[slice, :, :] = im
		pf.writeto(plots_path, cube, overwrite=True)

		plots_path = self.write_path.replace('SPEC2D', 'PLOTS/') + 'yShifted-negTrace-order'+str(onum)+'.fits'
		cube = np.zeros((self.stack2.shape[2], self.stack2.shape[0], self.stack2.shape[1]))
		for slice in range(self.stack2.shape[2]):
			im = self.stack2[:, :, slice]
			cube[slice, :, :] = im
		pf.writeto(plots_path, cube, overwrite=True)
		
		'''
		plots_path = self.write_path.replace('SPEC2D', 'PLOTS/') + 'yShifted-posTrace-unc-order'+str(onum)+'.fits'
		cube = np.zeros((self.ustack1.shape[2], self.ustack1.shape[0], self.ustack1.shape[1]))
		for slice in range(self.ustack1.shape[2]):
			im = self.ustack1[:, :, slice]
			cube[slice, :, :] = im
		pf.writeto(plots_path, cube, overwrite=True)

		plots_path = self.write_path.replace('SPEC2D', 'PLOTS/') + 'yShifted-negTrace-unc-order'+str(onum)+'.fits'
		cube = np.zeros((self.ustack2.shape[2], self.ustack2.shape[0], self.ustack2.shape[1]))
		for slice in range(self.ustack2.shape[2]):
			im = self.ustack2[:, :, slice]
			cube[slice, :, :] = im
		pf.writeto(plots_path, cube, overwrite=True)
		'''
		
		# Combine A-B images with weighted average
		self.image1, self.uimage1 = self._collapseOrder(stack=self.stack1, ustack=self.ustack1)
		self.image2, self.uimage2 = self._collapseOrder(stack=self.stack2, ustack=self.ustack2)
		self.lampimage, self.ulampimage = self._collapseOrder(stack=self.lampstack, ustack=self.ulampstack)
		

		# Fitting method for fitTrace function. Use the more sensitive Gauss method for faint companions
		if self.targetType == 'faint':
			fitMethod = 'Gauss'
		else:
			fitMethod = 'FFT'
		# Fit polynomial to rectify A-B image
		if trace is None:
			yr1, trace1 = self.fitTrace(self.image1, 1, fitMethod=fitMethod)
			yr2, trace2 = self.fitTrace(self.image2, 2, fitMethod=fitMethod)
			yrs, traces	 = [yr1, yr2], [trace1, trace2]
  
		images, uimages = [self.image1, self.image2], [self.uimage1, self.uimage2]
		lampimages, ulampimages = [self.lampimage, self.lampimage], [self.ulampimage, self.ulampimage]
		
		#plt.figure()
		#plt.imshow(self.image1, cmap='gray')
		#sh = self.image1.shape
		#plt.axhline(y=sh[0]/2-1)
		#plt.colorbar()
		#plt.figure()
		#plt.imshow(self.image2, cmap='gray')
		#sh = self.image2.shape
		#plt.axhline(y=sh[0]/2-1)
		#plt.colorbar()
		#plt.show()
		
		#plt.figure()
		#plt.imshow(self.uimage1, cmap='gray')
		#sh = self.uimage1.shape
		#plt.axhline(y=sh[0]/2-1)
		#plt.colorbar()
		#plt.figure()
		#plt.imshow(self.uimage2, cmap='gray')
		#sh = self.uimage2.shape
		#plt.axhline(y=sh[0]/2-1)
		#plt.colorbar()
		#plt.show()

		# Use trace to correct for non-linearity in combined order, i.e. rectifying
		self.image_rect, self.uimage_rect = self.yRectify(images, uimages, yrs, traces)
		self.lampimage_rect, self.ulampimage_rect = self.yRectify(lampimages, ulampimages, yrs, traces)
		
		
		self.sh = self.image_rect.shape
		#print self.sh

		# Subtract median of each column from the same column - SKY SUBTRACTION
		self._subMedian()

		# Repeat above for individual A and B images
		self.beamSkyStacks	= Nod.beamSkyStacks	 
		self.beamUSkyStacks = Nod.beamUSkyStacks

		beamSkyRect, beamUSkyRect = [], []
		for beamSkyStack, beamUSkyStack in zip(self.beamSkyStacks, self.beamUSkyStacks):
			# Crop sky images in spatial direction to extract a single order
			beamSkyStackO  = beamSkyStack[self.yrange[0]:self.yrange[1],:,:]
			beamUSkyStackO = beamUSkyStack[self.yrange[0]:self.yrange[1],:,:]
			
			beamSkyStackO1	= self._yShift(offsets1, beamSkyStackO)
			beamUSkyStackO1 = self._yShift(offsets1, beamUSkyStackO)
			beamSkyStackO2	= self._yShift(offsets2, beamSkyStackO)
			beamUSkyStackO2 = self._yShift(offsets2, beamUSkyStackO)

			beamSky1, beamUSky1 = self._collapseOrder(stack=beamSkyStackO1, ustack=beamUSkyStackO1)
			beamSky2, beamUSky2 = self._collapseOrder(stack=beamSkyStackO2, ustack=beamUSkyStackO2)

			beamSkies, beamUSkies = [beamSky1, beamSky2], [beamUSky1, beamUSky2]
			beam_sky_rect, beam_usky_rect = self.yRectify(beamSkies, beamUSkies, yrs, traces)
			beamSkyRect.append(beam_sky_rect), beamUSkyRect.append(beam_usky_rect)

		self.sky_rect  = beamSkyRect[0]
		self.usky_rect = beamUSkyRect[0]

		# print ('sky_rect shape')
		# print (self.sky_rect.shape)

		# plt.imshow(self.sky_rect, cmap='gray')
		# plt.colorbar()
		# plt.text(30, 30, 'sky rect initial')
		# plt.show()

		## These lines create the sky. This is done by using values from both the A and B image
		## Specifically, (B-A) image is made, and where the counts are < 2*error,
		## the counts from A are replace with the counts from B. This ensures the spectra in A is replaced by sky counts from B
		ssubs = np.where(beamSkyRect[1]-beamSkyRect[0]<2.*beamUSkyRect[0])
		self.sky_rect[ssubs]  = beamSkyRect[1][ssubs]
		self.usky_rect[ssubs] = beamUSkyRect[1][ssubs]

		# plt.imshow(self.sky_rect, cmap='gray')
		# plt.colorbar()
		# plt.text(30, 30, 'sky rect final')
		# plt.show()
		
		#### compare this to the one without bad pixels
		#plt.figure()
		#plt.imshow(self.image_rect, cmap='gray')
		#plt.colorbar()
		#plt.text(30, 30, 'SPEC 2D')
		#plt.figure()
		#plt.imshow(self.sky_rect, cmap='gray')
		#plt.colorbar()
		#plt.text(30, 30, 'Sky original')
		
		#plt.show()
		#assert(1==0)
		
		#plt.figure()
		#plt.imshow(self.image_rect, cmap='gray')
		#plt.colorbar()
		#plt.text(30, 30, 'SPEC 2D ')
		
		#plt.figure()
		#plt.imshow(self.lampimage_rect, cmap='gray')
		#plt.colorbar()
		#plt.text(30, 30, 'SPEC 2D ')
		#plt.show()
	
		'''
		### my thought for why this isn't working for NIRSPEC2 is because counts go up to 20,000-25,000 now, 
		### instead of 4000 and 
		### the top of the peak might just be more than 7sigma over the median, so it gets replaces with
		### values around it, and this cuts off the top of the PSF -- CB 8/14/19
		# MLB ADD - need to remove bad pixels in 2D SPEC
		img_medfilt = medfilt(self.image_rect, 3)
		sigma = np.median(self.uimage_rect)
		plt.figure()
		plt.imshow(np.abs((self.image_rect - img_medfilt)) / sigma, cmap='gray')
		plt.colorbar()
		badimg = np.abs((self.image_rect - img_medfilt)) / sigma > 7.0
		img_rectc = self.image_rect

		# Replace bad pixels with median filtered results
		img_rectc[badimg] = img_medfilt[badimg]

		self.image_rect = img_rectc
		
		skyimg_medfilt = medfilt(self.sky_rect, 3)
		sigma = np.median(self.usky_rect)
		badsky = np.abs((self.sky_rect - skyimg_medfilt)) / sigma > 7.0
		sky_rectc = self.sky_rect
		sky_rectc[badsky] = skyimg_medfilt[badsky]

		self.sky_rect = sky_rectc
		'''
		#img_medfilt = gaussian_filter1d(self.image_rect, 3, axis=0)
		##sigma = np.median(self.uimage_rect)
		#plt.figure()
		#plt.imshow(img_medfilt, cmap='gray')
		#plt.colorbar()
		#plt.figure()
		#plt.imshow(np.abs((self.image_rect - img_medfilt)) , cmap='gray')
		#plt.colorbar()
		#plt.show()
		#badimg = np.abs((self.image_rect - img_medfilt)) / sigma > 7.0
		#img_rectc = self.image_rect
		
		#plt.figure()
		#plt.imshow(self.image_rect, cmap='gray')
		#plt.colorbar()
		#plt.text(30, 30, 'SPEC 2D')
		#plt.show()
		#assert(1==0)
		
		if write_path:
			self.file = self.writeImage(path=write_path)
			self.filelamp = self.writeImage(filename=write_path+'/lamps.fits',lamp=True)
		
		
	def _collapseOrder(self, stack=None, ustack=None):

		sh = stack.shape
		int_stack = np.zeros(sh)
		int_ustack = np.zeros(sh)
		
		index = np.arange(sh[0])
   
		masked_stack = ma.masked_invalid(stack)
		masked_ustack = ma.masked_invalid(ustack)
		
		image = ma.average(masked_stack,2,weights=1./masked_ustack**2)
		uimage = np.sqrt(ma.mean(masked_ustack**2,2)/ma.count(masked_ustack,2))
		
		return image, uimage

	# Shifts images in y (spatial direction) so that they are aligned for combining
	def _yShift(self,offsets,stack):

		sh = stack.shape
		internal_stack = np.zeros(sh)

		index = np.arange(sh[0])

		for plane in np.arange(sh[2]):
			for i in np.arange(sh[1]):
				col = np.interp(index-offsets[plane],index,stack[:,i,plane])
				internal_stack[:,i,plane] = col

		return internal_stack

	# Compute offsets in y (spatial direction) between different images
	def _findYOffsets(self, kwidth=100):
		sh = self.stack.shape

		yr1 = (0, sh[0] / 2 - 1)  # Postive trace
		yr2 = (sh[0] / 2, sh[0] - 1)  # Negative trace
		yrs = [yr1,yr2]

		offsets1 = np.empty(0)
		offsets2 = np.empty(0)
		# Iterate twice, for positive trace and neg trace
		for i in range(0, len(yrs)):
			yr = yrs[i]
			yindex = np.arange(yr[0],yr[1]+1) #Top or Bottom of Order
			yindex = yindex.astype(np.int64)
			# Use the middle A-B image as the kernel by default
			kernel_image = self.stack[:, :, 0]

			# Median of values along the wavelength direction, over the central 200 pixels
			kernel_o = np.median(kernel_image[yindex,int(sh[1]/2)-kwidth:int(sh[1]/2)+kwidth],1)

			# Median of those medians
			kernel_med = np.median(kernel_o)
			kernel = np.subtract(kernel_o, kernel_med)

			# Iterate over each A-B image
			centroids = []
			for j in range(0,sh[2]):
				image = self.stack[:,:,j]
				# Median y values for a window of size kwidth
				profile_o = np.median(image[yindex,int(sh[1]/2)-kwidth:int(sh[1]/2)+kwidth],1)
				# Median of those medians
				profile_med = np.median(profile_o)
				# Normalize y profile by subtracting off the median
				profile = np.subtract(profile_o, profile_med)
				# Cross correlation between kernel and image - ???
				cc = fp.ifft(fp.fft(kernel)*np.conj(fp.fft(profile)))
				cc_sh = fp.fftshift(cc)
				cen = calc_centroid(cc_sh).real - yindex.shape[0]/2.
				centroids.append(cen)

				if (i == 0):
					offsets1 = np.append(offsets1,cen)
				if (i == 1):
					offsets2 = np.append(offsets2,cen)
				
		print ('')
		print ('Img	 A		B')
		print ('------')
		for k in range(0,sh[2]):
			print (('%i	 %5.2f %5.2f') % (k+1, offsets1[k], offsets2[k]))
		print ('')

		return offsets1, offsets2

	# Added option to tune kwidth by finding the smallest residuals to the polyfit
	def fitTrace(self, image, OneOrTwo, kwidth=100, fitMethod = 'FFT'):

		sh = image.shape
		print(sh)
		# A nod (1) or B nod (2) position
		if (OneOrTwo == 1):
			print ('positive trace')
			yr = (0,sh[0]/2-1)
			#yr = (0,150)
		if (OneOrTwo == 2):
			print ('negative trace')
			yr = (sh[0]/2,sh[0]-1)
			#yr = (130,sh[0]-1)


		yindex = np.arange(yr[0], yr[1] + 1)  # Top or Bottom of Order
		yindex = yindex.astype(np.int64)
		kernel = np.median(image[yindex, int(sh[1] / 2) - kwidth:int(sh[1] / 2) + kwidth], 1)
		centroids = []
		totals = []
		# Iterate over each col (wavelength direction)
		# Check
		for i in np.arange(sh[1]):
			col_med = np.median(image[yindex, i])
			total = np.abs((image[yindex, i] - col_med).sum())
			cc = fp.ifft(fp.fft(kernel) * np.conj(fp.fft(image[yindex, i] - col_med)))
			cc_sh = fp.fftshift(cc)
			centroid = calc_centroid(cc_sh).real - yindex.shape[0] / 2.
			centroids.append(centroid)
			totals.append(total)

		centroids = np.array(centroids)
		totals = np.array(totals)
		median_totals = np.nanmedian(totals)

		xindex = np.arange(sh[1])
		gsubs = np.where((np.isnan(centroids) == False) & (totals > median_totals * 0.25) & (totals < median_totals * 1.75))
		
		
		#### ADDED BY CB to cut in x range to get better fit to trace (12/17/2019)
		#print(np.shape(gsubs),gsubs,type(gsubs),gsubs[0])
		## cut off in x direction
		gsubs_cut = gsubs[0][np.where(gsubs[0]>400)]
		gsubs_cut = gsubs_cut[np.where(gsubs_cut<1500)]
		gsubs = tuple([np.array(gsubs_cut)],)
		#print(gsubs_cut,gsubs)
		
		
		#centroids[gsubs] = median_filter(centroids[gsubs], size=50)
		centroids[gsubs] = median_filter(centroids[gsubs], size=10)
		# Fit 3-order polynomial to order
		fit_params = np.polyfit(xindex[gsubs], centroids[gsubs], 3, full=True)
		# residuals
		res = fit_params[1][0]
		# coefficients of poly
		coeffs = fit_params[0]
		poly = np.poly1d(coeffs)
		print ('FFT polyfit residuals: ' + str(res))
		
		## show fit to trace
		#plt.figure()
		#plt.plot(xindex[gsubs], centroids[gsubs])
		#plt.plot(xindex[gsubs],poly(xindex[gsubs]))
		#plt.show()
		#assert(1==0)
		
		# Gaussian fit
		def gauss(x, a, x0, sigma):
			y = a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))
			return y

		if fitMethod == 'Gauss':

			# 2nd iteration using fft trace (from polyfit) as guesses for Gaussian fit
			gauss_centroids = []
			# Iterate over columns
			gauss_fail = 0
			model_centroids = []

			for i in np.arange(sh[1]):
				model_centroids.append(poly[i])

			for i in np.arange(sh[1]):
				counts = image[yindex, i]
				# Better guess?
				sigma = 1.5

				sigma_stats = sigmaclip(counts, low=6.0, high=6.0)
				upper_l = sigma_stats[2]
				lower_l = sigma_stats[1]

				# Positive trace
				if OneOrTwo == 1:
					# Five tallest peaks
					peaks = counts[counts.argsort()[-5:][::-1]]
				# Negative trace
				else:
					# Five lowest peals
					peaks = counts[counts.argsort()[:5]]

				peaks_clipped = []
				for p in peaks:
					if not (p > upper_l or p < lower_l):
						peaks_clipped.append(p)

				peaks_clipped = np.asarray(peaks_clipped)

				# If there are any valid peaks
				if len(peaks_clipped) > 0:
					# Fit a Gaussian for each of the peaks
					sigma_list = []
					a_list = []
					x0_list = []
					fit_list = []

					for peak_height in peaks_clipped:
						try:
							peak_index = np.argwhere(counts == peak_height)[0][0]
							mean = yindex[peak_index]

							initial_guess = [peak_height, mean, sigma]

						except:
							pass

						try:
							popt, pcov = curve_fit(gauss, yindex, counts, p0=initial_guess)
							#  Generate y-data based on the fit.
							count_fit = gauss(yindex, *popt)
							count_fit = count_fit.tolist()
							#print (count_fit)

							fit_list.append(count_fit)
							x0_list.append(popt[1])
							a_list.append(popt[0])
							sigma_list.append(popt[2])

						# Failed to fit Gaussian curve for this peak option at this pixel
						except:
							pass

					sigma_list = np.asarray(sigma_list)
					a_list = np.asarray(a_list)
					x0_list = np.asarray(x0_list)
					fit_list = np.asarray(fit_list)


					if len(sigma_list) > 0:
						bad_ind = []
						for j in range(len(sigma_list)):
							others = np.delete(x0_list, [j])
							offset = abs(x0_list[j] - np.median(others))

							## Double-check the peaks found by doing:
							# 1 Remove fits where sigma was either too big or too small
							# 2 Remove fits where height of Gauss fit is too large
							# 3 Remove fits where peak is far away from other peaks
							if (sigma_list[j] > 20.0 or sigma_list[j] < 1.0 or a_list[j] > upper_l or a_list[j] < lower_l or offset > (yindex[0]+yindex[-1])/4.0):
								bad_ind.append(j)

						# If there were bad indices
						if len(bad_ind) > 0:
							sigma_list = np.delete(sigma_list, bad_ind)
							a_list = np.delete(a_list, bad_ind)
							x0_list = np.delete(x0_list, bad_ind)

							fit_list = np.delete(fit_list, bad_ind, 0)

						# After deleting bad indices, check again if there are values left
						if len(sigma_list) > 0:
							# Pick the peak with the largest sigma
							final_ind = np.argmax(sigma_list)
							gauss_centroid = x0_list[final_ind]
							count_fit = fit_list[final_ind]
							gauss_centroids.append(gauss_centroid)

							# Plot gaussian fit
							if self.onum == 5 or self.onum == 6:
								if OneOrTwo == 2:
									if i % 50 == 0:
										# Create a plot of our work, showing both the data and the fit.
										fig, ax = plt.subplots()
										ax.plot(yindex, counts, 'b+:', label='data')
										ax.plot(yindex, count_fit, 'ro:', label='fit; final params: a = ' + str(
											a_list[final_ind]) + ' x0 = ' + str(
											gauss_centroid) + ' sigma = ' + str(sigma_list[final_ind]))

										ax.legend()
										ax.set_xlabel('y index')
										ax.set_ylabel('count')
										#plt.show()
										plt.savefig(
											self.write_path.replace('SPEC2D', 'PLOTS/') + 'gaussFit-order' + str(
												self.onum) + '-xpixel' + str(i + 1) + '.png')
										plt.close()

						else:
							gauss_fail = gauss_fail + 1
							gauss_centroids.append(np.nan)

					else:
						gauss_fail = gauss_fail+1
						gauss_centroids.append(np.nan)

				else:
					gauss_fail = gauss_fail + 1
					gauss_centroids.append(np.nan)

			print ('Gauss failed: ' + str(gauss_fail))
			gauss_centroids = np.array(gauss_centroids)
			print (gauss_centroids)
			print (len(gauss_centroids))

			gauss_gsubs = np.where(np.isnan(gauss_centroids) == False)
			print (len(gauss_gsubs[0]))

			gauss_centroids[gauss_gsubs] = median_filter(gauss_centroids[gauss_gsubs], size=50)
			gauss_cen_med = np.median(gauss_centroids[gauss_gsubs])

			# Fit 3-order polynomial to order
			g_fit_params = np.polyfit(xindex[gauss_gsubs],gauss_cen_med-gauss_centroids[gauss_gsubs], 3, full=True)
			# residuals
			min_g_res = g_fit_params[1][0]
			# Polynomial fit for Gauss fit
			g_coeffs = g_fit_params[0]
			g_poly = np.poly1d(g_coeffs)


			print ('Gauss polyfit residuals: '+str(min_g_res))
			fig = plt.figure(figsize=(18, 12))
			ax1 = fig.add_subplot(1, 1, 1)
			# FFT and Gauss fit models
			ax1.plot(xindex[gsubs], poly(xindex[gsubs]), label='FFT fit', ls=':', lw=4, color='black')
			ax1.plot(xindex[gauss_gsubs], g_poly(xindex[gauss_gsubs]), label='Gauss fit', ls='-.', lw=4, color='black')
			# Data points
			ax1.scatter(xindex[gsubs], centroids[gsubs], label='FFT centroids', color='red', s=2)
			#ax1.scatter(xindex[gauss_gsubs], gauss_centroids[gauss_gsubs], label='Gauss centroids', color='blue', s=2)
			#ax1.scatter(xindex[gauss_gsubs], ((yindex[0]+yindex[-1]) / 2.0) - gauss_centroids[gauss_gsubs], label='gauss centroids', color='blue')
			ax1.scatter(xindex[gauss_gsubs], gauss_cen_med - gauss_centroids[gauss_gsubs],label='gauss centroids', color='blue')
			ax1.set_ylim(-20, 15)
			ax1.set_xlabel('x pixel (wavelength direction)')
			ax1.set_ylabel('centroid position')
			ax1.legend(loc='upper left')
			#plt.show()

			ax1.set_title('fft residuals: ' + str(res) + ', kwidth: ' + str(kwidth))
			if OneOrTwo == 1:
				whichHalf = 'pos'
			else:
				whichHalf = 'neg'
			plt.savefig(self.write_path.replace('SPEC2D', 'PLOTS/') + 'fitTrace-order'+str(self.onum)+'-'+whichHalf+'.png')

			# print ('yr and g_poly: ')
			# print ('------------------------')
			# print yr
			# print g_poly

			## If using Gaussian method
			return yr, g_poly

		## With FFT method
		else:
			return yr,poly

	def yRectify(self,images,uimages,yrs,traces): 

		sh = images[0].shape
		image_rect = np.zeros(sh)
		uimage_rect = np.zeros(sh)
		
		for yr,trace,image,uimage in zip(yrs,traces,images,uimages):
			#yr,trace,image,uimage = yrs[1],traces[1],images[1],uimages[1] 
			index = np.arange(yr[0],yr[1]+1)
			index = index.astype(np.int64)
			# Iterate over columns
			for i in np.arange(sh[1]):
				#i = 1750
				col = np.interp(index-trace(i),index,image[index,i])
				#plt.plot(col)
				#plt.plot(image[index,i])
				#plt.show()
				image_rect[index,i] = col
				col = np.interp(index-trace(i),index,uimage[index,i])
				uimage_rect[index,i] = col

		return image_rect,uimage_rect
		
	# Subtract median of each column from that column
	def _subMedian(self):
		self.image_rect = self.image_rect-np.median(self.image_rect,axis=0)

	def writeImage(self,filename=None,path='.',lamp=False):

		if filename is None:
			time   = self.header['UTC']
			time   = time.replace(':','')
			time   = time[0:4]
			date   = self.header['DATE-OBS']
			date   = date.replace('-','')
			object = self.header['OBJECT']
			object = object.replace(' ','')
			filename = path+'/'+object+'_'+date+'_'+time+'_order'+str(self.onum)+'.fits'
		
		if not lamp:
			hdu	 = pf.PrimaryHDU(self.image_rect)
			uhdu = pf.ImageHDU(self.uimage_rect)
			sky_hdu = pf.ImageHDU(self.sky_rect)
			usky_hdu = pf.ImageHDU(self.usky_rect)
		else:
			hdu	 = pf.PrimaryHDU(self.lampimage_rect)
			uhdu = pf.ImageHDU(self.ulampimage_rect)


		#print (type(hdu.data))

		hdu.header['SETNAME'] = (self.setting, 'Setting name')
		hdu.header['ECHLPOS'] = (self.echelle, 'Echelle position')
		hdu.header['DISPPOS'] = (self.crossdisp, 'Cross disperser position')
		hdu.header['ORDER'] = (str(self.onum),'Order number')
		
		if not lamp:		
			hdulist = pf.HDUList([hdu,uhdu,sky_hdu,usky_hdu])
		else:
			hdulist = pf.HDUList([hdu,uhdu])

		hdulist.writeto(filename,overwrite=True)

		return filename

## Class to extract 1D spectra for positive and negative parts of an order respectively
## Returns spectra from positive and negative part, the associated errors, and the spectra of the sky
class Spec1D():
	def __init__(self,Order,sa=True,write_path=None,nirspec=1):
		self.Order	 = Order
		self.setting = Order.setting
		self.echelle = Order.echelle
		self.crossdisp = Order.crossdisp
		self.airmass = Order.airmass
		self.onum	 = Order.onum
		self.header	 = Order.header
		self.Envi	 = Order.Envi
		self.sh		 = Order.sh[1]
		self.sa		 = sa
		self.writepath = write_path

		# Normalized 1D PSF
		PSF = self.getPSF()
		self.disp = self.Envi.getDispersion(self.setting,self.onum)

		# Wavelength guesses for each pixel
		self.wave_pos = self.waveGuess(beam='pos')
		self.wave_neg = self.waveGuess(beam='neg')

		# print ('wave pos & wave neg')
		# print (self.wave_pos)
		# print (self.wave_neg)

		# Get 1D flux for pos and neg orders, as well as sky background
		self.flux_pos,self.uflux_pos,self.flux_neg,self.uflux_neg,self.sky_pos,self.sky_neg,self.usky_pos,self.usky_neg,self.lamp_pos,self.lamp_neg,self.ulamp_pos,self.ulamp_neg = self.extract(PSF,method='single')
		
		#plt.figure()
		#plt.plot(self.flux_pos)
		##plt.figure()
		#plt.plot(self.flux_neg)
		#plt.show()
		#assert(1==0)
		
		# print ('flux pos and flux neg')
		# print (self.flux_pos).shape
		# print (self.flux_neg).shape
		#
		# print (self.flux_pos)
		# print (self.flux_neg)

		# ???
		if sa:
			#self.sa_pos,self.usa_pos,self.sa_neg,self.usa_neg = self.SpecAst(PSF)
			self.sa_pos = self.flux_pos
			self.usa_pos = self.uflux_pos
			self.sa_neg = self.flux_neg
			self.usa_neg = self.uflux_neg

		if write_path:
			self.file = self.writeSpec(path=write_path)

	# Create normalized 1D PSF in a window from x=300 to x=900 pixels
	def getPSF(self,range=(300,900)):
		# PSF is the median of values in the wavelength direction
		PSF = np.median(self.Order.image_rect[:,range[0]:range[1]],1)
		
		## plot all PSFs
		#plt.figure()
		#print(np.shape(self.Order.image_rect))
		#ord = np.transpose(self.Order.image_rect)
		#print(np.shape(ord))
		#for i in np.arange(np.shape(ord)[1]):
		#	plt.plot(ord[i])
		
		#print (PSF.shape)
		npsf = PSF.size
		# Divide by half the sum of counts (because pos and neg flux each account for half the flux)
		PSF_norm = PSF/PSF[:int(npsf/2)-1].sum()
		#print (PSF_norm.shape)
		return PSF_norm 

	### fit Gaussian to PSF
	def Gaussian(self,x, a, x0, sigma):
		y = a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))
		return y
	
	def gaussresid(self,params, x, y):
		a, x0, sigma = params['A'].value,params['x0'].value,params['sigma'].value		
		fit = self.Gaussian(x, a, x0, sigma)
		resid = np.sum((fit - y)**2)
				
		return resid		
			
	def GaussScaling(self,x,gauss,A):
		return A*gauss
	
	def GaussScalingResid(self,params,x,mediangauss,individualgauss):
		a = params['A'].value
		fit = self.GaussScaling(x,mediangauss,a)
		resid = np.sum((fit-individualgauss)**2.)
		return resid
	
	def FitGaussiantoPSF(self,y,psf,PosOrNeg='pos'):
		## fit and plot gaussian on PSF
		params = Parameters()
		if PosOrNeg=='pos':
			params.add('A',value=np.max(psf), vary=True)
			params.add('x0',value=np.array(y)[np.argmax(psf)], vary=True)
		if PosOrNeg=='neg':
			params.add('A',value=np.min(psf), vary=True)
			params.add('x0',value=np.array(y)[np.argmin(psf)], vary=True)
		params.add('sigma',value=4, vary=True)
		
		result = minimize(self.gaussresid, params, args=(y, psf), method='powell') #'powell' 'least-sq'
		params = result.params	
		gauss = self.Gaussian(y, params["A"].value, params["x0"].value, params["sigma"].value)
		
		return gauss, params["A"].value, params["x0"].value, params["sigma"].value
	
	# Extract 1D spectra using method from Horne1986
	def extract(self,PSF,method='single'):
		# Placeholders
		pixsig = 1.
		# Frame size
		sh	 = self.sh
		npsf = PSF.size
		#plt.figure()
		#plt.plot(PSF)
		
		flux_pos  = np.zeros(sh)
		flux_neg  = np.zeros(sh)

		uflux_pos = np.zeros(sh)
		uflux_neg = np.zeros(sh)

		sky_pos	 = np.zeros(sh)
		sky_neg	 = np.zeros(sh)		   

		usky_pos  = np.zeros(sh)
		usky_neg  = np.zeros(sh)	
		
		lamp_pos  = np.zeros(sh)
		lamp_neg  = np.zeros(sh)
		
		ulamp_pos = np.zeros(sh)
		ulamp_neg = np.zeros(sh)	

		# Load rectified A-B images
		im = self.Order.image_rect
		uim = self.Order.uimage_rect
		sky = self.Order.sky_rect
		usky = self.Order.usky_rect
		lamp = self.Order.lampimage_rect
		ulamp = self.Order.ulampimage_rect		
				
		#plt.figure()
		#for i in allpsfs_pos:
		#	plt.plot(i)
		#plt.show()
		#assert(1==0)
		

		pos_psf = PSF[:int(npsf / 2) - 1]
		pos_y = np.linspace(0, npsf / 2 - 1, PSF[:int(npsf / 2) - 1].size, dtype=int)
		#
		# print (pos_psf.shape)
		# print pos_y.shape
		
		pos_gauss, A_pos, x0_pos, sigma_pos = self.FitGaussiantoPSF(pos_y,pos_psf,PosOrNeg='pos')
		
		neg_psf = PSF[int(npsf/2):-1]		
		neg_y = np.linspace(npsf/2, npsf-1, PSF[int(npsf/2):-1].size, dtype=int)

		neg_gauss, A_neg, x0_neg, sigma_neg = self.FitGaussiantoPSF(neg_y,neg_psf,PosOrNeg='neg')
		
		# Save PSF plots
		fig, ax = plt.subplots()
		ax.plot(pos_y,pos_gauss, 'gray',label=r'%.3f $\pm$ %.3f' % (x0_pos, sigma_pos))
		ax.plot(pos_y, pos_psf, 'b+:', label='PSF of positive trace')
		ax.legend()
		ax.set_xlabel('y index')
		ax.set_ylabel('normalized flux')
		plt.savefig(self.writepath.replace('SPEC1D', 'PLOTS/') + 'posPSF-order' + str(
				self.onum) + '.png')
		plt.close()
		
		fig, ax = plt.subplots()
		ax.plot(neg_y,neg_gauss, 'gray',label=r'%.3f $\pm$ %.3f' % (x0_neg, sigma_neg))
		ax.plot(neg_y, neg_psf, 'b+:', label='PSF of negative trace')
		ax.legend()
		ax.set_xlabel('y index')
		ax.set_ylabel('normalized flux')
		plt.savefig(self.writepath.replace('SPEC1D', 'PLOTS/') + 'negPSF-order' + str(
			self.onum) + '.png')
		plt.close()	

		### Identify bad pixels -- CB 8/2019
		#posinfo, neginfo = self.identifybadpixels_version2(im, PSF, pos_y, x0_pos, pos_psf, neg_y, x0_neg, neg_psf)
		#badpix_pos,badranges_pos,groupedposranges,posbounds = posinfo
		#badpix_neg,badranges_neg,groupednegranges,negbounds = neginfo
		
		# MLB(/CB) ADD - need to remove bad pixels in (rectified) 2D SPEC
		img_medfilt = []
		for i in im:
			img_medfilt.append(medfilt(i, 3))
		img_medfilt = np.array(img_medfilt)
		sigma = np.median(uim)
		#plt.figure()
		#plt.imshow(np.abs((self.image_rect - img_medfilt)) / sigma, cmap='gray')
		#plt.colorbar()
		badimg = np.abs((im - img_medfilt)) / sigma > 7.0
		# Replace bad pixels with median filtered results
		im[badimg] = img_medfilt[badimg]
	
		skyimg_medfilt = []
		for i in sky:
			skyimg_medfilt.append(medfilt(i, 3))
		skyimg_medfilt = np.array(skyimg_medfilt)	
		sigma = np.median(usky)
		badsky = np.abs((sky - skyimg_medfilt)) / sigma > 7.0
		sky[badsky] = skyimg_medfilt[badsky]

		
		# Iterate over wavelength direction
		for i in np.arange(sh):
			#There is something wrong with the optimal weights here. Will have to check on it later. For 
			#now, uniform weights make very little difference for high S/N spectra.
			
			'''
			### fix bad pixels identified above - CB 8/2019 
			## if multiple bad pixels in a row, interpolate from ends
			if np.any(badpix_pos==i) and np.any(badranges_pos==i): 
				for n,m in enumerate(groupedposranges):
					if np.any(m==i):
						ind = n
				bound1, bound2 = np.transpose(im[:int(npsf/2)-1,posbounds[ind][0]]), np.transpose(im[:int(npsf/2)-1,posbounds[ind][1]])
				skybound1, skybound2 = np.transpose(sky[:int(npsf/2)-1,posbounds[ind][0]]), np.transpose(sky[:int(npsf/2)-1,posbounds[ind][1]])
				imgs, skyimgs = [], []
				for row,row2, srow1, srow2 in zip(bound1,bound2,skybound1,skybound2):
					img = np.interp(i, posbounds[ind], [row,row2])
					imgs.append(img)
					skyimg = np.interp(i, posbounds[ind], [srow1,srow2])
					skyimgs.append(skyimg)
				im[:int(npsf/2)-1,i] = np.array(imgs)	
				sky[:int(npsf/2)-1,i] = np.array(skyimgs)	
			if np.any(badpix_pos==i) and not np.any(badranges_pos==i):	
				im[:int(npsf/2)-1,i] = np.mean([im[:int(npsf/2)-1,i-1],im[:int(npsf/2)-1,i+1]],axis=0)
				sky[:int(npsf/2)-1,i] = np.mean([sky[:int(npsf/2)-1,i-1],sky[:int(npsf/2)-1,i+1]],axis=0)
			
			if np.any(badpix_neg==i) and np.any(badranges_neg==i): 
				for n,m in enumerate(groupednegranges):
					if np.any(m==i):
						ind = n
				bound1, bound2 = np.transpose(im[int(npsf/2):-1,negbounds[ind][0]]), np.transpose(im[int(npsf/2):-1,negbounds[ind][1]])
				skybound1, skybound2 = np.transpose(sky[int(npsf/2):-1,negbounds[ind][0]]), np.transpose(sky[int(npsf/2):-1,negbounds[ind][1]])
				imgs, skyimgs = [], []
				for row,row2,srow, srow2 in zip(bound1,bound2,skybound1,skybound2):
					img = np.interp(i, negbounds[ind], [row,row2])
					imgs.append(img)		
					skyimg = np.interp(i, negbounds[ind], [srow,srow2])
					skyimgs.append(skyimg)		
				im[int(npsf/2):-1,i] = np.array(imgs)	
				sky[int(npsf/2):-1,i] = np.array(skyimgs)	
			# if single bad pix, replace with average of two surrounding columns
			if np.any(badpix_neg==i) and not np.any(badranges_neg==i):	
				im[int(npsf/2):-1,i] = np.mean([im[int(npsf/2):-1,i-1],im[int(npsf/2):-1,i+1]],axis=0)
				sky[int(npsf/2):-1,i] = np.mean([sky[int(npsf/2):-1,i-1],sky[int(npsf/2):-1,i+1]],axis=0)
			'''

			# Compute fluxes with optimal extraction algorithm (Horne1986)
			## optimal weights
			flux_pos[i] = (pos_psf*im[:int(npsf/2)-1,i]/uim[:int(npsf/2)-1,i]**2).sum() / (pos_psf**2/uim[:int(npsf/2)-1,i]**2).sum()
			flux_neg[i] = (neg_psf*im[int(npsf/2):-1,i]/uim[int(npsf/2):-1,i]**2).sum() / (neg_psf**2/uim[int(npsf/2):-1,i]**2).sum()
			## uniform weights
			#flux_pos[i] = (PSF[:int(npsf/2)-1]*im[:int(npsf/2)-1,i]).sum() / (PSF[:int(npsf/2)-1]**2).sum()
			#flux_neg[i] = (PSF[int(npsf/2):-1]*im[int(npsf/2):-1,i]).sum() / (PSF[int(npsf/2):-1]**2).sum()

			# Error in fluxes
			uflux_pos[i] = np.sqrt(1.0/(pos_psf**2/uim[:int(npsf/2)-1,i]**2.).sum())
			uflux_neg[i] = np.sqrt(1.0/(neg_psf**2/uim[int(npsf/2):-1,i]**2.).sum())

			# Same as above for the sky, but variance is zero (sky assumed to be noise-free)
			sky_pos[i] = (pos_psf*sky[:int(npsf/2)-1,i]).sum() / (pos_psf**2).sum()
			sky_neg[i] = -(neg_psf*sky[int(npsf/2):-1,i]).sum() / (neg_psf**2).sum()
			usky_pos[i] = np.sqrt(1.0/(pos_psf**2/usky[:int(npsf/2)-1,i]**2.).sum())
			usky_neg[i] = np.sqrt(1.0/(neg_psf**2/usky[int(npsf/2):-1,i]**2.).sum())
			
			# Lamps
			lamp_pos[i] = (pos_psf*lamp[:int(npsf/2)-1,i]/ulamp[:int(npsf/2)-1,i]**2).sum() / (pos_psf**2/ulamp[:int(npsf/2)-1,i]**2).sum()
			lamp_neg[i] = (neg_psf*lamp[int(npsf/2):-1,i]/ulamp[int(npsf/2):-1,i]**2).sum() / (neg_psf**2/ulamp[int(npsf/2):-1,i]**2).sum()
			ulamp_pos[i] = np.sqrt(1.0/(pos_psf**2/ulamp[:int(npsf/2)-1,i]**2.).sum())
			ulamp_neg[i] = np.sqrt(1.0/(neg_psf**2/ulamp[int(npsf/2):-1,i]**2.).sum())
			

		#plt.figure()
		#plt.imshow(im, cmap='gray')
		#plt.colorbar()
		#plt.text(30, 30, 'Bad Pixel Corrected')
		#
		#plt.figure()
		#plt.imshow(sky, cmap='gray')
		#plt.colorbar()
		#plt.text(30, 30, 'Sky')


		#print ('PSF and im')
		# print (PSF[:npsf/2-1])
		# print (im[:npsf/2-1,512])
		#
		# print (PSF[npsf / 2:-1])
		# print (im[npsf / 2:-1, 512])

		#### MICHAEL'S ADDITION TO THE CODE
		#print(np.shape(flux_pos))
		#print(np.shape(PSF))
		#model_image = PSF[:int(npsf/2)-1][:, np.newaxis] * flux_pos
		#residuals = im[:int(npsf/2)-1] - model_image
		#z_scores = residuals / uim[:int(npsf/2)-1]
		#np.save("/home/cbuzard/Pipeline/01_Reduction/reducs_test_michael/2019_04_08/HD187123/residuals_3.npy", residuals)
		#np.save("/home/cbuzard/Pipeline/01_Reduction/reducs_test_michael/2019_04_08/HD187123/z_scores_3.npy", z_scores)
		#np.save("/home/cbuzard/Pipeline/01_Reduction/reducs_test_michael/2019_04_08/HD187123/im3_3.npy", im[:int(npsf/2)-1])
		#assert(1==0)
        
		
		# Fill in faulty value with 1. for fluxes, 1000 for error, so they will have tiny weight
		flux_pos = ma.masked_invalid(flux_pos)
		flux_pos = ma.filled(flux_pos,1.)
		uflux_pos = ma.masked_invalid(uflux_pos)
		uflux_pos = ma.filled(uflux_pos,1000.)
		sky_pos = ma.masked_invalid(sky_pos)
		sky_pos = ma.filled(sky_pos,1.)
		usky_pos = ma.masked_invalid(usky_pos)
		usky_pos = ma.filled(usky_pos,1000.)
		lamp_pos = ma.masked_invalid(lamp_pos)
		lamp_pos = ma.filled(lamp_pos,1.)
		ulamp_pos = ma.masked_invalid(ulamp_pos)
		ulamp_pos = ma.filled(ulamp_pos,1000.)
		
		flux_neg = ma.masked_invalid(flux_neg)
		flux_neg = ma.filled(flux_neg,1.)
		uflux_neg = ma.masked_invalid(uflux_neg)
		uflux_neg = ma.filled(uflux_neg,1000.)
		sky_neg = ma.masked_invalid(sky_neg)
		sky_neg = ma.filled(sky_neg,1.)
		usky_neg = ma.masked_invalid(usky_neg)
		usky_neg = ma.filled(usky_neg,1000.)
		lamp_neg = ma.masked_invalid(lamp_neg)
		lamp_neg = ma.filled(lamp_neg,1.)
		ulamp_neg = ma.masked_invalid(ulamp_neg)
		ulamp_neg = ma.filled(ulamp_neg,1000.)

		# What is this???
		sky_pos_cont = self._fitCont(self.wave_pos,sky_pos)
		sky_neg_cont = self._fitCont(self.wave_neg,sky_neg)

		# print ('flux_pos')
		# print (flux_pos.shape)
		# print ('sky_pos')
		# print (sky_pos.shape)

		time = self.header['UTC']
		time = time.replace(':', '')
		time = time[0:4]
		date = self.header['DATE-OBS']
		date = date.replace('-', '')
		object = self.header['OBJECT']
		object = object.replace(' ', '')

		# Plot positive flux with initial wavelength guess
		fig, ax = plt.subplots()
		ax.plot(self.wave_pos, flux_pos, drawstyle ='steps-mid',label='Pos flux w/ initial wavelength fit')
		ax.plot(self.wave_pos, sky_pos, drawstyle='steps-mid', label='Pos sky w/ initial wavelength fit')
		ax.legend(loc='upper left')
		filename = self.writepath + '/' + object + '_' + date + '_' + time + '_spec1d' + str(self.onum) + '_pos' + '.png'
		plt.savefig(filename)
		plt.close()

		fig, ax = plt.subplots()
		ax.plot(self.wave_neg, flux_neg, drawstyle='steps-mid', label='Neg flux w/ initial wavelength fit')
		ax.plot(self.wave_neg, sky_neg, drawstyle='steps-mid', label='Neg sky w/ initial wavelength fit')
		ax.legend(loc='upper left')
		filename = self.writepath + '/' + object + '_' + date + '_' + time + '_spec1d' + str(
			self.onum) + '_neg' + '.png'
		#print (filename)
		plt.savefig(filename)
		plt.close()

		# Why is sky_pos_cont subtracted??? why is nothing subtracted for flux_pos?
		return flux_pos,uflux_pos,flux_neg,uflux_neg,sky_pos-sky_pos_cont,sky_neg-sky_neg_cont,usky_pos,usky_neg,lamp_pos,lamp_neg,ulamp_pos,ulamp_neg

	"""
	def identifybadpixels(self,im, PSF, pos_y, x0_pos, pos_psf, neg_y, x0_neg, neg_psf):
		### CB 08/2019
		
		sh	 = self.sh
		npsf = PSF.size
		### get individual psfs
		allpsfs = np.transpose(im)
		allpsfs_pos = [indpsf[:int(npsf / 2) - 1] for indpsf in allpsfs]
		allpsfs_neg = [indpsf[int(npsf/2):-1] for indpsf in allpsfs]
		
		### Identify bad pixels in the pos and neg traces - added by CB 08/2019
		pos_redchi = []
		#data_model_diff = []
		#### for order 0 
		##allpsfs_pos = allpsfs_pos[275:]
		for num,individualpsf in enumerate(allpsfs_pos):
			## find closest x value to x0
			miny = np.argmin(np.abs(pos_y-x0_pos))
			minind = miny - 10
			maxind = miny + 11
			# find scaling
			params = Parameters()
			params.add('A',value=np.max(individualpsf), vary=True)
			result = minimize(self.GaussScalingResid, params, args=(pos_y[minind:maxind],pos_psf[minind:maxind],individualpsf[minind:maxind]), method='powell') 
			params, redchi = result.params, result.redchi	
			pos_redchi.append(redchi)
			#good_individual_gauss = self.GaussScaling(pos_y, pos_psf, params["A"].value)
			#data_model_diff.append(individualpsf[minind:maxind]-good_individual_gauss[minind:maxind])
			#if num == 1035 or num == 1408 or num == 1621:
			#	good_individual_gauss = self.GaussScaling(pos_y, pos_psf, params["A"].value)
			#	plt.figure()
			#	plt.plot(pos_y,individualpsf)
			#	plt.plot(pos_y,good_individual_gauss)
		#print(np.argsort(pos_redchi)[::-1][:10],np.sort(pos_redchi)[::-1][:10])
		#badpix_pos, badval_pos = np.argsort(pos_redchi)[::-1]+275, np.sort(pos_redchi)[::-1]	
		badpix_pos, badval_pos = np.argsort(pos_redchi)[::-1], np.sort(pos_redchi)[::-1]
		
		neg_redchi = []
		### for order 0
		##allpsfs_neg = allpsfs_neg[:1900]
		for num,individualpsf in enumerate(allpsfs_neg):
			## find closest x value to x0
			miny = np.argmin(np.abs(neg_y-x0_neg))
			minind = miny - 10
			maxind = miny + 11
			# find scaling
			params = Parameters()
			params.add('A',value=np.min(individualpsf), vary=True)
			result = minimize(self.GaussScalingResid, params, args=(neg_y[minind:maxind],neg_psf[minind:maxind],individualpsf[minind:maxind]), method='powell') 
			params, redchi = result.params, result.redchi	
			neg_redchi.append(redchi)
			#if num == 1471 or num == 385 or num == 656:
			#	good_individual_gauss = self.GaussScaling(neg_y, neg_psf, params["A"].value)
			#	plt.figure()
			#	plt.plot(neg_y,individualpsf)
			#	plt.plot(neg_y,good_individual_gauss)
		#print(np.argsort(neg_redchi)[::-1][:10],np.sort(neg_redchi)[::-1][:10])
		badpix_neg, badval_neg = np.argsort(neg_redchi)[::-1], np.sort(neg_redchi)[::-1]		
		
		## difference is less than 1% of last difference
		for posind, ii in enumerate(badval_pos):
			if np.abs(badval_pos[posind+1] - ii) <= 0.01*ii and np.abs(badval_pos[posind+2] - ii) <= 0.01*ii: break
		for negind, ii in enumerate(badval_neg):
			if np.abs(badval_neg[negind+1] - ii) <= 0.01*ii and np.abs(badval_neg[negind+2] - ii) <= 0.01*ii: break
		
		## normalize then go until we reach 0.01
		#badval_pos = badval_pos/np.max(badval_pos)
		#badval_neg = badval_neg/np.max(badval_neg)
		#for posind, ii in enumerate(badval_pos):
		#	if ii <= 0.01: break
		#for negind, ii in enumerate(badval_neg):
		#	if ii <= 0.01: break
			
		badpix_pos = badpix_pos[:posind]
		badpix_neg = badpix_neg[:negind]	
		badpix_pos, badpix_neg = np.sort(badpix_pos), np.sort(badpix_neg)
		
		
		#plt.figure()
		#plt.plot(badval_pos)
		##plt.plot(scipy.signal.savgol_filter(badval_pos, 51, 3))
		#plt.yscale("log")
		#plt.axvline(x=posind)
		#plt.figure()
		#plt.plot(badval_neg)
		##plt.plot(scipy.signal.savgol_filter(badval_neg, 51, 3))
		#plt.yscale("log")
		#plt.axvline(x=negind)
		#plt.figure()
		#plt.plot(np.diff(np.diff(badval_pos)))
		##plt.plot(scipy.signal.savgol_filter(np.diff(badval_pos), 5, 2))
		#plt.axvline(x=posind)
		#plt.figure()
		#plt.plot(np.diff(np.diff(badval_neg)))
		##plt.plot(scipy.signal.savgol_filter(np.diff(badval_neg), 5, 3))
		#plt.axvline(x=negind)		
		
		#plt.figure()
		#for i in data_model_diff[:10]:
		#	plt.plot(i)
		#plt.figure()
		#for i in data_model_diff[500:600]:
		#	plt.plot(i)	
		#plt.figure()		
		#plt.plot(np.sort(np.sum(data_model_diff,axis=1))[::-1])		
				
		print(posind,negind)
		print(badpix_pos, badpix_neg)

		badranges_pos, badranges_neg = [], []

		for k, g in groupby(enumerate(badpix_pos), lambda ix : ix[0] - ix[1]):
			a=list(map(itemgetter(1), g))
			if len(a) > 1:
				badranges_pos += a
		for k, g in groupby(enumerate(badpix_neg), lambda ix : ix[0] - ix[1]):
			a=list(map(itemgetter(1), g))
			if len(a) > 1:
				badranges_neg += a

		#print(badranges_pos, badranges_neg)

		groupedposranges = np.array([np.array(list(g)) for k, g in groupby(badranges_pos, key=lambda i,j=count(): i-next(j))])
		groupednegranges = np.array([np.array(list(g)) for k, g in groupby(badranges_neg, key=lambda i,j=count(): i-next(j))])
		posbounds, negbounds = [], []
		for i in groupedposranges:
			posbounds.append([i[0]-1,i[-1]+1])
		for i in groupednegranges:
			negbounds.append([i[0]-1,i[-1]+1])
		
		#print(groupedposranges,groupednegranges)
		#print(posbounds,negbounds)
		
		posinfo = [badpix_pos,badranges_pos,groupedposranges,posbounds]
		neginfo = [badpix_neg,badranges_neg,groupednegranges,negbounds]
		return posinfo, neginfo	
		
	
	def identifybadpixels_version2(self,im, PSF, pos_y, x0_pos, pos_psf, neg_y, x0_neg, neg_psf):
		### CB 08/2019
		
		## find closest x value to x0
		miny_pos = np.argmin(np.abs(pos_y-x0_pos))
		minind_pos = miny_pos - 10
		maxind_pos = miny_pos + 11
		## find closest x value to x0
		miny_neg = np.argmin(np.abs(neg_y-x0_neg))
		minind_neg = miny_neg - 10
		maxind_neg = miny_neg + 11
		
		sh	 = self.sh
		npsf = PSF.size
		### get individual psfs
		allpsfs = np.transpose(im)
		allpsfs_pos = [indpsf[:int(npsf / 2) - 1][minind_pos:maxind_pos] for indpsf in allpsfs]
		allpsfs_neg = [indpsf[int(npsf/2):-1][minind_neg:maxind_neg] for indpsf in allpsfs]
		#allpsfs_pos = allpsfs_pos/np.max(allpsfs_pos)
		#allpsfs_neg = allpsfs_neg/np.min(allpsfs_neg)
		
		## Divide by half the sum of counts (because pos and neg flux each account for half the flux)
		#allpsfs_pos = [pos/PSF[:int(npsf/2)-1].sum() for pos in allpsfs_pos]
		#allpsfs_neg = [neg/PSF[:int(npsf/2)-1].sum() for neg in allpsfs_neg]
		
		#plt.figure()
		#for num, i in enumerate(allpsfs_pos):
		#	#if np.max(i)-np.median(i) > 3*np.std(i):
		#	outoftrace = np.delete(allpsfs[num][:int(npsf / 2) - 1],np.arange(minind_pos,maxind_pos))
		#	if np.max(i)-np.median(outoftrace) > 7*np.std(outoftrace):
		#		plt.plot(i)
		#plt.figure()
		#for i in allpsfs_neg:
		#	#if np.max(i)-np.median(i) > 3*np.std(i):
		#	outoftrace = np.delete(allpsfs[num][int(npsf/2):-1],np.arange(minind_neg,maxind_neg))
		#	if np.max(i)-np.median(outoftrace) > 7*np.std(outoftrace):
		#		plt.plot(i)
		
		### Identify bad pixels in the pos and neg traces - added by CB 08/2019
		#pos_redchi = []
		modeltodatadiff = []
		for num,individualpsf in enumerate(allpsfs_pos):
			#if np.any(np.abs(individualpsf) > 0.01):
			#if np.max(individualpsf)-np.median(individualpsf) > 3*np.std(individualpsf):
			#outoftrace = np.delete(allpsfs[num][:int(npsf / 2) - 1],np.arange(minind_pos,maxind_pos))
			#if np.max(individualpsf)-np.median(outoftrace) > 7*np.std(outoftrace):
			#individualpsf = individualpsf/np.max(individualpsf)
			#pos_psf_norm = pos_psf[minind_pos:maxind_pos]
			#pos_psf_norm = pos_psf_norm /np.max(pos_psf_norm)
			# find scaling
			params = Parameters()
			params.add('A',value=np.max(individualpsf), vary=True)
			result = minimize(self.GaussScalingResid, params, args=(pos_y[minind_pos:maxind_pos],pos_psf[minind_pos:maxind_pos],individualpsf), method='powell') 
			params, redchi = result.params, result.redchi	
			#pos_redchi.append(redchi)
			#if num == 1250 or num == 879 or num == 880:
			good_individual_gauss = self.GaussScaling(pos_y[minind_pos:maxind_pos], pos_psf[minind_pos:maxind_pos], params["A"].value)	
			#diff = np.abs(np.sum(individualpsf-pos_psf_norm))
			diff = np.abs(np.sum(individualpsf-good_individual_gauss))
			#if np.any(np.abs(individualpsf-pos_psf_norm) > 0.1):
			modeltodatadiff.append(diff)
			## plot
			#fig, axes = plt.subplots(2,1,sharex=True)
			#axes[0].plot(pos_y[minind_pos:maxind_pos],individualpsf,label='Data')
			#axes[0].plot(pos_y[minind_pos:maxind_pos],pos_psf_norm,label='Model')
			#axes[0].legend()
			#axes[1].plot(pos_y[minind_pos:maxind_pos],individualpsf-pos_psf_norm)
			#plt.title(diff)
			#plt.savefig(self.writepath.replace('SPEC1D', 'PLOTS/') + 'PSFs_order5/pos_{}.png'.format(num))
			#plt.close()
			#else:
			#	modeltodatadiff.append(0)
			#else:
			#	## we don't want to correct sat. tellurics as bad pixels
			#	modeltodatadiff.append(0)	
				
		#print(len(modeltodatadiff))
		#badpix_pos, badval_pos = np.argsort(pos_redchi)[::-1], np.sort(pos_redchi)[::-1]
		badpix_pos, badval_pos = np.argsort(modeltodatadiff)[::-1], np.sort(modeltodatadiff)[::-1]
		#badval_pos = badval_pos/np.max(badval_pos)
		#print(badpix_pos, badval_pos)
	

		#neg_redchi = []
		modeltodatadiff = []
		for num,individualpsf in enumerate(allpsfs_neg):
			#if np.any(np.abs(individualpsf) > 0.01):
			#if np.max(individualpsf)-np.median(individualpsf) > 3*np.std(individualpsf):
			#outoftrace = np.delete(allpsfs[num][int(npsf/2):-1],np.arange(minind_neg,maxind_neg))
			#if np.max(np.abs(individualpsf))-np.median(outoftrace) > 7*np.std(outoftrace):
			#	individualpsf = individualpsf/np.max(individualpsf)
			#	neg_psf_norm = neg_psf[minind_neg:maxind_neg]
			#	neg_psf_norm = neg_psf_norm /np.min(neg_psf_norm)
			#	diff = np.abs(np.sum(individualpsf-neg_psf_norm))
			#	if np.any(np.abs(individualpsf-neg_psf_norm) > 0.1):
			#		modeltodatadiff.append(diff)
			#	else:
			#		modeltodatadiff.append(0)	
			#else:
			#	modeltodatadiff.append(0)
			## find scaling
			params = Parameters()
			params.add('A',value=np.min(individualpsf), vary=True)
			result = minimize(self.GaussScalingResid, params, args=(neg_y[minind_neg:maxind_neg],neg_psf[minind_neg:maxind_neg],individualpsf), method='powell') 
			params, redchi = result.params, result.redchi	
			#neg_redchi.append(redchi)
			good_individual_gauss = self.GaussScaling(neg_y[minind_neg:maxind_neg],neg_psf[minind_neg:maxind_neg], params["A"].value)	
			diff = np.abs(np.sum(individualpsf-good_individual_gauss))
			modeltodatadiff.append(diff)
			## figure
			#fig, axes = plt.subplots(2,1,sharex=True)
			#axes[0].plot(neg_y[minind_neg:maxind_neg],individualpsf,label='Data')
			#axes[0].plot(neg_y[minind_neg:maxind_neg],neg_psf_norm,label='Model')
			#axes[0].legend()
			#axes[1].plot(neg_y[minind_neg:maxind_neg],individualpsf-neg_psf_norm)
			#plt.title(diff)
			#plt.savefig(self.writepath.replace('SPEC1D', 'PLOTS/') + 'PSFs_order5/neg_{}.png'.format(num))
			#plt.close()
			

		#print(len(modeltodatadiff))	
		badpix_neg, badval_neg = np.argsort(modeltodatadiff)[::-1], np.sort(modeltodatadiff)[::-1]		
		#print(badpix_neg, badval_neg)

		## difference is less than 1% of last difference
		#for posind, ii in enumerate(badval_pos):
		#	if np.abs(badval_pos[posind+1] - ii) <= 0.01*ii and np.abs(badval_pos[posind+2] - ii) <= 0.01*ii: break
		#	#if ii < 0.01*np.max(badval_pos): break
		#for negind, ii in enumerate(badval_neg):
		#	if np.abs(badval_neg[negind+1] - ii) <= 0.01*ii and np.abs(badval_neg[negind+2] - ii) <= 0.01*ii: break
		#	#if ii < 0.01*np.max(badval_neg): break
		
		#from scipy.signal import savgol_filter
		#yhat = savgol_filter(badval_pos, 15, 2)
		#
		#y_spl = UnivariateSpline(np.arange(len(badval_pos)),yhat,s=0,k=4)
		#y_spl_2d = y_spl.derivative(n=2)
		#secderiv_pos = y_spl_2d(np.arange(len(badval_pos)))

		plt.figure()
		plt.plot(badval_pos)
		plt.axhline(np.median(badval_pos))
		plt.axhline(np.median(badval_pos)+2*np.std(badval_pos))
		plt.figure()
		plt.plot(badval_neg)
		plt.axhline(np.median(badval_neg))
		plt.axhline(np.median(badval_neg)+2*np.std(badval_neg))

		posind  = np.argmin(np.abs(badval_pos-(np.median(badval_pos)+2*np.std(badval_pos))))
		negind  = np.argmin(np.abs(badval_neg-(np.median(badval_neg)+2*np.std(badval_neg))))

		badpix_pos = badpix_pos[:posind]
		badpix_neg = badpix_neg[:negind]	
		
		print(posind,negind)
		print(badpix_pos, badpix_neg)
		badpix_pos, badpix_neg = np.sort(badpix_pos), np.sort(badpix_neg)
		
		print(badpix_pos, badpix_neg)

		badranges_pos, badranges_neg = [], []

		for k, g in groupby(enumerate(badpix_pos), lambda ix : ix[0] - ix[1]):
			a=list(map(itemgetter(1), g))
			if len(a) > 1:
				badranges_pos += a
		for k, g in groupby(enumerate(badpix_neg), lambda ix : ix[0] - ix[1]):
			a=list(map(itemgetter(1), g))
			if len(a) > 1:
				badranges_neg += a
				
		groupedposranges = np.array([np.array(list(g)) for k, g in groupby(badranges_pos, key=lambda i,j=count(): i-next(j))])
		groupednegranges = np.array([np.array(list(g)) for k, g in groupby(badranges_neg, key=lambda i,j=count(): i-next(j))])
		posbounds, negbounds = [], []
		for i in groupedposranges:
			posbounds.append([i[0]-1,i[-1]+1])
		for i in groupednegranges:
			negbounds.append([i[0]-1,i[-1]+1])
		
		
		posinfo = [badpix_pos,badranges_pos,groupedposranges,posbounds]
		neginfo = [badpix_neg,badranges_neg,groupednegranges,negbounds]
		return posinfo, neginfo	
		
	"""
	
	def _fitCont(self,wave,spec):
		bg_temp = 210. #K
		niter = 2

		cont = self.bb(wave*1e-6,bg_temp)
		#print ('cont')
		#print (cont)
		gsubs = np.where(np.isfinite(spec))
		for i in range(niter):
			norm = np.median(spec[gsubs])
			# print ('norm')
			# print (norm)
			norm_cont = np.median(cont[gsubs])
			# normalize by specific intensity from a BB continuum?
			cont *= norm/norm_cont
			gsubs = np.where(spec<cont)

		return cont

	# Computes specific intensity: intensity at each wavelength
	def bb(self,wave,T):
		# speed of light
		cc = constants.c
		# Planck constant
		hh = constants.h
		# Boltzmann constant
		kk = constants.k
		
		blambda = 2.*hh*cc**2/(wave**5*(np.exp(hh*cc/(wave*kk*T))-1.))

		return blambda
		

	def SpecAst(self,PSF,method='centroid',width=5):
		'''
		The uncertainty on the centroid is:

				 SUM_j([j*SUM_i(F_i)-SUM_i(i*F_i)]^2 * s(F_j)^2)
		s(C)^2 = ------------------------------------------------
								[SUM_i(F_i)]^4

		
		'''
		
		#Guesstimated placeholder
		aper_corr = 1.4
		posloc = np.argmax(PSF)
		negloc = np.argmin(PSF)

		sa_pos = np.zeros(self.sh)
		sa_neg = np.zeros(self.sh)

		usa_pos = np.zeros(self.sh)
		usa_neg = np.zeros(self.sh)

		im = self.Order.image_rect
		uim = self.Order.uimage_rect

		for i in np.arange(self.sh):
			index = np.arange(width*2+1)-width

			# First calculate SUM_i(F_i)
			F_pos = (im[posloc-width:posloc+width+1,i]).sum() 
			F_neg = (im[negloc-width:negloc+width+1,i]).sum()
			#print index, negloc, width, i, np.shape(im)
			# then SUM_i(i*F_i)
			iF_pos = (index*im[posloc-width:posloc+width+1,i]).sum()
			iF_neg = (index*im[negloc-width:negloc+width+1,i]).sum()

			sa_pos[i] = iF_pos/F_pos
			sa_neg[i] = iF_neg/F_neg
	   
			# Now propagate the error
			uF_pos = uim[posloc-width:posloc+width+1,i]
			uF_neg = uim[negloc-width:negloc+width+1,i]
			usa_pos[i]	= np.sqrt(((index*F_pos - iF_pos)**2 * uF_pos**2).sum())/F_pos**2
			usa_neg[i]	= np.sqrt(((index*F_neg - iF_neg)**2 * uF_neg**2).sum())/F_neg**2

		#NIRSPEC flips the spectrum on the detector (as all echelles do).
		sa_pos[i] = -sa_pos[i]
		sa_neg[i] = -sa_neg[i]

		return sa_pos*aper_corr,usa_pos*aper_corr,sa_neg*aper_corr,usa_neg*aper_corr

	def plot(self):
		plt.plot(self.wave,self.flux_pos,drawstyle='steps-mid')
		plt.plot(self.wave,self.flux_neg,drawstyle='steps-mid')
		plt.show()

	def plotSA(self):
		plt.plot(self.wave,self.sa_pos,drawstyle='steps-mid')
		plt.plot(self.wave,self.sa_neg,drawstyle='steps-mid')
		plt.show()

	# what are the ABC ranges???
	def waveGuess(self,beam='pos'):
		wrange = self.Envi.getWaveRange(self.setting,self.onum)
		index = np.arange(self.sh)
		if beam is 'pos':
			j = 0
		elif beam is 'neg':
			j = 1
		else:
			raise AttributeError
	   
		wave = self.disp['A'][j]+self.disp['B'][j]*index+self.disp['C'][j]*index**2.
		return wave

	def writeSpec(self,filename=None,path='.'):
		
		c1	= pf.Column(name='wave_pos', format='D', array=self.wave_pos)
		c2	= pf.Column(name='flux_pos', format='D', array=self.flux_pos)
		c3	= pf.Column(name='uflux_pos', format='D', array=self.uflux_pos)
		c4	= pf.Column(name='sky_pos', format='D', array=self.sky_pos)
		c5	= pf.Column(name='usky_pos', format='D', array=self.usky_pos)		 
		if self.sa:
			c6	= pf.Column(name='sa_pos', format='D', array=self.sa_pos)
			c7	= pf.Column(name='usa_pos', format='D', array=self.usa_pos)
		c8	= pf.Column(name='wave_neg', format='D', array=self.wave_neg)
		c9	= pf.Column(name='flux_neg', format='D', array=self.flux_neg)
		c10	 = pf.Column(name='uflux_neg', format='D', array=self.uflux_neg)
		c11	 = pf.Column(name='sky_neg', format='D', array=self.sky_neg)
		c12	 = pf.Column(name='usky_neg', format='D', array=self.usky_neg)
		if self.sa:
			c13	 = pf.Column(name='sa_neg', format='D', array=self.sa_neg)
			c14 = pf.Column(name='usa_neg', format='D', array=self.usa_neg)
		c15 = pf.Column(name='lamp_pos', format='D', array=self.lamp_pos)
		c16 = pf.Column(name='lamp_neg', format='D', array=self.lamp_neg)
		c17 = pf.Column(name='ulamp_pos', format='D', array=self.ulamp_pos)
		c18 = pf.Column(name='ulamp_neg', format='D', array=self.ulamp_neg)


		if self.sa:
			coldefs = pf.ColDefs([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18])
		else:
			coldefs = pf.ColDefs([c1,c2,c3,c4,c7,c8,c9,c10,c11,c12])

		#tbhdu = pf.new_table(coldefs)
		tbhdu = pf.BinTableHDU.from_columns(coldefs)

		self.header['SETNAME'] = (self.setting, 'Setting name')
		self.header['ECHLPOS'] = (self.echelle, 'Echelle position')
		self.header['DISPPOS'] = (self.crossdisp, 'Cross disperser position')
		self.header['ORDER'] = (str(self.onum),'Order number')

		self.header['WAVA_POS'] = (self.disp['A'][0],'Dispersion coefficient, positive beam')
		self.header['WAVB_POS'] = (self.disp['B'][0],'Dispersion coefficient, positive beam')
		self.header['WAVC_POS'] = (self.disp['C'][0],'Dispersion coefficient, positive beam')

		self.header['WAVA_NEG'] = (self.disp['A'][1],'Dispersion coefficient, negative beam')
		self.header['WAVB_NEG'] = (self.disp['B'][1],'Dispersion coefficient, negative beam')
		self.header['WAVC_NEG'] = (self.disp['C'][1],'Dispersion coefficient, negative beam')

		self.header['RPOW_POS'] = (self.disp['R'][0],'Resolving power, positive beam')
		self.header['RPOW_NEG'] = (self.disp['R'][1],'Resolving power, negative beam')
		

		hdu = pf.PrimaryHDU(header=self.header)
		thdulist = pf.HDUList([hdu,tbhdu])

		if filename is None:
			time   = self.header['UTC']
			time   = time.replace(':','')
			time   = time[0:4]
			date   = self.header['DATE-OBS']
			date   = date.replace('-','')
			object = self.header['OBJECT']
			object = object.replace(' ','')
			filename = path+'/'+object+'_'+date+'_'+time+'_spec1d'+str(self.onum)+'.fits'

		
		thdulist.writeto(filename,overwrite=True)

		return filename

## Class to calibrate wavelength using rfm
class WaveCal():
	
	def __init__(self,specfile,path='./',am=1., hp=True):
		self.specfile = specfile
		self.path = path
		self.WavePos = rfm.Optimize(specfile,alt=4.145,rpower=29000.,beam='pos',cull=2,am=am)
		self.WaveNeg = rfm.Optimize(specfile,alt=4.145,rpower=29000.,beam='neg',cull=2,am=am) 
		self.file = self._updateWave(specfile,self.WavePos,self.WaveNeg)
		#self.plotWaveFit(hp)
		
	def _updateWave(self,specfile,WavePos,WaveNeg):
		spec1d = pf.open(specfile)
		spec1d[1].data.field('wave_pos')[:] = self.WavePos.getWave()
		spec1d[1].data.field('wave_neg')[:] = self.WaveNeg.getWave()
		oldpath, filename = os.path.split(specfile)
		length = len(filename)
		onum = filename[length - 6]
		filename = filename[0:length - 12]
		filename = filename + 'wave' + onum + '.fits'
		fullpath = self.path + '/' + filename
		spec1d.writeto(fullpath, overwrite=True)

		spec1d_ini_wfit_pos = self.WavePos.getModel()
		fig, ax = plt.subplots()
		ax.plot(spec1d_ini_wfit_pos['wave'], spec1d_ini_wfit_pos['flux'], label='pos flux')
		ax.plot(spec1d_ini_wfit_pos['wave'], spec1d_ini_wfit_pos['Radiance'], label='pos rad')
		ax.legend(loc='upper left')
		plt.savefig(self.path + '/wavePos' + onum + '.png')

		spec1d_ini_wfit_neg = self.WaveNeg.getModel()
		fig, ax = plt.subplots()
		ax.plot(spec1d_ini_wfit_neg['wave'], spec1d_ini_wfit_neg['flux'], label='neg flux')
		ax.plot(spec1d_ini_wfit_neg['wave'], spec1d_ini_wfit_neg['Radiance'], label='neg rad')
		ax.legend(loc='upper left')
		plt.savefig(self.path + '/aveNeg' + onum + '.png')
		
		return fullpath
		
	def plotWaveFit(self, hp=True):
		# NRC (19 JUNE 2014) - EDITED plotWaveFit SO THAT THE WAVE CALIBRATION PLOTS DON'T
		# HAVE TO STOP THE PIPELINE FROM RUNNING BY SETTING 'hp=False'.

		plt.close()
		plt.close()

		trans = self.WavePos.getModel()
		plt.figure(self.specfile+' (Pos)')
		plt.plot(trans['wave'],trans['Radiance'],label='radiance')
		plt.plot(trans['wave'],trans['sky'],label='sky')
		plt.legend(loc='upper left')
		if (hp == True):
			plt.show()
		#else:
			#plt.show(block=False)
			# oldpath,filename = os.path.split(self.specfile)
			# length = len(filename)
			# onum = filename[length-6]
			# plt.savefig(self.path + 'wavePosModel'+onum+'.png')

		trans = self.WaveNeg.getModel()
		plt.figure(self.specfile+' (Neg)')
		plt.plot(trans['wave'],trans['Radiance'],label='radiance')
		plt.plot(trans['wave'],trans['sky'],label='sky')
		plt.legend(loc='upper left')
		if (hp == True):
			plt.show()
		#else:
			# plt.show(block=False)
			# oldpath, filename = os.path.split(self.specfile)
			# length = len(filename)
			# onum = filename[length - 6]
			# plt.savefig(self.path + 'waveNegModel' + onum + '.png')

class CalSpec():
	def __init__(self,scifile,stdfile,shift=0.,dtau=0.0,dowave=True,write_path=None,order=0):

		self.scifile = scifile
		self.stdfile = stdfile

		self.Sci,self.header = readNIRSPEC(scifile)
		self.Std,self.header_std = readNIRSPEC(stdfile)

		self._normalize(self.Sci)
		self._normalize(self.Std)
		self._maskLowTrans()
		self._beers(self.Std,dtau)
		self._addFlux(shift)
		self._addSA(shift)

		if write_path:
			self.file = self.writeSpec(path=write_path,order=order)

	def _maskLowTrans(self,thres=0.25):
		bsubs = np.where(self.Std['flux_pos']<thres)
		self.Std.field('flux_pos')[bsubs] = np.nan
		self.Sci.field('flux_pos')[bsubs] = np.nan
		self.Sci.field('sa_pos')[bsubs] = np.nan
		
		bsubs = np.where(self.Std['flux_neg']<thres)
		self.Std.field('flux_neg')[bsubs] = np.nan
		self.Sci.field('flux_neg')[bsubs] = np.nan
		self.Sci.field('sa_neg')[bsubs] = np.nan

	def _specShift(self,spec,shift):
		nx = spec.size
		index = np.arange(nx)
		return np.interp(index-shift,index,spec)
	#
	def _normalize(self,Spec):
		niter = 3
		
		median = 0.
		gsubs = np.where(Spec['flux_pos']>median)
		for i in np.arange(niter):
			median = np.median(Spec['flux_pos'][gsubs])
			gsubs = np.where(Spec['flux_pos']>median)
		Spec.field('flux_pos')[:] /= median
		Spec.field('uflux_pos')[:] /= median

		median = 0.
		gsubs = np.where(Spec['flux_neg']>median)
		for i in np.arange(niter):
			median = np.median(Spec['flux_neg'][gsubs])
			gsubs = np.where(Spec['flux_neg']>median)
		Spec.field('flux_neg')[:] /= median
		Spec.field('uflux_neg')[:] /= median

	def _beers(self,Spec,dtau):
		Spec.field('flux_pos')[:] = np.exp((1.+dtau)*np.log(Spec['flux_pos']))
		Spec.field('flux_neg')[:] = np.exp((1.+dtau)*np.log(Spec['flux_neg']))

	def _addFlux(self,shift):
		
		sci_pos_resamp = np.interp(self.Std['wave_pos'],self.Sci['wave_pos'],self.Sci['flux_pos'])
		sci_neg_resamp = np.interp(self.Std['wave_pos'],self.Sci['wave_neg'],self.Sci['flux_neg'])
		usci_pos_resamp = np.interp(self.Std['wave_pos'],self.Sci['wave_pos'],self.Sci['uflux_pos'])
		usci_neg_resamp = np.interp(self.Std['wave_pos'],self.Sci['wave_neg'],self.Sci['uflux_neg'])

		std_pos_resamp = self.Std['flux_pos']
		std_neg_resamp = np.interp(self.Std['wave_pos'],self.Std['wave_neg'],self.Std['flux_neg'])
		ustd_pos_resamp = self.Std['uflux_pos']
		ustd_neg_resamp = np.interp(self.Std['wave_pos'],self.Std['wave_neg'],self.Std['uflux_neg'])

		self.flux_pos  = sci_pos_resamp/std_pos_resamp
		self.uflux_pos = np.sqrt((usci_pos_resamp/sci_pos_resamp)**2 + (ustd_pos_resamp/std_pos_resamp)**2)
		self.flux_neg  = sci_neg_resamp/std_neg_resamp
		self.uflux_neg = np.sqrt((usci_neg_resamp/sci_neg_resamp)**2 + (ustd_neg_resamp/std_neg_resamp)**2)

		self.flux = (self.flux_pos+self.flux_neg)/2.
		self.uflux = np.sqrt((self.uflux_pos/self.flux_pos)**2+(self.uflux_neg/self.flux_neg)**2)

		self.wave = self.Std['wave_pos']
		
	def _addSA(self,shift):
		self.sa	 = (self.Sci['sa_pos']+np.interp(self.wave,self.Sci['wave_neg'],self.Sci['sa_neg']))/2.
		self.usa = np.sqrt(self.Sci['usa_pos']**2+np.interp(self.wave,self.Sci['wave_neg'],self.Sci['usa_neg']**2))/2.
		self.sa	 = self._specShift(self.sa,shift)
		self.usa = self._specShift(self.usa,shift)

	def writeSpec(self, filename=None, path='.',order=1):
		c1	= pf.Column(name='wave', format='D', array=self.wave)
		c2	= pf.Column(name='flux', format='D', array=self.flux)
		c3	= pf.Column(name='uflux', format='D', array=self.uflux)
		c4	= pf.Column(name='sa', format='D', array=self.sa)
		c5	= pf.Column(name='usa', format='D', array=self.usa)

		coldefs = pf.ColDefs([c1,c2,c3,c4,c5])

		#tbhdu = pf.new_table(coldefs)
		tbhdu = pf.BinTableHDU.from_columns(coldefs)

		hdu = pf.PrimaryHDU(header=self.header)
		thdulist = pf.HDUList([hdu,tbhdu])

		if filename is None:
			basename = getBaseName(self.header)
			filename = path+'/'+basename+'_calspec'+str(order)+'.fits'

		thdulist.writeto(filename,overwrite=True)

		return filename

	def plotBeams(self):
		plt.plot(self.wave_pos,self.flux_pos)
		plt.plot(self.wave_neg,self.flux_neg)
		plt.show()

	def plotFlux(self):
		plt.plot(self.wave,self.flux,drawstyle='steps-mid')
		plt.show()

	def plotSA(self):
		plt.plot(self.wave,self.sa,drawstyle='steps-mid')
		plt.show()

class SASpec():
	def __init__(self, SA_ortho_file, SA_para_file, write_path=None,order=1):
		self.SA_ortho_file = SA_ortho_file
		self.SA_para_file  = SA_para_file

		self.Spec1,header1 = readNIRSPEC(SA_ortho_file)
		self.Spec2,header2 = readNIRSPEC(SA_para_file)

		self.header = header1
		self._subSA()
		self._addFlux()

		if write_path:
			self.file = self.writeSA(path=write_path,order=order)

	def _subSA(self):
		'''
		Note that we are assuming that the wavescale is close to identical, which
		we anyway need to be the case for SA to work. We could add a check and a warning at this point. 
		'''
		self.wave = (self.Spec1.wave+self.Spec2.wave)/2.
		self.sa = (self.Spec1.sa-self.Spec2.sa)/2.
		self.usa = np.sqrt((self.Spec1.usa**2+self.Spec2.usa**2))/2.

	def _addFlux(self):
		self.flux = (self.Spec1.flux+self.Spec2.flux)/2.
		self.uflux = np.sqrt((self.Spec1.uflux**2+self.Spec2.uflux**2))/2.

	def writeSA(self,filename=None,path='.',order=1):
		c1	= pf.Column(name='wave', format='D', array=self.wave)
		c2	= pf.Column(name='flux', format='D', array=self.flux)
		c3	= pf.Column(name='uflux', format='D', array=self.uflux)
		c4	= pf.Column(name='sa', format='D', array=self.sa)
		c5	= pf.Column(name='usa', format='D', array=self.usa)

		coldefs = pf.ColDefs([c1,c2,c3,c4,c5])

		#tbhdu = pf.new_table(coldefs)
		tbhdu = pf.BinTableHDU.from_columns(coldefs)
		hdu = pf.PrimaryHDU()
		thdulist = pf.HDUList([hdu,tbhdu])

		if filename is None:
			basename = getBaseName(self.header)
			filename = path+'/'+basename+'_saspec'+str(order)+'.fits'

		thdulist.writeto(filename,overwrite=True)

		return filename

# NRC (27 MAY 2014) - EDITED THE Reduction AND SAReduction CLASSES SO THERE IS A 'SettingsFile' ARGUMENT. THE USER
# CAN NOW MORE EASILY CHANGE THE PARAMETERS IN THIS FILE, WHICH DICTATE WHERE THE PIPELINE LOOKS FOR THE ECHELLE 
# ORDERS. THIS IS NECESSARY BECAUSE WHERE THE ORDERS APPEAR ON THE CHIP CAN SHIFT FOR DIFFERENT OBSERVATION EPOCHS,
# AND IF THESE VALUES ARE NOT SET PROPERLY IT CAN CAUSE THE PIPELINE TO FAIL. 

# NRC (19 JUNE 2014) - EDITED THE Reduction AND SAReduction CLASSES SO THAT TARGET AND STANDARD NAMES CAN BE EDITED
# BY THE USER WITH THE sci_tname AND std_tname ARGUMENTS, RESPECTIVELY. THESE NAMES ARE USED IN THE OUTPUT FILES.

# NRC (19 JUNE 2014) - ADDED AN ARGUMENT 'hold_plots' TO THE Reduction AND SAReduction CLASSES SO THAT WHEN SET TO
# True THE PIPELINE WILL NOT BE STOPPED BY THE WAVELENGTH CALIBRATION PLOTS.

class Reduction():
	'''
	Top level basic script for reducing an observation, consisting of a science target and a telluric standard, as well
	as associated calibration files. Output is saved in ../../NAME OF OUTPUT DIR/ut-date/target-name/
	'''

	# initialize file paths and configuration variables
	def __init__(self,flat_range=None, flat_dark_range=None, dark_range=None, lamp_range=None,
				 sci_range=None, std_range=None, path=None, output_path = None, base=None,level1=True,level2=True,
				 shift=0.0, dtau=0.0, save_dark=True, save_flat=True, SettingsFile=None,
				 ut_date = None, sci_tname=None, target_type='bright', std_tname=None, hold_plots=True, nirspec=2, hold=True, incont_fnum = False, **kwargs):

		if (hold == False):
			hold_plots = False

		# Save date and target name for later
		self.ut_date = ut_date
		self.sci_tname = sci_tname
		self.output_path = output_path
		# Whether to save processed darks and flats
		self.save_dark = save_dark
		self.save_flat = save_flat

		# Range of science images
		sci_range1 = sci_range
		self.targetType = target_type

		self.shift = shift
		self.dtau = dtau
		
		self.nirspec = nirspec

		# Path to save Level1 .json file: a summary of L1 output paths

		date_dir = self.output_path+self.ut_date+'/'
		if not os.path.exists(date_dir):
			os.mkdir(date_dir)
		out_dir = date_dir+self.sci_tname+'/'
		if not os.path.exists(out_dir):
			os.mkdir(out_dir)

		self.out_dir = out_dir

		print ('output directory for target')
		print (out_dir)

		self.level1_path = out_dir+'L1FILES'
		self.spec2d_path = out_dir+'SPEC2D'
		self.spec1d_path = out_dir+'SPEC1D'
		self.wave_path = out_dir+'WAVE'
		self.cal1d_path = out_dir+'CAL1D'
		self.plots_path = out_dir+'PLOTS'

		out_path_list = [self.level1_path, self.spec2d_path, self.spec1d_path,
						 self.wave_path, self.cal1d_path, self.plots_path]

		# Create output directories if they do not exist
		for out_path in out_path_list:
			if not os.path.exists(out_path):
				os.mkdir(out_path)

		# Path to .ini file containing instrument settings and details of the spectral orders etc.
		self.SettingsFile = SettingsFile

		# Create lists of file names for each file type
		if flat_dark_range:
			self.flat_dark_names = makeFilelist(base,flat_dark_range,path=path)
		if dark_range:
			self.obs_dark_names	 = makeFilelist(base,dark_range,path=path)
		self.flat_names = makeFilelist(base,flat_range,path=path)
		self.lamp_names = makeFilelist(base,lamp_range,path=path)
		self.sci_names1 = makeFilelist(base,sci_range1,path=path,incontinuous_fnum=incont_fnum)
		if std_range:
			self.std_names = makeFilelist(base,std_range,path=path)

		print (self.sci_names1)

		self.mode  = 'SciStd'

		# Dict for the data files
		self.tdict = {'science':self.sci_names1} #,'standard':self.std_names}
		# Dict for target names
		self.ndict = {'science':sci_tname}#, 'standard':std_tname}

		# To hold off showing plots in wavelength calibration so that code continues to run
		self.hold_plots = hold_plots

		# Whether to run Level1 steps
		if level1:
			self._level1()

		# Whether to run Level2 steps
		if level2:
			self._level2()

	## Level 1 processing, includes dark subtraction, bad pixel map, flat-fielding ...
	def _level1(self):

		try:
			# Create flat dark
			FDark = Dark(self.flat_dark_names,save=self.save_dark,sci_tname=self.sci_tname,out_dir=self.out_dir,nirspec=self.nirspec)
			# Create flats
			OFlat = Flat(self.flat_names, dark=FDark,save=self.save_flat, SettingsFile=self.SettingsFile,sci_tname=self.sci_tname,out_dir=self.out_dir,nirspec=self.nirspec)		
		except AttributeError:
			# Create flats
			OFlat = Flat(self.flat_names, dark=None,save=self.save_flat, SettingsFile=self.SettingsFile,sci_tname=self.sci_tname,out_dir=self.out_dir,nirspec=self.nirspec)		
		try:
			OLamp = Lamp(self.lamp_names, dark=FDark,save=self.save_flat, SettingsFile=self.SettingsFile,sci_tname=self.sci_tname,out_dir=self.out_dir,nirspec=self.nirspec)
		except AttributeError:
			None
		try:
			# Create observation dark
			ODark = Dark(self.obs_dark_names,sci_tname=self.sci_tname,out_dir=self.out_dir,nirspec=self.nirspec)
		except AttributeError:
			None

		
		img_flat = self.save_flat
		# Initialize dictionary for level1 files
		level1_files = {}
		#print (self.tdict)
		# Iterate over observation types (only science for now)
		for key in self.tdict.keys():
			#print (key)
			try:
				ONod	= Nod(self.tdict[key],flat=OFlat,dark=ODark,lamp=OLamp,tname=self.ndict[key], SettingsFile=self.SettingsFile,
							  plots_dir = self.plots_path, badpix=self.out_dir+'/badpix.dmp',nirspec=self.nirspec)
			except UnboundLocalError:
				ONod	= Nod(self.tdict[key],flat=OFlat,dark=None, tname=self.ndict[key], SettingsFile=self.SettingsFile,
							  plots_dir = self.plots_path, badpix=self.out_dir+'/badpix.dmp',nirspec=self.nirspec)
		  
			matrix = ONod.TargetStack

			# Number of spectral orders
			norders = ONod.getNOrders()
			target_files = []
			# Process each order separately
			for i in np.arange(norders):
				#if i != 2: continue					### if i == 0: continue, skips 0
				#if i == 1: continue
				#if i == 0: continue
				#if i == 3: continue
				
				print ('### Processing order',i+1)
				# Shift A-B images, combine them, and rectify combined A-B image
				OOrder	 = Order(ONod,onum=i+1,write_path=self.spec2d_path, targetType=self.targetType)
				print ('### 2D order extracted')

				OSpec1D	 = Spec1D(OOrder,sa=True,write_path=self.spec1d_path,nirspec=self.nirspec)
				
				# Only need to do WAVECAL for star, not faint companions
				if self.targetType == 'faint':
					OOrder_files = {'2d': OOrder.file, '1d': OSpec1D.file}
					target_files.append(OOrder_files)
				else:
					OWaveCal = WaveCal(OSpec1D.file, path=self.wave_path, am=OSpec1D.airmass, hp=self.hold_plots)
					OOrder_files = {'2d': OOrder.file, '1d': OSpec1D.file, 'wave': OWaveCal.file}
					target_files.append(OOrder_files)

			level1_files[key] = target_files

		filename = self._getLevel1File(self.ndict['science']) 
 
		f = open(self.level1_path+'/'+filename, 'w')
		json.dump(level1_files,f)
		f.close()

	def _level2(self):

		filename = self._getLevel1File(self.ndict['science'])
		f = open(self.level1_path+'/'+filename, 'r')
		level1_files = json.load(f)
		f.close()

		norders = len(level1_files['science'])

		for i in np.arange(norders):
			sci_file = level1_files['science'][i]['wave']
			std_file = level1_files['standard'][i]['wave']
			OCalSpec = CalSpec(sci_file,std_file,shift=self.shift,dtau=self.dtau,write_path=self.cal1d_path,order=i+1)


	def _getLevel1File(self, tname):
		warnings.resetwarnings()
		warnings.filterwarnings('ignore', category=UserWarning, append=True)
		header = pf.open(self.sci_names1[0],ignore_missing_end=True)[0].header
		basename = getBaseName(header, tname)
		filename = basename+'_files.json'
		return filename

def readFilelist(listfile):
	funit = open(listfile)
	flist = funit.readlines()
	return flist

## Create a list of file names for a specific file type (flats, darks etc.)
def makeFilelist(date,ranges,path='', incontinuous_fnum = False):
	# If the list of file numbers are not continuous (e.g.: 193, 194, 197, 198)

	#if not incontinuous_fnum:
	#if incontinuous_fnum:

	if not isinstance(ranges,list):
		ranges = [ranges]

	fnames = []
	for Range in ranges:
		fnumbers = np.arange(Range[0], Range[1]+1, 1)
		for number in fnumbers:
			#print(number)
			#fnames.append(path+date+'s'+str(number).zfill(4)+'.fits')
			fnames.append(path+'nspec'+date+'_'+str(number).zfill(4)+'.fits')

	#else:
	#	fnames = []
	#	fnumbers = ranges[0]
	#	for number in fnumbers:
	#		fnames.append(path + 'nspec' + date + '_' + str(number).zfill(4) + '.fits')
	
	return fnames

def readNIRSPEC(filename):
	header = pf.open(filename)[0].header
	data = pf.getdata(filename)
	return data, header

def getBaseName(header,tname=None):
	time   = header['UTC']
	time   = time.replace(':','')
	time   = time[0:4]
	date   = header['DATE-OBS']
	date   = date.replace('-','')
	if (tname is None):
		object = header['OBJECT']
	else:
		object = tname
	object = object.replace(' ','')
	basename = object+'_'+date+'_'+time
	return basename

def calc_centroid(cc,cwidth=15):
	maxind = np.argmax(cc)
	
	mini = max([0,maxind-cwidth])
	maxi = min([maxind+cwidth,cc.shape[0]])
	trunc = cc[mini:maxi]
	centroid = mini+(trunc*np.arange(trunc.shape[0])).sum()/trunc.sum()
	
	return centroid

def write_fits(array,filename='test.fits'):
	hdu = pf.PrimaryHDU(array)
	hdu.writeto(filename,overwrite=True)
