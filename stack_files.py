
from site import USER_BASE
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u
from PyAstronomy import pyasl
from spec_utils.spec_utils import vgeo

# -------------------------------------------------------------------------------------------
# Stack and shift spectra to heliocentric ref frame
# -------------------------------------------------------------------------------------------

class SpecFinal():
    def __init__(self,files,path=None,outfile=None,vr=0,del_v=0):
        '''
        Class to stack and shift two files
        or if only one item in files, shifts wavelengths to heliocentric ref frame

        Parameters
        ----------
        files: list of paths to two calspec files
        path: (str, opt)
            path to save new file in (CAL1D recommended)
        outfile: (str, opt)
            name of outfile (ex: 'FZTau_COice_2.fits')

        '''
        self.path=path
        self.outfile=outfile
            
        wave0,flux0,uflux0,sa0,usa0 = self.readfile(files[0])
        try:
            print('Stacking')
            wave1,flux1,uflux1,sa1,usa1 = self.readfile(files[1])

            self.wave0,self.flux0,self.sa0 = wave0,flux0,sa0
            self.uflux0,self.usa0 = uflux0,usa0
            self.flux1 = np.interp(wave0,wave1,flux1)
            self.uflux1 = np.interp(wave0,wave1,uflux1)
            self.sa1 = np.interp(wave0,wave1,sa1)
            self.usa1 = np.interp(wave0,wave1,usa1)

            # average spectra 
            fs = [flux0,flux1]
            self.flux = np.array(fs).mean(axis=0)
            sas = [sa0,sa1]
            self.sa = np.array(sas).mean(axis=0)

            # errors
            self.uflux = np.sqrt((self.uflux0/self.flux0)**2+(self.uflux1/self.flux1)**2)
            self.usa = np.sqrt((self.usa0/self.sa0)**2+(self.usa1/self.sa1)**2)

        except:
            self.wave0=wave0
            self.flux=flux0
            self.uflux=uflux0
            self.sa=sa0
            self.usa=usa0


        # doppler shift -> heliocentric ref frame
        if del_v == 0 and vr != None:
            del_v = CalShift(files[0],vr=vr).del_v
        self.wave = self.shift_wave(self.wave0,del_v)
        filename=self.writeSpec(outfile=self.outfile,path=self.path)
        #self.plotavg()
        self.plotshift()

    def readfile(self,file):
        # read wave, flux, sa data
        file = fits.open(file)
        data = file[1].data
        self.header = file[0].header

        wave = data['wave']
        flux = data['flux']
        uflux = data['uflux']
        sa = data['sa']
        usa = data['usa']
        return wave,flux,uflux,sa,usa

    def plotavg(self):
        fig,axs = plt.subplots(2)
        axs[0].plot(self.wave0,self.flux0,label='flux0')
        axs[0].plot(self.wave0,self.flux1,label='flux1')
        axs[0].plot(self.wave0,self.flux,label='flux')
        axs[0].legend()

        axs[1].plot(self.wave0,self.sa0,label='sa0')
        axs[1].plot(self.wave0,self.sa1,label='sa1')
        axs[1].plot(self.wave0,self.sa,label='sa')
        axs[1].legend()
        plt.show()

    def plotshift(self):
        fig,axs = plt.subplots(2)
        axs[0].plot(self.wave0,self.flux,label='telluric frame')
        axs[0].plot(self.wave,self.flux,label='heliocentric frame')
        axs[0].legend()

        axs[1].plot(self.wave0,self.sa,label='telluric frame')
        axs[1].plot(self.wave,self.sa,label='heliocentric frame')
        axs[1].legend()
        plt.show()

    def writeSpec(self, outfile=None, path='.'):
        c1	= fits.Column(name='wave', format='D', array=self.wave)
        c2	= fits.Column(name='flux', format='D', array=self.flux)
        c3	= fits.Column(name='uflux', format='D', array=self.uflux)
        c4	= fits.Column(name='sa', format='D', array=self.sa)
        c5	= fits.Column(name='usa', format='D', array=self.usa)
        coldefs = fits.ColDefs([c1,c2,c3,c4,c5])

        tbhdu = fits.BinTableHDU.from_columns(coldefs)

        hdu = fits.PrimaryHDU(header=self.header)
        thdulist = fits.HDUList([hdu,tbhdu])

        if outfile is None:
            filepath = self.getBaseName(self.header)
        else:
            filepath = self.path+outfile

        thdulist.writeto(filepath,overwrite=True)
        print('writing...',filepath)
        return filepath
    
    def getBaseName(self,header,tname=None):
        try:
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
            basename = object+'_'+date+'_'+time+'.fits'
        except:
            print('Error: file path')
        filepath = self.path+basename
        return filepath

    def shift_wave(self,waves,del_v):
        # shift to heliocentric reference frame
        wl = waves*1e-6 #[m] wavelength array to shift
        c = 3e8 #m/s
        v = del_v*1000 #m/s
        del_wl = (v/c)*wl #[m] shifts for each wavelength
        wl -= del_wl
        return wl*1e6 #microns

# -------------------------------------------------------------------------------------------
# Calculate Doppler Shift
# -------------------------------------------------------------------------------------------


class CalShift():
    def __init__(self,files,vr):
        '''
        Calculate heliocentric velocity
        files: list of paths to files
        '''
        # print outputs read from header
        hd     = ['File','Target Name','JD','UTC','RA (deg)','Dec (deg)','Bary_v (km/s)',
            'Doppler Shift (km/s)']
        string = '{:<5}{:<14}{:<16}{:<14}{:<12}{:<12}{:<15}{:<20}'
        hd_str = '{:<5}{:<14}{:<16}{:<14}{:<12}{:<12}{:<15}{:<20}'
        space  = '-' *110
        print(space + "\n" + hd_str.format(*hd) + "\n" + space)

        file1 = files
        entry = fits.getheader(file1, 0)
        jd,ra,dec,vb = self.calc_vb(file1)
        
        # Doppler shift with star's radial velocity - barycentric
        ds = self.calc_ds(vr,vb) #km/s

        mydate,mycoord,entry = self.readheader(files)

        # Print outputs
        var = [entry['FRAMENUM'], entry['OBJECT'], str(jd), entry['UTC'], 
        str(round(ra,6)), str(round(dec,6)),str(round(vb,8)),str(round(ds,8))]
        print(string.format(*var))

        # Calculate heliocentric velocity (Earth-induced + intrinsic)
        # + redshift, - blueshift
        earthv=vgeo(mydate, mycoord, vhel=0)
        myv=vgeo(mydate, mycoord, vhel=vr)
        print('Earth-induced shift',earthv)
        print('heliocentric velocity',myv)

        #return myv #doppler shift
        self.del_v = myv


    def readheader(self,filename):
        '''
        Get observation time and coordinates from header
        filename: (str) full path to file
        '''
        # read date, RA, Dec 
        entry = fits.getheader(filename,0)
        mydate = Time(entry['MJD'],format='mjd')

        ra = entry['RA']
        dec = entry['DEC']
        mycoord = SkyCoord(ra+dec, unit=(u.hourangle, u.deg))
        ra = mycoord.ra.degree
        dec = mycoord.dec.degree
        return mydate,mycoord,entry
        

        vh, vb = pyasl.baryCorr(jd, ra, dec, deq=2000.0)
        #print("Barycentric velocity of Earth toward star: ", vb,"km/s")
        return jd,ra,dec,vb

    
    def calc_ds(self,v_systemic,v_barycentric):
        # calculate doppler shift with respect to Earth
        v_primary = v_systemic - v_barycentric
        return v_primary

    def calc_vb(self,filename):
        entry = fits.getheader(filename,0)
        jd = float(entry['MJD'])+2400000.5
        ra = entry['RA']
        dec = entry['DEC']
        c = SkyCoord(ra+dec, unit=(u.hourangle, u.deg))
        ra = c.ra.degree
        dec = c.dec.degree



        vh, vb = pyasl.baryCorr(jd, ra, dec, deq=2000.0)
        #print("Barycentric velocity of Earth toward star: ", vb,"km/s")
        return jd,ra,dec,vb

    def calc_wave(self,del_v):
        wl = 5e-06 #m
        c = 3e8 #m/s
        v = del_v*1000 #m/s
        del_wl = (v/c)*wl #m
        return del_wl*(10**(-6)) #microns


#mydate=Time('2014-08-16T00:00:00.0', format='isot', scale='utc')
#mycoord=SkyCoord('16h31m33.46s', '-24d27m37.3s', frame='icrs')
#Calculate the heliocentric velocity (Earth-induced + intrinsic)
#myv=vgeo(mydate, mycoord, vhel=-7.93)

#(set vhel=0 if you just want the Earth-induced shift)


# -------------------------------------------------------------------------------------------
# Combine list of files
# -------------------------------------------------------------------------------------------

class Combine():
    def __init__(self,files,outfile=None,outpath=None):
        waves,fluxes,ufluxes,sas,usas = [],[],[],[],[]
        for f in files:
            wave,flux,uflux,sa,usa = self.readfile(f)
            waves = np.append(waves,wave)
            fluxes = np.append(fluxes,flux)
            ufluxes = np.append(ufluxes,uflux)
            sas = np.append(sas,sa)
            usas = np.append(usas,usa)
        self.wave = waves
        self.flux = fluxes
        self.uflux = ufluxes
        self.sa = sas
        self.usa = usas
        file = fits.open(files[0])
        self.header = file[0].header

        filepath = self.writeSpec(outfile=outfile,path=outpath)


    def readfile(self,file):
        # read wave, flux, sa data
        file = fits.open(file)
        data = file[1].data
        self.header = file[0].header

        wave = data['wave']
        flux = data['flux']
        uflux = data['uflux']
        sa = data['sa']
        usa = data['usa']
        return wave,flux,uflux,sa,usa
    
    def writeSpec(self, outfile=None, path='.'):
        c1	= fits.Column(name='wave', format='D', array=self.wave)
        c2	= fits.Column(name='flux', format='D', array=self.flux)
        c3	= fits.Column(name='uflux', format='D', array=self.uflux)
        c4	= fits.Column(name='sa', format='D', array=self.sa)
        c5	= fits.Column(name='usa', format='D', array=self.usa)
        coldefs = fits.ColDefs([c1,c2,c3,c4,c5])

        tbhdu = fits.BinTableHDU.from_columns(coldefs)

        hdu = fits.PrimaryHDU(header=self.header)
        thdulist = fits.HDUList([hdu,tbhdu])

        filepath = path+outfile

        thdulist.writeto(filepath,overwrite=True)
        print('writing...',filepath)
        return filepath
    