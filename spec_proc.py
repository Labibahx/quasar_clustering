""" Read a raw spectrum from the SDSS QSO DR10 catalog.
    correct for Galactic extinction and redshift and save as a new fits file
    """

import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy import units as u
from specutils import extinction

data= Table.read('dr10q.fits')

s= Table.read('sample_myflags.csv')

#s = data[(data['Z_PCA'] >1.6) & (data['Z_PCA'] <2.1) & (data['REWE_CIII'] >0) & (data['REWE_CIII'] <2000) & (data['REWE_CIV'] >0) & (data['REWE_CIV'] <2000) & (data['REWE_MGII'] >0) & (data['REWE_MGII'] <2000) & (data['BAL_FLAG_VI'] ==0)] #  subsample with: upper and lower redshift limits, reasonable measurements for EW, and BAL quasars excluded

for i in range(len(s)):
    print i
    z= s['Z_PCA'][i]
    mg= 3.1 * s['EXTINCTION'][i][1] #use the extinction in the g magnitude E(B-V)
    spec_name= './data/spec-'+str(s['PLATE'][i])+'-'+str(s['MJD'][i])+'-'+str(s['FIBERID'][i]).zfill(4)+'.fits' # fits file name
    spec= fits.open(spec_name) # read file
    flx= spec[1].data.field(0) # flux
    wlen= 10.**(spec[0].header['coeff0']+ spec[0].header['coeff1'] * np.arange(len(spec[1].data.field(1)))) # wavelength, coeff0: starting wavelength, coeff1: dispersion
    ext_lambda= extinction.extinction_ccm89(wlen * u.angstrom, a_v= mg, r_v= 3.1)
    tau_lambda= ext_lambda/((1 * u.angstrom) * 1.086) # for some reason the extinction has a unit of angstrom (should be magnitude i.e. unitless)
    dered_flx= flx * np.exp(tau_lambda) # extinction corrected flux
    rest_wlen= wlen /(1+z) # restframe wavelength
    x= np.arange(1100, 4000, 0.5) # the wavelength values for the rebinned spectrum with bin width = 0.5 A
    rebin_dered_flx= np.interp(x, rest_wlen, dered_flx) # resampled flux with delta lambda= 0.5 A
   # norm_flx= rebin_dered_flx/np.mean(rebin_dered_flx[2360:2390]) # normalize spectra
    proc_spec= np.vstack((x, rebin_dered_flx)) # 2D array. 1st row: restframe wavelength, 2nd row: corrected flux
    hdu= fits.PrimaryHDU(proc_spec) # HDU with flux and wavelngth to be saved as a fits file
    new_spec_name= './new_proc_data/spec-'+str(s['PLATE'][i])+'-'+str(s['MJD'][i])+'-'+str(s['FIBERID'][i]).zfill(4)+'_proc.fits' # name of new processed spectrum
    hdu.writeto(new_spec_name)
