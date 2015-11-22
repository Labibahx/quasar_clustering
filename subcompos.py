"""21 Nov 2015
    read FITS files of a given list of spectra and create a median composite spectrum from them
    """

import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.stats import sigmaclip

sample= "Main Sample"
sample_name= "main"
line= "c3"
line_name= "C III]"
k= 5
comments= "High Mi"

clstr= Table.read("./clusters/"+line+"_"+str(k)+"clstrs_"+sample_name+".fits")
#clstr= Table.read("./clusters/"+line+"_"+str(k)+"clstrs_"+sample_name+".fits")
#clstr= Table.read("./clusters/"+line+"_"+str(k)+"clstrs_"+sample_name+".fits")

data= Table.read("sample_myflags.fits")
#data= Table.read("sample_mixed_myflags.fits")
#data= Table.read("sample_bal_myflags.fits")

tt= join(clstr, data, keys= "SDSS_NAME")

#t= tt[(tt['label'] ==2) & (tt['MI'] > -25.5)]

for c in range(k):
    #t= tt[(tt['MI'] < -26.5) & (tt['label'] == c)] # low Mi
    t= tt[(tt['MI'] > -25.5) & (tt['label'] == c)] # high Mi

    compo_array= np.arange(1100, 4000, 0.5)

    for i in range(len(t)):
        spec_name="./new_proc_data/spec-"+str(t['PLATE'][i])+"-"+str(t['MJD'][i])+"-"+str(t['FIBERID'][i]).zfill(4)+"_proc.fits"
        spec=fits.open(spec_name)
        flx= spec[0].data[1]
        wlen= spec[0].data[0]
        norm_flx= flx/np.median(flx[2360:2390]) # normalize spectra
        compo_array= np.vstack((compo_array, norm_flx)) # 2D array. 1st row: restframe wavelength, other rows have corrected fluxes of spectra from clusters (one for each row)
        del spec

    n= len(compo_array[1:,])
    print "composite has", n, "spectra"

    clipped_compo=[]
    for j in range(compo_array.shape[1]):
        
        y= sigmaclip(compo_array[1:,j], 3, 3)
        m=median(y[0])
        clipped_compo.append(m)

    compo_name= "./composites/hiMI_"+line+"_"+str(c)+sample_name+".fits"
    spec_file= np.vstack((wlen,clipped_compo))
    hdu= fits.PrimaryHDU(spec_file)
    hdr= hdu.header
    hdr.set('SPEC_NUMBER', n)
    hdr.set('COMPOSITE', comments)
    hdu.writeto(compo_name)


