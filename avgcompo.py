"""13 Aug 2015 
read names of fits files from fits files and create a mean composite spectrum from them
"""

import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.stats import sigmaclip

def avg_compo(sample):

    """param:
    table: table file to read data from to create composite
    sample: name of the sample: main, mixed, bal"""
    
    if sample == "main":
        sample_name= ""
        comments= "Main sample, nonBALs with only EW >0"
    
    elif sample == "mixed":
        sample_name= "_mixed"
        comments= "Mixed sample, BAL and nonBAL"
    
    elif sample == "bal":
        sample_name= "_bal"
        comments= "BAL quasars only"

    table_name= "sample"+sample_name+"_myflags.fits"
    tab= Table.read(table_name)
    t= tab[tab['MY_FLAG'] ==0 ]
   
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

    avg_compo_name= "./composites/mean_compo"+sample_name+".fits"
    spec_file= np.vstack((wlen,clipped_compo))
    hdu= fits.PrimaryHDU(spec_file)
    hdr= hdu.header
    hdr.set('SPEC_NUMBER', n)
    hdr.set('COMPOSITE', comments)
    hdu.writeto(avg_compo_name)
