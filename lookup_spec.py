

import numpy as np
from astropy.io import fits

from astropy.table import Table



def spec_look_up(cluster_array, k):

    """ read a list of sdss names, corss-match list with table with mjd, plate, fiber IDs to generate spectra file name with format mjd-plate-fiber_proc.fits e.g., spec-3587-55182-0691_proc.fits
    the processed spectra (i.e., corrected for Galactic extinction and de-redshifted and normalized. See spec_proc.py)
    
    cluster_array: a 2D numpy array with the results of a clustering trial.  Each sample (row) has the values for the features (parameters) used in the clustering and the sdss names of the objects in each cluster. and a clolumn withe clusters labels
    k: cluster label. k=0 --> first cluster, k=1 --> second cluster... 
    
        """

    clstr= np.load(cluster_array)

    sdss_names= clstr[:,4][clstr[:,3].astype(int)== k]
    print sdss_names.shape

    data= Table.read('dr10q.fits') #DR10 catalog
    
    
   ## c1= ss[(ss['SDSS_NAME']== clstr_array[:,4]) & (clstr_array[:,3].astype(int) ==0)]
    
    mjd= data['MJD'][(data['SDSS_NAME'] == clstr[:,4]) & (clstr[:,3].astype(int) == 3)]
    #plate= data['PLATE'][data['SDSS_NAME'] == sdss_names]
    #fiber= data['FIBERID'][data['SDSS_NAME'] == sdss_names]

    print mjd.shape, sdss_names.shape

    spec_files_list=[]
    for s in range(len(sdss_names)):
        spectrum_name= "spec-"+str(plate[s])+"-"+str(mjd[s])+"-"+str(fiber[s]).zfill(4)+"_proc.fits"
        spec_files_list.append(spectrum_name)

    return(spec_files_list)




def spec_display(spec_list):

""" look up and display a list quasars. Read file names from a list
    """

    spectra= spec_list
    
    wavelen= np.arange(1100, 4000, 0.1)  #wavelength array

    for file in spectra:
        try:
            spec= fits.open(file)
            flx= spec[0].data
            plot(wavelen, flx)
            resume = input("Press Enter to plot next spectrum on list.")
        except SyntaxError:
            pass



