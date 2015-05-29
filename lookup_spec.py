

import numpy as np
from astropy.io import fits

from astropy.table import Table



def spec_look_up(cluster_array, k):

    """ read a list of sdss names, corss-match list with table with mjd, plate, fiber IDs to generate spectra file name with format mjd-plate-fiber_proc.fits e.g., spec-3587-55182-0691_proc.fits
    the processed spectra (i.e., corrected for Galactic extinction and de-redshifted and normalized. See spec_proc.py)
    
    cluster_array: a 2D numpy array with the results of a clustering trial.  Each sample (row) has the values for the features (parameters) used in the clustering and the sdss names of the objects in each cluster. and a clolumn withe clusters labels
    k: cluster label. k=0 --> first cluster, k=1 --> second cluster... 
    
        """

    all_clstrs= np.load(cluster_array)
    
    data= Table.read('dr10q.fits') #DR10 catalog
    
    ss = data[(data['Z_PCA'] >1.6) & (data['Z_PCA'] <2.1) & (data['REWE_CIII'] >0) & (data['REWE_CIII'] <2000) & (data['REWE_CIV'] >0) & (data['REWE_CIV'] <2000) & (data['REWE_MGII'] >0) & (data['REWE_MGII'] <2000) & (data['BAL_FLAG_VI'] ==0)]
    
    #corss-match the above two files
    
    clstr_k= ss[(ss['SDSS_NAME'] == all_clstrs[:,4]) & (all_clstrs[:, 3].astype(int) == k)] # only samples in cluster k


   ## c1= ss[(ss['SDSS_NAME']== clstr_array[:,4]) & (clstr_array[:,3].astype(int) ==0)]
    

    print len(clstr_k)

    spec_files_list=[]
    for s in range(len(clstr)):
        spectrum_name= "spec-"+str(clstr['PLATE'][s])+"-"+str(clstr['MJD'][s])+"-"+str(clstr['FIBERID'][s]).zfill(4)+"_proc.fits"
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



