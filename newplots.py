""" plots similar to what is in lstr_plots.py but ordered according to the number of objects in each cluster 
May 8 2015"""

import numpy as np
from astropy.io import fits
from astropy.table import Table
from sklearn.cluster import KMeans
from mpl_toolkits.mplot3d import Axes3D
from operator import itemgetter
from glob import glob


## the following is copied from quasar_cluster.py

data= Table.read('dr10q.fits')
ss = data[(data['Z_PCA'] >1.6) & (data['Z_PCA'] <2.1) & (data['REWE_CIII'] >0) & (data['REWE_CIII'] <2000) & (data['REWE_CIV'] >0) & (data['REWE_CIV'] <2000) & (data['REWE_MGII'] >0) & (data['REWE_MGII'] <2000) & (data['BAL_FLAG_VI'] ==0)] # subsample with: upper and lower redshift limits, reasonable measurements for EW, and BAL quasars excluded

### use only some of the parameters to do the clustering

# list of parameters to include in the clustering analysis (features)
redshift= ss['Z_PCA'] # PCA redshift
c4_fwhm= ss['FWHM_CIV'] # FWHM CIV emission line
c4_bhwhm= ss['BHWHM_CIV'] # blue HWHM of CIV emission line
c4_rhwhm= ss['RHWHM_CIV'] # red HWHM of CIV emission line
c4_amp = ss['AMP_CIV'] # amplitude of CIV emission line (median rms pixel noise)
c4_ew= ss['REWE_CIV'] # EW for the CIV emission line
c3_fwhm= ss['FWHM_CIII']
c3_bhwhm= ss['BHWHM_CIII']
c3_rhwhm= ss['RHWHM_CIII']
c3_amp= ss['AMP_CIII']
c3_ew= ss['REWE_CIII']
mg2_fwhm= ss['FWHM_MGII']
mg2_bhwhm= ss['BHWHM_MGII']
mg2_rhwhm= ss['RHWHM_MGII']
mg2_amp= ss['AMP_MGII']
mg2_ew= ss['REWE_MGII']
bal= ss['BAL_FLAG_VI'] # BAL flag from visual inspection
c4_tew= ss['REW_CIV'] # EW of the CIV trough
alpha_nu= ss['ALPHA_NU'] # spectra index between 1450−1500 A, 1700−1850 A and 1950−2750 A
sdss_name = ss['SDSS_NAME']
plate= ss['PLATE']
mjd= ss['MJD']
fiber= ss['FIBERID']

features= [[c4_ew, c4_bhwhm, c4_rhwhm, c4_amp, c4_fwhm],
           [c3_ew, c3_bhwhm, c3_rhwhm, c3_amp, c3_fwhm],
           [mg2_ew, mg2_bhwhm, mg2_rhwhm, mg2_amp, mg2_fwhm]]


## each cluster is saved into a numpy 2D array which includes the clustering parameters and the name of the sdss object. To get the other parameters for the same object that were not used in the clustering, I can cross match the cluster 2D array with the full subsample array ss (defined in line 14 here).



## each cluster is represented by a composite spectrum with the number of objects in each cluster given in the FITS header (keyword:'SPEC_NUMBER')

clstr_name= []
num_ls=[]
for c in clstr_name:
    c_array= np.load()
    num_ls.append(len(c_array))

for i,j in zip(clstr_name, num_ls):
    clstr_ls.append((i,j))

c1= fits.open()


line_clstrs= [[("mg2", 3), ("mg2", 4)], [("c3", 5), ("c3", 6)], [("c4", 3), ("c4", 4)]]

for i in range(3):
    clstr1= np.load(line_clstrs[0][i]+"_ew_hwhm_"+str(line_clstrs[1][i])+"clstrs.npy")
    clstr1= np.load(l_n[0]+"_ew_hwhm_"+str(l_n[1])+"clstrs.npy")

def profiles(line, k):
    """ plot line profiles for the clusters in 4 panels
        """
    compos= glob(line+"ew_hwhm_"+str(k)+"*.fits")
    
    
    
    
    clstrs= glob.glob(line+"_ew_hwhm*clstrs.npy")
    clstrs1= np.load(clstrs[0])
    clstrs2= np.load(clstrs[1])
    fig= figure()
    ax1= fig.add_sublot(241)
    for c in range(1,k1+1):
        ax1.plot()
    
    como_list.append(line)
    orderd_compos= sorted()
    for f in compos:
        spec_file= fits.open(f)

sorted(list, key=itemgetter(1))










