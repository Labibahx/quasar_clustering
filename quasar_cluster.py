""" Aesha Tammour 7 August 2014
Unsupervised clustering analysis on a sample of SDSS quasars.
Basically, imitating Mark's code!

"""

import numpy as np
from astropy.table import Table

from sklearn.cluster import KMeans

from scipy.cluster.hierarchy import linkage, dendrogram

### Read the SDSS DR10 Quasar catalog

all_data = Table.read('dr10q.fits')

""" extract the part with redshift 1.6 > z > 2.1
also remove the bad measurements and all the other junk (the catalog has EW of -1, -inf or ridiculously large negative numbers. The EW values are >0 for absorption lines so I am only keeping the >0 EWs)
I checked some of the objects with EWs < 0 and there seem to be something wrong with the measurements (the negative EW does not mean absorption line)
I also put an upper cut for the ew <2000 as there seems to be some outliers
"""


data = all_data[(all_data['Z_PCA'] >1.6) & (all_data['Z_PCA'] <2.1) & (all_data['REWE_CIII'] >0) & (all_data['REWE_CIII'] <2000) & (all_data['REWE_CIV'] >0) & (all_data['REWE_CIV'] <2000) & (all_data['REWE_MGII'] >0) & (all_data['REWE_MGII'] <2000) & (all_data['REW_CIV'] >0)] #good data

### use only some of the parameters to do the clustering

# list of parameters to include in the clustering analysis (features)
z= data['Z_PCA'] # PCA redshift
c4_fwhm= data['FWHM_CIV'] # FWHM CIV emission line
c4_bhwhm= data['BHWHM_CIV'] # blue HWHM of CIV emission line
c4_rhwhm= data['RHWHM_CIV'] # red HWHM of CIV emission line
c4_amp = data['AMP_CIV'] # amplitude of CIV emission line (median rms pixel noise)
c4_ew= data['REWE_CIV'] # EW for the CIV emission line
c3_fwhm= data['FWHM_CIII']
c3_bhwhm= data['BHWHM_CIII']
c3_rhwhm= data['RHWHM_CIII']
c3_amp= data['AMP_CIII']
c3_ew= data['REWE_CIII']
mg2_fwhm= data['FWHM_MGII']
mg2_bhwhm= data['BHWHM_MGII']
mg2_rhwhm= data['RHWHM_MGII']
mg2_amp= data['AMP_MGII']
mg2_ew= data['REWE_MGII']
bal= data['BAL_FLAG_VI'] # BAL flag from visual inspection
c4_tew= data['REW_CIV'] # EW of the CIV trough


#cols= ['Z_PCA', 'FWHM_CIV', 'BHWHM_CIV','RHWHM_CIV', 'AMP_CIV', 'REWE_CIV', 'FWHM_CIII', 'BHWHM_CIII', 'RHWHM_CIII', 'AMP_CIII', 'REWE_CIII', 'FWHM_MGII', 'BHWHM_MGII', 'RHWHM_MGII', 'AMP_MGII', 'REWE_MGII', 'BAL_FLAG_VI']

# list of names for the arrays
features= [z, c4_fwhm, c4_bhwhm, c4_rhwhm, c4_amp, c4_ew, c3_fwhm, c3_bhwhm, c3_rhwhm, c3_amp, c3_ew, mg2_fwhm, mg2_bhwhm, mg2_rhwhm, mg2_amp, mg2_ew, bal, c4_tew]

### combine 1D arrays to a 2D numpy array to perform the clustering on
#qs= np.column_stack((z, c4_fwhm, c4_bfwhm, c4_rfwhm, c4_amp, c4_ew, c3_fwhm, c3_bfwhm, c3_rfwhm, c3_amp, c3_ew, mg2_fwhm, mg2_bfwhm, mg2_rfwhm, mg2_amp, mg2_ew))

qs= np.column_stack(n for n in features)

kmeans= KMeans(init= 'k-means++', n_clusters= 3, n_init= 10)
kmeans.fit(qs)
labels= kmeans.predict(qs)

### visualize the results

clusc= kmeans.cluster_centers_  # returrns a list of tuples with the coordinates of the cluster center
sos= kmeans.inertia_ # Sum of distances of samples to their closest cluster center.

### choose k that gives the knee of the sos vs. k plot
num_c= np.arange(1,15) # number of clusters
sos_ls= [] # list of the sum of distances

for q in num_c:
    kmeans= KMeans(init= 'k-means++', n_clusters= q, n_init= 10)
    kmeans.fit(qs)
    labels= kmeans.predict(qs)
    sos_ls.append(kmeans.inertia_)

scatter(num_c, sos_ls)
show()



