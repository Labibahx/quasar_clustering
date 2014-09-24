""" Unsupervised clustering analysis on a sample of SDSS quasars.
"""

import numpy as np
from scipy.signal import resample
from astropy.table import Table
from astropy.io import fits
from astropy import units as u
from specutils import extinction
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import linkage, dendrogram

### Read the SDSS DR10 Quasar catalog

data = Table.read('dr10q.fits')

""" extract the part with redshift 1.6 > z > 2.1
also remove the bad measurements and all the other junk (the catalog has EW of -1, -inf or ridiculously large negative numbers. The EW values are >0 for absorption lines so I am only keeping the >0 EWs)
I checked some of the objects with EWs < 0 and there seem to be something wrong with the measurements (the negative EW does not mean absorption line)
I also put an upper cut for the ew <2000 as there seems to be some outliers
"""


sample = data[(data['Z_PCA'] >1.6) & (data['Z_PCA'] <2.1) & (data['REWE_CIII'] >0) & (data['REWE_CIII'] <2000) & (data['REWE_CIV'] >0) & (data['REWE_CIV'] <2000) & (data['REWE_MGII'] >0) & (data['REWE_MGII'] <2000) & (data['REW_CIV'] >0)] #good data with no zeros or inf EWs

### use only some of the parameters to do the clustering

# list of parameters to include in the clustering analysis (features)
z= sample['Z_PCA'] # PCA redshift
c4_fwhm= sample['FWHM_CIV'] # FWHM CIV emission line
c4_bhwhm= sample['BHWHM_CIV'] # blue HWHM of CIV emission line
c4_rhwhm= sample['RHWHM_CIV'] # red HWHM of CIV emission line
c4_amp = sample['AMP_CIV'] # amplitude of CIV emission line (median rms pixel noise)
c4_ew= sample['REWE_CIV'] # EW for the CIV emission line
c3_fwhm= sample['FWHM_CIII']
c3_bhwhm= sample['BHWHM_CIII']
c3_rhwhm= sample['RHWHM_CIII']
c3_amp= sample['AMP_CIII']
c3_ew= sample['REWE_CIII']
mg2_fwhm= sample['FWHM_MGII']
mg2_bhwhm= sample['BHWHM_MGII']
mg2_rhwhm= sample['RHWHM_MGII']
mg2_amp= sample['AMP_MGII']
mg2_ew= sample['REWE_MGII']
bal= sample['BAL_FLAG_VI'] # BAL flag from visual inspection
c4_tew= sample['REW_CIV'] # EW of the CIV trough
sdss_name = sample['SDSS_NAME']
plate= sample['PLATE']
mjd= sample['MJD']
fiber= sample['FIBERID']

# list features to be included in the clustering
features= [z, c4_fwhm, c4_bhwhm, c4_rhwhm, c4_amp, c4_ew, c3_fwhm, c3_bhwhm, c3_rhwhm, c3_amp, c3_ew, mg2_fwhm, mg2_bhwhm, mg2_rhwhm, mg2_amp, mg2_ew, bal, c4_tew]

### combine 1D arrays to a 2D numpy array to perform the clustering on (each row is one quasar, each column is one feature)
qs= np.column_stack(n for n in features)

### generate a list of sdss

### choose k that gives the knee of the sos vs. k plot
num_c= np.arange(1,15) # number of clusters
sos_ls= [] # list of the sum of distances squared

for q in num_c:
    kmeans= KMeans(init= 'k-means++', n_clusters= q, n_init= 10)
    kmeans.fit(qs)
    labels= kmeans.predict(qs)
    sos_ls.append(kmeans.inertia_)

scatter(num_c, sos_ls)
show()

### use k-means for the clustering analysis
nc= 4 #number of clusters
kmeans= KMeans(init= 'k-means++', n_clusters= nc, n_init= 10)
kmeans.fit(qs)
labels= kmeans.predict(qs)

### visualize the clusters -need to do more work on this
scatter(c3_fwhm, c4_fwhm, c=labels)

clusc= kmeans.cluster_centers_  # returrns a list of tuples with the coordinates of the cluster center
sos= kmeans.inertia_ # Sum of distances of samples to their closest cluster center.

### make composite spectra for each cluster
compos= []
for c in range(nc):
    cluster= sample[labels==c]
    clust_spec= np.arange(1100, 4000, 0.1)
    for q in range(len(cluster)):
        z= cluster['Z_PCA'][q]
        mg= cluster['EXTINCTION'][q][1] # using the extinction in g magnitude
        name='./data/spec-'+str(cluster['PLATE'][q])+'-'+str(cluster['MJD'][q])+'-'+str(cluster['FIBERID'][q]).zfill(4)+'.fits'
        obj=fits.open(name)
        flx= obj[1].data.field(0)
        wlen= 10.**(obj[0].header['coeff0']+ obj[0].header['coeff1'] * np.arange(len(obj[1].data.field(1)))) # coeff0: starting wavelength, coeff1: dispersion
        ext_lambda= extinction.extinction_ccm89(wlen * u.angstrom, a_v= mg, r_v=3.1)
        tau_lambda= ext_lambda/((1 * u.angstrom) * 1.086)  # for some reason the extinction has a unit of angstrom (should be magnitude i.e. unitless)
        dered_flx= flx * np.exp(tau_lambda) # extinction corrected flux
        rest_wlen= wlen / (1+ z) # restframe wavelength
        x= np.arange(1100, 4000, 0.1) # the wavelength values for the rebinned spectrum with bin width = 0.1 A
        rebin_dered_flx= np.interp(x, rest_wlen, dered_flx) # resampled flux with delta lambda= 0.1 A
##        rebin_dered_flx= resample(dered_flx, 4550)  # this does not seem to be correct
        norm_flx= rebin_dered_flx/mean(rebin_dered_flx[1810:1820]) # normalize spectra. need to check if this part of spec has no lines
        clust_spec= np.vstack((clust_spec, norm_flx)) # 2D array. 1st row: restframe wavelength, other rows have fluxes of spectra from single cluster
    compos.append(np.median(clust_spec, axis=0)) # list with the composites (compos[0] is composite from 1st cluster, compos[1] 2nd cluster,...)

## plot composites
fig= figure()
plot(x, compos[0])
plot(x, compos[1]+1)
plot(x, compos[2]+2)
plot(x, compos[3]+3)
axvline(1549, ls= '--')
text(1549, 5.5, r'C IV')
axvline(1908, ls= '--')
text(1908, 5.5, r'C III]')
axvline(2798, ls= '--')
text(2798, 5.5, r'Mg II')
xlabel(r'Restframe wavelength ($\AA$)')
ylabel(r'Normalized flux')






