""" Unsupervised clustering analysis on a sample of SDSS quasars.
link to the Paris et al. 2014 paper (SDSS-DR 10)
http://adslabs.org/adsabs/abs/2014A%26A...563A..54P/
This script is similar to spec_cluster.py but it reads corrected spectra (for extiction and redshift) then performs the clustering analysis and creates composites.
spec_cluster.py includes correcting the spectra, clustering and stacking.
I separated these processes to see if this is what is causing the script to be slow -turns out it is just the same.

The working directory for this script is the main quasar_clustering directory. The catalgue file is (dr10q.fits) is there. The script reads the proccessed fits spectra from the proc_data directory.

====================
Updated on 30 May 2015 to include a limit on the SNR in the selctions. This cut-off brought the size of the sample from 7754 to 4342 quasars.

====================
6 August 2015
This is a copy of quasar_cluster_balq.py to use for clustering on BAL quasars only sample

"""

import numpy as np
from astropy.table import Table
from astropy.io import fits
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.stats import sigmaclip
from sklearn.decomposition import PCA
from sklearn import metrics
import time
from mpl_toolkits.mplot3d import Axes3D

### Read the SDSS DR10 Quasar catalog

"""
data = Table.read('dr10q.fits')

##extract the part with redshift 1.6 > z > 2.1 with some cut-off on uncertainty and S/N

ss = data[(data['Z_PCA'] >1.6) & (data['Z_PCA'] <2.1)
          & (data['ERR_REWE_CIII'] < data['REWE_CIII']/10)
          & (data['ERR_REWE_CIV'] < data['REWE_CIV']/10)
          & (data['ERR_REWE_MGII'] < data['REWE_MGII']/10)
          & (data['SNR_1700'] > 3)]
          
          This is how the BAL+nonBAL sample was selected. File saved as dr10qsample_BAL_only.fits then I flagged a few objects with issues in their spectra. File saved as sample_BAL_myflags.fits
"""

data= Table.read('sample_mixed_myflags.fits')

tt= data[(data['BAL_FLAG_VI'] ==1) & (data['MY_FLAG']== 0)]

tt.write('sample_bal_myflags.fits')


# list of parameters to include in the clustering analysis (features)
redshift= tt['Z_PCA'] # PCA redshift
c4_fwhm= tt['FWHM_CIV'] # FWHM CIV emission line
c4_bhwhm= tt['BHWHM_CIV'] # blue HWHM of CIV emission line
c4_rhwhm= tt['RHWHM_CIV'] # red HWHM of CIV emission line
c4_amp = tt['AMP_CIV'] # amplitude of CIV emission line (median rms pixel noise)
c4_ew= tt['REWE_CIV'] # EW for the CIV emission line
c3_fwhm= tt['FWHM_CIII']
c3_bhwhm= tt['BHWHM_CIII']
c3_rhwhm= tt['RHWHM_CIII']
c3_amp= tt['AMP_CIII']
c3_ew= tt['REWE_CIII']
mg2_fwhm= tt['FWHM_MGII']
mg2_bhwhm= tt['BHWHM_MGII']
mg2_rhwhm= tt['RHWHM_MGII']
mg2_amp= tt['AMP_MGII']
mg2_ew= tt['REWE_MGII']
bal= tt['BAL_FLAG_VI'] # BAL flag from visual inspection
c4_tew= tt['REW_CIV'] # EW of the CIV trough
alpha_nu= tt['ALPHA_NU'] # spectra index between 1450−1500 A, 1700−1850 A and 1950−2750 A
sdss_name = tt['SDSS_NAME']
plate= tt['PLATE']
mjd= tt['MJD']
fiber= tt['FIBERID']

###
### a list of lists of the features entering the clustering analysis
###
features= [c3_ew, c3_bhwhm, c3_rhwhm]

### combine 1D arrays to create a 2D numpy array to perform the clustering analysis on (each row is one quasar, each column is one feature)

qs= np.column_stack(n for n in features) # one line only. 0 can be changed to use a different line

####
# use the following two lines to do dimensionality reduction on the features matrix using PCA
####

### apply PCA to the features matrix
pca = PCA(n_components=3, whiten=True).fit(features) ## Randomly chose 3 PC here. Number of PC can be changed.

### the new qs is now given by the pca.components_ attribute
pca_qs= pca.components_


'''choose k that gives the knee of the sos vs. k plot
    Two ways to do that:
    1- plot cost function (within cluster sum of squares) vs. number of clusters and look for elbow.
    2- silhouette score
    '''

num_c= np.arange(2,9) # number of clusters
sos_ls= [] # list of the sum of distances squared
sil_score_c3= []


for q in num_c:
    kmeans= KMeans(init= 'k-means++', n_clusters= q, n_init= 10)
    kmeans.fit(qs)
    labels= kmeans.predict(qs)
    sc= metrics.silhouette_score(qs, labels)
    print q, sc
    sos_ls.append(kmeans.inertia_)
    sil_score_c3.append(sc)


fig= figure(figsize=(8,6))
plot(num_c, sil_score_c3, marker= 'v', color= '0.3', ls='-.', label= 'C III]')


text(3.5,0.85, "Features: EW, RHWHM, BHWHM")

legend(numpoints= 1)
ylabel('Silhouette score')
xlabel(r'$K$')
xlim(1.9, 8.1)
ylim(0.3, 0.72)


"""plot sos vs K CIII]
"""

fig= figure(figsize=(8,6))

plot(num_c, sos_ls, marker= 'v', color= '0.3', ls='-.', label= 'C III]')

ylabel('Sum of squares')
xlabel(r'$K$')
xlim(1.9, 8.1)
legend(numpoints=1)


### test reproducibility -cluster centroids are the same for several runs of KMeans

lines = [('CIII', 3), ('CIII', 4), ('CIII', 5), ('CIII', 6)]

param_list = ['REWE_', 'BHWHM_', 'RHWHM_']


for l in lines:
    cntrs = open(l[0]+str(l[1])+"_bal_only.txt", 'wr')

    print l[0], ",K=", l[1]
    
    cntrs.write("#"+str(l)+'\n')
    qs= np.column_stack(tt[p+l[0]] for p in param_list)

    k=l[1] #number of clusters
    labels=[]
    clstr_cntrs=[]
    for r in range(50):
        #print r
        kmeans= KMeans(init= 'k-means++', n_clusters= k, n_init= 10)
        kmeans.fit(qs)
        labels.append(kmeans.predict(qs))
        clstr_cntrs = kmeans.cluster_centers_
        #print clstr_cntrs
        
        for i in range(k):
            cntrs.write(str(clstr_cntrs[i][0])+'\t'+ str(clstr_cntrs[i][1])+'\t' + str(clstr_cntrs[i][2])+ '\t')
        cntrs.write('\n')

    cntrs.close()
# see newplots.py for plotting the results

### Now do the clustering using K-Means

clstr_name= "c3_ew_hwhm_bal"
k=6 #number of clusters
kmeans= KMeans(init= 'k-means++', n_clusters= k, n_init= 10)
kmeans.fit(qs)
labels= kmeans.predict(qs)

## save the clustering results: the subset joined with the label (which point to the cluster the object belongs to) and the object name from the catalog.
clstr_with_names= np.column_stack((qs, labels, sdss_name))
save("./clusters/"+clstr_name+"_"+str(k)+"clstrs_name.npy", clstr_with_names) #save

clstr_with_labels= np.column_stack((qs, labels))
save("./clusters/"+clstr_name+"_"+str(k)+"clstrs.npy", clstr_with_labels)

### to read this array, use load(file name). All elements in the array will have dtype= S18.
### to use them I need to convert to floats. use new_array= old_array.astype(desired_dtype). dtype= float64


""" make median composite spectra for each cluster
"""

compos= [] # list of arrays, each array is a composite that represents one cluster
spec_num= [] # number of objects in each composite (cluster) to be used in the plotting

for c in range(k):
    cluster= tt[labels==c]
    clust_spec= np.arange(1100, 4000, 0.5) # wavelength -this is how it becomes after rebinning
    # t1= time.time()
    # print "t1=", t1
    
    for q in range(len(cluster)):
        name='./new_proc_data/spec-'+str(cluster['PLATE'][q])+'-'+str(cluster['MJD'][q])+'-'+str(cluster['FIBERID'][q]).zfill(4)+'_proc.fits'
        spec=fits.open(name)
        flx= spec[0].data[1]
        #spec.close()
        #   flxx= np.nan_to_num(flx) #making sure the flux array has no NAN or INF
        wlen= spec[0].data[0]
        norm_flx= flx/np.median(flx[2360:2390]) # normalize spectra
        clust_spec= np.vstack((clust_spec, norm_flx)) # 2D array. 1st row: restframe wavelength, other rows have corrected fluxes of spectra from clusters (one for each row)
        del spec
    
    print "cluster", c+1, "has", len(clust_spec[1:]), "objects"
    
    spec_num.append(len(clust_spec[1:]))
    
    clipped_compo=[]
    for i in range(clust_spec.shape[1]):
        
        y= sigmaclip(clust_spec[1:,i], 3, 3)
        m=median(y[0])
        clipped_compo.append(m)
    
    compos.append(clipped_compo) # list with the composites (compos[0] is composite from 1st cluster, compos[1] 2nd cluster,...)


#save the composites as fits files

for i,j in zip(range(1,k+1), spec_num):
    spec_name= "./composites/"+clstr_name+"_"+str(k)+"clstrs"+str(i)+".fits" #assumes there is a directory called composites in the working directory
    spec_file= np.vstack((wlen,compos[i-1]))
    hdu= fits.PrimaryHDU(spec_file)
    hdr= hdu.header
    hdr.set('SPEC_NUMBER', j)
    hdr.set('COMPOSITE', clstr_name)
    hdr.set('PARAMETERS USED', 'EW, RHWHM, BHWHM')
    hdu.writeto(spec_name)
    #hdu.close()

#generate tables with the number of objects in each cluster for each clustering run

clstr_num= [('c3', 3), ('c3', 4), ('c3', 5), ('c3', 6)]

tbl_file= open("clstrs_num_tbl_bal.txt", 'wr')
tbl_file.write("Line & k & 1 & 2 & 3 & 4 & 5 & 6\n")

for o in clstr_num:
    print o[0], o[1]
    spec_numbers=[]
    tbl_file.write(o[0] + "\t" )
    for pp in range(1, o[1]+1):
        sp= fits.open("./composites/"+o[0]+"_ew_hwhm_bal_"+str(o[1])+"clstrs"+str(pp)+".fits")
        spec_numbers.append(sp[0].header['SPEC_NUMBER'])
    spec_numbers_ordered= sorted(spec_numbers, reverse= True)
    print spec_numbers_ordered
    for s in spec_numbers_ordered:
        print s
        tbl_file.write(str(s) + " & ")
    tbl_file.write('\n')
tbl_file.close()

