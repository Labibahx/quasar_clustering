""" 
Unsupervised clustering analysis on a sample of SDSS quasars.
link to the Paris et al. 2014 paper (SDSS-DR 10)
http://adslabs.org/adsabs/abs/2014A%26A...563A..54P/
This script is similar to spec_cluster.py but it reads corrected spectra (for extiction and redshift) then performs the clustering analysis and creates composites.
spec_cluster.py includes correcting the spectra, clustering and stacking.
I separated these processes to see if this is what is causing the script to be slow -turns out it is just the same.

The working directory for this script is the main quasar_clustering directory. The catalgue file is (dr10q.fits) is there. The script reads the proccessed fits spectra from the proc_data directory.

====================
Updated on 30 May 2015 to include a limit on the SNR in the selctions. This cut-off brought the size of the sample from 7754 to 4342 quasars.

=======
This is a script that does the same stuff as quasar_cluster.py but written in a more friendly resuable way... hopefully

July 8 2015 - not using this script. Back to quasar_cluster.py

"""

import numpy as np
from astropy.table import Table
from astropy.io import fits
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.decomposition import PCA
from sklearn import metrics
import time
from mpl_toolkits.mplot3d import Axes3D



### Read the SDSS DR10 Quasar catalog

data = Table.read('dr10q.fits')

""" + extract the part with redshift 1.6 > z > 2.1
    + remove the bad measurements and all the other junk (the catalog has EW of -1, -inf or ridiculously large negative numbers. The EW values are >0 for absorption lines as well)
    + I checked some of the objects with EWs < 0 and there seem to be something wrong with the measurements (the negative EW does not mean absorption line).
    + I also put a limit on the ew <2000 as there seems to be some outliers
    + BAL quasars will not included in the main sample
    + limit on the SNR of >3 which reduced the size of the sample to ~1/3 if no SNR is used.
    
    """
ss = data[(data['Z_PCA'] >1.6) & (data['Z_PCA'] <2.1)
          & (data['REWE_CIII'] >0) & (data['ERR_REWE_CIII'] < data['REWE_CIII']/10)
          & (data['REWE_CIV'] >0) & (data['ERR_REWE_CIV'] < data['REWE_CIV']/10)
          & (data['REWE_MGII'] >0) & (data['ERR_REWE_MGII'] < data['REWE_MGII']/10)
          & (data['BAL_FLAG_VI'] ==0) & (data['SNR_1700'] > 3)]

def qclstr(line, k):

    """unsupervised clustering with KMeans
    """


    param_list= []



#var = raw_input("Please enter something: ")
#print "you entered", var








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

###
### a list of lists of the features entering the clustering analysis
###
features= [[c4_ew, c4_bhwhm, c4_rhwhm, c4_amp, c4_fwhm],
           [c3_ew, c3_bhwhm, c3_rhwhm, c3_amp, c3_fwhm],
           [mg2_ew, mg2_bhwhm, mg2_rhwhm, mg2_amp, mg2_fwhm]]

### combine 1D arrays to create a 2D numpy array to perform the clustering analysis on (each row is one quasar, each column is one feature)

# qs= np.column_stack(n for n in features[0][:3]+features[1][:3]+features[2][:3]) # all three lines. can specify which features to use by changing the range in the second []
qs= np.column_stack(n for n in features[2][:-2]) # one line only. 0 can be changed to use a different line


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
sil_score_c4= [] # list to store silhouette scores
sil_score_c3= []
sil_score_mg2= []
sil_score_all= []

for q in num_c:
    kmeans= KMeans(init= 'k-means++', n_clusters= q, n_init= 10)
    kmeans.fit(qs)
    labels= kmeans.predict(qs)
    sc= metrics.silhouette_score(qs, labels)
    print q, sc
    sos_ls.append(kmeans.inertia_)
    sil_score_mg2.append(sc)


#scatter(num_c, sos_ls)

plot(num_c, sil_score_mg2, marker= 'D', color= '0.1', ls='--', label= 'Mg II')
plot(num_c, sil_score_c3, marker= 'v', color= '0.3', ls='-.', label= 'C III]')
plot(num_c, sil_score_c4, marker= 'o', color= '0.5', ls=':', label= 'C IV')
#plot(num_c, sil_score_all, marker= 's', color= '0.7', ls='-', label='3 lines')

text(3.5,0.85, "Features: EW, RHWHM, BHWHM")

legend(numpoints= 1)
ylabel('Silhouette score')
xlabel(r'$K$')
xlim(1.9, 8.1)


"""plot sos for CIV, CIII], MgII, and 3lines clusters on one figure
"""

fig= figure()
feat_ls= [features[0]+features[1]+ features[2], features[0], features[1], features[2]]
lbl_ls= ['3 lines', 'C IV', 'C III', 'Mg II']
marker_ls= ['*', 'D', 'o', 's']
color_ls= ['purple', 'seagreen', 'navy', 'orange']

for g in range(4):
    quasars= np.column_stack(m for m in feat_ls[g])
    num_c= np.arange(1,15) # number of clusters
    sos_ls= [] # list of the sum of distances squared

    for q in num_c:
        kmeans= KMeans(init= 'k-means++', n_clusters= q, n_init= 10)
        kmeans.fit(quasars)
        labels= kmeans.predict(quasars)
        sos_ls.append(kmeans.inertia_)

    scatter(num_c, sos_ls, marker= marker_ls[g], color= color_ls[g] , label= lbl_ls[g])

ylabel('Sum of squares')
xlabel('Number of clusters')
legend(numpoints=1)
savefig('sos_all.pdf')


### Now do the clustering using K-Means

clstr_name= "c4_ew_hwhm"
k=3 #number of clusters
kmeans= KMeans(init= 'k-means++', n_clusters= k, n_init= 10)
kmeans.fit(qs)
labels= kmeans.predict(qs)

## save the clustering results: the subset joined with the label (which point to the cluster the object belongs to) and the object name from the catalog.
clstr_with_names= np.column_stack((qs, labels, sdss_name))
save(clstr_name+"_"+str(k)+"clstrs_name.npy", clstr_with_names) #save

clstr_with_labels= np.column_stack((qs, labels))
save(clstr_name+"_"+str(k)+"clstrs.npy", clstr_with_labels)

### to read this array, use load(file name). All elements in the array will have dtype= S18.
### to use them I need to convert to floats. use new_array= old_array.astype(desired_dtype). dtype= float64


""" make median composite spectra for each cluster
"""

compos= [] # list of arrays, each array is a composite that represents one cluster
spec_num= [] # number of objects in each composite (cluster) to be used in the plotting
for c in range(k):
    cluster= ss[labels==c]
    clust_spec= np.arange(1100, 4000, 0.1) # wavelength -this is how it becomes after rebinning
   # t1= time.time()
   # print "t1=", t1
   
    for q in range(len(cluster)):
        name='./proc_data/spec-'+str(cluster['PLATE'][q])+'-'+str(cluster['MJD'][q])+'-'+str(cluster['FIBERID'][q]).zfill(4)+'_proc.fits'
        spec=fits.open(name)
        flx= spec[0].data[1]
        spec.close()
     #   flxx= np.nan_to_num(flx) #making sure the flux array has no NAN or INF
        wlen= spec[0].data[0]
        clust_spec= np.vstack((clust_spec, flx)) # 2D array. 1st row: restframe wavelength, other rows have corrected fluxes of spectra from clusters (one for each row)
        del spec
    
    print "cluster", c+1, "has", len(clust_spec[1:]), "objects"
#    save(clstr_name+str(c+1)+'.npy', clust_spec)  #save spectra in cluster as 2D numpy array with wavelength in 1st row. Can be later read with load(file name)
    spec_num.append(len(clust_spec[1:]))
    compos.append(np.median(clust_spec[1:], axis=0)) # list with the composites (compos[0] is composite from 1st cluster, compos[1] 2nd cluster,...)

#save the composites as fits files

for i,j in zip(range(1,k+1), spec_num):
    spec_name= clstr_name+"_"+str(k)+"clstrs"+str(i)+".fits"
    hdu= fits.PrimaryHDU(compos[i-1])
    hdr= hdu.header
    hdr.set('SPEC_NUMBER', j)
    hdr.set('COMPOSITE', clstr_name)
    hdr.set('PARAMETERS USED', 'EW, RHWHM, BHWHM')
    hdu.writeto(spec_name)
    #hdu.close()

#compos.append(wlen)
#np.savetxt('3lines_5param.txt', compos, delimiter=',') #save as a 2D array with wavelength and a composite for each cluster



''' generate tables with the number of objects in each cluster for each clustering run (each with different number of clusters 4 to 8).
As these numbers might slightly vary for each run the numbers for the composites I currently have saved might be different than the ones I would get from repeating the clustering analysis again. 
So I will instead read the number of the objects in each composite from the FITS headers of the saved composites
'''

clstr_ls= ['c4', 'c3', 'mg2']

for cls in clstr_ls:
    f= open(cls+"_tbl.txt", 'wr')
    f.write("cluster"+"\t"+"1"+"\t"+"2"+"\t"+"3"+"\t"+"4"+"\t"+ "5"+"\t"+ "6"+"\t"+ "7"+"\t"+ "8" + "\n")
    for h in range(4,9):
        f.write(str(h) + "\t")
        for r in range(1,h+1):
            spec= fits.open(cls+"_5param_"+str(h)+"clstrs"+str(r)+".fits")
            num= spec[0].header['SPEC_NUMBER']
            f.write(str(num) + '\t')
        f.write('\n')
    f.close()


###
### format the tables a bit differently to match the figures below
###

clstr_num= [('mg2', 3), ('mg2', 4), ('c3', 4), ('c4', 3), ('c4', 4)]

tbl_file= open("clstrs_num_tbl.txt", 'wr')
tbl_file.write("Line & k & 1 & 2 & 3 & 4 \n")

for o in clstr_num:
    print o[0], o[1]
    spec_numbers=[]
    tbl_file.write(o[0] + "\t" )
    for pp in range(1, o[1]+1):
        sp= fits.open(o[0]+"_3param_"+str(o[1])+"clstrs"+str(pp)+".fits")
        spec_numbers.append(sp[0].header['SPEC_NUMBER'])
    spec_numbers_ordered= sorted(spec_numbers, reverse= True)
    print spec_numbers_ordered
    for s in spec_numbers_ordered:
        print s
        tbl_file.write(str(s) + " & ")
    tbl_file.write('\n')
tbl_file.close()


""" visualize the clusters -need to do more work on this
    
    make 3D plots to visualize the clusters
    """

### this is what I used to quickly look at clustering results.
### load the 2-D numpy array then
### scatter(c4_6param[:,2].astype(float64), c4_6param[:,5].astype(float64), c=c4_6param[:,6].astype(int))

fig= figure(figsize= (12,8))
ax= fig.add_subplot(111, projection = '3d')

ax.azim= 130
ax.elev = 0

ax.scatter(qs[:, 2], qs[:, 3], qs[:, 4], zdir= 'z', c=labels, marker='.')

ax.set_xlabel('FWHM')
ax.set_xlim3d(500, 8000)

ax.set_ylabel('BHWHM')
ax.set_ylim3d(500, 8000)

ax.set_zlabel('RHWHM')
ax.set_zlim3d(500, 8000)


clstr_cntrs= kmeans.cluster_centers_  # returrns a list of tuples with the coordinates of the cluster center

