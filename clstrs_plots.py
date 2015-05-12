""" various plots for the clusters and composites
"""

import numpy as np
from astropy.io import fits
from astropy.table import Table
from sklearn.cluster import KMeans
from mpl_toolkits.mplot3d import Axes3D


data= Tabl.read('dr10q.fits')

ss = data[(data['Z_PCA'] >1.6) & (data['Z_PCA'] <2.1) & (data['REWE_CIII'] >0) & (data['REWE_CIII'] <2000) & (data['REWE_CIV'] >0) & (data['REWE_CIV'] <2000) & (data['REWE_MGII'] >0) & (data['REWE_MGII'] <2000) & (data['BAL_FLAG_VI'] ==0)]

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

qs= np.column_stack(n for n in features[2])

clstr_name= "mg2_ew_hwhm"
k= 4 #number of clusters
kmeans= KMeans(init= 'k-means++', n_clusters= k, n_init= 10)
kmeans.fit(qs)
labels= kmeans.predict(qs)
clstr_cntrs= kmeans.cluster_centers_ # returrns a list of tuples with the coordinates of the cluster center

clstr_with_labels= np.column_stack((qs, labels))
save(clstr_name+str(k)+".npy", clstr_with_labels)


fig= figure(figsize= (12,8))
ax= fig.add_subplot(111, projection = '3d')

ax.azim= 100
ax.elev = 0

ax.scatter(qs[:,0], qs[:,1], qs[:,4], c=labels, marker='o', s=10)

#ax.scatter(qs[:,2], qs[:,3], zdir= 'z', c='grey', marker='+', s=10)
#ax.scatter(qs[:,3], qs[:,5], zdir= 'x', c= 'grey', marker='+', s=10)
#ax.scatter(qs[:,2], qs[:,5], zdir= 'y', c= 'grey', marker='+', s=10)

ax.scatter(clstr_cntrs[:,0], clstr_cntrs[:,1], clstr_cntrs[:,4], marker='x', s=100, c='orange')

ax.set_xlabel('EW')
ax.set_xlim3d(0, 200)

ax.set_ylabel('AMP')
ax.set_ylim3d(5, 50)

ax.set_zlabel('ALPHA_NU')
ax.set_zlim3d(-3,3)


# plot composites
#-----------------#

''' overplot line profiles in four panels (one for each line: CIV, HeII, CIII], MgII)
    '''

clstr_name= "c4_ew_hwhm"
x= 4 #number of clusters

## get the number of objects in each composite
flx_list=[]
num_obj=[]

for z in range(1,x+1):
    spec= fits.open(clstr_name+"_"+str(x)+"clstrs"+str(z)+".fits")
    flx_list.append(spec[0].data)
    hdr= spec[0].header
    num_obj.append(hdr['SPEC_NUMBER'])
    spec.close()

print num_obj

#fig_line= [[1,2,3,4], [5,6,7,8], [9,10,11,12]]

wlength= np.arange(1100, 4000, 0.1)
clr_ls= ['orange', 'navy', 'mediumvioletred','seagreen', 'brown' ,'cornflowerblue', 'khaki', 'olive', 'purple']
xlimit= [(1500,1600), (1600, 1700), (1835, 1950), (2740, 2850)]
ylimit= [(0.55,2.1), (0.49,0.9), (0.39,0.9), (0.19,0.5)]
label_x= [1549, 1640, 1908, 2800]
other_lines= [1663.5, 1857, 1892]
line_name=["C IV", "He II", "C III]", "Mg II"]

fig= figure()

### axes labels
fig1= fig.add_axes([0., 0., 1, 1])
fig1.set_axis_off()
fig1.set_xlim(0, 1)
fig1.set_ylim(0, 1)
fig1.text(.05, 0.5, r"Normalized Flux (erg s$^{-1}$ cm$^{-1}$ $\AA^{-1}$)", rotation='vertical', horizontalalignment='center', verticalalignment='center')
fig1.text(0.5, 0.05, r"Wavelength ($\AA$)", rotation='horizontal', horizontalalignment='center', verticalalignment='center')
###

fl=[1,2,3,4]

for (y,xl,yl,xlab, lname) in zip(fl, xlimit, ylimit, label_x, line_name):
    
    ax= fig.add_subplot(2,2,y)
    axvline(xlab, ls= ':', c= 'k')
    text(xlab+2, yl[1]-yl[1]/10, lname)
    axvline(1663.5, ls= ':', c= 'k')
    text(1666, yl[1]-yl[1]/10, "O III]")
    axvline(1857, ls= ':', c= 'k')
    text(1858, yl[1]-yl[1]/10, "Al III")
    axvline(1892, ls= ':', c= 'k')
    text(1892, yl[1]-yl[1]/10, "Si III")
    ax.axes.get_xaxis().set_ticks([1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2700, 2750, 2800, 2850])
    ax.axes.get_yaxis().set_ticks([.2, .4, .6, .8, 1, 1.2, 1.4, 1.6, 1.8])
    ax.tick_params(axis='both', which='major', labelsize=10)
    xlim(xl)
    ylim(yl)
    
    
    for z in range(x):
        if num_obj[z] >105:
            ax.plot(wlength, flx_list[z], c= clr_ls[z], lw=2, label=str(num_obj[z]))
            if y==1:
                ax.text(1510, 1.8-z/10., str(num_obj[z]), color= clr_ls[z])

savefig(clstr_name+"_"+str(x)+"clstrs"+".pdf")

""" cross match the clusters array with the data table to make histograms for the other lines and compare spectra of individual objects in each cluster with the composite spectrum
"""

clstr_name= "mg2_ew_hwhm"
k= 3 #number of clusters

clstr_array= np.load(clstr_name+str(k)+"clstrs_name.npy") #saved array

c1= ss[(ss['SDSS_NAME']== clstr_array[:,4]) & (clstr_array[:,3].astype(int) ==0)] #cluster 1
c2= ss[(ss['SDSS_NAME']== clstr_array[:,4]) & (clstr_array[:,3].astype(int) ==1)] #cluster 2
c3= ss[(ss['SDSS_NAME']== clstr_array[:,4]) & (clstr_array[:,3].astype(int) ==2)] #cluster 3
c4= ss[(ss['SDSS_NAME']== clstr_array[:,4]) & (clstr_array[:,3].astype(int) ==3)] #cluster 4

all_clstrs= ss[ss['SDSS_NAME']== clstr_array[:,4]] # cross match the array with the data table to get the columns for the other lines


### make histograms to compare distibutions of parameters in each cluster.

a= 'FWHM_CIII' # change this to make histogram for another parameter

fig= figure()
ax= fig.add_subplot(111)
ax.hist(c1[a], bins=12, histtype='step', normed=True, lw=1.5, ec='orange', label= str(len(c1)))
ax.hist(c2[a], bins=12, histtype='step', normed=True, lw=1.5, ec='navy', label= str(len(c2)))
ax.hist(c3[a], bins=12, histtype='step', normed=True, lw=1.5, ec='mediumvioletred', label= str(len(c3)))
ax.hist(c4[a], bins=12, histtype='step', normed=True, lw=1.5, ec='seagreen', label= str(len(c4)))
ylim(0, 0.001)
xlim(500, 5000)
xlabel(a)
ylabel('Normalized')
legend()

savefig(clstr_name+a+"_hist.pdf")


### create 5 panel figures with a histogram for each of the parameters in a panel

line="CIII"
clstr_name= "c3_ew_hwhm"
k=6

clstr_array= np.load(clstr_name+str(k)+"clstrs_name.npy") #saved array
c1= ss[(ss['SDSS_NAME']== clstr_array[:,4]) & (clstr_array[:,3].astype(int) ==0)] #cluster 1
c2= ss[(ss['SDSS_NAME']== clstr_array[:,4]) & (clstr_array[:,3].astype(int) ==1)] #cluster 2
c3= ss[(ss['SDSS_NAME']== clstr_array[:,4]) & (clstr_array[:,3].astype(int) ==2)] #cluster 3
c4= ss[(ss['SDSS_NAME']== clstr_array[:,4]) & (clstr_array[:,3].astype(int) ==3)] #cluster 4
c5= ss[(ss['SDSS_NAME']== clstr_array[:,4]) & (clstr_array[:,3].astype(int) ==4)] #cluster 5
c6= ss[(ss['SDSS_NAME']== clstr_array[:,4]) & (clstr_array[:,3].astype(int) ==5)] #cluster 6

all_clstrs= ss[ss['SDSS_NAME']== clstr_array[:,4]] # cross match the array with the data table to get the columns for the other lines

fig= figure(figsize=(10,10))
fig1= fig.add_axes([0., 0., 1, 1])
fig1.set_axis_off()
fig1.set_xlim(0, 1)
fig1.set_ylim(0, 1)
fig1.text(.05, 0.5, "Normalized Count", rotation='vertical', horizontalalignment='center', verticalalignment='center')

bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9) #used for labeling features in each panel

param_list=["REWE_"+line, "FWHM_"+line, "RHWHM_"+line, "BHWHM_"+line, "AMP_"+line] #edit this to change the line

ax1= fig.add_subplot(411)
ax1.hist(c1["REWE_"+line], bins=12, histtype='step', lw=1.5, ec='orange', label= str(len(c1)), normed=1)
ax1.hist(c2["REWE_"+line], bins=12, histtype='step', lw=1.5, ec='navy', label= str(len(c2)), normed=1)
ax1.hist(c3["REWE_"+line], bins=12, histtype='step', lw=1.5, ec='mediumvioletred', label= str(len(c3)), normed=1)
ax1.hist(c4["REWE_"+line], bins=12, histtype='step', lw=1.5, ec='seagreen', label= str(len(c4)), normed=1)
ax1.hist(c5["REWE_"+line], bins=12, histtype='step', lw=1.5, ec='brown', label= str(len(c5)), normed=1)
ax1.hist(c6["REWE_"+line], bins=12, histtype='step', lw=1.5, ec='cornflowerblue', label= str(len(c6)), normed=1)
xlim(min(all_clstrs["REWE_"+line]), max(all_clstrs["REWE_"+line])/2)
ylim(0, 0.045)
ax1.axes.get_yaxis().set_ticks([0.01, 0.03])
text(0.75, 0.8, r"EW ($\AA$)", transform= ax1.transAxes, bbox=bbox_props)
legend(loc=4)

ax2= fig.add_subplot(412)
ax2.hist(c1["RHWHM_"+line], bins=12, histtype='step', lw=1.5, ec='orange', label= str(len(c1)), normed=1)
ax2.hist(c2["RHWHM_"+line], bins=12, histtype='step', lw=1.5, ec='navy', label= str(len(c2)), normed=1)
ax2.hist(c3["RHWHM_"+line], bins=12, histtype='step', lw=1.5, ec='mediumvioletred', label= str(len(c3)), normed=1)
ax2.hist(c4["RHWHM_"+line], bins=12, histtype='step', lw=1.5, ec='seagreen', label= str(len(c4)), normed=1)
ax2.hist(c5["RHWHM_"+line], bins=12, histtype='step', lw=1.5, ec='brown', label= str(len(c5)), normed=1)
ax2.hist(c6["RHWHM_"+line], bins=12, histtype='step', lw=1.5, ec='cornflowerblue', label= str(len(c6)), normed=1)
xlim(min(all_clstrs["RHWHM_"+line]), max(all_clstrs["RHWHM_"+line]))
ylim(0, 0.00085)
ax2.axes.get_yaxis().set_ticks([0.0001, 0.0003, 0.0005, 0.0007])
text(0.75, 0.85, "RHWHM (km\s)", transform= ax2.transAxes, bbox=bbox_props)

ax3= fig.add_subplot(413)
ax3.hist(c1["BHWHM_"+line], bins=12, histtype='step', lw=1.5, ec='orange', label= str(len(c1)), normed=1)
ax3.hist(c2["BHWHM_"+line], bins=12, histtype='step', lw=1.5, ec='navy', label= str(len(c2)), normed=1)
ax3.hist(c3["BHWHM_"+line], bins=12, histtype='step', lw=1.5, ec='mediumvioletred', label= str(len(c3)), normed=1)
ax3.hist(c4["BHWHM_"+line], bins=12, histtype='step', lw=1.5, ec='seagreen', label= str(len(c4)), normed=1)
ax3.hist(c5["BHWHM_"+line], bins=12, histtype='step', lw=1.5, ec='brown', label= str(len(c5)), normed=1)
ax3.hist(c6["BHWHM_"+line], bins=12, histtype='step', lw=1.5, ec='cornflowerblue', label= str(len(c6)), normed=1)
xlim(min(all_clstrs["BHWHM_"+line]), max(all_clstrs["BHWHM_"+line]))
ylim(0, 0.00085)
ax3.axes.get_yaxis().set_ticks([0.0001, 0.0003, 0.0005, 0.0007])
text(0.75, 0.8, "BHWHM (km\s)", transform= ax3.transAxes, bbox=bbox_props)

ax4= fig.add_subplot(414)
ax4.hist(c1["FWHM_"+line], bins=12, histtype='step', lw=1.5, ec='orange', label= str(len(c1)), normed=1)
ax4.hist(c2["FWHM_"+line], bins=12, histtype='step', lw=1.5, ec='navy', label= str(len(c2)), normed=1)
ax4.hist(c3["FWHM_"+line], bins=12, histtype='step', lw=1.5, ec='mediumvioletred', label= str(len(c3)), normed=1)
ax4.hist(c4["FWHM_"+line], bins=12, histtype='step', lw=1.5, ec='seagreen', label= str(len(c4)), normed=1)
ax4.hist(c5["FWHM_"+line], bins=12, histtype='step', lw=1.5, ec='brown', label= str(len(c5)), normed=1)
ax4.hist(c6["FWHM_"+line], bins=12, histtype='step', lw=1.5, ec='cornflowerblue', label= str(len(c6)), normed=1)
xlim(min(all_clstrs["FWHM_"+line]), max(all_clstrs["FWHM_"+line]))
ylim(0, 0.0006)
ax4.axes.get_yaxis().set_ticks([0.0001, 0.0003, 0.0005, 0.0007])
text(0.75, 0.8, "FWHM (km\s)", transform= ax4.transAxes, bbox=bbox_props)

majorFormatter = FormatStrFormatter('%1.0e')
#majorLocator   = MultipleLocator()
ax1.yaxis.set_major_formatter(majorFormatter)
ax2.yaxis.set_major_formatter(majorFormatter)
ax3.yaxis.set_major_formatter(majorFormatter)
ax4.yaxis.set_major_formatter(majorFormatter)

savefig(clstr_name+str(k)+"clstrs_hist.pdf")



### overplot composites with individual spectra from the same cluster

#pick example spectrum from the original subsample table (called ss definded above) from each one of the clusters.

obj1='./proc_data/spec-'+str(c1['PLATE'][q])+'-'+str(c1['MJD'][q])+'-'+str(c1['FIBERID'][q]).zfill(4)+'_proc.fits'
spec1=fits.open(ob1)
flx1= spec[0].data[1]

