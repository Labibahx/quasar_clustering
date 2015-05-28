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


def profiles(line, k1, k2):
    """ plot line profiles for the clusters in 4 panels
        """
    comp_name1= line+"_ew_hwhm_"+str(k1)+"*.fits"
    comp_name2= line+"_ew_hwhm_"+str(k2)+"*.fits"
    compos1= glob(comp_name1)
    compos2= glob(comp_name2)
    
    compo_list_k1= []
    for obj in compos1:
        spec1= fits.open(obj)
        num_obj1= spec1[0].header['SPEC_NUMBER']
        compo_list_k1.append([obj, num_obj1])
    ordered_compos1= sorted(compo_list_k1, key= itemgetter(1), reverse= True)

    compo_list_k2= []
    for obj in compos2:
        spec2= fits.open(obj)
        num_obj2= spec2[0].header['SPEC_NUMBER']
        compo_list_k2.append([obj, num_obj2])
    ordered_compos2= sorted(compo_list_k2, key= itemgetter(1), reverse= True)

    print ordered_compos1, ordered_compos2

   # return ordered_compos1, ordered_compos2

    fig= figure(figsize=(12,8))
    fig1= fig.add_axes([0., 0., 1, 1])
    fig1.set_axis_off()
    fig1.set_xlim(0, 1)
    fig1.set_ylim(0, 1)
    fig1.text(.07, 0.5, r"Normalized Flux (erg s$^{-1}$ cm$^{-1}$ $\AA^{-1}$)", rotation='vertical', horizontalalignment='center', verticalalignment='center')
    fig1.text(0.5, 0.03, r"Wavelength ($\AA$)", rotation='horizontal', horizontalalignment='center', verticalalignment='center')
    fig1.text(0.92, 0.7, "K = "+str(k1), rotation= 'vertical', horizontalalignment= 'center')
    fig1.text(0.92, 0.3, "K = "+str(k2), rotation= 'vertical', horizontalalignment= 'center')
    
    w= np.arange(1100, 4000, 0.1) #wavelength array

    fl= range(1,5)

    clr_ls= ['orange', 'navy', 'mediumvioletred','seagreen', 'khaki', 'cornflowerblue', 'brown' , 'olive', 'purple']
    xlimit= [(1500,1600), (1600, 1700), (1835, 1950), (2740, 2850)]
    ylimit= [(0.55,2.1), (0.49,0.9), (0.39,0.9), (0.19,0.5)]
    line_mark= [[1549], [1640, 1663.5], [1857, 1892, 1908], [2800]]
    line_label= [['CIV'], ['HeII', 'OIII]'], ['AlIII', 'SiIII', 'CIII]'], ['MgII']]


    for (y,xl,yl, l_lab, l_mark) in zip(fl, xlimit, ylimit, line_label, line_mark):
        
        ax= fig.add_subplot(2,4,y)
        ax.axes.get_xaxis().set_ticks([1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2700, 2750, 2800, 2850])
        ax.axes.get_yaxis().set_ticks([.2, .4, .6, .8, 1, 1.2, 1.4, 1.6, 1.8])
        ax.tick_params(axis='both', which='major', labelsize=10)
        xlim(xl)
        ylim(yl)
        
        for ll in range(len(l_lab)):
            
            ax.axvline(l_mark[ll], ls= ':', c= 'k')
            ax.text(l_mark[ll]-4, yl[1]-yl[1]/20, l_lab[ll], rotation= 'vertical', horizontalalignment= 'center')


        ii=1.85
        for (o, clr) in zip(ordered_compos1, clr_ls):
            n= o[1]
            if n >105:
                spec= fits.open(o[0])
                flx= spec[0].data
                ax.plot(w, flx, c= clr, lw=2, label= str(n))
                ii-=0.15
                if y==1:
                    ax.text(1510, ii, str(n), color= clr)

                
    fl= range(5,9)


    clr_ls= ['orange', 'navy', 'mediumvioletred','seagreen', 'khaki', 'cornflowerblue', 'brown' , 'olive', 'purple']
    xlimit= [(1500,1600), (1600, 1700), (1835, 1950), (2740, 2850)]
    ylimit= [(0.55,2.1), (0.49,0.9), (0.39,0.9), (0.19,0.5)]
    line_mark= [[1549], [1640, 1663.5], [1857, 1892, 1908], [2800]]
    line_label= [['CIV'], ['HeII', 'OIII]'], ['AlIII', 'SiIII', 'CIII]'], ['MgII']]

    
    for (y,xl,yl, l_lab, l_mark) in zip(fl, xlimit, ylimit, line_label, line_mark):
    
        ax= fig.add_subplot(2,4,y)
        ax.axes.get_xaxis().set_ticks([1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2700, 2750, 2800, 2850])
        ax.axes.get_yaxis().set_ticks([.2, .4, .6, .8, 1, 1.2, 1.4, 1.6, 1.8])
        ax.tick_params(axis='both', which='major', labelsize=10)
        xlim(xl)
        ylim(yl)
        
        for ll in range(len(l_lab)):
            
            ax.axvline(l_mark[ll], ls= ':', c= 'k')
            ax.text(l_mark[ll]-4, yl[1]-yl[1]/20, l_lab[ll], rotation= 'vertical', horizontalalignment= 'center')


        ii= 1.85
        for (o, clr) in zip(ordered_compos2, clr_ls):
            n= o[1]
            if o[1] > 105:
                spec= fits.open(o[0])
                flx= spec[0].data
                ax.plot(w, flx, c=clr, lw=2, label= str(n))
                ii-=0.15
                if y==5:
                    ax.text(1510, ii, str(n), color= clr)


#####################


""" 2D scatter plots for clusters. read files from saved 2d numpy arrays
    the aray has the features in the "columns". the last column has the labels indicating which cluster each sample (row) belongs to.
    line: the emission line used in clustering: MgII, CIII], or CIV
    cluster: line_ew_hwhm -in this case. can refere to any other cluster name (eg, mg2_ew_hwhm)
    k: number of clusters used in creating the loaded numpy array. The labels column in the array goes from 0 to k-1
    feature 1: one of the features used in the clustering analysis. This will be plotted on the x-axis
    feature 2: to be plotted on the y-axis
    
    in the ew_hwhm clusters, the features are EW, BHWHM, and RHWHM given in the 0th, 1st and 2nd vectors on the array.
    
    """

def two_d_scatter(line, cluster, k, feature1, feature2, feature3):


    clstr_name= cluster+"_"+str(k)+"clstrs.npy"
    clstr= np.load(clstr_name)
    
    ew= clstr[:,0] #EW
    bhwhm= clstr[:,1] #BHWHM
    rhwhm= clstr[:,2] #RHWHM
    c_label= clstr[:,3] #label
    
    clstr_length=[]

    for c in range(0,k):
        clstr_length.append([c,len(clstr[c_label==c])])

    ordered_clstrs =sorted(clstr_length, key=itemgetter(1), reverse= True)
    print ordered_clstrs
    
    clr_ls= ['orange', 'navy', 'mediumvioletred','seagreen', 'khaki', 'cornflowerblue', 'brown' , 'olive', 'purple']

    fig= figure(figsize(8,8))
    ax= fig.add_subplot(111)

    text(0.05, 0.1, 'K= '+ str(k)+ ', Clustering features:'+' '+line+' ('+ feature1+','+feature2+','+ feature3+')',
         horizontalalignment='left', verticalalignment='center',
         transform=ax.transAxes, color= 'black')

    xlabel(line+" "+feature1+ " (km/s)")
    ylabel(line+" "+feature2+ " (km/s)")

    t=0
    for o in ordered_clstrs:
        print o[0], o[1]
        
        ax.scatter(bhwhm[c_label==o[0]], rhwhm[c_label==o[0]], c=clr_ls[t], label=str(o[1]))

        t+=1

    legend(scatterpoints=1)







