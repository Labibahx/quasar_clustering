""" plots similar to what is in lstr_plots.py but ordered according to the number of objects in each cluster 
May 8 2015"""

import numpy as np
from astropy.io import fits
from astropy.table import Table
from sklearn.cluster import KMeans
from mpl_toolkits.mplot3d import Axes3D
from operator import itemgetter
from glob import glob
import seaborn as sns


ss= np.load(dr10qsample.npy)


## each cluster is saved into a numpy 2D array which includes the clustering parameters and the name of the sdss object. To get the other parameters for the same object that were not used in the clustering, I can cross match the cluster 2D array with the full subsample array ss (defined in line 14 here).
## each cluster is represented by a composite spectrum with the number of objects in each cluster given in the FITS header (keyword:'SPEC_NUMBER')



def line_profile(line, line_name, k):

    """ plot profiles for lines in the same cluster (line+k, e.g, c4, k=4) in 6 panels for: Ly alpha, Si IV, C IV, He II, (Al III, Si III], C III], Mg II
        param: 
        line: c4, c3 or mg2 as str
        line_name: CIV, CIII, or MgII as str
        k: number of clusters (3, 4, ...)
        """

    compo_name= "./composites/"+line+"_ew_hwhm_"+str(k)+"*.fits"
    compos= glob(compo_name)


    compo_list= []
    for obj in compos:
        spec= fits.open(obj)
        num_obj= spec[0].header['SPEC_NUMBER']
        compo_list.append([obj, num_obj])
    ordered_compos= sorted(compo_list, key= itemgetter(1), reverse= True)

    print ordered_compos


    fig= figure(figsize=(14,8))
    sns.set_style("ticks")
    fig1= fig.add_axes([0., 0., 1, 1])
    fig1.set_axis_off()
    fig1.set_xlim(0, 1)
    fig1.set_ylim(0, 1)
    fig1.text(.07, 0.5, r"Normalized Flux (erg s$^{-1}$ cm$^{-1}$ $\AA^{-1}$)", rotation='vertical', horizontalalignment='center', verticalalignment='center', fontsize= 18)
    fig1.text(0.5, 0.01, r"Wavelength ($\AA$)", rotation='horizontal', horizontalalignment='center', verticalalignment='center', fontsize= 18)

    dx_list= [(1160, 1265), (1350, 1450), (1500, 1600), (1590, 1690), (1810, 1950), (2750, 2850)]
    dy_list= [(0.75, 2.3), (0.75, 1.5), (0.75, 3.0), (0.75, 1.2), (0.75, 1.8), (0.75, 1.6)]

    line_mark= [[1215.7, 1240], [1396.8], [1549], [1640, 1663.5], [1857, 1892, 1908], [2800]]
    line_label= [[r'Ly$\alpha$', 'N V'], ['Si IV'], ['C IV'], ['He II', 'O III]'], ['Al III', 'Si III]', 'C III]'], ['Mg II']]

    alphabet_list = ['a', 'b', 'c', 'd', 'e', 'f']
    compo_labels= [line_name+"-"+ a for a in alphabet_list]

    clr_ls= ['orange', 'navy', 'mediumvioletred','seagreen', '0.5', 'khaki', 'cornflowerblue', 'brown' , 'olive', 'purple']

    for (p,dx, dy, lm, lb) in zip(range(1,8), dx_list, dy_list, line_mark, line_label):
        ax= fig.add_subplot(2,3,p)
        ax.axes.get_xaxis().set_ticks([1150, 1175, 1200, 1225, 1250, 1275, 1300, 1325, 1350, 1375, 1400, 1425, 1450, 1475, 1500, 1525, 1550, 1575, 1600, 1625, 1650, 1675, 1700, 1725, 1750, 1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2700, 2725, 2750, 2775, 2800, 2825, 2850])
        ax.axes.get_yaxis().set_ticks([.2, .4, .6, .8, 1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8])
        
        xlim(dx)
        ylim(dy)

        for ll in range(len(line_mark[p-1])):
        
            axvline(lm[ll], ls=':', c='k')
            text(lm[ll]+0.7, dy[1]-dy[1]/20, lb[ll], rotation= 'vertical', fontsize= 12)

        ii= dy_list[0][1]-.2
        #props= dict(boxstyle='round', alpha=0.5)
        
        for (sp, clr, clab) in zip(ordered_compos, clr_ls, compo_labels):
            n= sp[1]
            if sp[1] >300:
                spec= fits.open(sp[0])
                wlen= spec[0].data[0]
                flx= spec[0].data[1]
                
                plot(wlen, flx/flx[(dx[0]-1100)*2], c= clr, lw= 2)
                ii-=0.15
                if p==1:
                    ax.text(1162, ii, clab+", N="+ str(n), color= clr) #bbox=props

#####
#####

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
    
    w= np.arange(1100, 4000, 0.5) #wavelength array

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
    
    """ line: CIV, CIII, MGII
        cluster: e.g., c4_ew_hwhm or mg2_ew_hwhm
        k: number of clusters
        features: EW, BHWHM, RHWHM
        """


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

    xlabel(line+" "+feature2+ " (km/s)")
    ylabel(line+" "+feature3+ " (km/s)")

    t=0
    for o in ordered_clstrs:
        print o[0], o[1]
        
        ax.scatter(bhwhm[c_label==o[0]], rhwhm[c_label==o[0]], c=clr_ls[t], label=str(o[1]), alpha= 0.7)

        t+=1

    legend(scatterpoints=1)

#############################

def plot_reprod(line, k):

    """ make plots with cluster centroids calculated 50 times to show reproducibility.
    read values from text files.
    """
    
<<<<<<< HEAD
    #lines = [('CIV', 3), ('CIV', 4), ('CIV', 5), ('CIII', 3), ('CIII', 4), ('CIII', 5), ('MGII', 3), ('MGII', 4), ('MGII', 5)]
    #for the BAL sample:lines = [('CIII', 3), ('CIII', 4), ('CIII', 5), ('CIII', 6), ('MGII', 3), ('MGII', 4), ('MGII', 5)]
    

    cntrs= loadtxt(line+str(k)+"_bal.txt") #sample with BALs
    #cntrs= loadtxt(line+str(k)+".txt") #use for the non-BAL sample
   
    fig= figure(figsize=(10,8))
    
    subplots_adjust(hspace = .05)
    ax1= fig.add_subplot(311)
    xlim(-4, 54)
    ylabel('EW '+line)
    ax1.set_xticks([10, 20, 30, 40, 50], [10, 20, 30, 40, 50])
    gca().yaxis.set_major_locator(MaxNLocator(nbins=7, prune= 'both'))
    
    for m in range(0, k*3, 3):
        ax1.scatter(range(50), cntrs[:, m], marker='s', edgecolor='k', facecolor='0.5')
    

    ax2= fig.add_subplot(312, sharex= ax1)
    xlim(-4, 54)
    ylabel("BHWHM "+line)
    gca().yaxis.set_major_locator(MaxNLocator(nbins=7, prune= 'both'))
    for m in range(1, k*3, 3):
        ax2.scatter(range(50), cntrs[:, m], marker='o', edgecolor='k', facecolor='0.5')

    ax3= fig.add_subplot(313, sharex= ax1)
    xlim(-4, 54)
    ylabel('RHWHM '+line)
    xlabel("Number of Repeats")
    gca().yaxis.set_major_locator(MaxNLocator(nbins=7, prune= 'both'))
    for m in range(2, k*3, 3):
        ax3.scatter(range(50), cntrs[:, m], marker='^', edgecolor='k', facecolor='0.5')

############
############

def twoD_cluster_kde(cluster_array, line):

    """ plot KDE for the clusters for the BHWHM, RHWHM
    still need to do a lot of work on it. but this should do for now
    cluster_array: a 2D numpy array with the clustering results
    line: a string with the line used for clustering. To be used as a title for the figure."""

    clstr= np.load(cluster_array)
    
  #  cmap_ls=['OrRd', 'PuBu', 'Purples', 'BuGn', 'RdPu', 'gray_r'] 'YlOrBr'
    cmap_ls= ['YlOrBr', 'Blues', 'RdPu', 'Greens', 'Greys', 'Reds']

    sns.set_style("ticks", {'font.family': u'sans-serif'})
   # sns.set(font_scale=1.5)
    
    fig= figure(figsize=(12, 12))
    ax= fig.add_subplot(111)
    
    #xlabel('BHWHM (km/s)', fontsize=18)
    xlabel(r'EW ($\AA$)', fontsize=18)
    ylabel('RHWHM (km/s)', fontsize=18)
    
    xlim(0,50)
    ylim(0,8000)
    
    x, y= [], []
    
    k_ls=[]
    
    clr_ls= ['orange', 'navy', 'mediumvioletred','seagreen', '0.5', 'red', 'cornflowerblue', 'brown' , 'olive', 'purple']
    
    for i in range(max(clstr[:,3].astype(int))+1):
    
        k_ls.append([i, len(clstr[clstr[:,3]==i])])
    
    ord_k_ls= sorted(k_ls, key= itemgetter(1), reverse= True)
    
    print ord_k_ls
    
    clstr_label= [ line+'-a', line+'-b', line+'-c', line+'-d', line+'-e', line+'-f']
    
    u=1
    for j in range(len(ord_k_ls)):
        k= ord_k_ls[j][0]
        print k
        
        u-=0.04
        x =mean(clstr[:,1][clstr[:,3]==k])
        y =mean(clstr[:,2][clstr[:,3]==k])
        n= len(clstr[:,2][clstr[:,3]==k])
        
        #sns.kdeplot(clstr[:,1][clstr[:,3]==k], clstr[:,2][clstr[:,3]==k], n_levels= 15, shade=True, shade_lowest=False, alpha= 0.5, cmap= cmap_ls[j])
        sns.kdeplot(clstr[:,0][clstr[:,3]==k], clstr[:,2][clstr[:,3]==k], n_levels= 10
                    , shade=True, shade_lowest=False, alpha= 0.5, cmap= cmap_ls[j])
        #scatter(x,y, marker= 'x', c='r', s=60)
        text(x, y, clstr_label[j], fontsize= 16)
        text(0.05, u,  clstr_label[j]+", N="+str(n), transform=ax.transAxes, color= clr_ls[j], fontsize= 16)



