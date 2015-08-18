""" 13 Aug 2015
a function to repeat the KMeans analysis many times and plot the mean and dispersion of the results of clustering (mean of coordiantes and dispersion as error bars).
the goal is to show that the clustering is stable and the clusters we find are not random.
"""

import numpy as np
from astropy.table import Table
from sklearn.cluster import KMeans
import seaborn as sns

def clstr_reprod(sample, k):
    """ show that K-Means is finding the same cluster centroids even after repeating hte clustering for 50 times.
    param:
    sample: ""main", "mixed", "bal"
    k: number of clusters
    """
    sns.set_style("ticks", {'font.family': u'serif', 'ytick.direction': u'in'})
    sns.set_palette("deep")

    if sample == "main":
        tab_name= "sample_myflags.fits"
        line_ls= [["MGII", 1], ["CIV", 4], ["CIII", 7]]
    else:
        tab_name= "sample_"+sample+"_myflags.fits"
        line_ls= [["CIII", 1]]

    t= Table.read(tab_name)
    tt= t[t['MY_FLAG'] ==0]

    fig= figure(figsize(14,9))
    subplots_adjust(hspace = .01)

    fig1= fig.add_axes([0., 0., 1, 1])
    fig1.set_axis_off()
    fig1.set_xlim(0, 1)
    fig1.set_ylim(0, 1)
    fig1.text(.07, 0.75, r"MgII", rotation='vertical', horizontalalignment='center', verticalalignment='center', fontsize= 14, family= 'serif')
    fig1.text(.07, 0.5, r"CIV", rotation='vertical', horizontalalignment='center', verticalalignment='center', fontsize= 14, family= 'serif')
    fig1.text(.07, 0.25, r"CIII]", rotation='vertical', horizontalalignment='center', verticalalignment='center', fontsize= 14, family= 'serif')

    fig1.text(.24, 0.92, r"EW ($\AA$)", rotation='horizontal', horizontalalignment='center', verticalalignment='center', fontsize= 14, family= 'serif')
    fig1.text(.51, 0.93, r"BHWHM (km/s)", rotation='horizontal', horizontalalignment='center', verticalalignment='center', fontsize= 14, family= 'serif')
    fig1.text(.78, 0.93, r"RHWHM (km/s)", rotation='horizontal', horizontalalignment='center', verticalalignment='center', fontsize= 14, family= 'serif')
    
    fig1.text(0.5, 0.035, r"Num of Repeats", rotation='horizontal', horizontalalignment='center', verticalalignment='center', fontsize= 18, family= 'serif')

    for l in line_ls:
        
        ew= tt["REWE_"+l[0]]
        bw= tt["BHWHM_"+l[0]]
        rw= tt["RHWHM_"+l[0]]

        features= [ew, bw, rw]

        qs= np.column_stack(n for n in features)

        clstr_cntrs=[]
        ew_c, bw_c, rw_c= [], [], []
        for r in range(50):
            kmeans= KMeans(init= 'k-means++', n_clusters= k, n_init= 10)
            kmeans.fit(qs)
            cc=kmeans.cluster_centers_
            for m in range(k):
                ew_c.append(cc[m][0])
                bw_c.append(cc[m][1])
                rw_c.append(cc[m][2])
        
        p =l[1]

        ax1= fig.add_subplot(3,3,p)
        xlim(-1, 51)
        gca().yaxis.set_major_locator(MaxNLocator(nbins=7, prune= 'both'))
        #ax1.ticklabel_format(axis='y', useOffset=False)
        if p in [1,2,3,4,5,6]:
            ax1.set_xticklabels([])
        scatter(range(50), ew_c[:50], marker='^',facecolor='0.7', edgecolor='r')
        scatter(range(50), ew_c[50:100], marker='^', facecolor='0.7', edgecolor='r')
        scatter(range(50), ew_c[100:150], marker='^', facecolor='0.7', edgecolor='r')

        ax2= fig.add_subplot(3,3,p+1)
        xlim(-1, 51)
        gca().yaxis.set_major_locator(MaxNLocator(nbins=7, prune= 'both'))
        if p in [1,2,3,4,5,6]:
            ax2.set_xticklabels([])
        scatter(range(50), bw_c[:50], marker='o', facecolor='0.7', edgecolor='b')
        scatter(range(50), bw_c[50:100], marker='o', facecolor='0.7', edgecolor='b')
        scatter(range(50), bw_c[100:150], marker='o', facecolor='0.7', edgecolor='b')

        ax3= fig.add_subplot(3,3,p+2)
        xlim(-1, 51)
        gca().yaxis.set_major_locator(MaxNLocator(nbins=7, prune= 'both'))
        if p in [1,2,3,4,5,6]:
            ax3.set_xticklabels([])
        scatter(range(50), rw_c[:50], marker='s', facecolor='0.7', edgecolor='')
        scatter(range(50), rw_c[50:100], marker='s',facecolor='0.7', edgecolor='g')
        scatter(range(50), rw_c[100:150], marker='s',facecolor='0.7', edgecolor='g')

    return
