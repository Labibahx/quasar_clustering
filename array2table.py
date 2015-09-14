""" 14 sept 2015
quick function to convert the clusters saved as numpy arrays to astropy fits tables.
"""

import numpy as np
from astropy.table import Table


def array2table():


    l_s= [["mg2", "_", "_main"], ["c4", "_", "_main"], ["c3", "_", "_main"], ["c3", "_mixed_", "_mixed"], ["c3", "_bal_", "_bal"]]
    k= [3,4,5,6]
    a_names= []
    t_names= []
    
    for l in l_s:
        
        for kk in k:
                
            a_names.append(l[0]+"_ew_hwhm"+l[1]+str(kk)+"clstrs_name.npy")
            t_names.append(l[0]+"_"+str(kk)+"clstrs"+l[2]+".fits")

    for a,t in zip(a_names, t_names):
        print a
        print t

        a_clstr= np.load("./clusters/"+a)

        t_clstr= Table([a_clstr[:,0], a_clstr[:,1], a_clstr[:,2], a_clstr[:,3], a_clstr[:,4]], \
               names= ('EW', 'BHWHM', 'RHWHM', 'label', 'SDSS_NAME'), \
               dtype= ('float64', 'float64', 'float64', 'int', 'S18'))

        t_clstr.write("./clusters/"+t, format= 'fits')

    return


