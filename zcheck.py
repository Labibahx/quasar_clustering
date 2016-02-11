## quick script to check redshift inconsistancy among different clusters using the Anderson-Darling test

from scipy import stats

from astropy.table import Table, join


def zcheck(line, sample, k):

    ''' param: line used for clustering: c3, c4, mg2
                sample: main, mixed, or bal
                k: number of clusters: 3,4,5,6
                
    '''

    if sample == "main":
        s= ""
    
    elif sample == "mixed":
        s= "mixed_"

    else:
        s= "bal_"

    t= Table.read("sample_"+s+"myflags.fits") # open table with catalog data


    c= Table.read("./clusters/"+line+"_"+str(k)+"clstrs_"+sample+".fits") # table with clustering results from a combination of K and line

    tt= join(t, c, keys= "SDSS_NAME") # join tables to have both clustering results and data from catalog

    # now calculate and print Anderson-Darling test for the redshift estimates using CIII], MgII, CIV and PCA. PCA is the one I used to shift spectra to restframe.

    f= open("z_match.txt", 'wr')
    f.write("Clstr"+'\t'+ "Num" + '\t'+ "Z_MgII"+ '\t' + "sig"+ '\t'+ \
            "Z_CIII"+ '\t' + "sig"+ '\t' +"Z_CIV"+ '\t' + "sig" + '\n')


    for l in range(k):
        
        #ss = stats.anderson_ksamp([tt['Z_PCA'][tt['label']==l], tt['Z_MGII'][tt['label'] ==l], \
                                 #tt['Z_CIII'][tt['label'] ==l], tt['Z_CIV'][tt['label'] ==l]])
        
        s_mg= stats.ks_2samp(tt['Z_PCA'][tt['label']==l], tt['Z_MGII'][tt['label'] ==l])
        s_c3= stats.ks_2samp(tt['Z_PCA'][tt['label']==l], tt['Z_CIII'][tt['label'] ==l])
        s_c4= stats.ks_2samp(tt['Z_PCA'][tt['label']==l], tt['Z_CIV'][tt['label'] ==l])

        f.write('\t' + '&' +str(len(tt[tt['label'] == l]))+ \
                '\t'+ '&'+'{:5.3f}'.format(s_mg[0])+ '\t' +'&'+'{:5.3f}'.format(s_mg[1])+ \
                '\t'+ '&'+'{:5.3f}'.format(s_c3[0])+ '\t' +'&'+'{:5.3f}'.format(s_c3[1])+ \
                '\t'+ '&'+'{:5.3f}'.format(s_c4[0])+'\t'+ '&'+'{:5.3f}'.format(s_c4[1])+'\n')



