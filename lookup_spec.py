

import numpy as np
from astropy.io import fits

from astropy.table import Table, join

from glob import glob

import seaborn as sns


def spec_display(plate, mjd, fiber):
    """ input Plate, MJD, FiberID. Look up the FITS file and plot the spectrum
        assumes we are in the quasar_clustering directory and FITS files are in the new_proc_data directory.
    """
    spec_name= "./new_proc_data/spec-"+str(plate)+"-"+str(mjd)+"-"+str(fiber).zfill(4)+"_proc.fits"

    spec= fits.open(spec_name)
    flx= spec[0].data[1]
    wlen= spec[0].data[0]

    fig= figure(figsize=(16,8))
    plot(wlen, flx, c='k')

    xlabel(r'Wavelength ($\AA$)')
    ylabel('Normalized flux (erg s$^{-1}$ cm$^{-1}$ $\AA^{-1}$)')

    return

###################


def spec_look_up(line, k, sample_name, label):

    """ read a list of sdss names, corss-match list with table with mjd, plate, fiber IDs to generate spectra file name with format mjd-plate-fiber_proc.fits e.g., spec-3587-55182-0691_proc.fits
    the processed spectra (i.e., corrected for Galactic extinction and de-redshifted and normalized. See spec_proc.py)
    
    param:
    line: "mg2", "c4", or "c3"
    k: total number of clusters used
    sample_name: "main", "mixed", or "bal"
    label: the cluster label you want to look at: k=0 --> first cluster, k=1 --> second cluster...
    
        """
    
    if sample_name== "main":
        sample= "_"
    elif sample_name== "mixed":
        sample= "_mixed_"
    elif sample_name == "bal":
        sample="_bal_"


    cluster_array = "./clusters/"+line+"_ew_hwhm"+sample+str(k)+"clstrs_name.npy"

    all_clstrs= np.load(cluster_array)

    compo_name= "./composites/"+line+"_ew_hwhm"+sample+str(k)+"clstrs"+str(label+1)+".fits"

    compo= fits.open(compo_name)
   
    t_name= "sample"+sample+"myflags.fits"
    
    data= Table.read(t_name) # BAL quasars sample
    #data= Table.read("sample_myflags.csv", format='ascii', delimiter=',') #sample (no BALs)
    
    ss = data[data['MY_FLAG'] ==0] # subsample. some objects were flagged if the had missing flux
    
    #corss-match the above two files
    
    clstr_k= ss[(ss['SDSS_NAME'] == all_clstrs[:,4]) & (all_clstrs[:, 3].astype(int) == label)] # only samples in cluster k

    print len(clstr_k)

    spec_files_ls=[]
    sdss_names_ls= []
    redshift_ls= []
    for s in range(len(clstr_k)):
        spectrum_name= "./new_proc_data/spec-"+str(clstr_k['PLATE'][s])+"-"+str(clstr_k['MJD'][s])+"-"+str(clstr_k['FIBERID'][s]).zfill(4)+"_proc.fits"
        spec_files_ls.append(spectrum_name)
        sdss_names_ls.append(clstr_k['SDSS_NAME'][s])
        redshift_ls.append(clstr_k['Z_PCA'][s])
    
    
        ## plot the spectra
    fig= figure(figsize=(16, 8))
    ax=fig.add_subplot(111)

    for (file, name, z) in zip(spec_files_ls, sdss_names_ls, redshift_ls):
        try:
            spec= fits.open(file)
            wavelen= spec[0].data[0]
            flx= spec[0].data[1]
            norm_flx= flx/np.median(flx[2360:2390])
            plot (wavelen, norm_flx, lw= 2, c='r')
            plot(compo[0].data[0], compo[0].data[1], lw= 2, c='k')
            xlim(1190, 3000)
            ylim(-1, 11)
            axvline(1549, ls= ':')
            axvline(2800, ls= ':')
            axvline(1908, ls=':')
            text(2400, 9, "SDSS "+name)
            text(2400, 8, "z= "+'{0:4.3f}'.format(z))
            #print('%.50f' % d)
            #print('{0:.50f}'.format(d))
            print str(file)
            
            resume = input("Press Enter to plot next spectrum on list.")
        
        except SyntaxError:
            pass
            clf()
    
################

def spec_flag(spec_ls, n1, n2):

    """ read a list of spectra and display them. Read input and use as flag (for either low SNR or BAL quasar).
        
        spec_ls: numpy array with the quasar sample as selected in quasar_cluster.py
        flag= 0 keep
        flag= 1 reject
        
        n1: start at line number n1
        n2: stop at line number n2
        
    """

    data= np.load(spec_ls)
    
    sample= data[n1:n2+1]
    print "Looking at lines", n1, "to", n2
    
    flag_ls=[]
    names=[]
    
    wavelen= np.arange(1100, 4000, 0.5)  #wavelength array
    fig= figure(figsize(20,8))
    
    for i in range(len(sample)):
        print "Looking at spectrum number", i+1
        try:
            spectrum_name= "./new_proc_data/spec-"+str(sample['PLATE'][i])+"-"+str(sample['MJD'][i])+"-"+str(sample['FIBERID'][i]).zfill(4)+"_proc.fits"
            spec= fits.open(spectrum_name)
            flx= spec[0].data[1]
            
            plot(wavelen, flx, c= 'k')
            xlim(1300, 3000)
            ylim(-5,20)
            axvline(1397, ls=':', c='r')
            text(1400, 3.5, 'Si IV', rotation= 'vertical')
            axvline(1549, ls=':', c='r')
            text(1551, 3.5, 'C IV', rotation= 'vertical')
            axvline(1640, ls= ':', c='r')
            text(1642, 3.5, 'He II', rotation='vertical')
            axvline(1908, ls=':', c='r')
            text(1910, 3.5, 'C III]', rotation= 'vertical')
            axvline(2800, ls=':', c='r')
            text(2802, 3.5, 'Mg II', rotation= 'vertical')
            text(2500, 18, "SDSS"+sample['SDSS_NAME'][i], color='r')
            text(2500, 16, "Z_PCA: %4.3f" % sample['Z_PCA'][i], color= 'r')
            print "Flags: 0= keep, 1= reject"
            flag= input()
            flag_ls.append(flag)
            names.append(sample['SDSS_NAME'][i])
            resume = input("Press Enter to plot next spectrum on list.")
        
        except SyntaxError:
            pass
            clf()

    close(fig)
    new_array= np.column_stack((names, flag_ls))
    save("myflags_"+str(n1)+"_to_"+str(n2)+".npy", new_array)

### concatenate the flag arrays

flag_arrays= glob('myflags*to*.npy')

f= np.load(flag_arrays[0])

for i in range(1, len(flag_arrays)):
    ff= np.load(flag_arrays[i])
    f= np.concatenate((f, ff), axis=0)

save("myflags.npy" ,f)


#join the sample file with the myflags list

mf= np.load("myflags.npy")
mf_newshape= np.transpose(mf)

mftable= Table([mf_newshape[0],mf_newshape[1]], names= ('SDSS_NAME', 'MY_FLAG') )
mftable.write("myflags.csv")

#join
sample= Table.read("dr10qsample.csv")
mf3= Table.read("myflags.csv")

x= join(sample, mf3, join_type='left').filled(0)

x.write("qsample_myflags.fits")




#---------------------

mf1= np.load("myflags1.npy")
mf2= np.load("myflags2.npy")

flags=[]
names=[]


for i in range(len(mf1)):
    for j in range(len(mf2)):
        if mf1[i][0]== mf2[j][0] :
            #print mf1[i][1].astype(int), mf2[j][1].astype(int)
            
            f= mf1[i][1].astype(int) * mf2[j][1].astype(int)
            
            names.append(mf1[i][0])
            
            flags.append(f)


new_array= np.column_stack((names, flags))
np.savetxt("myflags_full.csv", new_array, delimiter=",")











