"""

This code is outdated and no longer in use.

"""
import numpy as np
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.table import Table, hstack
import sys
from glob import glob


nfiles = len(sys.argv) - 1
month = 'FebMar'  # TEMP
match_criteria = 1.0
sdss_magtype = 'modelMag'

cattypes_to_merge = [sys.argv[x] for x in range(1, nfiles + 1)]
# clustername = filenames[0].split('_')[0]
filenames = glob('*mycatalog*gal.csv')
clustername = filenames[0].split('_')[0]

# Add these values to w1 and w2 mag to get AB mag
# allwise_ab_vega = [2.699, 3.339]  # From Seth
allwise_ab_vega = [2.661, 3.301]    # From EAZY filter files

panstarrs_filter = ['g','r','i','z','y']

sdss = None
panstarrs = None
allwise = None
clash = None

bsdss = None
bsdss_err = None
bpanstarrs = None
bpanstarrs_err = None

csdss = None
cpanstarrs = None
callwise = None
cclash = None

for f in range(nfiles):
    # filesplits = filenames[f].split('_')
    cattype = cattypes_to_merge[f]
    # t = ascii.read(filenames[f])


    # if filesplits[1] == 'mycatalog' and filesplits[2] == 'sdss':
    if cattype == 'sdss':
        t = ascii.read(clustername + '_mycatalog_' + cattype + '_gal.csv')
        zspec_sdss = t['zspec']
        csdss = SkyCoord(t['ra'], t['dec'], unit='deg')
        magstart = np.where(np.array(t.colnames) == 'dec')[0][0] + 1
        magerrend = np.where(np.array(t.colnames) == 'filter_count')[0][0]
        sdss = t[t.colnames[magstart:magerrend]]
        filts = sdss.colnames[::2]
        bsdss = []
        bsdss_err = []
        for i in range(len(filts[:-5])):
            bsdss.append(sdss.meta['b_' + filts[i].replace('MAG', '')])
            bsdss_err.append(sdss.meta['err_b_' + filts[i].replace('MAG', '')])

    # if filesplits[1] == 'mycatalog' and filesplits[2] == 'panstarrs':
    if cattype == 'panstarrs':
        t = ascii.read(clustername + '_mycatalog_' + cattype + '_gal.csv')
        zspec_panstarrs = t['zspec']
        cpanstarrs = SkyCoord(t['ra'], t['dec'], unit='deg')
        magstart = np.where(np.array(t.colnames) == 'dec')[0][0] + 1
        magerrend = np.where(np.array(t.colnames) == 'filter_count')[0][0]
        panstarrs = t[t.colnames[magstart:magerrend]]
        filts = panstarrs.colnames[::2]
        bpanstarrs = []
        bpanstarrs_err = []
        for i in range(len(filts[:-5])):
            bpanstarrs.append(panstarrs.meta['b_' + filts[i].replace('MAG', '')])
            bpanstarrs_err.append(panstarrs.meta['err_b_' + filts[i].replace('MAG', '')])


    if cattype == 'allwise':
        t = ascii.read(clustername + '_allwise_radec.txt')
        callwise = SkyCoord(t['ra'], t['dec'], unit='deg')
        t.rename_columns(('w1mpro', 'w1sigmpro', 'w2mpro', 'w2sigmpro'), ('allwise_w1MAG', 'allwise_w1ERR', 'allwise_w2MAG', 'allwise_w2ERR'))
        t['allwise_w1MAG'] = t['allwise_w1MAG'] + allwise_ab_vega[0]
        t['allwise_w2MAG'] = t['allwise_w2MAG'] + allwise_ab_vega[1]
        allwise = t['allwise_w1MAG', 'allwise_w1ERR', 'allwise_w2MAG', 'allwise_w2ERR']


bsdss = np.array(bsdss)
bsdss_err = np.array(bsdss_err)
bpanstarrs = np.array(bpanstarrs)
bpanstarrs_err = np.array(bpanstarrs_err)

# TODO figure out a way to match panstarrs to sdss
# b_diff = bsdss - bpanstarrs
# if sdss != None and panstarrs != None:



filters_all = []
cats0 = [sdss, panstarrs, allwise, clash]
if sdss is not None:
    filters_all.extend(sdss.colnames[::2])
    if panstarrs is not None:
        filters_all.extend(panstarrs.colnames[-10:][::2])

elif panstarrs is not None:
    filters_all.extend(panstarrs.colnames[::2])

if allwise is not None:
    filters_all.extend(allwise.colnames[::2])



# Set up table for filters information
eazy_filtname_all = ['F348', 'F349', 'F350', 'F351', 'F352', 'F353', 'F354', 'F355', 'F356', 'F357', 'F358', 'F359', 'F360', 'F361', 'F362', 'F363', 'F364', 'F365', \
        'F156', 'F157', 'F158', 'F159', 'F160', \
        'F244', 'F245', \
        'F334', 'F335', 'F336', 'F337', 'F338', \
        'F113', 'F114', 'F115', 'F285', 'F283', 'F117', 'F118', 'F222']

if month == 'FebMar':
    # Feb & Mar dictionary used for translating filters to eazy filter names
    filt_Dict_all = {'1swMAG': 'F350', '3nwMAG': 'F351', '1seMAG': 'F352', '2nwMAG': 'F353', \
        '3neMAG': 'F354', '2neMAG': 'F355', '3swMAG': 'F356', '2swMAG': 'F357', \
        '3seMAG': 'F358', '2seMAG': 'F359', '4nwMAG': 'F360', '4neMAG': 'F361', \
        '4swMAG': 'F362', '4seMAG': 'F363', '1neMAG': 'F364', '1nwMAG': 'F365', \
        'sdss_uMAG': 'F156', 'sdss_gMAG': 'F157', 'sdss_rMAG': 'F158', 'sdss_iMAG': 'F159', 'sdss_zMAG': 'F160', \
        'allwise_w1MAG': 'F244', 'allwise_w2MAG': 'F245', \
        'panstarrs_gMAG': 'F334', 'panstarrs_rMAG': 'F335', 'panstarrs_iMAG': 'F336', 'panstarrs_zMAG': 'F337', 'panstarrs_yMAG': 'F338'}
elif month == 'Jan':
    # Jan dictionary used for translating filters to eazy filter names
    filt_Dict_all = {'1swMAG': 'F348', '1neMAG': 'F349', '1nwMAG': 'F350', '3nwMAG': 'F351',
        '1seMAG': 'F352', '2nwMAG': 'F353', '3neMAG': 'F354', '2neMAG': 'F355',
        '3swMAG': 'F356', '2swMAG': 'F357', '3seMAG': 'F358', '2seMAG': 'F359',
        '4nwMAG': 'F360', '4neMAG': 'F361', '4swMAG': 'F362', '4seMAG': 'F363', \
        'sdss_uMAG': 'F156', 'sdss_gMAG': 'F157', 'sdss_rMAG': 'F158', 'sdss_iMAG': 'F159', 'sdss_zMAG': 'F160', \
        'allwise_w1MAG': 'F244', 'allwise_w2MAG': 'F245', \
        'panstarrs_gMAG': 'F334', 'panstarrs_rMAG': 'F335', 'panstarrs_iMAG': 'F336', 'panstarrs_zMAG': 'F337', 'panstarrs_yMAG': 'F338'}


# Translate filter names in my table columns to eazy column names, and add two column names id and z_spec
eazy_filts = [filt_Dict_all[filters_all[x]] for x in range(len(filters_all))]
hdrs = []
for x in range(len(eazy_filts)):
    hdrs = hdrs + [eazy_filts[x]] + [eazy_filts[x].replace('F','E')]


# Case 1: if both sdss and panstarrs are available
if sdss is not None and panstarrs is not None:

    idxp, d2dp, d3dp = csdss.match_to_catalog_sky(cpanstarrs)
    hp = d2dp.arcsec < match_criteria
    idxa, d2da, d3da = csdss.match_to_catalog_sky(callwise)
    ha = d2da.arcsec < match_criteria

    nobj1 = len(csdss) 
    tab = Table(np.full((nobj1, len(hdrs)), -99.0), names=hdrs)
    tab[hp] = hstack([sdss[hp], panstarrs[idxp[hp]], Table(np.full((len(csdss[hp]), 4), -99.0))])
    tab[~hp] = hstack([sdss[~hp], Table(np.full((len(csdss[~hp]), 4 + 2 * len(panstarrs_filter)), -99.0))])

    tab['F244'][ha] = allwise['allwise_w1MAG'][idxa[ha]]
    tab['E244'][ha] = allwise['allwise_w1ERR'][idxa[ha]]
    tab['F245'][ha] = allwise['allwise_w2MAG'][idxa[ha]]
    tab['E245'][ha] = allwise['allwise_w1ERR'][idxa[ha]]

    tab.add_columns([np.arange(len(tab)), zspec_sdss], names=['id', 'z_spec'], indexes = [0, 0])
    ascii.write(tab, output=clustername + '_merged_sdss_panstarrs_allwise_eazy.cat', format='commented_header', delimiter='\t', overwrite=True)
    

# Case 2: if only one of sdss/panstarrs is available
elif sdss is not None or panstarrs is not None:
    if sdss is not None:
        cat_name = 'sdss'
        c0 = csdss
        cat = sdss
        zspec = zspec_sdss
    elif panstarrs is not None:
        cat_name = 'panstarrs'
        c0 = cpanstarrs
        cat = panstarrs
        zspec = zspec_panstarrs

    idx, d2d, d3d = c0.match_to_catalog_sky(callwise)
    ha = d2d.arcsec < match_criteria

    nobj1 = len(cat)
    tab = Table(np.full((nobj1, len(hdrs)), -99.0), names=hdrs)
    tab[ha] = hstack([cat[ha], allwise[idx[ha]]])
    tab[~ha] = hstack([cat[~ha], Table(np.full((len(cat[~ha]), 4), -99.0))])

    tab.add_columns([np.arange(len(tab)), zspec], names=['id', 'z_spec'], indexes = [0, 0])

    ascii.write(tab, output=clustername + '_merged_' + cat_name + '_allwise_eazy.cat', format='commented_header', delimiter='\t', overwrite=True)




