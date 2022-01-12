"""

Usage:

    In python script directory:
    $ python combine_catalog_calibration.py [path to catalogs directory] [sdss/panstarrs] (optional arguments)

    In catalogs directory:
    $ python path_to_script/combine_catalog_calibration.py [sdss/panstarrs] (optional arguments)

Optional arguments and default values:

    match_criteria=1.0
    match=True
    plot=True
    x_cutoff=1.05   (or xcutoff=1.05)
    ransac=200      (or ransac_iterations=200)
    starmag=hybrid      (can be hybrid, auto, psf)
    sdssmag=modelmag    (can be modelmag, cmodelmag, petromag)
    psc=0.83
    my_gals_magtype_sdss=MAG_AUTO       (can be MAG_AUTO, MAG_PETRO, MAG_MODEL, MAG_HYBRID, MAG_DISK, MAG_SPHEROID)
    my_gals_magtype_panstarrs=MAG_AUTO  (can be MAG_AUTO, MAG_PETRO, MAG_MODEL, MAG_HYBRID, MAG_DISK, MAG_SPHEROID)

Required files:

    cluster_gal_sdss_radec.csv
    cluster_star_sdss_radec.csv
    cluster_panstarrs_radec.csv
    cluster_panstarrs_extinction.txt
    cluster_allwise_radec.txt                                       https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-scan?mission=irsa&submit=Select&projshort=WISE

This script search for 1ne.cat~4sw.cat and match them with SDSS/PanSTARRS catalog, create mycatalog then merged with ALLWISE catalog.
This script creates:
    clustername_mycatalog_sdss(panstarrs)_gal.csv                   calibrated and selected magnitudes, will change even if match=False
    clustername_mycatalog_sdss(panstarrs)_star_psf(auto).csv        calibrated psf and auto magnitudes for stars catalogs
    clustername_mycatalog_sdss(panstarrs)_gal(star)_eazy.cat        calibrated mycatalog as EAZY input format
    clustername_mycatalog_merged_sdss(panstarrs)_allwise_eazy.cat   calibrated mycatalog matched/merged with ALLWISE as EAZY input format
    clustername_mycatalog_sdss(panstarrs)_gal(star)_original.ecsv   uncalibrated magnitudes, just matched with SDSS/PanSTARRS, will be read by this script if match=False
    clustername_mycatalog_sdss(panstarrs)_gal(star)_ds9.csv         csv file with just ra/dec for DS9 catalog tool
    clustername_mycatalog_gal(star)_w_sdss_spec_radec.txt           txt file with ra/dec of objects that have SDSS spectra, can be upload to:
                                                                    https://dr16.sdss.org/optical/spectrum/search to get the spectra fits files


"""
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
import numpy as np
from astropy.table import Table, vstack, hstack
from astropy.coordinates import concatenate, SkyCoord
from astropy import units as u
from astropy.io import ascii, fits
import FoFCatalogMatching
import sys
import os
from glob import glob
import matplotlib.pyplot as plt
from sklearn import linear_model
import time
import multiprocessing as mp
import pathlib
import os.path

# plt.ion()
sys.dont_write_bytecode=True


# ========================================== PARAMETERS ================================================================

# Default values
# Parameters
barypeakdiff = 1.0
match_criteria = 1.0
sdss_calibration_gals_magtype = 'modelMag'
my_gals_magtype_sdss = 'MAG_AUTO'           # MAG_AUTO, MAG_PETRO, MAG_MODEL, MAG_HYBRID, MAG_DISK, MAG_SPHEROID
my_gals_magtype_panstarrs = 'MAG_AUTO'      # MAG_AUTO, MAG_PETRO, MAG_MODEL, MAG_HYBRID, MAG_DISK, MAG_SPHEROID
my_stars_magtype = 'hybrid'                 # mag type used for star calibration, auto/psf/hybrid
apply_extinction = True
read_only = False
make_plot = True
null_fillval = -99.0
apply_b_spectra_correction = True

# Calibration Criteria
x_cutoff = 1.05
ransac_iterations = 200
panstarrs_classification_psc = 0.83

# Plotting options
cmap = 'viridis'
linecolor = 'coral'
linealpha = 1.0             # fitted mx+b line in each filter alpha
markersize0 = 15            # all data points markersize 
markersize1 = 20            # points that are used to fit witn cmap
linewidth0 = 1              # r-i vs g-r plot linewidth
linewidth1 = 1              # fitted mx+b in each filter linewidth
small_title_fontsize = 8    # fontsize of filter wavelength
xlabelsize = 12             # labelsize for g-r
ylabelsize = 12             # labelsize for narrowband-r
titlesize = 12              # title size of the whole plot
ticksize = 6                # ticksize of each filter plot
xlim_l = -0.25
xlim_h = 2.25


# Add these values to w1 and w2 mag to get AB mag
# allwise_ab_vega = [2.699, 3.339]  # From Seth
allwise_ab_vega = [2.661, 3.301]    # From EAZY filter files


peakmax_c = 60000           # Not currently in use
peakmin_c = 2000            # Not currently in use
flag_c = 1
calibrationmagmax = 30
calibrationmagmin = 15
fwhmmax_c = 6               # Not currently in use
fwhmmin_c = 1               # Not currently in use


month_Dict = {'abell1576': 'FebMar', 'abell370': 'Jan', 'abell611': 'FebMar', 'macs0329': 'Jan', 'macs0717': \
                'FebMar', 'macs1115': 'Jan', 'macs1149': 'Jan', 'rxj1532': 'FebMar', 'zwicky1953': 'FebMar'}

# ======================================================================================================================


path = ''
script_dir = pathlib.Path(__file__).parent.absolute()
working_dir = pathlib.Path().absolute()


if script_dir == working_dir:
    if len(sys.argv) == 2:
        print(__doc__)
        sys.exit()
    else:
        path = sys.argv[1]
        calibration_cat = sys.argv[2]
        rest_arg_start = 3
        if path[-1] != '/':
            path = path + '/'
else:
    if len(sys.argv) == 1:
        print(__doc__)
        sys.exit()
    else:
        calibration_cat = sys.argv[1]
        rest_arg_start = 2

if calibration_cat != 'sdss' and calibration_cat != 'panstarrs':
    print('Incorrect calibration catalog name!')
    sys.exit()



nargs = len(sys.argv)
for x in range(rest_arg_start, nargs):
    if len(sys.argv[x].split('=')) == 1:
        print('Use "=" to separate argument keywords and their value')
        sys.exit()

    arg_type = sys.argv[x].split('=')[0].lower()
    arg_TF = sys.argv[x].split('=')[-1].lower()
    
    if arg_type == 'match':
        if arg_TF == 'true':
            read_only = False
        elif arg_TF == 'false':
            read_only = True

    if arg_type == 'plot':
        if arg_TF == 'true':
            make_plot = True
        elif arg_TF == 'false':
            make_plot = False

    if arg_type == 'match_criteria':
        match_criteria = float(arg_TF)

    if arg_type == 'x_cutoff' or arg_type == 'xcutoff':
        x_cutoff = float(arg_TF)

    if arg_type == 'starmag':
        my_stars_magtype = arg_TF
        if my_stars_magtype != 'psf' and my_stars_magtype != 'auto' and my_stars_magtype != 'hybrid':
            print('Incorrect star type')
            sys.exit()
    
    if arg_type == 'sdssmag':
        if arg_TF == 'modelmag':
            sdss_calibration_gals_magtype = 'modelMag'
        elif arg_TF == 'cmodelmag':
            sdss_calibration_gals_magtype = 'cmodelMag'
        elif arg_TF == 'petromag':
            sdss_calibration_gals_magtype = 'petroMag'
        else:
            print('Incorrect sdss mag type')
            sys.exit()

    if arg_type == 'ransac' or arg_type == 'ransac_iterations':
        ransac_iterations = int(arg_TF)

    if arg_type == 'psc':
        panstarrs_classification_psc = float(arg_TF)

    if arg_type == 'my_gals_magtype_sdss':
        my_gals_magtype_sdss = str(arg_TF).upper()

    if arg_type == 'my_gals_magtype_panstarrs':
        my_gals_magtype_panstarrs = str(arg_TF).upper()


globfilenames = glob(path + 'cutout_stack_*fits')
globfilenames.extend(glob(path + '*radec.csv'))
clustername = globfilenames[0].split('/')[-1].split('_')[2]

month = month_Dict[clustername]

if calibration_cat == 'sdss':
    if (not os.path.exists(path + clustername + '_gal_sdss_radec.csv')) or (not os.path.exists(path + clustername + '_star_sdss_radec.csv')):
        print("SDSS catalogs don't exist")
        sys.exit()
elif calibration_cat == 'panstarrs':
    if (not os.path.exists(path + clustername + '_panstarrs_radec.csv')) or (not os.path.exists(path + clustername + '_panstarrs_extinction.txt')):
        print("PanSTARRS catalog doesn't exist")
        sys.exit()

if not os.path.exists(path + clustername + '_allwise_radec.txt'):
    print("ALLWISE catalog doesn't exist")
    sys.exit()


print('calibration_cat={}\nmatch={}\nplot={}\nstarmag={}\nmatch_criteria={}\nx_cutoff={}\nransac_iterations={}'.format(calibration_cat, \
    not read_only, make_plot, my_stars_magtype, match_criteria, x_cutoff, ransac_iterations))
if calibration_cat == 'sdss':
    print('sdssmag={}\nmy_gals_magtype={}'.format(sdss_calibration_gals_magtype, my_gals_magtype_sdss))
elif calibration_cat == 'panstarrs':
    print('panstarrs_classification_psc={}\nmy_gals_magtype={}'.format(panstarrs_classification_psc, my_gals_magtype_panstarrs))



if month == 'FebMar':
    filters = ['1sw', '3nw', '1se', '2nw', '3ne', '2ne', '3sw', '2sw', '3se', '2se', '4nw', '4ne', '4sw', '4se', '1ne', '1nw']
# data with Jan
elif month == 'Jan':
    filters = ['1sw', '1ne', '1nw', '3nw', '1se', '2nw', '3ne', '2ne', '3sw', '2sw', '3se', '2se', '4nw', '4ne', '4sw', '4se']

actual_filters = []
filters_TFmask = []
for x in range(len(filters)):
    if os.path.isfile(path + 'cats_' + calibration_cat + '/' + filters[x] + '.cat'):
        actual_filters.append(filters[x])
        filters_TFmask.append(True)
    else:
        filters_TFmask.append(False)

filters = actual_filters


if calibration_cat == 'sdss':
    calibration_filters_original = ['u', 'g', 'r', 'i', 'z']
elif calibration_cat == 'panstarrs':
    calibration_filters_original = ['g', 'r', 'i', 'z', 'y']

calibration_filters = [ calibration_cat + '_' +  calibration_filters_original[x] for x in range(len(calibration_filters_original))]
calibration_extinction_colnames = [ calibration_filters[x] + 'Extinction' for x in range(len(calibration_filters)) ]


# Some useful name lists
calibration_mag_colnames = [calibration_filters[fs] + 'MAG' for fs in range(len(calibration_filters))]  # for filling sdss columns, calibration column names in my table
calibration_err_colnames = [calibration_filters[fs] + 'ERR' for fs in range(len(calibration_filters))]
calibration_catalog_extinction_names = ['extinction_' + calibration_filters_original[fs] for fs in range(len(calibration_filters))]

if calibration_cat == 'sdss':
    calibration_catalog_gal_magnames = [sdss_calibration_gals_magtype + '_' + calibration_filters_original[fs] for fs in range(len(calibration_filters))]
    calibration_catalog_gal_errnames = [sdss_calibration_gals_magtype + 'Err_' + calibration_filters_original[fs] for fs in range(len(calibration_filters))]
    calibration_catalog_star_magnames = ['psfMag_' + calibration_filters_original[fs] for fs in range(len(calibration_filters))]
    calibration_catalog_star_errnames = ['psfMagErr_' + calibration_filters_original[fs] for fs in range(len(calibration_filters))]
elif calibration_cat == 'panstarrs':
    calibration_catalog_gal_magnames = [calibration_filters_original[fs] + 'KronMag' for fs in range(len(calibration_filters))]
    calibration_catalog_gal_errnames = [calibration_filters_original[fs] + 'KronMagErr' for fs in range(len(calibration_filters))]
    calibration_catalog_star_magnames = [calibration_filters_original[fs] + 'MeanPSFMag' for fs in range(len(calibration_filters))]
    calibration_catalog_star_errnames = [calibration_filters_original[fs] + 'MeanPSFMagErr' for fs in range(len(calibration_filters))]



# Create column names
colnames = []
filters_all = filters.copy()
filters_all.extend(calibration_filters)


for x in range(len(filters)):
    colnames.append(filters[x] + 'MAG_AUTO')
    colnames.append(filters[x] + 'ERR_AUTO')
for x in range(len(filters)):
    colnames.append(filters[x] + 'MAG_PSF')
    colnames.append(filters[x] + 'ERR_PSF')
for x in range(len(calibration_filters)):
    colnames.append(calibration_filters[x] + 'MAG')
    colnames.append(calibration_filters[x] + 'ERR')

colnames.insert(0, 'dec')
colnames.insert(0, 'ra')

colnames.extend(['filter_count', 'multiple_matches', 'zspec', 'zspec_err', 'zphoto', 'zphoto_err', 'calibration_quality'])
colnames.extend(calibration_extinction_colnames)


panstarrs_A_EBV = {'panstarrs_g': 3.172, 'panstarrs_r': 2.271, 'panstarrs_i': 1.682, 'panstarrs_z': 1.322, 'panstarrs_y': 1.087}

for x in range(len(filters)):
    filt = filters[x]
    colnames.extend([filt + 'Flags', filt + 'Peak', filt + 'FWHM', filt + 'PSFFWHM',filt + '_chi2'])


if read_only == False:
    print('')
    tabs_Dict = {}
    actual_filters = []
    for x in range(len(filters)):
        tab = ascii.read(path + 'cats_' + calibration_cat + '/' + filters[x] + '.cat')
        c1 = SkyCoord(tab['ALPHA_J2000'], tab['DELTA_J2000'], unit='deg')
        c2 = SkyCoord(tab['ALPHAPEAK_J2000'], tab['DELTAPEAK_J2000'], unit='deg')
        diff = c1.separation(c2)
        k = diff.arcsecond < barypeakdiff
        tabs_Dict[filters[x]] = tab[k]
        del tab, c1, c2, diff, k

    if calibration_cat == 'sdss':
        oldnames = ('ra', 'dec')
        calibration_catalog_filenames = [path + clustername + '_gal_' + calibration_cat + '_radec.csv', path + clustername + '_star_' + calibration_cat + '_radec.csv']
        gals = ascii.read(calibration_catalog_filenames[0])
        stars = ascii.read(calibration_catalog_filenames[1])

    elif calibration_cat == 'panstarrs':
        oldnames = ('raMean', 'decMean')
        calibration_catalog_filenames = path + clustername + '_' + calibration_cat + '_radec.csv'
        panstarrs_tab = ascii.read(calibration_catalog_filenames)
        panstarrs_extinction_tab = ascii.read(path + clustername + '_' + calibration_cat + '_extinction.txt')
        EBV = panstarrs_extinction_tab['E_B_V_SFD']
        panstarrs_tab.add_columns([EBV * panstarrs_A_EBV[calibration_filters[x]] for x in range(len(calibration_filters))], names=calibration_catalog_extinction_names)

        # Convert all -999 entry to -99 for easier table handling in other scripts
        for c in range(len(panstarrs_tab.colnames)):
            pcolname = panstarrs_tab.colnames[c]
            h = panstarrs_tab[pcolname] == -999
            panstarrs_tab[pcolname][h] = -99.0

        psc0 = panstarrs_tab['ps_score']
        gals = panstarrs_tab[psc0 < panstarrs_classification_psc]
        stars = panstarrs_tab[psc0 > panstarrs_classification_psc]
        gals.meta['panstarrs_classification_psc'] = panstarrs_classification_psc
        stars.meta['panstarrs_classification_psc'] = panstarrs_classification_psc
        ascii.write(gals, output=path + clustername + '_gal_' + calibration_cat + '_radec.csv', format='ecsv', overwrite=True)
        ascii.write(stars, output=path + clustername + '_star_' + calibration_cat + '_radec.csv', format='ecsv', overwrite=True)
        # If panstarrs_stars xMeanPSFMag is not available, fill in xPSFMag
        for f in range(len(calibration_filters)):
            fn = calibration_catalog_star_magnames[f]
            en = calibration_catalog_star_errnames[f]
            new_magname = calibration_filters_original[f] + 'PSFMag'
            new_errname = calibration_filters_original[f] + 'PSFMagErr'
            hn = stars[fn] < 0
            stars[fn][hn] = stars[new_magname][hn]
            stars[en][hn] = stars[new_errname][hn]

        del fn, en, hn, panstarrs_tab, pcolname, psc0, new_errname, new_magname, c

    calibration_catalogs = [gals, stars]
    newnames = ('ALPHA_J2000', 'DELTA_J2000')

    # Prepare for multiprocessing    
    object_types = ['galaxy', 'star']
    calibration_catalog_magcolumes = [calibration_catalog_gal_magnames, calibration_catalog_star_magnames]
    calibration_catalog_errcolumes = [calibration_catalog_gal_errnames, calibration_catalog_star_errnames]

    # Define combine catalog function that can create either gal_table or star_table
    def combine_catalog(i, child):

        calibration_catalog = calibration_catalogs[i]
        calibration_catalog.rename_columns(oldnames, newnames)
        tabs_Dict[calibration_cat] = calibration_catalog

        match = FoFCatalogMatching.match(catalog_dict=tabs_Dict, ra_label='ALPHA_J2000', dec_label='DELTA_J2000', linking_lengths=match_criteria)
        if i == 0:
            print('Galaxies FOF matching done')
        elif i == 1:
            print('Stars FOF matching done')

        groupcount = match['group_id'][-1] + 1
        groupcount_temp = groupcount.copy() # Use when calibration objects have duplicate, this number will go up
        initial_array = np.full((groupcount, len(colnames)), -99.0)
        table = Table(initial_array, names=colnames)

        # Change some column values to other values depending on columns
        table['filter_count'][:] = 0
        table['multiple_matches'][:] = 0
        table['zspec'][:] = -1
        table['zspec_err'][:] = -1

        template_row = table[0] # a template row with initialization numbers for adding rows

        object_type = object_types[i]
        calibration_catalog_magcolumn = calibration_catalog_magcolumes[i]
        calibration_catalog_errcolumn = calibration_catalog_errcolumes[i]

        for x in range(groupcount):
            # current object
            obj = match[match['group_id'] == x]
            objs = []   # For consistency in coding, usually objs will just contain obj, in rare cases with duplicated SDSS or Panstarrs, it will contain multiple objects
            # If there are more than one object in this group, and unique catalog numbers is larger than 1, and has an SDSS object match
            if len(obj) > 2 and obj[-1]['catalog_key'] == calibration_cat and len(set(obj['catalog_key'])) > 1:
                # If the total number of objects in this group is larger than unique catalogs, that means some catalog has more than 1 object in this group
                if len(set(obj['catalog_key'])) != len(obj):

                    # If only 1 calibration object, but some filters have duplicate, it means either those filters successfully deblend objects while others don't, or the duplicate filter
                    # detect fake signal while others are correct; either way check separation to sdss object again and throw away anything with separation > 0.5 * match_criteria
                    if len(obj[obj['catalog_key'] == calibration_cat]) == 1:
                        calibration_row = obj[obj['catalog_key'] == calibration_cat]
                        calibration_row_idx = int(calibration_row['row_index'])
                        ra0 = tabs_Dict[calibration_cat][calibration_row_idx]['ALPHA_J2000']
                        dec0 = tabs_Dict[calibration_cat][calibration_row_idx]['DELTA_J2000']
                        c0 = SkyCoord(ra0, dec0, unit='deg')
                        select_mask = np.full(len(obj), True)

                        new_criteria = match_criteria / 2.0
                        
                        for k in range(len(obj) - 1):
                            cat_row = tabs_Dict[obj[k]['catalog_key']][int(obj[k]['row_index'])]
                            ra1 = cat_row['ALPHA_J2000']
                            dec1 = cat_row['DELTA_J2000']
                            c1 = SkyCoord(ra1, dec1, unit='deg')
                            sep = c1.separation(c0).arcsecond
                            if sep > new_criteria:
                                select_mask[k] = False

                        # If the above action still leaves some duplicated filters, choose the closer one in that filter
                        if len(obj[select_mask]) != set(obj[select_mask]['catalog_key']):
                            available_filts = np.array(np.unique(obj['catalog_key'])[:-1])
                            for l in range(len(available_filts)):
                                if len(obj[obj['catalog_key'] == available_filts[l]]) > 1:
                                    duplicated_objs_idxs = np.where(obj['catalog_key'] == available_filts[l])[0]
                                    seps = []
                                    for m in duplicated_objs_idxs:
                                        cat_row2 = tabs_Dict[obj[m]['catalog_key']][int(obj[m]['row_index'])]
                                        ra2 = cat_row2['ALPHA_J2000']
                                        dec2 = cat_row2['DELTA_J2000']
                                        c2 = SkyCoord(ra2, dec2, unit='deg')
                                        sep2 = c2.separation(c0).arcsecond
                                        seps.append(sep2)
                                    seps = np.array(seps)
                                    bad_obj_small_id = np.where(seps != min(seps))[0]
                                    select_mask[duplicated_objs_idxs[bad_obj_small_id]] = False

                        obj = obj[select_mask]
                        table[x]['multiple_matches'] = 1
                        objs.append(obj)

                    # Second case: if calibration catalog has duplicate, and some filters have duplicate, it means there are actually more objects, those duplicated filters correctly pick up
                    # all objects
                    elif len(obj[obj['catalog_key'] == calibration_cat]) > 1 and \
                        len(set(obj[obj['catalog_key'] != calibration_cat]['catalog_key'])) != len(obj[obj['catalog_key'] != calibration_cat]):
                        
                        calibration_rows = obj[obj['catalog_key'] == calibration_cat]
                        calibration_rows_idxs = calibration_rows['row_index']
                        ra0s = []
                        dec0s = []
                        # objs = []
                        # For each duplicated entries in the calibration catalog, create an obj
                        for k in range(len(calibration_rows_idxs)):
                            ra00 = tabs_Dict[calibration_cat][calibration_rows_idxs[k]]['ALPHA_J2000']
                            dec00 = tabs_Dict[calibration_cat][calibration_rows_idxs[k]]['DELTA_J2000']
                            ra0s.append(ra00)
                            dec0s.append(dec00)
                            nfilters = len(filters)
                            obj00 = Table([[0] * nfilters, filters, [x] * nfilters], names=['row_index', 'catalog_key', 'group_id'])
                            obj00.add_row([calibration_rows_idxs[k], calibration_cat,0])
                            objs.append(obj00)

                        c0s = SkyCoord(ra0s, dec0s, unit='deg')
                        available_filts = np.array(np.unique(obj['catalog_key'])[:-1])
                        new_criteria_single = 0.4 * match_criteria
                        new_criteria_multiple = 0.6 * match_criteria

                        for k in range(len(available_filts)):
                            obj1s = obj[obj['catalog_key'] == available_filts[k]]
                            cat1 = tabs_Dict[available_filts[k]]
                            if len(obj1s) > 1:
                                new_criteria = new_criteria_multiple
                            else:
                                new_criteria = new_criteria_single

                            for l in range(len(obj1s)):
                                ra1 = cat1[obj1s[l]['row_index']]['ALPHA_J2000']
                                dec1 = cat1[obj1s[l]['row_index']]['DELTA_J2000']
                                c1 = SkyCoord(ra1, dec1, unit='deg')
                                sep1 = c0s.separation(c1).arcsec
                                
                                while len(sep1[sep1 < new_criteria]) > 1:
                                    new_criteria = new_criteria - 0.02
                                # If the current entry satisfy criteria with respect to exactly one calibration object, write it into the corresponding objs
                                # Find out which calibration object it is
                                if len(sep1[sep1 < new_criteria]) == 1:
                                    objs_id = np.where(sep1 < new_criteria)[0][0]
                                    correct_obj = objs[objs_id]
                                    correct_obj[np.where(correct_obj['catalog_key'] == available_filts[k])[0][0]]['row_index'] = obj1s[l]['row_index']

                        table[x]['multiple_matches'] = 2
                        # Get rid of filters that have no entries
                        objs = [objs[m][objs[m]['row_index'] != 0] for m in range(len(objs))]

                    # If none of above satisfies, skip this loop
                    else:
                        continue

                # If none of above happen, no filters have duplicates, just put it in objs    
                else:
                    objs.append(obj)


                # ==========================================================================================
                # After above steps, the "obj" table should only have 1 object in each available filters, the list "objs" contains 1 or more "obj" table
                # Now write the magnitude and error and sdss's ra/dec to the x row of the final data table

                # First determine how many objects are in objs
                objs_len = len(objs)
                table_idx = [x] # The list of index for filling the table, usually when objs contain only one object, this will just be the number x
                if objs_len > 1:
                    for j in range(objs_len - 1):
                        print('Type: {}\tx={}\tDouble Calibration Objects, add row {}, groupcount={}'.format(object_type, x, groupcount_temp, groupcount))
                        table_idx.append(groupcount_temp)
                        table.add_row(template_row)
                        table[-1]['multiple_matches'] = 2
                        groupcount_temp = groupcount_temp + 1

                for y in range(objs_len):
                    x1 = table_idx[y]
                    obj1 = objs[y]

                    row_idx_calibration = int(obj1[obj1['catalog_key']==calibration_cat]['row_index'])
                    calibration = tabs_Dict[calibration_cat]
                    table[x1]['ra'] = calibration[row_idx_calibration]['ALPHA_J2000']
                    table[x1]['dec'] = calibration[row_idx_calibration]['DELTA_J2000']
                    # Number of detected filters in this object, not including SDSS

                    filter_count = len(obj) - 1
                    table[x1]['filter_count'] = filter_count
                    if calibration_cat == 'sdss':
                        table[x1]['calibration_quality'] = calibration[row_idx_calibration]['clean']
                        table[x1]['zspec'] = calibration[row_idx_calibration]['z']
                        table[x1]['zspec_err'] = calibration[row_idx_calibration]['zerr']
                        table[x1]['zphoto'] = calibration[row_idx_calibration]['photoz']
                        table[x1]['zphoto_err'] = calibration[row_idx_calibration]['photozerr']
                    elif calibration_cat == 'panstarrs':
                        table[x1]['calibration_quality'] = calibration[row_idx_calibration]['qualityFlag']
                        # TEMP need to find zphoto for panstarrs
                        # table[x1]['zspec'] = calibration[row_idx_calibration]['z']
                        # table[x1]['zspec_err'] = calibration[row_idx_calibration]['zerr']
                        # table[x1]['zphoto'] = calibration[row_idx_calibration]['photoz']
                        # table[x1]['zphoto_err'] = calibration[row_idx_calibration]['photozerr']


                    # Write the magnitudes and errors of all filter except sdss into the final data table
                    for j in range(len(obj1)-1):
                        row_idx = obj1[j]['row_index']
                        filtname = obj1[j]['catalog_key']
                        if i == 0:
                            if calibration_cat == 'sdss':
                                # For gal_table, grab the magtype that you choose (MAG_AUTO/MAG_PETRO/MAG_MODEL...), not necessarily MAG_AUTO even though the variable name is mag_auto
                                mag_auto = tabs_Dict[filtname][row_idx][my_gals_magtype_sdss]
                                err_auto = tabs_Dict[filtname][row_idx][my_gals_magtype_sdss.split('_')[0] + 'ERR_' + my_gals_magtype_sdss.split('_')[1]]
                            elif calibration_cat == 'panstarrs':
                                mag_auto = tabs_Dict[filtname][row_idx][my_gals_magtype_panstarrs]
                                err_auto = tabs_Dict[filtname][row_idx][my_gals_magtype_panstarrs.split('_')[0] + 'ERR_' + my_gals_magtype_panstarrs.split('_')[1]]
                        elif i == 1:
                            # For star_table, mag_auto will always be MAG_AUTO
                            mag_auto = tabs_Dict[filtname][row_idx]['MAG_AUTO']
                            err_auto = tabs_Dict[filtname][row_idx]['MAGERR_AUTO']
                        mag_psf = tabs_Dict[filtname][row_idx]['MAG_PSF']
                        err_psf = tabs_Dict[filtname][row_idx]['MAGERR_PSF']
                        flag = tabs_Dict[filtname][row_idx]['FLAGS']
                        peak = tabs_Dict[filtname][row_idx]['FLUX_MAX']
                        chi2 = tabs_Dict[filtname][row_idx]['CHI2_PSF']
                        fwhm = tabs_Dict[filtname][row_idx]['FWHM_IMAGE']
                        psffwhm = tabs_Dict[filtname][row_idx]['FWHMPSF_IMAGE']
                        if mag_auto == 99.0 and err_auto == 99.0:
                            mag_auto = -mag_auto
                            err_auto = -err_auto
                        if mag_psf == 99.0 and err_psf == 99.0:
                            mag_psf = -err_psf
                            err_psf = -err_psf
                        table[x1][filtname + 'MAG_AUTO'] = mag_auto
                        table[x1][filtname + 'ERR_AUTO'] = err_auto
                        table[x1][filtname + 'MAG_PSF'] = mag_psf
                        table[x1][filtname + 'ERR_PSF'] = err_psf
                        table[x1][filtname + 'Flags'] = flag
                        table[x1][filtname + '_chi2'] = chi2
                        table[x1][filtname + 'Peak'] = peak
                        table[x1][filtname + 'FWHM'] = fwhm
                        
                    # Get the sdss/panstarrs magnitudes, errors, and extinction
                    table[x1][calibration_mag_colnames] = calibration[row_idx_calibration][calibration_catalog_magcolumn]
                    table[x1][calibration_err_colnames] = calibration[row_idx_calibration][calibration_catalog_errcolumn]
                    table[x1][calibration_extinction_colnames] = calibration[row_idx_calibration][calibration_catalog_extinction_names]


        # Only keep the rows that have data
        table = table[table['ra'] != -99]

        # Write Metadata to Table
        table.meta['barypeakdiff'] = barypeakdiff
        table.meta['match_criteria'] = match_criteria
        table.meta['peakmax_c'] = peakmax_c
        table.meta['peakmin_c'] = peakmin_c
        table.meta['flag_c'] = flag_c
        table.meta['calibrationmagmax'] = calibrationmagmax
        table.meta['calibrationmagmin'] = calibrationmagmin
        table.meta['fwhmmax_c'] = fwhmmax_c
        table.meta['fwhmmin_c'] = fwhmmin_c
        table.meta['apply_extinction'] = apply_extinction
        table.meta['calibration_cat'] = calibration_cat
        table.meta['month'] = month
        table.meta['calibration_filters'] = calibration_filters
        if calibration_cat == 'panstarrs':
            table.meta['panstarrs_classification_psc'] = panstarrs_classification_psc
        elif calibration_cat == 'sdss':
            table.meta['calibration_gals_type'] = sdss_calibration_gals_magtype
            

        if i == 0:
            print('gal_table original done')
        elif i == 1:
            print('star_table original done')


        child.send(table)
        child.close()
        # return table, matching_log, match
        # End of function

    
    t1 = time.time()

    parent1, child1 = mp.Pipe(duplex=False)
    p1 = mp.get_context('fork').Process(target=combine_catalog, args=(0, child1))
    p1.start()

    parent2, child2 = mp.Pipe(duplex=False)
    p2 = mp.get_context('fork').Process(target=combine_catalog, args=(1, child2))
    p2.start()

    gal_table = parent1.recv()
    star_table = parent2.recv()
    p1.join()
    p2.join()

    t2 = time.time()
    print('\nTotal time to create both tables: {}'.format(t2 - t1))

    # For mags without psf mag, fill it with auto mag
    for f in range(len(filters)):
        filt = filters[f]
        psfmag = star_table[filt + 'MAG_PSF']
        automag = star_table[filt + 'MAG_AUTO']
        h = psfmag == -99
        star_table[filt + 'MAG_PSF'][h] = star_table[filt + 'MAG_AUTO'][h]


    ascii.write(gal_table, output=path+clustername+'_mycatalog_' + calibration_cat + '_gal_original.ecsv', format='ecsv', overwrite=True)
    ascii.write(star_table, output=path+clustername+'_mycatalog_' + calibration_cat + '_star_original.ecsv', format='ecsv', overwrite=True)

    del t1, t2, psfmag, automag, newnames, object_types, oldnames, stars, gals, calibration_catalog_magcolumes, calibration_catalog_filenames, calibration_catalog_errcolumes


else:
    star_table = ascii.read(path+clustername+'_mycatalog_' + calibration_cat + '_star_original.ecsv')
    gal_table = ascii.read(path+clustername+'_mycatalog_' + calibration_cat + '_gal_original.ecsv')


# Apply extinction if set to True
if apply_extinction == True:
    for x in range(len(calibration_filters)):
        star_table[calibration_filters[x] + 'MAG'] = star_table[calibration_filters[x] + 'MAG'] - star_table[calibration_filters[x] + 'Extinction']
        gal_table[calibration_filters[x] + 'MAG'] = gal_table[calibration_filters[x] + 'MAG'] - gal_table[calibration_filters[x] + 'Extinction']


star_table.meta['star_mag_type'] = my_stars_magtype
gal_table.meta['star_mag_type'] = my_stars_magtype

# copy gal and star table that doesn't have zero point calibration but have extinction (if set to True) for ipython
gtable = gal_table.copy()
stable = star_table.copy()




# ======================================================================================================  #
#                                               Calibration                                               #
# ======================================================================================================  #

sg = star_table[calibration_cat + '_gMAG']
dsg = star_table[calibration_cat + '_gERR']
sr = star_table[calibration_cat + '_rMAG']
dsr = star_table[calibration_cat + '_rERR']
si = star_table[calibration_cat + '_iMAG']
dsi = star_table[calibration_cat + '_iERR']
quality = star_table['calibration_quality']


x00 = sg - sr
yri = sr - si
hx = x00 < x_cutoff

hcalibration = (sg < calibrationmagmax) & (sg > calibrationmagmin) & (sr < calibrationmagmax) & (sr > calibrationmagmin)
if calibration_cat == 'sdss':
    hquality = quality == 1
elif calibration_cat == 'panstarrs':
    hquality = quality < 64

# First fit the line of Calibration catalog color-color diagram to make sure b~0, to make sure the x_cutoff value and other criterias are good
x01 = x00[hx]
y01 = yri[hx]
x11 = np.array([])
y11 = np.array([])
cts1 = np.zeros(len(x01))

for i in range(ransac_iterations):
    ransac = linear_model.RANSACRegressor()
    ransac.fit(np.array(x01).reshape(len(x01),1),y01)
    inlier_mask = ransac.inlier_mask_
    x11 = np.append(x11, x01[inlier_mask])
    y11 = np.append(y11, y01[inlier_mask])
    cts1 = cts1 + inlier_mask * 1

mt, bt = np.polyfit(x11, y11, 1)

if make_plot == True:
    fig1 = plt.figure(figsize=(10, 10))
    plt.scatter(x00, yri, marker='.', s=15, c='w', alpha=0.5, edgecolors='k', linewidth=1)
    plt.scatter(x01[cts1>0], y01[cts1>0], marker='.', s=30, c=cts1[cts1>0], edgecolors='gray', linewidth=0.1, cmap=cmap)
    plt.plot(x00, mt * x00 + bt, '-', color=linecolor, linewidth=linewidth0, alpha=linealpha)
    plt.grid()
    plt.xlabel('g - r')
    plt.ylabel('r - i')
    plt.xlim(-1, 3)
    plt.ylim(-1, 3)
    plt.title('g-r cutoff = {}, b_calibration={:.6f}'.format(x_cutoff, bt))

print('\n(r-i) vs (g-r) b = {:.6f}\n'.format(bt))

if make_plot == True:
    fig, axs = plt.subplots(4, 4, figsize=(15, 12), sharex=True, sharey=False)
    fig.tight_layout()
    fig.subplots_adjust(top=0.95, bottom=0.05, left=0.08,right=0.98)


orderx = [3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0]
ordery = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]
if month == 'FebMar':
    filter_lam = np.array(['5100', '5200', '5320', '5400', '5500', '5600', '5680', '5800', '5890', '6000', '6100', '6200', '6300', '6400', '6500', '6600'])
elif month == 'Jan':
    filter_lam = np.array(['4920', '5000', '5100', '5200', '5320', '5400', '5500', '5600', '5680', '5800', '5890', '6000', '6100', '6200', '6300', '6400'])


filter_lam = filter_lam[filters_TFmask]
xline = np.arange(-5,5,0.5)

filt_name_for_record = []
filt_lam_for_record = []
b_value = []
berr_value = []
fit_mags = []

# Dictionary for b correction term from actual SDSS spectra
if apply_b_spectra_correction == True:
    # Using 698 SEDs to calibrate
    b_spectra_correction_Dict_sdss = {'4920': 0.021081905200741817, '5000': -0.015733934540172863, '5100': -0.027419333943118483, '5200': -0.02953593608211138, \
                                '5320': -0.0053906074904824085, '5400': -0.0044880428835484395, '5500': -0.0049806868150285835, '5600': -0.006016068145396759, \
                                '5680': -0.002651410414121113, '5800': -0.0027220452878439997, '5890': -0.04159973136949347, '6000': -0.006519650156443971, \
                                '6100': -0.002927680956932478, '6200': 0.00035885025490195834, '6300': -0.004846132545664436, '6400': -0.010346126433205307, \
                                '6500': 0.010897822416220363, '6600': 0.04194828503459799,}

    b_spectra_correction_errDict_sdss = {'4920': 0.0015531780276028544, '5000': 0.001559402535974323, '5100': 0.001594591617123951, '5200': 0.00202382628478408, \
                                '5320': 0.002014247069868021, '5400': 0.0013678902544186381, '5500': 0.0012828750784658436, '5600': 0.0017947161977452768, \
                                '5680': 0.001843115463545148, '5800': 0.001559804724512924, '5890': 0.0016359537865492761, '6000': 0.0005693145485748402, \
                                '6100': 0.0005403462873010492, '6200': 0.000264381234352183, '6300': 0.0003023198499966118, '6400': 0.004994764239481349, \
                                '6500': 0.0009764994170727277, '6600': 0.0016989610845703188,}

    # Using ~20 SED from my fields
    # b_spectra_correction_Dict_sdss = {'4920': 0.029814306098656578, '5000': -0.0028190004933502286, '5100': -0.0090861137151594, '5200': 0.0014083853337358892, \
    #                             '5320': -3.491933037036188e-05, '5400': 0.0017521731522023949, '5500': -0.005681027029167839, '5600': -0.006912626095401223, \
    #                             '5680': -0.011787154323160288, '5800': -0.010208659392585119, '5890': -0.0035654361308840956, '6000': -0.007100657295703025, \
    #                             '6100': -0.00349335243576131, '6200': 0.0027922959870733885, '6300': 0.0012413967965923163, '6400': 0.0047055432862229225, \
    #                             '6500': 0.01774743228998505, '6600': 0.03892362187444068,}

    # b_spectra_correction_errDict_sdss = {'4920': 0.01441858272811864, '5000': 0.009490856461617382, '5100': 0.009831833775081756, '5200': 0.014944580719957717, \
    #                             '5320': 0.011920805349814928, '5400': 0.012267794293880406, '5500': 0.011352365261750196, '5600': 0.00969051406970445, \
    #                             '5680': 0.009564905096287737, '5800': 0.009359890453288501, '5890': 0.005547851126312819, '6000': 0.007519828700178888, \
    #                             '6100': 0.003344670308370953, '6200': 0.00365510650161222, '6300': 0.0036435704704655706, '6400': 0.005508197748384842, \
    #                             '6500': 0.008158579989755532, '6600': 0.00960722486811168,}

    b_spectra_correction_Dict_panstarrs = {'4920': 0.01562103301544379, '5000': -0.023583350426620985, '5100': -0.02703413492531374, '5200': -0.017907771491671068, \
                                '5320': -0.009142330116112343, '5400': -0.0038189304273664585, '5500': -0.0036810375213502063, '5600': -0.004073235768918441, \
                                '5680': -0.001228688854930752, '5800': -0.0020231022671291813, '5890': -0.03782788045402074, '6000': -0.008894007983504822, \
                                '6100': -0.001218747452394827, '6200': 0.0004789594611687742, '6300': -0.004670494930104791, '6400': -0.004230683795679477, \
                                '6500': 0.017756999880346926, '6600': 0.03711536933003344}

    b_spectra_correction_errDict_panstarrs = {'4920': 0.001235394628705121, '5000': 0.001517678610244854, '5100': 0.0010786743875286198, '5200': 0.0013973577524042647, \
                                '5320': 0.0012891508710228938, '5400': 0.0012006635304456293, '5500': 0.001198466456869297, '5600': 0.001457792714415801, \
                                '5680': 0.0014756939780594279, '5800': 0.0017105403468066548, '5890': 0.0016607389560113536, '6000': 0.0006000345086490612, \
                                '6100': 0.000673288491636023, '6200': 0.00015299898998969035, '6300': 0.0003145336680449157, '6400': 0.00038649067764463357, \
                                '6500': 0.001366910348873381, '6600': 0.002081687538197981}

else:
    b_spectra_correction_Dict_sdss = {'4920': 0.0, '5000': 0.0, '5100': 0.0, '5200': 0.0, \
                                '5320': 0.0, '5400': 0.0, '5500': 0.0, '5600': 0.0, \
                                '5680': 0.0, '5800': 0.0, '5890': 0.0, '6000': 0.0, \
                                '6100': 0.0, '6200': 0.0, '6300': 0.0, '6400': 0.0, \
                                '6500': 0.0, '6600': 0.0,}

    b_spectra_correction_errDict_sdss = {'4920': 0.0, '5000': 0.0, '5100': 0.0, '5200': 0.0, \
                                '5320': 0.0, '5400': 0.0, '5500': 0.0, '5600': 0.0, \
                                '5680': 0.0, '5800': 0.0, '5890': 0.0, '6000': 0.0, \
                                '6100': 0.0, '6200': 0.0, '6300': 0.0, '6400': 0.0, \
                                '6500': 0.0, '6600': 0.0,}

    b_spectra_correction_Dict_panstarrs = {'4920': 0.0, '5000': 0.0, '5100': 0.0, '5200': 0.0, \
                                '5320': 0.0, '5400': 0.0, '5500': 0.0, '5600': 0.0, \
                                '5680': 0.0, '5800': 0.0, '5890': 0.0, '6000': 0.0, \
                                '6100': 0.0, '6200': 0.0, '6300': 0.0, '6400': 0.0, \
                                '6500': 0.0, '6600': 0.0,}

    b_spectra_correction_errDict_panstarrs = {'4920': 0.0, '5000': 0.0, '5100': 0.0, '5200': 0.0, \
                                '5320': 0.0, '5400': 0.0, '5500': 0.0, '5600': 0.0, \
                                '5680': 0.0, '5800': 0.0, '5890': 0.0, '6000': 0.0, \
                                '6100': 0.0, '6200': 0.0, '6300': 0.0, '6400': 0.0, \
                                '6500': 0.0, '6600': 0.0,}


def ransac_iterations_fx(x0, y0, dy, likelihood, child):
    x1 = []
    y1 = []
    cts = np.zeros(len(x0))
    m1s = []
    b1s = []
    b1_errs = []
    mb_corrs = []
    weight = (y0 / dy)** 2
    for i in range(ransac_iterations):
        ransac = linear_model.RANSACRegressor()
        ransac.fit(np.array(x0).reshape(len(x0),1),y0, sample_weight=weight)
        inlier_mask = ransac.inlier_mask_
        x1 = np.append(x1, x0[inlier_mask])
        y1 = np.append(y1, y0[inlier_mask])
        cts = cts + inlier_mask * 1
        (m1, b1), cov1 = np.polyfit(x0[inlier_mask], y0[inlier_mask], 1, cov=True)
        m1s = np.append(m1s, m1)
        b1s = np.append(b1s, b1)
        b1_errs = np.append(b1_errs, cov1[1, 1]** 0.5)
        mb_corrs = np.append(mb_corrs, np.abs(cov1[0,1]))
    child.send([x1, y1, cts, m1s, b1s, b1_errs, mb_corrs])
    del x1, y1, cts, m1s, b1s, b1_errs, i, inlier_mask, cov1, mb_corrs
    child.close()


for f in range(len(filters)):
    filt = filters[f]
    mag_auto = star_table[filt + 'MAG_AUTO']
    dmag_auto = star_table[filt + 'ERR_AUTO']
    mag_psf = star_table[filt + 'MAG_PSF']
    dmag_psf = star_table[filt + 'ERR_PSF']
    peak = star_table[filt + 'Peak']
    fwhm = star_table[filt + 'FWHM']
    chi2 = star_table[filt + '_chi2']
    flag = star_table[filt + 'Flags']

    h0 = mag_psf > -90
    hpeak = (peak > peakmin_c) & (peak < peakmax_c)
    hfwhm = (fwhm > fwhmmin_c) & (fwhm < fwhmmax_c)
    hflag = flag < flag_c

    h = h0 & hx & hcalibration & hquality & hflag

    y00a = mag_auto - sr
    y00p = mag_psf - sr

    x0 = x00[h]
    y0a = y00a[h]
    y0p = y00p[h]
    likelihood = np.exp(-chi2[h])
    dya = np.sqrt(dmag_auto[h] ** 2 + dsr[h] ** 2)
    dyp = np.sqrt(dmag_psf[h] ** 2 + dsr[h] ** 2)

    parent1, child1 = mp.Pipe(duplex=False)
    p1 = mp.get_context('fork').Process(target=ransac_iterations_fx, args=(x0, y0a, dya, likelihood, child1))
    p1.start()

    parent2, child2 = mp.Pipe(duplex=False)
    p2 = mp.get_context('fork').Process(target=ransac_iterations_fx, args=(x0, y0p, dyp, likelihood, child2))
    p2.start()

    x1a, y1a, ctsa, m1sa, b1sa, b1_errsa, mb_corrsa = parent1.recv()
    x1p, y1p, ctsp, m1sp, b1sp, b1_errsp, mb_corrsp = parent2.recv()
    p1.join()
    p2.join()

    m_a = np.mean(m1sa)
    m_p = np.mean(m1sp)
    b_a = np.mean(b1sa)
    b_p = np.mean(b1sp)
    err_b_a = np.std(b1sa, ddof=1) / np.sqrt(ransac_iterations)  # standard deviation of the mean
    err_b_p = np.std(b1sp, ddof=1) / np.sqrt(ransac_iterations)  # standard deviation of the mean


    m_a2, cov_a2 = np.polyfit(x1a, y1a, 1, cov=True)
    m_p2, cov_p2 = np.polyfit(x1p, y1p, 1, cov=True)

    # m = ms[0]
    # b = ms[1]
    # err_b = np.sqrt(cov[1, 1])
    # print(m,b,err_b)
    

    if my_stars_magtype == 'auto':
        fit_mag = 'auto'
        m = m_a
        b = b_a
        err_b = err_b_a
        y00 = y00a
        y0 = y0a
        cts = ctsa
    elif my_stars_magtype == 'psf':
        fit_mag = 'psf'
        m = m_p
        b = b_p
        err_b = err_b_p
        y00 = y00p
        y0 = y0p
        cts = ctsp
    elif my_stars_magtype == 'hybrid':
        # diffa = y0a - (x0 * m_a + b_a)
        # diffp = y0p - (x0 * m_p + b_p)
        if cov_a2[1, 1] >= cov_p2[1, 1]:
            fit_mag = 'psf'
            m = m_p
            b = b_p
            err_b = err_b_p
            y00 = y00p
            y0 = y0p
            cts = ctsp
        else:
            fit_mag = 'auto'
            m = m_a
            b = b_a
            err_b = err_b_a
            y00 = y00a
            y0 = y0a
            cts = ctsa

    # Apply b spectra correction
    if calibration_cat == 'sdss':
        b_spectra_correction_Dict = b_spectra_correction_Dict_sdss
        b_spectra_correction_errDict = b_spectra_correction_errDict_sdss

    elif calibration_cat == 'panstarrs':
        b_spectra_correction_Dict = b_spectra_correction_Dict_panstarrs
        b_spectra_correction_errDict = b_spectra_correction_errDict_panstarrs

        
    b = b - b_spectra_correction_Dict[filter_lam[f]]
    err_b = np.sqrt(err_b ** 2 + b_spectra_correction_errDict[filter_lam[f]]** 2)

    b_value.append(b)
    berr_value.append(err_b)
    filt_name_for_record.append(filt)
    filt_lam_for_record.append(filter_lam[f])
    fit_mags.append(fit_mag)
    if b_spectra_correction_Dict[filter_lam[f]] > 0:
        b_corr_string = '-' + '{:f}'.format(np.abs(b_spectra_correction_Dict[filter_lam[f]]))[:8]
    else:
        b_corr_string = '+' + '{:f}'.format(np.abs(b_spectra_correction_Dict[filter_lam[f]]))[:8]


    print('{} {}  auto: {:.6f} +- {:.6f} | psf: {:.6f} +- {:.6f} ({})  select {}\t{:.6f} {:.6f} {:.6f} {:.6f}'.format(filter_lam[f], \
        filt, b_a, err_b_a, b_p, err_b_p, b_corr_string, fit_mag, np.sqrt(cov_a2[1, 1]), np.sqrt(cov_p2[1, 1]), cov_a2[0, 1] * 1e7, cov_p2[0, 1] * 1e7))


    # Record zero points and err_b
    star_table.meta['b_' + filt] = b
    star_table.meta['err_b_' + filt] = err_b
    gal_table.meta['b_' + filt] = b
    gal_table.meta['err_b_' + filt] = err_b

    h99s_auto = star_table[filt + 'MAG_AUTO'] > -90
    star_table[filt + 'MAG_AUTO'][h99s_auto] = star_table[filt + 'MAG_AUTO'][h99s_auto] - b
    star_table[filt + 'ERR_AUTO'][h99s_auto] = np.sqrt(star_table[filt + 'ERR_AUTO'][h99s_auto] ** 2 + err_b ** 2)
    h99g_auto = gal_table[filt + 'MAG_AUTO'] > -90
    gal_table[filt + 'MAG_AUTO'][h99g_auto] = gal_table[filt + 'MAG_AUTO'][h99g_auto] - b
    gal_table[filt + 'ERR_AUTO'][h99g_auto] = np.sqrt(gal_table[filt + 'ERR_AUTO'][h99g_auto] ** 2 + err_b ** 2)

    h99s_psf = star_table[filt + 'MAG_PSF'] > -90
    star_table[filt + 'MAG_PSF'][h99s_psf] = star_table[filt + 'MAG_PSF'][h99s_psf] - b
    star_table[filt + 'ERR_PSF'][h99s_psf] = np.sqrt(star_table[filt + 'ERR_PSF'][h99s_psf] ** 2 + err_b ** 2)

    if make_plot == True:
        axs[orderx[f], ordery[f]].scatter(x00[h0], y00[h0], c='k', marker='.', s=markersize0, alpha=0.25, edgecolor='w', linewidth=0)
        axs[orderx[f], ordery[f]].scatter(x0[cts>0], y0[cts>0], marker='.', s=markersize1, c=cts[cts>0], cmap=cmap, edgecolor='gray', linewidth=0.1)
        axs[orderx[f], ordery[f]].plot(xline, m * xline + b, '-', color=linecolor, linewidth=linewidth1, label=r'b=${:.6f}\pm{:.6f}$'.format(b, err_b), alpha=linealpha)
        axs[orderx[f], ordery[f]].grid()

        center_x = np.median(x0[cts > 0])
        center_y = np.median(y0[cts > 0])
        ylim_l = center_y - 1.0
        ylim_h = center_y + 1.5

        axs[orderx[f], ordery[f]].set_xlim(xlim_l, xlim_h)
        axs[orderx[f], ordery[f]].set_ylim(ylim_l, ylim_h)
        if fit_mag == 'auto' and my_stars_magtype != 'auto':
            axs[orderx[f], ordery[f]].set_title(filter_lam[f] + r'$\AA$' + ' (' + fit_mag + ')', fontsize=small_title_fontsize, pad=2)
        else:
            axs[orderx[f], ordery[f]].set_title(filter_lam[f] + r'$\AA$', fontsize=small_title_fontsize, pad=2)
        axs[orderx[f], ordery[f]].legend(loc=1, fontsize='x-small')
        axs[orderx[f], ordery[f]].tick_params(labelsize=ticksize)


    # Delete unnecessary variables in the loop
    del filt, mag_auto, dmag_auto, mag_psf, dmag_psf, peak, fwhm, chi2, flag, h0, hpeak, hfwhm, hflag, h, \
        y00, x0, y0, likelihood, x1a, x1p, y1a, y1p, y00a, y00p, y0a, y0p, ctsa, ctsp, m1sa, m1sp, b1sa, b1sp, \
        b1_errsa, b1_errsp, m, b, err_b, h99g_auto, h99s_auto, h99s_psf, m_a, m_p, b_a, b_p, err_b_a, err_b_p, \
        cov_a2, cov_p2, cts, dya, dyp, m_a2, m_p2, mb_corrsa, mb_corrsp#, diffa, diffp
    


# Drop psf mags from gal_table, drop unused mag type from star_table
gal_final_cols = [gal_table.colnames[x] for x in range(len(gal_table.colnames)) if gal_table.colnames[x].split('_')[-1] != 'PSF']
gal_table = gal_table[gal_final_cols]

star_final_cols_auto = [star_table.colnames[x] for x in range(len(star_table.colnames)) if star_table.colnames[x].split('_')[-1] != 'PSF']
star_table_auto = star_table[star_final_cols_auto]
star_final_cols_psf = [star_table.colnames[x] for x in range(len(star_table.colnames)) if star_table.colnames[x].split('_')[-1] != 'AUTO']
star_table_psf = star_table[star_final_cols_psf]

# Convert magnitude column names, get rid of "_AUTO" or "_PSF"
gal_mag_auto_cols = [gal_table.colnames[x] for x in range(len(gal_table.colnames)) if gal_table.colnames[x].split('_')[-1]=='AUTO']
gal_mag_auto_cols_new = [gal_mag_auto_cols[x].split('_')[0] for x in range(len(gal_mag_auto_cols))]
gal_table.rename_columns(gal_mag_auto_cols, gal_mag_auto_cols_new)

star_mag_cols_auto = [star_table_auto.colnames[x] for x in range(len(star_table_auto.colnames)) if star_table_auto.colnames[x].split('_')[-1] == 'AUTO']
star_mag_cols_new_auto = [star_mag_cols_auto[x].split('_')[0] for x in range(len(star_mag_cols_auto))]
star_table_auto.rename_columns(star_mag_cols_auto, star_mag_cols_new_auto)

star_mag_cols_psf = [star_table_psf.colnames[x] for x in range(len(star_table_psf.colnames)) if star_table_psf.colnames[x].split('_')[-1] == 'PSF']
star_mag_cols_new_psf = [star_mag_cols_psf[x].split('_')[0] for x in range(len(star_mag_cols_psf))]
star_table_psf.rename_columns(star_mag_cols_psf, star_mag_cols_new_psf)

if calibration_cat == 'sdss':
    calibration_cat_propername = 'SDSS'
elif calibration_cat == 'panstarrs':
    calibration_cat_propername = 'Pan-STARRS1'

if make_plot == True:
    fig.suptitle(clustername + ', Calibrated with ' + calibration_cat_propername + ', ' + my_stars_magtype.lower(), fontsize=titlesize)
    # fig.suptitle('ABELL 611 (Calibrated with ' + calibration_cat_propername + ')', fontsize=titlesize)# + ', ' + my_stars_magtype.lower())    # used for thesis plots
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.xlabel('g - r', labelpad=0, fontsize=xlabelsize)
    plt.ylabel('narrowband - r', labelpad=5, fontsize=ylabelsize)



# If using PanSTARRS, read SDSS to find z_specs if SDSS is available
if calibration_cat == 'panstarrs' and os.path.isfile(path + clustername + '_gal_sdss_radec.csv'):
    print('\nFill in SDSS zspecs')
    sdss_gtab = ascii.read(path + clustername + '_gal_sdss_radec.csv')
    sdss_zspec = sdss_gtab['z']
    h = sdss_zspec > 0
    sdss_zspec = sdss_zspec[h]
    sdss_zspecerr = sdss_gtab['zerr'][h]
    if len(sdss_zspec) > 0:
        sdss_ra = sdss_gtab['ra'][h]
        sdss_dec = sdss_gtab['dec'][h]
        c_zspec = SkyCoord(sdss_ra, sdss_dec, unit='deg')
        c_my = SkyCoord(gal_table['ra'], gal_table['dec'], unit='deg')
        idx, d2d, d3d = c_my.match_to_catalog_sky(c_zspec)
        h = d2d.arcsec < match_criteria
        gal_table['zspec'][h] = sdss_zspec[idx[h]]
        gal_table['zspec_err'][h] = sdss_zspecerr[idx[h]]

        del sdss_gtab, sdss_zspec, h, sdss_zspecerr, sdss_ra, sdss_dec, c_zspec, c_my, idx, d2d, d3d


ascii.write(gal_table, output=path + clustername + '_mycatalog_' + calibration_cat + '_gal.csv', format='ecsv', overwrite=True)
ascii.write(star_table_auto, output=path + clustername + '_mycatalog_' + calibration_cat + '_star_auto.csv', format='ecsv', overwrite=True)
ascii.write(star_table_psf, output=path + clustername + '_mycatalog_' + calibration_cat + '_star_psf.csv', format='ecsv', overwrite=True)


ascii.write([gal_table['ra'],gal_table['dec'],np.arange(len(gal_table['ra']))], output=path + clustername + '_mycatalog_' + calibration_cat + '_gal_ds9.csv', format='ecsv', overwrite=True)
ascii.write([star_table['ra'],star_table['dec'],np.arange(len(star_table['ra']))], output=path + clustername + '_mycatalog_' + calibration_cat + '_star_ds9.csv', format='ecsv', overwrite=True)


hspec_gal = (gal_table['zspec'] != -1) &(gal_table['zspec_err'] != -1)
hspec_star = (star_table['zspec'] != -1) &(star_table['zspec_err'] != -1)

if calibration_cat == 'sdss':
    ascii.write(gal_table['ra', 'dec'][hspec_gal], output=path+clustername+'_mycatalog_gal_w_' + calibration_cat + '_specs_radec.txt', format='no_header', delimiter=',', overwrite=True)
    ascii.write(star_table['ra', 'dec'][hspec_star], output=path+clustername+'_mycatalog_star_w_' + calibration_cat + '_specs_radec.txt', format='no_header', delimiter=',', overwrite=True)

# Write all b value and error to a table for future reference
b_table = Table([filt_name_for_record, filt_lam_for_record, b_value, berr_value, fit_mags], names=['name', 'wavelength', 'b', 'err_b', 'fit_mag'])
ascii.write(b_table, output=path + 'b_values_' + calibration_cat + '.ecsv', format='ecsv', overwrite=True)

plt.show()





# ======================================================================================================  #
#                                          Generates EAZY Format                                          #
# ======================================================================================================  #


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
eazy_filts = [filt_Dict_all[filters_all[x] + 'MAG'] for x in range(len(filters_all))]
hdrs = []
for x in range(len(eazy_filts)):
    hdrs = hdrs + [eazy_filts[x]] + [eazy_filts[x].replace('F', 'E')]
hdrs2 = hdrs.copy()
hdrs.insert(0, 'id')
hdrs.insert(1, 'z_spec')

mag_err_col_start = np.where(np.array(gal_table.colnames) == 'dec')[0][0] + 1
mag_err_col_end = np.where(np.array(gal_table.colnames) == 'filter_count')[0][0]
mag_err_colnames = gal_table.colnames[mag_err_col_start:mag_err_col_end]
gal_eazy = gal_table[mag_err_colnames]

# Add index column for EAZY
gal_eazy.add_column(np.arange(len(gal_eazy)), name='idx', index=0)
# Add z_spec column
gal_eazy.add_column(gal_table['zspec'], name='zspec', index=1)

# Add ra, dec, photoz, photoz_err from gal_table
gal_eazy.add_columns([gal_table['ra'], gal_table['dec'], gal_table['zphoto'], gal_table['zphoto_err']])

gal_eazy.rename_columns(gal_eazy.colnames[:-4], hdrs)

ascii.write(gal_eazy, output=path + clustername + '_mycatalog_' + calibration_cat + '_gal_eazy.cat', delimiter='\t', format='commented_header', overwrite=True)



# ======================================================================================================  #
#                                           Merges with ALLWISE                                           #
# ======================================================================================================  #




zspec = gal_table['zspec']
coords = SkyCoord(gal_table['ra'], gal_table['dec'], unit='deg')
magstart = np.where(np.array(gal_table.colnames) == 'dec')[0][0] + 1
magerrend = np.where(np.array(gal_table.colnames) == 'filter_count')[0][0]
mags = gal_table[gal_table.colnames[magstart:magerrend]]
bs= []
berrs = []
for i in range(len(filters)):
    bs.append(gal_table.meta['b_' + filters[i]])
    berrs.append(gal_table.meta['err_b_' + filters[i]])


# Read ALLWISE
t = ascii.read(path + clustername + '_allwise_radec.txt')
wise_A_EBV = {'WISE1': 0.189, 'WISE2': 0.146}
allwise_extinction_tab = ascii.read(path + clustername + '_allwise_extinction.txt')
aEBV = allwise_extinction_tab['E_B_V_SFD']
w1_extinction = aEBV * wise_A_EBV['WISE1']
w2_extinction = aEBV * wise_A_EBV['WISE2']
callwise = SkyCoord(t['ra'], t['dec'], unit='deg')
t.rename_columns(('w1mpro', 'w1sigmpro', 'w2mpro', 'w2sigmpro'), ('allwise_w1MAG', 'allwise_w1ERR', 'allwise_w2MAG', 'allwise_w2ERR'))
t['allwise_w1MAG'] = t['allwise_w1MAG'] + allwise_ab_vega[0]
t['allwise_w2MAG'] = t['allwise_w2MAG'] + allwise_ab_vega[1]
if apply_extinction == True:
    t['allwise_w1MAG'] = t['allwise_w1MAG'] - w1_extinction
    t['allwise_w2MAG'] = t['allwise_w2MAG'] - w2_extinction

allwise = t['allwise_w1MAG', 'allwise_w1ERR', 'allwise_w2MAG', 'allwise_w2ERR']
for x in range(len(allwise.colnames)):
    colname_temp = allwise.colnames[x]
    # hnull = allwise[colname_temp] == 'null'
    # print(hnull)
    for xi in range(len(allwise[colname_temp])):
        if allwise[colname_temp][xi] == 'null':
            allwise[colname_temp][xi] = null_fillval

filters_all.extend(allwise.colnames[::2])
hdrs2.extend([filt_Dict_all[allwise.colnames[0]], filt_Dict_all[allwise.colnames[0]].replace('F','E'), filt_Dict_all[allwise.colnames[2]], filt_Dict_all[allwise.colnames[2]].replace('F','E')])

idx, d2d, d3d = coords.match_to_catalog_sky(callwise)
ha = d2d.arcsec < match_criteria

nobj1 = len(mags)
tab = Table(np.full((nobj1, len(hdrs2)), -99.0), names=hdrs2)
tab[ha] = hstack([mags[ha], allwise[idx[ha]]])
tab[~ha] = hstack([mags[~ha], Table(np.full((len(mags[~ha]), 4), -99.0))])
tab.add_columns([np.arange(len(tab)), zspec], names=['id', 'z_spec'], indexes=[0, 0])
tab.add_columns([gal_table['ra'], gal_table['dec']], names=['ra', 'dec'], indexes=[0, 0])
ascii.write(tab, output=path + clustername + '_merged_' + calibration_cat + '_allwise_eazy.cat', format='commented_header', delimiter='\t', overwrite=True)




# # Delete unnecessary variables
# del x01, y01, x11, y11, cts1, mt, bt, orderx, ordery, f, xline, gal_final_cols, \
#     star_final_cols_auto, star_final_cols_psf, gal_mag_auto_cols, gal_mag_auto_cols_new, star_mag_cols_auto, star_mag_cols_psf, star_mag_cols_new_auto, star_mag_cols_new_psf, \
#     hspec_gal, hspec_star, apply_extinction, barypeakdiff, calibration_cat, calibration_catalog_extinction_names, \
#     calibration_catalog_gal_errnames, calibration_catalog_gal_magnames, calibration_catalog_star_errnames, calibration_catalog_star_magnames, \
#     calibration_err_colnames, calibration_extinction_colnames, calibration_files, calibration_filters, calibration_filters_original, sdss_calibration_gals_magtype, \
#     calibration_mag_colnames, calibrationmagmax, calibrationmagmin, cmap, colnames, flag_c, fwhmmax_c, fwhmmin_c, \
#     hcalibration, linewidth0, linewidth1, linewidth2, markersize0, markersize1, markersize2, match_criteria, path, peakmax_c, peakmin_c, \
#     ransac_iterations, read_only, small_title_fontsize, x_cutoff, x, x00, xlim_h, xlim_l, ylim_h, ylim_l, my_stars_magtype, hx, \
#     b_value, berr_value, filt_name_for_record, filt_lam_for_record, fit_mag, fit_mags, \
#     linealpha, linecolor, hquality, yri, p1, p2, parent1, parent2, child1, child2, i, dsg, dsi, dsr, month, panstarrs_classification_psc, quality, calibration_cat_propername, \
#     mag_err_col_start, mag_err_col_end, mag_err_colnames, gal_eazy

# if make_plot == True:
#     del axs, fig, fig1
# del make_plot

