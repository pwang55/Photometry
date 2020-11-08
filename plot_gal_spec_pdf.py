import numpy as np
import sys
from astropy.io import ascii,fits
from astropy.table import Table, Column, join
from astropy import units as u
from astropy.coordinates import match_coordinates_sky, SkyCoord
import subprocess
import pandas as pd
from glob import glob
import scipy.optimize as opt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from scipy.integrate import simps

# plt.ion()

path_file = sys.argv[1] # The path to the mycatalog gal file
filename = path_file.split('/')[-1]
clustername = filename.split('_')[0]
path = path_file[:-len(filename)]
calibration_cat = filename.split('_')[2]
if calibration_cat != 'sdss' and calibration_cat != 'panstarrs':
    print(calibration_cat)
    print('Incorrect calibration catalog type, check filename!')
    sys.exit()

cat = ascii.read(path_file)
# month = cat.meta['month']
# TEMP
month = 'FebMar'

script_dir = sys.path[0]
filterpath = '/Users/brianwang76/sources/90Prime/Filters/'

specs = np.array(glob(path + 'spec-*fits'))
# zout = ascii.read(path + clustername + '_' + calibration_cat + '.zout')
zout = ascii.read(path_file.replace('_eazy.cat', '.zout'))

c = 3.0 * 10 ** 18
scale0 = 10 ** 17

pdf = PdfPages(path + clustername + '_' + calibration_cat + '_flux_compare_spectra.pdf')


# function for fitting the arbritrary scale
def scalefunc(a, x):
    return a * x

# function that convert mag to f_lambda
def mag2flux(l,mag):
    return 10 ** ((mag + 48.6) / (-2.5)) * c / l ** 2
    



# Set up table for filters information
filt_name = ['F348', 'F349', 'F350', 'F351', 'F352', 'F353', 'F354', 'F355', 'F356', 'F357', 'F358', 'F359', 'F360', 'F361', 'F362', 'F363', 'F364', 'F365', \
        'F156', 'F157', 'F158', 'F159', 'F160', \
        'F244', 'F245', \
        'F334', 'F335', 'F336', 'F337', 'F338', \
        'F113', 'F114', 'F115', 'F285', 'F283', 'F117', 'F118', 'F222']

if month == 'FebMar':
    # Feb & Mar dictionary used for translating filters to eazy filter names
    filt_Dict = {'1swMAG': 'F350', '3nwMAG': 'F351', '1seMAG': 'F352', '2nwMAG': 'F353', \
        '3neMAG': 'F354', '2neMAG': 'F355', '3swMAG': 'F356', '2swMAG': 'F357', \
        '3seMAG': 'F358', '2seMAG': 'F359', '4nwMAG': 'F360', '4neMAG': 'F361', \
        '4swMAG': 'F362', '4seMAG': 'F363', '1neMAG': 'F364', '1nwMAG': 'F365', \
        'sdss_uMAG': 'F156', 'sdss_gMAG': 'F157', 'sdss_rMAG': 'F158', 'sdss_iMAG': 'F159', 'sdss_zMAG': 'F160', \
        'allwise_w1': 'F244', 'allwise_w2': 'F245', \
        'panstarrs_gMAG': 'F334', 'panstarrs_rMAG': 'F335', 'panstarrs_iMAG': 'F336', 'panstarrs_zMAG': 'F337', 'panstarrs_yMAG': 'F338'}
elif month == 'Jan':
    # Jan dictionary used for translating filters to eazy filter names
    filt_Dict = {'1swMAG': 'F348', '1neMAG': 'F349', '1nwMAG': 'F350', '3nwMAG': 'F351',
        '1seMAG': 'F352', '2nwMAG': 'F353', '3neMAG': 'F354', '2neMAG': 'F355',
        '3swMAG': 'F356', '2swMAG': 'F357', '3seMAG': 'F358', '2seMAG': 'F359',
        '4nwMAG': 'F360', '4neMAG': 'F361', '4swMAG': 'F362', '4seMAG': 'F363', \
        'sdss_uMAG': 'F156', 'sdss_gMAG': 'F157', 'sdss_rMAG': 'F158', 'sdss_iMAG': 'F159', 'sdss_zMAG': 'F160', \
        'allwise_w1': 'F244', 'allwise_w2': 'F245', \
        'panstarrs_gMAG': 'F334', 'panstarrs_rMAG': 'F335', 'panstarrs_iMAG': 'F336', 'panstarrs_zMAG': 'F337', 'panstarrs_yMAG': 'F338'}

filt_type = ['narrowband', 'narrowband', 'narrowband', 'narrowband', 'narrowband', 'narrowband', 'narrowband', 'narrowband', 'narrowband', 'narrowband', 'narrowband', 'narrowband', 'narrowband', 'narrowband', 'narrowband', 'narrowband', 'narrowband', 'narrowband', \
            'sdss', 'sdss', 'sdss', 'sdss', 'sdss', \
            'allwise', 'allwise', \
            'panstarrs', 'panstarrs', 'panstarrs', 'panstarrs', 'panstarrs', \
            'clash','clash','clash','clash','clash','clash','clash','clash']

filt_lam = [4920, 5000, 5100, 5200, 5320, 5400, 5500, 5600, 5680, 5800, 5890, 6000, 6100, 6200, 6300, 6400, 6500, 6600, \
        3543, 4770, 6231, 7625, 9134, \
        33682, 46179, \
        4849.11, 6201.19, 7534.96, 8674.18, 9627.77, \
        3817.2, 4448, 5470.2, 6505.4, 7960.5, 7671.2, 9028.2, 21574]


filt_table = Table()
filt_table.add_columns([filt_name, filt_lam, filt_type], names=['name', 'lam', 'type'])


# # Find where magnitudes end
# imagerr_end = np.where(np.array(cat.colnames) == 'filter_count')[0][0]
# fnames_original = cat.colnames[2:imagerr_end:2]
# # Translate original filter names from gal_table to eazy filter names using dictionary
# fnames = [filt_Dict[fnames_original[x]] for x in range(len(fnames_original))]
# hfilt = np.in1d(filt_table['name'], fnames)
# filt_table = filt_table[hfilt]


fnames = cat.colnames[4:-4:2]
# print(fnames)
hfilt = np.in1d(filt_table['name'], fnames)
filt_table = filt_table[hfilt]
# print(filt_table)

# Read FILTER.RES.latest.info anc convert info into index and number of lines in each filter
f = open(filterpath + 'FILTER.RES.latest.info')
lines = []
fidx = []
fnl = []
for l in f:
    lines.append(l.strip().split())
for x in range(len(lines)):
    fidx.append(int(lines[x][0]))
    fnl.append(int(lines[x][1]))
fidx = np.array(fidx)
fnl = np.array(fnl)
f.close()

# Now read the FILTER.RES.latest in order to get the filter curves of filters that are in filt_table[hfilt], save them into a dictionary filt_profiles
f = open(filterpath + 'FILTER.RES.latest')
lines = []
filt_profiles = {}
for l in f:
    lines.append(l.strip().split())
# for nth filter counting from 0, do: for x in range(sum(fnl[:n1])+fidx[n], sum(fnl[:n+1])+fidx[n]): filt.append(lines[x]), then Table(np.array(filt))
for x in range(len(filt_table['name'])):
    name = filt_table['name'][x]
    idx = int(name[1:])
    temp_filt = []
    n = np.where(fidx == idx)[0][0] # nth filter counting from zero
    for y in range(sum(fnl[:n]) + fidx[n], sum(fnl[:n + 1]) + fidx[n]):
        temp_filt.append(lines[y])
    filt_profiles[name] = Table(np.array(temp_filt), dtype=[int, float, float], names=['idx', 'lam', 'Response'])
f.close()




ra = []
dec = []
for x in range(len(specs)):
    r = subprocess.run(['gethead', 'plug_ra', 'plug_dec', path+specs[x]], capture_output=True, text=True)
    out = r.stdout.strip()
    ra.append(out.split(' ')[0])
    dec.append(out.split(' ')[1])

cspec = SkyCoord(ra, dec, unit='deg')



h = cat['z_spec'] != -1

cmy = SkyCoord(cat[h]['ra'], cat[h]['dec'], unit='deg')


idx, d2d, d3d = cmy.match_to_catalog_sky(cspec)
d2d = d2d.arcsecond

hd = d2d < 1.0

# my catalog that has zspec and the matched, corresponding spectra
cmy = cmy[hd]
cspec = cspec[idx[hd]]
specs = specs[idx[hd]]

my_mags = cat[h][hd][cat.colnames[4:-4:2]]
my_errs = cat[h][hd][cat.colnames[5:-4:2]]

zout_zspec = zout[h][hd]['z_spec']
zout_zphot = zout[h][hd]['z_m2']
zout_id = zout[h][hd]['id']
zout_1sigma = (zout[h][hd]['u68'] - zout[h][hd]['l68'])/2

ras = cat[h][hd]['ra']
decs = cat[h][hd]['dec']

# Loop over each spectra to make plots
for x in range(len(cmy)):

    if zout_zphot[x] > 0:
        spec = fits.open(specs[x])
        loglam = spec[1].data['loglam']
        lam = 10 ** loglam
        flux = spec[1].data['flux']
        best = spec[1].data['model']

        my_mag = my_mags[x]
        my_err = my_errs[x]

        # Loop over each magnitude to record them for comparison later
        final_lam = []
        final_mag_my = []
        final_magerr_my = []
        final_flux_my = []
        final_fluxerr_my = []
        final_type_my = []
        final_mag_spec = []
        final_flux_spec = []
        # TEMP other types of magnitude

        for y in range(len(my_mag)):
            if my_mag[y] > 0:
                filt_name_now = filt_table[y]['name']
                final_mag_my.append(my_mag[y])
                final_lam.append(filt_table[y]['lam'])
                final_magerr_my.append(my_err[y])
                final_type_my.append(filt_table[y]['type'])

                filt_profile = filt_profiles[filt_name_now]
                xl = filt_profile['lam']
                response = filt_profile['Response']
                # flux_interp = np.interp(xl, lam, best)
                flux_interp = np.interp(xl, lam, flux)

                ave_flux = simps(flux_interp * response, xl) / simps(response, xl)  # average flux from integration of spectra and response curve; theoretical value
                final_flux_spec.append(ave_flux)
                obd_flux = mag2flux(filt_table[y]['lam'], my_mag[y])                # Observed flux converted from magnitude
                final_flux_my.append(obd_flux)
                obd_flux_err = 0.4 * obd_flux * my_err[y]
                final_fluxerr_my.append(obd_flux_err)

        final_flux_my = np.array(final_flux_my)# * scale0
        final_fluxerr_my = np.array(final_fluxerr_my)# * scale0
        final_lam = np.array(final_lam)
        final_flux_spec = np.array(final_flux_spec) / scale0
        final_type_my = np.array(final_type_my)

        # print(final_lam)

        # h = (final_lam < max(lam)) & (final_lam > min(lam))
        h = (final_lam < max(lam)) & (final_lam > min(lam)) & (final_type_my == calibration_cat)
        final_scale, final_scale_err = opt.curve_fit(scalefunc, final_flux_my[h], final_flux_spec[h], sigma=final_fluxerr_my[h])
        # final_flux_my = final_flux_my * final_scale[0]
        # final_fluxerr_my = final_fluxerr_my * final_scale[0]
        final_flux_spec = final_flux_spec / final_scale[0]
        best = best / final_scale[0] / scale0
        flux = flux / final_scale[0] / scale0

        # if min(best) < 0:
        #     ybot = 1.1 * min(flux)
        # else:
        #     ybot = -0.5 / scale0
        ybot = -0.5 / scale0
        ytop = 1.85 * max(final_flux_my)
        
        fig = plt.figure(figsize=(15,10))
        # fig = plt.figure()
        plt.plot(lam, best, 'c--', alpha=0.3, label='SDSS spectra model fit')
        plt.plot(lam, flux, 'b:', alpha=0.1, label='SDSS spectra')

        if calibration_cat == 'sdss':
            calibration_label = 'SDSS Flux'
            calibration_flux_label = 'Flux from Spectra with SDSS filters'
        elif calibration_cat == 'panstarrs':
            calibration_label = 'PanSTARRS Flux'
            calibration_flux_label = 'Flux from Spectra with PanSTARRS filters'
        plt.errorbar(final_lam[final_type_my == 'narrowband'], final_flux_my[final_type_my == 'narrowband'], \
            final_fluxerr_my[final_type_my == 'narrowband'], fmt='ro', label='Narrowband Observed Flux')
        plt.errorbar(final_lam[final_type_my == calibration_cat], final_flux_my[final_type_my == calibration_cat], \
            final_fluxerr_my[final_type_my == calibration_cat], fmt='go', label=calibration_label)

        plt.errorbar(final_lam[final_type_my == 'narrowband'], final_flux_spec[final_type_my == 'narrowband'], fmt='bx', label='Flux from Spectra with Narrowband Filters')
        plt.errorbar(final_lam[final_type_my == calibration_cat], final_flux_spec[final_type_my == calibration_cat], fmt='kx', label=calibration_flux_label)

        plt.legend(loc=1, prop={'size': 10})
        plt.grid()
        plt.xlim(3000, 10000)
        plt.ylim(ybot, ytop)
        plt.xlabel(r'$\lambda  [\AA]$', fontsize=12)
        plt.ylabel(r'$f_\lambda  [erg/s/cm^2/\AA]$', fontsize=15)
        plt.title(r'$id: {}$  $z_s={}$,  $z_p={:.4f}\pm{:.4f}$,  $\sigma_z/(1+z_p)={:.4f}\%$,  $(z_s-z_p)/(1+z_s)={:.4f}\%$'.format(zout_id[x], zout_zspec[x], zout_zphot[x], \
            zout_1sigma[x], 100 * zout_1sigma[x] / (1 + zout_zphot[x]), 100 * (zout_zphot[x] - zout_zspec[x]) / (1 + zout_zspec[x])), fontsize=12)
        pdf.savefig(fig)
        plt.close()

pdf.close()
