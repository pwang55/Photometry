import numpy as np
from astropy.io import fits, ascii
from glob import glob
from scipy import integrate
import sys
from astropy.table import Table
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
import matplotlib.pyplot as plt

plt.ion()


path_file = sys.argv[1]     # a text list containing SEDs you want to use for corrections
filterpath = '/Users/brianwang76/sources/90Prime/Filters/'
int_type = 'simps'  # simps or trapz
bestTF = True   # use best fit SED or original data

h = 6.626 * 10 ** (-34)
c = 299792458 * 10 ** (10)


# The query used to get stars with spectrum
# SELECT p.ra, p.dec, p.psfMag_g, p.psfMag_r, p.psfMag_i, s.z, s. elodieTEff, s.subclass, s.rchi2
# FROM fGetNearbyObjEq(120, 36, 100) n, Star p 
# JOIN specObj s ON s.bestObjID=p.objID
# WHERE n.objID=p.objID AND (s.class = "STAR")
# AND p.psfMag_g-p.psfMag_r < 1.2
# AND p.psfMag_g-p.psfMag_r > -0.5




files = []
with open(path_file, 'r') as pf:
    for l in pf:
        files.append(l.strip())

narrowbands=['F348', 'F349', 'F350', 'F351', 'F352', 'F353', 'F354', 'F355', 'F356', 'F357', 'F358', 'F359', 'F360', 'F361', 'F362', 'F363', 'F364', 'F365']
sdss_g = 'F157'
sdss_r = 'F158'
panstarrs_g = 'F334'
panstarrs_r = 'F335'


filt_name = ['F348', 'F349', 'F350', 'F351', 'F352', 'F353', 'F354', 'F355', 'F356', 'F357', 'F358', 'F359', 'F360', 'F361', 'F362', 'F363', 'F364', 'F365', \
            'F156', 'F157', 'F158', 'F159', 'F160', 'F334', 'F335', 'F336', 'F337', 'F338']

filt_type = ['narrowband', 'narrowband', 'narrowband', 'narrowband', 'narrowband', 'narrowband', 'narrowband', 'narrowband', 'narrowband', 'narrowband', 'narrowband', 'narrowband', 'narrowband', 'narrowband', 'narrowband', 'narrowband', 'narrowband', 'narrowband', \
            'sdss', 'sdss', 'sdss', 'sdss', 'sdss', 'panstarrs', 'panstarrs', 'panstarrs', 'panstarrs', 'panstarrs']

filt_lam = [4920, 5000, 5100, 5200, 5320, 5400, 5500, 5600, 5680, 5800, 5890, 6000, 6100, 6200, 6300, 6400, 6500, 6600, \
            3543, 4770, 6231, 7625, 9134, 4849.11, 6201.19, 7534.96, 8674.18, 9627.77]

filt_table = Table()
filt_table.add_columns([filt_name, filt_lam, filt_type], names=['name', 'lam', 'type'])


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


def mag(fname, filtname, bestTF):
    f = fits.open(fname)
    loglam = f[1].data['loglam']
    lam = 10 ** loglam
    flux = f[1].data['flux'] * lam ** 2 / c
    best = f[1].data['model'] * lam ** 2 / c
    filt = filt_profiles[filtname]
    filtlam = filt['lam']
    response = filt['Response']
    resamp_r = np.interp(lam, filtlam, response)
    if int_type == 'trapz':
        integral_ubest = np.trapz(best * resamp_r * 10 ** (-17) * lam / h / c)
        integral_uflux = np.trapz(flux * resamp_r * 10 ** (-17) * lam / h / c)
    elif int_type == 'simps':
        integral_ubest = integrate.simps(best * resamp_r * 10 ** (-17) * lam / h / c)
        integral_uflux = integrate.simps(flux * resamp_r * 10 ** (-17) * lam / h / c)
    integral_l = np.trapz(10 ** (-23) * 3631 * resamp_r * lam / h / c)
    if bestTF == True:
        return - 2.5 * np.log10(integral_ubest / integral_l)
    else:
        return - 2.5 * np.log10(integral_uflux / integral_l)

bs = []
berrs = []

# TEMP
# for f in range(len(narrowbands)):
#     x0=[]
#     y0=[]
#     for x in range(len(files)):
#         gmag = mag(files[x], sdss_g, True)
#         rmag = mag(files[x], sdss_r, True)
#         x0.append(gmag - rmag)
#         nmag = mag(files[x], narrowbands[f], True)
#         y0.append(nmag - rmag)

#     (m, b), cov = np.polyfit(x0, y0, 1, cov=True)
#     bs.append(b)
#     berrs.append(cov[1, 1]** 0.5)
#     print(filt_lam[f], narrowbands[f], '\t', b, '\t', cov[1, 1]** 0.5)


x0=[]
y0 = []
f = 0
best_TF = True
for x in range(len(files)):
    gmag = mag(files[x], sdss_g, best_TF)
    rmag = mag(files[x], sdss_r, best_TF)
    x0.append(gmag - rmag)
    nmag = mag(files[x], narrowbands[f], best_TF)
    y0.append(nmag - rmag)

fig = plt.figure(figsize=(20, 13))

plt.plot(x0, y0, '.')
plt.grid()

x0 = np.array(x0)
y0 = np.array(y0)

h0 = (~np.isnan(x0)) & (~np.isnan(y0))
x0 = x0[h0]
y0 = y0[h0]

h1 = (x0 < 1.4) & (x0 > -0.5)
# h1 = x0 < 0.85

x1 = x0[h1]
y1 = y0[h1]
(m, b), cov = np.polyfit(x1, y1, 1, cov=True)
berr = cov[1, 1] ** 0.5

plt.plot(x0, y0, '.', alpha=0.3,color='k')
plt.plot(x1, y1, '.', alpha=0.6,color='teal')
plt.plot(x0, m * x0 + b, '-',color='salmon')
plt.title(r'$\lambda=${}$\AA$, y-intercept$={:.5f}\pm{:.5f}$'.format(filt_lam[f], b, berr))

print(filt_lam[f], '\t', b, berr)
