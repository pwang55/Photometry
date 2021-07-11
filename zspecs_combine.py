import warnings
warnings.filterwarnings("ignore", category=UserWarning)
import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
import sys
from glob import glob
from astropy.coordinates import SkyCoord
from astropy.table import Table

# matching distance
dist = 0.8
# std tolerance
std_tolerance = 0.01

path = ''
if len(sys.argv) > 1:
    path = sys.argv[1]
    if path[-1] != '/':
        path = path + '/'

zoutfiles = glob(path + 'eazy_photz/*zout')
clustername = zoutfiles[0].split('/')[-1].split('_')[0]
speczfiles = glob(path + '*_specz_*csv')
# clustername = speczfiles[0].split('_')[0]

s = ascii.read(path + clustername + '_merged_sdss_allwise_eazy.cat')
szout = ascii.read(path + 'eazy_photz/' + clustername + '_merged_sdss_allwise.zout')

p = ascii.read(path + clustername + '_merged_panstarrs_allwise_eazy.cat')
pzout = ascii.read(path + 'eazy_photz/' + clustername + '_merged_panstarrs_allwise.zout')

scoord = SkyCoord(s['ra'], s['dec'], unit='deg')
pcoord = SkyCoord(p['ra'], p['dec'], unit='deg')

# stab = Table()
# ptab = Table()

# stab.add_columns([s['ra'], s['dec'], szout['z_m2'], szout['l68'], szout['u68'], szout['l95'], szout['u95'], s['z_spec']], names=['ra', 'dec', 'z_m2', 'l68','u68','l95','u95', 'zspec_sdss'])
# ptab.add_columns([p['ra'], p['dec'], pzout['z_m2'], pzout['l68'], pzout['u68'], pzout['l95'], pzout['u95'], p['z_spec']], names=['ra', 'dec', 'z_m2', 'l68','u68','l95','u95', 'zspec_sdss'])

stab = szout
ptab = pzout


stab.add_columns([s['ra'], s['dec']], names=['ra', 'dec'], indexes=(0, 0))
ptab.add_columns([p['ra'], p['dec']], names=['ra', 'dec'], indexes=(0, 0))

# Determine how many filters are in each survey
stab.add_columns([[0] * len(stab), [0] * len(stab), [0] * len(stab)], names=['nfilts_narrowband', 'nfilts_calibration', 'nfilts_allwise'])
ptab.add_columns([[0] * len(ptab), [0] * len(ptab), [0] * len(ptab)], names=['nfilts_narrowband', 'nfilts_calibration', 'nfilts_allwise'])

s_narrowband = s.colnames[4::2][:-7]
p_narrowband = p.colnames[4::2][:-7]
for x in range(len(s_narrowband)):
    stab['nfilts_narrowband'] = stab['nfilts_narrowband'] + 1 * (s[s_narrowband[x]] > 0)
    ptab['nfilts_narrowband'] = ptab['nfilts_narrowband'] + 1 * (p[p_narrowband[x]] > 0)

s_calibration = s.colnames[4::2][-7:-2]
p_calibration = p.colnames[4::2][-7:-2]
for x in range(len(s_calibration)):
    stab['nfilts_calibration'] = stab['nfilts_calibration'] + 1 * (s[s_calibration[x]] > 0)
    ptab['nfilts_calibration'] = ptab['nfilts_calibration'] + 1 * (p[p_calibration[x]] > 0)

s_allwise = s.colnames[4::2][-2:]
p_allwise = p.colnames[4::2][-2:]
for x in range(len(s_allwise)):
    stab['nfilts_allwise'] = stab['nfilts_allwise'] + 1 * (s[s_allwise[x]] > 0)
    ptab['nfilts_allwise'] = ptab['nfilts_allwise'] + 1 * (p[p_allwise[x]] > 0)

if len(speczfiles) == 0:
    ascii.write(stab, output=path + 'eazy_photz/' + clustername + '_all_zspecs_sdss.zall', format='ecsv', overwrite=True)
    ascii.write(ptab, output=path + 'eazy_photz/' + clustername + '_all_zspecs_panstarrs.zall', format='ecsv', overwrite=True)

else:


    stab.rename_column('z_spec', 'zspec_sdss')
    ptab.rename_column('z_spec', 'zspec_sdss')


    stab_temp = Table()
    ptab_temp = Table()


    for i in range(len(speczfiles)):
        f = ascii.read(speczfiles[i])
        fra = f['RA']
        fdec = f['DEC']
        fz = f['Redshift']
        fcoord = SkyCoord(fra, fdec, unit='deg')
        idxs, d2ds, d3ds = fcoord.match_to_catalog_sky(scoord)
        idxp, d2dp, d3dp = fcoord.match_to_catalog_sky(pcoord)
        hs = d2ds.arcsec < dist
        hp = d2dp.arcsec < dist
        stab_temp.add_column([-1.0] * len(stab))
        ptab_temp.add_column([-1.0] * len(ptab))
        stab_temp[stab_temp.colnames[-1]][idxs[hs]] = fz[hs]
        ptab_temp[ptab_temp.colnames[-1]][idxp[hp]] = fz[hp]

    stab.add_column([-1.0] * len(stab), name='z_spec')
    ptab.add_column([-1.0] * len(ptab), name='z_spec')

    # zspec_sdss_col_idx = np.where(np.array(stab.colnames) == 'zspec_sdss')[0][0]


    for x in range(len(stab)):
        if stab[x]['zspec_sdss'] > 0:
            stab['z_spec'][x] = stab[x]['zspec_sdss']
        else:
            zspecs = []
            for i in range(len(stab_temp.colnames)):
                zspecs.append(stab_temp[x][stab_temp.colnames[i]])
            zspecs = np.array(zspecs)
            zspecs = zspecs[zspecs > 0]
            if len(zspecs) > 0 and np.std(zspecs) < std_tolerance:
                stab['z_spec'][x] = np.median(zspecs)
            else:
                stab['z_spec'][x] = -1.0

    for x in range(len(ptab)):
        if ptab[x]['zspec_sdss'] > 0:
            ptab['z_spec'][x] = ptab[x]['zspec_sdss']
        else:
            zspecs = []
            for i in range(len(ptab_temp.colnames)):
                zspecs.append(ptab_temp[x][ptab_temp.colnames[i]])
            zspecs = np.array(zspecs)
            zspecs = zspecs[zspecs > 0]
            if len(zspecs) > 0 and np.std(zspecs) < std_tolerance:
                ptab['z_spec'][x] = np.median(zspecs)
            else:
                ptab['z_spec'][x] = -1.0


    ascii.write(stab, output=path + 'eazy_photz/' + clustername + '_all_zspecs_sdss.zall', format='ecsv', overwrite=True)
    ascii.write(ptab, output=path + 'eazy_photz/' + clustername + '_all_zspecs_panstarrs.zall', format='ecsv', overwrite=True)




