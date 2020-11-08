"""

Usage:
    In script folder:
    $ python make_table_for_ps1_extinction.py path_to_file/cluster_panstarrs_radec.csv

    In data folder:
    $ path_to_script/python make_table_for_ps1_extinction.py cluster_panstarrs_radec.csv

This script generates a table of PanSTARRs ra/dec to upload to:
    https://irsa.ipac.caltech.edu/applications/DUST/
To get the extinction information.


"""
import numpy as np
from astropy.io import ascii
import sys

pathfile = sys.argv[1]
filename = pathfile.split('/')[-1]
filename_nosuffix = filename.split('.')[0]
path = pathfile[: - len(filename)]

t0 = ascii.read(pathfile)

t1 = t0['raMean', 'decMean']

t1.add_column(np.full(len(t1), 2), name='size')
t1.rename_columns(('raMean', 'decMean'), ('ra', 'dec'))


ascii.write(t1, output=path + filename_nosuffix + '_for_extinction_upload.csv', delimiter=',', format='basic', overwrite=True)

