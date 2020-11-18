"""

Usage:
    In script folder:
    $ python make_table_for_extinction.py path_to_file/

    In data folder:
    $ python path_to_script/make_table_for_extinction.py

This script generates a table of PanSTARRS and ALLWISE ra/dec to upload to:
    https://irsa.ipac.caltech.edu/applications/DUST/
To get the extinction information.


"""
import numpy as np
from astropy.io import ascii
import sys
from glob import glob
# pathfile = sys.argv[1]

path = ''
if len(sys.argv) == 2:
    path = sys.argv[1]
    if path[-1] != '/':
        path = path + '/'

filenames = glob(path + '*panstarrs_radec.csv')
clustername = filenames[0].split('_')[0]

# filename = pathfile.split('/')[-1]
# filename_nosuffix = filename.split('.')[0]
# path = pathfile[: - len(filename)]

panstarrs_file = path + clustername + '_panstarrs_radec.csv'
allwise_file = path + clustername + '_allwise_radec.txt'

panstarrs = ascii.read(panstarrs_file)
allwise = ascii.read(allwise_file)

# t0 = ascii.read(pathfile)

p1 = panstarrs['raMean', 'decMean']
a1 = allwise['ra', 'dec']

p1.add_column(np.full(len(p1), 2), name='size')
a1.add_column(np.full(len(a1), 2), name='size')
p1.rename_columns(('raMean', 'decMean'), ('ra', 'dec'))
a1.meta = {}

ascii.write(p1, output=path + clustername + '_panstarrs_radec_for_extinction_upload.csv', delimiter=',', format='basic', overwrite=True)
ascii.write(a1, output=path + clustername + '_allwise_radec_for_extinction_upload.csv', delimiter=',', format='basic', overwrite=True)

