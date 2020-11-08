from numpy import *
import sys
from glob import glob


# Finding the correct column for FWHM and MAG_AUTO
find_column_list=[]
file0=sys.argv[1]


with open(file0,'r') as find_column:
	find_column_list = [l.split() for l in find_column]
	
mag_c = [ci for ci in range(20) if find_column_list[ci][2]=='MAG_AUTO'][0]
fwhm_c = [ci for ci in range(20) if find_column_list[ci][2]=='FWHM_IMAGE'][0] 
peak_c = [ci for ci in range(20) if find_column_list[ci][2]=='FLUX_MAX'][0] 


f = loadtxt(file0,dtype=None)
h = (f[:,mag_c]>-20) & (f[:,mag_c]<-4) & (f[:,fwhm_c]<12) & (f[:,peak_c] < 65000.0) & (f[:,peak_c] > 4000.0)
dat = f[h]

fwhm = dat[:,fwhm_c]
hist1 = histogram(fwhm,bins=200)

m = argmax(hist1[0])
fwhm_star = hist1[1][m]

#print hist1[0][m-1],fwhm_star,hist1[0][m],hist1[0][m+1]

#print 'fwhm = ', fwhm_star

#int_fwhm = int(round(fwhm_star))
#print int_fwhm
filter_file=''

#gauss = glob('gauss_*conv')

#print gauss
#ng = len(gauss)

if fwhm_star < 1.75:
	filter_file = 'gauss_1.5_3x3.conv'
elif fwhm_star < 2.2:
	filter_file = 'gauss_2.0_3x3.conv'
elif fwhm_star < 2.7:
	filter_file = 'gauss_2.5_5x5.conv'
elif fwhm_star < 3.5:
	filter_file = 'gauss_3.0_7x7.conv'
elif fwhm_star < 4.5:
	filter_file = 'gauss_4.0_7x7.conv'
elif fwhm_star < 5.5:
	filter_file = 'gauss_5.0_9x9.conv'
elif fwhm_star < 7.5:
	filter_file = 'gauss_6.0_13x13.conv'
elif fwhm_star < 11.5:
	filter_file = 'gauss_9.0_19x19.conv'
elif fwhm_star < 15.5:
	filter_file = 'gauss_15.0_31x31.conv'
else:
	filter_file = 'gauss_16.0_33x33.conv'

#print filter_file

sys.exit(filter_file)

