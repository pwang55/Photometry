from numpy import *
import sys


psf = float(sys.argv[1])

def get_filter(fwhm_star):

	if fwhm_star < 1.75:
		filter_file = 'gauss_1.5_3x3.conv'
	elif fwhm_star < 2.25:
		filter_file = 'gauss_2.0_3x3.conv'
	elif fwhm_star < 2.75:
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

	return filter_file

sys.exit(get_filter(psf))

