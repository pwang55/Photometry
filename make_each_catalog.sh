#!/bin/bash


usage='

    In script folder:
    $ ./make_each_catalog.sh path_to_files/cutoutlist.txt

    In data folder:
    $ path_to_script/make_each_catalog.sh cutoutlist.txt


This script creates final photometry catalog for each filter. Also check if SDSS, PanSTARRS catalog and PanSTARRS extinction file exist.

'

# If no argument is given, print doc and exit
if (( "$#" == 0 )); then
	echo "$usage"
	exit 1
fi


# Get the absolute directory path of this script so that it can find extra/files
script_dir=$(cd `dirname $0` && pwd)

list_in=$1
path=""

# If argument 1 is a full path, path variable will be set to the path excluding the list name
len_list=`echo $list_in | awk '{n=split($1,a,"/"); print n}'`
if (( $len_list > 1 )); then
        path=`dirname $list_in`
        path=$path/
fi



for file in `cat $list_in`; do
    filename=`echo $file | sed -e 's/\.fits//g'`
    base2=`echo ${filename:(-3)}`
    # Get the initial catalog for fwhm estimation
    sex -c ${script_dir}/extra/find_fwhm.sex ${path}$file -CATALOG_NAME ${path}${filename}_fwhm.cat -PARAMETERS_NAME ${script_dir}/extra/fwhm.param -FILTER_NAME ${script_dir}/extra/gauss_4.0_7x7.conv -WEIGHT_IMAGE ${path}${filename}.wt.fits
    conv_filter=$(python ${script_dir}/extra/fwhm.py ${path}${filename}_fwhm.cat 2>&1)
    rm ${path}${filename}_fwhm.cat
    # Get the catalog for psf
    sex -c ${script_dir}/extra/psf.sex ${path}$file -CATALOG_NAME ${path}${filename}_psf.cat -CATALOG_TYPE FITS_LDAC -PARAMETERS_NAME ${script_dir}/extra/psf.param -FILTER_NAME ${script_dir}/extra/${conv_filter} -WEIGHT_IMAGE ${path}${filename}.wt.fits
    psfex -c ${script_dir}/extra/psf.psfex ${path}${filename}_psf.cat
    rm ${path}${filename}_psf.cat
    # Find the final photometry
    sex -c ${script_dir}/extra/photometry.sex ${path}$file -CATALOG_NAME ${path}${base2}.cat -FILTER_NAME ${script_dir}/extra/${conv_filter} -WEIGHT_IMAGE ${path}${filename}.wt.fits -PSF_NAME ${path}${filename}_psf.psf -PARAMETERS_NAME ${script_dir}/extra/photometry.param
done

clustername=`echo $filename | awk '{split($1,a,"_"); print a[3]}'`



if [ ! -e ${path}${clustername}_star_sdss_radec.csv ] && [ ! -e ${path}${clustername}_gal_sdss_radec.csv ]; then
    echo "SDSS catalogs doesn't exst!"
fi

if [ ! -e ${path}${clustername}_panstarrs_radec.csv ]; then
    echo "PanSTARRS catalogs doesn't exst!"
fi

if [ ! -e ${path}${clustername}_panstarrs_extinction.txt ]; then
    echo "PanSTARRS extinction doesn't exst!"
fi
