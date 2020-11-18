#!/bin/bash


if [ ! -d eazy_photz/ ]; then
    mkdir eazy_photz
fi


# Get the absolute directory path of this script so that it can find extra/files
script_dir=$(cd `dirname $0` && pwd)

~/sources/90Prime/Photometry/make_each_catalog.sh cutoutlist.txt

python ~/sources/90Prime/Photometry/combine_catalog_calibration.py sdss plot=false
python ~/sources/90Prime/Photometry/combine_catalog_calibration.py panstarrs plot=false

~/eazy-photoz/src/eazy -p merged_eazy_sdss.param
~/eazy-photoz/src/eazy -p merged_eazy_panstarrs.param


