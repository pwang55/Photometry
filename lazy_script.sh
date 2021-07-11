#!/bin/bash

# -------------------------------------------------------------------------------------------------------
# This script does all the photometry action at once:
# 1. Use SExtractor to create catalogs for each images
# 2. Combine catalogs together and combine with SDSS/PanSTARRS & ALLWISE, then calibrate zero points
# 3. Run EAZY to get photo-z
# 4. Combine EAZY .zout with spec-z and .cat to create .zall
# 5. Copy .zall files to a testing folder for plotting
#
# To disable any step simply comment out the lines. This script has to be run in the data folder.
# Requires:
# - cutoutlist.txt
# - cluster_gal_sdss_radec.csv
# - cluster_star_sdss_radec.csv
# - cluster_panstarrs_radec.csv
# - cluster_panstarrs_extinction.txt
# - merged_eazy_sdss.param
# - merged_eazy_panstarrs.param
# - A folder named eazy_photz
# -------------------------------------------------------------------------------------------------------

if [ ! -d eazy_photz/ ]; then
    mkdir eazy_photz
fi


# Get the absolute directory path of this script so that it can find extra/files
script_dir=$(cd `dirname $0` && pwd)

# Run SExtractor on images
# ~/sources/90Prime/Photometry/make_each_catalog.sh cutoutlist.txt sdss
# ~/sources/90Prime/Photometry/make_each_catalog.sh cutoutlist.txt panstarrs

# Combine and calibrate
python ~/sources/90Prime/Photometry/combine_catalog_calibration.py sdss plot=false
python ~/sources/90Prime/Photometry/combine_catalog_calibration.py panstarrs plot=false

# Run EAZY on combined catalogs
~/eazy-photoz/src/eazy -p merged_eazy_sdss.param
~/eazy-photoz/src/eazy -p merged_eazy_panstarrs.param

# Combine .zout with .cat and z-specs
python ~/sources/90Prime/Photometry/zspecs_combine.py

# TEMP Copy .zall files to plotting folder
cp eazy_photz/*zall ~/90PrimeData/testing/photometry/plots/

