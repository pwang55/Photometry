#!/bin/bash

eazy -p merged_eazy_sdss.param
eazy -p merged_eazy_panstarrs.param

python ~/sources/90Prime/Photometry/zspecs_combine.py

# TEMP
cp eazy_photz/*zall ~/90PrimeData/testing/photometry/plots/
