
**Create Catalogs for Each Cutout Stack**

In script directory:

`./make_each_catalog.sh path_to_files/cutoutlist.txt`

In data directory:

`~/path_to_script/make_each_catalog.sh cutoutlist.txt`

`cutoutlist.txt` should contain all the cutout stack images.
This script creates catalogs such as 1ne.cat.

------------------------------------

**Create Table for PanSTARRS Extinction Upload**

In script directory:

`$ python make_table_for_ps1_extinction.py path_to_file/cluster_panstarrs_radec.csv`

In data directory:

`$ path_to_script/python make_table_for_ps1_extinction.py cluster_panstarrs_radec.csv`

This script generates a table of PanSTARRs ra/dec to upload to:
https://irsa.ipac.caltech.edu/applications/DUST/ to get the extinction information.


------------------------------------

**Combine all Generated Catalogs/Zero points Calibration/Match with SDSS/PanSTARRS and ALLWISE**

In script directory:

`$ python combine_catalog_calibration.py [path to catalogs directory] [sdss/panstarrs] (optional arguments)`

In catalog directory:

`$ python path_to_script/combine_catalog_calibration.py [sdss/panstarrs] (optional arguments)`

Required files:

- cluster_gal_sdss_radec.csv
- cluster_star_sdss_radec.csv
- cluster_panstarrs_radec.csv
- cluster_panstarrs_extinction.txt
- cluster_allwise_radec.txt (https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-scan?mission=irsa&submit=Select&projshort=WISE)


This script search for 1ne.cat~4sw.cat and match them with SDSS/PanSTARRS catalog, create mycatalog then merged with ALLWISE catalog.

This script creates:

- clustername_mycatalog_sdss(panstarrs)_gal.csv                   calibrated and selected magnitudes, will change even if match=False
- clustername_mycatalog_sdss(panstarrs)_star_psf(auto).csv        calibrated psf and auto magnitudes for stars catalogs
- clustername_mycatalog_sdss(panstarrs)_gal(star)_eazy.cat        calibrated mycatalog as EAZY input format
- clustername_mycatalog_merged_sdss(panstarrs)_allwise_eazy.cat   calibrated mycatalog matched/merged with ALLWISE as EAZY input format
- clustername_mycatalog_sdss(panstarrs)_gal(star)_original.ecsv   uncalibrated magnitudes, just matched with SDSS/PanSTARRS, will be read by this script if match=False
- clustername_mycatalog_sdss(panstarrs)_gal(star)_ds9.csv         csv file with just ra/dec for DS9 catalog tool
- clustername_mycatalog_gal(star)_w_sdss_spec_radec.txt           txt file with ra/dec of objects that have SDSS spectra, can be upload to:
                                                                 https://dr16.sdss.org/optical/spectrum/search to get the spectra fits files

------------------------------------

**Plot SDSS Spectra**

Galaxies:

`$ python path_to_script/plot_gal_spec_pdf.py clustername_merged_sdss(panstarrs)_allwise_eazy.cat`

`$ python plot_gal_spec_pdf.py path_to_files/clustername_merged_sdss(panstarrs)_allwise_eazy.cat`

Requires:
- clustername_merged_sdss(panstarrs)_allwise_eazy.cat
- clustername_merged_sdss(panstarrs)_allwise.zout
- All sdss spectra in the same folder

Stars:

`$ python path_to_script/plot_star_spec_pdf.py clustername_mycatalog_sdss(panstarrs)_star_auto(psf).csv`

`$ python plot_star_spec_pdf.py path_to_files/clustername_mycatalog_sdss(panstarrs)_star_auto(psf).csv`

Requires:
- clustername_mycatalog_sdss(panstarrs)_star_auto(psf).csv
- All sdss spectra in the same folder

------------------------------------

From start to finish, in data folder:

`$ ~/path_to_script/lazy_script.sh`

Requires:
- cutoutlist.txt
- cluster_gal_sdss_radec.csv
- cluster_star_sdss_radec.csv
- cluster_panstarrs_radec.csv
- cluster_panstarrs_extinction.txt
- merged_eazy_sdss.param
- merged_eazy_panstarrs.param
- A folder named eazy_photz
