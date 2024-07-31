# Variable Catalog Xmatch
Code to cross-match between existing catalogs of variable stars in the Bulge
The code is modular, and designed to run in sequence as follows. 

## Identifying and cross-matching known variables within a region of interest
A number of pre-existing surveys have published catalogs of known variable stars discovered 
within the Galactic Plane and Bulge, including the OGLE and VVV microlensing surveys. 
Since these are very large source catalogs, served from different platforms, it is 
most convenient to extract the subset of these catalogs pertaining to the region of interest 
before cross-matching the entries.  

Two programs in this library query OGLE and VVV survey datasets respectively:
- identify_known_ogle_variables.py
- identify_known_vvv_variables.py

The OGLE and VVV catalogs used in this work are described in the data/README.md file.  
These two survey catalogs are queried first, producing the catalog subsets for the field:
data/rges_ogle4_variables.json 
data/rges_vvv_variables.json

The code combine_variable_catalogs.py is then used to cross-match these two catalogs and to
output a combined catalog data/rges_variable_catalog.json. 

## Identifying UKIRT data on known variables
The UKIRT microlensing survey observed the Galactic Bulge contemporaneously with OGLE and VVV. 
A catalog of microlensing events identified within the data was published, but no independent 
catalog of variables.  However, since the timeseries photometry of all stars in the field of view
was released, we can use the input catalog of known variables produced above to identify 
data on variable objects.  

The first step in this process is to cross-match the rges_variable_catalog.json with the UKIRT 
catalog to identify those stars for which photometry is available.  This is non-trival, due to 
the very large size of the full photometry dataset, so those data are not included in this repository 
but can be accessed from the NASA Exoplanet Database.

This library includes the program identify_variables_in_UKIRT_survey.py which is designed to 
query the UKIRT timeseries dataset held separately on local disk.  
The programs index_ukirt_lightcurves.py and build_ukirt_field_lut.py were written to accelerate
the process of searching the data holdings and retrieving the appropriate lightcurves. 
identify_variables_in_UKIRT_survey.py annotates the data/rges_variable_catalog.json with the 
information necessary to retrieve the appropriate lightcurve data.

## Retrieving optical and NIR timeseries photometry
The program compile_lightcurves.py reads in the rges_variable_catalog.json and compiles the 
available optical and NIR timeseries photometry from OGLE and UKIRT.  
The OGLE photometry is retrieved on a per-object basis from OGLE's online server, while 
the UKIRT photometry is retrieved from a copy of the survey data products on local disk. 

A single, combined file is produced for each star in FITS binary table format, containing 
multiple table extensions, one for each lightcurve in the available passbands.  
These extensions may include one or more of OGLE (I-band), OGLE (V-band), UKIRT (H-band), 
UKIRT (K-band).

## Utilities
The library also contains a number of utility functions in utils.py which are of general use.  
These include functions for reading and writing the various catalog and lightcurve files produced. 