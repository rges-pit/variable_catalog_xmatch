# Script to download TESS mission lightcurves for catalogs of known
# variable stars
#
# Based on MAST documentation:
# https://github.com/spacetelescope/mast_notebooks/blob/f552b976996c5b022d2676f886400a4b49485c9c/notebooks/TESS/interm_tasoc_lc/interm_tasoc_lc.ipynb

from os import path
import argparse
import json
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.mast import Tesscut, Observations, Catalogs

# Published catalogs of known variables based on TESS data
# Gunther et al. (2020), AJ, 159, 60
# https://ui.adsabs.harvard.edu/abs/2020AJ....159...60G/abstract
# Medina et al. (2020), ApJ, 905, 107
# https://ui.adsabs.harvard.edu/abs/2020ApJ...905..107M/abstract

def download_lc(args):

    # Load catalogs of published variable stars.  Naturally these are in different
    # formats, hence the separate load functions
    var_stars = {}
    var_stars = load_catalog(args, var_stars, 'Gunther')
    var_stars = load_catalog(args, var_stars, 'Medina')
    print('Loaded a total of ' + str(len(var_stars))
          + ' known flare stars from published catalogs')

    # Fetch lightcurves for the selected stars from MAST
    tic_id = '254231212'
    mission_res = Observations.query_criteria(obs_collection="TESS",
                                              target_name=tic_id)

    # Get list of corresponding data products
    mission_products = Observations.get_product_list(mission_res)
    #mission_products.pprint_all()

    # Download the lightcurve data products
    # Specify lightcurve only not target pixel products
    # _lc.fits, _fast-lc.fits indicate 2min and 20s cadence lightcurves
    # PDCSAP_FLUX = flux series that has the common instrumental systematics removed
    mission_manifest = Observations.download_products(mission_products, download_dir=args.output_dir)
    mission_manifest.pprint_all()

#def fetch_mast_lc(tic, ):

def load_catalog(args, var_stars, cat_name):

    if cat_name == 'Gunther':
        file_name = 'Gunther_TESS_flare_catalog_ajab5d3at1_mrt.txt'
    elif cat_name == 'Medina':
        file_name = 'Medina_TESS_var_catalog_apjabc686t1_mrt.txt'
    else:
        raise IOError('Unrecognized catalog name: ' + cat_name)

    file_path = path.join(args.data_dir, file_name)
    if not path.isfile(file_path):
        raise IOError('Cannot find input catalog file at ' + file_path)

    with open(file_path, 'r') as f:
        file_lines = f.readlines()

        for line in file_lines[47:]:
            entries = line.replace('\n','').split()

            if cat_name == 'Gunther':
                tic = entries[0]
                if tic not in var_stars.keys():
                    var_stars[tic] = {
                        "Type": "flare",
                        "RA": None,
                        "Dec": None
                    }
            else:
                tic = line[24:33]
                ra = line[34:45].replace(' ',':')
                dec = line[45:57].lstrip().replace(' ',':')
                s = SkyCoord(ra, dec, frame='icrs', unit=(u.hourangle, u.deg))
                if tic not in var_stars.keys():
                    var_stars[tic] = {
                        "Type": "flare",
                        "RA": s.ra.deg,
                        "Dec": s.dec.deg
                    }

    return var_stars

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('data_dir', help='Path to input data directory')
    parser.add_argument('output_dir', help='Path to output data directory')

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    download_lc(args)