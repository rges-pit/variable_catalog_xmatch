# Script to download TESS mission lightcurves for catalogs of known
# variable stars
#
# Based on MAST documentation:
# https://github.com/spacetelescope/mast_notebooks/blob/f552b976996c5b022d2676f886400a4b49485c9c/notebooks/TESS/interm_tasoc_lc/interm_tasoc_lc.ipynb

from os import path
import argparse
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.mast import Observations
import utils

def download_lc(args):

    # Load catalog of TESS published variable stars.
    tess_catalog = utils.load_json_catalog(args.tess_cat_file, decimal_degrees=True)
    print('Loaded a total of ' + str(len(tess_catalog))
          + ' known flare stars from published catalogs')

    # Query the MAST archive for TESS lightcurve data products for the stars
    # in this catalog
    for tic_id, params in tess_catalog.items():
        print('Downloading lightcurves for star ', tic_id)
        manifest = get_TESS_lc_from_MAST(args, tic_id)
        if manifest:
            params['lc_files'] = [path.basename(f) for f in manifest['Local Path']]
        else:
            params['lc_files'] = []
        tess_catalog[tic_id] = params

    # Output updated catalog
    utils.output_json_catalog_from_dict(tess_catalog, args.tess_cat_file)


def get_TESS_lc_from_MAST(args, tic_id):
    """
    Function to query the MAST catalog for FAST cadence lightcurves
    for a given TIC identifier.

    :param args: Argument parser object
    :param tic_id: integer TESS Input Catalog identitier
    :return: mission_manifest, Table, list of downloaded dataproducts for one star
    """

    # Fetch lightcurves for the selected stars from MAST
    # MAST doesn't seem to support queries to identify specific data products,
    # so this query serves a list of all related lightcurve and target pixel
    # data products of all cadences.
    mission_results = Observations.query_criteria(obs_collection="TESS",
                                                    target_name=tic_id)

    # Get list of corresponding data products
    print(mission_results)
    if mission_results:
        mission_products = Observations.get_product_list(mission_results)

        # Downselect from the table to identify the fast cadence lightcurves only
        # _lc.fits, _fast-lc.fits indicate 2min and 20s cadence lightcurves
        # FAST-LC selects the higher-cadence lightcurves as opposed to the target
        # pixel products
        # This often produces multiple individual lightcurves corresponding to
        # different observation sectors
        idx = mission_products['productSubGroupDescription'] == 'FAST-LC'

        # Download the lightcurve data products
        # PDCSAP_FLUX = flux series that has the common instrumental systematics removed
        mission_manifest = Observations.download_products(
            mission_products[idx],
            download_dir=args.output_dir
        )

        return mission_manifest

    else:
        None

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
    parser.add_argument('tess_cat_file', help='Path to TESS catalog of stars')
    parser.add_argument('output_dir', help='Path to output data directory')

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    download_lc(args)