# Script to build a complete JSON-format catalogs of flare stars
# identified from TESS data
#
# Based on MAST documentation:
# https://github.com/spacetelescope/mast_notebooks/blob/f552b976996c5b022d2676f886400a4b49485c9c/notebooks/TESS/interm_tasoc_lc/interm_tasoc_lc.ipynb

from os import path
import argparse
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table, Column
from astroquery.mast import Catalogs
import utils

# Published catalogs of known variables based on TESS data
# Gunther et al. (2020), AJ, 159, 60
# https://ui.adsabs.harvard.edu/abs/2020AJ....159...60G/abstract
# Medina et al. (2020), ApJ, 905, 107
# https://ui.adsabs.harvard.edu/abs/2020ApJ...905..107M/abstract

def build_catalog(args):

    # Load catalogs of published variable stars.  Naturally these are in different
    # formats, hence the separate load functions
    input_catalog = {}
    input_catalog = load_catalog(args, input_catalog, 'Gunther')
    input_catalog = load_catalog(args, input_catalog, 'Medina')
    print('Loaded a total of ' + str(len(input_catalog))
          + ' known flare stars from published catalogs')

    # Not all TESS stars from the variable catalogs have RA, Decs,
    # so we query for those if not available.  It is possible at this stage
    # to get no hit in the catalog.
    catalog = {
        'name': [],
        'type': [],
        'ra': [],
        'dec': []
    }
    for tic_id, star_data in input_catalog.items():

        if not star_data['RA'] or not star_data['Dec']:
            TICentry = get_TESS_star_from_MAST(tic_id)
            if TICentry:
                star_data['RA'] = TICentry['ra'][0]
                star_data['Dec'] = TICentry['dec'][0]
                catalog['name'].append(tic_id)
                catalog['type'].append('flare')
                catalog['ra'].append(star_data['RA'])
                catalog['dec'].append(star_data['Dec'])
        else:
            catalog['name'].append(tic_id)
            catalog['type'].append('flare')
            catalog['ra'].append(star_data['RA'])
            catalog['dec'].append(star_data['Dec'])

    tess_catalog = Table([
        Column(name='Name', data=catalog['name']),
        Column(name='Type', data=catalog['type']),
        Column(name='RA', data=catalog['ra']),
        Column(name='Dec', data=catalog['dec'])
    ])

    print('Found a total of ' + str(len(tess_catalog))
          + ' known flare stars with TIC entries')

    # Output catalog of known variables within the field
    utils.output_json_catalog_from_table(tess_catalog, args.output_file)

def get_TESS_star_from_MAST(tic_id):
    """
    Function to retrieve a TESS stars source catalog entry from MAST

    :param tic_id: integer TESS Input Catalog identitier
    :return: Table with TIC entry for exactly 1 star or None
    """

    catalog_entry = Catalogs.query_object(
        'TIC ' + str(tic_id),
        radius=0.0006,
        catalog="TIC"
    )

    if len(catalog_entry) == 1:
        return catalog_entry
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
    parser.add_argument('data_dir', help='Path to input data directory')
    parser.add_argument('output_file', help='Path to output file')

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    build_catalog(args)