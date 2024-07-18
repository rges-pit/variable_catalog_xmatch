# Program to compile a list of variable stars identified from the VVV survey
# within a specific set of pointings.

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, Column
from astropy.io import fits
import rges_survey_definition
import argparse
from os import path
import utils

# Configuration
# Location of data files
config = {
    'data_dir': './data/',
    'variable_catalogs': {
        'lpv': {'file': 'VVV_LPV_J_A+A_660_A35_tableb2.fits', 'columns': [0,2,3]},
        'miras': {'file': 'VVV_Miras_J_A+A_660_A35_tableb1.fits', 'columns': [0,1,2]},
    }
}

def find_vvv_variables(args, config):
    """
    Function to identify known variables from the VVV catalogs within the RGES survey
    region

    :param args: dict of user-provided configuration parameters
    :param config: dict of standard configurable key-value pairs
    :return: None, outputs JSON-format catalog
    """

    # Load and combine the lists of OGLE4 variables of different types
    vvv_catalog = load_vvv_variable_catalogs(config)

    # Identify those objects within the field of view of all given fields
    survey_catalog = utils.find_variables_in_fov(vvv_catalog, coord_type='degrees')

    # Output catalog of known variables within the field
    utils.output_json_catalog(survey_catalog, args.output_file)


def load_vvv_variable_catalogs(config):
    """
    Function to load the variable star catalogs from VVV.
    These are accessed via the FITS binary tables served from the CDS archive.
    See data/README.md for details.

    :param config:
    :return: dict combine catalog of variable stars
    """
    catalog = {
        'name': [],
        'type': [],
        'ra': [],
        'dec': []
    }
    for cat_type, cat_config in config['variable_catalogs'].items():
        with fits.open(path.join(config['data_dir'],cat_config['file'])) as hdul:
            data = hdul[1].data
            for entry in data:
                catalog['name'].append(str(entry[cat_config['columns'][0]]))
                catalog['type'].append(cat_type)
                catalog['ra'].append(float(entry[cat_config['columns'][1]]))
                catalog['dec'].append(float(entry[cat_config['columns'][2]]))

    vvv_catalog = Table([
        Column(name='Name', data=catalog['name']),
        Column(name='Type', data=catalog['type']),
        Column(name='RA', data=catalog['ra']),
        Column(name='Dec', data=catalog['dec'])
    ])

    print('Loaded ' + str(len(vvv_catalog)) + ' known variables from VVV')

    return vvv_catalog

def get_args():
    """
    Gather user-provided input

    :return: args dict Set of argument key-value pairs
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('output_file', help='Path to output file')
    args = parser.parse_args()

    return args

if __name__ == '__main__':
    args = get_args()
    find_vvv_variables(args, config)