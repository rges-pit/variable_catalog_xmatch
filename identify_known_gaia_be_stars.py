# Script crossmatches Kruzsynska's catalog of Be stars
# against the RGES field of view

# This code uses an input catalog of Be Stars produced by
# Katarzyna Kruszynska from the Gaia Alerts system

from os import path
import utils

from astropy.table import Table, Column
import argparse
from os import path
import numpy as np

# Configuration
# Location of data files
config = {
    'data_dir': './data/',
    'variable_catalogs': {
        'bestar': {'file': 'kruszynska_gaia_be_candidates.csv', 'columns': [0,2,3]},
    }
}

def find_gaia_variables(args, config):
    """
    Function to identify those stars from the Gaia catalog of Be star candidates
    that lie within the RGES field of view.

    :param args: dict of user-provided configuration parameters
    :param config: dict of standard configurable key-value pairs
    :return: None, outputs JSON-format catalog
    """

    # Load and combine the lists of Gaia candidate Be Stars
    gaia_catalog = load_gaia_variable_catalog(config)

    # Identify those objects within the field of view of all given fields
    survey_catalog = utils.find_variables_in_fov(gaia_catalog, coord_type='degrees')

    # Output catalog of known variables within the field
    utils.output_json_catalog_from_table(survey_catalog, args.output_file)

def load_gaia_variable_catalog(config):
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

    cat_type = 'bestar'
    cat_config = config['variable_catalogs'][cat_type]

    data = np.loadtxt(
        path.join(config['data_dir'], cat_config['file']),
        delimiter=',',
        dtype='U9,float,float',
        usecols=(cat_config['columns']),
        unpack=True
    )

    gaia_catalog = Table(
        [
            Column(name='Name', data=data[0]),
            Column(name='Type', data=[cat_type]*len(data[0])),
            Column(name='RA', data=data[1]),
            Column(name='Dec', data=data[2])
        ]
    )

    print('Loaded ' + str(len(gaia_catalog)) + ' known variables from VVV')

    return gaia_catalog


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
    find_gaia_variables(args, config)