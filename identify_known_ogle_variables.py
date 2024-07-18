# Program to compile a list of variable stars identified from the OGLE survey
# within a specific set of pointings.
# Program is a reimplementation of original code by Yiannis Tsapras
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, Column
import rges_survey_definition
import argparse
from os import path
import utils


# Configuration
# Location of data files of known OGLE variable catalogs
config = {
    'data_dir': './data/',
    'variable_catalogs': {
        'lpv': {'file': 'ogle4_LPV.dat', 'columns': {'name': 0, 'ra': 2, 'dec': 3}},
        'cephied': {'file': 'ogle4_cephieds.dat', 'columns': {'name': 0, 'ra': 2, 'dec': 3}},
        'cephied_type2': {'file': 'ogle4_cephieds_type2.dat', 'columns': {'name': 0, 'ra': 2, 'dec': 3}},
        'delta_scuti': {'file': 'ogle4_delta_scuti.dat', 'columns': {'name': 0, 'ra': 2, 'dec': 3}},
        'eclipsing_binary': {'file': 'ogle4_eclipsing_binaries.dat', 'columns': {'name': 0, 'ra': 2, 'dec': 3}},
        'hb': {'file': 'ogle4_hb.dat', 'columns': {
            'name': 0, 'ra': [2, 3, 4], 'dec': [5, 6, 7]
        }},
        'rrlyrae': {'file': 'ogle4_rrlyrae.dat', 'columns': {'name': 0, 'ra': 2, 'dec': 3}},
    }
}

def find_ogle_variables(args, config):
    """
    Function to search the OGLE catalogs of known variables within the list of
    field positions provided.
    Outputs dict ogle_variables catalog of OGLE variables output to file

    :param args: dict of user-provided configuration parameters
    :param config dict of standard configurable key-value pairs
    :return: None
    """

    # Load and combine the lists of OGLE4 variables of different types
    ogle_catalog = load_ogle_variable_catalogs(config)

    # Identify those objects within the field of view of all given fields
    survey_catalog = utils.find_variables_in_fov(ogle_catalog, coord_type='sexigesimal')

    # Output catalog of known variables within the field
    utils.output_json_catalog(survey_catalog, args.output_file)

def load_ogle_variable_catalogs(config):
    """
    Function to load the catalogs of known OGLE variables from the ASCII .dat
    format files downloaded from the OGLE 4 catalog website
    https://www.astrouw.edu.pl/ogle/ogle4/OCVS/blg/

    :param config:
    :return: ogle_catalog dict Names, types and sky positions of variables
    """

    # OGLE4 variable catalogs consist of a set of ASCII files of differing formats,
    # so here we extract the relevant columns in each case and compile them into
    # a single dictionary
    catalog = {
        'name': [],
        'type': [],
        'ra': [],
        'dec': []
    }
    for cat_type, cat_config in config['variable_catalogs'].items():
        with open(path.join(config['data_dir'],cat_config['file']), 'r') as f:
            lines = f.readlines()
            for entry in lines:
                if len(entry.replace('\n','')) > 0:
                    data = entry.replace('\n', '').split()
                    name = data[cat_config['columns']['name']]
                    if isinstance(cat_config['columns']['ra'], int):
                        ra = data[cat_config['columns']['ra']]
                        dec = data[cat_config['columns']['dec']]
                    else:
                        ra = data[cat_config['columns']['ra'][0]] \
                                   + ':' + data[cat_config['columns']['ra'][1]] \
                                   + ':' + data[cat_config['columns']['ra'][2]]
                        dec = data[cat_config['columns']['dec'][0]] \
                             + ':' + data[cat_config['columns']['dec'][1]] \
                             + ':' + data[cat_config['columns']['dec'][2]]

                    catalog['name'].append(name)
                    catalog['type'].append(cat_type)
                    catalog['ra'].append(ra)
                    catalog['dec'].append(dec)

    ogle_catalog = Table([
        Column(name='Name', data=catalog['name']),
        Column(name='Type', data=catalog['type']),
        Column(name='RA', data=catalog['ra']),
        Column(name='Dec', data=catalog['dec'])
    ])

    print('Loaded ' + str(len(ogle_catalog)) + ' known variables from OGLE 4')

    return ogle_catalog

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
    find_ogle_variables(args, config)