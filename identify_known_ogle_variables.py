# Program to compile a list of variable stars identified from the OGLE survey
# within a specific set of pointings.
# Program is based on original code by Yiannis Tsapras
from astropy import units as u
import rges_survey_definition
import argparse
from os import path


# Configuration
# Location of data files of known OGLE variable catalogs
config = {
    'data_dir': './data/',
    'variable_types': {
        'lpv': 'ogle4_LPV.dat',
        'cephieds': 'ogle4_cephieds.dat',
        'cephieds_type2': 'ogle4_cephieds_type2',
        'delta_scuti': 'ogle4_delta_scuti.dat',
        'eclipsing_binaries': 'ogle4_eclipsing_binaries.dat',
        'hb': 'ogle4_hb.dat',
        'rrlyrae': 'ogle4_rrlyrae.dat'
    }
}

def find_ogle_variables(args, field_centers):
    """
    Function to search the OGLE catalogs of known variables within the list of
    field positions provided.
    Outputs dict ogle_variables catalog of OGLE variables output to file

    :param args: dict of user-provided configuration parameters
    :param field_centers: list of tuples describing field pointings
    :return: None
    """

    # Load and combine the lists of OGLE4 variables of different types

    # Identify those objects within the field of view of all given fields

    # Output catalog of known variables within the field

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
    rges_fields = rges_survey_definition.get_rges_field_centers()
    find_ogle_variables(args, field_centers)