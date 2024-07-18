# Program to identify UKIRT lightcurves, if available, for variables identified by
# other Bulge surveys including OGLE, VVV
# While it is possible to search for individual stars in the UKIRT data through the
# NExSci archive website, bulk queries for large numbers of stars require this
# automated approach

import argparse
import utils

def find_ukirt_data_for_variables(args):
    """
    Function to search the source tables of the UKIRT microlensing survey for entries
    corresponding to stars from an input catalog.

    :param args: ArgumentParser object
    :return: None, outputs updated JSON catalog with lightcurve files where available
    """

    # Load the input catalog of stars to search for.  Add column(s) to store
    # UKIRT cross-matching information
    var_catalog = utils.load_json_catalog(args.catalog)
    var_catalog.add_column([None] * len(var_catalog), name='UKIRT_source_table')
    var_catalog.add_column([None] * len(var_catalog), name='UKIRT_lc_files')

    # Load UKIRT look-up table of field pointings
    lut = utils.load_json_catalog(args.lut)

    # For each catalog star, check whether it falls within a UKIRT field.
    # If it does, load the source table for that field and identify the corresponding
    # lightcurve file, if any.
    for star_id, params in var_catalog.items():
        source_tables = find_star_in_LUT(lut, params['RA'], params['Dec'])

        if len(source_tables) > 0:
            lc_files = []

    # Output updated catalog

def find_star_in_source_table(args, source_table_file, ra, dec):
    """
    Function to search the source table for a single field and CCD to determine
    whether at star of the given RA and Dec was detected within the field.
    If so, the function returns the year, field, CCD and star ID necessary
    to locate the corresponding lightcurve file.

    :param args: parsed arguments for this software
    :param source_table_file: Path to a single source table file
    :param ra: RA of star in decimal degrees
    :param dec: Dec of star in decimal degrees
    :return: dict of the parameters needed to identify the specific lightcurve file
    """


def find_star_in_LUT(lut, ra, dec):
    """
    Function to check whether a star is within the RA, Dec ranges of any field and CCD
    within the UKIRT survey.  If found, a list of filename(s) of the source table(s) of the
    corresponding field(s) and CCD(s) is returned.  Otherwise, an empty list is returned.

    :param ra: RA of star to search for
    :param dec: Dec of star to search for
    :return: source_tables list Filenames of source tables which may contain the star
    """

    raidx = np.where(ra >= lut['ra_min'] & ra <= lut['ra_max'])
    decidx = np.where(dec >= lut['dec_min'] & dec <= lut['dec_max'])
    idx = list(set(raidx).intersection(decidx))

    return lut['source_table'][idx]

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('catalog', help='Path to input star catalog')
    parser.add_argument('lut', help='Path to UKIRT LUT')
    parser.add_argument('ukirt_dir', help='Path to directory of UKIRT source tables')
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    find_ukirt_data_for_variables(args)