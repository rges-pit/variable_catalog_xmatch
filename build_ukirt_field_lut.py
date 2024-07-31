# Program to construct a lookup table of the range of on-sky pointings included
# within each field of the UKIRT microlensing survey.
# The output of this program is a LUT used to accelerate searches for individual
# objects within the UKIRT data

import argparse
import glob
from os import path
import utils
from astropy.table import Table, Column

def build_lut(args):
    """
    Function to construct a lookup table of the range of on-sky pointings included
    within each field of the UKIRT microlensing survey.

    :param args:
    :return: None, outputs LUT in JSON format
    """

    # The UKIRT dataset includes a large list of ASCII-format source tables of
    # all stars within each CCD for each field pointing.
    # There are four CCDs per field pointing.
    # Build a list of the files to iterate over
    source_table_files = glob.glob(path.join(args.ukirt_dir, 'UKIRT_year*md.tbl'))

    # Parse each source table file, and calculate the RA, Dec range of the stars listed.
    lut = Table([
        Column(name='source_table', data=[path.basename(f) for f in source_table_files]),
        Column(name='ra_min', data=[None]*len(source_table_files)),
        Column(name='ra_max', data=[None]*len(source_table_files)),
        Column(name='dec_min', data=[None]*len(source_table_files)),
        Column(name='dec_max', data=[None]*len(source_table_files))
    ])

    for i,source_filename in enumerate(source_table_files):
        print('Parsing source table ' + path.basename(source_filename) + ', '
              + str(i) + ' out of ' + str(len(source_table_files)))

        source_table = utils.load_ukirt_source_table(source_filename)

        lut['ra_min'][i] = source_table['ra'].min()
        lut['ra_max'][i] = source_table['ra'].max()
        lut['dec_min'][i] = source_table['dec'].min()
        lut['dec_max'][i] = source_table['dec'].max()

    # Output the LUT
    utils.output_json_catalog_from_table(lut, args.lut_file)

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('ukirt_dir', help='Path to directory of UKIRT source tables')
    parser.add_argument('lut_file', help='Path to output LUT file')
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    build_lut(args)