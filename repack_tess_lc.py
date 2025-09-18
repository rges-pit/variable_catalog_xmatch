# Utility to repackage the TESS LC into the standard format
# used throughout the rest of this package, for easier compatibility with
# software that uses its output.

import argparse
import utils
from astropy.io import fits

def run(args):
    """
    Function to run the repackaging of TESS lightcurves

    :param args: Arguments object
    :return: None
    """

    # Load the source catalog
    tess_catalog = utils.load_json_catalog(args.tess_cat_file, decimal_degrees=True)

    # Parse the TESS lightcurves into the standard format
    for tic_id, params in tess_catalog.items():

        # Each object may have multiple TESS lightcurves from different survey passes
        lc_data, hdr_info = utils.load_tess_lcs_for_star(args, tic_id, params)

        # Output to standard format lightcurve
        if lc_data:
            hdr = utils.make_lc_header(tic_id, params)

            photometry = {'TESS': lc_data}
            utils.output_multiband_lc(args, tic_id, params, hdr, photometry)

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('tess_cat_file', help='Path to source catalog file')
    parser.add_argument('input_dir', help='Path to input data directory')
    parser.add_argument('output_dir', help='Path to output data directory')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = get_args()
    run(args)