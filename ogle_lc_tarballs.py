from os import path
import glob
import argparse
import utils

def run(args):

    # List all existing multiband lightcurves
    lc_list = glob.glob(path.join(args.lc_dir, '*_multiband_lc.fits'))

    # Loop over all existing lightcurves.  If a corresponding OGLE lightcurve
    # exists, parse it and add the data to the multiband FITS table
    for lc_file in lc_list:
        star_id = path.basename(lc_file).split('_')[0]
        ogle_lc_file = path.join(args.ogle_lc_dir, star_id + '.dat')
        if path.isfile(ogle_lc_file):

            # Based on the paper "One Thousand New Dwarf Novae from the OGLE Survey"
            # by P. Mroz et al., all lightcurves in this dataset are I-band.
            ogleI = utils.parse_ogle_lc_file(ogle_lc_file)

            # Load the existing multiband lightcurve
            hdr, photometry = utils.load_multiband_lc(lc_file)

            # Add OGLE photometry
            photometry['LC_I'] = ogleI

            # Output updated multiband lightcurve
            utils.output_multiband_lc(args, star_id, hdr, photometry)

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('lc_dir', help='Path to directory of mulitband lightcurves')
    parser.add_argument('ogle_lc_dir', help='Path to directory of OGLE lightcurves from tarball')
    parser.add_argument('output_dir', help='Path to output directory for revised multiband lightcurves')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = get_args()
    run(args)