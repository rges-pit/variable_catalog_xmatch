# Program to gather and compile optical and NIR timeseries photometry for
# a catalog of known variable stars, using the OGLE and UKIRT photometry catalogs

import argparse
import utils

def gather_data(args):
    """
    Function to gather and compile optical and NIR timeseries photometry for
    a catalog of known variable stars, using the OGLE and UKIRT photometry catalogs

    :param args: Program arguments
    :return: None, outputs a set of lightcurve files
    """

    # Load catalog of variable stars
    var_catalog = utils.load_json_catalog(args.catalog)
    #print(var_catalog)

    # Load UKIRT lightcurve index
    #ukirt_index = utils.load_ukirt_index(args.ukirt_index)

    # For each star in the catalog
    for star_id, data in var_catalog.items():

        if data['UKIRT_lc_files']:
            print(star_id, data)

        # Harvest OGLE timeseries photometry from online archive
        #ogleI, ogleV = utils.fetch_ogle_photometry(star)

        # Fetch the UKIRT lightcurve for the star
        #ukirtH, ukirtK = utils.fetch_ukirt_photometry(star, ukirt_index)

        # Output combined lightcurve
        # photometry = {
        #     'I': ogleI,
        #     'V': ogleV,
        #     'H': ukirtH,
        #     'K': ukirtK
        # }
        # output_multiband_lc(args, star, photometry)



def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('catalog', help='Path to variable input catalog')
    parser.add_argument('ukirt_index', help='Path to UKIRT lightcurve index')
    parser.add_argument('data_dir', help='Path to output data directory')
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    gather_data(args)
