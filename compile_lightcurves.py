# Program to gather and compile optical and NIR timeseries photometry for
# a catalog of known variable stars, using the OGLE and UKIRT photometry catalogs

import argparse
import utils
from os import path

def gather_data(args):
    """
    Function to gather and compile optical and NIR timeseries photometry for
    a catalog of known variable stars, using the OGLE and UKIRT photometry catalogs

    :param args: Program arguments
    :return: None, outputs a set of lightcurve files
    """

    # Load catalog of variable stars
    var_catalog = utils.load_json_catalog(args.var_catalog_file)

    # The UKIRT catalogs are huge so loading the lightcurve index wholesale is
    # too slow to be practical.  Instead, we sort the var_catalog to group the
    # stars that are recorded in the same UKIRT catalog (i.e. same year, same field).
    # Then we extract the lightcurves for each set of stars.
    star_list = sort_stars_by_ukirt_catalog(args, var_catalog)
    #print(star_list)

    # Load UKIRT lightcurve index
    index_file = path.join(
        args.ukirt_index_dir,
        'UKIRT_year' + args.year + '_field' + args.field + '_ccd' + args.ccd + '.json'
    )
    ukirt_index = utils.load_ukirt_index(index_file)

    # For each selected star
    for j, star_id in enumerate(star_list):
        star_data = var_catalog[star_id]

        if star_data['UKIRT_lc_files']:
            print('Retrieving photometry for ' + star_id
                  + ', ' + str(j) + ' out of ' + str(len(star_list)))

            # Harvest OGLE timeseries photometry from online archive
            ogleI, ogleV = utils.fetch_ogle_photometry(star_id, star_data)

            # Fetch the UKIRT lightcurve for the star
            ukirtH, ukirtK = utils.fetch_ukirt_photometry(star_id, star_data, ukirt_index)

            # Output combined lightcurve
            photometry = {
                 'I': ogleI,
                 'V': ogleV
            }
            for k,lc in enumerate(ukirtH):
                photometry['H'+str(k)] = lc
            for k,lc in enumerate(ukirtK):
                photometry['K'+str(k)] = lc
            utils.output_multiband_lc(args, star_id, star_data, photometry)

def sort_stars_by_ukirt_catalog(args, var_catalog):
    """
    Function to group stars in the variable star catalog according to which
    UKIRT lightcurve index they fall into.  This is done to avoid loading very large
    UKIRT indices repeatedly into memory.

    :param args: program arguments
    :param var_catalog:
    :return:
    """

    # There is a lightcurve index file corresponding to each of the UKIRT source tables.
    # The var_catalog includes the source table ID for each star, so we can use
    # this to identify all stars from the UKIRT index requested
    source_table_id = 'UKIRT_year' + args.year + '_field' + args.field + '_ccd' + args.ccd + '_md.tbl'

    # Capture the identifiers of all stars which have the source_table_id listed in
    # their UKIRT_source_table parameter
    star_list = [key for key,data in var_catalog.items()
                 if data['UKIRT_source_table'] and source_table_id in data['UKIRT_source_table']]

    return star_list

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('var_catalog_file', help='Path to variable input catalog')
    parser.add_argument('ukirt_index_dir', help='Path to top-level UKIRT lightcurve index directory')
    parser.add_argument('year', help='Year of UKIRT survey')
    parser.add_argument('field', help='Field of survey')
    parser.add_argument('ccd', help='CCD in survey')
    parser.add_argument('output_dir', help='Path to output data directory')
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    gather_data(args)
