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
    src_table_id = 'UKIRT_year' + args.year + '_field' + args.field + '_ccd' + args.ccd
    index_file = path.join(args.ukirt_index_dir, src_table_id + '.json')
    ukirt_index = utils.load_ukirt_index(index_file)

    # For each selected star
    for j, star_id in enumerate(star_list):
        star_data = var_catalog[star_id]

        if star_data['UKIRT_lc_files']:
            print('Retrieving photometry for ' + star_id
                  + ', ' + str(j) + ' out of ' + str(len(star_list)) + ' stars')

            # Check to see if this star has a pre-existing lightcurve available.
            # If so, read it in; new data will be appended to it.
            # This is done to efficiently allow UKIRT data from different years to be
            # included.  As a result, it is assumed that the OGLE data will already
            # be in the file.
            if check_lightcurve_exists(args, star_id):
                hdr, photometry = utils.load_multiband_lc(utils.get_lc_path(args, star_id))

            else:
                hdr = utils.make_lc_header(star_id, star_data)

                # Harvest OGLE timeseries photometry from online archive
                ogleI, ogleV = utils.fetch_ogle_photometry(star_id, star_data)
                photometry = {
                     'LC_I': ogleI,
                     'LC_V': ogleV
                }

            # Fetch the UKIRT lightcurve for the star for the working index file
            # (Note this may or may not be in addition to UKIRT data from other years of
            # the survey if the lightcurve preexisted)
            if not check_ukirt_data_included(hdr, index_file):
                ukirtH, ukirtK = utils.fetch_ukirt_photometry(star_id, star_data, ukirt_index, src_table_id)

                # Record the UKIRT index file, identifying the extension number for
                # UKIRT data from this file
                uid, hdr = record_ukirt_index(hdr, index_file)
                photometry['LC_H'+str(uid)] = ukirtH
                photometry['LC_K'+str(uid)] = ukirtK

            # Output the (updated) lightcurve
            utils.output_multiband_lc(args, star_id, star_data, hdr, photometry)

def check_ukirt_data_included(hdr, index_file):
    """
    Function to check that data from the working UKIRT index file hasn't already
    been included in the lightcurve.

    :param hdr: FITS header object
    :param index_file: path to the UKIRT lightcurve index file
    :return: Boolean
    """
    status = False

    for key, value in hdr.items():
        if 'UKSRC' in key and path.basename(index_file) in value:
            status = True

    return status

def record_ukirt_index(hdr, index_file):
    """
    Function to record the working ukirt index in the lightcurve header.
    An existing lightcurve may already have UKIRT data in H and K bands,
    labelled with an index.  This function returns the highest index + 1,
    and records the corresponding UKIRT source file in the header.

    :param hdr: FITS header object
    :param ukirt_index: Index of UKIRT lightcurves
    :return:
    """

    # Find the highest index number of any existing UKIRT source table
    # header keywords
    uid = -1
    for key, value in hdr.items():
        if 'UKSRC' in key:
            if int(key[-1:]) > uid:
                uid = int(key[-1:])

    # Increment this index and store the current working index file
    uid += 1
    hdr['UKSRC'+str(uid)] = path.basename(index_file)

    return uid, hdr

def check_lightcurve_exists(args, star_id):
    """
    Function to check for a preexisting lightcurve file for a given star.

    :param args: Program commandline arguments
    :param star_id: Identifier name for the star, used to find the lightcurve file
    :return:
    :param new_lc: Boolean, True if existing lightcurve file is available
    """

    file_path = utils.get_lc_path(args, star_id)

    return path.isfile(file_path)

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
